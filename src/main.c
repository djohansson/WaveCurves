/* main.c -- SDL front-end for the Wave Curves reference implementation.
 *
 * Everything physical lives in wc_*.c; this file only owns the window, the
 * input handling, a band-parallel thread pool for the splatting stage, and a
 * headless mode that renders N frames to PNG for regression checking.
 */

#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "wc_base.h"
#include "wc_curve.h"
#include "wc_render.h"
#include "wc_png.h"

/* ------------------------------------------------------------------------- */

typedef struct {
    int         field_size;
    int         window_size;
    float       world;
    WcSceneId   scene;
    uint32_t    seed;
    int         max_points;
    bool        seeding;
    bool        headless;
    int         frames;
    int         dump_every;
    const char *out_prefix;
    float       dt;
    int         view;
    int         obstacles;
} Options;

static void options_default(Options *o)
{
    o->field_size  = 720;
    o->window_size = 900;
    o->world       = 4.0f;
    o->scene       = WC_SCENE_PADDLE;
    o->seed        = 12345u;
    o->max_points  = 0;          /* 0 = keep the wc_params_default value */
    o->seeding     = true;
    o->headless    = false;
    o->frames      = 0;
    o->dump_every  = 0;
    o->out_prefix  = "wavecurves";
    o->dt          = 1.0f / 60.0f;
    o->view        = WC_VIEW_SHADED;
    o->obstacles   = 0;
}

static void usage(const char *argv0)
{
    printf(
"Wave Curves -- reference implementation (Skrivan et al., SIGGRAPH 2020)\n"
"\n"
"usage: %s [options]\n"
"  --scene N          0=paddle 1=river 2=still        (default 0)\n"
"  --field N          displacement grid resolution    (default 720)\n"
"  --window N         window size in pixels           (default 900)\n"
"  --world W          domain edge length in metres    (default 4.0)\n"
"  --points N         wave curve control point budget (default 24000)\n"
"  --seed N           RNG seed                        (default 12345)\n"
"  --no-seeding       disable automatic wave seeding\n"
"  --dt S             fixed time step in seconds      (default 1/60)\n"
"  --view N           0=shaded 1=displacement 2=steepness 3=growth 4=flow\n"
"  --obstacles N      solid objects to reflect from:                (default 0)\n"
"                       0=none 1=pillars 2=harbour 3=wall with a gap\n"
"  --headless         no window; simulate and dump PNGs\n"
"  --frames N         headless: number of frames to simulate\n"
"  --dump-every N     headless: write a PNG every N frames (0 = last only)\n"
"  --out PREFIX       headless: output file prefix     (default wavecurves)\n"
"\n"
"keys: space pause | 1-5 view | c curves | s seeding | r reset | n scene\n"
"      left click ripple | right click single-wavelength ripple | f12 screenshot\n",
        argv0);
}

static bool parse_options(int argc, char **argv, Options *o)
{
    for (int i = 1; i < argc; ++i) {
        const char *a = argv[i];
        #define NEXT() (++i < argc ? argv[i] : NULL)
        if      (!strcmp(a, "--scene"))      { const char *v = NEXT(); if (v) o->scene = (WcSceneId)atoi(v); }
        else if (!strcmp(a, "--field"))      { const char *v = NEXT(); if (v) o->field_size = atoi(v); }
        else if (!strcmp(a, "--window"))     { const char *v = NEXT(); if (v) o->window_size = atoi(v); }
        else if (!strcmp(a, "--world"))      { const char *v = NEXT(); if (v) o->world = (float)atof(v); }
        else if (!strcmp(a, "--points"))     { const char *v = NEXT(); if (v) o->max_points = atoi(v); }
        else if (!strcmp(a, "--seed"))       { const char *v = NEXT(); if (v) o->seed = (uint32_t)strtoul(v, NULL, 10); }
        else if (!strcmp(a, "--dt"))         { const char *v = NEXT(); if (v) o->dt = (float)atof(v); }
        else if (!strcmp(a, "--frames"))     { const char *v = NEXT(); if (v) o->frames = atoi(v); }
        else if (!strcmp(a, "--dump-every")) { const char *v = NEXT(); if (v) o->dump_every = atoi(v); }
        else if (!strcmp(a, "--out"))        { const char *v = NEXT(); if (v) o->out_prefix = v; }
        else if (!strcmp(a, "--view"))       { const char *v = NEXT(); if (v) o->view = atoi(v); }
        else if (!strcmp(a, "--obstacles"))  { const char *v = NEXT(); if (v) o->obstacles = atoi(v); }
        else if (!strcmp(a, "--no-seeding")) { o->seeding = false; }
        else if (!strcmp(a, "--headless"))   { o->headless = true; }
        else if (!strcmp(a, "-h") || !strcmp(a, "--help")) { usage(argv[0]); return false; }
        else { fprintf(stderr, "unknown option: %s\n", a); usage(argv[0]); return false; }
        #undef NEXT
    }
    if ((int)o->scene < 0 || (int)o->scene >= WC_SCENE_COUNT) o->scene = WC_SCENE_PADDLE;
    if (o->field_size < 64) o->field_size = 64;
    if (o->headless && o->frames <= 0) o->frames = 300;
    return true;
}

/* ---------------------------------------------------------------------------
 * Band-parallel thread pool.
 *
 * Each stage of the render pipeline is split into horizontal bands that touch
 * disjoint rows, so no locking is needed anywhere -- workers just claim the
 * next band with an atomic increment.
 * ------------------------------------------------------------------------- */

#define WC_MAX_WORKERS 15
#define WC_BANDS_PER_WORKER 4

typedef enum { JOB_SPLAT, JOB_RESOLVE, JOB_SHADE } JobKind;

static struct {
    SDL_Thread   *thread[WC_MAX_WORKERS];
    SDL_Semaphore *start, *done;
    SDL_AtomicInt next_band;
    int           n_workers;
    int           n_bands;
    bool          quit;

    JobKind       job;
    WcField      *field;
    const WcSim  *sim;
    WcViewMode    mode;
} g_pool;

static void pool_run_bands(void)
{
    WcField *f = g_pool.field;
    for (;;) {
        int b = SDL_AddAtomicInt(&g_pool.next_band, 1);
        if (b >= g_pool.n_bands) break;
        int y0 = (int)((int64_t)b * f->h / g_pool.n_bands);
        int y1 = (int)((int64_t)(b + 1) * f->h / g_pool.n_bands);
        switch (g_pool.job) {
        case JOB_SPLAT:
            wc_field_clear_band(f, y0, y1);
            wc_field_splat_band(f, g_pool.sim, y0, y1);
            break;
        case JOB_RESOLVE:
            wc_field_resolve_band(f, g_pool.sim, y0, y1);
            break;
        case JOB_SHADE:
            wc_field_shade_band(f, g_pool.sim, g_pool.mode, y0, y1);
            break;
        }
    }
}

static int SDLCALL worker_main(void *ud)
{
    (void)ud;
    for (;;) {
        SDL_WaitSemaphore(g_pool.start);
        if (g_pool.quit) break;
        pool_run_bands();
        SDL_SignalSemaphore(g_pool.done);
    }
    return 0;
}

static void pool_init(void)
{
    memset(&g_pool, 0, sizeof(g_pool));
    int cores = SDL_GetNumLogicalCPUCores();
    int n = cores - 1;
    if (n < 0) n = 0;
    if (n > WC_MAX_WORKERS) n = WC_MAX_WORKERS;
    g_pool.n_workers = n;
    if (n == 0) return;

    g_pool.start = SDL_CreateSemaphore(0);
    g_pool.done  = SDL_CreateSemaphore(0);
    for (int i = 0; i < n; ++i) {
        char name[32];
        SDL_snprintf(name, sizeof(name), "wc_worker_%d", i);
        g_pool.thread[i] = SDL_CreateThread(worker_main, name, NULL);
    }
}

static void pool_shutdown(void)
{
    if (g_pool.n_workers == 0) return;
    g_pool.quit = true;
    for (int i = 0; i < g_pool.n_workers; ++i) SDL_SignalSemaphore(g_pool.start);
    for (int i = 0; i < g_pool.n_workers; ++i) SDL_WaitThread(g_pool.thread[i], NULL);
    SDL_DestroySemaphore(g_pool.start);
    SDL_DestroySemaphore(g_pool.done);
    g_pool.n_workers = 0;
}

static void pool_dispatch(JobKind job, WcField *f, const WcSim *s, WcViewMode mode)
{
    g_pool.job   = job;
    g_pool.field = f;
    g_pool.sim   = s;
    g_pool.mode  = mode;
    g_pool.n_bands = (g_pool.n_workers + 1) * WC_BANDS_PER_WORKER;
    SDL_SetAtomicInt(&g_pool.next_band, 0);

    for (int i = 0; i < g_pool.n_workers; ++i) SDL_SignalSemaphore(g_pool.start);
    pool_run_bands();                                   /* the main thread helps */
    for (int i = 0; i < g_pool.n_workers; ++i) SDL_WaitSemaphore(g_pool.done);
}

/* ------------------------------------------------------------------------- */

typedef struct {
    Options  opt;
    WcWater  water;
    WcParams params;
    WcBase   base;
    WcSim    sim;
    WcField  field;

    bool       paused;
    bool       show_curves;
    WcViewMode view;
    int        frame;
    int        shot_index;
    double     ms_sim, ms_render;
} App;

static void app_reset(App *a)
{
    wc_sim_free(&a->sim);
    wc_base_init(&a->base, a->opt.scene, a->opt.world, a->opt.world, a->water.gravity);
    wc_base_set_obstacles(&a->base, a->opt.obstacles);
    wc_sim_init(&a->sim, &a->base, &a->water, &a->params, a->opt.seed);
    wc_field_update_obstacles(&a->field, &a->base);
    a->frame = 0;
}

static void app_screenshot(App *a)
{
    char path[256];
    SDL_snprintf(path, sizeof(path), "%s_%04d.png", a->opt.out_prefix, a->shot_index++);
    if (wc_png_write(path, a->field.pixels, a->field.w, a->field.h))
        SDL_Log("wrote %s", path);
    else
        SDL_Log("failed to write %s", path);
}

/* ------------------------------------------------------------------------- */

/* Distinct hue per wavelength band, for the curve overlay. */
static void band_colour(const WcSim *s, float k, Uint8 *r, Uint8 *g, Uint8 *b)
{
    float lambda = WC_TAU / wc_maxf(k, 1e-6f);
    int   best = 0;
    float bd = 1e30f;
    for (int i = 0; i < s->p.n_bands; ++i) {
        float d = fabsf(logf(lambda / s->p.wavelength[i]));
        if (d < bd) { bd = d; best = i; }
    }
    static const Uint8 pal[6][3] = {
        { 255, 214,  92 }, { 120, 230, 180 }, { 130, 180, 255 },
        { 235, 130, 220 }, { 250, 150, 110 }, { 190, 190, 190 },
    };
    int i = best % 6;
    *r = pal[i][0]; *g = pal[i][1]; *b = pal[i][2];
}

static void draw_curve_overlay(SDL_Renderer *ren, const WcSim *s,
                               const SDL_FRect *dst, float world_w, float world_h)
{
    static SDL_FPoint *pts = NULL;
    static int pts_cap = 0;

    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
    for (int ci = 0; ci < s->n_curve; ++ci) {
        const WcCurve *c = &s->curve[ci];
        if (c->n < 2) continue;
        int n = c->n + (c->closed ? 1 : 0);
        if (n > pts_cap) {
            pts_cap = n * 2;
            pts = (SDL_FPoint *)realloc(pts, (size_t)pts_cap * sizeof(SDL_FPoint));
        }
        for (int i = 0; i < n; ++i) {
            const WcPoint *p = &c->pt[i % c->n];
            pts[i].x = dst->x + p->x.x / world_w * dst->w;
            pts[i].y = dst->y + p->x.y / world_h * dst->h;
        }
        Uint8 r, g, b;
        band_colour(s, v2_len(c->pt[0].k), &r, &g, &b);
        /* Fade with steepness so ramping-in seeds read as faint. */
        float st = c->pt[c->n / 2].amp * v2_len(c->pt[c->n / 2].k);
        Uint8 alpha = (Uint8)(60.0f + 195.0f * wc_clampf(st / 0.25f, 0.0f, 1.0f));
        SDL_SetRenderDrawColor(ren, r, g, b, alpha);
        SDL_RenderLines(ren, pts, n);
    }
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_NONE);
}

static void draw_hud(SDL_Renderer *ren, const App *a, int win_w)
{
    (void)win_w;
    SDL_SetRenderScale(ren, 2.0f, 2.0f);
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);

    char line[6][160];
    SDL_snprintf(line[0], sizeof(line[0]), "wave curves  t=%6.2fs  %s  obstacles:%s%s",
                 a->sim.time, wc_scene_name(a->base.scene),
                 wc_obstacle_layout_name(a->base.obstacle_layout),
                 a->paused ? "  [PAUSED]" : "");
    SDL_snprintf(line[1], sizeof(line[1]), "curves %5d   points %6d / %d   seeded/step %d",
                 a->sim.stat_curves, a->sim.stat_points, a->sim.p.max_points,
                 a->sim.stat_seeded);
    SDL_snprintf(line[2], sizeof(line[2]), "max steepness %.3f   seeding %s   view %s",
                 a->sim.stat_max_steep, a->sim.p.seeding_on ? "on" : "off",
                 wc_view_name(a->view));
    SDL_snprintf(line[3], sizeof(line[3]), "sim %5.1f ms   render %5.1f ms",
                 a->ms_sim, a->ms_render);
    SDL_snprintf(line[4], sizeof(line[4]), "lambda %.2f %.2f %.2f %.2f m   r=%.2f m",
                 a->sim.p.wavelength[0], a->sim.p.wavelength[1],
                 a->sim.p.wavelength[2], a->sim.p.wavelength[3], a->sim.p.radius0);
    SDL_snprintf(line[5], sizeof(line[5]),
                 "click ripple | 1-5 view | c curves | s seed | o obstacles | n scene | r reset");

    /* Drop shadow then text, so it reads over bright glints. */
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 0) SDL_SetRenderDrawColor(ren, 0, 0, 0, 190);
        else           SDL_SetRenderDrawColor(ren, 235, 240, 245, 255);
        float off = pass == 0 ? 1.0f : 0.0f;
        for (int i = 0; i < 6; ++i)
            SDL_RenderDebugText(ren, 6.0f + off, 6.0f + off + (float)i * 11.0f, line[i]);
    }

    SDL_SetRenderScale(ren, 1.0f, 1.0f);
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_NONE);
}

/* ------------------------------------------------------------------------- */

static int run_headless(App *a)
{
    const Options *o = &a->opt;
    printf("headless: scene=%s frames=%d field=%d world=%.2fm\n",
           wc_scene_name(a->base.scene), o->frames, a->field.w, o->world);

    /* Seed one visible ripple so a short run still shows something even with
     * automatic seeding disabled. With obstacles present, put it off to one
     * side so the wavefront actually meets them. */
    float sx = (o->obstacles != 0) ? 0.24f : 0.5f;
    wc_sim_spawn_ring_spectrum(&a->sim, v2(o->world * sx, o->world * 0.5f),
                               0.25f, 0.008f);

    for (int i = 0; i < o->frames; ++i) {
        wc_base_advance(&a->base, o->dt);
        wc_sim_step(&a->sim, o->dt);
        a->frame++;

        bool dump = (o->dump_every > 0 && (i % o->dump_every) == 0) ||
                    (i == o->frames - 1);
        if (dump) {
            wc_field_render(&a->field, &a->sim, a->view);
            char path[256];
            SDL_snprintf(path, sizeof(path), "%s_%04d.png", o->out_prefix, i);
            if (!wc_png_write(path, a->field.pixels, a->field.w, a->field.h)) {
                fprintf(stderr, "failed to write %s\n", path);
                return 1;
            }
            printf("frame %4d  curves %5d  points %6d  max steep %.3f  "
                   "peak |eta| %.4f m  -> %s\n",
                   i, a->sim.stat_curves, a->sim.stat_points, a->sim.stat_max_steep,
                   wc_field_peak(&a->field), path);
        } else if ((i % 60) == 0) {
            printf("frame %4d  curves %5d  points %6d  max steep %.3f\n",
                   i, a->sim.stat_curves, a->sim.stat_points, a->sim.stat_max_steep);
        }
        fflush(stdout);
    }
    return 0;
}

/* ------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
    static App app;
    options_default(&app.opt);
    if (!parse_options(argc, argv, &app.opt)) return 0;

    app.water = wc_water_default();
    wc_params_default(&app.params);
    if (app.opt.max_points > 0) app.params.max_points = app.opt.max_points;
    app.params.seeding_on = app.opt.seeding;
    if (app.opt.view < 0 || app.opt.view >= WC_VIEW_COUNT) app.opt.view = WC_VIEW_SHADED;
    if (app.opt.obstacles < 0 || app.opt.obstacles >= WC_OBSTACLE_LAYOUTS)
        app.opt.obstacles = 0;

    wc_base_init(&app.base, app.opt.scene, app.opt.world, app.opt.world, app.water.gravity);
    wc_base_set_obstacles(&app.base, app.opt.obstacles);
    wc_sim_init(&app.sim, &app.base, &app.water, &app.params, app.opt.seed);
    wc_field_init(&app.field, app.opt.field_size, app.opt.field_size,
                  app.opt.world, app.opt.world);
    wc_field_update_obstacles(&app.field, &app.base);
    app.view = (WcViewMode)app.opt.view;
    app.show_curves = false;

    if (app.opt.headless) {
        int rc = run_headless(&app);
        wc_field_free(&app.field);
        wc_sim_free(&app.sim);
        return rc;
    }

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow("Wave Curves (Skrivan et al. 2020) -- reference implementation",
                                       app.opt.window_size, app.opt.window_size,
                                       SDL_WINDOW_RESIZABLE);
    if (!win) { fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError()); return 1; }

    SDL_Renderer *ren = SDL_CreateRenderer(win, NULL);
    if (!ren) { fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError()); return 1; }
    SDL_SetRenderVSync(ren, 1);

    /* pack_rgb() writes 0xAABBGGRR, i.e. bytes R,G,B,A on a little-endian host,
     * which is SDL's ABGR8888 and PNG's natural RGBA byte order at once. */
    SDL_Texture *tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ABGR8888,
                                         SDL_TEXTUREACCESS_STREAMING,
                                         app.field.w, app.field.h);
    if (!tex) { fprintf(stderr, "SDL_CreateTexture: %s\n", SDL_GetError()); return 1; }
    SDL_SetTextureScaleMode(tex, SDL_SCALEMODE_LINEAR);

    pool_init();
    SDL_Log("worker threads: %d (+1 main)", g_pool.n_workers);

    bool running = true;
    bool debug_dirty = true;

    while (running) {
        int win_w, win_h;
        SDL_GetWindowSizeInPixels(win, &win_w, &win_h);
        float side = (float)((win_w < win_h) ? win_w : win_h);
        SDL_FRect dst = { ((float)win_w - side) * 0.5f, ((float)win_h - side) * 0.5f,
                          side, side };

        SDL_Event ev;
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_EVENT_QUIT) running = false;
            else if (ev.type == SDL_EVENT_KEY_DOWN) {
                switch (ev.key.key) {
                case SDLK_ESCAPE: case SDLK_Q: running = false; break;
                case SDLK_SPACE:  app.paused = !app.paused; break;
                case SDLK_1: app.view = WC_VIEW_SHADED;    break;
                case SDLK_2: app.view = WC_VIEW_HEIGHT;    break;
                case SDLK_3: app.view = WC_VIEW_STEEPNESS; break;
                case SDLK_4: app.view = WC_VIEW_GROWTH;    debug_dirty = true; break;
                case SDLK_5: app.view = WC_VIEW_FLOW;      debug_dirty = true; break;
                case SDLK_C: app.show_curves = !app.show_curves; break;
                case SDLK_S: app.sim.p.seeding_on = !app.sim.p.seeding_on; break;
                case SDLK_R: app_reset(&app); break;
                case SDLK_N:
                    app.opt.scene = (WcSceneId)((app.opt.scene + 1) % WC_SCENE_COUNT);
                    app_reset(&app);
                    debug_dirty = true;
                    break;
                case SDLK_O:
                    /* Cycle obstacle layouts in place: the base flow changes,
                     * but the existing wave curves keep going and start
                     * bouncing off whatever just appeared. */
                    app.opt.obstacles = (app.opt.obstacles + 1) % WC_OBSTACLE_LAYOUTS;
                    wc_base_set_obstacles(&app.base, app.opt.obstacles);
                    wc_field_update_obstacles(&app.field, &app.base);
                    debug_dirty = true;
                    break;
                case SDLK_F12: app_screenshot(&app); break;
                default: break;
                }
            }
            else if (ev.type == SDL_EVENT_MOUSE_BUTTON_DOWN) {
                float u = (ev.button.x - dst.x) / dst.w;
                float v = (ev.button.y - dst.y) / dst.h;
                if (u >= 0.0f && u <= 1.0f && v >= 0.0f && v <= 1.0f) {
                    Vec2 p = v2(u * app.opt.world, v * app.opt.world);
                    if (ev.button.button == SDL_BUTTON_LEFT)
                        wc_sim_spawn_ring_spectrum(&app.sim, p, 0.18f, 0.014f);
                    else
                        wc_sim_spawn_ring(&app.sim, p, 0.18f,
                                          app.sim.p.wavelength[app.sim.p.n_bands / 2], 0.010f);
                }
            }
        }

        if (!app.paused) {
            Uint64 t0 = SDL_GetPerformanceCounter();
            wc_base_advance(&app.base, app.opt.dt);
            wc_sim_step(&app.sim, app.opt.dt);
            Uint64 t1 = SDL_GetPerformanceCounter();
            app.ms_sim = 1000.0 * (double)(t1 - t0) / (double)SDL_GetPerformanceFrequency();
            app.frame++;
            debug_dirty = true;
        }

        Uint64 r0 = SDL_GetPerformanceCounter();
        pool_dispatch(JOB_SPLAT,   &app.field, &app.sim, app.view);
        pool_dispatch(JOB_RESOLVE, &app.field, &app.sim, app.view);
        if ((app.view == WC_VIEW_GROWTH || app.view == WC_VIEW_FLOW) && debug_dirty) {
            wc_field_update_debug(&app.field, &app.sim);
            debug_dirty = false;
        }
        pool_dispatch(JOB_SHADE,   &app.field, &app.sim, app.view);
        Uint64 r1 = SDL_GetPerformanceCounter();
        app.ms_render = 1000.0 * (double)(r1 - r0) / (double)SDL_GetPerformanceFrequency();

        SDL_UpdateTexture(tex, NULL, app.field.pixels, app.field.w * 4);

        SDL_SetRenderDrawColor(ren, 12, 14, 18, 255);
        SDL_RenderClear(ren);
        SDL_RenderTexture(ren, tex, NULL, &dst);
        if (app.show_curves)
            draw_curve_overlay(ren, &app.sim, &dst, app.opt.world, app.opt.world);
        draw_hud(ren, &app, win_w);
        SDL_RenderPresent(ren);
    }

    pool_shutdown();
    SDL_DestroyTexture(tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();

    wc_field_free(&app.field);
    wc_sim_free(&app.sim);
    return 0;
}
