/* wc_render.c -- wave stripe splatting and surface shading (paper Section 4.4). */

#include "wc_render.h"

#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------------------
 * Lookup tables.
 *
 * The inner loop runs tens of millions of times per frame, so the two
 * transcendentals it needs are tabulated. Both are smooth and periodic /
 * bounded, so linear interpolation between 4096 samples is far below the error
 * of the first-order phase reconstruction we are feeding them.
 * ------------------------------------------------------------------------- */
#define WC_LUT_N 4096

static float g_sin_lut[WC_LUT_N + 1];   /* sin over one period               */
static float g_win_lut[WC_LUT_N + 1];   /* Psi as a function of (dist/r)^2   */
static bool  g_lut_ready = false;

static void lut_init(void)
{
    if (g_lut_ready) return;
    for (int i = 0; i <= WC_LUT_N; ++i) {
        float t = (float)i / (float)WC_LUT_N;
        g_sin_lut[i] = sinf(WC_TAU * t);
        /* Eq. (22): a raised cosine that is 1 on the curve and falls to 0 at
         * the stripe edge. Tabulated against the *squared* normalised distance
         * so the splat loop never needs a square root. */
        g_win_lut[i] = 0.5f * (cosf(WC_PI * sqrtf(t)) + 1.0f);
    }
    g_lut_ready = true;
}

static inline float lut_sin(float phase)
{
    float t = phase * (1.0f / WC_TAU);
    t -= floorf(t);                       /* -> [0,1) */
    float f = t * (float)WC_LUT_N;
    int   i = (int)f;
    float a = f - (float)i;
    if (i < 0) i = 0; else if (i >= WC_LUT_N) i = WC_LUT_N - 1;
    return g_sin_lut[i] + (g_sin_lut[i + 1] - g_sin_lut[i]) * a;
}

static inline float lut_window(float u)   /* u = (dist/r)^2 in [0,1] */
{
    float f = u * (float)WC_LUT_N;
    int   i = (int)f;
    if (i < 0) return 1.0f;
    if (i >= WC_LUT_N) return 0.0f;
    float a = f - (float)i;
    return g_win_lut[i] + (g_win_lut[i + 1] - g_win_lut[i]) * a;
}

const char *wc_view_name(WcViewMode m)
{
    switch (m) {
        case WC_VIEW_SHADED:    return "shaded";
        case WC_VIEW_HEIGHT:    return "displacement";
        case WC_VIEW_STEEPNESS: return "steepness";
        case WC_VIEW_GROWTH:    return "growth rate";
        case WC_VIEW_FLOW:      return "base flow";
        default:                return "?";
    }
}

/* ------------------------------------------------------------------------- */

#define WC_DBG_N 96

void wc_field_init(WcField *f, int w, int h, float world_w, float world_h)
{
    lut_init();
    memset(f, 0, sizeof(*f));
    f->w = w; f->h = h;
    f->world_w = world_w; f->world_h = world_h;
    f->dx = world_w / (float)w;
    f->dy = world_h / (float)h;

    size_t n = (size_t)w * (size_t)h;
    f->disp     = (float *)calloc(n, sizeof(float));
    f->steep    = (float *)calloc(n, sizeof(float));
    f->height   = (float *)calloc(n, sizeof(float));
    f->pixels   = (uint32_t *)calloc(n, sizeof(uint32_t));
    f->obstacle = (uint8_t *)calloc(n, sizeof(uint8_t));

    f->dbg_n      = WC_DBG_N;
    f->dbg_growth = (float *)calloc((size_t)WC_DBG_N * WC_DBG_N, sizeof(float));
    f->dbg_speed  = (float *)calloc((size_t)WC_DBG_N * WC_DBG_N, sizeof(float));
}

void wc_field_free(WcField *f)
{
    free(f->disp); free(f->steep); free(f->height); free(f->pixels);
    free(f->obstacle);
    free(f->dbg_growth); free(f->dbg_speed);
    memset(f, 0, sizeof(*f));
}

/* Encode obstacle coverage: 0 = open water, otherwise 1..255 rising with depth
 * inside the solid over a short rim distance. */
void wc_field_update_obstacles(WcField *f, const WcBase *b)
{
    const float rim = 0.035f;   /* [m] */
    for (int py = 0; py < f->h; ++py) {
        const float wy = ((float)py + 0.5f) * f->dy;
        uint8_t *row = f->obstacle + (size_t)py * f->w;
        for (int px = 0; px < f->w; ++px) {
            const float wx = ((float)px + 0.5f) * f->dx;
            float sd = (b->n_obstacle > 0) ? wc_base_sdf(b, v2(wx, wy), NULL) : 1.0f;
            if (sd >= 0.0f) { row[px] = 0; continue; }
            float t = wc_clampf(-sd / rim, 0.0f, 1.0f);
            int v = 1 + (int)(t * 254.0f);
            row[px] = (uint8_t)v;
        }
    }
}

void wc_field_clear_band(WcField *f, int y0, int y1)
{
    size_t n = (size_t)(y1 - y0) * (size_t)f->w;
    memset(f->disp  + (size_t)y0 * f->w, 0, n * sizeof(float));
    memset(f->steep + (size_t)y0 * f->w, 0, n * sizeof(float));
}

/* ---------------------------------------------------------------------------
 * Splatting one wave stripe segment.
 *
 * For a pixel y inside the stripe we need the closest point x(s_y) on the
 * curve, which for a linear segment is the usual clamped projection. Then:
 *
 *   Eq. (21)  phi(y) ~ phi(s_y) + k(s_y) . (y - x(s_y))
 *   Eq. (23)  A(y)   ~ A(s_y) * Psi( dist(x(s_y), y), r(s_y) )
 *
 * and the segment contributes A(y) sin phi(y) to the displacement and
 * A(y) k(y) to the total steepness that Eq. (29) will use.
 * ------------------------------------------------------------------------- */
static void splat_segment(WcField *f, const WcPoint *pa, const WcPoint *pb,
                          int y0, int y1)
{
    const Vec2  a = pa->x;
    const Vec2  e = v2_sub(pb->x, a);
    const float len2 = v2_len2(e);
    if (len2 < 1e-12f) return;
    const float inv_len2 = 1.0f / len2;

    const float ra = pa->radius, rb = pb->radius;
    const float amp_a = pa->amp,  amp_b = pb->amp;
    if (amp_a <= 0.0f && amp_b <= 0.0f) return;

    const float rmax = wc_maxf(ra, rb);

    /* Segment bounding box, expanded by the stripe half-width. */
    float minx = wc_minf(a.x, pb->x.x) - rmax, maxx = wc_maxf(a.x, pb->x.x) + rmax;
    float miny = wc_minf(a.y, pb->x.y) - rmax, maxy = wc_maxf(a.y, pb->x.y) + rmax;

    int px0 = (int)(minx / f->dx);        if (px0 < 0) px0 = 0;
    int px1 = (int)(maxx / f->dx) + 1;    if (px1 > f->w) px1 = f->w;
    int py0 = (int)(miny / f->dy);        if (py0 < y0) py0 = y0;
    int py1 = (int)(maxy / f->dy) + 1;    if (py1 > y1) py1 = y1;
    if (px0 >= px1 || py0 >= py1) return;

    const Vec2  ka = pa->k, kb = pb->k;
    const float phase_a = pa->phase;
    /* Shortest-arc phase difference, resolved once per segment. */
    const float dphase = wc_wrap_pi(pb->phase - pa->phase);

    for (int py = py0; py < py1; ++py) {
        const float wy = ((float)py + 0.5f) * f->dy;
        float *drow = f->disp  + (size_t)py * f->w;
        float *srow = f->steep + (size_t)py * f->w;

        for (int px = px0; px < px1; ++px) {
            const float wx = ((float)px + 0.5f) * f->dx;

            /* Closest point on the segment: the curve parameter s_y. */
            float t = ((wx - a.x) * e.x + (wy - a.y) * e.y) * inv_len2;
            t = t < 0.0f ? 0.0f : (t > 1.0f ? 1.0f : t);

            const float ox = wx - (a.x + e.x * t);
            const float oy = wy - (a.y + e.y * t);
            const float d2 = ox * ox + oy * oy;

            const float r = ra + (rb - ra) * t;
            const float r2 = r * r;
            if (d2 >= r2) continue;

            const float win = lut_window(d2 / r2);
            const float amp = (amp_a + (amp_b - amp_a) * t) * win;      /* Eq. (23) */
            if (amp <= 0.0f) continue;

            const float kx = ka.x + (kb.x - ka.x) * t;
            const float ky = ka.y + (kb.y - ka.y) * t;

            const float phi = phase_a + dphase * t + kx * ox + ky * oy; /* Eq. (21) */

            drow[px] += amp * lut_sin(phi);
            srow[px] += amp * sqrtf(kx * kx + ky * ky);
        }
    }
}

void wc_field_splat_band(WcField *f, const WcSim *s, int y0, int y1)
{
    const float band_lo = (float)y0 * f->dy;
    const float band_hi = (float)y1 * f->dy;

    for (int ci = 0; ci < s->n_curve; ++ci) {
        const WcCurve *c = &s->curve[ci];
        if (c->n < 2) continue;
        int nseg = c->closed ? c->n : c->n - 1;
        for (int j = 0; j < nseg; ++j) {
            const WcPoint *pa = &c->pt[j];
            const WcPoint *pb = &c->pt[(j + 1) % c->n];
            /* Cheap band reject before the per-pixel work. */
            float rmax = wc_maxf(pa->radius, pb->radius);
            float lo = wc_minf(pa->x.y, pb->x.y) - rmax;
            float hi = wc_maxf(pa->x.y, pb->x.y) + rmax;
            if (hi < band_lo || lo > band_hi) continue;
            splat_segment(f, pa, pb, y0, y1);
        }
    }
}

/* ---------------------------------------------------------------------------
 * Eq. (29): soft-limit the accumulated displacement by the total steepness.
 *
 *   eta = (s_c / s_t) tanh(s_t / s_c) * sum_i A_i sin phi_i
 *
 * The prefactor is ~1 while the surface is gentle and rolls off as 1/s_t once
 * many curves pile up in one place, which is what keeps caustics and wave
 * blocking from blowing the surface apart.
 * ------------------------------------------------------------------------- */
void wc_field_resolve_band(WcField *f, const WcSim *s, int y0, int y1)
{
    const float sc = s->p.critical_steepness;
    const WcBase *b = s->base;

    for (int py = y0; py < y1; ++py) {
        const float wy = ((float)py + 0.5f) * f->dy;
        const float *drow = f->disp  + (size_t)py * f->w;
        const float *srow = f->steep + (size_t)py * f->w;
        float       *hrow = f->height + (size_t)py * f->w;

        const uint8_t *orow = f->obstacle + (size_t)py * f->w;

        for (int px = 0; px < f->w; ++px) {
            const float wx = ((float)px + 0.5f) * f->dx;
            float base = wc_base_height(b, v2(wx, wy), b->time);
            /* Flat inside solids, so the shading normal has a clean boundary. */
            if (orow[px]) { hrow[px] = base; continue; }
            const float st = srow[px];
            float scale = 1.0f;
            if (st > 1e-4f) scale = (sc / st) * tanhf(st / sc);
            hrow[px] = base + scale * drow[px];
        }
    }
}

/* ------------------------------------------------------------------------- */

void wc_field_update_debug(WcField *f, const WcSim *s)
{
    const int n = f->dbg_n;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            Vec2 x = v2(((float)i + 0.5f) / (float)n * f->world_w,
                        ((float)j + 0.5f) / (float)n * f->world_h);
            WcBaseSample bs;
            wc_base_sample(s->base, x, &bs);
            Mat2 D = m2_sym(bs.gradU);
            Eig2 e = m2_eigen_sym(D);
            /* Peak growth rate over wave direction, at the middle band. */
            float lambda = s->p.wavelength[s->p.n_bands / 2];
            float k0 = WC_TAU / lambda;
            float g = wc_growth_rate(&s->water, k0, bs.g_eff, e.v_lo, D, bs.dgeff_dt);
            f->dbg_growth[j * n + i] = g;
            f->dbg_speed [j * n + i] = v2_len(bs.U);
        }
    }
}

static inline float dbg_sample(const float *grid, int n, float u, float v)
{
    float fx = wc_clampf(u * (float)n - 0.5f, 0.0f, (float)(n - 1));
    float fy = wc_clampf(v * (float)n - 0.5f, 0.0f, (float)(n - 1));
    int x0 = (int)fx, y0 = (int)fy;
    int x1 = (x0 + 1 < n) ? x0 + 1 : x0;
    int y1 = (y0 + 1 < n) ? y0 + 1 : y0;
    float ax = fx - (float)x0, ay = fy - (float)y0;
    float a = wc_lerpf(grid[y0 * n + x0], grid[y0 * n + x1], ax);
    float b = wc_lerpf(grid[y1 * n + x0], grid[y1 * n + x1], ax);
    return wc_lerpf(a, b, ay);
}

/* ------------------------------------------------------------------------- */

static inline uint32_t pack_rgb(float r, float g, float b)
{
    /* sRGB-ish encode, then pack as ABGR (little-endian RGBA8888). */
    r = wc_clampf(r, 0.0f, 1.0f);
    g = wc_clampf(g, 0.0f, 1.0f);
    b = wc_clampf(b, 0.0f, 1.0f);
    r = powf(r, 1.0f / 2.2f);
    g = powf(g, 1.0f / 2.2f);
    b = powf(b, 1.0f / 2.2f);
    uint32_t R = (uint32_t)(r * 255.0f + 0.5f);
    uint32_t G = (uint32_t)(g * 255.0f + 0.5f);
    uint32_t B = (uint32_t)(b * 255.0f + 0.5f);
    return 0xFF000000u | (B << 16) | (G << 8) | R;
}

void wc_field_shade_band(WcField *f, const WcSim *s, WcViewMode mode, int y0, int y1)
{
    const int w = f->w, h = f->h;
    const float inv2dx = 1.0f / (2.0f * f->dx);
    const float inv2dy = 1.0f / (2.0f * f->dy);

    /* The field is drawn orthographically from straight above, but shading it
     * with a vertical view direction gives an almost constant Fresnel term and
     * the water reads as flat paint. We therefore light it as if the camera
     * were tilted ~30 degrees off vertical, which is what makes the glints
     * sweep across the ripples the way they do in a photograph. */
    const float Vx = 0.0f,   Vy = -0.527f, Vz = 0.850f;
    const float Lx = 0.427f, Ly =  0.527f, Lz = 0.736f;
    const float Hlen = sqrtf((Lx + Vx) * (Lx + Vx) + (Ly + Vy) * (Ly + Vy) +
                             (Lz + Vz) * (Lz + Vz));
    const float Hx = (Lx + Vx) / Hlen, Hy = (Ly + Vy) / Hlen, Hz = (Lz + Vz) / Hlen;

    /* The half-vector sits ~15 degrees off vertical: flat water stays dark and
     * only facets tilted by roughly that much light up, which is what turns the
     * wave field into legible glints instead of uniform glare. */
    const float inv_dx2 = 1.0f / (f->dx * f->dy);
    const float caustic_depth = 0.030f;    /* pseudo depth for the lens term [m] */

    for (int py = y0; py < y1; ++py) {
        const float *hrow = f->height + (size_t)py * w;
        const float *hup  = f->height + (size_t)(py > 0 ? py - 1 : py) * w;
        const float *hdn  = f->height + (size_t)(py < h - 1 ? py + 1 : py) * w;
        const float *srow = f->steep  + (size_t)py * w;
        const uint8_t *orow = f->obstacle + (size_t)py * w;
        uint32_t    *out  = f->pixels + (size_t)py * w;

        for (int px = 0; px < w; ++px) {
            int xm = px > 0 ? px - 1 : px;
            int xp = px < w - 1 ? px + 1 : px;

            /* Solids are opaque: a lit top face with a darker wet rim, so the
             * reflecting boundary is unambiguous in every view mode. */
            if (orow[px]) {
                float d = (float)orow[px] * (1.0f / 255.0f);
                float shade = 0.30f + 0.70f * d;
                float ao = wc_lerpf(0.35f, 1.0f, d);
                out[px] = pack_rgb(0.115f * shade * ao,
                                   0.125f * shade * ao,
                                   0.140f * shade * ao);
                continue;
            }

            switch (mode) {
            case WC_VIEW_SHADED: {
                float hc   = hrow[px];
                float dhdx = (hrow[xp] - hrow[xm]) * inv2dx;
                float dhdy = (hdn[px]  - hup[px])  * inv2dy;
                float nl = 1.0f / sqrtf(dhdx * dhdx + dhdy * dhdy + 1.0f);
                float nx = -dhdx * nl, ny = -dhdy * nl, nz = nl;

                /* Schlick Fresnel against the tilted view direction. */
                float ndv = nx * Vx + ny * Vy + nz * Vz;
                if (ndv < 0.0f) ndv = 0.0f;
                float c = 1.0f - ndv, c2 = c * c;
                float fres = 0.02f + 0.98f * c2 * c2 * c;

                /* Mirror the view about the normal and look the reflected ray
                 * up in a two-stop sky gradient. */
                float rx = 2.0f * ndv * nx - Vx;
                float ry = 2.0f * ndv * ny - Vy;
                float rz = 2.0f * ndv * nz - Vz;
                float t = wc_clampf(rz, 0.0f, 1.0f);
                t = t * t * (3.0f - 2.0f * t);
                float sr = wc_lerpf(0.72f, 0.24f, t);
                float sg = wc_lerpf(0.80f, 0.44f, t);
                float sb = wc_lerpf(0.92f, 0.78f, t);

                /* Water body: dark, slightly translucent, plus a cheap caustic
                 * term. Light refracting through the surface focuses in concave
                 * patches and defocuses in convex ones, which to first order is
                 * -laplacian(h); it is what gives real top-down water its
                 * bright filaments between the glints. */
                float lap = (hrow[xp] + hrow[xm] + hup[px] + hdn[px] - 4.0f * hc)
                            * inv_dx2;
                float caustic = wc_clampf(1.0f / (1.0f + caustic_depth * lap),
                                          0.35f, 2.6f);
                float diff = wc_maxf(nx * Lx + ny * Ly + nz * Lz, 0.0f);
                float br = (0.0015f + 0.010f * diff) * caustic;
                float bg = (0.0290f + 0.045f * diff) * caustic;
                float bb = (0.0620f + 0.080f * diff) * caustic;

                float ndh = nx * Hx + ny * Hy + nz * Hz;
                ndh = ndh > 0.0f ? ndh : 0.0f;
                float p2 = ndh * ndh, p4 = p2 * p2, p8 = p4 * p4;
                float p16 = p8 * p8, p32 = p16 * p16, p64 = p32 * p32;
                float p128 = p64 * p64, p256 = p128 * p128;
                float spec = p256 * 5.0f       /* tight sun disc  */
                           + p32  * 0.05f;     /* broad sheen     */

                /* Sun visible directly in the reflected ray adds the sparkle
                 * that a pure Blinn lobe misses at grazing slopes. */
                float rdl = rx * Lx + ry * Ly + rz * Lz;
                if (rdl > 0.985f) spec += (rdl - 0.985f) * (1.0f / 0.015f) * 6.0f;

                float r = wc_lerpf(br, sr * 0.9f, fres) + spec;
                float g = wc_lerpf(bg, sg * 0.9f, fres) + spec;
                float b = wc_lerpf(bb, sb * 0.9f, fres) + spec;

                /* Reinhard-ish highlight rolloff so glints do not clip flat. */
                r = r / (1.0f + 0.55f * r);
                g = g / (1.0f + 0.55f * g);
                b = b / (1.0f + 0.55f * b);
                out[px] = pack_rgb(r * 1.55f, g * 1.55f, b * 1.55f);
                break;
            }
            case WC_VIEW_HEIGHT: {
                float v = wc_clampf(f->disp[(size_t)py * w + px] * 90.0f, -1.0f, 1.0f);
                float p = wc_maxf(v, 0.0f), n = wc_maxf(-v, 0.0f);
                out[px] = pack_rgb(0.05f + 0.95f * p, 0.06f + 0.35f * (p + n),
                                   0.10f + 0.90f * n);
                break;
            }
            case WC_VIEW_STEEPNESS: {
                float v = wc_clampf(srow[px] / (3.0f * s->p.critical_steepness), 0.0f, 1.0f);
                out[px] = pack_rgb(v * 1.15f, v * v * 0.9f, 0.08f + 0.25f * v * v * v);
                break;
            }
            case WC_VIEW_GROWTH: {
                float u = ((float)px + 0.5f) / (float)w;
                float vv = ((float)py + 0.5f) / (float)h;
                float g = dbg_sample(f->dbg_growth, f->dbg_n, u, vv);
                float p = wc_clampf(g * 0.5f, 0.0f, 1.0f);
                float n = wc_clampf(-g * 0.5f, 0.0f, 1.0f);
                out[px] = pack_rgb(0.04f + 0.96f * p, 0.05f + 0.30f * (p + n),
                                   0.10f + 0.75f * n);
                break;
            }
            case WC_VIEW_FLOW:
            default: {
                float u = ((float)px + 0.5f) / (float)w;
                float vv = ((float)py + 0.5f) / (float)h;
                float sp = dbg_sample(f->dbg_speed, f->dbg_n, u, vv);
                float v = wc_clampf(sp / 1.6f, 0.0f, 1.0f);
                out[px] = pack_rgb(0.05f + 0.55f * v * v, 0.10f + 0.75f * v,
                                   0.22f + 0.65f * sqrtf(v));
                break;
            }
            }
        }
    }
}

void wc_field_render(WcField *f, const WcSim *s, WcViewMode mode)
{
    wc_field_update_obstacles(f, s->base);
    wc_field_clear_band(f, 0, f->h);
    wc_field_splat_band(f, s, 0, f->h);
    wc_field_resolve_band(f, s, 0, f->h);
    if (mode == WC_VIEW_GROWTH || mode == WC_VIEW_FLOW) wc_field_update_debug(f, s);
    wc_field_shade_band(f, s, mode, 0, f->h);
}

float wc_field_peak(const WcField *f)
{
    float m = 0.0f;
    size_t n = (size_t)f->w * (size_t)f->h;
    for (size_t i = 0; i < n; ++i) {
        float a = f->disp[i] < 0.0f ? -f->disp[i] : f->disp[i];
        if (a > m) m = a;
    }
    return m;
}
