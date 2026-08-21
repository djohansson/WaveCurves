/* wc_curve.c -- wave curve discretisation and time stepping (paper Section 4). */

#include "wc_curve.h"

#include <stdlib.h>
#include <string.h>

#define WC_MAX_CURVE_POINTS 4096

/* ------------------------------------------------------------------------- */

void wc_params_default(WcParams *p)
{
    memset(p, 0, sizeof(*p));

    p->spacing    = 0.045f;
    p->radius0    = 0.20f;          /* the paper's value for all its examples */
    p->radius_min = 0.06f;
    p->radius_max = 0.55f;
    p->substeps   = 1;

    p->min_steepness = 0.01f;       /* Section 4.2 */
    p->max_steepness = 0.70f;       /* Appendix C  */
    p->max_points    = 24000;       /* the paper's production runs use ~1e6    */

    p->seeding_on      = true;
    p->seed_candidates = 900;
    p->seed_gain       = 0.50f;
    p->seed_ramp_time  = 0.5f;      /* Section 4.3 */
    /* beta(k) = seed_spectrum * rho * gstar / k. The paper leaves the constant
     * as a tuning knob. It matters more than it looks: linear wave theory only
     * holds while steepness is small, and because a few dozen stripes overlap
     * at any point the *total* steepness of Eq. (29) grows with the curve
     * budget. This value puts a typical strain-driven seed near a steepness of
     * 0.03 after its 0.5 s ramp, which keeps the summed steepness comfortably
     * under the critical value s_c rather than pinned against it. */
    p->seed_spectrum   = 6.0e-6f;
    /* Undirected seeding, off by default: see the long note in seed_curves().
     * Raise it (~0.05) together with max_points for a livelier background. */
    p->seed_isotropic  = 0.0f;
    p->seed_length     = 1.1f;
    p->seed_align      = 2.0f;      /* alpha of Eq. (28) */

    p->n_bands      = 4;
    p->wavelength[0] = 0.05f;
    p->wavelength[1] = 0.10f;
    p->wavelength[2] = 0.20f;
    p->wavelength[3] = 0.40f;

    p->critical_steepness = 3.0f;   /* s_c of Eq. (29) */
}

/* ---------------------------------------------------------------------------
 * Scratch buffers. Single-threaded by design: the solver runs on one thread and
 * only the splatting stage (wc_render.c) is parallelised.
 * ------------------------------------------------------------------------- */
static WcPoint *g_pt_scratch;   static int g_pt_scratch_cap;
static float   *g_f_scratch;    static int g_f_scratch_cap;

/* Kept separate so growing one never invalidates a live pointer into the other. */
static void scratch_pt(int n)
{
    if (n <= g_pt_scratch_cap) return;
    int cap = g_pt_scratch_cap ? g_pt_scratch_cap : 256;
    while (cap < n) cap *= 2;
    g_pt_scratch = (WcPoint *)realloc(g_pt_scratch, (size_t)cap * sizeof(WcPoint));
    g_pt_scratch_cap = cap;
}

static void scratch_f(int n)
{
    if (n <= g_f_scratch_cap) return;
    int cap = g_f_scratch_cap ? g_f_scratch_cap : 256;
    while (cap < n) cap *= 2;
    g_f_scratch = (float *)realloc(g_f_scratch, (size_t)cap * sizeof(float));
    g_f_scratch_cap = cap;
}

/* ---------------------------------------------------------------------------
 * Curve storage. Removed curves are swapped to the tail rather than freed so
 * their point arrays get recycled by the next spawn; over a few seconds the
 * simulation stops calling the allocator entirely.
 * ------------------------------------------------------------------------- */

static void curve_reserve(WcCurve *c, int n)
{
    if (n <= c->cap) return;
    int cap = c->cap ? c->cap : 32;
    while (cap < n) cap *= 2;
    c->pt = (WcPoint *)realloc(c->pt, (size_t)cap * sizeof(WcPoint));
    c->cap = cap;
}

static WcPoint *curve_push(WcCurve *c)
{
    curve_reserve(c, c->n + 1);
    WcPoint *p = &c->pt[c->n++];
    memset(p, 0, sizeof(*p));
    p->alive  = true;
    p->fade   = 1.0f;
    p->radius = 0.2f;
    return p;
}

static WcCurve *sim_new_curve(WcSim *s)
{
    if (s->n_curve >= s->cap_curve) {
        int cap = s->cap_curve ? s->cap_curve * 2 : 128;
        s->curve = (WcCurve *)realloc(s->curve, (size_t)cap * sizeof(WcCurve));
        memset(&s->curve[s->cap_curve], 0,
               (size_t)(cap - s->cap_curve) * sizeof(WcCurve));
        s->cap_curve = cap;
    }
    WcCurve *c = &s->curve[s->n_curve++];   /* may still own a recycled buffer */
    c->n      = 0;
    c->closed = false;
    c->id     = s->next_id++;
    return c;
}

/* Sign of cross(tangent, khat) at point i, used to seed and to test p->orient. */
static float point_handedness(const WcCurve *c, int i)
{
    int im = (i > 0) ? i - 1 : (c->closed ? c->n - 1 : i);
    int ip = (i < c->n - 1) ? i + 1 : (c->closed ? 0 : i);
    Vec2 T = v2_sub(c->pt[ip].x, c->pt[im].x);
    if (v2_len2(T) < 1e-14f) return 0.0f;
    return v2_cross(v2_norm(T), v2_norm(c->pt[i].k));
}

static void sim_remove_curve(WcSim *s, int i)
{
    /* Swap-with-last, preserving both allocations. */
    WcCurve tmp = s->curve[i];
    s->curve[i] = s->curve[s->n_curve - 1];
    s->curve[s->n_curve - 1] = tmp;
    s->curve[s->n_curve - 1].n = 0;
    s->n_curve--;
}

/* ------------------------------------------------------------------------- */

void wc_sim_init(WcSim *s, WcBase *base, const WcWater *w, const WcParams *p, uint32_t seed)
{
    memset(s, 0, sizeof(*s));
    s->base  = base;
    s->water = *w;
    s->p     = *p;
    wc_rng_seed(&s->rng, seed);
}

void wc_sim_free(WcSim *s)
{
    for (int i = 0; i < s->cap_curve; ++i) free(s->curve[i].pt);
    free(s->curve);
    memset(s, 0, sizeof(*s));
}

void wc_sim_clear(WcSim *s)
{
    for (int i = 0; i < s->cap_curve; ++i) s->curve[i].n = 0;
    s->n_curve = 0;
    s->stat_points = s->stat_curves = 0;
}

/* ---------------------------------------------------------------------------
 * Area patches, Fig. 6.
 *
 * Each segment owns a trapezoid of width 2r at either end and length equal to
 * the segment length; half of it is credited to each endpoint. Open-curve
 * endpoints therefore get roughly half the area of an interior point, which is
 * exactly what the figure shows.
 * ------------------------------------------------------------------------- */
static void compute_areas(WcCurve *c, float *out)
{
    for (int i = 0; i < c->n; ++i) out[i] = 0.0f;
    if (c->n < 2) return;

    int nseg = c->closed ? c->n : c->n - 1;
    for (int j = 0; j < nseg; ++j) {
        int a = j, b = (j + 1) % c->n;
        float L = v2_len(v2_sub(c->pt[b].x, c->pt[a].x));
        float half = 0.5f * L * (c->pt[a].radius + c->pt[b].radius);
        out[a] += half;
        out[b] += half;
    }
}

/* ---------------------------------------------------------------------------
 * One forward-Euler step of Eqs. (15), (16), (17), (24), followed by the wave
 * action update of Eq. (25).
 * ------------------------------------------------------------------------- */
static void evolve(WcSim *s, float dt)
{
    const WcWater *w = &s->water;
    const WcParams *P = &s->p;
    const float margin = P->radius_max;
    const float wx = s->base->world_w, wy = s->base->world_h;

    for (int ci = 0; ci < s->n_curve; ++ci) {
        WcCurve *c = &s->curve[ci];
        if (c->n < 2) { c->n = 0; continue; }

        scratch_f(c->n);
        float *area_before = g_f_scratch;
        compute_areas(c, area_before);
        /* stash: compute_areas will be called again into the same buffer, so
         * move the "before" areas onto the points themselves. */
        for (int i = 0; i < c->n; ++i) c->pt[i].area = area_before[i];

        for (int i = 0; i < c->n; ++i) {
            WcPoint *p = &c->pt[i];
            if (!p->alive) continue;

            float k = v2_len(p->k);
            if (k < 1e-3f || !(k < 1e6f)) { p->alive = false; continue; }
            Vec2 khat = v2_mul(p->k, 1.0f / k);

            WcBaseSample bs;
            wc_base_sample(s->base, p->x, &bs);

            const float geff = bs.g_eff;
            const float om   = wc_omega(w, k, geff);
            const float cg   = wc_group_speed(w, k, geff);
            const Vec2  cgv  = v2_mul(khat, cg);

            /* --- Eq. (15): xdot = U + c_g ------------------------------- */
            p->x = v2_mad(p->x, v2_add(bs.U, cgv), dt);

            /* --- Eq. (16): kdot = -(domega/dg*) grad g* - [grad U]^T k --- */
            Vec2 kdot = v2_sub(v2_mul(bs.grad_geff, -wc_domega_dgeff(w, k, geff)),
                               m2_apply_t(bs.gradU, p->k));
            p->k = v2_mad(p->k, kdot, dt);

            /* --- Eq. (17): phidot = -omega + c_g . k -------------------- */
            p->phase = wc_wrap_pi(p->phase + (-om + cg * k) * dt);

            /* --- Eq. (24): radial stretching by the background flow ------
             * rdot = r * ( n . grad U . n ) with n = khat, the in-plane curve
             * normal. The paper drops the grad c_g contribution because
             * evaluating it needs grad k, i.e. a neighbourhood search. */
            float stretch = m2_quad(bs.gradU, khat);
            p->radius = wc_clampf(p->radius * (1.0f + stretch * dt),
                                  P->radius_min, P->radius_max);

            /* --- Section 4.3: ramp energy into freshly seeded points ----
             * gamma is re-evaluated every step, so a seed that drifts out of a
             * growing region simply stops gaining energy. */
            if (p->seed_ramp > 0.0f) {
                Mat2 D = m2_sym(bs.gradU);
                float gamma = wc_growth_rate(w, k, geff, khat, D, bs.dgeff_dt);
                gamma = wc_maxf(gamma, 0.0f) + P->seed_isotropic;
                float beta = wc_seed_spectrum(w, k, geff, P->seed_spectrum);
                float E = wc_energy_from_action(w, p->action, k, geff) + beta * gamma * dt;
                p->action = wc_action_from_energy(w, E, k, geff);
                p->seed_ramp -= dt;
            }

            /* --- Section 4.2: wavenumber-dependent exponential damping --- */
            float decay = wc_amp_decay(w, k, dt);
            p->action *= decay * decay;          /* action ~ amplitude^2 */

            if (p->fade < 1.0f) p->action *= p->fade;

            /* --- Reflection off solid obstacles -------------------------
             * Ray theory reflects each control point independently: push it
             * back onto the surface and mirror its wavevector about the normal.
             * A wavefront that is only partly past the obstacle therefore forms
             * a moving kink at the contact point, and the reflected arm
             * reassembles into the correct mirror image over the next few steps.
             *
             * Mirroring khat also mirrors the handedness of (tangent, khat), so
             * the fold test's reference sign has to flip with it -- otherwise
             * the whole reflected arm reads as folded and gets deleted.
             *
             * This is specular reflection only. Ray tracing does not model
             * diffraction, so waves leave a hard-edged geometric shadow behind
             * an obstacle rather than bending around it. */
            if (s->base->n_obstacle > 0) {
                Vec2 n;
                float sd = wc_base_sdf(s->base, p->x, &n);
                if (sd < 0.0f) {
                    p->x = v2_mad(p->x, n, -sd + 1.0e-4f);
                    float kn = v2_dot(p->k, n);
                    if (kn < 0.0f) {
                        p->k = v2_mad(p->k, n, -2.0f * kn);
                        p->orient = -p->orient;
                    }
                }
            }

            p->g_eff = geff;
            p->age  += dt;

            if (p->x.x < -margin || p->x.x > wx + margin ||
                p->x.y < -margin || p->x.y > wy + margin) {
                p->alive = false;
            }
        }

        /* --- Eq. (25) / errata (E2): action density times advected area is
         * invariant, so rescale by the area ratio across the step. --------- */
        scratch_f(c->n);
        float *area_after = g_f_scratch;
        compute_areas(c, area_after);
        for (int i = 0; i < c->n; ++i) {
            WcPoint *p = &c->pt[i];
            if (!p->alive) continue;
            float a0 = p->area, a1 = area_after[i];
            if (a1 > 1e-12f && a0 > 1e-12f) p->action *= a0 / a1;
            p->area = a1;
        }

        /* --- Fold detection (Section 4.2) ---------------------------------
         * The stripe is generated by offsetting the curve along khat. That map
         * inverts exactly where cross(tangent, khat) changes sign, which is
         * where the paper's "curve normal opposes the surface normal" test
         * fires. Inverted control points get their amplitude zeroed. */
        for (int i = 0; i < c->n; ++i) {
            WcPoint *p = &c->pt[i];
            if (!p->alive) continue;
            /* A point next to a reflection kink has a meaningless tangent for
             * one step; skip it rather than deleting a healthy wavefront. */
            int im = (i > 0) ? i - 1 : (c->closed ? c->n - 1 : i);
            int ip = (i < c->n - 1) ? i + 1 : (c->closed ? 0 : i);
            if (c->pt[im].orient != p->orient || c->pt[ip].orient != p->orient) continue;
            float cr = point_handedness(c, i);
            if (cr == 0.0f) continue;
            if (cr * p->orient <= 0.0f) { p->action = 0.0f; p->alive = false; }
        }

        /* --- Amplitude reconstruction, clamping and culling --------------- */
        for (int i = 0; i < c->n; ++i) {
            WcPoint *p = &c->pt[i];
            if (!p->alive) { p->amp = 0.0f; continue; }

            float k = v2_len(p->k);
            float amp = wc_amp_from_action(w, p->action, k, p->g_eff);

            /* Appendix C: cap steepness at 0.7. */
            float smax = P->max_steepness / wc_maxf(k, 1e-6f);
            if (amp > smax) {
                amp = smax;
                p->action = wc_action_from_amp(w, amp, k, p->g_eff);
            }
            p->amp = amp;

            /* Section 4.2: delete control points whose steepness falls below
             * the visibility threshold. Points still ramping in are exempt --
             * they legitimately start at zero amplitude. */
            if (p->seed_ramp <= 0.0f && amp * k < P->min_steepness) p->alive = false;
        }
    }
}

/* ---------------------------------------------------------------------------
 * Remove dead control points, splitting a curve into several where a gap opens
 * up in its middle.
 * ------------------------------------------------------------------------- */
static void split_and_compact(WcSim *s)
{
    int killed = 0;
    int n_curve_before = s->n_curve;

    for (int ci = 0; ci < n_curve_before; ++ci) {
        WcCurve *c = &s->curve[ci];
        if (c->n == 0) continue;

        int n_dead = 0;
        for (int i = 0; i < c->n; ++i) if (!c->pt[i].alive) n_dead++;
        killed += n_dead;
        if (n_dead == 0) continue;

        if (n_dead == c->n) { c->n = 0; continue; }

        /* A closed curve that loses a point becomes open. Rotate so index 0 is
         * dead, then the alive runs never straddle the seam. */
        if (c->closed) {
            int first_dead = 0;
            while (c->pt[first_dead].alive) first_dead++;
            if (first_dead > 0) {
                scratch_pt(c->n);
                memcpy(g_pt_scratch, c->pt, (size_t)c->n * sizeof(WcPoint));
                for (int i = 0; i < c->n; ++i)
                    c->pt[i] = g_pt_scratch[(i + first_dead) % c->n];
            }
            c->closed = false;
        }

        /* Copy out, then rebuild this curve from the first run and spawn new
         * curves for the rest. */
        scratch_pt(c->n);
        int n = c->n;
        memcpy(g_pt_scratch, c->pt, (size_t)n * sizeof(WcPoint));
        c->n = 0;

        int i = 0;
        bool first_run = true;
        while (i < n) {
            while (i < n && !g_pt_scratch[i].alive) i++;
            int run0 = i;
            while (i < n && g_pt_scratch[i].alive) i++;
            int run_len = i - run0;
            if (run_len < 2) continue;      /* a lone point cannot be a stripe */

            WcCurve *dst;
            if (first_run) {
                dst = c;
                first_run = false;
            } else {
                dst = sim_new_curve(s);     /* may realloc s->curve ... */
                c   = &s->curve[ci];        /* ... so refresh the pointer */
                dst = &s->curve[s->n_curve - 1];
            }
            curve_reserve(dst, run_len);
            memcpy(dst->pt, &g_pt_scratch[run0], (size_t)run_len * sizeof(WcPoint));
            dst->n = run_len;
            dst->closed = false;
        }
    }

    /* Drop empty curves. */
    for (int ci = s->n_curve - 1; ci >= 0; --ci)
        if (s->curve[ci].n < 2) sim_remove_curve(s, ci);

    s->stat_killed = killed;
}

/* ---------------------------------------------------------------------------
 * Re-sampling to a uniform control-point spacing (Section 4.2).
 *
 * The wavevector is interpolated as magnitude + direction and recomposed, per
 * the paper; phase uses shortest-arc interpolation; action and radius are
 * densities defined along the curve and interpolate linearly.
 * ------------------------------------------------------------------------- */
static WcPoint lerp_point(const WcPoint *a, const WcPoint *b, float t)
{
    WcPoint r = *a;
    r.x = v2_lerp(a->x, b->x, t);

    float ka = v2_len(a->k), kb = v2_len(b->k);
    Vec2  da = v2_norm(a->k), db = v2_norm(b->k);
    Vec2  dd = v2_lerp(da, db, t);
    if (v2_len2(dd) < 1e-8f) dd = da;      /* near-antiparallel: keep the left one */
    r.k = v2_mul(v2_norm(dd), wc_lerpf(ka, kb, t));

    r.phase     = wc_wrap_pi(wc_lerp_phase(a->phase, b->phase, t));
    r.action    = wc_lerpf(a->action, b->action, t);
    r.radius    = wc_lerpf(a->radius, b->radius, t);
    r.g_eff     = wc_lerpf(a->g_eff,  b->g_eff,  t);
    r.amp       = wc_lerpf(a->amp,    b->amp,    t);
    r.seed_ramp = wc_lerpf(a->seed_ramp, b->seed_ramp, t);
    r.age       = wc_lerpf(a->age,    b->age,    t);
    r.fade      = wc_minf(a->fade,    b->fade);
    r.orient    = (t < 0.5f) ? a->orient : b->orient;   /* never blend a sign */
    r.alive     = true;
    return r;
}

static void resample_curve(WcSim *s, WcCurve *c)
{
    const float spacing = s->p.spacing;
    if (c->n < 2) { c->n = 0; return; }

    int n = c->n;
    int nseg = c->closed ? n : n - 1;

    scratch_f(nseg + 1);
    float *seg = g_f_scratch;
    float L = 0.0f;
    for (int j = 0; j < nseg; ++j) {
        seg[j] = v2_len(v2_sub(c->pt[(j + 1) % n].x, c->pt[j].x));
        L += seg[j];
    }

    if (L < spacing * 0.75f) {
        /* Too short to carry a wavefront any more. */
        if (c->closed || L < spacing * 0.25f) { c->n = 0; return; }
        return;                                   /* leave a short open stub alone */
    }

    int m;
    float step;
    if (c->closed) {
        m = (int)(L / spacing + 0.5f);
        if (m < 8) { c->n = 0; return; }
        step = L / (float)m;
    } else {
        m = (int)(L / spacing + 0.5f) + 1;
        if (m < 2) m = 2;
        step = L / (float)(m - 1);
    }
    if (m > WC_MAX_CURVE_POINTS) { m = WC_MAX_CURVE_POINTS; step = L / (float)(m - 1); }

    /* scratch_pt() cannot invalidate `seg`, which lives in the float scratch. */
    scratch_pt(m > n ? m : n);
    memcpy(g_pt_scratch, c->pt, (size_t)n * sizeof(WcPoint));
    const WcPoint *src = g_pt_scratch;

    curve_reserve(c, m);

    int j = 0;
    float acc = 0.0f;                 /* arclength at the start of segment j */
    for (int i = 0; i < m; ++i) {
        float target = (float)i * step;
        while (j < nseg - 1 && acc + seg[j] < target) { acc += seg[j]; j++; }
        float t = (seg[j] > 1e-9f) ? (target - acc) / seg[j] : 0.0f;
        t = wc_clampf(t, 0.0f, 1.0f);
        c->pt[i] = lerp_point(&src[j], &src[(j + 1) % n], t);
    }
    c->n = m;
}

/* ---------------------------------------------------------------------------
 * Seeding (Section 4.3).
 * ------------------------------------------------------------------------- */

int wc_sim_grow_curve(WcSim *s, Vec2 seed, Vec2 khat, float wavelength)
{
    const WcParams *P = &s->p;
    const float ds    = P->spacing;
    const float half  = 0.5f * P->seed_length;
    const int   nhalf = (int)(half / ds) + 1;
    if (nhalf < 2) return 0;

    const float k0 = WC_TAU / wavelength;
    const float wx = s->base->world_w, wy = s->base->world_h;

    /* March outward from the seed in both directions along Eq. (27),
     * xdot = N x khat  (= perp(khat) on a planar domain), turning khat as we go
     * with Eq. (28), khat_dot = -alpha (I - khat khat^T) D khat. */
    scratch_pt(2 * nhalf + 1);
    WcPoint *tmp = g_pt_scratch;
    int lo = nhalf, hi = nhalf;                     /* fill outward from the middle */

    for (int dir = -1; dir <= 1; dir += 2) {
        Vec2  x  = seed;
        Vec2  kh = khat;
        int   idx = nhalf;

        for (int step = 0; step < nhalf; ++step) {
            if (step > 0) {
                WcBaseSample bs;
                wc_base_sample(s->base, x, &bs);
                Mat2 D = m2_sym(bs.gradU);
                Vec2 Dk = m2_apply(D, kh);
                Vec2 transverse = v2_sub(Dk, v2_mul(kh, v2_dot(Dk, kh)));  /* (I - kk^T) D k */
                kh = v2_norm(v2_mad(kh, transverse, -P->seed_align * ds));
                x  = v2_mad(x, v2_mul(v2_perp(kh), (float)dir), ds);
                idx += dir;
            }
            if (x.x < 0.0f || x.x > wx || x.y < 0.0f || x.y > wy) break;
            if (idx < 0 || idx > 2 * nhalf) break;
            /* Stop the march at a solid rather than threading a curve through it. */
            if (s->base->n_obstacle > 0 && wc_base_sdf(s->base, x, NULL) < 0.0f) break;

            WcPoint *p = &tmp[idx];
            memset(p, 0, sizeof(*p));
            p->alive     = true;
            p->fade      = 1.0f;
            p->x         = x;
            p->k         = v2_mul(kh, k0);
            p->phase     = 0.0f;          /* a wavefront: constant phase along s */
            p->action    = 0.0f;          /* energy ramps in over seed_ramp_time */
            p->radius    = P->radius0;
            p->seed_ramp = P->seed_ramp_time;

            if (dir < 0 && idx < lo) lo = idx;
            if (dir > 0 && idx > hi) hi = idx;
        }
    }

    int count = hi - lo + 1;
    if (count < 3) return 0;

    WcCurve *c = sim_new_curve(s);
    curve_reserve(c, count);
    memcpy(c->pt, &tmp[lo], (size_t)count * sizeof(WcPoint));
    c->n = count;
    c->closed = false;

    /* Record the reference handedness used by the fold test. */
    Vec2 T = v2_sub(c->pt[count - 1].x, c->pt[0].x);
    float cr = v2_cross(v2_norm(T), v2_norm(c->pt[0].k));
    float orient = (cr >= 0.0f) ? 1.0f : -1.0f;
    for (int i = 0; i < count; ++i) c->pt[i].orient = orient;

    return count;
}

static void seed_curves(WcSim *s, float dt)
{
    const WcParams *P = &s->p;
    s->stat_seeded = 0;
    if (!P->seeding_on || P->n_bands <= 0) return;
    if (s->stat_points > P->max_points) return;

    const WcWater *w = &s->water;
    const float wx = s->base->world_w, wy = s->base->world_h;

    int seeded = 0;
    for (int i = 0; i < P->seed_candidates; ++i) {
        Vec2 x = v2(wc_rng_range(&s->rng, 0.0f, wx), wc_rng_range(&s->rng, 0.0f, wy));
        if (s->base->n_obstacle > 0 && wc_base_sdf(s->base, x, NULL) < 0.0f) continue;

        WcBaseSample bs;
        wc_base_sample(s->base, x, &bs);
        Mat2 D = m2_sym(bs.gradU);

        int band = (int)(wc_rng_f(&s->rng) * (float)P->n_bands);
        if (band >= P->n_bands) band = P->n_bands - 1;
        float lambda = P->wavelength[band];
        float k0 = WC_TAU / lambda;

        /* The paper maximises gamma over wave direction; for the strain term
         * the optimum is the eigenvector of D with the most negative eigenvalue
         * (the most compressive direction). */
        Eig2 e = m2_eigen_sym(D);

        float strain, gravity;
        wc_growth_terms(w, k0, bs.g_eff, e.v_lo, D, bs.dgeff_dt, &strain, &gravity);

        /* A wave curve is a long, coherent wavefront, so seeding one is only
         * meaningful where the flow actually prefers a direction to align it
         * with. Two parts of Eq. (20) offer no such preference: the
         * effective-gravity term, and the isotropic constant of Section 4.3.
         * Spawning curves from those produces isolated straight ripples drifting
         * across otherwise calm water at whatever angle the sampler happened to
         * pick -- the undersampling failure the paper itself warns about in
         * Section 5 ("if a perfectly isotropic wave spectrum is undersampled by
         * our method, then the spectrum will consist of few waves in a few
         * randomly chosen directions, which may also appear unnatural"). At the
         * paper's ~1e6 control points those blend into texture; at the budgets
         * here each one reads as a line.
         *
         * So undirected growth only *modulates the rate* of direction-aligned
         * seeds; it never spawns on its own. Set seed_isotropic > 0 to opt back
         * into undirected seeding, which is worth doing once max_points is high
         * enough for the background to be densely sampled.
         *
         * The direction is meaningful only when D's eigenvalues are genuinely
         * separated. Testing the compressive eigenvalue alone is not enough:
         * pure radial convergence gives D = -c*I, where every direction is an
         * eigenvector and m2_eigen_sym must fall back to a coordinate axis. */
        const float strain_eps = 1.0e-3f;      /* [1/s] */
        bool directed = (e.hi - e.lo > strain_eps) && (strain > 0.0f);

        float g_directed  = directed ? wc_maxf(strain + gravity, 0.0f) : 0.0f;
        float g_undirected = P->seed_isotropic;

        float total = g_directed + g_undirected;
        if (total <= 0.0f) continue;
        float prob = total * P->seed_gain * dt;
        if (wc_rng_f(&s->rng) >= prob) continue;

        Vec2 khat;
        if (wc_rng_f(&s->rng) * total < g_directed) {
            khat = e.v_lo;
        } else {
            float th = wc_rng_range(&s->rng, 0.0f, WC_TAU);
            khat = v2(cosf(th), sinf(th));
        }
        if (wc_rng_u32(&s->rng) & 1u) khat = v2_mul(khat, -1.0f);  /* either sense */

        if (wc_sim_grow_curve(s, x, khat, lambda) > 0) seeded++;
    }
    s->stat_seeded = seeded;
}

/* ---------------------------------------------------------------------------
 * Global point budget (Appendix C): fade out the least steep control points
 * rather than deleting them outright, and never touch a point that was created
 * in the last fraction of a second.
 * ------------------------------------------------------------------------- */
#define WC_HIST_BINS 128

static void enforce_budget(WcSim *s)
{
    const WcParams *P = &s->p;

    /* Reset any fade from a previous step. */
    for (int ci = 0; ci < s->n_curve; ++ci)
        for (int i = 0; i < s->curve[ci].n; ++i) s->curve[ci].pt[i].fade = 1.0f;

    if (s->stat_points <= P->max_points) return;

    int hist[WC_HIST_BINS];
    memset(hist, 0, sizeof(hist));
    const float inv = (float)WC_HIST_BINS / P->max_steepness;

    for (int ci = 0; ci < s->n_curve; ++ci) {
        WcCurve *c = &s->curve[ci];
        for (int i = 0; i < c->n; ++i) {
            float st = c->pt[i].amp * v2_len(c->pt[i].k);
            int b = (int)(st * inv);
            if (b < 0) b = 0;
            if (b >= WC_HIST_BINS) b = WC_HIST_BINS - 1;
            hist[b]++;
        }
    }

    int excess = s->stat_points - P->max_points;
    excess += excess / 4 + 1;                 /* aim slightly under the budget */

    int acc = 0, bin = 0;
    for (; bin < WC_HIST_BINS && acc < excess; ++bin) acc += hist[bin];
    float threshold = (float)bin / inv;

    for (int ci = 0; ci < s->n_curve; ++ci) {
        WcCurve *c = &s->curve[ci];
        for (int i = 0; i < c->n; ++i) {
            WcPoint *p = &c->pt[i];
            if (p->age < 0.2f) continue;      /* protect fresh seeds */
            float st = p->amp * v2_len(p->k);
            if (st < threshold) p->fade = 0.55f;   /* ~5 steps to invisibility */
        }
    }
}

/* ------------------------------------------------------------------------- */

void wc_sim_step(WcSim *s, float dt)
{
    int sub = s->p.substeps > 0 ? s->p.substeps : 1;
    float h = dt / (float)sub;

    for (int it = 0; it < sub; ++it) {
        evolve(s, h);
        split_and_compact(s);

        for (int ci = 0; ci < s->n_curve; ++ci) resample_curve(s, &s->curve[ci]);
        for (int ci = s->n_curve - 1; ci >= 0; --ci)
            if (s->curve[ci].n < 2) sim_remove_curve(s, ci);

        /* Refresh statistics before seeding so the budget check is honest. */
        s->stat_points = 0;
        s->stat_max_steep = 0.0f;
        for (int ci = 0; ci < s->n_curve; ++ci) {
            WcCurve *c = &s->curve[ci];
            s->stat_points += c->n;
            for (int i = 0; i < c->n; ++i) {
                float st = c->pt[i].amp * v2_len(c->pt[i].k);
                if (st > s->stat_max_steep) s->stat_max_steep = st;
            }
        }

        enforce_budget(s);
        seed_curves(s, h);

        s->time += h;
    }

    s->stat_points = 0;
    for (int ci = 0; ci < s->n_curve; ++ci) s->stat_points += s->curve[ci].n;
    s->stat_curves = s->n_curve;
}

/* ---------------------------------------------------------------------------
 * Interactive seeding.
 * ------------------------------------------------------------------------- */

void wc_sim_spawn_ring(WcSim *s, Vec2 center, float radius, float wavelength, float amplitude)
{
    const WcParams *P = &s->p;
    float circ = WC_TAU * radius;
    int m = (int)(circ / P->spacing + 0.5f);
    if (m < 12) m = 12;
    if (m > WC_MAX_CURVE_POINTS) m = WC_MAX_CURVE_POINTS;

    WcBaseSample bs;
    wc_base_sample(s->base, center, &bs);

    float k0 = WC_TAU / wavelength;
    /* Cap the launch amplitude at the steepness limit so a fat click does not
     * immediately violate Appendix C. */
    float amp = wc_minf(amplitude, P->max_steepness / k0);

    WcCurve *c = sim_new_curve(s);
    curve_reserve(c, m);
    c->n = m;
    c->closed = true;

    for (int i = 0; i < m; ++i) {
        float th = WC_TAU * (float)i / (float)m;
        Vec2 khat = v2(cosf(th), sinf(th));      /* radially outward */
        WcPoint *p = &c->pt[i];
        memset(p, 0, sizeof(*p));
        p->alive     = true;
        p->fade      = 1.0f;
        p->x         = v2_mad(center, khat, radius);
        p->k         = v2_mul(khat, k0);
        p->phase     = 0.0f;
        p->radius    = P->radius0;
        p->g_eff     = bs.g_eff;
        p->amp       = amp;
        p->action    = wc_action_from_amp(&s->water, amp, k0, bs.g_eff);
        p->seed_ramp = 0.0f;
        /* Points launched inside a solid never existed. */
        if (s->base->n_obstacle > 0 && wc_base_sdf(s->base, p->x, NULL) < 0.0f)
            p->alive = false;
    }

    Vec2 T = v2_sub(c->pt[1].x, c->pt[m - 1].x);
    float cr = v2_cross(v2_norm(T), v2_norm(c->pt[0].k));
    float orient = (cr >= 0.0f) ? 1.0f : -1.0f;
    for (int i = 0; i < m; ++i) c->pt[i].orient = orient;
}

void wc_sim_spawn_ring_spectrum(WcSim *s, Vec2 center, float radius, float amplitude)
{
    for (int b = 0; b < s->p.n_bands; ++b) {
        float lambda = s->p.wavelength[b];
        /* Give each band the same steepness rather than the same amplitude, so
         * no single band dominates the splash. */
        float amp = amplitude * lambda / s->p.wavelength[s->p.n_bands - 1];
        wc_sim_spawn_ring(s, center, radius + 0.02f * (float)b, lambda, amp);
    }
}
