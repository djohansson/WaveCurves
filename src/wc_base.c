#include "wc_base.h"

#include <string.h>

/* Finite-difference steps used to build grad U, grad gstar and dgstar/dt. Real base
 * simulations would hand these over directly; differencing an analytic field is
 * the honest stand-in and exercises the same code path. */
#define WC_FD_DX 2.0e-3f   /* [m] */
#define WC_FD_DT 2.0e-3f   /* [s] */

const char *wc_scene_name(WcSceneId s)
{
    switch (s) {
        case WC_SCENE_PADDLE: return "paddle";
        case WC_SCENE_RIVER:  return "river";
        case WC_SCENE_STILL:  return "still";
        default:              return "?";
    }
}

void wc_base_init(WcBase *b, WcSceneId scene, float world_w, float world_h, float gravity)
{
    memset(b, 0, sizeof(*b));
    b->scene    = scene;
    b->world_w  = world_w;
    b->world_h  = world_h;
    b->gravity  = gravity;
    b->time     = 0.0f;

    float cx = 0.5f * world_w, cy = 0.5f * world_h;

    switch (scene) {
    case WC_SCENE_PADDLE:
        b->current = v2(0.0f, 0.0f);

        b->paddle_on    = true;
        b->paddle_pos   = v2(cx, cy);
        b->paddle_len   = 0.22f * world_w;
        b->paddle_width = 0.10f * world_w;
        b->paddle_omega = 1.1f;

        b->n_vortex = 4;
        b->vortex[0] = (WcVortex){ v2(0.25f * world_w, 0.28f * world_h),  0.55f, 0.30f, v2( 0.05f,  0.03f) };
        b->vortex[1] = (WcVortex){ v2(0.74f * world_w, 0.30f * world_h), -0.45f, 0.26f, v2(-0.04f,  0.05f) };
        b->vortex[2] = (WcVortex){ v2(0.28f * world_w, 0.75f * world_h), -0.50f, 0.28f, v2( 0.03f, -0.05f) };
        b->vortex[3] = (WcVortex){ v2(0.76f * world_w, 0.72f * world_h),  0.40f, 0.24f, v2(-0.05f, -0.02f) };

        b->n_cell = 3;
        b->cell[0] = (WcCell){ v2(0.18f * world_w, 0.52f * world_h), -0.55f, 0.45f, 0.13f };
        b->cell[1] = (WcCell){ v2(0.84f * world_w, 0.50f * world_h), -0.45f, 0.40f, 0.17f };
        b->cell[2] = (WcCell){ v2(0.50f * world_w, 0.88f * world_h),  0.40f, 0.42f, 0.11f };

        b->swell_amp[0]   = 0.030f;
        b->swell_k[0]     = v2(2.10f, 0.60f);
        b->swell_amp[1]   = 0.018f;
        b->swell_k[1]     = v2(-0.90f, 1.70f);
        break;

    case WC_SCENE_RIVER:
        b->current = v2(0.85f, 0.0f);

        b->n_vortex = 6;
        b->vortex[0] = (WcVortex){ v2(0.20f * world_w, 0.35f * world_h),  0.60f, 0.22f, v2(0.55f, 0.02f) };
        b->vortex[1] = (WcVortex){ v2(0.42f * world_w, 0.66f * world_h), -0.60f, 0.22f, v2(0.55f,-0.02f) };
        b->vortex[2] = (WcVortex){ v2(0.64f * world_w, 0.33f * world_h),  0.50f, 0.20f, v2(0.60f, 0.03f) };
        b->vortex[3] = (WcVortex){ v2(0.86f * world_w, 0.62f * world_h), -0.50f, 0.20f, v2(0.60f,-0.03f) };
        b->vortex[4] = (WcVortex){ v2(0.08f * world_w, 0.72f * world_h), -0.35f, 0.18f, v2(0.65f, 0.00f) };
        b->vortex[5] = (WcVortex){ v2(0.55f * world_w, 0.14f * world_h),  0.35f, 0.18f, v2(0.65f, 0.00f) };

        b->n_cell = 4;
        b->cell[0] = (WcCell){ v2(0.30f * world_w, 0.50f * world_h), -0.70f, 0.35f, 0.23f };
        b->cell[1] = (WcCell){ v2(0.58f * world_w, 0.42f * world_h), -0.55f, 0.30f, 0.19f };
        b->cell[2] = (WcCell){ v2(0.75f * world_w, 0.58f * world_h),  0.50f, 0.32f, 0.15f };
        b->cell[3] = (WcCell){ v2(0.12f * world_w, 0.28f * world_h), -0.45f, 0.28f, 0.27f };

        b->swell_amp[0]   = 0.022f;
        b->swell_k[0]     = v2(1.40f, 0.30f);
        b->swell_amp[1]   = 0.014f;
        b->swell_k[1]     = v2(0.40f, -1.90f);
        break;

    case WC_SCENE_STILL:
    default:
        /* Deliberately empty: a flat, motionless surface. Everything that moves
         * in this scene is a wave curve, which makes it the right place to look
         * at dispersion in isolation. */
        break;
    }

    /* Deep-water swell: omega = sqrt(g |k|). */
    for (int i = 0; i < 2; ++i) {
        float k = v2_len(b->swell_k[i]);
        b->swell_omega[i] = (k > 0.0f) ? sqrtf(gravity * k) : 0.0f;
    }
}

/* ---------------------------------------------------------------------------
 * Obstacles.
 * ------------------------------------------------------------------------- */

const char *wc_obstacle_layout_name(int layout)
{
    switch (layout) {
        case 0:  return "none";
        case 1:  return "pillars";
        case 2:  return "harbour";
        case 3:  return "slit";
        default: return "?";
    }
}

/* Exact SDF + outward normal of an axis-aligned box in its own frame. */
static float sdf_box_local(Vec2 p, Vec2 half, Vec2 *n)
{
    float sx = (p.x < 0.0f) ? -1.0f : 1.0f;
    float sy = (p.y < 0.0f) ? -1.0f : 1.0f;
    Vec2  q  = v2(fabsf(p.x) - half.x, fabsf(p.y) - half.y);

    if (q.x > 0.0f || q.y > 0.0f) {                 /* outside */
        Vec2 m = v2(wc_maxf(q.x, 0.0f), wc_maxf(q.y, 0.0f));
        float l = v2_len(m);
        if (l > 1e-9f) *n = v2(m.x / l * sx, m.y / l * sy);
        else           *n = v2(sx, 0.0f);
        return l + wc_minf(wc_maxf(q.x, q.y), 0.0f);
    }
    if (q.x > q.y) { *n = v2(sx, 0.0f); return q.x; }   /* inside: nearest face */
    *n = v2(0.0f, sy);
    return q.y;
}

float wc_base_sdf(const WcBase *b, Vec2 x, Vec2 *out_normal)
{
    float best = 1.0e30f;
    Vec2  bn   = v2(0.0f, 1.0f);

    for (int i = 0; i < b->n_obstacle; ++i) {
        const WcObstacle *o = &b->obstacle[i];
        float d;
        Vec2  n;
        if (o->kind == WC_OBS_CIRCLE) {
            Vec2 r = v2_sub(x, o->pos);
            float l = v2_len(r);
            d = l - o->radius;
            n = (l > 1e-9f) ? v2_mul(r, 1.0f / l) : v2(0.0f, 1.0f);
        } else {
            float c = cosf(o->angle), s = sinf(o->angle);
            Vec2 r = v2_sub(x, o->pos);
            Vec2 lp = v2(c * r.x + s * r.y, -s * r.x + c * r.y);   /* world -> local */
            Vec2 ln;
            d = sdf_box_local(lp, o->half, &ln);
            n = v2(c * ln.x - s * ln.y, s * ln.x + c * ln.y);      /* local -> world */
        }
        if (d < best) { best = d; bn = n; }
    }

    if (out_normal) *out_normal = bn;
    return best;
}

void wc_base_set_obstacles(WcBase *b, int layout)
{
    const float W = b->world_w, H = b->world_h;
    b->obstacle_layout = layout;
    b->n_obstacle = 0;

    #define ADD_CIRCLE(px, py, rr) do {                                        \
        WcObstacle *o = &b->obstacle[b->n_obstacle++];                         \
        o->kind = WC_OBS_CIRCLE; o->pos = v2(px, py); o->radius = (rr);        \
        o->half = v2(0.0f, 0.0f); o->angle = 0.0f;                             \
    } while (0)
    #define ADD_BOX(px, py, hx, hy, ang) do {                                   \
        WcObstacle *o = &b->obstacle[b->n_obstacle++];                          \
        o->kind = WC_OBS_BOX; o->pos = v2(px, py); o->radius = 0.0f;            \
        o->half = v2(hx, hy); o->angle = (ang);                                 \
    } while (0)

    switch (layout) {
    case 1:  /* three pillars of different size: circular reflection + shadowing */
        ADD_CIRCLE(0.30f * W, 0.34f * H, 0.11f * W);
        ADD_CIRCLE(0.70f * W, 0.62f * H, 0.075f * W);
        ADD_CIRCLE(0.36f * W, 0.76f * H, 0.045f * W);
        break;

    case 2:  /* harbour: an angled sea wall plus two round piles */
        ADD_BOX(0.68f * W, 0.30f * H, 0.30f * W, 0.030f * W, -0.55f);
        ADD_BOX(0.22f * W, 0.68f * H, 0.20f * W, 0.030f * W,  0.42f);
        ADD_CIRCLE(0.78f * W, 0.80f * H, 0.070f * W);
        ADD_CIRCLE(0.14f * W, 0.24f * H, 0.055f * W);
        break;

    case 3:  /* a wall with a gap: the classic two-slit / shadow-edge test */
        ADD_BOX(0.50f * W, 0.50f * H, 0.021f * W, 0.190f * H, 0.0f);
        ADD_BOX(0.50f * W, 0.94f * H, 0.021f * W, 0.190f * H, 0.0f);
        ADD_BOX(0.50f * W, 0.06f * H, 0.021f * W, 0.190f * H, 0.0f);
        break;

    default:
        break;
    }

    #undef ADD_CIRCLE
    #undef ADD_BOX
}

void wc_base_advance(WcBase *b, float dt)
{
    b->time += dt;
    b->paddle_angle += b->paddle_omega * dt;

    for (int i = 0; i < b->n_vortex; ++i) {
        WcVortex *v = &b->vortex[i];
        v->pos = v2_mad(v->pos, v->drift, dt);
        /* Wrap so the scenes run indefinitely. */
        if (v->pos.x < 0.0f)          v->pos.x += b->world_w;
        if (v->pos.x > b->world_w)    v->pos.x -= b->world_w;
        if (v->pos.y < 0.0f)          v->pos.y += b->world_h;
        if (v->pos.y > b->world_h)    v->pos.y -= b->world_h;
    }
}

/* ------------------------------------------------------------------------- */

/* Closest point on segment [a,b] to p, plus the parameter along it. */
static Vec2 closest_on_segment(Vec2 a, Vec2 bb, Vec2 p, float *out_t)
{
    Vec2 e = v2_sub(bb, a);
    float len2 = v2_len2(e);
    float t = (len2 > 1e-12f) ? wc_clampf(v2_dot(v2_sub(p, a), e) / len2, 0.0f, 1.0f) : 0.0f;
    if (out_t) *out_t = t;
    return v2_mad(a, e, t);
}

Vec2 wc_base_velocity(const WcBase *b, Vec2 x, float t)
{
    Vec2 U = b->current;

    /* Lamb-Oseen style vortices: divergence free, finite at the core. */
    for (int i = 0; i < b->n_vortex; ++i) {
        const WcVortex *v = &b->vortex[i];
        /* vortex positions are integrated in wc_base_advance; for the tiny
         * finite-difference-in-time offsets we extrapolate along the drift. */
        Vec2 vp = v2_mad(v->pos, v->drift, t - b->time);
        Vec2 d  = v2_sub(x, vp);
        float r2 = v2_len2(d);
        float a2 = v->core * v->core;
        if (r2 < 1e-10f) continue;
        float f = v->strength * (1.0f - expf(-r2 / a2)) / r2;
        U = v2_mad(U, v2_perp(d), f);
    }

    /* Converging / diverging surface cells. */
    for (int i = 0; i < b->n_cell; ++i) {
        const WcCell *c = &b->cell[i];
        Vec2 d = v2_sub(x, c->pos);
        float r2 = v2_len2(d);
        float s2 = c->radius * c->radius;
        float amp = c->strength * (0.6f + 0.4f * sinf(WC_TAU * c->pulse_hz * t));
        U = v2_mad(U, d, amp * expf(-r2 / s2));
    }

    /* Rotating paddle: the surface is dragged with the bar, falling off over a
     * boundary-layer width. Strongly compressive on the leading face and
     * strongly extensional on the trailing one -- exactly the anisotropic strain
     * that Eq. (20) rewards. */
    if (b->paddle_on) {
        float ang = b->paddle_angle + b->paddle_omega * (t - b->time);
        Vec2 dir = v2(cosf(ang), sinf(ang));
        Vec2 a   = v2_mad(b->paddle_pos, dir, -b->paddle_len);
        Vec2 bb  = v2_mad(b->paddle_pos, dir,  b->paddle_len);
        float seg_t;
        Vec2 cp  = closest_on_segment(a, bb, x, &seg_t);
        Vec2 d   = v2_sub(x, cp);
        float w2 = b->paddle_width * b->paddle_width;
        float fall = expf(-v2_len2(d) / w2);
        /* Rigid-body velocity of the paddle at that point. */
        Vec2 rad = v2_sub(cp, b->paddle_pos);
        Vec2 vel = v2_mul(v2_perp(rad), b->paddle_omega);
        U = v2_mad(U, vel, fall);
    }

    /* Free-slip against obstacles: smoothly remove the normal component of the
     * flow as we approach a solid, so the current slides around it instead of
     * dragging wave curves through it. This is a kinematic fudge, not a
     * projection of the whole field -- a real base simulation would supply a
     * velocity that already satisfies the boundary condition. */
    if (b->n_obstacle > 0) {
        Vec2 n;
        float sd = wc_base_sdf(b, x, &n);
        const float skin = 0.16f;
        if (sd < 3.0f * skin) {
            float w = expf(-(wc_maxf(sd, 0.0f) * wc_maxf(sd, 0.0f)) / (skin * skin));
            U = v2_mad(U, n, -w * v2_dot(U, n));
        }
    }

    return U;
}

float wc_base_height(const WcBase *b, Vec2 x, float t)
{
    float h = 0.0f;
    for (int i = 0; i < 2; ++i) {
        if (b->swell_amp[i] == 0.0f) continue;
        h += b->swell_amp[i] * sinf(v2_dot(b->swell_k[i], x) - b->swell_omega[i] * t);
    }
    return h;
}

/* Vertical acceleration of the base surface. For a height field this is the
 * dominant component of the surface acceleration `a` in Eq. (13). */
static float base_accel_z(const WcBase *b, Vec2 x, float t)
{
    float az = 0.0f;
    for (int i = 0; i < 2; ++i) {
        if (b->swell_amp[i] == 0.0f) continue;
        float w = b->swell_omega[i];
        az += -b->swell_amp[i] * w * w * sinf(v2_dot(b->swell_k[i], x) - w * t);
    }
    return az;
}

/* Effective gravity, Eq. (13):   g* = -N . (g - a).
 *
 * With N the unit normal of the height field z = h(x,y,t), g = (0,0,-g0) and
 * a = (0,0,az), this reduces to  g* = Nz * (g0 + az):  the surface tilts away
 * from gravity (Nz < 1) and its vertical acceleration adds to or subtracts from
 * what a wave riding on it feels. */
static float base_geff(const WcBase *b, Vec2 x, float t)
{
    float hx = (wc_base_height(b, v2(x.x + WC_FD_DX, x.y), t) -
                wc_base_height(b, v2(x.x - WC_FD_DX, x.y), t)) / (2.0f * WC_FD_DX);
    float hy = (wc_base_height(b, v2(x.x, x.y + WC_FD_DX), t) -
                wc_base_height(b, v2(x.x, x.y - WC_FD_DX), t)) / (2.0f * WC_FD_DX);
    float nz = 1.0f / sqrtf(1.0f + hx * hx + hy * hy);
    return nz * (b->gravity + base_accel_z(b, x, t));
}

void wc_base_sample(const WcBase *b, Vec2 x, WcBaseSample *out)
{
    const float t = b->time;
    const float h = WC_FD_DX;

    out->U = wc_base_velocity(b, x, t);

    Vec2 uxp = wc_base_velocity(b, v2(x.x + h, x.y), t);
    Vec2 uxm = wc_base_velocity(b, v2(x.x - h, x.y), t);
    Vec2 uyp = wc_base_velocity(b, v2(x.x, x.y + h), t);
    Vec2 uym = wc_base_velocity(b, v2(x.x, x.y - h), t);

    float inv = 1.0f / (2.0f * h);
    out->gradU.a = (uxp.x - uxm.x) * inv;   /* dUx/dx */
    out->gradU.b = (uyp.x - uym.x) * inv;   /* dUx/dy */
    out->gradU.c = (uxp.y - uxm.y) * inv;   /* dUy/dx */
    out->gradU.d = (uyp.y - uym.y) * inv;   /* dUy/dy */

    out->height = wc_base_height(b, x, t);

    float g0  = base_geff(b, x, t);
    float gxp = base_geff(b, v2(x.x + h, x.y), t);
    float gxm = base_geff(b, v2(x.x - h, x.y), t);
    float gyp = base_geff(b, v2(x.x, x.y + h), t);
    float gym = base_geff(b, v2(x.x, x.y - h), t);

    out->g_eff       = g0;
    out->grad_geff.x = (gxp - gxm) * inv;
    out->grad_geff.y = (gyp - gym) * inv;
    out->dgeff_dt    = (base_geff(b, x, t + WC_FD_DT) -
                        base_geff(b, x, t - WC_FD_DT)) / (2.0f * WC_FD_DT);
}
