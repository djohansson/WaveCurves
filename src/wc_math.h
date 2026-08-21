/* wc_math.h -- small 2D linear algebra + deterministic RNG.
 *
 * Part of the Wave Curves reference implementation.
 * Header-only, no dependencies beyond <math.h>.
 */
#ifndef WC_MATH_H
#define WC_MATH_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>

#define WC_PI  3.14159265358979323846f
#define WC_TAU 6.28318530717958647692f

typedef struct { float x, y; } Vec2;

/* Row-major 2x2 matrix:   | a  b |
 *                         | c  d |
 * For a velocity gradient we use the convention (grad U)_ij = dU_i/dx_j, i.e.
 *   a = dUx/dx,  b = dUx/dy,
 *   c = dUy/dx,  d = dUy/dy. */
typedef struct { float a, b, c, d; } Mat2;

static inline Vec2 v2(float x, float y)            { Vec2 r = { x, y }; return r; }
static inline Vec2 v2_add(Vec2 p, Vec2 q)          { return v2(p.x + q.x, p.y + q.y); }
static inline Vec2 v2_sub(Vec2 p, Vec2 q)          { return v2(p.x - q.x, p.y - q.y); }
static inline Vec2 v2_mul(Vec2 p, float s)         { return v2(p.x * s, p.y * s); }
static inline Vec2 v2_mad(Vec2 p, Vec2 q, float s) { return v2(p.x + q.x * s, p.y + q.y * s); }
static inline float v2_dot(Vec2 p, Vec2 q)         { return p.x * q.x + p.y * q.y; }
static inline float v2_cross(Vec2 p, Vec2 q)       { return p.x * q.y - p.y * q.x; }
static inline float v2_len2(Vec2 p)                { return p.x * p.x + p.y * p.y; }
static inline float v2_len(Vec2 p)                 { return sqrtf(v2_len2(p)); }
static inline Vec2 v2_perp(Vec2 p)                 { return v2(-p.y, p.x); }
static inline Vec2 v2_lerp(Vec2 p, Vec2 q, float t){ return v2(p.x + (q.x - p.x) * t,
                                                               p.y + (q.y - p.y) * t); }
static inline Vec2 v2_norm(Vec2 p)
{
    float l = v2_len(p);
    return (l > 1e-20f) ? v2(p.x / l, p.y / l) : v2(1.0f, 0.0f);
}

static inline float wc_clampf(float v, float lo, float hi)
{
    return v < lo ? lo : (v > hi ? hi : v);
}
static inline float wc_lerpf(float a, float b, float t) { return a + (b - a) * t; }
static inline float wc_minf(float a, float b) { return a < b ? a : b; }
static inline float wc_maxf(float a, float b) { return a > b ? a : b; }

/* Wrap an angle into [-pi, pi). Phases are stored wrapped so that they keep
 * full float precision no matter how long the simulation runs. */
static inline float wc_wrap_pi(float a)
{
    a = fmodf(a + WC_PI, WC_TAU);
    if (a < 0.0f) a += WC_TAU;
    return a - WC_PI;
}

/* Shortest-arc interpolation between two wrapped phases. */
static inline float wc_lerp_phase(float a, float b, float t)
{
    return a + wc_wrap_pi(b - a) * t;
}

static inline Vec2 m2_apply(Mat2 m, Vec2 v)
{
    return v2(m.a * v.x + m.b * v.y, m.c * v.x + m.d * v.y);
}

/* m^T * v -- needed for the [grad U]^T k term of Eq. (16). */
static inline Vec2 m2_apply_t(Mat2 m, Vec2 v)
{
    return v2(m.a * v.x + m.c * v.y, m.b * v.x + m.d * v.y);
}

/* D = 1/2 (grad U + grad U^T): the strain-rate tensor of Eq. (19). */
static inline Mat2 m2_sym(Mat2 m)
{
    float off = 0.5f * (m.b + m.c);
    Mat2 r = { m.a, off, off, m.d };
    return r;
}

static inline float m2_trace(Mat2 m) { return m.a + m.d; }

/* v^T M v, with M assumed symmetric. */
static inline float m2_quad(Mat2 m, Vec2 v)
{
    return v.x * (m.a * v.x + m.b * v.y) + v.y * (m.c * v.x + m.d * v.y);
}

/* Eigen-decomposition of a symmetric 2x2 matrix.
 * Returns eigenvalues sorted lo <= hi with their unit eigenvectors.
 * Used by the seeding step to find the most compressive flow direction.
 *
 * CAUTION: when hi == lo the eigenvectors are degenerate -- every direction is
 * an eigenvector -- and this returns an arbitrary basis, which for a zero
 * matrix is the coordinate axes. Callers must not treat that as a meaningful
 * direction; check that the eigenvalue separation is real first. */
typedef struct { float lo, hi; Vec2 v_lo, v_hi; } Eig2;

static inline Eig2 m2_eigen_sym(Mat2 m)
{
    Eig2 e;
    float mean = 0.5f * (m.a + m.d);
    float diff = 0.5f * (m.a - m.d);
    float off  = 0.5f * (m.b + m.c);           /* symmetrise defensively */
    float disc = sqrtf(diff * diff + off * off);

    e.hi = mean + disc;
    e.lo = mean - disc;

    if (fabsf(off) > 1e-20f) {
        e.v_hi = v2_norm(v2(off, e.hi - m.a));
    } else {
        e.v_hi = (m.a >= m.d) ? v2(1.0f, 0.0f) : v2(0.0f, 1.0f);
    }
    e.v_lo = v2_perp(e.v_hi);                  /* symmetric => eigenvectors orthogonal */
    return e;
}

/* ---------------------------------------------------------------------------
 * Deterministic RNG (xorshift32). The whole simulation is reproducible from a
 * seed, which matters for the headless regression mode.
 * ------------------------------------------------------------------------- */
typedef struct { uint32_t s; } WcRng;

static inline void wc_rng_seed(WcRng *r, uint32_t seed)
{
    r->s = seed ? seed : 0x9e3779b9u;
}

static inline uint32_t wc_rng_u32(WcRng *r)
{
    uint32_t x = r->s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return (r->s = x);
}

/* Uniform in [0,1). */
static inline float wc_rng_f(WcRng *r)
{
    return (float)(wc_rng_u32(r) >> 8) * (1.0f / 16777216.0f);
}

/* Uniform in [lo,hi). */
static inline float wc_rng_range(WcRng *r, float lo, float hi)
{
    return lo + (hi - lo) * wc_rng_f(r);
}

#endif /* WC_MATH_H */
