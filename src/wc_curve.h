/* wc_curve.h -- the wave curve primitive and its solver (paper Section 4).
 *
 * A *wave curve* is a polyline of control points lying on the base surface.
 * Each control point is a Lagrangian wave packet: it carries a position, a
 * wavevector, a phase, a wave action density and a radius. The polyline itself
 * is a wavefront -- a curve of (almost) constant phase -- which is what lets a
 * single primitive represent a long, coherent, connected ripple instead of an
 * isolated blob.
 *
 * Nothing here is hidden behind accessors: this is a reference implementation
 * and the renderer, the HUD and any future analysis code are all expected to
 * walk the data structures directly.
 */
#ifndef WC_CURVE_H
#define WC_CURVE_H

#include "wc_math.h"
#include "wc_wave.h"
#include "wc_base.h"

/* ------------------------------------------------------------------------- */

typedef struct {
    Vec2  x;          /* position on the surface                    [m]        */
    Vec2  k;          /* wavevector                                 [rad/m]    */
    float phase;      /* phi, wrapped to [-pi,pi)                   [rad]      */
    float action;     /* wave action *density*, Eq. (10)            [J s/m^2]  */
    float radius;     /* stripe half-width r                        [m]        */

    /* Derived / bookkeeping, refreshed every step. */
    float area;       /* advected area patch A_i, Fig. 6 & errata E2 [m^2]     */
    float g_eff;      /* effective gravity sampled here             [m/s^2]    */
    float amp;        /* amplitude reconstructed from `action`      [m]        */
    float seed_ramp;  /* seconds of energy ramp-in remaining        [s]        */
    float age;        /* [s] */
    float fade;       /* 1 = normal, decays to 0 when over budget              */
    /* Reference sign of cross(tangent, khat). A point whose actual sign no
     * longer matches has folded and gets killed. It is stored per point rather
     * than per curve because reflecting off an obstacle mirrors the wavefront,
     * which flips the handedness of only the part that has already bounced. */
    float orient;
    bool  alive;
} WcPoint;

typedef struct {
    WcPoint *pt;
    int   n, cap;
    bool  closed;
    int   id;
} WcCurve;

/* ------------------------------------------------------------------------- */

#define WC_MAX_BANDS 6

typedef struct {
    /* Discretisation (Section 4.1 / 4.2) */
    float spacing;          /* target control-point spacing              [m]   */
    float radius0;          /* radius r given to new curves (paper: 0.20) [m]  */
    float radius_min;       /* clamp range for r                         [m]   */
    float radius_max;
    int   substeps;         /* forward-Euler substeps per wc_sim_step()        */

    /* Practical heuristics (Appendix C) */
    float min_steepness;    /* delete control points below this (paper: 0.01)  */
    float max_steepness;    /* clamp amplitude to this (paper: 0.70)           */
    int   max_points;       /* global budget (paper: ~1e6)                     */

    /* Seeding (Section 4.3) */
    bool  seeding_on;
    int   seed_candidates;  /* growth-rate samples evaluated per step          */
    float seed_gain;        /* candidate -> acceptance probability scale       */
    float seed_ramp_time;   /* energy ramps in over this long (paper: 0.5 s)   */
    float seed_spectrum;    /* scale of beta(k) in Eq. (20)-driven growth      */
    float seed_isotropic;   /* small constant added to gamma, Section 4.3      */
    float seed_length;      /* arclength of a newly grown curve          [m]   */
    float seed_align;       /* alpha of Eq. (28); 1 = optimal, 0 = geodesic    */
    int   n_bands;
    float wavelength[WC_MAX_BANDS];   /* seeded wavelengths               [m]  */

    /* Rendering-side steepness limiter, Eq. (29) */
    float critical_steepness;         /* s_c (paper: 3)                        */
} WcParams;

void wc_params_default(WcParams *p);

/* ------------------------------------------------------------------------- */

typedef struct {
    WcWater  water;
    WcParams p;
    WcBase  *base;

    WcCurve *curve;
    int      n_curve, cap_curve;
    int      next_id;

    WcRng    rng;

    /* Statistics for the HUD. */
    int   stat_points;
    int   stat_curves;
    int   stat_seeded;      /* curves grown in the last step */
    int   stat_killed;      /* control points removed in the last step */
    float stat_max_steep;
    float time;
} WcSim;

void wc_sim_init (WcSim *s, WcBase *base, const WcWater *w, const WcParams *p, uint32_t seed);
void wc_sim_free (WcSim *s);
void wc_sim_clear(WcSim *s);
void wc_sim_step (WcSim *s, float dt);

/* Interactive seeding: a closed circular wavefront. Because every wavelength
 * band gets its own ring, dropping one shows dispersion directly -- the long
 * waves outrun the short ones and the ring separates into a train. */
void wc_sim_spawn_ring(WcSim *s, Vec2 center, float radius, float wavelength, float amplitude);
void wc_sim_spawn_ring_spectrum(WcSim *s, Vec2 center, float radius, float amplitude);

/* Grow one curve from a seed point, marching along Eq. (27)/(28). Returns the
 * number of control points created. */
int  wc_sim_grow_curve(WcSim *s, Vec2 seed, Vec2 khat, float wavelength);

#endif /* WC_CURVE_H */
