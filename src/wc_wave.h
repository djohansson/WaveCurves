/* wc_wave.h -- linear water wave theory on a moving surface.
 *
 * This file is the "formula sheet" of the reference implementation: it holds
 * nothing but the closed-form relations of Sections 3.1 and 3.3 of
 *
 *   Skrivan, Soderstrom, Johansson, Sprenger, Museth, Wojtan.
 *   "Wave Curves: Simulating Lagrangian water waves on dynamically deforming
 *    surfaces." ACM TOG 39(4), SIGGRAPH 2020.
 *
 * Equation numbers in the comments refer to that paper. Everything is SI:
 * metres, seconds, kilograms, radians.
 */
#ifndef WC_WAVE_H
#define WC_WAVE_H

#include "wc_math.h"

typedef struct {
    float gravity;       /* g      [m/s^2]   nominal gravity                  */
    float surf_tension;  /* sigma  [N/m]     air/water surface tension        */
    float density;       /* rho    [kg/m^3]  water density                    */
    float viscosity;     /* nu     [m^2/s]   damping coefficient (see below)  */
} WcWater;

static inline WcWater wc_water_default(void)
{
    WcWater w;
    w.gravity      = 9.81f;
    w.surf_tension = 0.0728f;   /* clean water at 20 C */
    w.density      = 1000.0f;
    /* Molecular viscosity (1e-6) damps far too slowly to be visible at these
     * scales, so like Jeschke & Wojtan [2017] we use an artistic value. At 3e-5
     * the amplitude e-folding time is ~1 s for a 5 cm ripple and ~1 min for a
     * 40 cm wave, which gives the spectrum a natural short-wave cutoff. */
    w.viscosity    = 3.0e-5f;
    return w;
}

/* ---------------------------------------------------------------------------
 * Dispersion relation, Eq. (14):
 *
 *      omega(k, g*) = sqrt( ( g* + (sigma/rho) k^2 ) k )
 *
 * where g* is the *effective* gravity of Eq. (13),
 *
 *      g*(x,t) = -N(x,t) . ( g - a(x,t) ),
 *
 * i.e. the component of gravity minus surface acceleration along the surface
 * normal. On a flat surface with no vertical acceleration this collapses to
 * the classic Airy relation omega = sqrt((g + sigma k^2 / rho) k).
 *
 * Note that g* can go negative under strong downward acceleration; the radicand
 * then flips sign and the wave is Rayleigh-Taylor unstable. Linear theory has
 * nothing useful to say there, so we clamp to a small positive floor and let
 * the amplitude limiter of Appendix C carry the load.
 * ------------------------------------------------------------------------- */

#define WC_GEFF_FLOOR 0.05f

static inline float wc_geff_safe(float g_eff)
{
    return g_eff > WC_GEFF_FLOOR ? g_eff : WC_GEFF_FLOOR;
}

static inline float wc_omega(const WcWater *w, float k, float g_eff)
{
    float g = wc_geff_safe(g_eff);
    return sqrtf((g + (w->surf_tension / w->density) * k * k) * k);
}

/* Group speed magnitude, Eq. (4):  c_g = d omega / d k.
 *
 *   omega^2 = g* k + (sigma/rho) k^3
 *   => 2 omega domega/dk = g* + 3 (sigma/rho) k^2
 *   => domega/dk = ( g* + 3 (sigma/rho) k^2 ) / ( 2 omega )
 *
 * The group velocity vector is c_g = (domega/dk) * khat. */
static inline float wc_group_speed(const WcWater *w, float k, float g_eff)
{
    float g  = wc_geff_safe(g_eff);
    float s  = w->surf_tension / w->density;
    float om = sqrtf((g + s * k * k) * k);
    if (om < 1e-9f) return 0.0f;
    return (g + 3.0f * s * k * k) / (2.0f * om);
}

/* Phase speed  c_p = omega / k. */
static inline float wc_phase_speed(const WcWater *w, float k, float g_eff)
{
    if (k < 1e-9f) return 0.0f;
    return wc_omega(w, k, g_eff) / k;
}

/* d omega / d g*, needed for the refraction term of Eq. (16):
 *   2 omega domega/dg* = k   =>   domega/dg* = k / (2 omega). */
static inline float wc_domega_dgeff(const WcWater *w, float k, float g_eff)
{
    float om = wc_omega(w, k, g_eff);
    if (om < 1e-9f) return 0.0f;
    return k / (2.0f * om);
}

/* ---------------------------------------------------------------------------
 * Energy, action and amplitude.
 *
 * Energy density, Eq. (3):     E = 1/2 ( rho g* + sigma k^2 ) A^2
 * Wave action,    Eq. (10):    Action = E / omega
 *
 * Action -- not energy -- is the quantity conserved by Eq. (11) when the
 * environment varies in time, so the solver stores action per control point and
 * converts to amplitude only when it needs to draw something.
 * ------------------------------------------------------------------------- */
static inline float wc_energy_coeff(const WcWater *w, float k, float g_eff)
{
    /* the ( rho g* + sigma k^2 ) factor */
    return w->density * wc_geff_safe(g_eff) + w->surf_tension * k * k;
}

static inline float wc_action_from_amp(const WcWater *w, float amp, float k, float g_eff)
{
    float om = wc_omega(w, k, g_eff);
    if (om < 1e-9f) return 0.0f;
    return 0.5f * wc_energy_coeff(w, k, g_eff) * amp * amp / om;
}

static inline float wc_amp_from_action(const WcWater *w, float action, float k, float g_eff)
{
    float coeff = wc_energy_coeff(w, k, g_eff);
    if (action <= 0.0f || coeff < 1e-9f) return 0.0f;
    return sqrtf(2.0f * action * wc_omega(w, k, g_eff) / coeff);
}

static inline float wc_energy_from_action(const WcWater *w, float action, float k, float g_eff)
{
    return action * wc_omega(w, k, g_eff);
}

static inline float wc_action_from_energy(const WcWater *w, float energy, float k, float g_eff)
{
    float om = wc_omega(w, k, g_eff);
    return (om < 1e-9f) ? 0.0f : energy / om;
}

/* ---------------------------------------------------------------------------
 * Damping.
 *
 * Section 4.2 borrows the Lagrangian damping model of Jeschke & Wojtan [2017]:
 * "a gradual exponential decay of the packet's amplitude with a rate dependent
 * on its wavenumber". For viscous gravity-capillary waves the linear result is
 * an amplitude decay rate of 2 nu k^2 (so energy decays at 4 nu k^2), which is
 * what we use. Short waves die quickly, long waves persist.
 * ------------------------------------------------------------------------- */
static inline float wc_amp_decay(const WcWater *w, float k, float dt)
{
    return expf(-2.0f * w->viscosity * k * k * dt);
}

/* ---------------------------------------------------------------------------
 * Wave energy growth rate, Eq. (20) (with the erratum's effective-gravity term):
 *
 *      gamma(x,k) = -(c_g / c_p) khat . D . khat
 *                   + (1/omega) (domega/dgstar) (dgstar/dt)
 *
 * D = 1/2 (grad U + grad U^T) is the strain rate of the base flow. The first
 * term is positive where the flow *compresses* along the wave direction (-D
 * positive definite), which is the paper's criterion for where small waves grow
 * fastest; the second injects energy where a surface stops accelerating
 * downward (e.g. ballistic water landing in a still pool).
 *
 * Only the seeding step uses gamma; the evolution of existing curves is driven
 * by action conservation instead.
 * ------------------------------------------------------------------------- */
/* Eq. (20) has two structurally different halves:
 *
 *   gamma = -(c_g/c_p) khat.D.khat  +  (1/omega)(domega/dgstar)(dgstar/dt)
 *           \_____ strain _______/     \________ gravity _________/
 *
 * Only the strain term depends on the wave direction. Seeding needs the two
 * separately, because only the strain term is entitled to pick a direction --
 * see the note in seed_curves(). */
static inline void wc_growth_terms(const WcWater *w, float k, float g_eff,
                                   Vec2 khat, Mat2 D, float dgeff_dt,
                                   float *out_strain, float *out_gravity)
{
    float om = wc_omega(w, k, g_eff);
    if (om < 1e-9f) { *out_strain = 0.0f; *out_gravity = 0.0f; return; }
    float cg = wc_group_speed(w, k, g_eff);
    float cp = om / k;
    *out_strain  = -(cg / cp) * m2_quad(D, khat);      /* khat . D . khat */
    *out_gravity = (wc_domega_dgeff(w, k, g_eff) / om) * dgeff_dt;
}

static inline float wc_growth_rate(const WcWater *w, float k, float g_eff,
                                   Vec2 khat, Mat2 D, float dgeff_dt)
{
    float strain, gravity;
    wc_growth_terms(w, k, g_eff, khat, D, dgeff_dt, &strain, &gravity);
    return strain + gravity;
}

/* Seeding spectrum beta(k) of Section 4.3: the (unobservably small) background
 * wave energy that gamma amplifies. The paper uses a form proportional to
 * gstar/k; the constant is a tuning knob, exposed here as `scale`. */
static inline float wc_seed_spectrum(const WcWater *w, float k, float g_eff, float scale)
{
    if (k < 1e-9f) return 0.0f;
    return scale * w->density * wc_geff_safe(g_eff) / k;
}

#endif /* WC_WAVE_H */
