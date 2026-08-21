/* wc_base.h -- the "base simulation" that wave curves ride on.
 *
 * Wave Curves is a post-process: it consumes an existing water surface and adds
 * detail on top. Everything the solver needs from that surface is funnelled
 * through WcBaseSample, so swapping the procedural scene below for real FLIP /
 * level-set / mesh data means implementing exactly one function.
 *
 * The reference scenes keep the surface a height field over a planar domain.
 * That is not a limitation of the method -- the paper runs it on level sets and
 * triangle meshes -- but it keeps geodesics trivial and lets the demo stay
 * focused on the wave physics rather than on surface parameterisation.
 */
#ifndef WC_BASE_H
#define WC_BASE_H

#include "wc_math.h"

typedef struct {
    Vec2  U;           /* surface (tangential) velocity                  [m/s]   */
    Mat2  gradU;       /* (grad U)_ij = dU_i/dx_j                        [1/s]   */
    float height;      /* vertical displacement of the base surface      [m]     */
    float g_eff;       /* effective gravity, Eq. (13)                    [m/s^2] */
    Vec2  grad_geff;   /* spatial gradient of g_eff, for Eq. (16)        [1/s^2] */
    float dgeff_dt;    /* d g_eff / dt, for the seeding rate Eq. (20)    [m/s^3] */
} WcBaseSample;

typedef enum {
    WC_SCENE_PADDLE = 0,   /* stirring paddle + vortices: the paper's Fig. 4     */
    WC_SCENE_RIVER  = 1,   /* strong current, shear layers and obstacles         */
    WC_SCENE_STILL  = 2,   /* no flow, flat surface: pure dispersion sandbox     */
    WC_SCENE_COUNT  = 3
} WcSceneId;

#define WC_MAX_VORTICES 6
#define WC_MAX_CELLS    6
#define WC_MAX_OBSTACLES 8
#define WC_OBSTACLE_LAYOUTS 4

/* Static solid objects the waves reflect from. Represented as signed distance
 * functions so the solver only ever needs "how deep am I, and which way is out",
 * which makes the reflection step shape-agnostic. */
typedef enum { WC_OBS_CIRCLE = 0, WC_OBS_BOX = 1 } WcObstacleKind;

typedef struct {
    WcObstacleKind kind;
    Vec2  pos;
    float radius;    /* circle */
    Vec2  half;      /* box half-extents in its own frame */
    float angle;     /* box rotation [rad] */
} WcObstacle;

typedef struct {
    Vec2  pos;
    float strength;   /* circulation-like [m^2/s]; sign sets rotation direction */
    float core;       /* core radius [m] */
    Vec2  drift;      /* [m/s] */
} WcVortex;

/* A radially converging or diverging patch. Real 3D incompressible flow still
 * has non-zero *surface* divergence wherever fluid rises and spreads (or sinks),
 * and that divergence is what feeds the seeding criterion of Eq. (20). */
typedef struct {
    Vec2  pos;
    float strength;   /* <0 converging (waves grow), >0 diverging */
    float radius;     /* [m] */
    float pulse_hz;
} WcCell;

typedef struct {
    WcSceneId scene;
    float world_w, world_h;    /* domain extent [m] */
    float time;                /* [s] */

    Vec2  current;             /* uniform background current [m/s] */

    WcVortex vortex[WC_MAX_VORTICES];
    int      n_vortex;
    WcCell   cell[WC_MAX_CELLS];
    int      n_cell;

    WcObstacle obstacle[WC_MAX_OBSTACLES];
    int        n_obstacle;
    int        obstacle_layout;

    /* Rotating paddle: a bar of half-length `paddle_len` spinning about
     * `paddle_pos`, dragging the surface along with it. */
    bool  paddle_on;
    Vec2  paddle_pos;
    float paddle_len;
    float paddle_width;
    float paddle_omega;        /* [rad/s] */
    float paddle_angle;        /* [rad]   */

    /* Long background swell, sum of two travelling sinusoids. Its vertical
     * acceleration is what makes the effective gravity of Eq. (13) interesting:
     * without it g_eff would be a constant and half the paper would be moot. */
    float swell_amp[2];        /* [m]     */
    Vec2  swell_k[2];          /* [rad/m] */
    float swell_omega[2];      /* [rad/s] */

    float gravity;             /* nominal g, mirrored from WcWater */
} WcBase;

void  wc_base_init(WcBase *b, WcSceneId scene, float world_w, float world_h, float gravity);
void  wc_base_advance(WcBase *b, float dt);

/* Primitive queries. */
Vec2  wc_base_velocity(const WcBase *b, Vec2 x, float t);
float wc_base_height(const WcBase *b, Vec2 x, float t);

/* Everything the wave-curve solver needs, in one call. */
void  wc_base_sample(const WcBase *b, Vec2 x, WcBaseSample *out);

/* Install one of the built-in obstacle layouts (0 = none). */
void  wc_base_set_obstacles(WcBase *b, int layout);

/* Signed distance to the nearest obstacle: negative inside, and `out_normal`
 * (may be NULL) receives the outward unit normal there. Returns a large
 * positive number when the scene has no obstacles. */
float wc_base_sdf(const WcBase *b, Vec2 x, Vec2 *out_normal);

const char *wc_scene_name(WcSceneId s);
const char *wc_obstacle_layout_name(int layout);

#endif /* WC_BASE_H */
