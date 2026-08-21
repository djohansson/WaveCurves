/* wc_render.h -- turning wave curves into a displaced, shaded water surface
 * (paper Section 4.4).
 *
 * The paper thickens every wave curve into a *wave stripe* of half-width r and
 * then ray-casts each vertex of a high-resolution surface against those stripes,
 * evaluating Eqs. (21) and (23) at every hit. Here the high-resolution surface
 * is a regular grid, so we invert the loop and scatter each stripe segment into
 * the pixels it covers. The arithmetic per (segment, pixel) pair is identical;
 * scattering just skips the acceleration structure.
 *
 * Every band function takes an explicit row range [y0,y1) and only writes rows
 * inside it, so a caller can hand disjoint bands to different threads without
 * any synchronisation. That is the "trivially parallelizable" claim of the
 * paper, made literal.
 */
#ifndef WC_RENDER_H
#define WC_RENDER_H

#include "wc_curve.h"
#include <stdint.h>

typedef enum {
    WC_VIEW_SHADED = 0,   /* lit water surface                                 */
    WC_VIEW_HEIGHT,       /* signed displacement, blue/red                     */
    WC_VIEW_STEEPNESS,    /* total steepness s_t before the Eq. (29) limiter   */
    WC_VIEW_GROWTH,       /* seeding criterion gamma of Eq. (20)               */
    WC_VIEW_FLOW,         /* base flow speed                                   */
    WC_VIEW_COUNT
} WcViewMode;

const char *wc_view_name(WcViewMode m);

typedef struct {
    int   w, h;
    float world_w, world_h;
    float dx, dy;              /* world units per pixel */

    float *disp;               /* sum_i A_i sin(phi_i)                  [m]   */
    float *steep;              /* sum_i s_i = sum_i A_i k_i             [-]   */
    float *height;             /* base height + limited displacement    [m]   */
    uint32_t *pixels;          /* RGBA8888, little-endian ABGR word          */

    /* Obstacle coverage, rebuilt only when the layout changes: 0 outside, and
     * inside a solid it encodes depth so the edge can be rimmed. Obstacles are
     * static, so this costs nothing per frame. */
    uint8_t *obstacle;

    /* Coarse debug fields, refreshed on demand (they need base-simulation
     * samples, which are far too expensive to evaluate per pixel). */
    int    dbg_n;
    float *dbg_growth;
    float *dbg_speed;
} WcField;

void wc_field_init(WcField *f, int w, int h, float world_w, float world_h);
void wc_field_free(WcField *f);

/* Pipeline, in order. Each may be called concurrently on disjoint bands. */
void wc_field_clear_band  (WcField *f, int y0, int y1);
void wc_field_splat_band  (WcField *f, const WcSim *s, int y0, int y1);
void wc_field_resolve_band(WcField *f, const WcSim *s, int y0, int y1);
void wc_field_shade_band  (WcField *f, const WcSim *s, WcViewMode mode, int y0, int y1);

/* Single-threaded convenience wrapper around the whole pipeline. */
void wc_field_render(WcField *f, const WcSim *s, WcViewMode mode);

/* Recompute the coarse debug fields. Call before shading a debug view. */
void wc_field_update_debug(WcField *f, const WcSim *s);

/* Rasterise the obstacle mask. Call once after changing the obstacle layout. */
void wc_field_update_obstacles(WcField *f, const WcBase *b);

/* Peak |displacement| over the field, for the HUD. */
float wc_field_peak(const WcField *f);

#endif /* WC_RENDER_H */
