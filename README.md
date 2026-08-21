# Wave Curves — a reference implementation

A self-contained C99 + SDL3 implementation of

> **Wave Curves: Simulating Lagrangian water waves on dynamically deforming surfaces**
> Tomáš Skřivan, Andreas Söderström, John Johansson, Christoph Sprenger, Ken Museth, Chris Wojtan
> *ACM Transactions on Graphics 39(4)* — SIGGRAPH 2020
> <https://visualcomputing.ist.ac.at/publications/2020/WaveCurves/>
> ([paper PDF](https://pub.ista.ac.at/group_wojtan/projects/2020_Skrivan_WaveCurves/wave_curves_2020.pdf) ·
> [errata](https://pub.ista.ac.at/group_wojtan/projects/2020_Skrivan_WaveCurves/wave_curves_errata.pdf))

The goal is legibility, not speed or production features: every equation in the paper that
the method actually needs appears in the code with its number attached, and the numerical
choices are visible rather than buried. SDL3 provides the window, the input and a thread
pool; nothing physical depends on it.

The errata's corrections are both incorporated: the effective-gravity term in the growth
equation (E1) is present, and area is called `area` rather than overloading `a` (E2).

---

## Building

Requires CMake ≥ 3.16 and a C11 compiler. If SDL3 is not installed, the build fetches and
builds it from <https://github.com/libsdl-org/SDL> automatically.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

On Windows with Visual Studio:

```
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --target wavecurves
```

Pass `-DWAVECURVES_FETCH_SDL=OFF` to require a system SDL3 instead of downloading one.

## Running

```sh
./build/wavecurves                          # interactive, paddle scene
./build/wavecurves --scene 1                # river: strong current + vortex street
./build/wavecurves --scene 2 --no-seeding   # still water: pure dispersion sandbox
./build/wavecurves --scene 2 --no-seeding --obstacles 3   # reflect off a wall
```

| key | |
|---|---|
| `space` | pause |
| `1`–`5` | shaded / displacement / total steepness / growth rate γ / base flow |
| `c` | overlay the wave curves themselves, coloured by wavelength band |
| `s` | toggle automatic seeding |
| `o` | cycle obstacle layouts — none · pillars · harbour · wall with a gap |
| `n` | next scene · `r` reset · `f12` screenshot |
| left click | drop a ripple in every wavelength band |
| right click | drop a single-wavelength ripple |

Headless mode simulates without a window and writes PNGs — useful for regressions, since
the whole simulation is deterministic given `--seed`:

```sh
./build/wavecurves --headless --scene 0 --frames 300 --dump-every 100 --out out
```

Full option list: `wavecurves --help`.

---

## Code map

| file | contents |
|---|---|
| `src/wc_math.h` | 2D vectors, 2×2 matrices, symmetric eigen-decomposition, xorshift RNG |
| `src/wc_wave.h` | **the theory**: dispersion, group speed, action↔amplitude, damping, growth rate. Header-only and formula-shaped on purpose |
| `src/wc_base.{h,c}` | the base simulation the waves ride on, behind a one-struct interface |
| `src/wc_curve.{h,c}` | the wave curve primitive: time stepping, action transport, resampling, folding, seeding, budget |
| `src/wc_render.{h,c}` | wave stripes → displacement field → shaded image |
| `src/wc_png.c` | dependency-free PNG writer (stored deflate blocks) |
| `src/main.c` | SDL window, input, band-parallel thread pool, headless driver |

### Where each equation lives

| paper | code |
|---|---|
| (3) energy density, (10) wave action | `wc_action_from_amp` / `wc_amp_from_action` |
| (4) group speed `c_g = ∂ω/∂k` | `wc_group_speed` |
| (13) effective gravity `g* = −N·(g − a)` | `base_geff` in `wc_base.c` |
| (14) dispersion `ω = √((g* + σk²/ρ)k)` | `wc_omega` |
| (15) `ẋ = U + c_g` | `evolve` in `wc_curve.c` |
| (16) `k̇ = −(∂ω/∂g*)∇g* − [∇U]ᵀk` | `evolve` |
| (17) `φ̇ = −ω + c_g·k` | `evolve` |
| (20) growth rate γ | `wc_growth_rate` |
| (21) phase Taylor expansion | `splat_segment` in `wc_render.c` |
| (22) falloff kernel Ψ, (23) amplitude | `lut_window`, `splat_segment` |
| (24) radius stretching `ṙ = r (n·∇U·n)` | `evolve` |
| (25)/errata (E2) action per advected area | `compute_areas` + the area-ratio rescale in `evolve` |
| (27)(28) marching a new curve | `wc_sim_grow_curve` |
| (29) total-steepness limiter | `wc_field_resolve_band` |
| Fig. 6 trapezoidal area patches | `compute_areas` |
| App. C steepness clamp, budget fade-out | `evolve`, `enforce_budget` |

### The one function to replace

Everything the solver knows about the underlying fluid arrives through

```c
void wc_base_sample(const WcBase *b, Vec2 x, WcBaseSample *out);
```

which returns the surface velocity `U`, its gradient `∇U`, the surface height, the effective
gravity `g*` with its spatial gradient and time derivative. Point that at real FLIP or
level-set data and the rest of the method is unchanged — that is the whole premise of the
paper, that this is a post-process.

The bundled scenes build `U` from drifting Lamb–Oseen vortices, converging/diverging surface
cells and a rotating paddle, and difference it numerically for `∇U`. The long background
swell exists specifically so that `∂g*/∂t ≠ 0` and the effective-gravity machinery of
Section 3.3 has something to do.

---

## What is and is not here

**Implemented.** Lagrangian curve advection at `U + c_g`; wavevector refraction by both
`∇U` and `∇g*`; phase transport; radius stretching; wave-action transport over advected
trapezoidal area patches; arclength resampling with magnitude/direction wavevector
interpolation and shortest-arc phase interpolation; fold detection and amplitude zeroing;
minimum- and maximum-steepness heuristics; a global point budget with gradual fade-out;
growth-rate-driven seeding with the 0.5 s energy ramp and curve marching along (27)/(28);
wave-stripe splatting with the Eq. (29) steepness limiter; a band-parallel renderer.

**Deliberately out of scope.**

- **Non-planar surfaces.** The theory of Section 3.3 is written for arbitrary moving
  surfaces, and the solver never assumes otherwise, but the demo's base surface is a height
  field over a plane. Geodesics are therefore trivial and the surface-projection step of
  Section 4.2 is a no-op. Level sets and triangle meshes would need a projection operator
  and a tangent-plane frame; nothing else changes.
- **Domain edges.** Control points that leave the domain are deleted, not reflected. (Solid
  obstacles *inside* the domain do reflect — see below.)
- **Curvature clamps.** Appendix C also clamps amplitude and radius against the base
  surface's maximal curvature, to stop large waves living on droplets. A planar base has no
  droplets, so this is omitted.
- **Scale.** The paper's production scenes carry ~10⁶ control points; the default budget
  here is 24 000, which saturates a 4 m domain about 25 times over.

**Where this deviates from the paper, and why.**

1. **Rendering is a scatter, not a gather.** Section 4.4 ray-casts every vertex of a
   high-resolution surface against the wave stripes. Here the high-resolution surface is a
   regular grid, so each stripe segment is rasterised into the pixels it covers instead.
   The per-(segment, pixel) arithmetic — Eqs. (21) and (23) — is identical; only the
   traversal order differs.
2. **The falloff kernel.** Eq. (22) is printed as `Ψ(x,r) = ½(cos(2πx/r) + 1)` for `|x| ≤ r`,
   which is not monotone over that support — it returns to 1 at the stripe edge. The
   intent is clearly a raised cosine that is 1 on the curve and 0 at the edge, so
   `½(cos(πx/r) + 1)` is used (equivalently, the printed form with a full width of `2r`).
3. **Eq. (20)'s sign convention** was reconstructed from the paper's own prose ("`−D` is
   positive definite where the velocity field squishes the surface together", and the
   fastest-growing direction is the eigenvector of the most negative eigenvalue) plus
   dimensional analysis, because the equation's typesetting is ambiguous. The implemented
   form is `γ = −(c_g/c_p)(k̂·D·k̂) + (1/ω)(∂ω/∂g*)(∂g*/∂t)`.
4. **β(k)'s constant is a free parameter.** The paper gives the shape `β(k) ∝ g*/k` but the
   constant is a tuning knob. It is exposed as `seed_spectrum`, and its default matters more
   than it looks: because a few dozen stripes overlap at any point, the *total* steepness in
   Eq. (29) scales with the curve budget, and if it is pinned against `s_c` the limiter
   saturates and the surface turns into corrugated metal regardless of how correct the
   dynamics are.
5. **α in Eq. (28) carries units of s/m.** The paper's `α = 2` is dimensionally implicit;
   it works out here because the scenes have comparable strain magnitudes.
6. **Seeding uses rejection sampling** on γ over uniformly drawn candidates rather than
   Houdini's density-driven point scatter, and **only seeds where the flow prefers a
   direction** — see the next section.
7. **Forward Euler, one substep** (as in Section 4.2), and our own arclength resampler in
   place of Houdini's polyline tools.

---

## Seeding only where the flow prefers a direction

A wave curve is a long, coherent wavefront, so spawning one commits to a direction. Eq. (20)
has two parts that offer no direction to commit to: the effective-gravity term
`(1/ω)(∂ω/∂g*)(∂g*/∂t)`, and the small isotropic constant of Section 4.3. Seeding from those
produces isolated straight ripples drifting across otherwise calm water at whatever angle the
sampler happened to draw — visually, a background of wandering lines.

That is the failure mode the paper itself names in Section 5:

> If a perfectly isotropic wave spectrum is undersampled by our method, then the spectrum
> will consist of few waves in a few randomly chosen directions, which may also appear
> unnatural.

At the paper's ~10⁶ control points those blend into texture. At the budgets here, each one
reads as a line. So `seed_curves()` lets undirected growth **modulate the rate** of
direction-aligned seeds but never spawn on its own, and `seed_isotropic` defaults to 0.
Raise it (~0.05) together with `max_points` if you want a livelier background.

A direction only counts as preferred when `D`'s eigenvalues are genuinely *separated*.
Testing the compressive eigenvalue alone is not enough — pure radial convergence gives
`D = −c·I`, where every direction is an eigenvector and any eigensolver must fall back to an
arbitrary basis, in practice a coordinate axis. Calm water (`D ≈ 0`) is the same trap. Both
produce axis-aligned wavefronts marching across the domain, which is a convincing-looking
artifact rather than an obviously broken one.

---

## Obstacles and reflection *(an extension — not in the paper)*

`--obstacles N` (or `o` at runtime) drops solid objects into the scene for the waves to
bounce off. Layouts: `1` three pillars, `2` a harbour of angled sea walls and piles,
`3` a wall with two narrow gaps. Objects are signed distance functions (`wc_base_sdf`), so
the solver only ever asks *how deep am I and which way is out* and the reflection code is
shape-agnostic.

Reflection follows ray theory, which is what the rest of the method already assumes. Each
control point is handled independently: if it ends a step inside a solid it is pushed back
onto the surface and its wavevector is mirrored about the outward normal,
`k ← k − 2(k·n̂)n̂`. Phase, action and radius are untouched, so the wall is a perfect
reflector. A wavefront only partly past an obstacle therefore develops a travelling kink at
the contact point, and the reflected arm reassembles into the correct mirror image over the
following few steps — the same way a bundle of rays does.

One thing this required changing in the solver: the fold test of Section 4.2 compares the
sign of `cross(tangent, k̂)` against a reference recorded when the curve was born. Mirroring
`k̂` flips that handedness, so the reference sign had to move from the curve to the
individual control point and flip along with the reflection — otherwise the entire reflected
arm reads as folded and gets deleted the moment it bounces.

Known limits, all inherited from ray theory:

- **No diffraction.** Waves leave a hard-edged geometric shadow behind an obstacle instead
  of bending around it, and the gaps in layout 3 do not spread the wave into a circular
  front the way a real slit would. Section 6 of the paper notes the same limitation.
- **Total reflection.** There is no reflection coefficient and no transmission; a wall
  returns all the energy that reaches it.
- **Corners focus.** The normal turns discontinuously at a sharp corner, so control points
  reflecting near one converge and pile up. Eq. (29) keeps the result bounded but there is a
  visible bright spot at wall ends.
- **The base flow is faked around solids.** `wc_base_velocity` smoothly removes the normal
  component of the current near an obstacle so the flow slides along it rather than dragging
  wave curves through it. That is a kinematic patch, not a projection — a real base
  simulation would already satisfy the boundary condition, and this is exactly the sort of
  thing that disappears when you point `wc_base_sample` at real data.

---

## Checking that it works

Run the still-water scene with automatic seeding off and drop one click. The four
wavelength bands are launched as coincident rings and immediately separate, because each
travels at its own group speed `c_g = ½√(g/k)` for deep-water gravity waves:

```sh
./build/wavecurves --headless --scene 2 --no-seeding --frames 150 --dump-every 149 --out disp
```

At `t = 2.5 s` the 40 cm ring has reached `r ≈ 1.2 m` and the 5 cm ring `r ≈ 0.5 m`,
against predictions of `0.40 m/s` and `0.14 m/s` — the ring radii read straight off the
image match. The 5 cm ring is also visibly the faintest, which is the `exp(−2νk²t)` damping
of Section 4.2 doing its job.

For reflection, the flat wall of layout 3 is the cleanest check — the reflected wavefront
should be a circle centred on the *mirror image* of the source:

```sh
./build/wavecurves --headless --scene 2 --no-seeding --obstacles 3 \
                   --frames 200 --dump-every 199 --out slit
```

The other two scenes exercise the parts of the theory that need a moving surface: in
`--scene 0` the paddle's strain field seeds new curves along its compressive axis and
stretches everything else into flow-aligned streaks (the paper's Figs. 4 and 12); in
`--scene 1` a vortex street winds wavefronts into spirals while the mean current carries
them downstream. Press `4` in either to see the raw γ field that drives the seeding.

---

## Parameters worth knowing

All defaults live in `wc_params_default()` (`wc_curve.c`) and `wc_water_default()`
(`wc_wave.h`).

| name | default | meaning |
|---|---|---|
| `spacing` | 0.045 m | target control-point spacing after resampling |
| `radius0` | 0.20 m | stripe half-width `r` given to new curves (the paper's value) |
| `min_steepness` | 0.01 | below this a control point is deleted (Section 4.2) |
| `max_steepness` | 0.70 | per-wave amplitude clamp (Appendix C) |
| `critical_steepness` | 3.0 | `s_c` in the Eq. (29) limiter |
| `max_points` | 24 000 | global control-point budget |
| `seed_ramp_time` | 0.5 s | energy ramp for a newly seeded curve (Section 4.3) |
| `seed_spectrum` | 6e-6 | constant in `β(k) = seed_spectrum · ρ g*/k` |
| `seed_isotropic` | 0 (off) | undirected background seeding; see the section above |
| `seed_align` | 2.0 | `α` in Eq. (28); 1 → optimal directions, 0 → geodesic curves |
| `viscosity` | 3e-5 m²/s | amplitude damping `exp(−2νk²Δt)`; artistic, not molecular |
| wavelength bands | 5, 10, 20, 40 cm | the paper seeds 5, 10, 20, 25 and 125 cm |
