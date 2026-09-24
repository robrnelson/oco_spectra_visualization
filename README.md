# OCO-2 Spectra Explorer — extended

An interactive, single-file visualisation of what the three OCO-2 bands see.

Built on the original
[**OCO-2 Spectra Explorer**](https://rsfralab.github.io/oco/widgets/spectra-explorer.html)
created by **Christian Frankenberg** (Caltech/JPL) — this is an extended rewrite of
that widget, and it reads the same pre-computed line-by-line `data/` JSON files
produced for it. The original is mirrored here as `orig.html` for comparison.

```
python3 validate.py                 # validate data, write the selftest reference
python3 -m http.server 8000         # then open http://localhost:8000/
```

Open `http://localhost:8000/?selftest=1` to have the page prove its forward model
matches the NumPy reference. It must be served over HTTP — opening `index.html`
directly from the filesystem will fail on `fetch` due to CORS.

## Forward model

Line-by-line, then convolved with the real OCO-2 ILS:

```
τ_gas = psurf/P₀ · (XCO₂/420 · τ_co2 + τ_o2) + h2o · τ_h2o
τ_ray = 0.008569 µm⁻⁴ (1 + 0.0113 µm⁻² + 0.00013 µm⁻⁴) · psurf/P₀     [Bodhaine 1999]
τ_aer = AOD₅₅₀ · (λ/550)^(−å)
m     = 1/cos(SZA) + 1/cos(VZA)
a(λ)  = albScale · ECOSTRESS(λ) + albSlope · (ν − ν_centre),   ν = 10⁷/λ cm⁻¹

L  = (E₀ cosSZA / π) [ R·exp(−m·f·τ_gas) + T·a·T_diff·exp(−m·τ_gas) / (1 − a·s) ]
L_c = (E₀ cosSZA / π) [ R              + T·a·T_diff              / (1 − a·s) ]
```

`R`, `T` (beam) and `s`, `T_diff` (diffuse) are the reflectance and transmittance of
a single scattering layer holding the mixed Rayleigh + aerosol optical depth — see
**Scattering** below. With no aerosol and Rayleigh off, `R = 0`, `T = T_diff = 1`,
`s = 0`, and this collapses *bit-identically* to the plain
`L = E₀ cosSZA · a/π · exp(−m τ_gas)`.

`τ_h2o` is deliberately **not** scaled by surface pressure — it is a tropospheric
profile, and the data documentation says to scale only `τ_co2` and `τ_o2`. The H₂O
slider is a separate column scale factor.

### Scattering

Aerosol and Rayleigh are mixed into **one two-stream layer coupled to the Lambertian
surface**, not treated as pure extinction. That distinction matters enormously:
photons taken out of the direct beam are not lost, they still reach the ground
diffusely and still come back up. At the bundled sounding's AOD of 0.071, pure
extinction dims the A-band **10.2 %**; the two-stream treatment dims it **0.35 %**.

- Optical properties mix by the standard τ-weighted rule: τ adds, and ω₀ and g are
  τ-weighted between Rayleigh (ω₀ = 1, g = 0) and aerosol.
- The layer reflectance/transmittance use **delta-Eddington scaling** with the
  Meador & Weaver (1980) Eddington-closure beam solution, plus the classical
  two-stream forms for the diffuse upward path off the surface.
- `f` is the fraction of the well-mixed gas column **above** the layer. The
  atmospheric-scatter term crosses only that much gas, which is precisely the
  mechanism by which aerosol biases retrieved XCO₂. Raise the layer and watch the
  lines fill in. Rayleigh scattering is distributed ∝ pressure, so its mean
  scattering level sits at psurf/2, i.e. `f = 0.5`.
- **Accuracy is measured, not asserted.** `validate_scattering.py` compares the
  system albedo against an independent Monte Carlo solution (Henyey–Greenstein
  phase function, Lambertian surface, collimated beam) over the full slider range:
  **worst absolute error 0.010**, and 0.6 % for the actual sounding case.

Turn the **two-stream scattering** toggle off to get extinction-only back and see
the difference directly.

Still approximate, and worth knowing: one homogeneous layer rather than a profile;
the atmospheric path radiance is treated as isotropic, so there is no
phase-function angular dependence; and the diffuse surface path is charged the
direct-beam gas slant path. It is a two-stream approximation, not a benchmark RT
code.

### Chlorophyll fluorescence (SIF)

SIF is an **additive surface source**, not a reflectance. It emits in the far red only, so
it is present in the A-band and **identically zero** in the two SWIR CO₂ bands. It leaves
the surface, so it crosses the gas **once** rather than twice:

```
L   += SIF · T_diff/(1 − a·s) · exp(−τ_gas / cos VZA)
L_c += SIF · T_diff/(1 − a·s)
```

In a line with τ_gas = 3 that is 32× less attenuation than reflected sunlight at SZA 30°,
and 403× at SZA 60°. So SIF **fills in deep lines** far more than it lifts the continuum:
measured here, SIF = 2 raises the continuum 0.92 % but the deep line cores 2.60 % —
**2.8× stronger in the cores**. That asymmetry is the whole basis of SIF retrieval. Because
the solar spectrum here carries the Toon Fraunhofer lines while SIF is flat, in-filling of
**solar Fraunhofer lines** also appears — the signal OCO-2 actually uses at 757 and 771 nm.
Range 0–4 W m⁻² µm⁻¹ sr⁻¹ (= mW m⁻² nm⁻¹ sr⁻¹); typical vegetated scenes are 0.5–2.5.
Treated as spectrally flat, whereas real retrievals fit a low-order polynomial.

### Pressure broadening — what you can and cannot see

The *consequences* of pressure broadening are in the data and visible; its *pressure
dependence* is not.

- **Visible.** The lines are genuinely pressure-broadened. An isolated A-band O₂ line at
  769.13 nm has FWHM 0.0040 nm = 0.068 cm⁻¹ (HWHM 0.034 cm⁻¹) against a Doppler HWHM of
  0.013 cm⁻¹ at 250 K — ~2.6× wider than Doppler-limited. But the OCO-2 ILS FWHM is 10.4×
  the line FWHM, so you must narrow the **ILS width** to ~**0.1–0.25×** before the shape
  emerges at all.
- **Not visible.** τ is a *single fixed vertical profile* scaled by psurf/P₀. The surface
  pressure slider changes absorber **amount**, never line **width** — the τ FWHM stays at
  0.00400 nm at 1013, 500 and 200 hPa. Real physics would narrow the lines too.
- **A trap.** The *apparent* width in transmittance does shrink with pressure (0.0100 nm at
  1013 hPa → 0.0070 nm at 200 hPa), but only because the line is less saturated. That is
  not broadening.

So the tool captures the *amount* half of why the A-band constrains surface pressure; the
line-shape half needs per-layer τ, which this data set does not provide.

### The AOD–albedo degeneracy

With the aerosol at the surface (`f = 1`) both radiance terms share `exp(−m·τ_gas)`, so the
model factorises into (system albedo) × `exp(−m·τ_gas)`. AOD and surface albedo then enter
through a single scalar and **neither can change the spectrum's shape** — `L/L_c` becomes
independent of AOD, albedo, ω₀ and g (verified to machine precision for albedo). Shape
residual after best-fit scaling: 0.003–0.014 % of continuum. Geometry does not break it.

What does: **lifting the layer** (0.008 % at the surface → 0.11 % at 1 km → 0.35 % at 3 km,
and geometry then helps — 0.35 % → 1.19 % going to SZA 65°), **Rayleigh** (pinned at f = 0.5,
worth ~0.25 % by itself), and **cross-band leverage** (τ_aer A-band : strong CO₂ = 2.7 : 1 at
å = 1, so all three bands together separate them if the surface spectral shape is constrained).

### Other approximations you should know about

- **Spectral albedo is sampled at 5 nm** (the ECOSTRESS grid): 5 points across the
  A-band, 10 across the strong-CO₂ band. Enough for a continuum tilt, not for fine
  spectral structure.
- **The albedo slope control is global**, while L2 reports a slope per band. It is a
  manual exploration knob; `load L2 state` sets it to 0 rather than guessing a
  convention (see "About load L2 state" below).
- **Noise is a display effect**, σ = continuum/SNR, hashed deterministically on
  sample index so the curve does not shimmer while you drag other sliders.

## What is different from the original

**A real bug fixed.** The original's `convolveSame` assumes the ILS kernel is
centred on its array midpoint. It is not: `delta_lambda_nm` shows the A-band centre
sits 0.0090 nm away — 22 % of that band's 0.0415 nm FWHM. Placing the kernel using
`delta_lambda_nm` cuts the A-band residual against the real sounding:

| band | original centring | corrected | of continuum |
|---|---|---|---|
| O₂ A-band | 6.599 | **3.883** | 5.13 % → 3.02 % |
| weak CO₂ | 1.152 | 1.134 | 3.84 % → 3.78 % |
| strong CO₂ | 0.330 | 0.330 | already centred |

(RMS of sounding − model, W m⁻² µm⁻¹ sr⁻¹, at the L2 retrieval state.)
`validate.py` asserts all six numbers, so the fix cannot silently regress.

**Faster.** The original convolves both `L` and `L_c` across the whole fine grid,
then plots every 4th point. Evaluating the convolution only at the output positions
is bit-identical (verified: max abs diff 0.000e+00) and much cheaper. Measured in
gjs/mozjs60 at 2000 output samples per band: **98.6 ms → 13.9 ms per redraw
(7.1×)**. A continuum cache removes roughly half of what remains while dragging
XCO₂, surface pressure, H₂O or aerosol, since none of those change `L_c`.

**More physics.** Independent viewing zenith angle, H₂O column scale, a two-stream
aerosol/Rayleigh scattering layer with adjustable AOD, Ångström exponent, single-scatter
albedo, asymmetry parameter and layer height, spectrally-resolved ECOSTRESS albedo with
scale and slope, ILS width and spectral-shift controls, and a noise/SNR model.

**Data the original discarded.** `albedo_surfaces.json` carries full 400–2200 nm
reflectance spectra but only its scalar band means were used — snow varies 34 %
across the strong-CO₂ band and conifer 8 % across the weak-CO₂ band. The `spectral
albedo` toggle shows the difference. From `sounding_oco2.json`, the per-band
`albedo_l2`, `vza` and `xco2_uncertainty_ppm` are now used, and `aod_total_l2` is
reported in the footnote — though deliberately not loaded into the model, for the
reason given under "About load L2 state" below.

**More to explore.** Four view modes (radiance, transmittance, optical depth,
per-absorber breakdown), single-band focus, brush-to-zoom and pan per band,
snapshot ghost curves for A/B comparison, permalinks, and CSV export alongside the
2× PNG export.

## Spectral sorting (the "sorted" view)

Follows the method of **Zeng et al. (2018)**. Channels are ordered by ascending
**clear-sky** radiance — the same scene with AOD = 0 — and that single ordering is
then applied to every curve drawn, so scenarios stay directly comparable. Saturated
line cores land on the left, intermediate lines in the middle, continuum on the right.

```
perm   = argsort(L_clear)      (ascending; AOD = 0, everything else unchanged)
x-axis = sorted channel index
y-axis = radiance, reordered by perm
```

The physical basis, in the paper's words: *"atmospheric aerosols scatter photons back
to space and therefore reduce the chance of the photons being absorbed by the oxygen…
Photons scattered by higher aerosol layers undergo shorter absorption paths, thereby
reducing the O₂ absorption depths."* That is exactly the `f` factor in the scattering
model above, so the tool reproduces the method's discriminants. Measured on the
A-band (desert surface, SZA 45°):

| change | continuum (right) | intermediate (middle) | saturated (left) |
|---|---|---|---|
| AOD 0 → 1.0 at 800 hPa | **+4.5 %** | +7.6 % | +21 % |
| layer 1000 → 150 hPa at AOD 0.5 | +0.8 % | **+22.6 %** | +240 % |
| albedo × 0.6 → × 1.6 | ×2.3 | ×2.3 | ×2.3 |

So the **continuum tracks total AOD**, the **intermediate lines track aerosol layer
height** (28× more strongly than the continuum does), and **surface albedo scales the
whole curve without changing its shape** — the intermediate/continuum ratio holds at
0.3895 / 0.3883 / 0.3876 across that albedo range. That last row is the AOD–albedo
degeneracy the paper works around, and it is why ALH is separable while AOD is not
without knowing the surface.

The three regions are shaded and labelled. The cuts used here are clear-sky
transmittance < 0.05, 0.05–0.80 and ≥ 0.80; those are a **display convention for this
tool**, not the paper's exact retrieval windows, which are in its supporting
information. Note also that the paper works in the 1.27 µm O₂ band with a mountaintop
FTS and reports reflectance, whereas this applies the same construction to the OCO-2
bands in radiance. One effect it cannot reproduce: the paper's small sorted-curve
wiggles come from O₂ absorption coefficients varying with pressure and temperature
across altitudes, and this tool has a single fixed vertical τ profile.

Drag on the sorted axis to zoom into the intermediate region (the paper's Figure 2d
does the same), and use **snapshot** to freeze one scenario and compare.

## Controls

| Group | Controls |
|---|---|
| Geometry | solar zenith angle, viewing zenith angle |
| Atmosphere | XCO₂ (0–800 ppm), surface pressure (0–1030 hPa), H₂O column scale (0–10×) |
| Aerosol & scattering | AOD₅₅₀ (0–1.5), Ångström exponent, single-scatter albedo ω₀ (0.5–1), asymmetry g (0–0.9), layer pressure (50–1030 hPa), two-stream scattering toggle, Rayleigh |
| Surface | 5 ECOSTRESS surfaces + L2, spectral-albedo toggle, albedo scale, albedo slope, SIF (0–4 W m⁻² µm⁻¹ sr⁻¹) |
| Instrument | ILS width (0.1–3× FWHM), spectral shift (±0.05 nm), SNR, noise |
| Sounding | overlay the real sounding, load its L2 state |

In **radiance** and **transmittance** the model line takes its band's own colour — O₂ A-band
blue, weak CO₂ green, strong CO₂ red. The other views keep one model colour, since
**absorbers** already colours by species and **optical depth** and **sorted** are for
comparing bands against each other.

Views: **radiance** · **transmittance** (L/L_c, comparable across bands) ·
**optical depth** (−ln L/L_c, log axis) · **absorbers** (apparent τ per species) ·
**sorted** (spectral sorting after Zeng et al. 2018, see below).

A **docs** button in the header opens an in-app page with the full model description,
all equations, the approximations, and the citation list — the same material as this
README, available without leaving the tool.

Plot interaction: drag to zoom a band, shift-drag to pan, double-click to reset.

## Things to try

- Move **XCO₂**: lines deepen in the two CO₂ panels, the A-band holds still.
- Raise **AOD** to 0.5, then toggle **two-stream scattering** off and on. Off, the
  A-band collapses; on, it barely moves. That gap is the whole reason extinction-only
  is not good enough.
- With AOD up, drag the **aerosol layer** from 1000 hPa to 100 hPa: absorption lines
  fill in as scattered photons short-circuit more of the gas column. This is the
  aerosol-induced XCO₂ bias, live.
- Pull **ω₀** down to 0.5 (soot rather than sulphate): now the aerosol really does
  darken the scene, because the photons are absorbed rather than redirected.
- Pick **snow**, switch to the strong-CO₂ band, toggle **spectral albedo**: the
  continuum tilts, because snow's reflectance falls 34 % across that band.
- Switch to **absorbers**: O₂ owns the A-band; H₂O contaminates the strong-CO₂ band.
- Raise **SIF** to 2 over **conifer forest** and watch the A-band line cores lift while the
  continuum barely moves — then check the two CO₂ bands, which do not move at all.
- Narrow the **ILS to 0.1×** and zoom into a single O₂ line to see the pressure-broadened
  line shape the ILS normally hides.
- Switch to **sorted**, focus the **O₂ A** band, set AOD to 0.5, then drag the
  **aerosol layer** up and down. The middle of the curve lifts strongly while the
  right-hand continuum barely moves — that separation is what makes aerosol layer
  height retrievable from an O₂ band. Then change **albedo scale** instead: the whole
  curve scales and the shape is untouched.
- Overlay the **sounding**, hit **load L2 state**, then detune XCO₂ and watch the
  residual grow.
- Drag **XCO₂** across its full 0–800 ppm range: the two CO₂ bands deepen steadily
  (mean transmittance 0.995 → 0.898 weak, 0.933 → 0.410 strong) and the A-band does
  not move by a single digit — `tau_co2` is zero there.
- Push **H₂O to 10×**: the strong-CO₂ band suffers most (0.58 → 0.40), which is why
  it is the band most contaminated by water vapour.
- Set **surface pressure to 0**: the A-band goes transparent, because O₂ scales with
  pressure. H₂O stays until you zero its own slider — it is a fixed tropospheric
  column in the data, not pressure-scaled. Zero both for a true vacuum, where
  transmittance is exactly 1.
- Narrow the **ILS to 0.1×**: absorption lines deepen towards their line-by-line
  depth. The weak-CO₂ band changes most (deepest line 0.23 → 0.95 in absorption),
  because its lines are narrowest relative to its ILS.

## Files

| File | Purpose |
|---|---|
| `index.html` | The whole app — physics, convolution, rendering, UI. No dependencies. |
| `validate.py` | NumPy reference model. Validates data, asserts the residual regression, writes `data/selftest_reference.json`. |
| `validate_scattering.py` | Monte Carlo check of the two-stream scattering model, over the full slider range. |
| `selftest_headless.js` | Runs the real `index.html` script under `gjs` with a DOM shim, for CI without a browser. |
| `data/*.json` | Inputs, mirrored unmodified from the upstream site. |
| `orig.html` | The original widget, for comparison. |

## About "load L2 state"

It sets XCO₂, SZA, VZA, surface pressure, the per-band retrieved albedo, and the L2
aerosol optical depth. Residuals against the sounding:

| band | rms | % of continuum | bias |
|---|---|---|---|
| O₂ A-band | 4.014 | 3.11 % | −1.13 |
| weak CO₂ | 1.071 | 3.58 % | −1.03 |
| strong CO₂ | 0.327 | 2.71 % | +0.01 |

**Do not expect these to go to zero, and the weak-CO₂ band least of all.** The L2
Full Physics retrieval carries far more machinery than this tool — a multi-layer
scattering atmosphere, retrieved aerosol types and profiles, temperature and water
profiles, per-footprint dispersion and ILS, zero-level offsets, polarisation. Loading
a handful of its posterior values into a single-layer forward model with one
Lambertian albedo per band cannot reproduce its radiances, and the residual is
dominated by that structural difference rather than by anything adjustable here.
The point of the button is to put you in a physically sensible neighbourhood so the
sliders do something meaningful, not to claim a fit.

Two specifics worth recording:

- **The L2 aerosol is loaded, but only because scattering is two-stream.** As pure
  extinction the same AOD drove the A-band 3× worse — rms 3.88 → 11.51, with the model
  11 W m⁻² µm⁻¹ sr⁻¹ too dark. With scattering it is roughly neutral (3.88 → 4.01) and
  slightly *helps* the weak-CO₂ band. The aerosol layer *height* is not in the L2 file,
  so 800 hPa is just a default.
- **`albedo_slope_l2_per_wn` is not applied.** The reference wavenumber it is defined
  against is not recorded in the file, and measured against the sounding it does not
  help: the weak-CO₂ residual is already flat (0.23 % trend across the band) with no
  slope and degrades to 1.44 % with it. Use the albedo slope slider to tilt by hand.

`validate.py` asserts the A-band fit quality and that two-stream beats extinction-only
by >2×, so neither can regress quietly.

## Attribution

The original **OCO-2 Spectra Explorer** was created by **Christian Frankenberg**
(Caltech/JPL). This project is an extended rewrite of that widget and would not exist
without it; the forward-model formulation, band selection, fixed-axis presentation and
all of the pre-computed line-by-line data come from that work.

None of the data in `data/` originates here. It is mirrored unmodified from
[`rsfralab.github.io/oco/data`](https://rsfralab.github.io/oco/), where the JSON
files were assembled; the underlying sources are:

| Input | Source |
|---|---|
| CO₂ / O₂ / H₂O optical depths | ABSCO v5.2 (JPL/AER tables, incl. line mixing and continua), on a US Standard Atmosphere 1976 profile |
| Solar irradiance | TSIS-1 Hybrid Solar Reference Spectrum (Coddington et al. 2021) × G. Toon's disk-integrated solar line transmittance |
| Surface reflectance | ECOSTRESS Spectral Library 1.0 (Meerdink et al. 2019; JHU Becknic and UCSB ASD spectra) |
| ILS and the example sounding | OCO-2 L1b `oco2_L1bScND_53642a_240801_B11205r` and L2 `oco2_LtCO2_240801_B11211Ar` (NASA GES DISC) |

Rayleigh optical depth uses the Bodhaine et al. (1999) parameterisation and is added
here, not present upstream.

Methods and references added by this version:

- Zeng, Z.-C., Natraj, V., Xu, F., Pongetti, T. J., Shia, R.-L., Kort, E. A., Toon,
  G. C., Sander, S. P., & Yung, Y. L. (2018). Constraining aerosol vertical profile in
  the boundary layer using hyperspectral measurements of oxygen absorption.
  *Geophysical Research Letters*, 45(19), 10772–10780.
  [doi:10.1029/2018GL079286](https://doi.org/10.1029/2018GL079286) — the spectral
  sorting method used by the **sorted** view.
- Meador, W. E., & Weaver, W. R. (1980). Two-stream approximations to radiative
  transfer in planetary atmospheres. *J. Atmos. Sci.*, 37, 630–643 — the beam solution
  used for the scattering layer.
- Joseph, J. H., Wiscombe, W. J., & Weinman, J. A. (1976). The delta-Eddington
  approximation for radiative flux transfer. *J. Atmos. Sci.*, 33, 2452–2459.
- Bodhaine, B. A., Wood, N. B., Dutton, E. G., & Slusser, J. R. (1999). On Rayleigh
  optical depth calculations. *J. Atmos. Oceanic Technol.*, 16, 1854–1861.

`orig.html` is a verbatim copy of the upstream widget, kept only for side-by-side
comparison. It carries no licence notice, so consider removing it from a public
repository — the table above and the link at the top of this file preserve the
attribution without redistributing the code.

## Verification

```
python3 validate.py             # data checks, ILS + L2 regressions, writes the reference
python3 validate_scattering.py  # Monte Carlo check of the scattering model
gjs selftest_headless.js        # runs the shipped JS against that reference
```

`validate.py` checks grid uniformity, array lengths, kernel normalisation, finiteness
and sounding overlap, then asserts the six residual numbers in the table above.

`selftest_headless.js` extracts the actual `<script>` from `index.html`, boots it
against a shimmed DOM with `?selftest=1`, and exits non-zero on failure. Current
status: **768 values across 16 parameter sets match to 4.9e-15** (float round-off;
the gate is 1e-9), and **44 render configurations** — every view × focus × overlay
combination, plus snapshot, zoom, and the slider extremes — draw without throwing.
The parameter sets cover the defaults, the L2 state, heavy aerosol, humid/low-pressure,
broadened and narrowed ILS with spectral shifts, a dark ocean scene, zero XCO₂, 800 ppm
XCO₂, 10× H₂O, zero surface pressure, a full vacuum, a 0.1× razor ILS, extinction-only
scattering, thick-high and strongly-absorbing aerosol layers, and three SIF cases.

The same check runs in the browser via `?selftest=1`, reported in a banner.

Not verified here: **the visual appearance**. This machine has no browser, so the
numerics and every render code path are tested, but the layout, colours and
interaction feel have not been seen rendered. Expect to want small cosmetic
adjustments on first look.
