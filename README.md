# OCO-2 Spectra Explorer — extended

An interactive, single-file visualisation of what the three OCO-2 bands see. It is a
substantially extended rewrite of
[`rsfralab.github.io/oco/widgets/spectra-explorer.html`](https://rsfralab.github.io/oco/widgets/spectra-explorer.html)
(mirrored here as `orig.html` for comparison), reading the same `data/` JSON files.

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
τ    = psurf/P₀ · (XCO₂/420 · τ_co2 + τ_o2) + h2o · τ_h2o + τ_ray + τ_aer
τ_ray = 0.008569 µm⁻⁴ (1 + 0.0113 µm⁻² + 0.00013 µm⁻⁴) · psurf/P₀     [Bodhaine 1999]
τ_aer = AOD₅₅₀ · (λ/550)^(−å)
m    = 1/cos(SZA) + 1/cos(VZA)
a(λ) = albScale · ECOSTRESS(λ) + albSlope · (ν − ν_centre),   ν = 10⁷/λ cm⁻¹
L    = E₀ · cos(SZA) · a(λ)/π · exp(−m τ)
L_c  = E₀ · cos(SZA) · a(λ)/π                                (continuum)
```

`τ_h2o` is deliberately **not** scaled by surface pressure — it is a tropospheric
profile, and the data documentation says to scale only `τ_co2` and `τ_o2`. The H₂O
slider is a separate column scale factor.

### Approximations you should know about

- **Aerosol and Rayleigh are pure extinction along the slant path.** There is no
  scattering source term, so dimming is overestimated — real scattering redirects
  photons back into the beam. This is not a substitute for a radiative-transfer
  solve; it is here to show the *wavelength dependence* of extinction.
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

**More physics.** Independent viewing zenith angle, H₂O column scale, aerosol AOD
with Ångström exponent, Rayleigh scattering, spectrally-resolved ECOSTRESS albedo
with scale and slope, ILS width and spectral-shift controls, and a noise/SNR model.
At the bundled sounding's L2 AOD of 0.071 and airmass 2.10, aerosol dims the A-band
10.2 % but the strong-CO₂ band only 3.9 %; Rayleigh dims them 4.4 % and 0.08 %.

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

## Controls

| Group | Controls |
|---|---|
| Geometry | solar zenith angle, viewing zenith angle |
| Atmosphere | XCO₂ (0–440 ppm), surface pressure (0–1030 hPa), H₂O column scale, aerosol AOD₅₅₀, Ångström exponent, Rayleigh |
| Surface | 5 ECOSTRESS surfaces + L2, spectral-albedo toggle, albedo scale, albedo slope |
| Instrument | ILS width (0.1–3× FWHM), spectral shift (±0.05 nm), SNR, noise |
| Sounding | overlay the real sounding, load its L2 state |

Views: **radiance** · **transmittance** (L/L_c, comparable across bands) ·
**optical depth** (−ln L/L_c, log axis) · **absorbers** (apparent τ per species).

Plot interaction: drag to zoom a band, shift-drag to pan, double-click to reset.

## Things to try

- Move **XCO₂**: lines deepen in the two CO₂ panels, the A-band holds still.
- Raise **AOD** to 0.5: the A-band dims far more than the strong-CO₂ band.
- Pick **snow**, switch to the strong-CO₂ band, toggle **spectral albedo**: the
  continuum tilts, because snow's reflectance falls 34 % across that band.
- Switch to **absorbers**: O₂ owns the A-band; H₂O contaminates the strong-CO₂ band.
- Overlay the **sounding**, hit **load L2 state**, then detune XCO₂ and watch the
  residual grow.
- Drag **XCO₂ to 0**: the two CO₂ bands lose their CO₂ lines entirely (what is left
  is H₂O), and the A-band does not move by a single digit — `tau_co2` is zero there.
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
| `selftest_headless.js` | Runs the real `index.html` script under `gjs` with a DOM shim, for CI without a browser. |
| `data/*.json` | Inputs, mirrored unmodified from the upstream site. |
| `orig.html` | The original widget, for comparison. |

## About "load L2 state"

It sets XCO₂, SZA, VZA, surface pressure and the per-band retrieved albedo, and
gives an A-band residual of 3.88 (3.0 % of continuum). Two things it deliberately
does **not** do:

- **It does not load `aod_total_l2` or enable Rayleigh.** The L2 albedo was retrieved
  with a full multiple-scattering model, in which aerosol does not simply attenuate
  the beam. Adding it here as pure extinction double-counts the loss and drives the
  A-band 3× worse — rms 3.88 → 11.51, with the model 11 W m⁻² µm⁻¹ sr⁻¹ too dark.
  Raise the aerosol slider yourself to see the effect; just don't expect a better fit.
- **It does not apply `albedo_slope_l2_per_wn`.** The reference wavenumber that slope
  is defined against is not recorded in the file, and measured against the sounding
  it does not help: the weak-CO₂ residual is already flat (0.23 % trend across the
  band) with no slope, and degrades to 1.44 % with it. Use the albedo slope slider to
  tilt by hand.

An earlier version did both, and also applied one band's slope to all three (a −4.8 %
spurious tilt across the A-band). `validate.py` now asserts the A-band fit quality so
this cannot regress.

## Attribution

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

`orig.html` is a verbatim copy of the upstream widget, kept only for side-by-side
comparison. It carries no licence notice, so consider removing it from a public
repository — the table above and the link at the top of this file preserve the
attribution without redistributing the code.

## Verification

```
python3 validate.py          # data checks + ILS regression + writes the reference
gjs selftest_headless.js     # runs the shipped JS, compares against that reference
```

`validate.py` checks grid uniformity, array lengths, kernel normalisation, finiteness
and sounding overlap, then asserts the six residual numbers in the table above.

`selftest_headless.js` extracts the actual `<script>` from `index.html`, boots it
against a shimmed DOM with `?selftest=1`, and exits non-zero on failure. Current
status: **528 values across 11 parameter sets match to 4.9e-15** (float round-off;
the gate is 1e-9), and **39 render configurations** — every view × focus × overlay
combination, plus snapshot, zoom, and the slider extremes — draw without throwing.
The parameter sets cover the defaults, the L2 state, heavy aerosol, humid/low-pressure,
broadened and narrowed ILS with spectral shifts, a dark ocean scene, zero XCO₂, zero
surface pressure, a full vacuum, and a 0.1× razor ILS.

The same check runs in the browser via `?selftest=1`, reported in a banner.

Not verified here: **the visual appearance**. This machine has no browser, so the
numerics and every render code path are tested, but the layout, colours and
interaction feel have not been seen rendered. Expect to want small cosmetic
adjustments on first look.
