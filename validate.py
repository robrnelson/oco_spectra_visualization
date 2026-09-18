#!/usr/bin/env python3
"""
Reference implementation and validation for the OCO-2 Spectra Explorer (extended).

This file is the numerics contract. index.html implements the same algorithms in
JavaScript; `python3 validate.py` checks the input data, reproduces the residual
regression numbers, and writes data/selftest_reference.json so that
index.html?selftest=1 can prove the JS matches this reference to <1e-9 relative.

Usage:
    python3 validate.py            # validate data, run regressions, write reference
    python3 validate.py --quiet    # only report failures
"""

import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
BANDS = ["aband", "wco2", "sco2"]
P0 = 1013.25          # reference surface pressure, hPa
XCO2_REF = 420.0      # ppm, the XCO2 the tau_co2 tables were computed at
AOD_REF_NM = 550.0    # wavelength at which the aerosol optical depth is specified

_cache = {}


def load(name):
    if name not in _cache:
        with open(os.path.join(DATA, name)) as fh:
            _cache[name] = json.load(fh)
    return _cache[name]


def band_data(band):
    """Grid + per-absorber optical depths + solar irradiance + ILS for one band."""
    key = f"_band_{band}"
    if key in _cache:
        return _cache[key]
    tau = load(f"tau_{band}.json")
    sol = load(f"solar_{band}.json")
    ils = load(f"ils_{band}.json")
    m = tau["meta"]
    n = m["n_points"]
    d = {
        "n": n,
        "step": m["wl_step_nm"],
        "wl0": m["wl_start_nm"],
        "wl": m["wl_start_nm"] + np.arange(n) * m["wl_step_nm"],
        "tau_co2": np.asarray(tau["tau_co2"], float),
        "tau_o2": np.asarray(tau["tau_o2"], float),
        "tau_h2o": np.asarray(tau["tau_h2o"], float),
        "e0": np.asarray(sol["e0"], float),
        "ils_dl": np.asarray(ils["delta_lambda_nm"], float),
        "ils_r": np.asarray(ils["response"], float),
        "meta": m,
    }
    _cache[key] = d
    return d


# ---------------------------------------------------------------------------
# Numerics core -- must match index.html exactly
# ---------------------------------------------------------------------------

def make_kernel(dl, r, step, stretch=1.0, shift_nm=0.0, legacy_centre=False):
    """Resample an ILS onto the tau grid spacing.

    The ILS is a function of dl = lambda' - lambda. Stretching by `stretch`
    scales its width; `shift_nm` moves the band centre the instrument reports.
    The returned kernel is indexed on offsets (j - centre) * step, so placing
    index `centre` at the output wavelength puts dl = 0 there -- which is what
    the original widget got wrong by assuming centre == (len-1)>>1.

    legacy_centre=True reproduces that original behaviour for regression tests.
    """
    if legacy_centre:
        k = r.copy()
        return k / k.sum(), (len(k) - 1) >> 1

    # integer offsets (in grid steps) whose dl falls inside the tabulated support
    jmin = math.ceil((stretch * dl[0] + shift_nm) / step)
    jmax = math.floor((stretch * dl[-1] + shift_nm) / step)
    j = np.arange(jmin, jmax + 1)
    # sample the stretched, shifted ILS: R((offset - shift) / stretch)
    k = np.interp((j * step - shift_nm) / stretch, dl, r, left=0.0, right=0.0)
    s = k.sum()
    if s <= 0:
        raise ValueError("degenerate ILS kernel")
    return k / s, -jmin


def conv_at(x, k, centre, idx):
    """Edge-normalized correlation of x with kernel k, evaluated only at idx.

    Evaluating at just the output positions is bit-identical to convolving the
    whole grid and then subsampling, but costs far less (measured ~7x in JS).
    """
    n = len(x)
    K = len(k)
    out = np.empty(len(idx), float)
    for o, i in enumerate(idx):
        lo = i - centre
        j0 = lo if lo > 0 else 0
        j1 = lo + K - 1
        if j1 > n - 1:
            j1 = n - 1
        kk = k[j0 - lo:j1 - lo + 1]
        out[o] = float(np.dot(x[j0:j1 + 1], kk) / kk.sum())
    return out


def rayleigh_tau(wl_nm):
    """Bodhaine et al. (1999) sea-level Rayleigh optical depth, per wavelength."""
    um = wl_nm / 1000.0
    return 0.008569 * um ** -4 * (1 + 0.0113 * um ** -2 + 0.00013 * um ** -4)


# ---------------------------------------------------------------------------
# Scattering: a single two-stream layer coupled to the Lambertian surface.
#
# Pure extinction is badly wrong for a scattering layer, because photons removed
# from the direct beam are not lost -- they still reach the surface diffusely and
# still come back up. Treating the aerosol as extinction only dimmed the A-band
# 10.2 % at the sounding's AOD; the two-stream treatment below dims it 0.35 %,
# and agrees with a Monte Carlo reference to <0.009 in system albedo.
# ---------------------------------------------------------------------------

def delta_scale(tau, w, g):
    """Delta-Eddington scaling: strip the forward-scattering spike."""
    f = g * g
    return (1 - w * f) * tau, np.clip(w * (1 - f) / (1 - w * f), 0.0, 1 - 1e-7), g / (1 + g)


def two_stream_beam(tau, w, g, mu0):
    """Meador & Weaver (1980), Eddington closure: collimated beam on a layer with
    a black lower boundary. Returns (R_beam, T_beam_total incl. direct)."""
    t, w, g = delta_scale(np.asarray(tau, float), np.asarray(w, float),
                          np.asarray(g, float))
    g1 = (7 - w * (4 + 3 * g)) / 4
    g2 = -(1 - w * (4 - 3 * g)) / 4
    lam = np.sqrt(np.maximum(g1 * g1 - g2 * g2, 1e-30))
    # lam*mu0 == 1 is a removable singularity of this closed form; nudge mu0 off
    # it where it bites (standard practice in RT codes). 0.3 % in mu0 is nothing.
    mu = np.where(np.abs(1 - (lam * mu0) ** 2) < 1e-4, mu0 * (1 - 3e-3), mu0)
    g3 = (2 - 3 * g * mu) / 4
    g4 = 1 - g3
    a1 = g1 * g4 + g2 * g3
    a2 = g1 * g3 + g2 * g4
    el, eml = np.exp(np.minimum(lam * t, 700)), np.exp(-lam * t)
    e0 = np.exp(-np.minimum(t / mu, 700))
    den = (lam + g1) * el + (lam - g1) * eml
    q = 1 - (lam * mu) ** 2
    R = (w / (q * den)) * ((1 - lam * mu) * (a2 + lam * g3) * el
                           - (1 + lam * mu) * (a2 - lam * g3) * eml
                           - 2 * lam * (g3 - a2 * mu) * e0)
    Td = -(w / (q * den)) * ((1 + lam * mu) * (a1 + lam * g4) * el * e0
                             - (1 - lam * mu) * (a1 - lam * g4) * eml * e0
                             - 2 * lam * (g4 + a1 * mu))
    return np.clip(R, 0.0, 1.0), np.clip(Td + e0, 0.0, 1.0)


def two_stream_diffuse(tau, w, g):
    """Classical two-stream for diffuse illumination: spherical albedo and
    diffuse transmittance (used for the upward path off a Lambertian surface)."""
    t, w, g = delta_scale(np.asarray(tau, float), np.asarray(w, float),
                          np.asarray(g, float))
    with np.errstate(all="ignore"):
        u = np.sqrt((1 - w) / (1 - w * g))
        rinf = (1 - u) / (1 + u)
        k = np.sqrt(3 * (1 - w) * (1 - w * g))
        e = np.exp(-2 * k * t)
        R = rinf * (1 - e) / (1 - rinf ** 2 * e)
        T = (1 - rinf ** 2) * np.exp(-k * t) / (1 - rinf ** 2 * e)
    # Conservative-scattering limit (w -> 1), where the closed form above is 0/0.
    # The threshold must be LOOSER than delta_scale's 1-1e-7 clip, or this branch
    # can never trigger -- which matters because Rayleigh alone gives w = 1 exactly.
    qc = 0.75 * (1 - g) * t
    cons = w >= 1 - 1e-6
    R = np.where(cons, qc / (1 + qc), R)
    T = np.where(cons, 1 / (1 + qc), T)
    return np.clip(R, 0.0, 1.0), np.clip(T, 0.0, 1.0)


def scattering_layer(band, d, state):
    """Mix Rayleigh and aerosol into one effective scattering layer.

    Standard optical-property mixing: tau adds, and omega/g are tau-weighted.
    `f` is the fraction of the well-mixed gas column lying ABOVE the layer, which
    sets how much gas the atmospherically-scattered photons actually traverse --
    the mechanism behind the aerosol-induced XCO2 bias. Rayleigh scattering is
    distributed proportional to pressure, so its mean scattering level sits at
    psurf/2, i.e. f = 0.5.
    """
    ps = state["psurf"] / P0
    t_ray = rayleigh_tau(d["wl"]) * ps if state["rayleigh"] else np.zeros(d["n"])
    t_aer = (state["aod"] * (d["wl"] / AOD_REF_NM) ** -state["angstrom"]
             if state["aod"] > 0 else np.zeros(d["n"]))
    t_tot = t_ray + t_aer
    f_aer = np.clip(state["aerP"] / max(state["psurf"], 1e-9), 0.0, 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = np.where(t_tot > 0, (t_ray * 1.0 + t_aer * state["ssa"]) / t_tot, 1.0)
        g = np.where(t_tot > 0, (t_ray * 0.0 + t_aer * state["asym"]) / t_tot, 0.0)
        f = np.where(t_tot > 0, (t_ray * 0.5 + t_aer * f_aer) / t_tot, 1.0)
    return t_tot, w, g, f, t_ray, t_aer


def albedo_spectrum(band, d, state):
    """Surface albedo across the band: spectral, band-mean, or L2 constant."""
    wl = d["wl"]
    surf = state["surface"]
    if surf == "l2":
        # per-band L2 albedo as a constant; albedo_slope_l2_per_wn is deliberately
        # not applied (undocumented reference wavenumber, and it does not improve
        # the fit) -- see albedoBase() in index.html for the measurements
        base = np.full(len(wl), load("sounding_oco2.json")["meta"]["albedo_l2"][band])
    else:
        alb = load("albedo_surfaces.json")
        if state["spectralAlbedo"]:
            base = np.interp(wl, np.asarray(alb["wl"], float),
                             np.asarray(alb[surf], float))
        else:
            base = np.full(len(wl), alb["band_means"][surf][band])
    # slope is per cm-1, referenced to the band centre, matching albedo_slope_l2_per_wn
    nu = 1e7 / wl
    nu_c = 1e7 / (0.5 * (wl[0] + wl[-1]))
    a = state["albScale"] * base + state["albSlope"] * (nu - nu_c)
    return np.clip(a, 0.0, None)


def optical_depths(band, d, state):
    """Per-absorber vertical optical depths for the current state."""
    ps = state["psurf"] / P0
    t_tot, _, _, _, t_ray, t_aer = scattering_layer(band, d, state)
    return {
        # tau_co2/tau_o2 are tabulated at P0 and scale with surface pressure;
        # tau_h2o is a tropospheric profile and deliberately does NOT (see data meta)
        "co2": ps * (state["xco2"] / XCO2_REF) * d["tau_co2"],
        "o2": ps * d["tau_o2"],
        "h2o": state["h2o"] * d["tau_h2o"],
        "ray": t_ray,
        "aer": t_aer,
    }


def forward(band, state):
    """Line-by-line radiance and continuum on the fine grid (unconvolved).

    With scattering on, the signal splits into two paths with *different* gas
    path lengths, so it no longer factorises as continuum x exp(-m tau):

        L = (E0 mu0 / pi) [ R_beam         . exp(-m f tau_gas)      <- atmospheric
                          + T_beam a Tdiff exp(-m tau_gas) / (1-a s) ] <- surface
                            \\_______________________________/
                                       multiple surface-layer reflections

    The first term only crosses the gas *above* the scattering layer, which is
    exactly why aerosol biases retrieved XCO2. With aod = 0 and Rayleigh off this
    reduces identically to the old L = Lc exp(-m tau).
    """
    d = band_data(band)
    mu0 = math.cos(math.radians(state["sza"]))
    m = 1.0 / mu0 + 1.0 / math.cos(math.radians(state["vza"]))
    taus = optical_depths(band, d, state)
    tau_gas = taus["co2"] + taus["o2"] + taus["h2o"]
    a = albedo_spectrum(band, d, state)
    solar = d["e0"] * mu0 / math.pi

    t_sc, w_sc, g_sc, f_sc, _, _ = scattering_layer(band, d, state)
    # scalar condition so index.html can mirror the branch exactly
    have_sc = state["aod"] > 0 or (state["rayleigh"] and state["psurf"] > 0)
    if not state.get("scatter", True) or not have_sc:
        # extinction only: the old behaviour, kept so the difference is visible
        Lc = solar * a
        L = Lc * np.exp(-m * (tau_gas + t_sc))
        return L, Lc, m, taus

    Rb, Tb = two_stream_beam(t_sc, w_sc, g_sc, mu0)
    s, Td = two_stream_diffuse(t_sc, w_sc, g_sc)
    surf = Tb * a * Td / np.maximum(1 - a * s, 1e-12)
    Lc = solar * (Rb + surf)
    L = solar * (Rb * np.exp(-m * f_sc * tau_gas) + surf * np.exp(-m * tau_gas))
    return L, Lc, m, taus


def model_at(band, state, idx, legacy_centre=False):
    """Convolved radiance and continuum at integer grid indices `idx`.

    legacy_centre is Python-only, used by regression() to show what the original
    widget's array-midpoint assumption costs. The widget itself has no such mode.
    """
    d = band_data(band)
    L, Lc, m, _ = forward(band, state)
    k, centre = make_kernel(d["ils_dl"], d["ils_r"], d["step"],
                            state["ilsStretch"], state["shiftNm"], legacy_centre)
    return conv_at(L, k, centre, idx), conv_at(Lc, k, centre, idx), m


def default_state():
    return dict(xco2=420.0, sza=30.0, vza=0.0, psurf=1013.25, h2o=1.0,
                aod=0.0, angstrom=1.0, rayleigh=False, surface="dry_grass",
                albScale=1.0, albSlope=0.0, spectralAlbedo=False,
                ilsStretch=1.0, shiftNm=0.0,
                scatter=True, ssa=0.95, asym=0.65, aerP=800.0)


def l2_state():
    """State matching the L2 retrieval of the bundled sounding."""
    M = load("sounding_oco2.json")["meta"]
    s = default_state()
    s.update(xco2=M["xco2_l2_ppm"], sza=M["sza"], vza=M["vza"],
             psurf=M["surface_pressure_hpa"], surface="l2",
             # with two-stream scattering the L2 aerosol is ~neutral on the fit,
             # so loading it is both faithful and harmless (it was not, as pure
             # extinction: that drove the A-band 3x worse)
             aod=M["aod_total_l2"], angstrom=1.0, rayleigh=True, scatter=True)
    return s


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

class Report:
    def __init__(self, quiet=False):
        self.quiet, self.failures = quiet, []

    def ok(self, msg):
        if not self.quiet:
            print(f"  ok    {msg}")

    def fail(self, msg):
        self.failures.append(msg)
        print(f"  FAIL  {msg}")

    def check(self, cond, msg):
        self.ok(msg) if cond else self.fail(msg)


def validate_data(rep):
    print("\n[1] data files")
    for band in BANDS:
        d = band_data(band)
        m = d["meta"]
        rep.check(len(d["tau_co2"]) == d["n"] and len(d["tau_o2"]) == d["n"]
                  and len(d["tau_h2o"]) == d["n"] and len(d["e0"]) == d["n"],
                  f"{band}: all arrays have n_points={d['n']}")
        lo, hi = m["band_range_nm"]
        rep.check(abs(d["wl"][0] - lo) < 1e-9 and abs(d["wl"][-1] - hi) < 1e-6,
                  f"{band}: uniform grid spans {lo}-{hi} nm")
        finite = all(np.isfinite(d[k]).all()
                     for k in ("tau_co2", "tau_o2", "tau_h2o", "e0"))
        rep.check(finite, f"{band}: no NaN/Inf")
        rep.check((d["tau_co2"] >= 0).all() and (d["tau_o2"] >= 0).all()
                  and (d["tau_h2o"] >= 0).all(), f"{band}: optical depths non-negative")
        rep.check((d["e0"] > 0).all(), f"{band}: solar irradiance positive")

        ils = load(f"ils_{band}.json")
        r = d["ils_r"]
        rep.check(len(r) == len(d["ils_dl"]), f"{band}: ILS response/delta_lambda same length")
        rep.check((r >= 0).all(), f"{band}: ILS response non-negative")
        step_ils = ils["meta"]["spacing_nm"]
        rep.check(abs(step_ils - d["step"]) < 1e-12,
                  f"{band}: ILS spacing {step_ils} matches tau grid step {d['step']}")
        # kernel built at stretch=1 must be unit-sum and centred on the ILS peak
        k, centre = make_kernel(d["ils_dl"], d["ils_r"], d["step"])
        rep.check(abs(k.sum() - 1.0) < 1e-12, f"{band}: kernel sums to 1")
        rep.check(abs(int(k.argmax()) - centre) <= 2,
                  f"{band}: kernel peak within 2 samples of centre index {centre}")

        # the mis-centering the original suffers from
        legacy = (len(d["ils_r"]) - 1) >> 1
        off = (centre - legacy) * d["step"]
        fwhm = ils["meta"]["fwhm_nm"]
        print(f"        {band}: true centre idx {centre}, array midpoint {legacy}"
              f" -> legacy offset {off:+.4f} nm ({100*abs(off)/fwhm:.0f}% of FWHM {fwhm:.4f})")

    snd = load("sounding_oco2.json")
    for band in BANDS:
        sd, d = snd[band], band_data(band)
        swl = np.asarray(sd["wl"], float)
        rep.check(len(swl) == len(sd["radiance"]) == len(sd["good"]),
                  f"{band}: sounding arrays same length ({len(swl)})")
        rep.check(swl.max() > d["wl"][0] and swl.min() < d["wl"][-1],
                  f"{band}: sounding {swl.min():.1f}-{swl.max():.1f} nm overlaps tau grid")
        rep.check(sum(sd["good"]) > 0, f"{band}: sounding has {sum(sd['good'])} good samples")

    alb = load("albedo_surfaces.json")
    awl = np.asarray(alb["wl"], float)
    for s in ("ocean", "desert", "snow", "conifer", "dry_grass"):
        v = np.asarray(alb[s], float)
        rep.check(len(v) == len(awl) and (v >= 0).all() and (v <= 1.2).all(),
                  f"albedo {s}: {len(v)} points in [0,1.2]")
    for band in BANDS:
        lo, hi = alb["meta"]["band_ranges_nm"][band]
        rep.check(awl.min() <= lo and awl.max() >= hi,
                  f"albedo spectrum covers {band} ({lo}-{hi} nm)")


def sounding_residuals(state, legacy_centre=False):
    """RMS of (sounding - model) per band, at the sounding's good samples."""
    snd = load("sounding_oco2.json")
    M = snd["meta"]
    pf = M.get("polarization_factor") or 0.5
    out = {}
    for band in BANDS:
        d = band_data(band)
        sd = snd[band]
        swl = np.asarray(sd["wl"], float)
        rad = np.asarray(sd["radiance"], float) / pf
        good = np.asarray(sd["good"], bool)
        # keep a margin from the grid edges so the kernel is not truncated
        g = good & (swl >= d["wl"][0] + 0.4) & (swl <= d["wl"][-1] - 0.4)
        idx = np.clip(np.round((swl[g] - d["wl0"]) / d["step"]).astype(int), 0, d["n"] - 1)
        Lm, Lcm, _ = model_at(band, state, idx, legacy_centre)
        r = rad[g] - Lm
        out[band] = dict(rms=float(np.sqrt((r ** 2).mean())), bias=float(r.mean()),
                         cont=float(Lcm.max()), n=int(g.sum()))
    return out


def regression(rep):
    """The ILS-centering fix must reproduce the improvement measured while planning."""
    print("\n[2] ILS centering regression (model vs real sounding, L2 parameters)")
    # aerosol held out so this isolates the ILS-centring effect; l2_state() itself
    # now carries the L2 aerosol, which would move these constants
    st = dict(l2_state(), aod=0.0, rayleigh=False)
    legacy = sounding_residuals(st, legacy_centre=True)
    fixed = sounding_residuals(st, legacy_centre=False)
    # (band, expected legacy rms, expected fixed rms)
    expect = [("aband", 6.599, 3.883), ("wco2", 1.152, 1.134), ("sco2", 0.330, 0.330)]
    for band, e_leg, e_fix in expect:
        L, F = legacy[band]["rms"], fixed[band]["rms"]
        print(f"        {band}: legacy rms {L:.3f} -> fixed {F:.3f}"
              f"  ({100*L/legacy[band]['cont']:.2f}% -> {100*F/fixed[band]['cont']:.2f}% of continuum)")
        rep.check(abs(L - e_leg) < 0.02, f"{band}: legacy rms {L:.3f} ~ expected {e_leg}")
        rep.check(abs(F - e_fix) < 0.02, f"{band}: fixed  rms {F:.3f} ~ expected {e_fix}")
    rep.check(fixed["aband"]["rms"] < 0.65 * legacy["aband"]["rms"],
              "A-band residual improves by >35% with correct ILS centering")
    # Where the ILS is already centred (sco2) the two paths agree to round-off: the
    # tabulated delta_lambda_nm values are stored rounded, so resampling them onto the
    # exact grid lattice perturbs the kernel by ~2e-5 relative. Allow that, catch more.
    for band in BANDS:
        rel = (fixed[band]["rms"] - legacy[band]["rms"]) / legacy[band]["rms"]
        rep.check(rel <= 1e-3,
                  f"{band}: correct centering is not worse than legacy "
                  f"(rel change {rel:+.2e})")


def scattering_checks(rep):
    """Cheap invariants for the two-stream layer. The quantitative accuracy check
    against Monte Carlo lives in validate_scattering.py."""
    print("\n[3] scattering model")

    # 1. with no scattering layer, the two paths must collapse to the old model
    for band in BANDS:
        st = default_state()
        L1, Lc1, _, _ = forward(band, dict(st, scatter=True))
        L2, Lc2, _, _ = forward(band, dict(st, scatter=False))
        rep.check(np.array_equal(L1, L2) and np.array_equal(Lc1, Lc2),
                  f"{band}: aod=0 + Rayleigh off is bit-identical with/without scattering")

    # 2. conservative scattering conserves energy for a black surface
    for tau in (0.05, 0.5, 2.0):
        R, T = two_stream_diffuse(tau, 1.0, 0.65)
        rep.check(abs(float(R) + float(T) - 1.0) < 1e-12,
                  f"diffuse R+T = 1 at tau={tau} for omega=1 (energy conserved)")
        Rb, Tb = two_stream_beam(tau, 1.0, 0.65, 0.9)
        rep.check(float(Rb) + float(Tb) <= 1.0 + 1e-9,
                  f"beam R+T <= 1 at tau={tau} (no energy created)")

    # 3. brightening over a dark surface must be monotonic in AOD
    d = band_data("aband")
    prev = -1.0
    mono = True
    for aod in (0.0, 0.1, 0.3, 0.6, 1.0):
        st = dict(default_state(), surface="ocean", aod=aod, ssa=1.0, rayleigh=False)
        _, Lc, _, _ = forward("aband", st)
        v = float(Lc.mean())
        if v < prev - 1e-12:
            mono = False
        prev = v
    rep.check(mono, "continuum over dark ocean rises monotonically with scattering AOD")

    # 4. the path-shortening that biases XCO2: higher layer -> less gas traversed
    fs = []
    for p in (1000.0, 700.0, 400.0, 150.0):
        _, _, _, f, _, _ = scattering_layer(
            "sco2", band_data("sco2"),
            dict(default_state(), aod=0.5, rayleigh=False, aerP=p))
        fs.append(float(np.mean(f)))
    rep.check(all(fs[i] > fs[i+1] for i in range(len(fs)-1)),
              f"gas fraction above the layer falls as it rises: "
              f"{' > '.join(f'{v:.2f}' for v in fs)}")

    # 5. an absorbing aerosol must darken, a conservative one must not (dark scene)
    st = dict(default_state(), surface="ocean", aod=0.8, rayleigh=False)
    _, Lc_w1, _, _ = forward("aband", dict(st, ssa=1.0))
    _, Lc_w5, _, _ = forward("aband", dict(st, ssa=0.5))
    rep.check(float(Lc_w5.mean()) < float(Lc_w1.mean()),
              "absorbing aerosol (omega=0.5) is darker than conservative (omega=1)")


SORT_SAT, SORT_CONT = 0.05, 0.80   # clear-sky transmittance region cuts (display)


def sorted_regions(band, state, idx):
    """Spectral sorting after Zeng et al. (2018): order channels by ascending
    CLEAR-SKY radiance (same state, AOD = 0) and apply that order to the scene.
    Returns mean radiance in the continuum / intermediate / saturated regions.
    Mirrors the "sorted" view in index.html."""
    L, _, _ = model_at(band, state, idx)
    L0, Lc0, _ = model_at(band, dict(state, aod=0.0), idx)
    L, L0, Lc0 = np.asarray(L), np.asarray(L0), np.asarray(Lc0)
    perm = np.argsort(L0, kind="stable")
    t0 = np.where(Lc0 > 0, L0 / Lc0, 0.0)[perm]
    Ls = L[perm]
    pick = lambda msk: float(Ls[msk].mean()) if msk.any() else float("nan")
    return (pick(t0 >= SORT_CONT),
            pick((t0 >= SORT_SAT) & (t0 < SORT_CONT)),
            pick(t0 < SORT_SAT))


def sorting_checks(rep):
    """The sorted view must reproduce the discriminants the method relies on:
    the continuum tracks AOD, the intermediate lines track aerosol layer height,
    and surface albedo scales the level without changing the shape."""
    print("\n[5] spectral sorting (Zeng et al. 2018)")
    band = "aband"
    d = band_data(band)
    idx = np.arange(int(0.05 * d["n"]), int(0.95 * d["n"]), 4)
    base = dict(default_state(), surface="desert", sza=45.0)

    # 1. continuum rises monotonically with AOD
    cont = [sorted_regions(band, dict(base, aod=a, rayleigh=True), idx)[0]
            for a in (0.0, 0.1, 0.3, 0.6, 1.0)]
    rep.check(all(cont[i] < cont[i+1] for i in range(len(cont)-1)),
              "continuum mean rises monotonically with AOD: "
              + " < ".join(f"{v:.2f}" for v in cont))

    # 2. raising the layer lifts the intermediate lines, and lifts them far more
    #    than the continuum -- this separation is what makes ALH retrievable
    mids, conts = [], []
    for p in (1000.0, 850.0, 700.0, 500.0, 300.0, 150.0):
        c, m, _ = sorted_regions(band, dict(base, aod=0.5, rayleigh=True, aerP=p), idx)
        conts.append(c); mids.append(m)
    rep.check(all(mids[i] < mids[i+1] for i in range(len(mids)-1)),
              "intermediate mean rises monotonically as the aerosol layer rises: "
              + " < ".join(f"{v:.2f}" for v in mids))
    d_mid = (mids[-1] - mids[0]) / mids[0]
    d_cont = (conts[-1] - conts[0]) / conts[0]
    print(f"        1000 -> 150 hPa at AOD 0.5: intermediate {100*d_mid:+.1f} %, "
          f"continuum {100*d_cont:+.1f} %")
    rep.check(d_mid > 8 * d_cont,
              f"layer height moves the intermediate lines >8x more than the continuum "
              f"({100*d_mid:+.1f} % vs {100*d_cont:+.1f} %)")

    # 3. surface albedo scales the level but not the shape (the AOD/albedo degeneracy)
    ratios = []
    for sc in (0.6, 1.0, 1.6):
        c, m, _ = sorted_regions(band, dict(base, aod=0.3, rayleigh=True, albScale=sc), idx)
        ratios.append(m / c)
    spread = (max(ratios) - min(ratios)) / np.mean(ratios)
    print(f"        intermediate/continuum ratio across albScale 0.6-1.6: "
          + ", ".join(f"{r:.4f}" for r in ratios))
    rep.check(spread < 0.01,
              f"albedo scales the sorted curve without changing its shape "
              f"(ratio spread {100*spread:.2f} % < 1 %)")


def l2_match(rep):
    """The widget's "load L2 state" must actually fit the sounding.

    An earlier version also switched on Rayleigh and the L2 aerosol optical depth,
    and applied one band's albedo slope to all three. Because aerosol here is pure
    extinction while the L2 albedo was retrieved with full multiple scattering,
    that double-counted the loss and drove the A-band 3x worse.
    """
    print("\n[6] \"load L2 state\" fit quality")
    st = l2_state()                                    # what the button now sets
    good = sounding_residuals(st)
    for band in BANDS:
        g = good[band]
        print(f"        {band:6s} rms {g['rms']:6.3f}  ({100*g['rms']/g['cont']:4.2f} % of "
              f"continuum, bias {g['bias']:+.3f})")
    rep.check(good["aband"]["rms"] < 4.2,
              f"A-band rms {good['aband']['rms']:.3f} < 4.2 at the L2 state")
    rep.check(abs(good["aband"]["bias"]) < 1.5,
              f"A-band bias {good['aband']['bias']:+.3f} within +/-1.5 (not systematically dark)")

    # show why aerosol/Rayleigh are deliberately left off
    M = load("sounding_oco2.json")["meta"]
    # the same aerosol treated as pure extinction is what used to wreck the A-band
    bad = sounding_residuals(dict(st, scatter=False))
    print(f"        same state with extinction-only scattering: A-band rms "
          f"{bad['aband']['rms']:.3f} (bias {bad['aband']['bias']:+.3f}) — "
          f"{bad['aband']['rms']/good['aband']['rms']:.1f}x worse")
    rep.check(bad["aband"]["rms"] > 2.0 * good["aband"]["rms"],
              "two-stream scattering beats extinction-only by >2x on the A-band")


def write_selftest(rep):
    """Sample the model over several states so the JS can be checked against it."""
    print("\n[7] selftest reference")
    M = load("sounding_oco2.json")["meta"]
    cases = []

    def add(name, **over):
        st = default_state()
        st.update(over)
        cases.append((name, st))

    add("default")
    cases.append(("l2", l2_state()))
    add("aerosol_heavy", aod=0.5, angstrom=1.5, rayleigh=True)
    add("humid_lowp", h2o=3.0, psurf=700.0, xco2=440.0, sza=65.0, vza=40.0)
    add("ils_broad_shift", ilsStretch=3.0, shiftNm=0.05, surface="snow",
        spectralAlbedo=True)
    add("ils_narrow_negshift", ilsStretch=0.5, shiftNm=-0.05, surface="conifer",
        spectralAlbedo=True, albScale=1.4, albSlope=-1.2e-4)
    add("ocean_dark", surface="ocean", spectralAlbedo=True, aod=0.25, rayleigh=True,
        xco2=380.0, psurf=1030.0)
    # the extremes the sliders now reach
    add("xco2_zero", xco2=0.0)
    add("vacuum", psurf=0.0, xco2=0.0, h2o=0.0)
    add("psurf_zero_humid", psurf=0.0)
    add("ils_razor", ilsStretch=0.1, surface="desert", spectralAlbedo=True)
    # newly reachable slider space + both scattering treatments
    add("xco2_800", xco2=800.0)
    add("h2o_10x", h2o=10.0, xco2=800.0)
    add("extinction_only", aod=0.4, angstrom=1.2, rayleigh=True, scatter=False)
    add("scatter_thick_high", aod=1.5, angstrom=1.2, rayleigh=True, scatter=True,
        ssa=0.99, asym=0.80, aerP=200.0, surface="conifer", spectralAlbedo=True)
    add("scatter_absorbing_low", aod=0.8, rayleigh=True, scatter=True,
        ssa=0.50, asym=0.10, aerP=1000.0, sza=65.0, vza=35.0)

    ref = {
        "generated_by": "validate.py",
        "note": ("Reference model values for index.html?selftest=1. Radiance/continuum "
                 "are convolved values at the given integer grid indices. Noise is "
                 "excluded: it is a display feature and uses a PRNG, not part of the "
                 "numerics contract."),
        "constants": {"P0": P0, "XCO2_REF": XCO2_REF, "AOD_REF_NM": AOD_REF_NM},
        "tolerance_rel": 1e-9,
        "cases": [],
    }
    for name, st in cases:
        entry = {"name": name, "state": st, "bands": {}}
        for band in BANDS:
            d = band_data(band)
            # 8 indices spread across the band, kept clear of the grid edges
            idx = np.linspace(int(0.08 * d["n"]), int(0.92 * d["n"]), 8).astype(int)
            L, Lc, m = model_at(band, st, idx)
            # same flags model_at used, or the recorded kernel metadata would
            # describe a different kernel than the radiances above
            k, centre = make_kernel(d["ils_dl"], d["ils_r"], d["step"],
                                    st["ilsStretch"], st["shiftNm"])
            entry["bands"][band] = {
                "idx": [int(i) for i in idx],
                "wl_nm": [float(d["wl0"] + i * d["step"]) for i in idx],
                "L": [float(v) for v in L],
                "Lc": [float(v) for v in Lc],
                "airmass": float(m),
                "kernel_len": int(len(k)),
                "kernel_centre": int(centre),
            }
        ref["cases"].append(entry)

    out = os.path.join(DATA, "selftest_reference.json")
    with open(out, "w") as fh:
        json.dump(ref, fh, separators=(",", ":"))
    size = os.path.getsize(out)
    rep.check(size > 0, f"wrote {os.path.relpath(out, HERE)} "
                        f"({len(ref['cases'])} cases, {size/1024:.1f} kB)")

    # sanity: the cases must actually differ somewhere, or the selftest proves
    # nothing. Not per-band: XCO2 cannot move the A-band at all (tau_co2 is zero
    # there), which is itself the point of the O2 band.
    base = ref["cases"][0]["bands"]
    for c in ref["cases"][1:]:
        moved = [b for b in BANDS
                 if any(abs(x - y) > 1e-6
                        for x, y in zip(base[b]["L"], c["bands"][b]["L"]))]
        rep.check(bool(moved),
                  f"case '{c['name']}' differs from 'default' (bands: "
                  f"{','.join(moved) or 'NONE'})")


def physics_summary():
    """Magnitudes of the newly added terms, for the record."""
    print("\n[8] magnitudes of the added physics (airmass 2.10, AOD 0.071, alpha 1)")
    M = load("sounding_oco2.json")["meta"]
    m = 1 / math.cos(math.radians(M["sza"])) + 1 / math.cos(math.radians(M["vza"]))
    for band, lam in (("aband", 765.0), ("wco2", 1605.0), ("sco2", 2062.0)):
        tr = rayleigh_tau(lam) * M["surface_pressure_hpa"] / P0
        ta = M["aod_total_l2"] * (lam / AOD_REF_NM) ** -1.0
        print(f"        {band:6s} {lam:6.0f} nm: Rayleigh dims {100*(1-math.exp(-m*tr)):5.2f}%,"
              f" aerosol dims {100*(1-math.exp(-m*ta)):5.2f}%")
    alb = load("albedo_surfaces.json")
    awl = np.asarray(alb["wl"], float)
    print("        spectral albedo variation across each band (why band means are lossy):")
    for surf in ("snow", "conifer", "dry_grass"):
        v = np.asarray(alb[surf], float)
        parts = []
        for band in BANDS:
            lo, hi = alb["meta"]["band_ranges_nm"][band]
            s = v[(awl >= lo) & (awl <= hi)]
            parts.append(f"{band} {100*(s.max()-s.min())/s.mean():5.1f}%")
        print(f"          {surf:10s} " + "  ".join(parts))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    print("OCO-2 Spectra Explorer -- reference validation")
    rep = Report(args.quiet)
    validate_data(rep)
    regression(rep)
    scattering_checks(rep)
    sorting_checks(rep)
    l2_match(rep)
    write_selftest(rep)
    physics_summary()

    print()
    if rep.failures:
        print(f"FAILED: {len(rep.failures)} check(s)")
        for f in rep.failures:
            print(f"  - {f}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
