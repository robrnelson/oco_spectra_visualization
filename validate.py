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
    taus = {
        # tau_co2/tau_o2 are tabulated at P0 and scale with surface pressure;
        # tau_h2o is a tropospheric profile and deliberately does NOT (see data meta)
        "co2": ps * (state["xco2"] / XCO2_REF) * d["tau_co2"],
        "o2": ps * d["tau_o2"],
        "h2o": state["h2o"] * d["tau_h2o"],
    }
    taus["ray"] = (rayleigh_tau(d["wl"]) * ps) if state["rayleigh"] else np.zeros(d["n"])
    taus["aer"] = (state["aod"] * (d["wl"] / AOD_REF_NM) ** -state["angstrom"]
                   if state["aod"] > 0 else np.zeros(d["n"]))
    return taus


def forward(band, state):
    """Line-by-line radiance and continuum on the fine grid (unconvolved)."""
    d = band_data(band)
    mu0 = math.cos(math.radians(state["sza"]))
    m = 1.0 / mu0 + 1.0 / math.cos(math.radians(state["vza"]))
    taus = optical_depths(band, d, state)
    tau_tot = taus["co2"] + taus["o2"] + taus["h2o"] + taus["ray"] + taus["aer"]
    a = albedo_spectrum(band, d, state)
    Lc = d["e0"] * mu0 * a / math.pi
    L = Lc * np.exp(-m * tau_tot)
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
                ilsStretch=1.0, shiftNm=0.0)


def l2_state():
    """State matching the L2 retrieval of the bundled sounding."""
    M = load("sounding_oco2.json")["meta"]
    s = default_state()
    s.update(xco2=M["xco2_l2_ppm"], sza=M["sza"], vza=M["vza"],
             psurf=M["surface_pressure_hpa"], surface="l2")
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
    st = l2_state()
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


def l2_match(rep):
    """The widget's "load L2 state" must actually fit the sounding.

    An earlier version also switched on Rayleigh and the L2 aerosol optical depth,
    and applied one band's albedo slope to all three. Because aerosol here is pure
    extinction while the L2 albedo was retrieved with full multiple scattering,
    that double-counted the loss and drove the A-band 3x worse.
    """
    print("\n[3] \"load L2 state\" fit quality")
    st = l2_state()                                    # what the button now sets
    good = sounding_residuals(st)
    for band in BANDS:
        g = good[band]
        print(f"        {band:6s} rms {g['rms']:6.3f}  ({100*g['rms']/g['cont']:4.2f} % of "
              f"continuum, bias {g['bias']:+.3f})")
    rep.check(good["aband"]["rms"] < 4.0,
              f"A-band rms {good['aband']['rms']:.3f} < 4.0 at the L2 state")
    rep.check(abs(good["aband"]["bias"]) < 1.5,
              f"A-band bias {good['aband']['bias']:+.3f} within +/-1.5 (not systematically dark)")

    # show why aerosol/Rayleigh are deliberately left off
    M = load("sounding_oco2.json")["meta"]
    bad = sounding_residuals(dict(st, rayleigh=True, aod=M["aod_total_l2"], angstrom=1.0))
    print(f"        if aerosol+Rayleigh were loaded too: A-band rms would be "
          f"{bad['aband']['rms']:.3f} (bias {bad['aband']['bias']:+.3f}) — "
          f"{bad['aband']['rms']/good['aband']['rms']:.1f}x worse, hence excluded")
    rep.check(bad["aband"]["rms"] > 1.5 * good["aband"]["rms"],
              "extinction-only aerosol demonstrably degrades the A-band fit")


def write_selftest(rep):
    """Sample the model over several states so the JS can be checked against it."""
    print("\n[4] selftest reference")
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
    print("\n[5] magnitudes of the added physics (airmass 2.10, AOD 0.071, alpha 1)")
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
