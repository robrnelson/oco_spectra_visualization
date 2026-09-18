#!/usr/bin/env python3
"""
Monte Carlo validation of the two-stream scattering model used by index.html.

The widget replaces pure extinction with a delta-Eddington two-stream layer coupled
to the Lambertian surface. That is an approximation, so this script checks it
against an independent Monte Carlo solution of the same problem: a plane-parallel
homogeneous scattering layer with a Henyey-Greenstein phase function over a
Lambertian surface, illuminated by a collimated beam.

    python3 validate_scattering.py [--photons N]

It reports the system albedo from both methods. The two-stream functions are
imported from validate.py, so this validates the shipped implementation rather
than a copy of it.
"""

import argparse
import sys

import numpy as np

from validate import two_stream_beam, two_stream_diffuse


def system_albedo(tau, w, g, a, mu0):
    """What the widget computes: beam reflectance + surface path through the layer."""
    Rb, Tb = two_stream_beam(tau, w, g, mu0)
    s, Td = two_stream_diffuse(tau, w, g)
    return float(Rb + Tb * a * Td / (1 - a * s)), float(Rb), float(Tb)


def monte_carlo(tau, w, g, a, mu0, n, rng):
    """Independent reference. Returns the system albedo (upward flux / incident)."""
    y = np.zeros(n)                 # optical depth from layer top
    mu = np.full(n, float(mu0))     # direction cosine, >0 means travelling downward
    weight = np.ones(n)
    alive = np.ones(n, bool)
    up = 0.0

    for _ in range(10000):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        y_new = y[idx] + (-np.log(rng.random(idx.size))) * mu[idx]

        # escaped through the top
        esc = y_new < 0
        if esc.any():
            up += weight[idx[esc]].sum()
            alive[idx[esc]] = False

        # reached the surface
        hit = y_new > tau
        ih = idx[hit]
        if ih.size:
            if a <= 0:
                alive[ih] = False
            else:
                weight[ih] *= a
                y[ih] = tau
                mu[ih] = -np.sqrt(rng.random(ih.size))   # cosine-weighted, upward
                dead = weight[ih] < 1e-4                 # russian-roulette-free cutoff
                alive[ih[dead]] = False

        # scattering or absorption inside the layer
        ins = idx[~esc & ~hit]
        if ins.size:
            y[ins] = y_new[~esc & ~hit]
            surv = rng.random(ins.size) < w
            alive[ins[~surv]] = False
            iv = ins[surv]
            if iv.size:
                u = rng.random(iv.size)
                if abs(g) < 1e-6:
                    ct = 2 * u - 1
                else:
                    ct = (1 + g * g - ((1 - g * g) / (1 + g - 2 * g * u)) ** 2) / (2 * g)
                ct = np.clip(ct, -1, 1)
                phi = 2 * np.pi * rng.random(iv.size)
                m0 = mu[iv]
                st = np.sqrt(np.maximum(0.0, 1 - ct * ct))
                sm = np.sqrt(np.maximum(0.0, 1 - m0 * m0))
                mu[iv] = np.clip(m0 * ct + sm * st * np.cos(phi), -1, 1)
    return up / n


CASES = [
    # tau,   w0,   g,     a,     mu0     note
    (0.051, 0.95, 0.65, 0.353, 0.907, "the bundled sounding's A-band aerosol"),
    (0.051, 1.00, 0.65, 0.353, 0.907, "conservative scattering"),
    (0.200, 0.95, 0.65, 0.353, 0.907, ""),
    (0.500, 0.95, 0.65, 0.353, 0.907, ""),
    (0.500, 0.90, 0.70, 0.100, 0.707, "dark surface, 45 deg sun"),
    (0.200, 0.80, 0.60, 0.600, 1.000, "absorbing aerosol, bright surface"),
    (1.000, 0.95, 0.65, 0.300, 0.866, ""),
    (0.100, 0.99, 0.75, 0.950, 0.500, "snow, low sun"),
    (0.300, 0.92, 0.70, 0.020, 0.940, "near-black surface (ocean)"),
    (1.500, 0.99, 0.80, 0.150, 0.800, "thick high aerosol"),
    (0.800, 0.50, 0.10, 0.400, 0.500, "strongly absorbing, near-isotropic"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--photons", type=int, default=400000)
    args = ap.parse_args()
    rng = np.random.default_rng(20240801)

    print("Two-stream (delta-Eddington, Meador & Weaver beam) + Lambertian surface")
    print(f"vs Monte Carlo, {args.photons} photons per case\n")
    print(f"{'tau':>6} {'w0':>5} {'g':>5} {'alb':>5} {'mu0':>5} | "
          f"{'A_2str':>7} {'A_MC':>7} {'diff':>8} {'rel%':>6}  note")
    print("-" * 78)

    worst_abs = worst_rel = 0.0
    for tau, w, g, a, mu0, note in CASES:
        A, _, _ = system_albedo(tau, w, g, a, mu0)
        Amc = monte_carlo(tau, w, g, a, mu0, args.photons, rng)
        d = A - Amc
        rel = 100 * d / Amc if Amc > 0 else 0.0
        worst_abs = max(worst_abs, abs(d))
        worst_rel = max(worst_rel, abs(rel))
        print(f"{tau:6.3f} {w:5.2f} {g:5.2f} {a:5.3f} {mu0:5.3f} | "
              f"{A:7.4f} {Amc:7.4f} {d:+8.4f} {rel:+6.1f}  {note}")

    print(f"\nworst absolute error in system albedo: {worst_abs:.4f}")
    print(f"worst relative error:                  {worst_rel:.1f} %")
    print("\nMonte Carlo statistical noise is roughly 1/sqrt(N) ~ "
          f"{100/np.sqrt(args.photons):.2f} % of the albedo, so differences below")
    print("that are not meaningful. Large *relative* errors occur only where the")
    print("albedo itself is tiny (a near-black surface), where the absolute error")
    print("stays small. This is a two-stream approximation, not a benchmark RT code.")

    ok = worst_abs < 0.02
    print("\n" + ("PASS" if ok else "FAIL")
          + f": worst absolute error {worst_abs:.4f} "
          + ("<" if ok else ">=") + " 0.02 tolerance")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
