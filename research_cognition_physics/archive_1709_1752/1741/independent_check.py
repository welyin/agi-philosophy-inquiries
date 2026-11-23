"""Independent scalar/fermion spectra and vacuum-orientation response at the 1739 point."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
NC = 3
LOOP = NC / (8 * math.pi**2)


def spectrum(y, loss, p, v2, soft_left, soft_right, phase=.91):
    """Masses in pole units. Soft terms are gauge invariant in the original row basis."""
    lines = np.zeros((2, 3), dtype=complex)
    a, b = math.sqrt(1-p), math.sqrt(p)*np.exp(1j*phase)
    lines[:, [0, 2]] = [[a, b], [-b.conjugate(), a]]
    heavy = np.column_stack((np.diag(np.sqrt(loss)) @ lines, np.diag(np.sqrt(1-loss))))
    top = np.column_stack((np.diag(y) * math.sqrt(v2), np.zeros((3, 2))))
    matrix = np.vstack((top, heavy))
    fermions = np.linalg.svd(matrix, compute_uv=False)**2
    left = np.linalg.eigvalsh(matrix @ matrix.conj().T + np.diag(soft_left))
    right = np.linalg.eigvalsh(matrix.conj().T @ matrix + np.diag(soft_right))
    assert min(left) > 0 and min(right) > 0
    return fermions, left, right


def cw(spectra, log_scale=0.):
    def trace(vals):
        return np.sum(vals**2*(np.log(vals)-log_scale-1.5))
    fermions, left, right = spectra
    return NC/(32*math.pi**2)*(trace(left) + trace(right) - 2*trace(fermions))


def orientation_integral(y, loss, v2, soft, order):
    """Slopes and free-energy gaps. The positive integrands resolve a very narrow interval."""
    nodes, gw = np.polynomial.legendre.leggauss(order)
    splits = sorted({-84., 84., 0., math.log(soft), *(math.log(v2*yi**2) for yi in y)})
    weights = np.concatenate([gw*(b-a)/2 for a, b in zip(splits, splits[1:])])
    points = np.concatenate([(a+b)/2 + (b-a)*nodes/2 for a, b in zip(splits, splits[1:])])
    q = np.exp(points)
    diag = [q + v2*yi**2*(1-r/(1+q)) for yi, r in zip(y[[0, 2]], loss)]
    slope = q/(1+q)/(diag[0]*diag[1])
    contrast = v2 * (y[2]**2-y[0]**2) * (loss[0]-loss[1])
    h = slope*contrast
    assert np.all((h > 0) & (h < 1))
    measure = LOOP*weights*q*np.minimum(q, soft)
    c0 = measure @ slope
    width = measure @ (slope*h/(1-h))
    halfshift = measure @ (slope*.5*h/(1-.5*h))
    gaps = {p: float(measure @ (-np.log1p(-p*h))) for p in (0., .25, .5, .75, 1.)}
    return c0, width, halfshift, contrast, gaps


def matched_gauge_gap(y, loss, v2, soft, h1, reference_gap, order):
    """Left heavy partner masses may split, with their spectral product held fixed."""
    hbar = 1 + soft
    h2 = hbar**2/h1
    nodes, weights = np.polynomial.legendre.leggauss(order)
    cuts = sorted({-82., 82., 0., math.log(soft), math.log(hbar), math.log(h1), math.log(h2)})
    zd, zb = v2*y[[0, 2]]**2
    total_change = 0.
    for lo, hi in zip(cuts, cuts[1:]):
        q = np.exp((hi+lo)/2 + (hi-lo)*nodes/2)
        common = q+soft
        def det(z, h, r):
            return (common+z)*(q+h)-z*r
        oriented = np.zeros_like(q)
        for sign, h, r in zip((1,-1), (h1,h2), loss):
            cross = (h-hbar)*r*(zb-zd)*common
            denominator = det(zd,h,r)*det(zb,hbar,r)
            oriented += sign*np.log1p(cross/denominator)
        total_change += (hi-lo)/2 * (weights @ (q*q*oriented))
    return reference_gap + LOOP/2*total_change


def run():
    old = json.loads((HERE.parent/"1739"/"higgs_results.json").read_bytes())["absolute_and_oscillations"]
    y = np.array(old["down_yukawas_at_M"])
    m, t, c = old["M_GeV"], (old["soft_GeV"]/old["M_GeV"])**2, old["c"]
    loss = np.array([1-c*c, (1-c)/2])
    v2 = (246.22/m)**2 / 2
    common = np.full(5, t)
    eps = .001
    shifts = np.array([.1, 1.])
    split_left = np.array([eps]*3 + (eps + shifts).tolist())
    split_right = np.full(5, eps)

    def divided(high, low):
        return (high*math.log(high)-low*math.log(low))/(high-low)
    coeff = np.array([divided(1+eps+s, eps)+divided(1+eps, eps) for s in shifts])
    split_small = NC/(16*math.pi**2) * (y[2]**2-y[0]**2) * (loss[0]*coeff[0]-loss[1]*coeff[1])
    assert split_small < 0
    small_bg = LOOP*((1+t)*math.log1p(t) - t*math.log(t))
    rows = []
    for multiplier in (1., 100., 1e6):
        bg = v2*multiplier
        c0, wid, hs, contrast, gaps = orientation_integral(y, loss, bg, t, 72)
        low_order = orientation_integral(y, loss, bg, t, 48)
        np.testing.assert_allclose([c0,wid,hs],low_order[:3],rtol=5e-9,atol=0)
        assert 0 < hs < wid and c0 < small_bg
        baseline = spectrum(y, loss, 0., bg, common, common)
        reference = cw(baseline)
        max_residual = 0.
        for p in (.25, .5, .75, 1.):
            current = spectrum(y, loss, p, bg, common, common)
            observed = cw(current) - reference
            max_residual = max(max_residual, abs(observed/gaps[p]-1))
            assert abs(observed/gaps[p]-1) < 5e-4
            assert abs(cw(current,1.)-cw(baseline,1.)-gaps[p]) < max(1e-11,2e-7*gaps[p])
        split_baseline = spectrum(y, loss, 0., bg, split_left, split_right)
        split_swap = spectrum(y, loss, 1., bg, split_left, split_right)
        split_gap = cw(split_swap)-cw(split_baseline)
        assert split_gap < 0
        if multiplier == 1:
            assert abs(split_gap/(bg*split_small)-1) < 2e-4
        rows.append(dict(background_over_physical=multiplier,
                         common_CW_swap_gap_GeV4=gaps[1.]*m**4,
                         counterterm_onset=-c0, width_to_full_swap=wid,
                         distance_onset_to_half_rotation=hs,
                         weakest_gap_over_slope_estimate=gaps[1.]/(contrast*small_bg),
                         full_svd_max_relative_error=max_residual,
                         split_swap_gap_GeV4=split_gap*m**4))
    ref_gap = rows[0]["common_CW_swap_gap_GeV4"]/m**4
    threshold_rows = []
    for h1 in (1., 1+t, (1+t)**2):
        gap = matched_gauge_gap(y, loss, v2, t, h1, ref_gap, 64)*m**4
        fine = matched_gauge_gap(y, loss, v2, t, h1, ref_gap, 96)*m**4
        np.testing.assert_allclose(gap, fine, rtol=2e-10)
        h2 = (1+t)**2/h1
        part = common.copy()
        part[3:] = [h1-1,h2-1]
        aligned = spectrum(y, loss, 0., v2, part, common)
        swapped = spectrum(y, loss, 1., v2, part, common)
        observed = (cw(swapped)-cw(aligned))*m**4
        assert abs(observed/gap-1) < .014
        threshold_rows.append(dict(left_heavy_scalar_mass_over_M=[math.sqrt(h1),math.sqrt(h2)],
                                   swap_gap_GeV4=gap))
    assert 0 < threshold_rows[0]["swap_gap_GeV4"] < threshold_rows[1]["swap_gap_GeV4"] < threshold_rows[2]["swap_gap_GeV4"]
    return dict(source="1739/higgs_results.json absolute_and_oscillations",
                pole_GeV=m, higgs_background_GeV2=v2*m**2, common_soft_over_pole_sq=t,
                small_background_CW_coefficient=small_bg, rows=rows,
                alternate_common_soft_mass_squared_over_M2=eps,
                alternate_extra_left_heavy_mass_squared_over_M2=shifts.tolist(),
                alternate_gate_coefficients=coeff.tolist(),
                alternate_swap_slope_GeV2=split_small*m**2,
                fixed_one_loop_threshold_extremes_and_uniform=threshold_rows)


def compare(actual, saved):
    if isinstance(actual, dict):
        assert actual.keys() == saved.keys()
        for k in actual:
            compare(actual[k], saved[k])
    elif isinstance(actual, list):
        assert len(actual) == len(saved)
        for a,b in zip(actual,saved):
            compare(a,b)
    elif isinstance(actual, (float, np.floating)):
        np.testing.assert_allclose(actual,saved,rtol=3e-9,atol=1e-18)
    else:
        assert actual == saved


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    out = HERE/"independent_results.json"
    if args.check:
        compare(result, json.loads(out.read_bytes()))
        print("1741 independent spectrum, finite matching and soft-splitting tests passed")
    else:
        with out.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
            f.write("\n")
