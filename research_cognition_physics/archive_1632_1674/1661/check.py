"""Bound every frequency cell for single-flavon matching; test a neutral two-path link."""
import argparse
import json
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EPS = 2e-8


def read(name):
    return json.loads((HERE / name).read_bytes())


def root(R):
    B = 6 * (R + 1 / R) + 16
    return 8 / (B + np.sqrt(B**2 - 1936))


def coefficients(E, phases):
    Q = np.exp(1j * np.array([0, *phases]))[:, None] * E
    return (Q.conj().T @ [0, 1, 1], Q.conj().T @ [1, 3, 1],
            Q.conj().T @ [2, -1, 1] / np.sqrt(6))


def projector(coefs, phase, R, r, i, j):
    a, b, v = coefs
    return ((2/r)*a[i]*a[j].conj() + (11*r)*b[i]*b[j].conj()
            + 4*(phase*a[i]*b[j].conj() + phase.conj()*b[i]*a[j].conj())
            - (6/R)*(float(i == j)-v[i]*v[j].conj())) / (6*(R-1/R))


def enclose(coefs, phase, R, r, h, windows):
    a, b, v = coefs
    # h = cell half-width; each projector entry has absolute derivative below 2.
    error = 2*h + EPS
    pe = projector(coefs, phase, R, r, 0, 0).real
    pm = projector(coefs, phase, R, r, 1, 1).real
    pl, pu = np.maximum(0, pe-error), np.minimum(1, pe+error)
    assert max(pu) < .5
    p1 = abs(v[0])**2
    low = np.array([4*pl*(1-pl), 1-(p1+EPS)/(1-pu), (pm-error)/(1-pl)])
    high = np.array([4*pu*(1-pu), 1-(p1-EPS)/(1-pl), (pm+error)/(1-pu)])
    c = windows.mean(axis=1)[:, None]
    hw = (windows[:, 1:] - windows[:, :1])/2
    score = np.maximum((low-c)/hw, (c-high)/hw).max(axis=0)
    return score <= 1, float(min(score))


def stable(x):
    if isinstance(x, dict):
        return {k: stable(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [stable(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(f"{x:.10g}")
    return x


def endpoints(row):
    return row[0]-2*row[1], row[0]+2*row[2]


def solve():
    cfg = read("inputs.json")
    old = runpy.run_path(str(HERE.parent / "1659/check.py"))
    charged_cfg = read(cfg["charged_and_two_sources"])
    E, _ = old["charged_basis"](charged_cfg)
    # Eigenvector accuracy of the frozen charged input matters for the phase exclusion.
    charged = runpy.run_path(str(HERE.parent / "1657/independent_check.py"))
    point = old["read"](charged_cfg["charged"])["extended_box_GJ_test"]
    M_e = charged["mtx"](charged["texture"](point["marginal_coordinates"]), [1,1,3], point["t"]).T
    H_e = M_e @ M_e.conj().T
    lam = np.diag(E.conj().T @ H_e @ E).real
    residual = np.linalg.norm(H_e @ E - E * lam)
    assert residual / np.min(np.diff(lam)) < 1e-9
    osc = old["read"](charged_cfg["oscillation"])
    checked = read(cfg["charged_reference"])
    windows = np.array(checked["joint_windows"])
    ds = endpoints(osc["KamLAND_dm21_eV2"])
    da = endpoints(osc["T2K_dm32_NO_eV2"])
    B0, B1 = cfg["frequency_outer"]
    assert B0 < np.sqrt(1+da[0]/ds[1]) < np.sqrt(1+da[1]/ds[0]) < B1
    grid = np.linspace(B0, B1, cfg["cells"]+1)
    R = (grid[1:]+grid[:-1])/2
    h = (B1-B0)/(2*cfg["cells"])
    r = root(R)
    assert .08 < root(B1) <= root(B0) < .13
    poly_min = 6*(B0+1/B0)+16
    assert 6*root(B0)/np.sqrt(poly_min**2-1936) < .05
    cfg4 = read(cfg["flavour"])
    assert cfg4["N"] == 4 and cfg4["right"] == [2, 0, 1]
    phi = np.pi * cfg["vacuum_branch_phi_over_pi"]
    lifts = [(k, l) for k in range(12) for l in range(12)
             if k % 4 == (cfg4["right"][1]-cfg4["right"][0]) % 4
             and l % 4 == (cfg4["right"][2]-cfg4["right"][0]) % 4]
    signed = {((k+i) % 12, (l+j) % 12) for k, l in lifts for i in (0,6) for j in (0,6)}
    classes = {}
    for sign in (-1, 1):
        scores = {}
        for k in range(12):
            for l in range(12):
                coefs = coefficients(E, [k*phi, l*phi])
                ok, score = enclose(coefs, np.exp(sign*2j*np.pi/3), R, r, h, windows)
                assert not np.any(ok)
                scores[k, l] = score
        def worst(keys):
            key = min(keys, key=scores.get)
            return {"pairs": len(keys), "least_excluded_pair": list(key),
                    "minimax_lower_bound": scores[key]}
        classes[str(sign)] = {"all_integer": worst(scores), "charge_lifts": worst(lifts),
                              "charge_lifts_with_real_signs": worst(sorted(signed))}
    ex = cfg["example"]
    assert ex["port"] == 3 and (ex["k"], ex["l"]) in lifts
    z = ex["constant_amplitude"] + ex["wound_amplitude"]*np.exp(4j*phi)
    angles = [ex["k"]*phi, ex["l"]*phi + np.angle(z)]
    coefs = coefficients(E, angles)
    a, b, v = coefs
    ranges = {}
    for sign in (-1, 1):
        phase = np.exp(sign*2j*np.pi/3)
        ok, score = enclose(coefs, phase, R, r, h, windows)
        if sign < 0:
            assert not ok.any()
            ranges[str(sign)] = {"minimax_lower_bound": score}
            continue
        assert ok.any()
        x = R[ok]
        rx = r[ok]
        lower = np.maximum(ds[0], da[0]/((x+h)**2-1))
        upper = np.minimum(ds[1], da[1]/((x-h)**2-1))
        J = -(v[0]*v[1].conj()*projector(coefs, phase, x, rx, 1, 0)).imag
        elem = abs(a[0].conj()**2 + rx*phase*b[0].conj()**2)
        elem_error = .05*h*abs(b[0])**2 + EPS
        m0lo = np.sqrt(lower*(x-h)/(6*root(x-h)))
        m0hi = np.sqrt(upper*(x+h)/(6*root(x+h)))
        beta_beta = [min(m0lo*(elem-elem_error))-EPS, max(m0hi*(elem+elem_error))+EPS]
        ranges[str(sign)] = {"possible_R_outer": [min(x)-h, max(x)+h],
                             "J_outer": [min(J)-2*h-EPS, max(J)+2*h+EPS],
                             "m_beta_beta_eV_outer": beta_beta}
        assert ranges["1"]["J_outer"][1] < 0 and beta_beta[0] > 0
    d21, d32 = ex["ds"], ex["ds"]*(ex["R"]**2-1)
    assert ds[0] <= d21 <= ds[1] and da[0] <= d32 <= da[1]
    M, _, rs = old["rank_two"](d21, d32, 2*np.pi/3, "small")
    witness = old["predict"](E, M, *angles)
    D = np.diag(np.exp(1j*np.array([.3, -.2, .7])))
    physical = np.diag(np.exp(1j*np.array([0, *angles])))
    physical = physical @ M @ physical
    same_frame = old["predict"](D @ E, D.conj() @ physical @ D.conj().T, 0, 0)
    for key in witness:
        np.testing.assert_allclose(same_frame[key], witness[key], atol=2e-14)
    values = witness["A13_solar_atmospheric"]
    assert np.all((windows[:, 0] < values) & (values < windows[:, 1]))
    rx = ex["R"]
    phase = np.exp(2j*np.pi/3)
    s13 = projector(coefs, phase, rx, rs, 0, 0).real
    qmu = projector(coefs, phase, rx, rs, 1, 1).real
    np.testing.assert_allclose([4*s13*(1-s13), 1-abs(v[0])**2/(1-s13), qmu/(1-s13)], values)
    jcp = -(v[0]*v[1].conj()*projector(coefs, phase, rx, rs, 1, 0)).imag
    mbb = np.sqrt(d21*rx/(6*rs))*abs(a[0].conj()**2 + rs*phase*b[0].conj()**2)
    np.testing.assert_allclose([jcp, mbb], [witness["J_dirac"], witness["m_beta_beta_eV"]])
    large = old["large_bound"](E, osc)["A13_lower_bound"]
    assert abs(large-checked["large_root_full_frequency_box"]["A13_lower_bound"]) < EPS
    return stable({"frequency_cover": [B0, B1], "half_cell": h, "monomial_by_eta_sign": classes,
                   "large_root_A13_min_any_phases": large, "two_path_phase_deg": np.degrees(np.angle(z)),
                   "two_path_heldouts": ranges, "two_path_witness": {"dm21": d21, "dm32": d32, **witness}})


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    cli.add_argument("--write-new", action="store_true")
    args = cli.parse_args()
    out = solve()
    if args.write_new:
        with (HERE/"results.json").open("x", encoding="utf8") as f:
            f.write(json.dumps(out, indent=2) + "\n")
    else:
        assert out == read("results.json")
        print("1661 PASS: 144 exponent pairs x both signs, full frequency enclosure and two-path holdouts")
