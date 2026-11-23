"""Positive bottom-tagged norm and a triangular down-sector meson blind spot."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def rotation(i, j, sine, phase=0.0):
    result = np.eye(3, dtype=complex)
    result[i, i] = result[j, j] = np.sqrt(1 - sine**2)
    result[i, j] = sine * np.exp(-1j * phase)
    result[j, i] = -result[i, j].conjugate()
    return result


def secular_floor(y, p):
    squares = y**2
    left, right = 0.0, max(squares)
    for _ in range(100):
        mid = (left + right) / 2
        if 2 * np.dot(p, squares / (squares + mid)) > 1:
            left = mid
        else:
            right = mid
    k = (left + right) / 2
    assert abs(2 * np.dot(p, squares / (squares + k)) - 1) < 1e-13
    loss = np.dot(p, (squares - k)**2 / (squares + k))
    return float(k), float(loss)


def tagged_strength(s):
    return float(np.sum(np.abs(s[2, :2])**2) + np.sum(np.abs(s[:2, 2])**2))


def blind_example(a, order):
    # Orthogonalize the row functionals in the specified order. Restore named d,s,b.
    perm = np.eye(3)[order, :]
    q, _ = np.linalg.qr((perm @ a).conj().T)
    right = q @ perm
    np.testing.assert_allclose(right @ right.conj().T, np.eye(3), atol=2e-15)
    s = a @ right
    ordered = perm @ s @ perm.T
    assert np.linalg.norm(np.triu(ordered, 1)) < 1e-14
    mh = 1000.0
    out = {"bottom_tag_squared": tagged_strength(s)}
    for ma in (mh, 1100.0):
        rates = []
        wilsons = []
        for i, j in ((0, 1), (0, 2), (1, 2)):
            x, z = s[i, j], s[j, i].conjugate()
            ih, ia = mh**-2, ma**-2
            coeff = [x*x*(ih-ia)/4, z*z*(ih-ia)/4, x*z*(ih+ia)/2]
            # One real color-singlet pole; count both charge-conjugate flavor decays.
            rate = 3*(abs(x)**2 + abs(z)**2)/(16*np.pi)
            rates.append(float(rate))
            wilsons.append(float(np.linalg.norm(coeff)))
            if ma != mh:
                np.testing.assert_allclose(max(map(abs, coeff)),
                    (4*np.pi/3)*rate*abs(ih-ia), atol=1e-22)
        if ma == mh:
            assert max(wilsons) < 1e-20
            out["single_real_pole_width_over_mass_K_Bd_Bs_flavors"] = rates
            out["degenerate_tree_abs_Cnorm_GeV_minus2_numerical_residual"] = max(wilsons)
        else:
            out["split_1000_1100_abs_Cnorm_GeV_minus2_K_Bd_Bs"] = wilsons
    return out


def calculate():
    ref = json.loads((HERE.parent / "1386" / "inputs.json").read_text(encoding="utf-8"))
    angles = np.asarray(ref["quark_sines_and_phase"])
    ydata = np.asarray(ref["yukawa"]["up"])
    s12, s13, s23, phase = angles[:, 0]
    v = rotation(1, 2, s23) @ rotation(0, 2, s13, phase) @ rotation(0, 1, s12)
    y = ydata[:, 0]
    p = abs(v[:, 2])**2
    k, floor = secular_floor(y, p)
    a = v.conj().T @ np.diag(y)
    w = a.conj().T[:, 2]
    form = np.diag(y**2) - 2*np.outer(w, w.conjugate())
    values, axes = np.linalg.eigh(form)
    np.testing.assert_allclose(values[0], -k, rtol=2e-12)
    right = axes[:, [2, 1, 0]]
    np.testing.assert_allclose(tagged_strength(a @ right), floor, rtol=2e-12)
    # An entirely independent uniform (deliberately conservative) two-width bound.
    low, high = angles[:, 0] - 2*angles[:, 1], angles[:, 0] + 2*angles[:, 1]
    tlo, chi, thi = y[2]-2*ydata[2, 1], y[1]+2*ydata[1, 1], y[2]+2*ydata[2, 1]
    pc_low = low[2]**2*(1-high[1]**2)
    pt_low = (1-high[2]**2)*(1-high[1]**2)
    box = pc_low*pt_low*(tlo**2-chi**2)**2/(2*thi**2)
    cases = {name: blind_example(a, order)
             for name, order in (("light_axes_first", [0,1,2]), ("bottom_axis_first", [2,1,0]))}
    assert all(case["bottom_tag_squared"] >= floor-1e-12 for case in cases.values())
    return {"central_secular_k": k, "central_positive_sum": floor,
            "two_width_variance_floor": float(box),
            "tree_blind_with_positive_partial_width": cases,
            "split_to_elastic_same_chirality_coefficient_ratio": abs(1000**2-1100**2)/(2*(1000**2+1100**2))}


def verify(a,b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key, value in a.items():
            verify(value,b[key])
    else:
        np.testing.assert_allclose(a,b, rtol=1e-9, atol=1e-19)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    path = HERE / "independent_results.json"
    data = calculate()
    if args.write_new:
        with path.open("x", encoding="utf-8") as stream:
            json.dump(data, stream, indent=2)
            stream.write("\n")
    else:
        verify(data, json.loads(path.read_text(encoding="utf-8")))
    print("PASS 1549 positive spectral sum, triangular charged-current completion, and both meson-matching branches")
