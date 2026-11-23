"""Independent signed-mass, conjugation, and finite-momentum seesaw checks."""
import argparse
import json
from math import cos, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8"))


def flavour_map(x, y, z, axis, alpha, direction=1):
    cx, sx, cy, sy = sqrt(1-x), sqrt(x), sqrt(1-y), direction*sqrt(y)
    cz, sz = sqrt(1-z), sqrt(z)
    right = np.array([[cx, sx, 0], [-sx, cx, 0], [0, 0, 1]])
    left = np.array([[cy, 0, sy], [-sz*sy, cz, sz*cy], [-cz*sy, -sz, cz*cy]])
    phases = np.ones(3, dtype=complex)
    phases[axis] = np.exp(-1j*alpha)
    return (left*phases) @ right


def signed_response(m, heavy, sign, q):
    """Mass numerator of the spacelike pole sum, including its active residues."""
    block = np.array([[0., sqrt(m*heavy)], [sqrt(m*heavy), sign*(heavy-m)]])
    poles, vectors = np.linalg.eigh(block)
    assert np.allclose(poles, sorted([-sign*m, sign*heavy]), rtol=2e-12, atol=1e-15)
    residues = vectors[0]**2
    direct = np.sum(residues*poles / (1+(poles/q)**2))
    closed = -sign*m*(heavy/(heavy+m))*(1/(1+(m/q)**2)-1/(1+(heavy/q)**2))
    assert abs(direct-closed) < max(1e-16, 2e-13*m)
    assert abs(np.dot(residues, poles)) < 2e-13*m
    return float(direct)


def verify():
    cfg = read("inputs.json")
    data = read(cfg["oscillation"])
    alpha, error = read(cfg["quark_phase"])["quark_sines_and_phase"][3]
    z = cfg["representative_s23sq"]
    output = {"spectrum": {}, "inverted_conjugation": {}}
    for name in ("NO", "IO"):
        x, y, solar, atm = [data[name][k][0] for k in ("s12", "s13", "solar", "atm_abs")]
        masses = np.array([0., sqrt(solar), sqrt(atm)] if name == "NO"
                          else [sqrt(atm-solar), sqrt(atm), 0.])
        idx = np.flatnonzero(masses)
        for parity in (1, -1):
            signed = masses.copy()
            signed[idx[-1]] *= parity
            for rev in (1, -1):
                matrices = []
                for axis in range(2):
                    u = flavour_map(x, y, z, axis, alpha, rev)
                    matrix = (u*signed) @ u.T
                    matrices.append(matrix)
                    assert np.allclose(np.linalg.svd(matrix, compute_uv=False), sorted(masses, reverse=True))
                    val = 1000*abs(matrix[0, 0])
                    key = f"{name}_axis{axis+1}_parity{parity}"
                    if rev == 1:
                        output["spectrum"][key] = float(val)
                    else:
                        assert np.isclose(val, output["spectrum"][key])
                if name == "IO":
                    first, second = matrices
                    assert np.allclose(second, np.exp(-2j*alpha)*first.conj(), atol=5e-17)
                    def cp(m):
                        return float((m[0,0]*m[1,1]*m[0,1].conjugate()**2).imag)
                    assert abs(cp(first)+cp(second)) < 1e-20
                    assert abs(cp(first)) > 1e-12
                    output["inverted_conjugation"][f"parity{parity}_direction{rev}"] = {
                        "CP_odd_e_mu_axis1_and_axis2_eV4": [cp(first), cp(second)],
                        "max_abs_entry_difference_eV": float(np.max(abs(abs(first)-abs(second))))}
    no = data["NO"]
    m = [sqrt(no["solar"][0]), sqrt(no["atm_abs"][0])]
    def response(axis, hs, q):
        coeff = flavour_map(no["s12"][0], no["s13"][0], z, axis, alpha)[0,1:]
        return float(1000*abs(sum(c*c*signed_response(mass, h, 1, q)
                                  for c, mass, h in zip(coeff, m, hs))))
    large = np.array(cfg["heavy_GeV"])*1e9
    first_q, next_q = 1e8, 2e8
    target = response(0, large, first_q)
    lo, hi = 1e6, 1e9
    assert response(1, [lo]*2, first_q) < target < response(1, [hi]*2, first_q)
    for _ in range(70):
        mid = (lo+hi)/2
        if response(1, [mid]*2, first_q) < target:
            lo = mid
        else:
            hi = mid
    common = (lo+hi)/2
    rows = [[response(0, large, q), response(1, [common]*2, q)] for q in (first_q, next_q)]
    assert abs(rows[0][0]-rows[0][1]) < 1e-12 and abs(rows[1][0]-rows[1][1]) > 1.
    output["NO_full_kernel_counterexample"] = {
        "axis1_heavy_GeV": cfg["heavy_GeV"], "axis2_equal_heavy_MeV": common/1e6,
        "spacelike_Q_MeV": [100, 200], "absolute_kernel_axis1_axis2_meV_by_Q": rows,
        "max_axis2_active_probability_deficit_per_mass": max(m)/(max(m)+common)}
    output["NO_min_separation_abs_X"] = 1-max(abs(cos(2*(alpha+i*error)))
                                                for i in (-2, 2))
    return output


def equivalent(a, b):
    if isinstance(a, dict):
        return a.keys()==b.keys() and all(equivalent(a[k], b[k]) for k in a)
    return np.allclose(a,b,rtol=2e-11,atol=1e-18)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = verify()
    dst = HERE / "independent_results.json"
    if args.write:
        with dst.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert equivalent(result, read(dst.name))
        print("1426 independent signed spectra, IO conjugation and exact pole check passed")
