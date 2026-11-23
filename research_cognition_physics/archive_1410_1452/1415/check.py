"""Flavor-changing neutral scalar strength forced by an observed weak quartet."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def rotation(i, j, sine):
    real = np.eye(3)
    real[i, i] = real[j, j] = np.sqrt(1-sine**2)
    real[i, j] = sine
    real[j, i] = -sine
    return real


def phase(n, angle):
    return np.eye(3) + np.expm1(1j*angle)*np.outer(n, n)


def distance_bound(j, moduli):
    def growth(d):
        upper = 0.
        base = 1.
        for x in moduli:
            upper = upper*(x+d)+base*d
            base *= x
        return upper
    lo, hi = 0., 1.
    for _ in range(66):
        mid = (lo+hi)/2
        if growth(mid) < j:
            lo = mid
        else:
            hi = mid
    return lo


def make_result():
    cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    data = json.loads((HERE/cfg["flavor_inputs"]).read_text(encoding="utf-8"))
    values = np.array(data["quark_sines_and_phase"])
    x, y, z, delta = values[:, 0]
    e = np.array([1., 0., 0.])
    up_left_t = rotation(1, 2, z) @ rotation(0, 2, y)
    down_left = rotation(0, 1, x)
    basis_phase = phase(e, -delta)
    mixing = basis_phase @ up_left_t @ phase(e, delta) @ down_left
    quartet = mixing[0,1]*mixing[1,2]*np.conj(mixing[0,2]*mixing[1,1])
    signed_j = float(quartet.imag)
    assert np.allclose([abs(mixing[0,1]), abs(mixing[0,2]), abs(mixing[1,2])],
                       [x*np.sqrt(1-y*y), y, z*np.sqrt(1-y*y)])
    assert signed_j > 0 and np.isclose(signed_j, x*y*z*np.cos(np.arcsin(x))*np.cos(np.arcsin(z))*(1-y*y)*np.sin(delta))
    bounds = []
    for width in cfg["marginal_widths"]:
        low = values[:, 0] - width*values[:, 1]
        high = values[:, 0] + width*values[:, 1]
        if width == 0:
            j = signed_j
            mods = np.abs(mixing[[0, 1, 0, 1],[1, 2, 2, 1]])
        else:
            j = np.prod(low[:3])*np.sqrt((1-high[0]**2)*(1-high[2]**2))*(1-high[1]**2)*min(np.sin(low[3]), np.sin(high[3]))
            mods = [high[0], high[2], high[1], 1.]
        d = distance_bound(j, mods)
        eps = d*d/2
        total = 8*eps*(1-eps)
        bounds.append({"width": width, "J_lower": float(j), "distance_lower": float(d),
                       "each_sector_inclusive_strength_lower": float(total),
                       "some_pair_rms_lower": float(np.sqrt(total/3))})
    beta = cfg["benchmark"]["tan_beta"]
    pref = beta+1/beta
    potentials_rotation = rotation(0, 1, .6) @ rotation(1, 2, .3)
    answers = {}
    determinant = 1.+0j
    recovered = []
    for sector, n, left, sign in [("up", up_left_t@e, up_left_t.T, -1), ("down", down_left.T@e, down_left, 1)]:
        masses = np.array(data["yukawa"][sector])[:,0]
        raw = phase(e, sign*delta/2) @ left @ np.diag(masses)
        determinant *= np.linalg.det(raw)
        unitary, singular, _ = np.linalg.svd(raw)
        recovered.append(unitary[:, ::-1])
        assert np.allclose(singular[::-1], masses)
        pairs = [(0,1), (0,2), (1,2)]
        cs = pref*np.array([n[i]*n[j] for i, j in pairs])
        coefficients = np.outer(potentials_rotation[:,1] + 1j*sign*potentials_rotation[:,2], cs)
        strength = float(np.sum(abs(coefficients)**2))
        expected = pref**2*(1-np.sum(n**4))
        assert np.isclose(strength, expected) and strength >= bounds[0]["each_sector_inclusive_strength_lower"]
        answers[sector] = {"n": n.tolist(), "pair_rms_12_13_23": np.linalg.norm(coefficients,axis=0).tolist(),
                           "inclusive_strength": strength}
    rec = recovered[0].conj().T @ recovered[1]
    assert np.allclose(abs(rec), abs(mixing)) and abs(np.angle(determinant)) < 1e-11
    assert np.isclose((rec[0,1]*rec[1,2]*np.conj(rec[0,2]*rec[1,1])).imag, signed_j)
    assert max(answers["down"]["pair_rms_12_13_23"][1:]) == 0.
    b = cfg["benchmark"]
    top, higgs, vev, width = [b[k] for k in ("m_top_GeV","m_h_GeV","v_GeV","top_width_GeV")]
    norm_ut = pref*np.prod(np.array(answers["up"]["n"])[[0,2]])
    branching = top**3/(32*np.pi*vev**2*width)*(1-(higgs/top)**2)**2*norm_ut**2
    return {"round": 1415, "signed_J_input": signed_j, "quartet_bounds": bounds, "real_seed_witness": answers,
            "vacuum_theta_rad": float(delta/2), "tree_mass_det_phase_rad": float(np.angle(determinant)),
            "witness_t_to_uh_branching_ceiling_LO_massless_u": float(branching),
            "higgs_and_heavy_spectrum_not_fitted": True, "new_independent_SM_parameter_predictions": 0}


def close(actual, stored):
    if isinstance(actual, dict):
        return actual.keys()==stored.keys() and all(close(actual[k], stored[k]) for k in actual)
    if isinstance(actual, list):
        return len(actual)==len(stored) and all(close(x,y) for x,y in zip(actual,stored))
    if isinstance(actual, bool):
        return actual is stored
    return bool(np.isclose(actual, stored, rtol=2e-8, atol=5e-14))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = make_result()
    if args.check:
        assert close(result, json.loads((HERE/"results.json").read_text(encoding="utf-8")))
        print("1415 charged-current, neutral-scalar sum and real-seed witness verified")
    else:
        print(json.dumps(result, indent=2))
