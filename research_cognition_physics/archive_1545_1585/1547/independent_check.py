"""Bottom-tagged colour-octet widths under independent left and right identifications."""
import argparse
import json
from itertools import product
from math import pi, sqrt
from pathlib import Path
import numpy as np


def bottom_counts(y):
    bb = abs(y[-1, -1])**2
    mixed = np.sum(abs(y[:-1, -1])**2) + np.sum(abs(y[-1, :-1])**2)
    return float(2*bb+mixed), float(mixed)


def exact_limits(b, mu, tau):
    tag = 2*min(abs(b-mu), abs(b-tau))**2
    mixed_if_bb_zero = 2*(b-mu)*(tau-b)
    assert mu < b < tau and abs(b-mu) < b-1e-5
    return np.array([tag, mixed_if_bb_zero])


def spectra():
    here = Path(__file__).resolve().parent
    source = json.loads((here.parent/"1402"/"inputs.json").read_text(encoding="utf-8"))
    rows = np.asarray(source["scales"]["1e16"])
    chosen = rows[[2, 4, 5]]
    (b, mu, tau) = chosen[:, 0]
    central = exact_limits(b, mu, tau)
    corners = [exact_limits(*(chosen[:, 0]+2*chosen[:, 1]*sign))
               for sign in product((-1, 1), repeat=3)]
    low = np.min(corners, axis=0)
    # Both extrema for the gap-product on the b interval occur at endpoints
    # when minimizing: it is concave in b and monotone in each bracketing lepton.
    weight = (b-mu)/(tau-mu)
    rot = np.eye(3)
    rot[1:, 1:] = [[sqrt(weight), -sqrt(1-weight)],
                   [sqrt(1-weight), sqrt(weight)]]
    d = np.diag(rows[:3, 0])
    e = np.diag(rows[3:, 0])
    y = d-rot@e@rot.T
    assert abs(y[-1, -1]) < 1.e-16
    assert np.isclose(bottom_counts(y)[1], central[1], atol=1.e-17)
    no_mix = bottom_counts(d-e)
    assert abs(no_mix[0]-central[0]) < 1e-17 and no_mix[1] == 0.
    rng = np.random.default_rng(1547)
    for _ in range(250):
        matrices = [np.linalg.qr(rng.normal(size=(3, 3))+1j*rng.normal(size=(3, 3)))[0]
                    for _ in (0, 1)]
        delta = d-matrices[0]@e@matrices[1].conj().T
        assert bottom_counts(delta)[0] >= central[0]-1.e-15
        aa = abs(delta[-1, -1])
        finite_bound = central[1]-2*aa*abs(mu+tau-2*b)-2*aa**2
        assert bottom_counts(delta)[1] >= finite_bound-1.e-15
    c = 3/(32*pi)
    return {"input_y_b_mu_tau": [b, mu, tau],
            "central_delta_squared_tag_and_mixed_if_bb_zero": central.tolist(),
            "two_width_minimum_delta_squared": low.tolist(),
            "central_width_floors_if_adjoint_exhausts_EW_VEV": (c*central).tolist(),
            "two_width_minimum_width_floors_if_adjoint_exhausts_EW_VEV": (c*low).tolist(),
            "zero_bb_witness_tau_overlap_squared": weight,
            "zero_bb_witness_D_minus_E": y.tolist()}


def partial_spectrum():
    # One accessible real mass eigenstate from two complex neutral components.
    # A unit vector of a real four-dimensional neutral mixing matrix
    # is represented by a complex two-component row of norm 1.
    measured = np.array([[1, 1j]])/sqrt(2)
    visible = measured.conj().T@measured
    witness = np.array([[1], [1j]])/sqrt(2)
    assert np.allclose(measured@witness, 0)
    c_aligned, c_blind = measured, witness.conj().T
    assert np.allclose(c_aligned@np.linalg.pinv(visible)@c_aligned.conj().T, 1)
    assert np.allclose(c_blind@witness, 1)
    return {"aligned_single_pole_norm_squared_gain": 1,
            "independent_mass_direction_with_one_pole": 1,
            "same_pole_rate_on_mass_carrying_null_vector": 0}


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    out = {"physical_b": spectra(), "calibrated_subset": partial_spectrum()}
    target = Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with target.open("x", encoding="utf-8") as f:
            json.dump(out, f, indent=2)
            f.write("\n")
    else:
        assert out == json.loads(target.read_text(encoding="utf-8"))
        print("PASS 1547 tagged widths, sharp compatible alignments and partial spectrum")


if __name__ == "__main__":
    main()
