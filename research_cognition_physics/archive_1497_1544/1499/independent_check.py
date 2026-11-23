"""Rayleigh-deficit stability and on-shell FCNC versus scalar exchange, independent of the author code."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def read(number, filename="inputs.json"):
    return json.loads((BASE.parent / str(number) / filename).read_text(encoding="utf-8-sig"))


def kinematic_weight(m1, m2, pole, cos2phase):
    """Charge-summed width is Nc M |P c|^2 times this expression / (8 pi)."""
    a, b = m1/pole, m2/pole
    assert a+b < 1 and abs(cos2phase) <= 1
    momentum = math.sqrt((1-(a+b)**2)*(1-(a-b)**2))
    return momentum*((m1*m1+m2*m2)*(1-a*a-b*b)-4*m1*m2*a*b*cos2phase)


def scalar_grams(couplings, masses):
    y = np.array(couplings, dtype=complex)/masses
    return sum(y*y), np.vdot(y, y).real


def audit():
    inputs = read(1386)
    y = inputs["yukawa"]
    up = np.square([row[0] for row in y["up"]])
    gaps = up[2]-up[:2]
    relative_gap_values = {}
    for gap_fraction in (0., 0.001, 0.01, 0.05):
        deficit = up[2]*gap_fraction
        probabilities = deficit/(2*gaps)
        rays = np.sqrt([*probabilities, 1-sum(probabilities)])
        measured = up[2] - np.dot(rays**2, up)
        cap = deficit/(2*math.sqrt(np.prod(gaps)))
        assert abs(deficit-measured) <= 1e-15
        assert math.isclose(rays[0]*rays[1], cap)
        relative_gap_values[str(gap_fraction)] = {"abs_P_uc_cap": cap, "abs_C_D_over_F0_cap": cap**2/2}

    # The following are one-pole tree examples. Weak-scale Yukawa to masses is an explicitly
    # chosen illustration for the finite phase-space kernel, not a pole-mass extraction.
    masses = np.array([row[0] for row in y["down"]])*246.22/math.sqrt(2)
    md, ms, mb = masses
    pole = 125.
    transfer = [kinematic_weight(md, mb, pole, k)/kinematic_weight(ms, mb, pole, k)
                for k in (-1., 1.)]
    angles = read(1490, "results.json")["tree_and_mixing"]
    width_ratios = [[row["td_over_ts"][0]**2*min(transfer),
                     row["td_over_ts"][1]**2*max(transfer)] for row in angles]
    # Algebra verification by explicit spin traces in the rest frame.
    unit = np.eye(2)
    z2 = np.zeros((2, 2))
    sigma3 = np.diag([1., -1.])
    gamma0 = np.diag([1., 1., -1., -1.])
    gamma3 = np.block([[z2, sigma3], [-sigma3, z2]])
    gamma5 = np.block([[z2, unit], [unit, z2]])
    pr, pl = (np.eye(4)+gamma5)/2, (np.eye(4)-gamma5)/2
    phase_errors = []
    pz = math.sqrt((pole**2-(ms+mb)**2)*(pole**2-(ms-mb)**2))/(2*pole)
    for phi in (0., 0.73, 1.14, np.pi/2):
        coupling = np.exp(1j*phi)
        vertex = coupling*mb*pr + coupling.conjugate()*ms*pl
        adjoint = gamma0@vertex.conj().T@gamma0
        sij = math.hypot(ms, pz)*gamma0-pz*gamma3 + ms*np.eye(4)
        sjb = math.hypot(mb, pz)*gamma0+pz*gamma3 - mb*np.eye(4)
        measured = np.trace(sij@vertex@sjb@adjoint).real/pole**2*(2*pz/pole)
        expected = kinematic_weight(ms, mb, pole, math.cos(2*phi))
        phase_errors.append(abs(measured-expected))
    assert max(phase_errors) < 2e-13
    # Saturation of the phase-independent lower bound is at cos(2phi)=-1 for a single pole.
    cap = kinematic_weight(ms, mb, pole, -1)
    ratio = (8*math.pi/3)*ms*mb/(pole**3*cap)
    accepted_width = .001
    fixed_bs = .04
    coupling_mag_sq = accepted_width/((3*pole/(8*np.pi))*fixed_bs**2*cap)
    c_lr = fixed_bs**2*coupling_mag_sq*ms*mb/pole**2
    assert abs(ratio*accepted_width-c_lr) < 1e-23
    square, positive = scalar_grams([1., 1j], np.array([500., 500.]))
    assert abs(square) < 1e-14 and positive > 0
    rng = np.random.default_rng(1499)
    for _ in range(24):
        mat = rng.normal(size=(4,))+1j*rng.normal(size=(4,))
        amplitude, norm = scalar_grams(mat, rng.uniform(110., 600., 4))
        assert abs(amplitude) <= norm*(1+1e-14)
    return {
        "normalization": "P_LR basis Q=(bar i PL j)(bar i PR j), positive residue at tree matching",
        "top_y_squared_central": up.tolist(),
        "relative_Rayleigh_deficit": relative_gap_values,
        "illustrative_tree_mass_GeV_down_MZ_yukawa_times_246_22_over_sqrt2": masses.tolist(),
        "pole_GeV": pole, "finite_bd_over_bs_kinematic_factor": transfer,
        "two_tree_gamma_windows_for_charge_summed_bd_over_bs": width_ratios,
        "unit_example_only_width_GeV": accepted_width,
        "min_abs_C_LR_bs_GeV_neg2_at_example_width": ratio*accepted_width,
        "max_width_GeV_per_abs_C_LR_in_GeV_neg2": 1/ratio,
        "two_equal_mass_opposite_CP_scalar_example": {"holomorphic_sum": abs(square),
                                                       "hermitian_sum": positive},
        "largest_spin_trace_check_discrepancy": max(phase_errors)}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
