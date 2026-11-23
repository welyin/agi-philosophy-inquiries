"""Charged spectrum and independently measurable vector couplings of general Higgs sectors."""
import argparse
import hashlib
import json
import math
import runpy
from fractions import Fraction as F
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE / "inputs.json").read_bytes())


def multiplet(j, y):
    w = F(j*(j+2)-y*y, 2)
    z = y*y
    diff = w-z
    # Cheapest field with these weights is real for Y=0 and complex otherwise.
    divisor = 1 if y else 2
    charges = [(F(m+y, 2), F(m, 2)) for m in range(-j, j+1, 2)]
    em = sum(q*q for q, _ in charges)/divisor
    pole = tuple(sum((term(q, m) for q, m in charges if abs(q) >= 2), F(0))/divisor
                 for term in (lambda q, m: m*m, lambda q, m: -2*q*m, lambda q, m: q*q))
    return {"j": j, "y": y, "w": w, "z": z, "d": diff/w,
            "v": diff**2/w, "em": em, "pole": pole, "label": [j, y]}


def minimum_quadratic(polynomial, low=F(22,100), high=F(24,100)):
    a, b, c = polynomial
    vertex = min(high, max(low, -b/(2*c))) if c else low
    return min(a+b*x+c*x*x for x in (low, vertex, high))


def frontier(states, mean_d=F(0)):
    # An optimal fraction vector at fixed sum and fixed first moment needs <=2 vertices.
    rows = []
    for i, p in enumerate(states):
        if p["d"] == mean_d:
            rows.append((p["v"],p["em"],minimum_quadratic(p["pole"]),[p["label"]]))
        for q in states[i+1:]:
            if (p["d"]-mean_d)*(q["d"]-mean_d) >= 0:
                continue
            fraction = (mean_d-q["d"])/(p["d"]-q["d"])
            v = fraction*p["v"] + (1-fraction)*q["v"]
            coeffs = [u+t for u,t in zip(p["pole"],q["pole"])]
            rows.append((v,p["em"]+q["em"],minimum_quadratic(coeffs),[p["label"],q["label"]]))
    return rows


def verify_representation_identities():
    # Use actual raising/lowering matrices separately from the polynomial proof.
    tested = 0
    for j in range(2, 28):
        t = j/2
        ms = np.arange(-t, t+1)
        raising = np.diag(np.sqrt((t-ms[:-1])*(t+ms[:-1]+1)), -1)
        lowering = raising.T
        for y in range(j%2, j+1, 2):
            v = np.eye(j+1)[:, (j-y)//2]
            g2 = 2*(np.linalg.norm(raising@raising@v)**2+np.linalg.norm(lowering@lowering@v)**2)
            p = multiplet(j, y)
            w, z = p["w"], p["z"]
            assert abs(g2-float(w*w-4*w+3*z)) < 1e-8
            assert w-1 >= 2+p["d"]
            assert 3*w*(w-16)+13*(w-z)**2-16*(w-z) >= 0
            tested += 1
    return tested


def max_mass_from_width(weight, width_cap, prefactor, mz):
    fraction = width_cap/(weight*prefactor)
    if fraction >= 1:
        return None
    # Required lower bound on the maximum physical mass in the contributing charge sectors.
    return mz/2*math.sqrt(1-fraction**(2/3))


def run():
    states = [multiplet(j,y) for j in range(2,11) for y in range(j%2,j+1,2) if j+y<=10]
    rows = frontier(states)
    # For fixed maximum charge all possible non-doublet representations are finite.
    min_mismatch_for_charge = {}
    for charge_cap, expected in enumerate((None, F(3), F(3), F(24,17), F(0)), start=1):
        allowed = [p for p in states if p["j"] + p["y"] <= 2 * charge_cap]
        candidates = frontier(allowed)
        best = min((row[0] for row in candidates), default=None)
        assert best == expected
        min_mismatch_for_charge[str(charge_cap)] = str(best) if best is not None else None
    limits = [(None, F(6), F(2704,10000)), (F(3), F(36), F(32508,10000)),
              (F(24,17), F(56), F(46304,10000))]
    # j+y>10 implies q_max >=6 and a consecutive ladder containing 2..q_max.
    tail_em = sum(i*i for i in range(1,7))
    tail_pole = F(76,100)**2*F(5*(25-1),12)
    assert tail_em > 56 and tail_pole > F(46304,10000)
    input_pole = runpy.run_path(str(HERE/DATA["pole_reference"]))
    mz = input_pole["MZ"]-2*input_pole["SIGMZ"]
    gf = input_pole["GF"]-2*input_pole["SIGGF"]
    tree_factor = input_pole["scale"](mz,gf)
    quark_floor = tree_factor*input_pole["extrema"]((3,0))[0]*input_pole["MATCH"][0]
    leptons_floor = sum(a-2*b for a,b in DATA["partial_width_leptons_MeV"])/1000
    total_cap = DATA["Z_width_GeV"][0]+2*DATA["Z_width_GeV"][1]
    room = total_cap-quark_floor-leptons_floor
    # One complex scalar-pair coefficient in the same normalization as the FIVE-quark width.
    prefactor = gf*mz**3/(6*math.pi*math.sqrt(2))
    beta = math.sqrt(1-4*DATA["examples"]["common_mass_cap_over_cm_energy"]**2)
    tables = []
    for vcap, em, pole in limits:
        admitted = [row for row in rows if vcap is None or row[0]<vcap]
        assert min(row[1] for row in admitted) == em and min(row[2] for row in admitted) == pole
        witness = next(row for row in admitted if row[1]==em and row[2]==pole)
        tables.append({"mismatch_per_fraction_strict_upper": str(vcap), "minimum_photon_weight": int(em),
                       "minimum_q_abs_ge2_Z_weight": float(pole), "attained_by_j_y": witness[3],
                       "photon_R_floor_for_given_mass_cap":float(em)/4*beta**3,
                       "massless_multicharge_Z_width_floor_GeV":prefactor*float(pole),
                       "necessary_max_mass_lower_GeV":max_mass_from_width(float(pole),room,prefactor,mz)})
    rho = DATA["rho0_fit"]
    dlo, dhi = (1-1/(rho["central"]+i*rho["sigma"]) for i in (-2,2))
    f0 = DATA["examples"]["conditional_non_doublet_fraction_lower"]
    mean_range = (F(str(dlo/f0)),F(str(dhi/f0)))
    # All d for the finite list outside this interval, except the septet at d=0.
    assert all(p["d"]==0 or p["d"] < mean_range[0] or p["d"] > mean_range[1] for p in states)
    robust = {}
    for budget in (36,56):
        values = [row[0] for d in (*mean_range,F(0)) for row in frontier(states,d) if row[1]<budget]
        robust[str(budget)] = float(min(values))
    s = DATA["examples"]["upper_double_charge_WW_strength"]
    v = DATA["examples"]["upper_W_minus_Z_neutral_strength"]
    return {"representation_matrix_checks":verify_representation_identities(),"rho_one_frontier":tables,
            "minimum_mismatch_over_f_for_charge_cap":min_mismatch_for_charge,
            "rho_one_example_fraction_max":min(1,s/2,(s+13*v/3)/15),
            "total_width_upper_GeV":total_cap,"ordinary_quark_width_floor_GeV":quark_floor,
            "ordinary_charged_leptons_width_floor_GeV":leptons_floor,"remaining_width_GeV":room,
            "rho_fit_conditional_D_range":[dlo,dhi],
            "necessary_mismatch_over_f_below_photon_budgets_at_f_ge_0p1":robust,
            "large_charge_tail_min_em":tail_em, "large_charge_tail_min_Z_weight":float(tail_pole)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(run(),indent=2,allow_nan=False)+"\n"
    path = HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==json.loads(payload)
    else:
        path.write_text(payload, encoding="utf-8")
    print("1789 electroweak charge and condensation bounds",hashlib.sha256(payload.encode()).hexdigest()[:12])
