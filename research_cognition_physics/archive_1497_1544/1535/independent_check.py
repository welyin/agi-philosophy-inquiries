"""Spin-trace and three-flavour checks independent of the author rate code."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
cfg = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
cal = json.loads((HERE/cfg["mass_calibration"]).read_text(encoding="utf-8"))
YUK = np.array(cal["one_TeV_DR"]["up"])[:,0]
T = cfg["kinetic_top_mass_GeV"][1]
Q = [0, cfg["kinetic_charm_mass_GeV"][1]]
H = cfg["mH_GeV"]
TOTAL = cfg["total_width_GeV"]["reference"]
RUN1_BR = 0.0079
SINGLE_AXIS_ATLAS = [0.00026, 0.00034]


def spin_trace(mass, right, left_star):
    ident = np.eye(4)
    zero = np.zeros((2,2))
    pauli1 = np.array([[0,1],[1,0]])
    pauli2 = np.array([[0,-1j],[1j,0]])
    pauli3 = np.diag([1,-1])
    g0 = np.diag([1,1,-1,-1])
    g1,g2,g3 = (np.block([[zero,s],[-s,zero]]) for s in [pauli1,pauli2,pauli3])
    g5 = 1j*g0@g1@g2@g3
    amp = (ident+g5)*right/2 + (ident-g5)*left_star/2
    energy = (T*T+mass*mass-H*H)/(2*T)
    p = math.sqrt(energy**2-mass**2)
    final = energy*g0 - p*g3 + mass*ident
    initial = T*(g0+ident)
    trace = np.trace(final @ amp @ initial @ g0 @ amp.conj().T @ g0)/2
    assert abs(trace.imag) < 1e-8
    return float(trace.real * p/(8*math.pi*T*T))


def example(x, real_shifts=None):
    odd = np.array([-x, x-1, 1.])
    shifts = np.zeros(3) if real_shifts is None else np.array(real_shifts)
    diag = YUK*(shifts+1j*odd)/math.sqrt(2)
    left = np.sqrt(abs(diag))
    right = np.divide(diag, left, out=np.zeros(3,dtype=complex), where=left>0)
    z = np.outer(left,right)
    assert np.linalg.matrix_rank(z,tol=1e-12) == 1
    for r in (.8,1.,1.3):
        m = np.diag(YUK/math.sqrt(2))
        ref = np.linalg.det(r*m)
        val = np.linalg.det(r*(m+(r*r-1)*z/2))/ref
        assert abs(val.imag)<1e-10
    return [spin_trace(Q[i], z[i,2], z[2,i].conjugate()) for i in range(2)]


def calc():
    ku, kc = example(1)[0],example(0)[1]
    k = np.array([ku,kc])
    for x in (0,.22,.5,.999,1):
        np.testing.assert_allclose(example(x), k*[x,1-x], atol=1e-13)
    for x in (-.24,0,.3,.87,1,1.54):
        signed_weights = np.array([x,1-x])
        floors = np.array([spin_trace(Q[i], math.sqrt(YUK[i]*YUK[2]/2),
                                 math.copysign(math.sqrt(YUK[i]*YUK[2]/2), signed_weights[i]))
                           for i in range(2)]) * abs(signed_weights)
        for shifts in ([0,0,0],[-4,1,3],[.3,.7,-2]):
            assert np.all(np.array(example(x,shifts))+1e-14 >= floors)
    branching = k/TOTAL
    hypothetical_axes = np.array(SINGLE_AXIS_ATLAS)/branching
    return {"individual_floor_GeV_at_unit_slope": k.tolist(),
            "individual_floor_branching_at_unit_slope": branching.tolist(),
            "unit_slope_80percent_charm_compensation_GeV": example(.2),
            "run1_flavour_inclusive_limit_on_abs_top_odd_if_template_unchanged":
                math.sqrt(RUN1_BR/min(branching)),
            "separate_ATLAS_2024_single_channel_axis_on_top_odd": np.sqrt(hypothetical_axes).tolist(),
            "hypothetical_rectangular_joint_sensitivity": math.sqrt(sum(hypothetical_axes)),
            "adversarial_u_share_in_hypothetical_joint_rectangle":
                hypothetical_axes[0]/sum(hypothetical_axes)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    flag = parser.add_mutually_exclusive_group(required=True)
    flag.add_argument("--check",action="store_true")
    flag.add_argument("--write-new",action="store_true")
    opts = parser.parse_args()
    out = HERE/"independent_results.json"
    obj = calc()
    if opts.write_new:
        with out.open("x",encoding="utf-8") as f:
            json.dump(obj,f,indent=2)
            f.write("\n")
    else:
        old=json.loads(out.read_text(encoding="utf-8"))
        assert obj.keys()==old.keys()
        for key in obj:
            np.testing.assert_allclose(obj[key],old[key],rtol=1e-12,atol=1e-14)
    print("1535 independent OK",obj)
