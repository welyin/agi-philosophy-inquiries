"""Charge reconstruction from commuting gauge currents of heavy Dirac matter."""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
# Colour Cartan generator diag(1,1,-2), with its normalization absorbed into
# z. This rational convention makes every extraction of a coefficient exact.
T = [1, 1, -2]
S2, S3, S4 = [sum(t**p for t in T) for p in (2, 3, 4)]


def extract(coloured, uncoloured):
    coeff = [F(0) for _ in range(5)]  # z-power in a homogeneous quartic in x,z
    pol = F(0)
    for species, charges in [(coloured, T), (uncoloured, [0])]:
        for q, mass in species:
            for t in charges:
                pol += q**2/mass**4
                for k in range(5):
                    coeff[k] += math.comb(4,k)*q**(4-k)*t**k/mass**4
    s0 = coeff[4]/S4
    s1 = coeff[3]/(4*S3)
    s2 = coeff[2]/(6*S2)
    wc = s2 - s1/3 - 2*s0/9
    wall = coeff[0]-pol + F(4,3)*s2 + F(4,27)*s0
    assert coeff[1] == 0
    assert wc == sum((q+F(1,3))*(q-F(2,3))/m**4 for q,m in coloured)
    assert wall == (3*sum((q**2-F(1,9))*(q**2-F(4,9))/m**4 for q,m in coloured)
                    + sum(q**2*(q**2-1)/m**4 for q,m in uncoloured))
    return dict(coeff_x4_to_z4=coeff, current_curvature=pol, colour_odd_test=wc,
                all_charge_even_test=wall)


def current_curvature():
    # pi(t) = e^2/(2 pi^2) int w log(1/(1-t w/m^2)) for time-like t.
    x, weights=np.polynomial.legendre.leggauss(16)
    w=(1-x**2)/4
    return float(np.dot(weights/2, w**3)/2)


def pair_kernel(energy, mass):
    if energy <= 2*mass:
        return 0.0
    beta = math.sqrt(1-(2*mass/energy)**2)
    return (3-beta**4)*math.atanh(beta)-beta*(2-beta**2)


def finite_cut(states, energy):
    coefficients = [0.0]*5
    for q, mass in states:
        kernel = pair_kernel(energy, float(mass))
        for t in T:
            for k in range(5):
                coefficients[k] += (kernel * math.comb(4,k) *
                                    float(q**(4-k)) * t**k)
    a2 = coefficients[2]/(6*S2)
    a1 = coefficients[3]/(4*S3)
    a0 = coefficients[4]/S4
    return a2-a1/3-2*a0/9


def mismatch_cost():
    # Beyond these three labels, a +3 lepton or at least a 7/3 coloured
    # component on its own already costs more than the proposed minimum.
    ceiling=F(92,3)
    assert 3*(F(49,9)-F(1,9))*(F(49,9)-F(4,9))>ceiling
    assert 9*(9-1)>ceiling
    weights=[F(20,3),F(56,3),F(12)]
    candidates=[]
    for a in range(-4,5):
        for b in range(-1,2):
            for c in range(-2,3):
                if 5*a-7*b+3*c==1:
                    cost=sum(abs(i)*w for i,w in zip((a,b,c),weights))
                    candidates.append((cost,(a,b,c)))
    cost,labels=min(candidates)
    assert cost==ceiling and labels==(-1,0,2)
    saturated=extract([(F(-4,3),F(1)),(F(-1,3),F(1))],
                      [(F(-2),F(1)),(F(-1),F(1))]*2)
    assert saturated["all_charge_even_test"]==cost
    return cost,labels


def stringify(obj):
    if isinstance(obj, F):
        return str(obj)
    if isinstance(obj, dict):
        return {k:stringify(v) for k,v in obj.items()}
    if isinstance(obj, list):
        return [stringify(v) for v in obj]
    return obj


def run():
    assert abs(current_curvature()-1/280)<1e-15
    m = F(1)  # TeV
    ordinary=[(F(2,3),m)]*5
    exotic=[(F(-4,3),m)]+[(F(-1,3),m)]*4
    ref, comp = extract(ordinary,[]), extract(exotic,[])
    assert ref["current_curvature"] == comp["current_curvature"]
    assert ref["all_charge_even_test"] == 0 and comp["all_charge_even_test"] == F(20,3)
    assert ref["colour_odd_test"] == 0 and comp["colour_odd_test"] == 2
    assert ref["coeff_x4_to_z4"][2] == comp["coeff_x4_to_z4"][2]
    assert ref["coeff_x4_to_z4"][4] == comp["coeff_x4_to_z4"][4]
    # 1779's |Q|<2 branch with three coloured -5/6 and two +7/6 doublets.
    shifted=[(F(-4,3),m),(F(-1,3),m)]*3+[(F(2,3),m),(F(5,3),m)]*2
    paired = extract(shifted,[])
    assert paired["colour_odd_test"] == 10 and paired["all_charge_even_test"] == F(172,3)
    irreps = [F(k)-F(1,3) for k in range(-5,6)]
    for i, q in enumerate(irreps):
        for j in range(-5,6):
            h = extract([(q,F(i+2,3))],[(F(j),F(abs(j)+3,5))])
            assert h["colour_odd_test"] >= 0 and h["all_charge_even_test"] >= 0
            assert (h["colour_odd_test"] == 0) == (q in (F(-1,3),F(2,3)))
            assert (h["all_charge_even_test"] == 0) == (abs(j)<=1 and q in (F(-1,3),F(2,3)))
    outside = extract([(F(0),m)]*9 + [(F(-4,3),m)], [])
    assert outside["colour_odd_test"] == 0
    # Use a saved threshold purely as an alternative common heavy mass.
    old=json.loads((HERE/"../1776/results.json").read_bytes())
    gate=max(map(float, old["required_projected_Y_if_only_three_below"]))/1000
    vev=json.loads((HERE/"../1754/inputs.json").read_bytes())["v_GeV"]/1000
    mass_cap=2*vev/math.sqrt(2)
    energy=1.0
    pattern=[(q,F(str(mass_cap))) for q,_ in shifted]
    floor=10*math.sqrt(1-(2*mass_cap/energy)**2)
    sharp,labels=mismatch_cost()
    assert abs(finite_cut(pattern,energy)-10*pair_kernel(energy,mass_cap))<1e-12
    assert finite_cut(pattern,energy) >= floor
    assert finite_cut(pattern,1.9*mass_cap)==0
    for beta in np.linspace(.002,.9999,120):
        m_test=math.sqrt(1-beta**2)/2
        assert pair_kernel(1.,m_test) >= beta - 1e-12
    # Compare arbitrary genuinely split masses, both below and above cut.
    split=[(q,F(i+1,7)) for i,(q,_) in enumerate(shifted)]
    direct=sum(float((q+F(1,3))*(q-F(2,3)))*pair_kernel(1.5,float(m)) for q,m in split)
    assert abs(finite_cut(split,1.5)-direct)<1e-12
    return stringify(dict(mass_unit="TeV", colour_generator=T, response_coefficients=ref,
        equal_two_point_competitor=comp, shifted_doublet_control=paired,
        off_lattice_competitor=outside, split_mass_charge_tests=len(irreps)*11,
        old_threshold_TeV=gate,
        sharp_unit_mismatch=dict(labels_q_minus1_q_plus1_l1=list(labels),
            minimum_at_one_TeV=sharp, minimum_at_old_cap=round(float(sharp)/gate**4,12)),
        timelike_cut=dict(yukawa_norm_cap=2, mass_cap_TeV=round(mass_cap,9),
            energy_TeV=energy, bound=round(floor,9), equal_cap_value=round(finite_cut(pattern,energy),9)),
        odd_exotic_signal_at_old_threshold=round(2/gate**4,12),
        even_exotic_signal_at_old_threshold=round((20/3)/gate**4,12)))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    payload=json.dumps(result,indent=2)+"\n"
    path=HERE/"independent_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
    else:
        path.write_text(payload,encoding="utf-8")
    print("1780 heavy charge moments", hashlib.sha256(payload.encode()).hexdigest()[:12])
