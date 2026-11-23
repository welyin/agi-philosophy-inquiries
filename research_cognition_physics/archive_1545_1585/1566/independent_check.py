"""Solve the two odd-return equations for the actual three mass tags."""
import argparse
import json
from fractions import Fraction as F
from math import acos, cos, sin, sqrt, pi
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def independent_velocity(high, cfg):
    def position(reading, key):
        return reading[0]+reading[1]*cfg["width_coordinate"][key]
    aa,bb,cc,delta = high["angles_theta12_theta13_theta23_delta"]
    angles = [position(aa,"theta12"),position(bb,"theta13"),cc[0]]
    a,b,c = angles
    def rotation(i,j,angle,phase=0.):
        m = np.eye(3,dtype=complex)
        m[i,i] = m[j,j] = cos(angle)
        m[i,j] = sin(angle)*np.exp(-1j*phase)
        m[j,i] = -m[i,j].conj()
        return m
    up = np.diag([v[0] for v in high["up"]])
    signs = np.asarray(cfg["signs"])
    velocities = []
    for phase in [0.,pi,delta[0],pi/3]:
        v = rotation(1,2,c)@rotation(0,2,b,phase)@rotation(0,1,a)
        weights = abs(v[0])**2
        ms = position(high["down"][1],"strange")
        mb = position(high["down"][2],"bottom")
        first = -(weights[1]*ms*signs[1]+weights[2]*mb*signs[2])/(signs[0]*weights[0])
        y = v@np.diag([first,ms,mb])
        diff = 1.5*(y@y.conj().T-up@up.T)/(16*pi**2)
        def read_F(step):
            # Identity pieces cancel at F=0. Re-diagonalize both physical Yukawas.
            axes, _, _ = np.linalg.svd((np.eye(3)-step*diff)@up)
            left, masses, _ = np.linalg.svd((np.eye(3)+step*diff)@y)
            return (abs(axes[:,-1].conj()@left[:,::-1])**2)@(signs*masses[::-1])
        def central(step):
            return (read_F(step)-read_F(-step))/(2*step)
        estimate = (4*central(.02)-central(.04))/3
        velocities.append(float(estimate))
    saved = read(HERE/"results.json")
    assert abs(velocities[2]-saved["initial"]["derivative_total"]) < 5e-15
    expected = (velocities[0]+velocities[1]+cos(pi/3)*(velocities[0]-velocities[1]))/2
    assert abs(velocities[3]-expected) < 5e-15
    return velocities


def certificate(down,angles,angles_are_sines):
    masses = [F(str(x[0])) for x in down]
    d,s,b = masses
    # Exact common null of the one- and three-step amplitudes.
    kernel = [b**2-s**2, -(b**2-d**2), s**2-d**2]
    weights = [abs(k/m) for k,m in zip(kernel,masses)]
    z = sum(weights)
    p = [w/z for w in weights]
    sig = [1,-1,1]
    odd = [sum(w*e*m**n for w,e,m in zip(p,sig,masses)) for n in [1,3,5]]
    assert odd[0]==odd[1]==0 and odd[2]!=0
    assert p[2]/p[1] == s*(s**2-d**2)/(b*(b**2-d**2))
    assert p[2]/p[1] < (s/b)**3
    c12,c13 = angles[:2]
    def trig(x):
        return x if angles_are_sines else sin(x)
    s12,s13 = trig(c12[0]),trig(c13[0])
    data = [(1-s12**2)*(1-s13**2),s12**2*(1-s13**2),s13**2]
    overlap = sum(sqrt(float(q)*r) for q,r in zip(p,data))
    # In the next bounds no value for the smallest mass or its sign is used.
    s12hi = trig(c12[0]+2*c12[1])
    s13lo = trig(c13[0]-2*c13[1])
    rate_lo = s13lo/(s12hi*sqrt(1-s13lo**2))
    sb_max = (down[1][0]+2*down[1][1])/(down[2][0]-2*down[2][1])
    upper = sb_max**1.5
    leakage = 1-sb_max**2-sb_max**3/rate_lo**2
    assert rate_lo>upper and leakage>0
    return {"opposite_middle_sign_needed":True,"exact_double_null_probabilities":[float(q) for q in p],
            "nonzero_five_step_return_abs":float(abs(odd[2])),
            "minimum_ray_angle_to_double_null_center_deg":acos(min(1,overlap))*180/3.141592653589793,
            "heavy_to_middle_weak_modulus_ratio_data_lower":rate_lo,
            "double_null_modulus_ratio_upper_for_every_first_mass":upper,
            "minimum_ratio_discrepancy_factor":rate_lo/upper,
            "if_one_step_zero_three_step_abs_over_pb_b_cubed_lower":leakage}


def verify(x,y):
    if isinstance(x,dict):
        assert set(x)==set(y)
        for name in x:
            verify(x[name],y[name])
    elif isinstance(x,list):
        assert len(x)==len(y)
        for a,b in zip(x,y):
            verify(a,b)
    else:
        assert abs(x-y) <= 2e-10*max(abs(x),abs(y),1e-35)


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    flag=parser.parse_args().check
    high=read(HERE.parent/"1397"/"inputs.json")["high_scale"]
    low=read(HERE.parent/"1386"/"inputs.json")
    result={"high":certificate(high["down"],high["angles_theta12_theta13_theta23_delta"],False),
            "direct_yukawa_derivatives_at_zero_pi_observed_pi_over_three":independent_velocity(high,read(HERE/"inputs.json")),
            "low":certificate(low["yukawa"]["down"],low["quark_sines_and_phase"],True)}
    path=HERE/"independent_results.json"
    if flag:
        verify(result,read(path))
        print("1566 signed moments, two-scale current and nonzero three-step return checked")
    else:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
