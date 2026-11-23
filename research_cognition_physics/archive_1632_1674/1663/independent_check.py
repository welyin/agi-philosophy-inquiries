"""Independent one-loop joint-threshold test with sharp mass-order bounds."""
import argparse
from fractions import Fraction
from itertools import product
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
MATCH_SCALE = 1e16
Q = Fraction
SM = (Q(41,10), Q(-19,6), -7)
D = (Q(4,15), 0, Q(2,3))
L = (Q(2,5), Q(2,3), 0)
TRIP = (Q(3,5), Q(2,3), 0)
BI = (Q(1,30), Q(1,2), Q(1,3))
SEXT = (Q(8,15), 0, Q(5,6))


def read(rel):
    return json.loads((HERE.parent / rel).read_bytes())


def projections(v):
    return v[0]-v[1], v[2]-v[1]


def solve(inverse_couplings, k, running):
    a, weak, color = inverse_couplings
    r = 2*math.pi*(a-weak)-float(SM[0]-SM[1])*running+4*math.log(k)/5
    t = 2*math.pi*(color-weak)-float(SM[2]-SM[1])*running+2*math.log(k)
    return (-25*r+16*t)/9, (-20*r+2*t)/9


def ordered_bounds(mz, x):
    slack = 4*math.log(MATCH_SCALE/mz)-x
    return [mz*math.exp(slack/(5-j)) for j in range(1,5)]


def integrate(inverse, k, b_masses, x, z, mz):
    # Within-copy scalar mass splittings need not be identical; projected sums use log means.
    means = math.prod(b_masses)**0.25
    m_delta, m_sext = means*math.exp(x/4), means*math.exp(z/4)
    assert min(b_masses) >= mz and max(*b_masses,m_delta,m_sext) <= MATCH_SCALE*(1+1e-12)
    logs = [math.log(MATCH_SCALE/mz), 3*math.log(MATCH_SCALE/(1e8*k)), 3*math.log(MATCH_SCALE/1e8),
            4*math.log(MATCH_SCALE/m_delta), sum(math.log(MATCH_SCALE/m) for m in b_masses),4*math.log(MATCH_SCALE/m_sext)]
    out = np.asarray(inverse)-np.asarray([SM,D,L,TRIP,BI,SEXT],dtype=float).T@logs/(2*math.pi)
    assert np.ptp(out) < 3e-12 and min(out)>0
    return {"B_GeV":b_masses,"Delta_GeV":m_delta,"S_GeV":m_sext,"alpha_U_inverse":float(np.mean(out))}


def results():
    src = read("1369/inputs.json")
    mz = src["mu_ref_GeV"]
    run = math.log(MATCH_SCALE/mz)
    input_gain = read("1662/independent_results.json")["two_width_flavour_gain"]
    endpoints = [1.30003,15.701]
    assert endpoints[0]<input_gain[0]<input_gain[1]<endpoints[1]
    tableau = [src[key] for key in ("g1_GUT","g2","g3")]
    w = np.array([5,-3,-2])
    cols = (SM,D,L,TRIP,BI,SEXT)
    assert list(w@np.asarray(cols,dtype=object).T)==[44,0,0,1,-2,1]
    assert tuple(sum(c[j] for c in (TRIP,BI,SEXT)) for j in range(3))==(Q(7,6),)*3
    p1,p2 = projections(TRIP), projections(SEXT)
    prod_matrix = np.array([p1,p2],dtype=object).T@np.array([[-25,16],[-20,2]],dtype=object)
    assert np.array_equal(prod_matrix,-9*np.eye(2,dtype=int))
    a0 = [4*math.pi/g**2 for g,_ in tableau]
    witness = read("1657/independent_results.json")["extended_box_GJ_test"]
    k0 = witness["k_times_yb"]/witness["marginal_coordinates"][2]
    x0,z0 = solve(a0,k0,run)
    wx, wz, sum_ratio = [], [], []
    for signs in product((-2,2), repeat=3):
        a = [4*math.pi/(g+sign*e)**2 for (g,e),sign in zip(tableau,signs)]
        for gain in endpoints:
            x,z = solve(a,gain,run)
            assert math.isclose(x+z,44*run-2*math.pi*np.dot(w,a),rel_tol=3e-14)
            wx.append(x)
            wz.append(z)
            sum_ratio.append((x+z)/(2*run))
            integrate(a,gain,[mz]*4,x,z,mz)
    assert 0 < min(wz) <= max(wz) < min(wx) <= max(wx) < 4*run
    assert 3 < min(sum_ratio) <= max(sum_ratio) < 4
    sharp = ordered_bounds(mz,min(wx))
    w_only = MATCH_SCALE*math.exp(-2*run*min(sum_ratio)/8)
    assert sharp[0] < 520 and w_only < 2910
    # Saturate each sorted order statistic: j-1 copies at MZ, rest share the log budget.
    for j,ceiling in enumerate(ordered_bounds(mz,x0),start=1):
        spectrum = [mz]*(j-1)+[ceiling]*(5-j)
        demo = integrate(a0,k0,spectrum,x0,z0,mz)
        assert math.isclose(demo["Delta_GeV"],MATCH_SCALE,rel_tol=3e-14)
    example = integrate(a0,k0,[100.0,100.0,100.0,8000.0],x0,z0,mz)
    # A shared summed type-II mass need not fix normalized flavour in each multiplet.
    u = np.array([[0.8,0.6,0],[-0.6,0.8,0],[0,0,1.]])
    tensors = [m*np.outer(u[:,j],u[:,j]) for j,m in enumerate([.01,.02,.05])]
    def electron(t):
        gram = t.conj().T@t
        return gram[0,0]/gram.trace()
    electronic = [float(electron(t)) for t in tensors]
    total = float(electron(sum(tensors)))
    assert all(abs(f-total)>0.01 for f in electronic)
    return {"threshold_free_n_lower_bound":[min(sum_ratio),max(sum_ratio)],"whole_box_X":[min(wx),max(wx)],
            "whole_box_Z":[min(wz),max(wz)],"four_copy_capacity":4*run,
            "w_only_bound_on_B_geometric_mean_GeV":w_only,
            "four_B_order_statistic_upper_bounds_GeV":sharp,"physical_central_gain":k0,
            "central_XZ":[x0,z0],"central_B_order_upper_bounds_GeV":ordered_bounds(mz,x0),
            "explicit_distinct_copy_spectrum":example,
            "three_rank_one_tensors":{"electron_shares":electronic,"electron_share_of_sum":total}}


def clean(obj):
    if isinstance(obj,dict): return {key:clean(val) for key,val in obj.items()}
    if isinstance(obj,list): return list(map(clean,obj))
    if isinstance(obj,float): return float(f"{obj:.11g}")
    return obj


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    output=HERE/"independent_results.json"
    data=clean(results())
    if args.write_new:
        with output.open("x",encoding="utf-8") as out:
            out.write(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    else:
        assert json.loads(output.read_bytes())==data
        print("PASS 1663 independent double projection, complete gauge boxes, four ordered mass bounds and flavour witness")
