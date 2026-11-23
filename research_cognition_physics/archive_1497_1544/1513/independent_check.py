"""Positive four-simplex weights: sharp axis tilt and independent Bs clock."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE=Path(__file__).resolve().parent


def read(round_id):
    return json.loads((BASE.parent/str(round_id)/"inputs.json").read_text(encoding="utf-8-sig"))


def theta(r):
    if not 0 <= r < .5:
        raise ValueError("the isolated heavy branch needs 0 <= m2/m3 < 1/2")
    v=math.sqrt(2)*(1-r)
    return math.atan(r/(v+math.sqrt(v*v-r)))


def frame():
    q=[np.array([0.,0.,1.])]
    for k in range(3):
        t=2*math.pi*k/3
        q.append(np.array([math.sqrt(8)*math.cos(t)/3,math.sqrt(8)*math.sin(t)/3,-1/3]))
    q=np.array(q).T
    assert np.allclose(q.T@q,4*np.eye(4)/3-np.ones((4,4))/3)
    return q


def positive_weights(m):
    m1,m2,m3=m
    s=3*m1/4
    t=m2+m3-2*m1
    product=9*(m2-m1)*(m3-m1)/8
    disc=math.sqrt(t*t-4*product)
    return np.array([s+(t+disc)/2,s+2*product/(t+disc),s,s])


def axis_test():
    rays=frame()
    tests=0
    for b in (0.,.0007,.1,.35,.6,.95):
        for c in (0.,.001,.2,.45):
            for d in (0.,.003,.25):
                weights=np.array([1.,b,c,d])
                m,v=np.linalg.eigh((rays*weights)@rays.T)
                r=max(0.,m[-2]/m[-1])
                if r >= .5:
                    continue
                angle=math.atan2(np.linalg.norm(v[:2,-1]),abs(v[-1,-1]))
                assert angle <= theta(r)+3e-14
                if c==d==0:
                    assert abs(angle-theta(r))<4e-14
                tests+=1
    return tests


def audit():
    low=read(1386)["yukawa"]
    high=read(1397)["high_scale"]
    tree=read(1412)["BelleII_Vcb_BGL"]
    min_cb=tree[0]-2*sum(tree[1:])
    weak=read(1404)
    ca=weak["const_GeV"]
    lattice=read(1407)
    ft,ef=lattice["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    s0,es=weak["short_distance"]["S0_published"]
    eta,ee=weak["short_distance"]["etaB"]
    tau=read(1468)["omega_ps_inverse"]
    obs=tau[0]-2*sum(tau[1:])
    scenarios={}
    for label,data in (("MZ",low),("common_SM_1e16",high)):
        maxratios=[(data[k][1][0]+2*data[k][1][1])/(data[k][2][0]-2*data[k][2][1])
                   for k in ("up","down")]
        cone=sum(map(theta,maxratios))
        entry={"mass_ratio_upper":maxratios, "max_heavy_angle_rad":cone,
               "Vcb_upper":math.sin(cone), "J_upper_any_phase":math.sin(cone)**2*math.cos(cone)/4}
        if label=="MZ":
            lam=math.sin(cone)*math.cos(cone)
            coef=ca["GF"]**2*ca["MW"]**2*ca["M_Bs"]/(6*math.pi**2*lattice["hbar_GeV_ps"])
            top=coef*(ft+2*ef)**2*(eta+2*ee)*(s0+2*es)*lam*lam
            entry.update({"direct_Vcb_lower":min_cb,"minimum_principal_source_angle_deg":
                          math.degrees(math.asin(min_cb)-cone), "lambda_t_upper":lam,
                          "top_box_frequency_upper_ps_inverse":top,
                          "observed_clock_lower_ps_inverse":obs,
                          "min_other_mixing_in_frequency_units":obs-top,
                          "min_other_fraction_of_observed":1-top/obs})
        scenarios[label]=entry
    rays=frame().astype(complex)
    weights=[positive_weights([m for m,e in low[k]]) for k in ("up","down")]
    alpha=.23
    t=np.array([[math.cos(alpha),math.sin(alpha)],[-math.sin(alpha),math.cos(alpha)]],complex)
    rotate=np.eye(3,dtype=complex)
    rotate[:2,:2]=t@np.diag([1j,np.exp(.71j)])
    hz=(rays*weights[0])@rays.T
    # Up Takagi tensor in a real tetrahedral basis; down Hermitian left mass.
    hd=rotate@((rays*weights[1])@rays.T)@rotate.conj().T
    mu,uu=np.linalg.eigh(hz)
    md,ud=np.linalg.eigh(hd)
    assert np.allclose(mu,[m for m,e in low["up"]])
    assert np.allclose(md,[m for m,e in low["down"]])
    assert min(map(min,weights))>0 and np.allclose(rotate[:,2],rays[:,0])
    v=uu.conj().T@ud
    j=(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj()).imag
    assert abs(j)>1e-10 and abs(v[1,2])<scenarios["MZ"]["Vcb_upper"]
    return {"grid_and_equality_tests":axis_test(),"two_width_bounds":scenarios,
            "positive_common_axis_exact_central_mass_example":{
            "up_weights":weights[0].tolist(),"down_weights":weights[1].tolist(),
            "masses_up":mu.real.tolist(),"masses_down":md.real.tolist(),"J":float(j),
            "Vcb":float(abs(v[1,2])),"Vub":float(abs(v[0,2])),
            "same_eigenvalues_conjugate_J":float(-j)}}


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
