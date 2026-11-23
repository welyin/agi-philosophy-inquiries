"""Gauge and right-electron tests on the split five-plet spectrum of round 1744."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
def read(p):
    return json.loads((HERE/p).read_bytes())

def load(p):
    return runpy.run_path(str(HERE/p))

SOFT=30000.
TWO_PI=2*math.pi


def gauge_and_wilson(mw):
    pars=read("../1369/inputs.json")
    base=read("../1698/results.json")["shared_witness"]
    mc=read("../1706/results.json")["cluster_reference_scales_GeV"][1]
    low=np.array([pars[k] for k in ["g1_GUT","g2","g3"]])
    gg=low[:,0]+np.array([0,0,-.67])*low[:,1]
    # Infer the crossing from a1=a2 directly; alpha3 is an independent residual.
    a=4*math.pi/gg**2
    sm=np.array([41/10,-19/6,-7.])
    mssm=np.array([33/5,1.,-3.])
    color=np.array([.4,0,1.])
    weak=np.array([.6,1.,0])
    events=[(SOFT,mssm-sm,-np.array([0.,2.,3.])/(12*math.pi))]
    for label, charge, newmass in [("colored",color,mc),("weak",weak,mw)]:
        for kind,weight in [("fermions",2/3),("scalars",1/3)]:
            extra=newmass if kind=="fermions" else math.hypot(SOFT,newmass)
            for m in [*base[f"{label}_{kind}_GeV"],extra,extra]:
                events.append((m,weight*charge,np.zeros(3)))
    old=pars["mu_ref_GeV"]
    b=sm.copy()
    gamma=np.array([[23/20,9/4,2],[23/30,3/2,4/3]])
    accum=np.zeros((2,3))
    nodes,w=np.polynomial.legendre.leggauss(20)
    def integrate(a,b,lo,hi):
        t=math.log(hi/lo)/TWO_PI
        return t/2*(w@(1/(a-(nodes[:,None]+1)*t*b/2)))
    for m,db,matching in sorted(events,key=lambda e:e[0]):
        accum[int(old>=SOFT)]+=integrate(a,b,old,m)
        a=a-b*math.log(m/old)/TWO_PI+matching
        b+=db
        old=m
    l=(a[0]-a[1])/(b[0]-b[1])
    cross=old*math.exp(TWO_PI*l)
    accum[1]+=integrate(a,b,old,cross)
    at=a-b*l
    return dict(cross_GeV=cross,a_at_cross=at.tolist(),colour_minus_weak=at[2]-at[1],
                qque_below_above_soft=np.exp(np.sum(gamma*accum,axis=1)).tolist(),
                down_to_2GeV_gain=1.25*np.exp(np.sum(gamma*accum)))


def projective_chord(v,w):
    # Stable for nearly identical normalized vectors, including arbitrary column phases.
    a=np.vdot(v,w)
    return np.linalg.norm(w-v*a/abs(a))


def at_low(bg, filter_weak):
    y=bg["yW"].copy()
    y[2]=filter_weak@y[2]
    y[3]=filter_weak@y[3]
    flow=bg["pt"]["flow"]
    tags=bg["model"]
    svd=lambda mat:np.linalg.svd(mat[:,:3])
    leg=np.eye(3,dtype=complex)
    for start,_,dt,grid in reversed(bg["segments"]):
        y,leg=tags["advance"](y,leg,-dt,grid[::-1],flow)
        if start==bg["Mcolor"]:
            up=svd(y[0])[0][:,-1]
            left,sv,vh=svd(y[1])
            c=abs(np.vdot(up,left[:,-1]))
            factors=np.array([math.sqrt((1+c)/2),1.,c])
            y[1,:,:3]=y[1,:,:3]@(vh.conj().T*factors)@vh
        for j,pole in enumerate(flow["fields"]["neutral_poles_GeV"]):
            if start==pole:
                y[3,:,j]=0
    return y


def flavour():
    bg=load("check.py")["mother"](.125)
    pt=bg["pt"]
    high=bg["high"]
    E,D=pt["E"],pt["D"]
    swapped=E.T
    leftq=bg["left_quark_to_GJ"]
    dre=np.linalg.svd(D)[0][:,-1]
    ere,es,_=np.linalg.svd(E)
    e=ere[:,-1]
    v=np.array([1,0,0],complex)
    tv=(np.eye(3)+np.outer(v,v.conj()))@dre
    norms=abs(tv*e)
    fmin=2*max(norms)-sum(norms)
    # Rotate canonical right fields so that their high tensors agree, including all phases.
    re=np.linalg.solve(high[2,:,:3],swapped)
    aa=np.linalg.svd(high[0,:,:3],compute_uv=False)[::-1]
    ru=np.linalg.solve(high[0,:,:3],leftq*aa)
    for r in [re,ru]:
        np.testing.assert_allclose(r.conj().T@r,np.eye(3),atol=2e-6)
    rows={}
    for key,F in [("no_weak_gate",np.eye(3)),("registered_common_gate",bg["F"])]:
        low=at_low(bg,F)
        dl=np.linalg.svd(leftq.conj().T@low[1,:,:3])[0][:,-1]
        el=np.linalg.svd((low[2,:,:3]@re).T)[0][:,-1]
        ul,_,ur=np.linalg.svd(leftq.conj().T@low[0,:,:3]@ru)
        # qque uses the conjugate right-handed spinor. Its distance to e1 is phase independent.
        chords=[projective_chord(dre,dl),projective_chord(e,el),projective_chord(v,ul[:,-1]),
                projective_chord(v,ur[-1,:])]
        # Two normalized contractions ((w†u)d+(w†d)u). Triangle bound uniform over phases.
        allowance=2*sum(chords)
        rows[key]=dict(chord_d_e_upL_upC=chords,possible_amplitude_shift_upper=allowance,
                       flavor_amplitude_lower_with_external_transport=fmin-allowance)
    # Radius under any single left gate; SM/MSSM RGE optional conditional spectral-norm cap.
    r=np.linalg.svd(bg["yW"][2,:,:3],compute_uv=False)[-2:]
    r=r[1]/r[0]
    c=bg["cW"]
    def angle(r,c):
        return math.degrees(math.atan(r*(1-c*c)/(2*(c*c-r*r))))
    msh=pt["flow"]["MZ"]
    cap=(3*math.log(SOFT/msh)+4*math.log(1e7/SOFT))/(16*math.pi**2)
    ctotal=c*math.exp(-cap)
    return dict(mw=bg["Mweak"],mc=bg["Mcolor"],terms=norms.tolist(),fmin=fmin,
                threshold_right_e_gate_anystage_angle_max_deg=angle(r,c),
                full_right_e_if_all_Ye_Ynu_operator_norms_le_1_angle_deg=angle(es[-1]/es[-2],ctotal),
                lower_condition_number_inverse=ctotal,external_rows=rows)


def calculate():
    xx=flavour()
    mw=xx["mw"]
    spectrum=gauge_and_wilson(mw)
    ref=read("../1742/results.json")
    target=ref["U_GeV"]
    rel=spectrum["cross_GeV"]/target-1
    assert abs(rel)<1e-8 and abs(spectrum["colour_minus_weak"])<2e-8
    # An entirely free new weak-pole hypothesis would not guarantee the same alpha3 reading.
    fail=gauge_and_wilson(xx["mc"])
    assert abs(fail["colour_minus_weak"])>.01
    canonical=4*math.pi/(spectrum["a_at_cross"][0]*spectrum["cross_GeV"]**2)
    old=load("../1725/independent_check.py")
    summary={}
    for label,amp in [("high_family_only",xx["fmin"]),("registered_common_gate_with_transport",
                       xx["external_rows"]["registered_common_gate"]["flavor_amplitude_lower_with_external_transport"])]:
        entry={}
        for lower,key in [(0,"central_hadron"),(1,"rectangular_hadron_floor")]:
            cap=old["coefficients"]("pi_e",lower)
            gain=spectrum["down_to_2GeV_gain"]
            entry[key]=dict(U_coefficient_limit=cap/(amp*gain),
                            rate_to_single_mode_limit_lower=(amp*gain*canonical/cap)**2)
        summary[label]=entry
    return dict(new_spectrum=spectrum,
                wrong_equal_weak_colour_Dirac_control=fail,initial_family=xx,
                qque_increment_over_round1742_equal_weak=spectrum["down_to_2GeV_gain"]/ref["degenerate_weak"]["total_qque"]-1,
                finite_e_bounds=summary)


def tidy(o):
    if isinstance(o,dict):return {k:tidy(v) for k,v in o.items()}
    if isinstance(o,list):return [tidy(v) for v in o]
    if isinstance(o,float):return float(f"{o:.9g}")
    return o


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("--check",action="store_true")
    args=a.parse_args()
    out=tidy(calculate())
    path=HERE/"independent_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1744 independent OK: split thresholds, qque transport, right electron and finite final state.")
    else:
        path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")

