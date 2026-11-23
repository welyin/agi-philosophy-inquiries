"""Two-heavy-flavour loss in a common five-bar gauge current."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial

HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())
MAT = runpy.run_path(str(HERE/CFG["matter"]))
OLD = runpy.run_path(str(HERE/CFG["proton_rows"]))


def tagged_frame(row, name):
    values = np.array([row["moduli"][k] for k in ("Vus","Vub","Vcb")])
    trial=values.copy()
    for _ in range(3):
        actual=MAT["eval_point"](np.log(trial),name,.25)
        v=np.array([actual["moduli"][k] for k in ("Vus","Vub","Vcb")])
        if max(abs(v-values))<2e-10:
            break
        trial+=values-v
    assert max(abs(v-values))<2e-10
    bg,filt,_=MAT["load_mother"](*trial,.25)
    y=bg["yW"].copy()
    y[2]=actual["tau_gain"]*filt.conj().T@y[2]
    y[3]=filt.conj().T@y[3]
    transport=np.eye(3,dtype=complex)
    for start,susy,dt,grid in reversed(bg["segments"]):
        y,transport=bg["model"]["advance"](y,transport,-dt,grid[::-1],bg["pt"]["flow"])
        if start==bg["Mcolor"]:
            break
        for j,m in enumerate(bg["pt"]["flow"]["fields"]["neutral_poles_GeV"]):
            if start==m:
                y[3,:,j]=0.
    assert start==bg["Mcolor"]
    u,_,vh=np.linalg.svd(y[1,:,:3])
    lu=np.linalg.svd(y[0,:,:3])[0]
    ckabs=abs(lu[:,-1].conj()@u[:,::-1])
    c=ckabs[0]
    gain=np.array([{"c":c,"h":math.sqrt((1+c)/2),"1":1}[k] for k in name])
    # The ordinary right down quark is conjugated in the common left-handed five-bar.
    q=(bg["right_quark_to_GJ"].conj().T@vh.conj().T[:,::-1]).conj()
    np.testing.assert_allclose(q.conj().T@q,np.eye(3),atol=3e-11)
    vals,vecs=np.linalg.eigh(filt)
    assert abs(vals[-1]-1)<5e-12
    return q,gain,vals[:2],vecs[:,:2],ckabs


def bounding_loss(abs_d,s,f):
    # sqrt-coordinate magnitudes (s cos theta,s sin theta); phases extremize separately.
    if s<1e-14:
        return [float(abs_d@abs_d)]*2
    ranges=[]
    x=Polynomial([0,1])
    for sign in (-1,1):
        p,q=abs_d*f*s
        diff=s*s*(f[1]**2-f[0]**2)
        deriv=2*diff*x*(1-x*x)+sign*(-2*p*x+q*(1-x*x))*(1+x*x)
        candidates=[0.,1.]
        for root in deriv.roots():
            if abs(root.imag)<2e-9 and 0<root.real<1:
                candidates.append(float(root.real))
        scores=[]
        for t in candidates:
            cs=np.array([1-t*t,2*t])/(1+t*t)
            scores.append(float(np.sum((abs_d+sign*s*f*cs)**2)))
        ranges.append(min(scores) if sign<0 else max(scores))
    return ranges


def modes(q,g,f,weak):
    d=np.sqrt(1-f*f)[:,None]*(weak.conj().T@q)*g
    s=np.sqrt(1-g*g)
    exact={}
    for j,name in enumerate("dsb"):
        lo,hi=bounding_loss(abs(d[:,j]),s[j],f)
        exact[name]=[1-hi,1-lo]
    support=np.flatnonzero(s>0)
    j0=int(np.flatnonzero(s==0)[0])
    h=np.zeros((2,3),complex)
    h[np.arange(2),support]=s[support]
    traces=[]
    for U in (np.eye(2),np.diag([-1,1]),np.array([[0,1],[-1,0]]),np.array([[1,1j],[1j,1]])/math.sqrt(2)):
        loss=d-f[:,None]*(U@h)
        H=np.eye(3)-loss.conj().T@loss
        assert np.linalg.eigvalsh(H)[0]>=-1.e-12
        assert abs(np.linalg.eigvalsh(H)[-1]-1)<2e-12
        for j,name in enumerate("dsb"):
            lo,hi=exact[name]
            assert lo-1e-12 <= H[j,j].real <= hi+1e-12
        assert abs(H[j0,j0]-exact["dsb"[j0]][0])<3e-12
        traces.append(float(np.trace(H).real))
    return exact,traces


def calculate():
    refits=json.loads((HERE/CFG["six_refits"]).read_bytes())["refits"]
    proton=json.loads((HERE/CFG["proton_inputs"]).read_bytes())
    ref=np.array(OLD["rows"](proton,0,False))
    angles=OLD["load"](proton["angles"])["high_scale"]["angles_theta12_theta13_theta23_delta"]
    orig=np.cos(angles[1][0])*np.array([math.cos(angles[0][0]),math.sin(angles[0][0]),math.cos(angles[0][0])])
    results={}
    for name,row in refits.items():
        q,g,f,weak,ckabs=tagged_frame(row,name)
        rows,test_traces=modes(q,g,f,weak)
        A,B,C=ref*ckabs[[0,1,0]]/orig
        floor=math.cos(math.acos(min(g))+math.acos(min(f)))**2
        norms=np.array([A*A,B*B+C*C])
        min_bounds=np.array([rows["d"][0],floor])*norms
        max_bounds=np.array([rows["d"][1],1.])*norms
        # Pion dominates even the most favorable kaon rate throughout the common U(2).
        assert min_bounds[0] > max_bounds[1]
        lower=max(min_bounds)**.25*1e15
        upper=max(max_bounds)**.25*1e15
        results[name]={"tagged_visible_ranges":rows, "unitary_test_traces":test_traces,
          "any_coherent_flavour_floor":floor, "threshold_first_CKM_moduli":ckabs.tolist(),
          "unit_response_rates_at_1e15":norms.tolist(),
          "lower_upper_bracket_on_required_scale_GeV":[lower,upper],
          "pi_only_required_scale_range_GeV":(A**.5*np.array(rows["d"])**.25*1e15).tolist(),
          "kaon_lifetime_gain_necessary_for_possible_dominance":min_bounds[0]/max_bounds[1],
          "kaon_lifetime_gain_sufficient_for_universal_dominance":max_bounds[0]/min_bounds[1],
          "possible_K_to_pi_flavour_factor_envelope":[floor/rows["d"][1],1/rows["d"][0]]}
    return results


def rounded(value):
    if isinstance(value,dict):
        return {k:rounded(v) for k,v in value.items()}
    if isinstance(value,list):
        return [rounded(v) for v in value]
    return float(f"{value:.10g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=rounded(calculate())
    path=HERE/"results.json"
    if args.check:
        if result!=json.loads(path.read_bytes()):
            raise SystemExit("1768 recalculation differs")
        print("1768 five-bar and two-search checks passed")
    else:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
