"""Finite-band charged-current index and fourth-family Yukawa budget at the measured-input mother."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE = Path(__file__).resolve().parent
PORT = runpy.run_path(str(HERE.parent / "1748/check.py"))


def svd(m):
    l,s,rh = np.linalg.svd(m)
    return l[:,::-1],s[::-1],rh.conj().T[:,::-1]


def lower_fraction(eu,ed):
    amp=max(0.,math.sqrt((1-eu**2)*(1-ed**2))-eu*ed)
    return amp**2-(eu*ed)**2


def chiral_measure(mu,md,wl,wr,n):
    lu,su,ru=svd(mu)
    ld,sd,rd=svd(md)
    gl=lu.conj().T@wl@ld
    gr=ru.conj().T@wr@rd
    def strength(left,right):
        return np.linalg.norm(left,"fro")**2-np.linalg.norm(right,"fro")**2
    whole=strength(gl,gr)
    np.testing.assert_allclose(whole,n,atol=1e-12)
    return strength(gl[:n,:n],gr[:n,:n]),su,sd


def test_gap():
    rng=np.random.default_rng(1754)
    checks=0
    smallest=1.
    for n in (3,4,6):
        q,s=2,2
        nd=n+q+s
        ldouble=np.arange(nd)<n+q
        rdouble=(np.arange(nd)>=n)&(np.arange(nd)<n+q)
        wl=np.diag(ldouble.astype(float))
        wr=np.diag(rdouble.astype(float))
        allowed=ldouble[:,None]!=rdouble[None,:]
        for t in (.02,.10,.20,.28,.39):
            for rep in range(8):
                mats=[]
                for pole_scale in (1.,1.2):
                    m0=np.diag(np.r_[np.zeros(n),[1.,1.7],pole_scale*np.array([1.2,1.9])])
                    y=(rng.normal(size=(nd,nd))+1j*rng.normal(size=(nd,nd)))*allowed
                    mats.append(m0+t*y/np.linalg.norm(y,2))
                result,su,sd=chiral_measure(*mats,wl,wr,n)
                assert all(v[n-1]<=t+1e-14 and v[n]>=1-t-1e-14 for v in (su,sd))
                e=t/(1-t)
                assert n*lower_fraction(e,e)-1e-12<=result<=n+1e-12
                smallest=min(smallest,result/n)
                checks+=1
    return checks,smallest


def same_light_index():
    # Four independently mixed left doublets, vectorlike singlets only. Dimensionful unit M=1.
    # Tiny unequal ordinary masses split the four light eigenvalues.
    n=4
    wl=np.diag(np.r_[np.ones(n),np.zeros(n)])
    wr=np.zeros_like(wl)
    def build(t,m):
        return np.block([[np.diag(m),t*np.eye(n)],[np.zeros((n,n)),np.eye(n)]])
    u=np.array([.001,.002,.003,.007])
    d=np.array([.0005,.0011,.0033,.005])
    def evaluate(t):
        return chiral_measure(build(t,u),build(t,d),wl,wr,n)[0]
    lo,hi=0.,.49
    for _ in range(60):
        mid=(lo+hi)/2
        if evaluate(mid)>3:
            lo=mid
        else:
            hi=mid
    t=(lo+hi)/2
    measured=evaluate(t)
    np.testing.assert_allclose(measured,3.,atol=2e-12)
    return dict(t=t,finite_band_residue=measured,full_residue=4.,
                largest_breaking_to_unbroken_gap=math.hypot(t,max(u)))


def solve():
    b=PORT["OLD"]["mother"](.125)
    v=json.loads((HERE.parent/"1749/inputs.json").read_bytes())["v_GeV"]
    prep=json.loads((HERE.parent/"1753/independent_results.json").read_bytes())
    r=prep["fourth_response_R"]
    weighted=prep["two_sources"]["quark"]["exact_large_R_crossing_for_rate_preparation"]
    data={}
    for f,ys,nc in zip(("up","down","charged_leptons"),b["yW"][:3],(3,3,1)):
        y=np.linalg.svd(ys[:,:3],compute_uv=False)[::-1]
        def read(r):
            return math.sqrt(y[0]**2+r*(y[-1]**2-y[0]**2))
        extra=read(r)
        # This is a declared perturbative power-counting budget, not a universal unitarity bound.
        maxr=(16*math.pi**2/nc-y[0]**2)/(y[-1]**2-y[0]**2)
        data[f]=dict(reference_y_min_max=y[[0,-1]].tolist(),fourth_y_if_universal_R=extra,
                     vev_times_y_over_sqrt2_GeV=extra*v/math.sqrt(2),
                     minimal_y_at_down_reversal_threshold=read(weighted),
                     color_weighted_y_squared_over_16pi2=nc*extra**2/(16*math.pi**2),
                     maximal_R_in_unit_loop_budget=maxr)
    assert data["up"]["maximal_R_in_unit_loop_budget"]<weighted<r
    assert r<min(data[k]["maximal_R_in_unit_loop_budget"] for k in ("down","charged_leptons"))
    eps=math.sqrt((4-math.sqrt(13))/6)
    cases,ratio=test_gap()
    return dict(matching_scale_GeV=b["Mweak"],R=r,down_sign_threshold_R=weighted,
                separated_sectors=data,random_mass_checks=cases,
                smallest_sampled_low_index_divided_by_true_count=ratio,
                maximal_common_deltaM_ratio_sufficient_for_four_to_exceed_three=eps/(1+eps),
                equal_low_sum_with_different_full_sum=same_light_index())


def clean(o):
    if isinstance(o,dict):
        return {k:clean(v) for k,v in o.items()}
    if isinstance(o,list):
        return [clean(x) for x in o]
    if isinstance(o,(np.floating,float)):
        return float(f"{o:.11g}")
    return o


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=clean(solve())
    path=HERE/"independent_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1754 complete and restricted weak indices, gap witnesses and common fourth coordinate verified")
    else:
        with path.open("x",encoding="utf8") as out_file:
            json.dump(out,out_file,indent=2,allow_nan=False)
            out_file.write("\n")
