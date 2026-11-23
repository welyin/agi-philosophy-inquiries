"""A third isospin observable from charged/neutral K semileptonic data."""
import argparse
import cmath
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def epsilon(S, r):
    return 0.5*math.atan(math.sqrt(3)*(1-r)/((S-1)*(1+r)))


def delta_lo(S, r):
    e=epsilon(S,r)
    return math.cos(e)+math.sqrt(3)*math.sin(e)-1


def delta_nlo(S, r, mass_ratio_sq):
    Q2=(S*S-1)*(1+r)/(4*(1-r))
    return 0.75*(mass_ratio_sq+0.219*(1+S)/2)/Q2


def ratio_for_lo(S, delta):
    # On the physical (small-mixing) branch, cos e + sqrt(3) sin e = 1+delta.
    e=math.pi/3-math.acos((1+delta)/2)
    x=(S-1)*math.tan(2*e)/math.sqrt(3)
    return (1-x)/(1+x)


def fit_cofactor(seed, target):
    def apply(k):
        u,d,s=seed
        return (u+k*d*s, d+k*u*s, s+k*u*d)
    left,right=0.,5.
    for _ in range(65):
        mid=(left+right)/2
        u,d,s=apply(mid)
        if 2*s/(u+d)>target:
            left=mid
        else:
            right=mid
    return (left+right)/2,apply((left+right)/2)


def run():
    S,dS1,dS2=load("../1440/inputs.json")["ms_over_average"]
    sratio_hi=S+2*math.hypot(dS1,dS2)
    ms,dms=load("../1480/inputs.json")["lattice"]["ms"]
    mshi=(ms+2*dms)*.001
    beta,db1,db2=load("../1500/inputs.json")["sea_slope_GeV_inv"]
    rhi=(beta+2*(db1+db2))*mshi
    data=.027
    sigma=.004
    upper=data+2*sigma
    # BMW pure-QCD isoscalar tuning inputs in MeV. Same convention as their S.
    mK,emK=494.2,.5
    mpi,empi=134.8,.3
    A=((mK-2*emK)/(mpi+2*empi))**2
    B=((mK+2*emK)/(mpi-2*empi))**2
    seed_y=load("../1440/results.json")["per_scale"]["MZ"]["0"]["ordered_positive_root_ranges"][0]
    down=load("../1386/inputs.json")["yukawa"]["down"]
    u0=sum(seed_y)/2*ms*.001/down[1][0]
    d0=down[0][0]*ms*.001/down[1][0]
    seed=[u0,d0,ms*.001]
    k,eff=fit_cofactor(seed,S)
    r=eff[0]/eff[1]
    zero_LO=delta_lo(sratio_hi,rhi)
    zero_NLO=delta_nlo(sratio_hi,rhi,A)
    witness_LO=delta_lo(S,r)
    witness_NLO=[delta_nlo(S,r,a) for a in (A,B)]
    assert data-2*sigma<witness_LO<upper
    assert data-2*sigma<min(witness_NLO)<max(witness_NLO)<upper
    assert delta_lo(sratio_hi,0)>zero_LO>upper
    assert zero_NLO>upper
    inverse={}
    for rem in (0,.005,.01):
        allowed=ratio_for_lo(sratio_hi,upper+rem)
        assert math.isclose(delta_lo(sratio_hi,allowed),upper+rem,abs_tol=1e-12)
        inverse[str(rem)]={"minimum_effective_up_down":allowed,
                           "minimum_kappa_if_no_direct_up_GeV_inv":allowed/mshi}
    # Identical real effective masses at theta=0, different physical theta source.
    for theta in (.19,-.93,2.71):
        soft=cmath.exp(-1j*theta)*k*seed[1]*seed[2]
        hard=k*seed[1]*seed[2]
        assert abs(cmath.phase(soft)+theta)<1e-12 and cmath.phase(hard)==0
    return {"K_ratio_experiment":data,"K_ratio_sigma":sigma,"S_upper":sratio_hi,
            "effective_r_max_for_all_induced_positive_sea":rhi,"isospin_pure_QCD_K2_pi2_band":[A,B],
            "lo_zero_without_induction_min":delta_lo(sratio_hi,0),
            "all_induced":{"LO_min":zero_LO,"LO_required_signed_residual_max":upper-zero_LO,
                           "NLO_min_quoted_truncation":zero_NLO,"NLO_required_signed_residual_max":upper-zero_NLO},
            "LO_inverse_for_bounded_remainder":inverse,
            "partial_witness":{"kappa_GeV_inv":k,"effective_r":r,"K_LO":witness_LO,"K_NLO_range":witness_NLO,
                               "effective_Q":math.sqrt((S*S-1)*(1+r)/(4*(1-r)))}}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[key],b[key]) for key in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    out=run()
    path=HERE/"independent_results.json"
    if args.check:
        assert same(load(path.name),out)
        print("1620 independent OK: charged/neutral kaons, finite mixing, all-induced box, partial witness")
    else:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
