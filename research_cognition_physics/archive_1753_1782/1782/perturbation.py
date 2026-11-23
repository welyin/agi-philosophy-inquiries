"""Conditional operator-norm stability test for the electronic reactor spectrum."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
# Conservative vacuum spectral windows; these are not a four-state fit.
SOLAR = (6.9e-5, 8.1e-5)
ATM = (2.2e-3, 2.8e-3)
BIG, MIDDLE, TINY = (.64,.72), (.27,.33), (.018,.026)


def root(lower_weight, min_gap, higher_weight, max_gap, max_zero):
    def deficit(eps):
        return (lower_weight*math.sqrt(min_gap) - higher_weight*math.sqrt(max_gap+eps**2)
                - (1+max_zero)*eps)
    lo, hi=0., .1
    assert deficit(lo)>0>deficit(hi)
    for _ in range(80):
        mid=(lo+hi)/2
        if deficit(mid)>0:
            lo=mid
        else:
            hi=mid
    return lo


def matrix_witness(m_a, m_b, wa, wb):
    z=wa+wb
    w0=1-z
    d=(m_a+m_b)*math.sqrt(wa*wb)/z
    direct=(wb*m_b-wa*m_a)/z
    bare=(wa*m_b-wb*m_a)/z
    x=np.array([math.sqrt(z),math.sqrt(w0),0.0])
    matrix=np.zeros((4,4))
    matrix[:3,:3]=direct*np.outer(x,x)
    matrix[3,:3]=d*x
    matrix[:3,3]=d*x
    matrix[3,3]=bare
    val, vec=np.linalg.eigh(matrix)
    weights=vec[0,:]**2
    actual=[sum(weights[abs(val)<1e-12]),sum(weights[abs(val+m_a)<1e-12]),
            sum(weights[abs(val-m_b)<1e-12])]
    np.testing.assert_allclose(actual,[w0,wa,wb],atol=1e-12)
    eps=np.linalg.norm(matrix[:3,:3],2)
    np.testing.assert_allclose(eps,abs(direct))
    return dict(direct_active_norm_eV=eps, abs_ee_entry_eV=abs(matrix[0,0]),
                portal_dirac_eV=d, singlet_bare_signed_eV=bare, verified_weights=actual)


def run():
    # If the low unresolved pair is at a member of the close solar pair,
    # the middle mass carries either the big or middle electron weight.
    lows=[root(pair[0],SOLAR[0],TINY[1],ATM[1],other[1])
          for pair,other in [(MIDDLE,BIG),(BIG,MIDDLE)]]
    # If it is isolated, both nonzero roots have atmospheric-scale masses,
    # regardless of which member of the solar pair carries the largest weight.
    inverted=root(BIG[0],ATM[0],MIDDLE[1],ATM[1],TINY[1])
    central_s12=.3092
    tiny=(1-math.sqrt(1-.0851))/2
    middle,big=(1-tiny)*central_s12,(1-tiny)*(1-central_s12)
    ds=7.50e-5
    # Distances from the bottom cluster in eV^2; these are examples within
    # the registered rectangles, not a common best-fit likelihood.
    norm=matrix_witness(math.sqrt(ds),math.sqrt(.002466+ds),middle,tiny)
    inv=matrix_witness(math.sqrt(.002571-ds),math.sqrt(.002571),big,middle)
    # Explicitly check the Weyl perturbation claim including complex phases.
    rng=np.random.default_rng(1782)
    for _ in range(96):
        d=(rng.normal(size=3)+1j*rng.normal(size=3))*.02
        m0=np.zeros((4,4),complex)
        m0[3,:3]=m0[:3,3]=d
        m0[3,3]=rng.uniform(-1,1)*.1
        e=(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))*.001
        e=(e+e.T)/2
        m=m0+e
        eps=np.linalg.norm(e,2)
        assert np.sort(np.linalg.svd(m,compute_uv=False))[1] <= eps+1e-13
        assert abs(m[0,0])<=eps+1e-13
    return dict(units="eV, spectral gaps in eV^2",
        windows=dict(solar=SOLAR,atmospheric=ATM,big=BIG,middle=MIDDLE,tiny=TINY),
        close_pair_at_low_end_eV=lows,isolated_at_low_end_given_two_state_low_cluster_eV=inverted,
        smallest_branch_distance_eV=min(*lows,inverted,math.sqrt(SOLAR[0])),
        constructive_normal=norm,constructive_inverted=inv,complex_perturbation_examples=96)


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [clean(v) for v in x]
    return float(f"{x:.11g}") if isinstance(x,float) else x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=clean(run())
    payload=json.dumps(result,indent=2)+"\n"
    path=HERE/"perturbation_results.json"
    if args.check:
        assert json.loads(path.read_bytes())==result
    else:
        path.write_text(payload,encoding="utf-8")
    print("1782 reactor perturbation window",hashlib.sha256(payload.encode()).hexdigest()[:12])
