"""Tensor projection, full messenger spectrum and a tuned counterexample."""
import argparse
import json
import math
from itertools import combinations, combinations_with_replacement
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
VAC = np.diag([2.,2.,2.,-3.,-3.])


def weights(h):
    raw = np.einsum("ik,j->ijk",VAC,h)
    unit = np.eye(5)
    traceless = []
    for sign in (-1,1):
        part = (raw + sign*raw.swapaxes(0,1))/2
        contracted = np.einsum("iji->j",part)
        source = contracted/(5+sign)
        trace = np.einsum("ik,j->ijk",unit,source) + sign*np.einsum("jk,i->ijk",unit,source)
        traceless.append(part-trace)
        assert np.linalg.norm(np.einsum("iji->j",part-trace))<1e-12
    five = raw - sum(traceless)
    irreps = (five,*traceless)
    projections = np.array([a.ravel() for a in irreps])
    gram = projections @ projections.conj().T
    np.testing.assert_allclose(gram,np.diag(np.diag(gram)),atol=1e-13)
    return np.diag(gram).real


def dynkin_pair():
    # Trace one normalized Lie generator squared on each projected irrep.
    t = np.array([.5,-.5,0.,0.,0.])
    indices=[]
    for bag in (combinations(range(5),2),combinations_with_replacement(range(5),2)):
        terms=[(t[i]+t[j]-t[k])**2 for i,j in bag for k in range(5)]
        indices.append(sum(terms)-sum(t**2))
    return [2*sum(t**2),*(2*x for x in indices)]


def exact_mass(wt,wd):
    left,right,heavy=.18,.29,1.
    mu=wd*left*right/heavy
    t=np.array([[mu,left*math.sqrt(wt)],[right*math.sqrt(wt),heavy]])
    d=np.array([[mu,left*math.sqrt(wd)],[right*math.sqrt(wd),heavy]])
    st,sd=np.linalg.svd(t,compute_uv=False),np.linalg.svd(d,compute_uv=False)
    assert sd[-1]<2e-15
    schur=abs(mu-wt*left*right/heavy)
    exact=float(np.prod(st)/sd[0])
    wave=math.sqrt((1+wd*(left/heavy)**2)*(1+wd*(right/heavy)**2))
    np.testing.assert_allclose(exact,schur/wave,rtol=3e-14)
    return {"mu_over_M":mu,"triplet_masses_over_M":st.tolist(),
            "heavy_doublet_over_M":float(sd[0]),"effective_triplet_over_M":exact,
            "schur_over_M":schur,"ratio_full_to_schur":exact/schur}


def audit():
    wt,wd=weights(np.eye(5)[0]),weights(np.eye(5)[3])
    for k in range(5):
        np.testing.assert_allclose(weights(np.eye(5)[k]),wt if k<3 else wd,atol=4e-14)
    frac=wt/wd
    np.testing.assert_allclose(frac,[4/9,4/3,8/9])
    real_cost=dynkin_pair()
    assert real_cost==[1.,24.,49.]
    # Independently check the operator bounds for nondiagonal complex portals.
    rng = np.random.default_rng(73217)
    margins=[]
    for i in range(40):
        matrices=[rng.normal(size=(4,3))+1j*rng.normal(size=(4,3)) for _ in range(3)]
        forms=[a.conj().T@a for a in matrices]
        dou=sum(x*g for x,g in zip(wd,forms))
        tri=sum(x*g for x,g in zip(wt,forms))
        low=np.linalg.eigvalsh(tri-min(frac)*dou)[0]
        high=np.linalg.eigvalsh(max(frac)*dou-tri)[0]
        assert min(low,high)>-1e-11
        margins.append(min(low,high))
    # Two inequivalent holomorphic products with relative minus sign cancel only the doublet.
    coeff=np.array([1.,-wd[0]/wd[1],0.])
    assert abs(wd@coeff)<2e-14 and abs(wt@coeff)>0.5
    return {"channel_order":[5,45,70],"triplet_weights":wt.tolist(),"doublet_weights":wd.tolist(),
            "ratios":frac.tolist(),"positive_gram_test_count":len(margins),
            "SUSY_delta_b5":real_cost,"signed_cancellation":{"coefficients":coeff.tolist(),
                "triplet_exchange":float(wt@coeff),"doublet_exchange":float(wd@coeff)},
            "one_tree_messenger_exact_spectra":{str(k):exact_mass(float(t),float(d))
                for k,t,d in zip((5,45,70),wt,wd)}}

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
