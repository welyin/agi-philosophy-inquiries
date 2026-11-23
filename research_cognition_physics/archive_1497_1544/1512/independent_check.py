"""Weighted one-loop Ward identity and explicit finite-model matter totals."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path

BASE = Path(__file__).resolve().parent
DIMS = (10,5,5,5)
TWICE_INDEX = (3,1,1,1)


def block(xs, u2, d2):
    # One family and its dedicated Higgs pair. In this block g=1,
    # gamma denotes 16 pi^2 times the one-loop anomalous dimension.
    gamma = (3*u2+2*d2-F(36,5),4*d2-F(24,5),3*u2-F(24,5),4*d2-F(24,5))
    k5 = sum(x*i for x,i in zip(xs,TWICE_INDEX))
    trace = sum(dim*x*gamma_i for dim,x,gamma_i in zip(DIMS,xs,gamma))
    t,f,h,hbar=xs
    cubic = 15*(2*t+h)*u2+20*(t+f+hbar)*d2
    assert trace==cubic-24*k5
    if cubic==0:
        trace_norm = sum(dim*abs(x) for dim,x in zip(DIMS,xs))
        assert 24*abs(k5) <= trace_norm*max(abs(z) for z in gamma)
    return {"K5":k5, "gamma":list(map(str,gamma)), "trace_X_gamma":str(trace),
            "weighted_cubic_breaking":str(cubic)}


def sm_index(c,w,y):
    return (w*int(c==3), c*int(w==2), c*w*y*y)


def add(*arrays):
    return tuple(map(sum,zip(*arrays)))


def scaled(x,a):
    return tuple(x*t for t in a)


def audit():
    # 10 = Q+u^c+e^c, 5 = colour triplet + weak doublet.
    q = sm_index(3,2,1)
    uc,ec = sm_index(3,1,-4),sm_index(1,1,6)
    five = add(sm_index(3,1,-2),sm_index(1,2,3))
    ten = add(q,uc,ec)
    assert five==(1,1,30) and ten==(3,3,90)
    phase_line=scaled(-1,q)
    other_members=add(scaled(-1,uc),scaled(-1,ec),scaled(3,five))
    assert add(other_members,phase_line)==(0,0,0)
    u2,d2=F(8,5),F(6,5)
    exact=block((-1,0,2,1),u2,d2)
    broken=block((1,0,0,0),u2,d2)
    nonfinite=block((-1,0,2,-1),u2,F(0))
    assert exact["K5"]==0 and broken["weighted_cubic_breaking"]=="72"
    assert nonfinite["K5"]==-2 and nonfinite["trace_X_gamma"]=="48"
    for x in ((-1,4,2,-3),(0,2,0,-2),(1,-1,0,2),(1,0,-2,-1)):
        for a,d in ((u2,d2),(F(11,7),F(8,3)),(F(0),F(2))):
            block(x,a,d)
    # Full three-family, four Higgs pairs and the adjoint: gaugino included.
    total_index=3*(F(3,2)+F(1,2))+4*(F(1,2)+F(1,2))+5
    b0=3*5-total_index
    r_anomaly=2*5+2*(F(2,3)-1)*total_index
    assert r_anomaly==F(2,3)*b0==0
    mid=(F(8,3)+F(5,3))/2
    qcd_hi=F(2035,1000)
    return {"normalized_current_on_family_and_higgs":["10","bar5","5_H","bar5_H"],
            "unified_5_indices":five,"unified_10_indices":ten,"Q1_only":phase_line,
            "missing_ordinary_and_higgs_index_in_example":other_members,
            "finite_exact_charge_example":exact,"finite_explicit_breaking_example":broken,
            "symmetry_preserved_but_down_vertices_removed":nonfinite,
            "necessary_max_gamma_for_last_example":str(F(48,25)),
            "b0_finite_assignment":str(b0),"total_gaugino_plus_matter_K5_canonical_R":str(r_anomaly),
            "unsigned_common_QCD_mirror_offset":str(mid),
            "additional_shift_above_old_QCD_max_needed_for_exact_mirror":str(mid-qcd_hi)}


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
