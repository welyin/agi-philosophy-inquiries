"""Inverse lifetime bound from a positive charged-current loop weight."""
import argparse
import json
from math import sqrt,pi
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
MW=80.4
ALPHA=1/137.036
BR_LIMIT=1.5e-13
T_MIN=(1.0018-.0014)**2
TINY_MASS_GEV=.1e-9


def loop(x):
    z,w=np.polynomial.legendre.leggauss(48)
    z=(z+1)/2
    return .75*x*np.dot(w,z*(1-z)**2/(1-z+x*z))


def min_overlap(tau_squared):
    return (tau_squared-1)/(2*(2-sqrt(4-tau_squared)))


def work(M):
    gl=loop((M/MW)**2)
    # sqrt(det X_plane)*|n_tau| <= m_max/(sqrt(6)*M)
    loss=(.5+gl)*TINY_MASS_GEV/(2*sqrt(6)*M) + (TINY_MASS_GEV/MW)**2/4
    return gl,loss


def bound(tau_squared,M):
    gl,loss=work(M)
    amp=max(0.,gl*min_overlap(tau_squared)-loss)
    return 3*ALPHA/(2*pi)*amp*amp


def bisect(predicate,lo,hi):
    assert not predicate(lo) and predicate(hi)
    for _ in range(52):
        mid=(lo+hi)/2
        if predicate(mid):
            hi=mid
        else:
            lo=mid
    return [lo,hi]


def calculate():
    # For x>=1 monotonicity is analytic. Check the integral against the closed form.
    for x in (.2,.4,1.55,2.,5.):
        expected=-x*(2*x*x+5*x-1)/(4*(1-x)**3)-3*x**3*np.log(x)/(2*(1-x)**4)
        assert np.isclose(expected,loop(x),rtol=1e-11)
    assert np.isclose(loop(1),.125)
    a=min_overlap(T_MIN)
    r=2-sqrt(4-T_MIN)
    assert np.isclose(a,(T_MIN-1)*(2-r)/(T_MIN-r*r))
    mcrit=bisect(lambda m: bound(T_MIN,m)>BR_LIMIT,4,60)
    gcrit=bisect(lambda g: bound(g*g,100)>BR_LIMIT,1,1.001)
    # Rigorous rounded 20 GeV sufficient cutoff, without quadrature for G:
    # x/(8 max(1,x)) lower; bound all subtractions with gl<=1/2.
    x=(20/MW)**2
    amp20=x/8*a-TINY_MASS_GEV/(2*sqrt(6)*20)-(TINY_MASS_GEV/MW)**2/4
    br20=3*ALPHA/(2*pi)*amp20**2
    assert amp20>0 and br20>BR_LIMIT
    return dict(minimum_overlap_norm=a,optimal_tau_to_e_loss_sqrt=r,
                min_loop_100GeV=work(100)[0],branch_floor_100GeV=bound(T_MIN,100),
                conservative_branch_floor_20GeV=br20,MEGII_limit=BR_LIMIT,
                critical_lower_heavy_mass_GeV=mcrit,upper_tau_over_e_at_100GeV=gcrit)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    res=calculate()
    if p.parse_args().check:
        old=json.loads((HERE/"independent_results.json").read_text(encoding="utf8"))
        assert old.keys()==res.keys() and all(np.allclose(v,old[k],atol=1e-25,rtol=1e-10) for k,v in res.items())
        print("1409 sharp row-geometry envelope, inverse MEG bound and conservative mass cutoff verified")
    else:
        print(json.dumps(res,indent=2))
