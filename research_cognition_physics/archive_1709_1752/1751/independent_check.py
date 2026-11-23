"""Contractive polynomial alternatives for a middle-family reflection on two actual spectra."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
BG = runpy.run_path(str(HERE/"../1744/check.py"))
MARG = json.loads((HERE/"../1747/independent_results.json").read_bytes())["fine"]


def theta(x):
    return 2*math.asin(math.sqrt(x))


def angles_and_weights(a):
    t=theta(a)
    lower=math.ceil(math.pi/t)
    N=math.ceil(math.pi/(2*t))
    powers=np.array([N-2,N])
    values=np.cos(powers*t)
    weights=np.array([-values[1],values[0]])/(values[0]-values[1])
    assert min(weights)>0 and 2*N-lower in (0,1)
    return lower,powers,weights


def port(x,powers,weights,full=False):
    t=theta(x)
    signal=float(weights@np.cos(t*powers))
    response=2*signal**2-1
    if full:
        # Program register plus two path states; a single shared reflection is
        # applied after the controlled evolutions and before their inverses.
        rot=np.array([[math.cos(t),-math.sin(t)],[math.sin(t),math.cos(t)]])
        out=np.zeros((2*len(powers),)*2)
        for j,n in enumerate(powers):
            out[2*j:2*j+2,2*j:2*j+2]=np.linalg.matrix_power(rot,int(n))
        ref=np.zeros(len(out))
        ref[::2]=np.sqrt(weights)
        actual=out.T@(2*np.outer(ref,ref)-np.eye(len(out)))@(out@ref)
        assert abs(ref@actual-response)<2e-13
        np.testing.assert_allclose(np.linalg.norm(actual),1,atol=1e-13)
        return response,np.linalg.norm(actual-response*ref)
    return response


def separating_certificate(a,b,n):
    # A single real hyperplane with strictly positive readings separates
    # every normalized nonnegative program with these powers from (0,0).
    numbers=np.arange(n%2,n+1,2)
    points=np.array([np.cos(numbers*theta(x)) for x in (a,b)]).T
    pairs=[]
    for x,y in itertools.combinations(points,2):
        vec=y-x
        lam=np.clip(-x@vec/(vec@vec),0,1)
        q=x+lam*vec
        pairs.append(q)
    normal=min(pairs,key=np.linalg.norm)
    normal/=np.linalg.norm(normal)
    margin=np.min(points@normal)
    assert margin>1e-3
    return dict(powers_parity=int(n%2),highest_power=n,normal=normal,margin=margin)


def interval_effect(alpha,window,cap):
    dtheta=max(abs(theta(q)-theta(alpha)) for q in window)
    return dict(a_min_max=window,
                minimum_degree_range_if_error_zero=[math.ceil(math.pi/theta(v)) for v in window[::-1]],
                max_radian_error=dtheta,same_calibration_max_middle_error=min(2.,2*(cap*dtheta)**2))


def solve(step):
    b=BG["mother"](step)
    lepton=b["yW"][2,:,:3]
    down=b["yW"][1,:,:3]@b["right_quark_to_GJ"]
    spectra=[(np.linalg.svd(Y,compute_uv=False)[::-1]/np.linalg.norm(Y,2))**2 for Y in (down,lepton)]
    reference=json.loads((HERE/"../1750/results.json").read_bytes())
    out={}
    a_list=[]
    for name,h in zip(("quark","lepton"),spectra):
        np.testing.assert_allclose(h,reference[name]["normalized_eigenvalues"],atol=5e-13)
        a=(h[1]-h[0])/(h[2]-h[0])
        a_list.append(a)
        lower,powers,weights=angles_and_weights(a)
        cap=int(max(powers))
        audited=[port(x,powers,weights,True) for x in (0.,a,1.)]
        np.testing.assert_allclose(np.array(audited)[:,0],[1.,-1.,1.],atol=1e-14)
        assert max(x[1] for x in audited)<3e-14
        # Under conjugation, the output reflection rotates twice as far as its
        # middle projector. This does not depend on intermediate eigenvalue probes.
        out[name]=dict(squared_eigenvalues_in_heaviest_units=h,middle_fraction=a,
            minimal_polynomial_degree_lower=lower,explicit_degree_upper=2*cap,powers=powers,weights=weights,
            required_degree_if_error_at_most_one_percent=math.ceil((math.pi-2*math.acos(.99))/theta(a)),
            exact_quadratic_at_half=1-.5/(a*(1-a)),contractive_polynomial_at_quarter=port(.25,powers,weights),
            uniform_orbit_Lipschitz_floor=2/a,explicit_midpoint_return_loss=audited[1][1])
    # Common coefficients in each sector's *normalized* min-max scanner coordinate.
    powers=np.array([33,35,53])
    weights=np.linalg.solve(np.array([np.ones(3),*[np.cos(powers*theta(a)) for a in a_list]]),[1.,0.,0.])
    checks=[port(x,powers,weights,True) for x in [0.,*a_list,1.]]
    assert min(weights)>.26
    np.testing.assert_allclose(np.array(checks)[:,0],[1.,-1.,-1.,1.],atol=2e-14)
    assert max(x[1] for x in checks)<5e-14
    probe=.1
    original=port(probe,powers,weights)
    changed=(original+1)*(1-4*probe*(1-probe))**2-1
    out["common"]=dict(powers=powers,weights=weights,degree=2*int(max(powers)),
        exact_bounds_for_subclass_of_nonnegative_same_parity_programs=[separating_certificate(*a_list,n) for n in (51,52)],
        quarter_half_three_quarters=[port(x,powers,weights) for x in (.25,.5,.75)],
        largest_fiducial_remainder=max(x[1] for x in checks),probe_decimal_one_and_alternative=[original,changed])
    # Both charged ratios (relative to tau) may move anywhere inside the envelope
    # used in 1747. This is a deterministic stress test of the *frozen* calibration.
    h=spectra[1]
    lo,hi=np.array(MARG["charged_gain_ratio_intervals"])
    al=(h[1]*lo[1]**2-h[0]*hi[0]**2)/(1-h[0]*hi[0]**2)
    ah=(h[1]*hi[1]**2-h[0]*lo[0]**2)/(1-h[0]*lo[0]**2)
    out["lepton"]["original_charged_ratio_box"]=interval_effect(a_list[1],(al,ah),max(out["lepton"]["powers"]))
    out["common"]["charged_box"]=interval_effect(a_list[1],(al,ah),max(powers))
    return out


def clean(o):
    if isinstance(o,dict):
        return {k:clean(v) for k,v in o.items()}
    if isinstance(o,(tuple,list,np.ndarray)):
        return [clean(v) for v in o]
    if isinstance(o,(int,np.integer)):
        return int(o)
    if isinstance(o,(float,np.floating)):
        return float(f"{o:.10g}")
    return o


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    fine,coarse=solve(.125),solve(.25)
    for key in ("quark","lepton"):
        assert fine[key]["minimal_polynomial_degree_lower"]==coarse[key]["minimal_polynomial_degree_lower"]
    result=clean(dict(fine=fine,largest_step_halving_shift_in_common_prepare=max(abs(fine["common"]["weights"]-coarse["common"]["weights"]))))
    path=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1751 passive spectral reflection bounds, two-sector unitary witness and public row resolution verified")
    else:
        with path.open("w",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")

