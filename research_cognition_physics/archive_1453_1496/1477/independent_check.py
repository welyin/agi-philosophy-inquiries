"""Check the signed-response identity, including an arbitrary CP-odd tau coupling."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def window(n):
    # Direct transcription of the sample inputs, independent of the author's verifier.
    rho = ((80.3602-n*.0099)/(91.1876+n*.0021))**2/(1-.23097+n*.00027)/1.02
    amin = math.sqrt((1.77-n*.46)/((2.36+n*.05)*1.05))
    amax = math.sqrt((1.77+n*.59)/((2.36-n*.05)*.95))
    b = math.sqrt((.95-n*.11)/1.05)
    required_negative_w = 2*b+1
    assert required_negative_w > amax
    scalar_min = (2*rho*b-1)/(2*rho-1)
    assert 0 < scalar_min < amax
    alpha_max = math.degrees(math.acos(scalar_min/amax))
    pure_odd_transfer = (.95-n*.11)*(2*rho)**2
    p_abs_max = math.sqrt(2*rho*(1-rho)/(2*rho-1))
    # Unit scalar direction and Goldstone-orthogonal pseudoscalar at the loose
    # beta=1 angular boundary, including a spectator singlet for absolute strength.
    d,c = math.sqrt(2*rho-1),math.sqrt((1-rho)/2)
    p = math.sqrt(amax*amax - scalar_min**2)
    vec = np.array([d*scalar_min, (1-d*d*scalar_min)/(4*c), d*p,-d*d*p/(2*c)])
    scale = .5/np.linalg.norm(vec)
    h = np.append(scale*vec, math.sqrt(.75))
    assert abs(h@h-1) < 1e-12 and abs(d*h[2]+2*c*h[3]) < 1e-12
    cZ = math.sqrt(rho)*(d*h[0]+4*c*h[1])
    cW = (d*h[0]+2*c*h[1])/math.sqrt(rho)
    tau = math.sqrt(rho)*h[[0,2]]/d
    assert np.allclose([cW/cZ, *list(tau/cZ)],[b, scalar_min,p])
    assert abs(tau[1]) <= p_abs_max
    return dict(rho_lower=rho,tau_amplitude_radius=[amin,amax],min_abs_W=b,
                minimum_negative_W_required_abs_scalar=required_negative_w,
                min_tau_CP_even_over_Z=scalar_min,angle_bound_degrees_for_beta_at_most_one=alpha_max,
                required_W_rate_transfer_for_pure_CP_odd=pure_odd_transfer,
                absolute_odd_tau_kappa_bound=p_abs_max,
                normalized_direction_singlet_example=h.tolist(),example_absolute_kZ=cZ)


def run():
    ans={}
    for n in (1,2):
        ans[f"Run1_tau_Run2_W_{n}_width"]=window(n)
    vev=np.array([.72,.31,.25])
    h=np.array([-.33,.47,.18])
    # Finite differences of log pole mass squares and log tau mass.
    def squared(d):
        return np.array([d[0]**2+2*d[1]**2+4*d[2]**2,d[0]**2+4*d[1]**2,d[0]**2])
    eps=.001
    bare=squared(vev)
    deriv=(squared(vev+eps*h)-squared(vev-eps*h))/(2*eps)/bare
    ratio,tau=deriv[[0,2]]/deriv[1]
    rho=bare[0]/bare[1]
    lhs=2*rho*ratio-1-(2*rho-1)*tau
    rhs=8*vev[2]*(h[2]-vev[2]*h[0]/vev[0])/(vev[0]*h[0]+4*vev[1]*h[1])
    assert math.isclose(lhs,rhs,rel_tol=1e-10)
    ans["finite_difference_check_three_representations"]=[lhs,rhs]
    ans["independent_ATLAS_tau_polarimetry_mod_180_deg_68"]=[-7,25]
    return ans


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,abs_tol=2e-10,rel_tol=2e-10)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    choices=cli.add_mutually_exclusive_group()
    choices.add_argument("--write",action="store_true")
    choices.add_argument("--check",action="store_true")
    opt=cli.parse_args()
    path=ROOT/"independent_results.json"
    ans=run()
    if opt.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(ans,f,indent=2)
            f.write("\n")
    else:
        assert equal(ans,json.loads(path.read_text(encoding="utf-8")))
        print("1477 independent: signed mass derivatives and CP-generalized threshold checked.")
