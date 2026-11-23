"""Exact winding counts, angular potential and strictly holomorphic neutral-insertion test."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def json_at(path):
    return json.loads((HERE/path).read_bytes())


def count_windings(N):
    count=[0,0,0]
    for l2,l3,r2,r3 in itertools.product(range(N),repeat=4):
        c,r,e,f=[x%N for x in [l2+r2,l2+r3,l3+r2,l3+r3]]
        w=r+e-c-f
        assert w%N==0
        m=w//N
        assert m in (-1,0,1)
        count[m+1]+=1
    outer=N*N*(N*N-1)//6
    assert count==[outer,N**4-2*outer,outer]
    return count


def potential_test(phi):
    # Dimensionless single-flavon model: r = sqrt(2)|S|/f, theta = arg S.
    eta=4*phi
    orbit=np.array([(sign*abs(eta)+2*math.pi*k)/4 for sign in (1,-1) for k in range(4)])
    def energy(r,theta):
        return (r*r-1)**2+r**8*(math.cos(4*theta)-math.cos(eta))**2
    assert max(energy(1,t) for t in orbit)<1e-25
    assert len(set(np.round(np.mod(orbit,2*np.pi),9)))==8
    def kernel(r,theta):
        return np.array([r*r-1,r**4*(math.cos(4*theta)-math.cos(eta))])
    step=1e-6
    J=np.stack([(kernel(1+step,phi)-kernel(1-step,phi))/(2*step),
                (kernel(1,phi+step)-kernel(1,phi-step))/(2*step)],axis=1)
    np.testing.assert_allclose(2*J.T@J,np.diag([8,32*math.sin(eta)**2]),atol=4e-8)
    return {"two_discrete_orbits_deg_mod360":np.round(np.degrees(orbit)%360,8).tolist(),
            "dimensionless_radial_angular_hessian_diagonal":[8,32*math.sin(eta)**2]}


def physical_moments(extra,config):
    N=config["N"];eps=config["epsilon"]
    old=json_at("../1489/results.json")["reconstruction"]
    phase=math.radians(old["eta_example_deg"])/N
    xi=config["v_GeV"]*config["sin_theta"]/config["f_GeV"]
    hproj=math.sqrt(1-config["sin_theta"]**2)
    real=np.array(old["normalized_real_edges_A_B_C_E_F"]+old["two_real_coefficients"])
    edges=[(0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2)]
    k=np.array([(config["left"][i]+config["right"][j])%N for i,j in edges])
    k[-1]+=N
    def action(hov):
        x=eps*np.exp(1j*phase)*(1+xi*hov)
        vals=real*(x/eps)**k
        if extra=="nonholomorphic":
            vals[2]*=abs(x/eps)**2
        elif extra=="holomorphic":
            vals[2]*=1-eps**(3*N)+x**(3*N)
        Y=np.zeros((3,3),complex)
        for (i,j),val in zip(edges,vals):
            Y[i,j]+=val
        return (1+hproj*hov)*Y
    # Retain raw weak-basis phases. SVD gives rephasing-independent physical vertices.
    mass=action(0)
    step=1e-5
    dY=(action(step)-action(-step))/(2*step)
    left,ms,vt=np.linalg.svd(mass)
    idx=np.argsort(ms);ms=ms[idx]
    p=np.imag(np.diag(left[:,idx].conj().T@dY@vt.conj().T[:,idx])/ms)
    lam=ms**2
    H=mass@mass.conj().T
    D=real[-1]*np.exp(1j*N*phase)
    norm1=real[-2]*D.imag
    norm2=(H[0,1]*H[1,2]*H[2,0]).imag/H[0,0].real
    T=H[1,1].real+real[4]**2
    out=[float(lam@p/norm1),float((lam*(lam-T))@p/norm2)]
    kappa={"minimal":0,"nonholomorphic":2,"holomorphic":3*N*eps**(3*N)}[extra]
    np.testing.assert_allclose(out,[N*xi,(N-kappa)*xi],atol=3e-8)
    return {"mass":mass,"measured":out,"shift_in_C":kappa,"J_mass_weak_ring":norm2*H[0,0].real}


def solve():
    cfg=json_at("inputs.json")
    assert cfg["N"]==4
    old=json_at("../1489/results.json")["reconstruction"]
    eta=math.radians(old["eta_example_deg"])
    runs={name:physical_moments(name,cfg) for name in ("minimal","nonholomorphic","holomorphic")}
    for run in runs.values():
        np.testing.assert_allclose(run["mass"],runs["minimal"]["mass"],atol=1e-14)
    for run in runs.values():
        del run["mass"]
        run["measured"]=[round(v,8) for v in run["measured"]]
    hist={str(N):count_windings(N) for N in (2,3,4,5,6,7,8,9)}
    return {"winding_exact_counts_minus_zero_plus":hist,"phase_potential":potential_test(eta/4),
            "physical_moments_by_C_insertion":runs}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=solve()
    path=HERE/"independent_results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
        print("1653 independent PASS: exact charge counts, scalar minima, physical moments and two neutral insertions")
    else:
        with path.open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")

