"""Exact-isospectral two-scalar/D quark family; tree benchmark and charge check."""
import argparse
import json
from math import acos, cos, pi, sin, sqrt
from pathlib import Path
import numpy as np


def mass(phi,delta,rs,rt,iso):
    s,tt=rs*np.exp(1j*phi),rt*np.exp(1j*(phi+delta))
    v=np.array([1.,2.,2.])/3
    o=2*np.outer(v,v)-np.eye(3)
    if iso:
        m=o@np.diag([1.,2.,3.])
        c=o[:,0]/2
        b=[s,tt,s]
    else:
        m=np.array([[1.,0.,0.],[1.,2.,0.],[1.,1.,3.]])
        c=m[:,1]/2
        b=[s,tt,s+tt]
    out=np.column_stack([m,c])
    return np.vstack([out,[*b,4*s]]),o


def read(m):
    h=m@m.conj().T
    values,l=np.linalg.eigh(h)
    cp=(h[0,1]*h[1,3]*h[3,0]).imag
    j=(l[0,0]*l[1,1]*l[0,1].conj()*l[1,0].conj()).imag
    return values,cp,j


def benchmark():
    d=acos(1/3)
    arg=np.angle(4-.5*np.exp(1j*d))
    phi=-arg
    m,_=mass(phi,d,1,1,False)
    vals,cp,j=read(m)
    assert np.isclose(cp,-2*sin(d)) and abs(j)>.0001
    assert abs(np.angle(np.linalg.det(m)))<5e-15
    # V=(r_s²−1)²+(r_t²−1)²+(r_s r_t cosδ−1/3)²+1−cos(theta).
    # At the common zero residual the Hessian is the positive weighted Gram of its gradients.
    z=4-.5*np.exp(1j*d)
    dtheta_s=np.imag(4/z)
    dtheta_t=np.imag(-.5*np.exp(1j*d)/z)
    dtheta_delta=np.imag(-.5j*np.exp(1j*d)/z)
    jac=np.array([[2.,0.,0.,0.],[0.,2.,0.,0.],[cos(d),cos(d),0.,-sin(d)],
                  [dtheta_s,dtheta_t,1.,dtheta_delta]])
    hess=jac.T@np.diag([2.,2.,2.,1.])@jac
    evals=np.linalg.eigvalsh(hess)
    assert min(evals)>0 and np.linalg.matrix_rank(jac)==4
    return {"relative_phase":d,"common_phase":phi,"physical_left_invariant":cp,
            "light_quartet":j,"full_tree_potential_coordinate_hessian_eigenvalues":evals.tolist(),
            "positive_down_squared_masses":vals.tolist()}


def isospectral():
    points=[]
    for rs,rt in [(1.,1.),(1.3,.8)]:
        zero,o=mass(0,0,rs,rt,True)
        for phi,d in [(0,0),(.43,.7),(-.21,-.83),(0,pi/2),(0,-pi/2)]:
            m,_=mass(phi,d,rs,rt,True)
            phase=np.diag([1.,np.exp(1j*d),1.])
            l=np.eye(4,dtype=complex); rr=l.copy()
            l[:3,:3]=o@phase.conj().T@o.T
            l[3,3]=np.exp(1j*phi)
            rr[:3,:3]=phase
            assert np.allclose(m,l@zero@rr)
            assert abs(np.linalg.det(m)-21*rs*np.exp(1j*phi))<4e-13
            v,cp,j=read(m)
            assert np.allclose(v,read(zero)[0])
            assert np.isclose(cp,(790/81)*rs*rt*sin(d))
            if rs==rt==1 and abs(d)==pi/2:
                assert abs(j)>.003
                points.append({"delta":d,"light_quartet":j,"physical_left_invariant":cp,
                               "squared_masses":v.tolist()})
    return {"weak_CP_at_minima_with_identical_pure_QCD":points}


def compute():
    # Primitive convention: N=2 T(3) Δq=1; E=2*3*Q² Δq=2/3.
    return {"primitive_N":1,"E_over_N":2/3,"author_benchmark":benchmark(),
            "fixed_spectrum_competition":isospectral()}


def rounded(x):
    if isinstance(x,dict): return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list): return [rounded(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--check",action="store_true")
    dst=Path(__file__).with_name("independent_results.json")
    out=rounded(compute())
    if ap.parse_args().write:
        with dst.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2,ensure_ascii=False); f.write("\n")
    else: assert out==json.loads(dst.read_text(encoding="utf-8"))
    print("1391 common-anomaly minimum, Hessian and exact QCD isospectral weak-CP family verified.")
