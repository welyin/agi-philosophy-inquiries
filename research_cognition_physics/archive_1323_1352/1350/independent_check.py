"""Independent matrix moment checks and explicit single-clock instruments, round 1350."""
import argparse
import json
import math
from pathlib import Path
import numpy as np


def spectral_fun(A,f):
    eig,vec=np.linalg.eigh(A)
    return (vec*f(eig))@vec.conj().T


def kron_filter(X):
    # Conditional rail tau=1, U=exp(-iH), input control |+>, observed control |+>.
    H=spectral_fun(X,lambda x:2*np.arccos(np.sqrt(np.clip(x,0,1))))
    U=spectral_fun(H,lambda x:np.exp(-1j*x))
    return (np.eye(len(U))+U)/2


def records(X,V,k):
    K=kron_filter(X)
    output=np.linalg.matrix_power(K,k)@V
    return output.conj().T@output


def run():
    eye=np.eye(2)
    # Noncommuting A,B and matrix-valued feasible extensions.
    A=np.array([[.42,.09j],[-.09j,.57]])
    G=np.array([[.14,.03],[.03,.11]])
    K=G@G
    B=A@A+K
    low=B@np.linalg.solve(A,B)
    high=B-(A-B)@np.linalg.solve(eye-A,A-B)
    DD=G@np.linalg.solve(A,G)
    DU=eye-G@np.linalg.solve(eye-A,G)
    V=np.vstack([eye,np.zeros((2,2))])
    out=[]
    for s in [0.,.13,.61,1.]:
        X=np.block([[A,G],[G,(1-s)*DD+s*DU]])
        assert np.linalg.eigvalsh(X).min()>=-1e-14
        assert np.linalg.eigvalsh(X).max()<=1+1e-14
        for n,M in enumerate([eye,A,B,(1-s)*low+s*high]):
            assert np.allclose(records(X,V,n),M,atol=1e-12)
        out.append(np.linalg.eigvalsh(records(X,V,3)).tolist())
    # Full-rank source is not assumed to be an isometry in the physical tomography.
    scale=np.array([[.71,.15],[0.,.39]])
    raw=[records(X,V@scale,k) for k in range(3)]
    invsqrt=spectral_fun(raw[0],lambda x:1/np.sqrt(x))
    U=scale@invsqrt
    assert np.allclose(U.conj().T@U,eye)
    for k,pop in enumerate([eye,A,B]):
        assert np.allclose(invsqrt@raw[k]@invsqrt,U.conj().T@pop@U)
    # M3-A^3 need not be positive semidefinite, despite X positive.
    X=np.array([[.5,.1,.1],[.1,.5,0],[.1,0,.5]])
    V3=np.eye(3)[:,:2]
    M1=V3.T@X@V3
    third_defect=records(X,V3,3)-np.linalg.matrix_power(M1,3)
    negative=np.linalg.eigvalsh(third_defect)
    assert negative[0]<-1e-6
    # A common phase-sensitive pair of weak inputs; all real input readouts agree.
    sy=np.array([[0.,-1j],[1j,0.]])
    targets=[.5*eye+.25*sy,.5*eye-.25*sy]
    real_inputs=[np.array([np.cos(t),np.sin(t)]) for t in np.linspace(0,np.pi,19)]
    phase_input=np.array([1.,1j])/np.sqrt(2)
    coherent=[]
    for X in targets:
        coherent.append(float(np.vdot(phase_input,records(X,eye,3)@phase_input).real))
    for i in range(1,6):
        M,N=[records(X,eye,i) for X in targets]
        assert all(np.isclose(np.vdot(r,M@r),np.vdot(r,N@r)) for r in real_inputs)
    assert np.allclose(coherent,[27/64,1/64])
    # Exact effect moments vs. a successful state update escaping the source subspace.
    E=np.pi*np.array([1/3,5/3,1/2])
    Vmix=np.array([[1.,0],[1.,0],[0.,np.sqrt(2)]])/np.sqrt(2)
    X=np.diag(np.cos(E/2)**2)
    S=np.diag((1+np.exp(-1j*E))/2)
    A0=Vmix.T@X@Vmix
    assert np.allclose(Vmix.T@X@X@Vmix,A0@A0)
    success_output=S@Vmix[:,0]
    leak=float(np.linalg.norm((np.eye(3)-Vmix@Vmix.T)@success_output)**2/np.linalg.norm(success_output)**2)
    assert np.isclose(leak,.25)
    # Two irrationally related real-filter clocks do not even distinguish these energies.
    alias=np.pi*np.array([np.sqrt(2)+1,np.sqrt(2)-1])
    alias_vals=np.cos(np.outer([1,np.sqrt(2),np.sqrt(3)],alias)/2)**2
    assert np.allclose(alias_vals[:2,0],alias_vals[:2,1])
    assert abs(alias_vals[2,0]-alias_vals[2,1])>.1
    # Best uniform fixed-k bound, and its long-time limitation.
    sharp={}
    for k in [3,4,7]:
        h=.001
        err=(1+(1-2*h)**k)/2-(1-h)**k
        ratio=err/h**2
        assert .99*math.comb(k,2)<ratio<=math.comb(k,2)+1e-9
        sharp[str(k)]=ratio
    return dict(ordered_third_extrema_interpolation=out,
                third_interval_width_opnorm=float(np.linalg.norm(high-low,2)),
                nonpositive_third_defect=negative.tolist(),
                coherent_triple_pass=coherent, normalized_postselection_escape=leak,
                two_clock_alias_and_third_clock=alias_vals.tolist(),
                near_sharp_normalized_errors=sharp)


def compare(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for k in x:
            compare(x[k],y[k])
    else:
        assert np.allclose(x,y,rtol=2e-9,atol=2e-11),(x,y)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    dst=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(result,json.loads(dst.read_text(encoding="utf-8")))
    else:
        with dst.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False);f.write("\n")
    print("PASS",result["third_interval_width_opnorm"],result["coherent_triple_pass"])
