"""Exact multiplicities tested with the full SU(15) kinetic Gram matrix."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(filename):
    return json.loads((HERE/filename).read_text(encoding="utf-8-sig"))


def generators(dim):
    # Real orthonormal basis of traceless Hermitian matrices, tr(Ta Tb)=delta_ab.
    series = []
    for i in range(dim):
        for j in range(i):
            t = np.zeros((dim,dim),dtype=complex)
            t[i,j] = t[j,i] = 1/math.sqrt(2)
            series.append(t)
            u = np.zeros_like(t)
            u[i,j], u[j,i] = 1j/math.sqrt(2),-1j/math.sqrt(2)
            series.append(u)
    for k in range(1,dim):
        t = np.zeros((dim,dim),dtype=complex)
        t[np.arange(k),np.arange(k)] = 1
        t[k,k] = -k
        series.append(t/math.sqrt(k*(k+1)))
    return np.array(series)


def vacuum(n,a,b1,b2):
    diag = np.diag(np.repeat([-a,0,a],n))
    edge = np.array([[0,b1,0],[b1,0,b2],[0,b2,0]])
    return diag, np.kron(edge,np.eye(n))


def gram_modes(T,A,S):
    dA = T@A-A@T
    dS = T@S+S@T.transpose(0,2,1)
    deriv = np.concatenate([dA.reshape(len(T),-1),dS.reshape(len(T),-1)],axis=1)
    # Gauge fields have real coefficients in this Hermitian basis.
    return np.linalg.eigvalsh((deriv@deriv.conj().T).real)


def secular(trace,det):
    largest = (trace+math.sqrt(trace**2-4*det))/2
    return [det/largest,largest]


def analytic(n,a,b):
    aa,bb = a*a,b*b
    lam = secular(4*aa+8*bb,24*aa*bb)
    cond = secular(4*aa+8*bb/3,8*aa*bb/3)
    return [(0,n*n-1),(2*bb,n*n),(4*aa+2*bb,n*n),
            (aa,n*(n-1)),(aa+2*bb,2*n*n),(aa+8*bb,n*(n+1)),
            (lam[0],n*n-1),(lam[1],n*n-1),(cond[0],1),(cond[1],1)]


def potential(A,S,a,b):
    eye=np.eye(len(A))
    x=A/a
    y=S/b
    pm=(x@x-x)/2
    pp=(x@x+x)/2
    p0=eye-x@x
    bm,bp=pm@y@p0.T,pp@y@p0.T
    constraints=[x@x@x-x,y-bm-bm.T-bp-bp.T]
    for p,block in ((pm,bm),(pp,bp)):
        constraints.extend((block@block.conj().T-p,block.conj().T@block-p0.conj()))
    return sum(np.linalg.norm(m)**2 for m in constraints)


def values():
    cfg=load("inputs.json")
    n=cfg["block_n"]
    dim=cfg["N"]
    assert dim==3*n
    T=generators(dim)
    rows=[]
    for sample in cfg["vacuum_tests"]:
        a,b=sample["a"],sample["b"]
        A,S=vacuum(n,a,b,b)
        spectrum=analytic(n,a,b)
        modes=gram_modes(T,A,S)
        full=np.sort(np.concatenate([np.full(count,mass) for mass,count in spectrum]))
        assert len(modes)==len(full)==dim**2-1
        assert np.allclose(modes,full,atol=cfg["relative_numeric_tolerance"],rtol=0)
        assert abs(potential(A,S,a,b))<1e-20
        rows.append({"a":a,"b":b,"analytic_mass_squared_and_counts":
                     [[round(m,10),k] for m,k in sorted(spectrum)],
                     "numeric_kernel":int(np.count_nonzero(abs(modes)<1e-9)),
                     "full_gram_dimension":len(modes),
                     "relative_squared_mass_50_over_20":round((a*a+2*b*b)/(a*a),10)})
    A,S=vacuum(n,1,1,0)
    reduced=gram_modes(T,A,S)
    missing=int(np.count_nonzero(abs(reduced)<1e-9))
    assert missing==2*n*n-1 and potential(A,S,1,1)>0
    two_A=np.diag(np.repeat([-1,1],n))
    two_S=np.kron(np.array([[0,1],[1,0]]),np.eye(n))
    two=gram_modes(generators(2*n),two_A,two_S)
    two_count=int(np.count_nonzero(abs(two)<1e-9))
    assert two_count==n*n
    # Test a complex unitary change of variables, both right actions distinct.
    rng=np.random.default_rng(1577)
    unitary,_=np.linalg.qr(rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim)))
    unitary[:,0]/=np.linalg.det(unitary)
    A,S=vacuum(n,1,1,1)
    Anew,Snew=unitary@A@unitary.conj().T,unitary@S@unitary.T
    assert potential(Anew,Snew,1,1)<1e-20
    return {"three_scale_points":rows,"two_block_kernel":two_count,
            "three_blocks_one_missing_edge_kernel":missing,"complex_gauge_covariance":True,
            "zero_mode_coupling_g5_over_g15":round(1/math.sqrt(3),10)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    flags=parser.parse_args()
    result=values()
    if flags.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as output:
            json.dump(result,output,indent=2,allow_nan=False)
            output.write("\n")
    elif result!=load("results.json"):
        raise SystemExit("Saved mass spectrum differs")
    print(json.dumps(result))
