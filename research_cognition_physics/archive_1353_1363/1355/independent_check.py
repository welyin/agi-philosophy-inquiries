"""Round 1355: direct forward and reverse coherently selected preparations."""
import argparse
import json
from pathlib import Path
import numpy as np


def make_routes(name, theta):
    mats=[]
    for j,k in [(0,1),(1,2),(2,0)]:
        e=np.zeros((3,3));e[j,k]=1
        if name=="SU3":
            mats.extend([(e+e.T)/2,(e-e.T)*.5j])
        else:
            mats.append(1j*(e-e.T))
    if name=="SU3":
        mats.extend([np.diag([1.,-1.,0.])/2,np.diag([1.,1.,-2.])/np.sqrt(12)])
    c=sum(np.trace(T@T).real for T in mats)/3
    off=np.concatenate([(np.kron(np.eye(3),T)+np.kron(T,np.eye(3)))/np.sqrt(c) for T in mats])
    H=np.zeros((9+len(off),)*2,complex)
    H[:9,9:]=off.conj().T;H[9:,:9]=off
    e,U=np.linalg.eigh(H)
    prop=(U*np.exp(-1j*theta*e))@U.conj().T
    swap=np.eye(9).reshape(3,3,3,3).transpose(1,0,2,3).reshape(9,9)
    swap_all=np.kron(np.eye(1+len(mats)),swap)
    return [np.eye(len(H)),swap_all,prop]


def data(name):
    routes=make_routes(name,np.pi/np.sqrt(3))
    full_anchor=np.zeros((len(routes[0]),9),complex)
    full_anchor[:9]=np.eye(9)/3
    source0=[r@full_anchor for r in routes]
    X=np.column_stack([raw[:9].ravel() for raw in source0])
    G=X.conj().T@X
    eig,ev=np.linalg.eigh(G)
    B=(ev*(.5/np.sqrt(np.maximum(eig,.25))))@ev.conj().T
    Bfail=(ev*np.sqrt(np.maximum(0,1-(.25/np.maximum(eig,.25)))))@ev.conj().T
    assert np.allclose(B.conj().T@B+Bfail.conj().T@Bfail,np.eye(3))
    dft=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)

    def forward(logical):
        coeff=B@logical
        rail=np.array([z*r for z,r in zip(coeff,source0)])
        measured=np.einsum("ij,jkl->ikl",dft,rail)
        success=measured[0,:9].ravel()
        all_flags=(float(np.linalg.norm(Bfail@logical)**2),
                   float(np.linalg.norm(measured[1:])**2),
                   float(np.linalg.norm(measured[0,9:])**2),
                   float(np.linalg.norm(success)**2))
        assert np.isclose(sum(all_flags),np.linalg.norm(logical)**2)
        return success

    def decode(physical):
        start=np.zeros_like(full_anchor)
        start[:9]=physical.reshape(9,9)
        rails=[U.conj().T@start/np.sqrt(3) for U in routes]
        anchor_flags=np.array([np.vdot(full_anchor,r) for r in rails])
        rejected=sum(np.linalg.norm(r)**2 for r in rails)-np.linalg.norm(anchor_flags)**2
        result=B@anchor_flags
        assert np.isclose(rejected+np.linalg.norm(Bfail@anchor_flags)**2+
                          np.linalg.norm(result)**2,np.linalg.norm(physical)**2)
        return result

    K=np.column_stack([forward(q) for q in np.eye(3)])
    threshold=1e-10
    nonzero=eig>threshold
    support=ev[:,nonzero]@ev[:,nonzero].conj().T
    assert np.min(eig[nonzero])>.25
    assert np.allclose(12*K.conj().T@K,support)
    assert np.allclose(K@(12*K.conj().T@K),K)
    ar=make_routes(name,np.pi/(2*np.sqrt(3)))[2]@full_anchor
    probes=[(X[:,0]+1j*X[:,1])/np.sqrt(2), ar[:9].ravel()/np.linalg.norm(ar[:9])]
    probe=(np.array([1,1j,0]))/np.sqrt(2)
    probability=[]
    label_example=[]
    for v in probes:
        u=decode(v)
        assert np.allclose(u,K.conj().T@v)
        assert np.isclose(np.linalg.norm(u)**2,1/12)
        assert np.allclose(forward(u),v/12)
        probability.append(float(12*abs(np.vdot(probe,u))**2))
        # A missing relative phase calibration, despite unchanged individual route norms.
        rephased=u*np.array([1,-1,1])
        label_example.append(float(12*abs(np.vdot(probe,rephased))**2))
    # An external E with rank-dimensional Schmidt support is never acted upon.
    Q=ev[:,nonzero]
    bell=Q/np.sqrt(sum(nonzero))
    encoded=np.column_stack([forward(col) for col in bell.T])
    recovered=np.column_stack([decode(col) for col in encoded.T])
    assert np.allclose(recovered,bell/12)
    # Leaving a copy of the orthogonal input label in the environment destroys Bell coherences.
    # In the Schmidt basis the dephased Bell state has fidelity 1/r.
    rank=int(sum(nonzero))
    # Local-only counterexample: symmetric vs. antisymmetric retained partner supports.
    swap=X[:,1].reshape(9,9)*3
    ps=(np.eye(9)+swap)/2;pa=(np.eye(9)-swap)/2
    cp,cm=ps/np.sqrt(6),pa/np.sqrt(3)
    assert np.allclose(cp@cm.conj().T,0)
    assert np.allclose(cp.conj().T@cm,0)
    return dict(rank=rank,gram_eigenvalues=eig.tolist(),neutral_readout_given_decode=probability,
                uncontrolled_second_rail_sign=label_example,
                bell_round_trip_unconditional=float(np.linalg.norm(recovered)**2),
                optimal_herald_rate_using_raw_erasure_and_a_single_prefilter=float(min(eig[nonzero])/3),
                raw_erasure_three_basis_weights=(np.diag(G).real/3).tolist(),
                unreturned_which_label_Bell_fidelity=1/rank,
                passive_internal_partners_Bell_fidelity_upper=.5)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    else:
        assert np.allclose(a,b,rtol=1e-11,atol=1e-11),(a,b)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    res={name:data(name) for name in ["SU3","SO3","mu3SO3"]}
    path=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(res,json.loads(path.read_text(encoding="utf8")))
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(res,f,indent=2,allow_nan=False);f.write("\n")
    print("PASS",res["SU3"]["neutral_readout_given_decode"],res["SO3"]["neutral_readout_given_decode"])
