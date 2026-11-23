"""Role reassignment in the two remaining SU(5) antifive equations."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def grams(D,L,F):
    matter=np.einsum("iAa,iBa->AB", D.conj(),D)+L.conj().T@L
    higgs=np.einsum("iAa,iAb->ab",D.conj(),D)+F.conj().T@F
    return matter-np.eye(len(matter)),higgs-np.eye(len(higgs))


def norm(M):
    return float(np.linalg.norm(M,ord=2))


def saturated(leak):
    # N=3,r=1,h=3. Disjoint cyclic entries keep the three Higgs columns orthogonal.
    D=np.zeros((3,4,3))
    for i in range(3):
        D[i,i,i]=math.sqrt(1-leak)
        D[(i+1)%3,3,i]=math.sqrt(2*leak)
    L=np.array([[0,0,0,math.sqrt(1-6*leak)]])
    F=np.zeros((3,3))
    M,H=grams(D,L,F)
    spent=np.linalg.norm(D[:,3,:])**2
    np.testing.assert_allclose(spent,6*leak,atol=1e-14)
    np.testing.assert_allclose(max(norm(M),norm(H)),leak,atol=1e-14)
    return D,L,F


def two_partners():
    D=np.zeros((3,5,2))
    for i in range(3):
        angle=2*math.pi*i/3
        D[i,i,:]=np.sqrt(.8)*np.array([math.cos(angle),math.sin(angle)])
    L=np.eye(5)[-2:]
    M,H=grams(D,L,np.zeros((2,2)))
    np.testing.assert_allclose(np.diag(M),[-.2,-.2,-.2,0,0],atol=2e-15)
    np.testing.assert_allclose(H,.2*np.eye(2),atol=2e-15)
    return max(norm(M),norm(H))


def traces_only():
    # Embed the actual parent singular spectrum chosen in 1515, without imposing
    # matrix cancellation on the two antifives. Choose the active Higgs to be slot 0.
    vals=np.array(json.loads((HERE/"../1402/inputs.json").read_text(encoding="utf8"))["scales"]["1e16"])
    d,e=np.split(vals[:,0],2)
    s=np.r_[np.sqrt(np.maximum(d,e)[:2]*np.minimum(d,e)[1:]),1.2*max(d[2],e[2])]
    g=.7
    d0=math.sqrt(6/5)*g
    Y=np.column_stack((np.zeros(3),np.diag(s)))
    D=np.zeros((3,4,3))
    D[:,:,0]=Y/d0
    D[0,0,1]=math.sqrt(3-np.linalg.norm(D[:,:,0])**2)
    L=np.array([[0.,1.,1.,1.]])/math.sqrt(3)
    M,H=grams(D,L,np.zeros((3,3)))
    np.testing.assert_allclose([np.trace(M),np.trace(H),(L@L.T)[0,0]-1],0.,atol=2e-15)
    assert np.linalg.norm(Y@L.T)>0
    return {"g_example":g,"parent_singular_values":s.tolist(),
            "normalized_tensor_frobenius_squared":float(np.linalg.norm(D)**2),
            "three_normalized_traces":[float(np.trace(M)),float(np.trace(H)),float((L@L.T)[0,0]-1)],
            "max_antifive_matrix_residual":max(norm(M),norm(H))}


def calculate():
    D,L,F=saturated(0)
    Y=np.einsum("iAa,a->iA",D,np.array([1,2j,3])/math.sqrt(14))
    ell=L.conj().T[:,0]
    np.testing.assert_allclose(Y@ell,0,atol=1e-15)
    m=np.array([.7,1.1,1.3,1j])
    sig=[]
    for q in [2,-3]:
        h=m+q*ell
        h=h/np.linalg.norm(h)
        G=Y@(np.eye(4)-np.outer(h,h.conj()))@Y.conj().T
        sig.append(np.sqrt(np.linalg.eigvalsh(G)))
    assert np.all(sig[0]<=sig[1]+1e-14) or np.all(sig[1]<=sig[0]+1e-14)
    # A genuine difference of ordered spectra: d>e and s<mu in every input corner.
    tab=np.array(json.loads((HERE/"../1402/inputs.json").read_text(encoding="utf8"))["scales"]["1e16"])
    lo,hi=tab[:,0]-2*tab[:,1],tab[:,0]+2*tab[:,1]
    assert lo[0]>hi[3] and hi[1]<lo[4]
    cases={}
    for n in range(1,5):
        feasible=[]
        for r in range(11):
            h=10-2*n-r
            if r>=0 and h>=n:
                feasible.append([r,h])
        cases[str(n)]=feasible
    for eps in (0.0001,.01,.1):
        saturated(eps)
    return {"not_excluded_by_counting_r_h":cases,
            "N3_r1_kernel_read_by_Y":float(np.linalg.norm(Y@ell)),
            "N3_r1_common_direction_spectra_example":[x.tolist() for x in sig],
            "mass_signs_over_two_widths":{"yd_lower":lo[0],"ye_upper":hi[3],"ys_upper":hi[1],"ymu_lower":lo[4]},
            "N3_r2_necessary_normalized_gamma_error_and_witness":two_partners(),
            "r1_max_squared_leakage_per_equal_small_gamma_error_sharp":6,
            "three_trace_alternative":traces_only()}


def check(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            check(a[key],b[key])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=2e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=calculate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        check(result,json.loads(path.read_text(encoding="utf8")))
    print("Verified family count, saturated kernel, sign obstruction and sharp finite tolerances.")
