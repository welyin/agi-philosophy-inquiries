"""Apply the inverse test to the physical angular witness and recover its heavy charged couplings."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def observable_gram(y,f,charged_daughter):
    # Charged Yukawa has no independent direct coupling to the two vectorlike copies.
    _,ell,vh=np.linalg.svd(charged_daughter)
    right=vh.conj().T[:,::-1]
    gram=right.conj().T @ y.conj().T @ (np.eye(3)-f@f) @ y @ right
    return gram/np.outer(ell[::-1],ell[::-1])


def solve():
    m=runpy.run_path(str(HERE/"check.py"))
    run=runpy.run_path(str(HERE/m["CFG"]["physical_mother"]))
    b=run["mother"](.125)
    pt=b["pt"]
    leg=np.array([b["cW"],b["fW"],1.])
    archived=json.loads((HERE/"independent_results.json").read_bytes())["fine"]
    approx=np.array(archived["F_real"])+1j*np.array(archived["F_imag"])
    _,u=np.linalg.eigh(approx)
    f=(u*leg)@u.conj().T
    assert np.linalg.norm(f-approx,2)<1e-9
    y=b["yW"][2,:,:3]
    t=np.linalg.inv(b["S"]@pt["low_E"])
    n=t.T@(pt["C"]+b["nu"]["FROZEN"][0]*np.outer(pt["b"],pt["b"]))@t
    e0,n0,le=m["charged_frame"](y,n)
    e1,n1,_=m["charged_frame"](f@y,f.T@n@f)
    v,ratio=m["visible"](n1)
    inv=m["inverse_case"](e0,n0,e1,v,ratio,leg,recover=True)
    actual=abs(observable_gram(y,f,f@y))
    named=abs(observable_gram(y,b["F"],b["F"]@y))
    gj=np.linalg.eigh(b["F"])[1]
    accepted=[]
    rejected=[]
    for branch in inv["phase_sign_branches"]:
        for trial in branch["trials"]:
            if trial["match"]:
                root=le@trial["filter_src"]@le.conj().T
                err=np.linalg.norm(root-f,2)
                gram=abs(trial["heavy_gram"])
                rephasing_error=np.max(abs(gram-actual))
                evec=np.linalg.eigh(root)[1]
                accepted.append(dict(error_to_injected_root=err,
                                     max_absolute_heavy_gram_difference=rephasing_error,
                                     sorted_old_axes_in_new_eigenlines=(abs(gj.conj().T@evec)**2).tolist(),
                                     normalized_abs_heavy_Higgs_Gram=gram.tolist()))
            else:
                rejected.append(trial["max_singular_error"])
    assert accepted and min(row["error_to_injected_root"] for row in accepted)<3e-7
    assert min(row["max_absolute_heavy_gram_difference"] for row in accepted)<5e-7
    np.testing.assert_allclose(np.linalg.eigvalsh(observable_gram(y,f,f@y)),
                               np.sort(leg**-2-1),atol=1e-13)
    return dict(scale_GeV=b["Mweak"], matched=accepted,
                reference_GJ_abs_heavy_Higgs_Gram=named.tolist(),
                direct_from_fitted_witness_abs_heavy_Gram=actual.tolist(),
                singular_failures_on_trace_crossings=rejected,
                branch_min_trace_excess=[r["min_trace_minus_required"] for r in inv["phase_sign_branches"]])


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    else:
        np.testing.assert_allclose(a,b,rtol=2e-7,atol=2e-8)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=solve()
    out=HERE/"witness_inverse_results.json"
    if args.check:
        compare(result,json.loads(out.read_bytes()))
        print("1746 reconstructed charged heavy Gram and physical spectral witness verified")
    else:
        with out.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
