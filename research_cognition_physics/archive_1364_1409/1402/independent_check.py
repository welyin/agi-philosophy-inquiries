"""Three-dimensional common positive factors, determinant pencil, and compact resource bound."""
import argparse
import json
from math import cos, pi, sin, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def eig_trace_test(dmat, emat):
    # det(zE-D)/det(E) from four sample evaluations, without diagonalizing D/E.
    grid = np.array([1.,-1.,1j,-1j])
    values = [np.linalg.det(z*emat-dmat)/np.linalg.det(emat) for z in grid]
    cubic = np.linalg.solve(np.vander(grid,4), values)
    np.testing.assert_allclose(cubic[0],1,atol=1e-12)
    a,b,d = -cubic[1],cubic[2],-cubic[3]
    t = a/np.exp(1j*np.angle(d)/3)
    disc = abs(t)**4-8*(t**3).real+18*abs(t)**2-27
    return {"det_modulus":float(abs(d)), "self_inverse_residual":float(abs(b-d*a.conjugate())),
            "normalized_discriminant":float(disc),
            "circle_prerequisites":bool(abs(abs(d)-1)<1e-10 and abs(b-d*a.conjugate())<1e-10)}


def causal_pair(q, scale):
    base = np.diag([sqrt(q),1/sqrt(q),1.])
    u = np.array([[cos(pi/6),1j*sin(pi/6),0],[1j*sin(pi/6),cos(pi/6),0],[0,0,1.]])
    matrix = np.linalg.solve(base.T,(base @ u).T).T
    emat = scale*np.eye(3)
    dmat = scale*matrix
    three = np.linalg.matrix_power(matrix,3)
    attained = np.linalg.norm(three,2)
    assert np.isclose(attained,q)
    assert np.isclose(np.linalg.cond(base),q)
    f = base/sqrt(q)
    b = scale*np.linalg.inv(f)
    np.testing.assert_allclose(f @ b,emat,atol=1e-14)
    np.testing.assert_allclose(f @ u @ b,dmat,atol=1e-14)
    # Perturb each full normalized matrix independently in operator norm.
    pair_error = .001
    ratio_norm = np.linalg.norm(matrix,2)
    err = (1+ratio_norm)*pair_error/(scale-pair_error)
    lower= (scale-pair_error)*(attained-err*(3*ratio_norm**2+3*ratio_norm*err+err**2))
    return {"q":q,"a":scale,"individual_D_norm":float(np.linalg.norm(dmat,2)),
            "individual_E_norm":scale, "norm_Q_cubed":float(attained),
            "condition_of_F_witness":float(np.linalg.cond(f)),
            "B_norm_with_max_F_one":float(np.linalg.norm(b,2)),
            "common_contractions_possible":bool(scale*q<=1),
            "both_operator_norm_error":pair_error,
            "perturbed_min_sigma_E_times_norm_Q3_lower":float(lower),
            "pencil":eig_trace_test(dmat,emat)}


def calculate():
    good=causal_pair(2.,.5)
    bad=causal_pair(2.,.6)
    assert good["common_contractions_possible"] and not bad["common_contractions_possible"]
    assert bad["individual_D_norm"]<1 and bad["perturbed_min_sigma_E_times_norm_Q3_lower"]>1
    assert good["pencil"]["normalized_discriminant"]<0 and bad["pencil"]["circle_prerequisites"]
    spec=np.diag([.8,.4,.2])
    orient=np.diag([.2,.4,.8])
    identical=eig_trace_test(spec,spec)
    swap=eig_trace_test(spec,orient)
    assert np.allclose(np.linalg.svd(spec,compute_uv=False),np.linalg.svd(orient,compute_uv=False))
    assert swap["circle_prerequisites"] and swap["normalized_discriminant"]>0
    # Exactly the same polynomial and root multiplicities in regular and defective pencils.
    jordan=np.diag([1.,1.,-1.])
    ok=jordan.copy()
    jordan[0,1]=.3
    ej=.5*np.eye(3)
    j1,j2=eig_trace_test(.5*ok,ej),eig_trace_test(.5*jordan,ej)
    assert np.isclose(j1["normalized_discriminant"],0,atol=1e-10)
    assert np.isclose(j2["normalized_discriminant"],0,atol=1e-10)
    dims=[3-int(np.linalg.matrix_rank(.5*m-ej,tol=1e-13)) for m in (ok,jordan)]
    assert dims==[2,1]
    # |det Q|=1 can also fail the self-inversive coefficient condition.
    unequal=np.diag([2.,3.,1/6])
    nonself=eig_trace_test(unequal,np.eye(3))
    assert not nonself["circle_prerequisites"]
    # Nontrivial determinant phase and a complex pair-independent SU(3) trace.
    phases=np.array([.22,.74,1.32])
    complex_example=eig_trace_test(np.diag(np.exp(1j*phases)),np.eye(3))
    assert complex_example["circle_prerequisites"] and complex_example["normalized_discriminant"]<0
    # Recover free orientation for any two equal-volume spectra by a cyclic permutation.
    dspec=np.diag([.9,.2,.16])
    espec=np.diag([.8,.4,.09])
    perm=np.eye(3)[[1,2,0]]
    cycle=dspec @ np.linalg.inv(espec @ perm)
    np.testing.assert_allclose(np.linalg.matrix_power(cycle,3),np.eye(3),atol=1e-13)
    # Reuse independent old low data, and independently transcribe the published
    # charged-lepton Table 2 column for 10^16 GeV (same 2024-data SM trajectory).
    low=json.loads((HERE/"../1386/inputs.json").read_text(encoding="utf-8"))["yukawa"]
    high=json.loads((HERE/"../1397/inputs.json").read_text(encoding="utf-8"))["high_scale"]
    charged_high=[[2.6935e-6,.0091e-6],[.00056745,.00000192],[.009639,.000033]]
    masses={}
    for label,down,leptons in [("MZ",low["down"],low["charged_leptons"]),
                                ("1e16_GeV",high["down"],charged_high)]:
        dvals,derr=np.array(down).T
        evals,eerr=np.array(leptons).T
        rad=np.exp(np.log(dvals).sum()-np.log(evals).sum())
        tau_bounds=[float(np.prod(dvals-i*derr)/np.prod((evals+i*eerr)[:2])) for i in (2,-2)]
        assert (tau_bounds[0]>evals[-1]+2*eerr[-1] if label=="MZ" else
                tau_bounds[1]<evals[-1]-2*eerr[-1])
        masses[label]={"volume_ratio":float(rad),"central_strength_ratio":float(rad**(1/3)),
                       "two_width_predicted_tau_for_equal_strengths":tau_bounds}
    return {"common_gain_sharp_threshold":[good,bad],
            "physical_input_endpoint_cross_checks":masses,
            "same_singular_spectra_free_and_fixed_orientation":{"same":identical,"exchanged":swap},
            "same_polynomial_semisimplicity":{"regular":j1,"defective":j2,"double_root_nullities":dims},
            "equal_determinant_only":nonself,"complex_determinant_case":complex_example,
            "free_cyclic_identification":eig_trace_test(dspec,espec @ perm)}


def verify(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for k in x:
            verify(x[k],y[k])
    elif isinstance(x,list):
        assert len(x)==len(y)
        for a,b in zip(x,y):
            verify(a,b)
    else:
        np.testing.assert_allclose(x,y,atol=2e-11,rtol=2e-10)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    out=HERE/"independent_results.json"
    res=calculate()
    if args.write:
        with out.open("w",encoding="utf-8") as stream:
            stream.write(json.dumps(res,indent=2,allow_nan=False)+"\n")
    else:
        verify(res,json.loads(out.read_text(encoding="utf-8")))
        print("1402 pencil, repeated roots and sharp common contraction threshold verified")
