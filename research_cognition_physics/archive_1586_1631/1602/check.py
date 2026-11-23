"""Two mass tensors, one CP-even scalar: exact rate identities and a cross-experiment holdout."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def width_test(cfg, n):
    wh, wh_minus, _ = cfg["ATLAS_WHbb"]
    bsm, bsm_err = cfg["ATLAS_Bbb_SM"]
    bb, _, bb_plus = cfg["CMS_ratio_bb_ZZ_normalized"]
    ww, ww_minus, _ = cfg["CMS_ratio_WW_ZZ_normalized"]
    b_ratio_sm, db = cfg["CMS_prediction_bb_ZZ"]
    w_ratio_sm, dw = cfg["CMS_prediction_WW_ZZ"]
    transfer = cfg["WH_additive_surplus"]
    wh_min = wh - n*wh_minus
    bsm_min = bsm-n*bsm_err
    required_b = (wh_min-transfer)*bsm_min
    # Includes WW/ZZ interference and transfer to the physical four-fermion width.
    four_over_bb = cfg["four_fermion_rate_transfer_lower"]*(
        1+(ww-n*ww_minus)*(w_ratio_sm-n*dw))/((bb+n*bb_plus)*(b_ratio_sm+n*db))
    k = cfg["matched_example_k"]
    c_lep = 3*cfg["lepton_width_kernel_lower"]/cfg["bottom_width_kernel_upper"]
    allowed_lep = 1/required_b - 1-four_over_bb
    pure_bmax = 1/(1+four_over_bb+c_lep*k*k)
    return {"fourfermion_over_bb_lower":four_over_bb,
            "BR_bb_required_lower_from_WH":required_b,
            "BR_bb_pure45_upper":pure_bmax,
            "WH_mu_pure45_upper":pure_bmax/bsm_min+transfer,
            "WH_mu_data_lower":wh_min,
            "residual_lower_at_matched_k1":max(0,1-math.sqrt(allowed_lep/c_lep)/k),
            "residual_lower_without_CMS_at_k1":max(0,1-math.sqrt((1/required_b-1)/c_lep)/k),
            "pure45_necessary_k_upper":math.sqrt(allowed_lep/c_lep),
            "pure45_necessary_normalized_production_lower":required_b/pure_bmax,
            "critical_allowed_WH_additive_surplus_at_k1":wh_min-pure_bmax/bsm_min}


def common_rows(singular, v1, v2):
    # v1,v2 are eigenvalue-squared targets for the shared two-row principal block.
    lam = singular**2
    w1 = np.prod(v1-lam)/(v1*(v2-v1))
    w2 = -np.prod(v2-lam)/(v2*(v2-v1))
    w3 = np.prod(lam)/(v1*v2)
    assert min(w1,w2,w3)>0
    return np.array([[math.sqrt(v1),0.,0.],[0.,math.sqrt(v2),0.],np.sqrt([w1,w2,w3])])


def ckm(angles):
    s12,s13,s23,phase = [row[0] for row in angles]
    a,b,c = (math.sqrt(1-v*v) for v in (s12,s13,s23))
    z = np.exp(1j*phase)
    return np.array([[a*b,s12*b,s13/z],[-s12*c-a*s23*s13*z,a*c-s12*s23*s13*z,s23*b],
                     [s12*s23-a*c*s13*z,-a*s23-s12*c*s13*z,c*b]])


def norm2(mat):
    return float(np.vdot(mat,mat).real)


def mass_example(cfg):
    data=load(cfg["charged_masses_at_MZ"])
    u,d,e = (np.array(data["yukawa"][label])[:,0] for label in ["up","down","charged_leptons"])
    d2,e2=d**2,e**2
    first=math.sqrt(max(d2[0],e2[0])*min(d2[1],e2[1]))
    second=math.sqrt(max(d2[1],e2[1])*min(d2[2],e2[2]))
    Yd, base_e = common_rows(d,first,second),common_rows(e,first,second)
    UL, singular_d, _ = np.linalg.svd(Yd)
    UL=UL[:,::-1]
    V=ckm(data["quark_sines_and_phase"])
    Uu=UL@V.conj().T
    Yu=Uu@np.diag(u)
    assert np.allclose(Uu.conj().T@UL,V) and np.allclose(np.linalg.svd(Yu,compute_uv=False)[::-1],u)
    out={}
    for angle in [0,90,180]:
        E=base_e.astype(complex)
        E[2]*=np.exp(1j*math.radians(angle))
        A,B=(3*Yd+E)/4,(Yd-E)/4
        assert np.allclose(np.linalg.svd(E,compute_uv=False)[::-1],e,rtol=1e-10,atol=1e-15)
        assert np.allclose(singular_d[::-1],d,rtol=1e-10,atol=1e-15)
        assert np.linalg.matrix_rank(B,tol=1e-13)==1
        assert np.allclose(A+B,Yd) and np.allclose(A-3*B,E)
        beta=math.radians(cfg["positive_norm_vev_angle_witness_degrees"])
        # Second CP-even unit eigenvector is orthogonal to the 45-degree vacuum.
        r,s=-math.tan(beta),1/math.tan(beta)
        Sd,Sl=r*A+s*B,r*A-3*s*B
        a,b,c,dcoef=(3*r+s)/4,(r-s)/4,3*(r-s)/4,(r+3*s)/4
        lhs=c*dcoef*norm2(Sd)-a*b*norm2(Sl)
        rhs=r*s*(a*c*norm2(Yd)-b*dcoef*norm2(E))
        assert math.isclose(lhs,rhs,rel_tol=1e-11,abs_tol=1e-18)
        assert math.isclose(norm2(-3*B)/norm2(B),9.)
        out[str(angle)]={"correction_Frobenius":math.sqrt(norm2(B)),
                         "orthogonal_scalar_lepton_over_down_squared_norm":norm2(Sl)/norm2(Sd),
                         "radial_lepton_over_down_squared_norm":norm2(E)/norm2(Yd)}
    return out


def calculate():
    cfg=load("inputs.json")
    return {"observations":{str(n):width_test(cfg,n) for n in cfg["working_widths"]},
            "common_rank_one_examples_at_MZ":mass_example(cfg)}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(v,b[k]) for k,v in a.items())
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-15)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    modes=cli.add_mutually_exclusive_group(required=True)
    modes.add_argument("--check",action="store_true")
    modes.add_argument("--write-new",action="store_true")
    args=cli.parse_args()
    data=calculate()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(data,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        if not same(data,load("results.json")):
            raise SystemExit("1602 certificate mismatch")
        print("1602 verified: two width budgets, independent WH holdout and three rank-one charged-mass witnesses.")
