"""Exact tree-level standard-doublet mixing, independent weak and scalar ports."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def diag(mass):
    left, singular, right_h = np.linalg.svd(mass)
    order = np.argsort(singular)
    return singular[order], left[:, order], right_h.conj().T[:, order]


def width(parent, daughter, scalar, left, right):
    if parent <= daughter + scalar:
        return 0.
    phase = math.sqrt((parent**2 - (daughter+scalar)**2)*(parent**2 - (daughter-scalar)**2))
    amplitude = ((abs(left)**2+abs(right)**2)*(parent**2+daughter**2-scalar**2)
                 + 4*parent*daughter*(left*right.conjugate()).real)
    return float(phase*amplitude/(32*math.pi*parent**3))


def one_example(B, config):
    mt = config["m_top_reference_GeV"]
    v = config["v_GeV"]
    ell = mt**2
    loss = config["left_W_loss"]
    kappa = 1 - B**2/ell * loss
    assert 0 < kappa < 1
    T2 = ell+(B**2-ell)/kappa
    T = math.sqrt(T2)
    x = math.sqrt(ell+T2-B**2)
    f = ell*T2/(B**2*x*x)
    sc = math.sqrt(f)
    cc = math.sqrt(1-f)
    basis = np.array([[-cc,sc],[sc,cc]])
    upper_y = np.array([[0.,0.],[x,0.]])
    top = np.array([[0.,B*sc],[x,B*cc]])
    assert np.allclose(basis.T@top,[[x*sc,0],[x*cc,B]])
    bottom = basis @ np.diag([config["m_bottom_example_GeV"],B])
    masses,L,R = diag(top)
    down_masses,Lb,Rb = diag(bottom)
    assert np.allclose(masses,[mt,T]) and np.allclose(down_masses,[config["m_bottom_example_GeV"],B])
    Yh = L.T @ upper_y @ R / v
    # Both LH states are weak doublets; only the heavy RH weak state has nonzero isospin.
    WL = L.T@Lb
    WR = R.T@np.diag([0.,1.])@Rb
    singular_projection = R.T@np.diag([1.,0.])@R
    assert np.allclose(WL@WL.T,np.eye(2)) and np.allclose(WL.T@WL,np.eye(2))
    assert np.allclose(WR@WR.T,np.eye(2)-singular_projection)
    assert abs(WR[0,0])<1e-12 and np.isclose(abs(WL[0,0])**2,1-loss)
    assert np.allclose(v*Yh,np.diag(masses)@singular_projection)
    assert np.allclose(np.diag(v*Yh)/masses,[kappa,1-kappa])
    off = float(abs(Yh[0,1])**2+abs(Yh[1,0])**2)
    assert np.isclose(off,(ell+T2)*kappa*(1-kappa)/v**2)
    pole_slope = []
    for factor in (1-1e-5,1+1e-5):
        varied,_,_ = diag(top+(factor-1)*upper_y)
        pole_slope.append(np.log(varied))
    assert np.allclose((pole_slope[1]-pole_slope[0])/(2e-5),[kappa,1-kappa])
    return {"m_B":B,"m_T":T,"site_fraction":f,"underlying_Y":math.sqrt(2)*x/v,
            "kappa_t":kappa,"kappa_T":1-kappa,"light_left_W_squared":1-loss,
            "light_right_Z_doublet_fraction":1-kappa,"ttH_norm":kappa**2,
            "mixed_higgs_strength_squared":off,
            "T_to_th_width_GeV":width(T,mt,config["m_h_GeV"],Yh[1,0].conjugate(),Yh[0,1])}


def minimum_f(mt,B,kappa):
    r=(mt/B)**2
    deficit=1-kappa
    return r*(1-r*deficit)/(r+(1-2*r)*deficit)


def calculate(data):
    config=data["benchmarks"]
    mt=config["m_top_reference_GeV"]
    examples=[one_example(B,config) for B in config["m_B_GeV"]]
    s=data["signal"]
    allowed=[math.sqrt(s["mu"]-s["minus"]),min(1,math.sqrt(s["mu"]+s["plus"]))]
    # Unobserved heavy masses are hypothetical. The mass input is not a rate calibration.
    bridge=data["charm_bridge"]
    old=json.loads((ROOT/bridge["prev_results"]).read_text(encoding="utf-8-sig"))
    epsilon=1-old["MZ"]["outer_bounds"]["2"][bridge["key"]]
    charm=bridge["yc_MZ"]*bridge["v_MZ_GeV"]/math.sqrt(2)
    B=bridge["m_B_working_GeV"]
    if not epsilon*B*B>charm*charm:
        raise ValueError("positive finite threshold split is impossible")
    split=charm**2*B**2*(1-epsilon)/(epsilon*B**2-charm**2)
    mass_gap=split/(math.sqrt(B**2+split)+B)
    assert np.isclose(epsilon,charm**2*(B**2+split)/(B**2*(charm**2+split)))

    # A named additional scalar keeps its charge-site fixed while all EW matrices
    # are co-rotated. Remove basis ambiguities by comparing in the same mass basis.
    singlet=data["additional_scalar"]
    ref=examples[0]
    x=ref["underlying_Y"]*config["v_GeV"]/math.sqrt(2)
    shared_y=np.array([x*math.sqrt(ref["site_fraction"]),x*math.sqrt(1-ref["site_fraction"])])
    B=ref["m_B"]
    universal=np.column_stack([shared_y,[0.,B]])
    tagged=[]
    for f in singlet["test_f"]:
        s0,c0=math.sqrt(f),math.sqrt(1-f)
        frame=np.array([[-c0,s0],[s0,c0]])
        qm=frame@universal
        down=frame@np.diag([config["m_bottom_example_GeV"],B])
        assert np.allclose(frame.T@qm,universal)
        assert np.allclose(frame.T@down,np.diag([config["m_bottom_example_GeV"],B]))
        scalar=frame.T@np.array([0.,singlet["unit_yukawa"]])
        value=width(B,config["m_bottom_example_GeV"],singlet["m_S_GeV"],scalar[0],0.)
        tagged.append({"site_fraction":f,"S_bB":scalar[0],"S_BB":scalar[1],"B_to_bS_width_GeV":value})
    return {"benchmarks":examples,"ttH_kappa_window_if_rescaled":allowed,
            "lower_f_from_ttH_at_registered_masses":[minimum_f(mt,b,allowed[0]) for b in config["m_B_GeV"]],
            "charm_conditional_bridge":{"MZ_mass":charm,"allowed_site_share":epsilon,
                "minimum_square_mass_split_GeV2":split,"minimum_mass_gap_GeV":mass_gap},
            "tagged_fixed_weak_and_higgs":tagged,
            "max_per_model_equal_error_th_GeV":abs(examples[0]["T_to_th_width_GeV"]-examples[1]["T_to_th_width_GeV"])/2,
            "max_per_model_equal_error_bS_GeV":abs(tagged[0]["B_to_bS_width_GeV"]-tagged[1]["B_to_bS_width_GeV"])/2}


def close(got,expected):
    if isinstance(got,dict):
        return got.keys()==expected.keys() and all(close(v,expected[k]) for k,v in got.items())
    if isinstance(got,list):
        return len(got)==len(expected) and all(close(x,y) for x,y in zip(got,expected))
    return bool(np.isclose(got,expected,rtol=2e-10,atol=1e-12))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=calculate(json.loads((ROOT/"inputs.json").read_text(encoding="utf-8-sig")))
    path=ROOT/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        assert close(result,json.loads(path.read_text(encoding="utf-8-sig")))
        print("1399 masses, weak and scalar currents, transfer and tagged comparison verified")
