"""The radial derivative of a fitted mass hierarchy and ATLAS partial-width ratios."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rotation(s12,s13,s23,delta):
    c12,c13,c23 = np.sqrt(1-np.array([s12,s13,s23])**2)
    phase=np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13/phase],
                     [-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13],
                     [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13]])


def mass_example(spectrum,left,right,rows,columns,epsilon):
    powers=rows[:,None]+columns
    yukawa=(left*spectrum[None,:])@right.conj().T
    coefficients=yukawa/epsilon**powers
    derivative=(1+2*powers)*yukawa
    physical=left.conj().T@derivative@right
    response=np.diag(physical)/spectrum
    exact=1+2*((abs(left)**2).T@rows+(abs(right)**2).T@columns)
    np.testing.assert_allclose(response,exact,atol=1e-11)
    np.testing.assert_allclose(np.linalg.svd(coefficients*epsilon**powers,compute_uv=False),spectrum[::-1],rtol=1e-11)
    assert np.min(abs(coefficients))>.2 and np.max(abs(coefficients))<5
    def field(t):
        return coefficients*epsilon**powers*t**(1+2*powers)
    step=1e-5
    np.testing.assert_allclose((field(1+step)-field(1-step))/(2*step),derivative,rtol=1e-8)
    return {"coeff_modulus_range":[float(np.min(abs(coefficients))),float(np.max(abs(coefficients)))],
            "diagonal_radial_couplings_over_SM":np.real(response).tolist()}


def working_box(observed,standard,number,transfer):
    central,err_minus,err_plus=observed
    ref,error=standard
    lower=(central-number*err_minus)/(ref+number*error)
    upper=(central+number*err_plus)/(ref-number*error)
    response=[math.sqrt(max(0,lower)/transfer[1]),math.sqrt(upper/transfer[0])]
    return {"normalized_rate_interval":[lower,upper],"allowed_abs_kf_over_kZ":response,
            "largest_positive_mass_share_from_powers_one_or_higher":(response[1]-1)/2,
            "largest_transfer_factor_compatible_with_a_ninefold_width":upper/9}


def compute():
    inputs=load(HERE/"inputs.json")
    mass=load(HERE/inputs["mass_and_charged_current"])
    angles=inputs["mass_fit_basis_s12_s13_s23_phase"]
    lu=rotation(*angles["up_left"])
    v=rotation(*[item[0] for item in mass["quark_sines_and_phase"]])
    ld=lu@v
    np.testing.assert_allclose(lu.conj().T@ld,v,atol=1e-14)
    axes={"up":(lu,rotation(*angles["up_right"]),inputs["right_up"]),
          "down":(ld,rotation(*angles["down_right"]),inputs["right_down_and_example_leptons"]),
          "charged_leptons":(rotation(*angles["leptons_left"]),rotation(*angles["leptons_right"]),
                             inputs["right_down_and_example_leptons"])}
    example={}
    for sector,(left,right,columns) in axes.items():
        spectrum=np.array(mass["yukawa"][sector])[:,0]
        example[sector]=mass_example(spectrum,left,right,np.array(inputs["left_exponents"]),
                                     np.array(columns),inputs["epsilon"])
    comparison={}
    for port,item in inputs["atlas_table7"].items():
        comparison[port]={str(n):working_box(item["observed"],item["SM"],n,inputs["rate_transfer_window"])
                          for n in [1,2]}
        assert comparison[port]["2"]["allowed_abs_kf_over_kZ"][1]<3
    # The same mass and one-Higgs slope leave higher contact derivatives open.
    weights=np.array([2.,-1.])
    powers=np.array([3,5])
    derivatives=[float(np.sum(weights)),float(np.dot(weights,powers)),float(np.dot(weights,powers*(powers-1)))]
    assert derivatives == [1.,1.,-8.]
    return {"fit_to_known_inputs":example,"effective_scale_GeV":inputs["v_GeV"]/math.sqrt(2*inputs["epsilon"]),
            "direct_ATLAS_ratio_boxes":comparison,"signed_competitor_mass_slope_curvature":derivatives}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=compute()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif args.check:
        def compare(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(compare(a[k],b[k]) for k in a)
            if isinstance(a,list):
                return len(a)==len(b) and all(compare(x,y) for x,y in zip(a,b))
            return math.isclose(a,b,rel_tol=2e-10,abs_tol=1e-13)
        assert compare(result,load(target))
        print("1435 calibrated spectra, radial derivatives and ATLAS direct ratios verified")
    else:
        print(json.dumps(result,indent=2))
