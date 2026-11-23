"""Round 1613: residual reflection axes, finite-flavor mixing, and independent tree decays."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DATA = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
OLD = json.loads((HERE.parent/"1381"/"inputs.json").read_text(encoding="utf-8"))
TH = math.pi / DATA["group"]["polygon_order"]


def window(value, err, widths):
    return [value-widths*err, value+widths*err]


def divide(pos_a, pos_b):
    assert min(pos_b)>0 and min(pos_a)>0
    return [pos_a[0]/pos_b[1],pos_a[1]/pos_b[0]]


def angle_floor(value, typ):
    if typ=="ud":
        return max(0., math.acos(value)-TH)
    return max(0., math.atan(value)-TH)


def block(axis,theta,phase):
    co,si=math.cos(axis),math.sin(axis)
    c,s=math.cos(theta),math.sin(theta)
    a=np.array([[co,-si,0],[si,co,0],[0,0,1.]])
    b=np.array([[1,0,0],[0,c,s*np.exp(1j*phase)],
                [0,-s*np.exp(-1j*phase),c]],complex)
    return a@b


def reflection(axis):
    v=np.array([math.cos(axis),math.sin(axis),0.])
    return 2*np.outer(v,v)-np.eye(3)


def example(theta, bu, complex_phase):
    u=block(0.,.045,complex_phase)
    d=block(theta,math.asin(bu/math.sin(theta)),0.)
    for f,axis,values in [(u,0.,[.000007,.0035,.96]),(d,theta,[.000015,.00031,.016])]:
        h=f@np.diag(np.square(values))@f.conj().T
        ref=reflection(axis)
        assert np.max(abs(ref@h@ref-h))<2e-15
    mix=u.conj().T@d
    assert np.max(abs(mix.conj().T@mix-np.eye(3)))<1e-13
    assert math.isclose(abs(mix[0,2]),bu,rel_tol=1e-13)
    j=(mix[0,0]*mix[1,1]*mix[0,1].conjugate()*mix[1,0].conjugate()).imag
    return {"moduli":abs(mix).tolist(),"J":float(j)}


def analyze():
    cosine,sine,tan=math.cos(TH),math.sin(TH),math.tan(TH)
    ru,rd=reflection(0.),reflection(TH)
    rr=rd@ru
    assert np.allclose(np.linalg.matrix_power(rr,14),np.eye(3))
    assert all(not np.allclose(np.linalg.matrix_power(rr,k),np.eye(3)) for k in range(1,14))
    assert np.allclose(ru@rr@ru,np.linalg.inv(rr))
    names=OLD["ub_over_cb"]
    bc=OLD["F_times_cb"]
    bu=bc[0]/OLD["F_example"]*names[0]
    results={}
    axial=DATA["leptonic"]
    for n in [1,2]:
        bmax=(bc[0]+n*sum(bc[1:]))/OLD["F_lower"]*(names[0]+n*sum(names[1:]))
        lowus=math.sqrt(sine**2-bmax**2)
        f=window(*OLD["K_form"],n)
        young=DATA["kse3"]
        vec={
            "old":divide(window(*OLD["K_product"],n),f),
            "KSe3_body":divide(window(young["product"],young["body_error"],n),f),
            "KSe3_abstract":divide(window(young["product"],young["abstract_error"],n),f)}
        kr=divide(window(axial["combination"],axial["branching_error"]+axial["radiative_error"],n),
                  [axial["decay_constant_ratio"]-n*axial["lower"],
                   axial["decay_constant_ratio"]+n*axial["upper"]])
        nuclear={}
        for key,case in [("2020",OLD["ud_beta"]),("2025",DATA["new_nuclear"])]:
            ud=window(*case,n)
            theta_floor=angle_floor(ud[1],"ud")
            usfloor=math.sqrt(1-ud[1]**2-bmax**2)
            chosen_us=math.sqrt(1-ud[1]**2-bu**2)
            nuclear[key]={"ud":ud,"angle_floor_deg":math.degrees(theta_floor),
                          "minimum_endpoint_sum_projective_rad":theta_floor,
                          "positive_rate_excess_at_fixed_nuclear_response":(cosine/ud[1])**2-1,
                          "unitary_us_floor_at_bmax":usfloor,
                          "historic_vector_upper_minus_unitary_floor":vec["old"][1]-usfloor,
                          "chosen_us_at_upper_ud":chosen_us,
                          "chosen_us_over_ud":chosen_us/ud[1]}
        results[str(n)]={"beta":nuclear,"ub_upper":bmax,"exact_us_range":[lowus,sine],
                         "semileptonic":vec,"axial_ratio":kr,
                         "axial_angle_floor_deg":math.degrees(angle_floor(kr[0],"ratio")),
                         "minimum_axial_amplitude_transfer":kr[0]/tan,
                         "axial_reference_rate_factor":(kr[0]/tan)**2,
                         "KSe3_vs_predicted_products":[f[0]*lowus,f[1]*sine]}
    two=results["2"]
    ex=two["beta"]["2025"]
    witness=example(math.acos(ex["ud"][1]),bu,1.2)
    assert two["semileptonic"]["KSe3_body"][0]<witness["moduli"][0][1]<two["semileptonic"]["KSe3_body"][1]
    assert two["axial_ratio"][0]<ex["chosen_us_over_ud"]<two["axial_ratio"][1]
    assert ex["historic_vector_upper_minus_unitary_floor"]<0
    cases={str(x):example(TH,bu,x) for x in (0.,1.2,-1.2)}
    assert abs(cases["0.0"]["J"])<1e-14 and abs(cases["1.2"]["J"])>1e-6
    assert abs(cases["1.2"]["J"]+cases["-1.2"]["J"])<1e-14
    return {"axis_angle_rad":TH,"exact_ud":cosine,"exact_max_us":sine,"exact_max_us_over_ud":tan,
            "center_ub_illustrative":bu,"width_tests":results,"exact_unitary_examples":cases,
            "deformed_witness_at_2025_beta_upper":witness}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [rounded(v) for v in x]
    return float(f"{x:.10g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=rounded(analyze())
    out=HERE/"results.json"
    if args.check:
        assert result==json.loads(out.read_text(encoding="utf-8"))
    else:
        with out.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
            f.write("\n")
    print("1613: two reflections, three-flavor completions, three weak-current tests verified")
