"""Local weak-current completions from independent charged-current strengths."""
import argparse
import json
from math import acos, cos, pi, radians, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def limits(values, n=1):
    center, *parts = values
    width = n*sum(parts)
    return center-width, center+width


def third_dimension_c(a, b, d, angle):
    linear = a*b*d*cos(radians(angle))
    positive = (1-a*a)*(1-d*d)-b*b
    return (linear + sqrt(linear**2+(1-b*b)*positive))/(1-b*b)


def block(a, b, c, d, angle):
    # gamma = arg[-Vud Vub* / (Vcd Vcb*)] in the stated phase convention.
    return np.array([[a, -b*np.exp(-1j*radians(angle))], [c, d]])


def extend_rows(rows):
    _, _, vh = np.linalg.svd(rows,full_matrices=True)
    unitary = np.vstack([rows, vh[rows.shape[0]:]])
    np.testing.assert_allclose(unitary@unitary.conj().T,np.eye(len(unitary)),atol=2e-12)
    return unitary


def three_witness(a,b,d,angle):
    c = third_dimension_c(a,b,d,angle)
    K = block(a,b,c,d,angle)
    D = np.eye(2)-K@K.conj().T
    col_u = sqrt(D[0,0].real)
    col = np.array([col_u, D[1,0]/col_u])
    U = extend_rows(np.column_stack([K,col]))
    return U


def four_witness(a,b,c,d,angle,us,cs):
    K = block(a,b,c,d,angle)
    D = np.eye(2)-K@K.conj().T
    p, q, overlap = D[0,0].real,D[1,1].real,D[1,0]
    first = np.array([us,sqrt(p-us*us)])
    perpendicular = np.array([-first[1],first[0]]) / sqrt(p)
    fixed = overlap/p*first
    free = sqrt(q-abs(overlap)**2/p)*perpendicular
    # Unknown phase in the two missing columns fixes |Vcs|, if inside its annulus.
    lo, hi = abs(abs(fixed[0])-abs(free[0])), abs(fixed[0])+abs(free[0])
    assert lo <= cs <= hi
    cosine = (cs*cs-abs(fixed[0])**2-abs(free[0])**2)/(2*abs(fixed[0]*free[0]))
    eta = np.angle(fixed[0])-np.angle(free[0])+acos(cosine)
    second = fixed+np.exp(1j*eta)*free
    U = extend_rows(np.column_stack([K,np.stack([first,second])]))
    np.testing.assert_allclose(abs(U[1,2]),cs,atol=2e-12)
    return U, [lo,hi], np.linalg.eigvalsh(D).tolist()


def calc(x):
    a=x["ud_beta"][0]
    d=x["F_times_cb"][0]/x["F_example"]
    b=d*x["ub_over_cb"][0]
    c=x["charm_d"][0]
    g=x["gamma_conditional_deg"][0]
    us=x["K_product"][0]/x["K_form"][0]
    cs=x["charm_s"][0]
    three=three_witness(a,b,d,g)
    four,cs_range,defects=four_witness(a,b,c,d,g,us,cs)
    angle_data=[]
    for dg in [0,180]:
        angle_data.append([third_dimension_c(a,b,d,v+dg) for v in x["gamma_conditional_deg"]])
    data_boxes={}
    for width in [1,2]:
        al,ah=limits(x["ud_beta"],width)
        dmax=limits(x["F_times_cb"],width)[1]/x["F_lower"]
        bmax=dmax*limits(x["ub_over_cb"],width)[1]
        kl,kh=limits(x["K_product"],width)
        fl,fh=limits(x["K_form"],width)
        k_lo,k_hi=kl/fh,kh/fl
        min_k_three=sqrt(1-ah**2-bmax**2)
        max_k_three=sqrt(1-al**2)
        loss_min=1-ah**2-bmax**2-k_hi**2
        loss_max=1-al**2-k_lo**2
        c_lo=(sqrt((1-ah**2-bmax**2)*(1-dmax**2-bmax**2))-ah*bmax*dmax)
        c_hi=(sqrt(1-al**2)+ah*bmax*dmax)/(1-bmax**2)
        data_boxes[str(width)]={"Vus_three": [min_k_three,max_k_three],
                                "Vus_measured":[k_lo,k_hi],"row_loss":[loss_min,loss_max],
                                "Vcd_three_outer_any_angle":[c_lo,c_hi]}
    assert data_boxes["1"]["row_loss"][0]>0>data_boxes["2"]["row_loss"][0]
    assert min(third_dimension_c(a,b,d,0),third_dimension_c(a,b,d,180))>c
    assert limits(x["charm_d"])[0] < abs(three[1,0]) < limits(x["charm_d"])[1]
    # A concrete three-dimensional example in the enlarged input boxes.
    a2=limits(x["ud_beta"],2)[1]
    d2=.0495
    b2=d2*limits(x["ub_over_cb"],2)[1]
    witness2=three_witness(a2,b2,d2,g)
    k2=abs(witness2[0,2])
    assert data_boxes["2"]["Vus_measured"][0] <= k2 <= data_boxes["2"]["Vus_measured"][1]
    for observed,value in [("charm_d",witness2[1,0]),("charm_s",witness2[1,2])]:
        lower,upper=limits(x[observed],2)
        assert lower <= abs(value) <= upper
    return {"Vcd_three_any_phase_central_other_moduli":[third_dimension_c(a,b,d,180),third_dimension_c(a,b,d,0)],
            "conditional_Vcd_gamma_center_left_right":angle_data,
            "nuclear_and_kaon_boxes":data_boxes,"central_four_defect_eigenvalues":defects,
            "central_four_allowed_Vcs":cs_range,"central_missing_u_amplitude":abs(four[0,3]),
            "three_central_inputs_predicted_Vcd_Vus":[abs(three[1,0]),abs(three[0,2])],
            "explicit_four_unitary_real":four.real.tolist(),"explicit_four_unitary_imag":four.imag.tolist(),
            "doubled_box_three_abs_matrix":np.abs(witness2).tolist()}


def cmp(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k,v in a.items():
            cmp(v,b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-9,atol=1e-12)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=calc(json.loads((HERE/"inputs.json").read_text(encoding="utf-8")))
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        cmp(out,json.loads(path.read_text(encoding="utf-8")))
        print("1381: contraction ranks, two held-out charm/strange reads and real-data witnesses verified.")
