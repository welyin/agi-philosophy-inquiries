"""Bounds on the radial angular velocity, with independent top and VBF readouts."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def row_weight(ratio):
    return ratio * ratio / (1 + ratio * ratio)


def gap_factor(y_heavy, y_light):
    # The physical scalar coupling has m/v = y/sqrt(2).
    return (y_heavy**2 - y_light**2)**2 / (2*(y_heavy**2 + y_light**2))


def unitary(z, angle):
    c, s = math.cos(angle), math.sin(angle)
    out = np.eye(3, dtype=complex)
    out[:2,:2] -= (1-c)*np.outer(z, z.conj())
    out[:2,2] = s*z
    out[2,:2] = -s*z.conj()
    out[2,2] = c
    return out


def endpoint_example(weak, inp, yuk):
    # A long-arc CKM that already passes the independent weak readouts of 1614.
    r, t, vcb = .096, .208624126345, .04424
    p, q = row_weight(r), row_weight(t)
    angle = math.pi/weak["angle_denominator"]
    C, S = math.cos(angle), math.sin(angle)
    phase = -math.acos((C*C*p+S*S*(1-p)-q)/(2*C*S*math.sqrt(p*(1-p))))
    z = np.array([math.sqrt(p),math.sqrt(1-p)*np.exp(1j*phase)])
    R = np.eye(3)
    R[:2,:2] = [[C,S],[-S,C]]
    delta = math.pi - math.asin(vcb * math.sqrt(1+r*r))
    cu, cd = -.45*delta, .55*delta
    mats = [np.diag([x[0] for x in yuk[sector]])/math.sqrt(2) for sector in ("up","down")]
    def matrix(rho, a, start, mass, flat):
        x = rho*rho
        f = (2*x-x*x) if flat else x
        return rho*unitary(z, a*f)@start@mass
    du = unitary(z,cu)
    dd = unitary(z,cd)@R
    V = du.conj().T@dd
    h = 1e-5
    flat_errors = []
    direct_errors = []
    for a, start, mass in zip((cu,cd),(np.eye(3),R),mats):
        U = unitary(z,a)@start
        for flat, errors in [(True,flat_errors),(False,direct_errors)]:
            Y = U.conj().T@(matrix(1+h,a,start,mass,flat)-matrix(1-h,a,start,mass,flat))/(2*h)
            errors.append(float(np.linalg.norm(Y-np.diag(np.diag(Y)))))
        assert np.allclose(matrix(1,a,start,mass,False),matrix(1,a,start,mass,True))
    assert max(flat_errors)<2e-8 and min(direct_errors)>1e-3
    assert abs(abs(V[0,0])-.97434994)<3e-7
    # Put the whole displacement in the down sector, keep the right frame fixed.
    # For the direct-lineshape comparison, set the kinematic final masses to zero,
    # which here increases the conservative upper estimate of the extra width.
    max_down_width = (3*inp["kinematics_GeV"]["h_center"]/(8*math.pi)
                      *(2*delta)**2*(yuk["down"][2][0]**2+yuk["down"][1][0]**2)/2
                      *(1+(5.5/inp["kinematics_GeV"]["h_center"])**2))
    return {"delta_rad":delta,"Vud":abs(V[0,0]),"Vus":abs(V[0,1]),
            "flat_numerical_offdiagonal_norms":flat_errors,
            "linear_numerical_offdiagonal_norms":direct_errors,
            "all_down_linear_upper_h_width_with_reference_GeV":max_down_width+inp["h_reference_width_GeV"]}


def check_minimum():
    # Arbitrary complex left velocity; check the attained right-rotation minimum.
    a = .49+.89j
    di,dj=.37,3.21
    best = 2*di*dj*a/(di*di+dj*dj)
    floor=(dj*dj-di*di)**2/(dj*dj+di*di)*abs(a)**2
    for shift in (0,.61-.32j,-1.3+.14j):
        b=best+shift
        norm=abs(dj*a-di*b)**2+abs(di*a-dj*b)**2
        assert math.isclose(norm,floor+(di*di+dj*dj)*abs(shift)**2,rel_tol=1e-13)


def run():
    inputs=read(HERE/"inputs.json")
    weak=read(HERE/inputs["weak"])
    yuk=read(HERE/inputs["yukawa"])["yukawa"]
    n=inputs["marginal_widths"]
    def hi(name):
        val=weak[name]
        return val["central"]+n*sum(val["errors"])
    def yend(sector,index,side):
        y,epsilon=yuk[sector][index]
        return y+side*n*epsilon
    hs=hi("cb")*math.sqrt(1+hi("ub_over_cb")**2)
    displacement_floor=math.pi-math.asin(hs)
    kin=inputs["kinematics_GeV"]
    hmin=kin["h_center"]-n*kin["h_sigma"]
    hmax=kin["h_center"]+n*kin["h_sigma"]
    mtlo,mthi=kin["top"]
    rh,rc=hmax/mtlo,kin["charm_upper"]/mtlo
    top_phase=math.sqrt((1-(rh+rc)**2)*(1-rh**2))*(1-rh**2-2*rc)
    bottom_phase=(1-(kin["b_plus_s_upper"]/hmin)**2)**2
    ku=mtlo/(32*math.pi)*top_phase*gap_factor(yend("up",2,-1),yend("up",1,1))/(1+hi("ub_over_cb")**2)
    kd=3*hmin/(8*math.pi)*bottom_phase*gap_factor(yend("down",2,-1),yend("down",1,1))/(1+hi("td_over_ts")**2)
    transfer=inputs["transfer"]
    rate=transfer["tree_FCNC_width_min"]
    vbf=inputs["ATLAS_VBF_WW_pb"]
    sigmin=vbf["value"]-n*vbf["lower_error"]
    smmax=vbf["SM"]+n*vbf["SM_error"]
    dilution = sigmin/(smmax*transfer["VBF_WW_numerator_and_acceptance_max"])
    width_h=inputs["h_reference_width_GeV"]/dilution
    top=inputs["CMS_top_same_sign"]
    width_t=top["top_width_upper_GeV"]
    response=transfer["top_search_response_min_excluding_h_dilution"]
    up_direct=math.sqrt(top["B_ch_upper"]*width_t/(rate*ku*response*dilution))
    down=math.sqrt(width_h/(rate*kd))
    combined=up_direct+down
    Wlo,Whi=kin["W"]
    bmax=kin["b_plus_s_upper"]
    wb_phase=math.sqrt((1-((Whi+bmax)/mtlo)**2)*(1-(Whi/mtlo)**2))
    wb_shape=(1-(bmax/mtlo)**2)**2+(Wlo/mthi)**2-2*(Whi/mtlo)**4
    wb=inputs["GF_lower_GeV_inverse2"]*mtlo**3/(8*math.pi*math.sqrt(2))*(1-hs**2)*wb_phase*wb_shape
    wb_phys=wb*transfer["tree_Wb_width_min"]
    up_width=math.sqrt((width_t-wb_phys)/(rate*ku))
    positive_min=2*displacement_floor
    assert combined<up_width+down<positive_min
    check_minimum()
    return {"working_h_sin_upper":hs,"long_displacement_min_rad":displacement_floor,
            "same_sign_endpoint_velocity_min":positive_min,"tree_GeV_per_unit_velocity_squared":{"t_ch":ku,"h_bs":kd},
            "VBF_dilution_min":dilution,"VBF_total_width_max_GeV":width_h,
            "with_FCNC_search":{"max_up":up_direct,"max_down":down,"max_difference":combined,
                "max_normalized_endpoint_Fprime":combined/positive_min,
                "minimum_opposite_weighted_n_minus_one":1-combined/positive_min,
                "minimum_width_h_GeV_if_direct_up_bound_fixed":rate*kd*(positive_min-up_direct)**2},
            "without_FCNC_search":{"Wb_tree_min_GeV":wb,"Wb_physical_min_GeV":wb_phys,
                "max_up":up_width,"max_difference":up_width+down,
                "critical_top_total_width_GeV":wb_phys+rate*ku*(positive_min-down)**2},
            "endpoint_example":endpoint_example(weak,inputs,yuk)}


def compare(old,new):
    if isinstance(old,dict):
        assert old.keys()==new.keys()
        for k in old:
            compare(old[k],new[k])
    elif isinstance(old,list):
        assert len(old)==len(new)
        for a,b in zip(old,new):
            compare(a,b)
    else:
        assert math.isclose(old,new,rel_tol=1e-9,abs_tol=1e-10)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    path=HERE/"results.json"
    if args.check:
        compare(read(path),result)
        print("1618 OK: finite chirality, radial velocities, two width tests and equal-endpoint profiles")
    else:
        with path.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,ensure_ascii=False)
            f.write("\n")
