"""Round 1614: two short/long branches of one shared coherent heavy-light rotation."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def marginal(data, widths):
    err = widths * sum(data["errors"])
    return data["central"]-err, data["central"]+err


def probability(ratio):
    return ratio**2/(1+ratio**2)


def light_rotation(theta):
    c,s = math.cos(theta), math.sin(theta)
    return np.array([[c,s,0.],[-s,c,0.],[0.,0.,1.]],dtype=complex)


def rotate(z,angle):
    c,s=math.cos(angle),math.sin(angle)
    out=np.eye(3,dtype=complex)
    out[:2,:2]-=(1-c)*np.outer(z,z.conj())
    out[:2,2]=s*z
    out[2,:2]=-s*z.conj()
    out[2,2]=c
    return out


def scalar_prediction(theta,p,q,k):
    C,S=math.cos(theta),math.sin(theta)
    return C*C-k*(p+q-S*S)+k*k*p*q


def build(theta,r,t,v,branch,handedness,inputs):
    C,S=math.cos(theta),math.sin(theta)
    p,q=probability(r),probability(t)
    h=v*math.sqrt(1+r*r)
    c=branch*math.sqrt(1-h*h)
    # z=(sqrt(p),sqrt(1-p)*e^{i phi}); physical d axis=(C,-S).
    cos_phi=(C*C*p+S*S*(1-p)-q)/(2*C*S*math.sqrt(p*(1-p)))
    if abs(cos_phi)>1+1e-12:
        raise ValueError("inequivalent heavy projections for this angle")
    phi=handedness*math.acos(max(-1,min(1,cos_phi)))
    z=np.array([math.sqrt(p),math.sqrt(1-p)*cmath.exp(1j*phi)])
    delta=math.atan2(h,c)
    # Nonzero physical actions in both sectors. Their difference alone enters charged currents.
    up=rotate(z,.12)
    down=rotate(z,.12+delta) @ light_rotation(theta)
    V=up.conj().T @ down
    for basis,values in ((up,inputs["parent_mass_square_example_u"]),
                         (down,inputs["parent_mass_square_example_d"])):
        H=(basis*np.array(values))@basis.conj().T
        assert np.allclose(np.linalg.eigvalsh(H),values)
    assert np.allclose(V@V.conj().T,np.eye(3))
    assert abs(abs(V[0,2]/V[1,2])-r)<3e-12
    assert abs(abs(V[2,0]/V[2,1])-t)<3e-12
    assert abs(abs(V[1,2])-v)<3e-12
    ud=abs(V[0,0])
    assert abs(ud**2-scalar_prediction(theta,p,q,1-c))<3e-12
    roots=[]
    vh=V[2,2]
    for eta in (+1,-1):
        # In the negative branch the "denominator" is 1-|Vtb|; use 1-hypot identity.
        denominator=1+abs(vh) if eta==1 else h*h/(1+abs(vh))
        S0=V[:2,:2]-eta*(vh.conjugate()/abs(vh))*np.outer(V[:2,2],V[2,:2])/denominator
        assert np.allclose(S0@S0.conj().T,np.eye(2),atol=1e-11)
        roots.append(float(abs(S0[0,0])))
    assert abs(roots[0 if branch==1 else 1]-C)<2e-11
    gamma=math.degrees(cmath.phase(-V[0,0]*V[0,2].conjugate()/(V[1,0]*V[1,2].conjugate())))
    J=(V[0,0]*V[1,1]*V[0,1].conjugate()*V[1,0].conjugate()).imag
    return {"r_ub_cb":r,"t_td_ts":t,"Vcb":v,"signed_cosine":c,"action_difference_degree":math.degrees(delta),
            "relative_z_phase_degree":math.degrees(phi),"abs_Vud":float(ud),"abs_Vus":float(abs(V[0,1])),
            "gamma_mod_180_deg":gamma%180,"J":float(J),"two_reconstructed_cosines":roots}


def guaranteed_window(theta,inputs,widths,branch):
    rlo,rhi=marginal(inputs["ub_over_cb"],widths)
    tlo,thi=marginal(inputs["td_over_ts"],widths)
    vlo,vhi=marginal(inputs["cb"],widths)
    hmin=vlo*math.sqrt(1+rlo*rlo)
    hmax=vhi*math.sqrt(1+rhi*rhi)
    ks=sorted([1-branch*math.sqrt(1-h*h) for h in (hmin,hmax)])
    plo,phi=map(probability,(rlo,rhi))
    qlo,qhi=map(probability,(tlo,thi))
    assert ks[1]*max(phi,qhi)<1
    # For fixed k the squared modulus decreases with p,q. In k it is convex.
    zero=(phi+qhi-math.sin(theta)**2)/(2*phi*qhi)
    minimizer=max(ks[0],min(ks[1],zero))
    low=math.sqrt(scalar_prediction(theta,phi,qhi,minimizer))-2e-13
    high=math.sqrt(max(scalar_prediction(theta,plo,qlo,k) for k in ks))+2e-13
    ud_lo,ud_hi=marginal(inputs["ud"],widths)
    tree_shift=(1-math.sqrt(1-hmax*hmax))*math.sqrt(phi)*(math.cos(theta)*math.sqrt(phi)+math.sin(theta))
    C,S=math.cos(theta),math.sin(theta)
    cos_phi_lo=(S*S+(C*C-S*S)*plo-qhi)/(2*C*S*math.sqrt(phi*(1-phi)))
    cos_phi_hi=(S*S+(C*C-S*S)*phi-qlo)/(2*C*S*math.sqrt(plo*(1-plo)))
    assert 0 < cos_phi_lo < cos_phi_hi < 1 and phi < .5
    j_floor=math.sqrt(1-hmax*hmax)*hmin**2*C*S*math.sqrt(plo*(1-plo))*math.sqrt(1-cos_phi_hi**2)
    return {"cos_phi_outer":[cos_phi_lo-2e-13,cos_phi_hi+2e-13],
            "actual_abs_J_lower":j_floor-2e-13,
            "min_phase_distance_from_real_axis_degree":math.degrees(math.acos(cos_phi_hi))-2e-10,
            "abs_Vud_outer":[low,high],"beta_own_window":[ud_lo,ud_hi],
            "positive_gap_above_beta_upper":max(low-ud_hi,0),
            "short_arc_tree_only_abs_shift_upper":tree_shift+2e-13}


def run(data):
    theta=math.pi/data["angle_denominator"]
    r,t,v=(data[k]["central"] for k in ("ub_over_cb","td_over_ts","cb"))
    models={}
    for branch in (1,-1):
        models[str(branch)]={"central_conjugates":[build(theta,r,t,v,branch,h,data) for h in (1,-1)],
                             "one_width":guaranteed_window(theta,data,1,branch),
                             "two_widths":guaranteed_window(theta,data,2,branch)}
    example=data["ordinary_long_example"]
    ordinary=build(theta,example["r"],example["t"],example["v"],-1,-1,data)
    other=data["all_readouts_long_example"]
    r,v,ud=other["r"],other["v"],other["ud"]
    k=1+math.sqrt(1-v*v*(1+r*r))
    p=probability(r)
    q=(math.cos(theta)**2-k*(p-math.sin(theta)**2)-ud**2)/(k-k*k*p)
    complete=build(theta,r,math.sqrt(q/(1-q)),v,-1,-1,data)
    prod=marginal(data["kaon_product"],2)
    ff=marginal(data["kaon_form"],2)
    krange=(prod[0]/ff[1],prod[1]/ff[0])
    for label,val in (("ub_over_cb",r),("cb",v),("td_over_ts",complete["t_td_ts"]),("ud",complete["abs_Vud"])):
        assert marginal(data[label],2)[0] < val < marginal(data[label],2)[1]
    assert krange[0]<complete["abs_Vus"]<krange[1]
    gamma=data["charged_B_gamma"]
    gamma_window=[gamma["central_deg"]-gamma["minus"]-gamma["transport_deg"],
                  gamma["central_deg"]+gamma["plus"]+gamma["transport_deg"]]
    assert all(gamma_window[0]<m["gamma_mod_180_deg"]<gamma_window[1] for m in [complete,ordinary])
    r,t,v=(data[key]["central"] for key in ("ub_over_cb","td_over_ts","cb"))
    k=1-math.sqrt(1-v*v*(1+r*r))
    p,q=probability(r),probability(t)
    C2=(data["ud"]["central"]**2-k*(1-p-q)-k*k*p*q)/(1-k)
    tilted=math.acos(math.sqrt(C2))
    repair=build(tilted,r,t,v,1,1,data)
    A=models["1"]["two_widths"]["short_arc_tree_only_abs_shift_upper"]
    max_ud=marginal(data["ud"],2)[1]
    min_angle=math.degrees(math.acos(max_ud + (1-math.sqrt(1-marginal(data["cb"],2)[1]**2*
                      (1+marginal(data["ub_over_cb"],2)[1]**2)))*math.sqrt(probability(marginal(data["ub_over_cb"],2)[1])))-theta)
    # Independently transported nuclear extractions using the same historical kaon box.
    beta_rows={}
    for name,ud_input in (("2020",data["ud"]),("2025",data["reanalysis_2025_ud"])):
        ud_max=marginal(ud_input,2)[1]
        beta_rows[name]=math.sqrt(max(0.,1-ud_max**2-krange[1]**2))
    return {"general_unitary_minimum_Vub_from_nucleus_kaon":beta_rows,
            "max_Vub_from_two_heavy_tree_inputs":marginal(data["ub_over_cb"],2)[1]*marginal(data["cb"],2)[1],
            "D14_theta_deg":math.degrees(theta),"models":models,"central_extra_light_rotation_deg":math.degrees(tilted-theta),
            "central_short_with_extra_light_rotation":repair,"complete_two_width_long_example":complete,
            "ordinary_two_width_long_example":ordinary,"Kloe_two_width_Vus":krange,"gamma_working_window":gamma_window,
            "short_tree_only_minimum_beta_gap_two_width":math.cos(theta)-A-max_ud,
            "necessary_light_rotation_degrees_short_tree_only":min_angle-1e-10}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for key in a:
            compare(a[key],b[key])
    elif isinstance(a,(tuple,list)):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            compare(x,y)
    else:
        assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-12),(a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    inputs=json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig"))
    result=run(inputs)
    path=HERE/"results.json"
    if args.check:
        compare(result,json.loads(path.read_text(encoding="utf-8-sig")))
        print("1614 OK: both arcs, exact reconstruction, nuclear, kaon and charged-B witnesses")
    else:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,indent=2)
            handle.write("\n")
