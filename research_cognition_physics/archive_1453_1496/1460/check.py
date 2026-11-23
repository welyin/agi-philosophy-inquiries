"""Rate-only sharp row contraction and exact one-doublet UV construction."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(p):return json.loads(p.read_text(encoding="utf-8"))


def number(k):
    val=float(k["mean"])
    err=float(k["sigma"]) if "sigma" in k else math.hypot(*map(float,k["sigma_components"]))
    return val,err,float(k["rounding_half_unit"])


def axial(x,y,r,S):
    den=1+r*r; b=x+r*y
    return (-b+math.sqrt(b*b+den*(4*S-x*x-y*y)))/den


def product_floor(x,y,r0,r1,S):
    points=[]
    for r in (r0,r1):
        t=axial(x,y,r,S)
        points.append(((t-x)/2,(r*t-y)/2))
    if min(z for pair in points for z in pair)<=0:return 0.
    return min(a*b for a,b in points)


def basis(bottom):
    # Columns adjacent to the supplied last row span its kernel.
    a=np.asarray(bottom,dtype=float)
    q=np.linalg.svd(a[None,:],full_matrices=True)[2][1:,:]
    return np.vstack([q,a])


def witness(x,y,r,z,w):
    fu=w["Fu"]; mu,mc,mt=w["up_light_masses"]; md,ms,mb=w["down_light_masses"]; B=w["B_GeV"]
    a=np.array([math.sqrt(fu),0.,0.,math.sqrt(1-fu)])
    def directions(t):
        bd,bs=(x-t)/(2*a[0]),(y-r*t)/(2*a[0])
        b=np.array([bd,bs,0.,math.sqrt(1-bd*bd-bs*bs)])
        ld,ls=(x+t)/2,(y+r*t)/2
        lB=(mu*a[0]-md*ld*bd-ms*ls*bs)/(B*b[3])
        return b,np.array([ld,ls,z,lB])
    ref=axial(x,y,r,1-z*z)
    lo,hi=ref-1e-7,ref+1e-7
    for _ in range(65):
        mid=(lo+hi)/2
        if np.linalg.norm(directions(mid)[1])<1:lo=mid
        else:hi=mid
    b,row=directions((lo+hi)/2)
    d=np.array([md,ms,mb,B])
    va=d*b
    other=va-mu*a[0]*row
    T=np.linalg.norm(other)/a[3]
    rowT=other/np.linalg.norm(other)
    inner=np.linalg.svd(np.vstack([row,rowT]),full_matrices=True)[2][2:,:]
    vl=np.vstack([row,inner,rowT]); vr=np.outer(a,b)
    massu=np.diag([mu,mc,mt,T])@basis(a).T
    massd=vl@np.diag(d)@basis(b).T
    assert np.allclose(vl@vl.T,np.eye(4),atol=3e-12)
    assert np.max(abs(massu[:,3]-massd[:,3]))<3e-10
    assert np.allclose([abs(row[0]+vr[0,0]),abs(row[1]+vr[0,1]),abs((row[1]-vr[0,1])/(row[0]-vr[0,0]))],[x,y,r])
    fd=b[0]**2
    return {"V_left_u_row":row.tolist(),"right_doublet_components_up":a.tolist(),"right_doublet_components_down":b.tolist(),
       "T_B_masses_GeV":[T,B],"same_mass_column_residual_GeV":float(np.max(abs(massu[:,3]-massd[:,3]))),
       "delta_Qp":2*fu-fd,"delta_QCs":188*fu-211*fd,"abs_Fds":abs(b[0]*b[1]),"abs_RdRs":abs(vr[0,0]*vr[0,1])}


def run():
    p=read(HERE/"inputs.json"); c=read(HERE/p["charged"]); pv=read(HERE/p["parity"]); rare=read(HERE/p["rare"])
    x,y,r,z=[number(c[k]) for k in ("x","y","r","z")]
    A=p["mass_scope_GeV"]; bB=p["heavy_right_fraction_amplitude_lower"]
    def row(w,rounding=True):
        def box(t):
            v,e,ro=t;spread=w*e+(ro if rounding else 0)
            return v-spread,v+spread
        xu=box(x)[1];yu=box(y)[1];rr=box(r);zu=box(z)[1]/math.sqrt(p["B_inclusive_min_left_rate_response"])
        leaks=(A["up_max"]+math.sqrt(1-bB*bB)*math.hypot(A["down_max"],A["strange_max"],A["bottom_max"]*zu))/(A["B_min"]*bB)
        S=1-zu*zu-leaks*leaks
        return {"work_width":w,"Vub_left_upper":zu,"VuB_left_upper":leaks,"light_ds_left_square_lower":S,
                "abs_RdRs_lower":product_floor(xu,yu,*rr,S)}
    rows=[row(w,rounding=w>0) for w in p["charged_work_width"]]
    proton=pv["proton"];cs=pv["caesium"]; ww=p["parity_width"]
    qp=proton["pv_electron_scattering_only"];bp=proton["SM"]
    qc=cs["value_and_errors_experiment_atom"];bc=cs["SM"]
    p_hi=qp[0]-bp[0]+ww*math.hypot(qp[1],bp[1])
    c_lo=qc[0]-bc[0]-ww*math.hypot(qc[1],qc[2],bc[1])
    fu=(211*p_hi-c_lo)/234
    wrare=p["rare_width"];be=rare["B_e3"]; bn=rare["B_rare"];other=rare["other_amplitude"]
    bu=bn["central"]+wrare*bn["plus"]
    be_low=be["central"]-wrare*(be["stat"]+be["sys"])
    yn=y[0]+2.5*y[1]+y[2]
    C=3*be_low*min(rare["r_conversion"]["sensitivity"])/(2*yn**2)
    b0=other["SM_rate_anchor"]+wrare*other["SM_width"]
    Pmax=(math.sqrt(bu)+math.sqrt(b0))/math.sqrt(C)
    for v in rows:
        pmin=v["abs_RdRs_lower"]/fu
        v.update({"abs_Fds_lower_using_PV":pmin,"min_other_equivalent_BR":max(0.,math.sqrt(C)*pmin-math.sqrt(bu))**2})
    ex=witness(x[0],y[0],r[0],z[0],p["full_witness"])
    ex["min_other_equivalent_BR"]=max(0.,math.sqrt(C)*ex["abs_Fds"]-math.sqrt(bu))**2
    assert rows[1]["abs_Fds_lower_using_PV"]>Pmax and rows[2]["abs_RdRs_lower"]==0
    # Rigorous existence of a zero-Rd rate solution at two widths (take t=x).
    def zero_rd(w):
        xx=x[0]+w*x[1]; yy=y[0]+w*y[1]; rr=r[0]+w*r[1];zz=z[0]
        ls=math.sqrt(1-xx**2-zz**2)
        cos=(4*ls**2-yy**2-(rr*xx)**2)/(2*yy*rr*xx)
        assert -1<=cos<=1
        return {"x":xx,"y":yy,"r":rr,"abs_Ls":ls,"relative_vector_axial_cos":cos}
    pure_d=.974471186518; pure_s=.224475158195; pure_z=math.sqrt(1-pure_d**2-pure_s**2)
    assert all(abs(val-param[0])<=2.5*param[1]+param[2] for val,param in zip([pure_d,pure_s,pure_s/pure_d,pure_z],[x,y,r,z]))
    return {"charged_windows":rows,"weak_charges_Fu_upper":fu,"rare_C_lower":C,"rare_SM_other_conditional_budget":b0,
        "rare_abs_Fds_upper":Pmax,"joint_max_abs_RdRs":fu*Pmax,"full_UV_example":ex,"two_width_current_only_zero_Rd":zero_rd(2),
        "two_point_five_width_exact_decoupled_UV_left_row":[pure_d,pure_s,pure_z,0.]}


if __name__=="__main__":
    p=argparse.ArgumentParser(); g=p.add_mutually_exclusive_group()
    g.add_argument("--write",action="store_true");g.add_argument("--check",action="store_true")
    args=p.parse_args();out=json.dumps(run(),indent=2,allow_nan=False)+"\n";dest=HERE/"results.json"
    if args.write:
        with dest.open("x",encoding="utf-8") as f:f.write(out)
    elif dest.read_text(encoding="utf-8")!=out:raise SystemExit("1460 mismatch")
    else:print("1460 OK")
