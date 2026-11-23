"""Two charged masses determine the third mass and CKM phase with a measured port rule."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
read=lambda key:json.loads((HERE/CFG[key]).read_bytes())
MODEL=runpy.run_path(str(HERE/CFG["flavor_graph"]))
FLOW=MODEL["fl"]
PHI=math.pi*read("vacuum")["scalar_angle_over_pi"]
LEPTONS=np.array(read("weak_masses")["yukawa"]["charged_leptons"])
DOWN=np.array(read("weak_masses")["yukawa"]["down"])
DIRECT=read("weak_moduli")["direct_semileptonic_moduli"]
PRIOR=read("quark_points")["points"]["interior_all_one_width"]
NPOLES=MODEL["spectra"]["neutral_poles_GeV"]
MOMENT=MODEL["spectra"]["neutral_only_Higgs_mass_squared_GeV2"]


def boundary(q, e, mu, gamma, sign, rule):
    edges,_=FLOW["texture"]["inputs_at_gamma"](gamma,q)
    a,b,c,d32,d33,re,im=np.array(edges)*q[2]
    d23=re+1j*im
    ud=math.cos(q[3])*math.cos(q[4])
    z={"open":ud,"unit":1.,"closed":ud*ud}[rule]*np.exp(1j*PHI)
    K=a*a+b*b+9*c*c
    L=d32*d32+d33*d33+abs(d23)**2
    k=K*e*mu/(abs(z)*abs(a*b*d33))
    c0=e*e+mu*mu-abs(z)**2*L
    assert c0<0
    tau=(k+math.sqrt(k*k-4*c0))/2
    x=sign*math.sqrt(e*mu*tau/(abs(z)*abs(a*b*d33)))
    E=np.array([[0,x*a,0],[x*b,3*x*c,z*d23],[0,z*d32,z*d33]])
    target,vals=FLOW["basis"](E.T)
    # Trace and absolute determinant are satisfied before the remaining residual.
    assert abs(sum(vals**2)/(e*e+mu*mu+tau*tau)-1)<2e-13
    assert abs(np.prod(vals)/(e*mu*tau)-1)<2e-10
    D=np.array([[0,a,0],[b,c,d23],[0,d32,d33]])
    return math.log(vals[1]/mu), tau, x, z, target, np.linalg.svd(E/z-D,compute_uv=False)


def case(offsets, zero_return, seed, sign=1, step=.25, rule="open"):
    for f,name in (("Vus","us"),("Vub","ub"),("Vcb","cb")):
        FLOW["cfg"]["direct_semileptonic_moduli"][f][0]=DIRECT[f][0]+offsets.get(name,0)*np.linalg.norm(DIRECT[f][1:])
    ydown=DOWN[:,0]+np.array([offsets.get(n,0) for n in ("down","strange","bottom")])*DOWN[:,1]
    for i,d in enumerate(ydown):
        FLOW["low"]["yukawa"]["down"][i]=[float(d),0.]
    theta=FLOW["angles"]("direct")
    gamma=seed
    tau=CFG["blind_seeds"]["tau_MZ"]
    tan=4.7
    gauge_soft=[g[-1] for a,susy,_,g in FLOW["segments"](step) if not susy][-1]
    gy,gw=math.sqrt(.6*gauge_soft[0]),math.sqrt(gauge_soft[1])
    quartic_params={"gY":gy,"g2":gw,"D":(gy*gy+gw*gw)/8}
    gain=CFG["blind_seeds"]["neutral_gain"]
    frame=np.eye(3,dtype=complex)
    low_insert=None
    us=math.sin(theta[0])*math.cos(theta[1])
    phase0=np.exp(1j*PHI*np.array([2,0,1]))
    phase=np.exp(1j*np.array([0.,2*PHI,3*PHI+np.angle(3+2*np.exp(4j*PHI))]))/phase0
    upper=MODEL["sim"]["U"]
    v=MODEL["sim"]["vev"]
    for it in range(32):
        FLOW["low"]["yukawa"]["charged_leptons"][2]=[float(tau),0.]
        shape=MODEL["new_dirac"](2*us*us)
        yl,yh=((gain*phase[:,None]*shape)@upper.conj()).conj().T/v
        if low_insert is None:
            low_insert=yl
        FLOW["physics"]["Y0"]=np.column_stack((low_insert,0*yl,yh,0*yh))
        # Integrate below the high cluster again to close the two physical
        # insertions and optional u-return without any measured tau as input.
        y=np.zeros((4,3,4),complex)
        y[0,:,:3]=np.diag(FLOW["central"]("up"))
        y[1,:,:3]=FLOW["low_ckm"](gamma,theta)@np.diag(FLOW["central"]("down"))
        y[2,:,:3]=frame@np.diag(FLOW["central"]("charged_leptons"))
        beta=FLOW["beta"]
        for start,susy,dt,gauge in FLOW["segments"](step):
            if start==NPOLES[2]:
                break
            if start==NPOLES[0]:
                y[3,:,0]=low_insert
            for i in range(0,len(gauge)-1,2):
                g=gauge[i+1]
                k1=beta(y,gauge[i],False)
                k2=beta(y+dt*k1/2,g,False)
                k3=beta(y+dt*k2/2,g,False)
                k4=beta(y+dt*k3,gauge[i+2],False)
                y+=dt*(k1+2*k2+2*k3+k4)/6
        lu,_=FLOW["basis"](y[0])
        ld,sd=FLOW["basis"](y[1])
        p=abs(lu[:,0].conj()@ld)**2
        new_us=math.sqrt(p[1])
        new_d_factor=(p[1]*sd[1]-p[2]*sd[2])/(p[0]*sd[0]) if zero_return else 1.
        low_err=yl-y[3,:,0]
        m=MODEL["sim"]["m"]
        moment=(m[0]**2*np.linalg.norm(low_insert)**2+m[1]**2*np.linalg.norm(yh)**2)/(8*math.pi**2)
        gg=math.sqrt(MOMENT/moment)
        # Do not import the old fitted tau through the soft Higgs projection.
        FLOW["shoot"].__globals__["TAN"]=tan
        high,soft,lam=FLOW["shoot"](gamma,theta,frame,step)
        quartic_params["yt"]=FLOW["basis"](soft[0])[1][-1]
        phys=FLOW["physics"]
        next_tan=phys["bisect"](lambda t:phys["finite_quartic"](quartic_params,t,0.)-lam,2.,9.)
        q,l,ga,actual,_=FLOW["extraction"](high)
        residual,hi_tau,x,z,target,arm_values=boundary(q,l[0],l[1],ga,sign,rule)
        mapped=FLOW["mapped_frame"](frame,actual,target)
        tau_factor=hi_tau/l[2]
        error=max(abs(residual),abs(tau_factor-1),abs(next_tan/tan-1),np.linalg.norm(mapped-frame),
                  abs(new_us-us),np.linalg.norm(low_err)/np.linalg.norm(yl),abs(gg-1),abs(new_d_factor-1))
        if error<2e-10:
            break
        deriv=(boundary(q,l[0],l[1],ga+.002,sign,rule)[0]-boundary(q,l[0],l[1],ga-.002,sign,rule)[0])/.004
        gamma-=np.clip(residual/deriv,-16.,16.)
        tau*=tau_factor
        frame=mapped
        low_insert=(low_insert+low_err)*gg
        gain*=gg
        us=new_us
        tan=next_tan
        FLOW["low"]["yukawa"]["down"][0][0]*=new_d_factor
    else:
        raise RuntimeError(("three-source shooting did not close",offsets,sign,seed,error,gamma))
    return {"gamma_MZ_deg":gamma,"gamma_10PeV_deg":ga,"tau_MZ":tau,"tau_10PeV":hi_tau,
            "Vud_10PeV":math.cos(q[3])*math.cos(q[4]),"abs_z":abs(z),"x":x,"Vus_at_neutral_upper":us,"moment_GeV2":moment,
            "MZ_down_yukawas":FLOW["central"]("down").tolist(),"high_down_yukawas":q[:3].tolist(),
            "high_e_mu":l[:2].tolist(),"return_relative":(p[1]*sd[1]-p[2]*sd[2]-p[0]*sd[0])/(p[0]*sd[0]),
            "positive_zero_stop_tan_beta":tan,"extra_eff_vertex_singular_values":arm_values.tolist(),
            "iterations":it+1,"shoot_residual":error}


def observed(row):
    gamma=read("direct_gamma")["time_integrated_gamma_deg"]
    tau,sigma=LEPTONS[2]
    return {"tau_marginal_widths":(row["tau_MZ"]-tau)/sigma,
            "in_time_integrated_gamma_68":gamma["68"][0]<=row["gamma_MZ_deg"]<=gamma["68"][1],
            "in_time_integrated_gamma_95":gamma["95"][0]<=row["gamma_MZ_deg"]<=gamma["95"][1]}


def exchange_tau_for_bottom(offsets,zero_return,seed,step,rule="open"):
    # Choose the initial bracket in absolute Yukawa units. The measured bottom's
    # center and standard deviation are used solely to encode input/output units.
    low,high=0.014,0.019
    def trace(yb):
        shift=(yb-DOWN[2,0])/DOWN[2,1]
        return case({**offsets,"bottom":shift},zero_return,seed,step=step,rule=rule)
    tau=LEPTONS[2,0]
    r0,r1=trace(low),trace(high)
    assert r0["tau_MZ"]<tau<r1["tau_MZ"]
    for _ in range(8):
        yb=low+(high-low)*(tau-r0["tau_MZ"])/(r1["tau_MZ"]-r0["tau_MZ"])
        out=trace(yb)
        if abs(out["tau_MZ"]/tau-1)<2e-10:
            break
        if out["tau_MZ"]<tau:
            low,r0=yb,out
        else:
            high,r1=yb,out
    else:
        raise RuntimeError("bottom from physical tau not closed")
    out["measured_bottom_widths"]=(yb-DOWN[2,0])/DOWN[2,1]
    out["measured_down_widths"]=(out["MZ_down_yukawas"][0]-DOWN[0,0])/DOWN[0,1]
    return out


def main():
    seed=CFG["blind_seeds"]["gamma_deg"]
    step=CFG["log_step"]
    cases={}
    for name,p,zero in (("all_quark_central",{},False),
                        ("registered_moduli_zero_return",PRIOR,True),
                        ("central_moduli_zero_return",{},True)):
        cases[name]=case(p,zero,seed,step=step)
    high=case({},False,seed,step=step/2)
    finer={k:abs(cases["all_quark_central"][k]-high[k]) for k in ("gamma_MZ_deg","tau_MZ","x")}
    assert finer["gamma_MZ_deg"]<1e-6 and finer["tau_MZ"]<2e-10
    alternatives=[]
    for sign,other in ((1,-110.),(-1,-135.),(-1,150.)):
        row=case({},False,other,sign,step)
        alternatives.append({"sign_x":sign,"seed_deg":other,"gamma_deg":row["gamma_MZ_deg"],
                             "tau_MZ":row["tau_MZ"]})
    reverse={k:exchange_tau_for_bottom(p,zero,seed,step)
             for k,p,zero in (("central_except_bottom",{},False),
                              ("registered_moduli_zero_return",PRIOR,True))}
    assert all(abs(row["measured_bottom_widths"])<1 and observed(row)["in_time_integrated_gamma_68"]
               for row in reverse.values())
    b_variation=[]
    for shift in (-1.,1.):
        row=case({"bottom":shift},False,seed,step=step)
        b_variation.append({"bottom_widths":shift,"tau_MZ":row["tau_MZ"],"gamma_MZ_deg":row["gamma_MZ_deg"]})
    assert b_variation[0]["tau_MZ"]<LEPTONS[2,0]<b_variation[1]["tau_MZ"]
    competitors={}
    for key,points,zero in (("central",{},False),("registered_zero_return",PRIOR,True)):
        competitors[key]={}
        for rule in ("unit","closed"):
            row=exchange_tau_for_bottom(points,zero,seed,step,rule)
            competitors[key][rule]={**{k:row[k] for k in ("gamma_MZ_deg","tau_MZ","abs_z","measured_bottom_widths")},
                                    "bottom_MZ":row["MZ_down_yukawas"][-1]}
    # Original high point only, not an input to any of the new mass solutions.
    frozen=json.loads((HERE/"../1723/results.json").read_bytes())["original_10PeV_masses_down_tau_desc"]
    angles=runpy.run_path(str(HERE/"../1715/check.py"))["charged_background"](.10)[1]["CKM_angles_at_10PeV_rad"]
    vtb=math.cos(angles[1])*math.cos(angles[2])
    ceiling=frozen[1][0]/(vtb*frozen[0][0])
    assert ceiling<1
    return {"cases":cases,"holdouts":{k:observed(v) for k,v in cases.items()},
            "calibrate_tau_hold_out_bottom_and_gamma":reverse,
            "same_inputs_competing_port_rules":competitors,
            "frozen_quark_spectrum_third_row_bound_for_all_x_and_gamma":{"Vtb_high":vtb,"abs_z_upper":ceiling},
            "independent_bottom_endpoints_at_other_quark_centers":b_variation,
            "other_signed_or_CP_branches_central":alternatives,"halfstep_differences":finer}


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [tidy(a) for a in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.11g}")
    if isinstance(x,np.bool_): return bool(x)
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    results=tidy(main())
    output=HERE/"results.json"
    if check:
        assert results==json.loads(output.read_bytes())
        print("1727 verified two-mass shooting, recomputed neutral columns and blind tau/gamma holdouts")
    else:
        with output.open("w",encoding="utf8") as out:
            json.dump(results,out,indent=2,allow_nan=False)
            out.write("\n")
