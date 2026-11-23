"""Transfer physical quark and charged-lepton solutions to the two-source Majorana problem."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
def read(p): return json.loads((HERE/p).read_bytes())
def load(p): return runpy.run_path(str(HERE/p))
CFG=read("inputs.json")
TRAN=load(CFG["transport"])
NU=load(CFG["neutral"])
FREQ=read(CFG["frequency"])
SOL=FREQ["juno_NO_dm21"][0]
ANG=np.array(read(CFG["angles"])["joint_windows"])
OLD=read(CFG["old_charged_output"])
FROZEN=sorted(r["raw_real_b_coefficient"] for r in
              read(CFG["old_neutral"])["original_basis_verification"]["direct_forward_both_roots"])

def endpoint(name,step=.25,explicit_offsets=None):
    """Return charged physical data and an unscaled complex-symmetric MZ pencil C+rho*b*b.T."""
    charged=load(CFG["charged"])
    flow=charged["FLOW"]
    offsets=charged["PRIOR"].copy() if name.startswith("registered") else {}
    zero=name.startswith("registered")
    forward=name.startswith("bottom_")
    rule="open" if forward else name.split("_")[1]
    # Reuse a known scalar root as a starting input; reconstruct the complete frame and
    # verify the physical tau equation to the original shooting tolerance.
    if explicit_offsets is not None:
        offsets={key:explicit_offsets[key] for key in ("us","ub","cb","strange","bottom")}
        forward,zero,rule=False,True,"open"
    elif forward:
        offsets["bottom"]=-1 if "minus" in name else 1
    else:
        if rule=="open":
            key="registered_moduli_zero_return" if zero else "central_except_bottom"
            old=OLD["calibrate_tau_hold_out_bottom_and_gamma"][key]
        else:
            key="registered_zero_return" if zero else "central"
            old=OLD["same_inputs_competing_port_rules"][key][rule]
        offsets["bottom"]=old["measured_bottom_widths"]
    original=flow["mapped_frame"]
    cache={}
    def remember(frame,actual,target):
        cache["low"]=frame.copy()
        cache["high"]=target.copy()
        return original(frame,actual,target)
    flow["mapped_frame"]=remember
    row=charged["case"](offsets,zero,59.,step=step,rule=rule)
    if not forward:
        assert abs(row["tau_MZ"]/charged["LEPTONS"][2,0]-1)<1e-8
    flow["mapped_frame"]=original
    frame=cache["low"]
    gamma=row["gamma_MZ_deg"]
    flow["TAN"]=row["positive_zero_stop_tan_beta"]
    high=flow["shoot"](gamma,flow["angles"]("direct"),frame,step)[0]
    q,l,g,_,A=flow["extraction"](high)
    edges=flow["texture"]["inputs_at_gamma"](g,q)[0]
    u,d,c,h,k,re,im=q[2]*np.array(edges)
    D=np.array([[0,u,0],[d,c,re+1j*im],[0,h,k]])
    z=row["abs_z"]*np.exp(1j*charged["PHI"])
    x=row["x"]
    E=np.array([[0,x*u,0],[x*d,3*x*c,z*(re+1j*im)],[0,z*h,z*k]])
    pack=flow["pack"]
    rep={"charged_left_MZ":pack(frame),"charged_left_10PeV":pack(cache["high"]),
         "gamma_MZ_deg":gamma,"gamma_10PeV_deg":g,"Y_down_at_10PeV":q[:3],
         "Y_charged_at_10PeV":l,"CKM_angles_at_10PeV_rad":q[3:],
         "neutral_columns_10PeV":pack(high[3]),"up_texture_diag_10PeV":A}
    T,cert=TRAN["FLOW"]["congruence"](flow,rep,step)
    high_N=high[3]
    v=np.linalg.svd(high_N)[0][:,-1]
    n=np.linalg.solve(T,v)
    n/=np.linalg.norm(n)
    labels=TRAN["read"](TRAN["CFG"]["family_labels"])
    ab=np.transpose([labels["source_a"],labels["source_b"]])
    phi=charged["PHI"]
    phase=np.exp(1j*np.array([0,2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))]))
    dc=TRAN["at_cluster"](flow,rep,step).conj()
    coeff=np.linalg.lstsq(dc,phase[:,None]*ab,rcond=None)[0]
    lift=np.linalg.norm(dc@coeff-phase[:,None]*ab)
    assert lift<1e-9
    b=(T.T@high_N[:,[0,2]].conj()@coeff)[:,1]
    W=D@T
    def base(powers):
        ir=np.exp(4j*phi*np.array(powers))
        c=W.T@np.diag(1/(ir*A))@W
        return c-np.outer(c@n,c@n)/(n@c@n)
    C=base(CFG["up_monomials"])
    assert np.linalg.norm(C@n)<1e-11 and abs(n@b)<1e-10
    color=TRAN["read"](TRAN["CFG"]["flavor_vacuum"])
    s=color["scalar_magnitude_over_cutoff"]*np.exp(1j*phi)
    P10=np.diag(np.exp(1j*phi*np.array([0,0,1])))
    P5=np.diag(np.exp(1j*phi*np.array([2,0,1])))
    detheavy=np.linalg.det([[0,s**3],[-s**3,1]])
    strong=np.angle(np.linalg.det(P10@D@P5)*np.linalg.det(P10@np.diag(A*np.exp(4j*phi*np.array(CFG["up_monomials"])))@P10)*detheavy)
    assert abs(strong)<2e-12
    return dict(C=C,b=b,n=n,E=E,D=D,A=A,T=T,high_N=high_N,low_E=frame,rep=rep,flow=flow,
                original_C=base([0,0,0]),strong=strong,row=row,zero=zero,rule=rule,
                dyb=float(offsets["bottom"]),tau_width=(row["tau_MZ"]-charged["LEPTONS"][2,0])/charged["LEPTONS"][2,1],
                lift=lift,cert=cert)

def inside_box(value,table):
    x,d,u=table
    return bool(x-d<=value<=x+u)

def clocks(C,b,n):
    Q=np.linalg.svd(n[:,None],full_matrices=True)[0][:,1:]
    block=Q.T@C@Q
    arm=Q.T@b
    norm=np.linalg.norm(block,2)
    w=np.linalg.norm(arm)**2
    O=block/norm
    V=np.outer(arm,arm)/w
    r=math.sqrt(FREQ["juno_NO_dm31"][0]/SOL)
    F=np.polynomial.Polynomial([np.vdot(O,O).real,2*np.vdot(O,V).real,np.vdot(V,V).real])
    det=np.polynomial.Polynomial([np.linalg.det(O),O[0,0]*V[1,1]+O[1,1]*V[0,0]-2*O[0,1]*V[1,0]])
    det_squared=det*np.polynomial.Polynomial(det.coef.conj())
    eq=F*F-(r+1/r)**2*np.polynomial.Polynomial(det_squared.coef.real)
    roots=sorted(x*norm/w for x in NU["real_roots"](eq))
    return roots

def observables(C,b,rho):
    M=C+rho*np.outer(b,b)
    U,s,_=np.linalg.svd(M)
    U=U[:,::-1].conj()
    p=abs(U)**2
    z=p[0,2]
    a,t=p[0,1]/(1-z),p[1,2]/(1-z)
    J=np.imag(U[0,0]*U[1,1]*np.conj(U[0,1]*U[1,0]))
    den=math.sqrt(a*(1-a)*t*(1-t)*z)
    delta=math.atan2(J/((1-z)*den),(p[1,0]-a*(1-t)-(1-a)*t*z)/(2*den))
    scale=math.sqrt(SOL)/s[1]
    masses=s[::-1]*scale
    dm31=masses[2]**2
    ang=np.array([4*z*(1-z),a,t])
    assert s[2]<1e-10*s[1]
    return {"rho":rho,"scale_eV":scale,"dm31_eV2":dm31,"dm32_eV2":dm31-SOL,"angle_readings":ang,
            "angles_in_box":[bool(x) for x in (ANG[:,0]<=ang)&(ang<=ANG[:,1])],"J":J,"delta_deg":np.degrees(delta),
            "SK_T2K_FC68_90":NU["fc_read"](t,delta),"dm31_and_32_in_box_at_central_solar":[
                inside_box(dm31,FREQ["juno_NO_dm31"]),inside_box(dm31-SOL,FREQ["accelerator_NO_dm32"])],
            "mbb_meV":abs(M[0,0])*scale*1000,"mbeta_meV":math.sqrt(p[0]@(masses**2))*1000,"mass_sum_meV":sum(masses)*1000}

def absolute_window(C,b,intervals):
    """Outer bounds over all allowed raw gains AND all three independent frequency windows."""
    if not intervals:
        return None
    boxes=[NU["bounds"](FREQ[k]) for k in ("juno_NO_dm21","juno_NO_dm31","accelerator_NO_dm32")]
    B=np.outer(b,b)
    hulls=[[],[],[]]
    for start,end in intervals:
        grid=np.linspace(start,end,129)
        for l,r in zip(grid,grid[1:]):
            M=C+(.5*l+.5*r)*B
            s3,s2,_=np.linalg.svd(M,compute_uv=False)
            d=.5*(r-l)
            err=d*np.linalg.norm(B,2)+5e-15
            rlo,rhi=(s3-err)/(s2+err),(s3+err)/(s2-err)
            zlo=max(boxes[0][0],boxes[1][0]/rhi**2,boxes[2][0]/(rhi**2-1))
            zhi=min(boxes[0][1],boxes[1][1]/rlo**2,boxes[2][1]/(rlo**2-1))
            if zhi<zlo:
                continue
            scl=1000*math.sqrt(zlo)/(s2+err)
            sch=1000*math.sqrt(zhi)/(s2-err)
            for h,norm,rate in ((hulls[0],abs(M[0,0]),abs(B[0,0])),
                                (hulls[1],np.linalg.norm(M[:,0]),np.linalg.norm(B[:,0]))):
                h.append((max(0,norm-d*rate-5e-15)*scl,(norm+d*rate+5e-15)*sch))
            hulls[2].append((1000*math.sqrt(zlo)*(1+rlo),1000*math.sqrt(zhi)*(1+rhi)))
    return {key:[min(x[0] for x in h),max(x[1] for x in h)] for key,h in
            zip(("mbb_meV","mbeta_meV","sum_meV"),hulls)}

def compute(name,step):
    pt=endpoint(name,step)
    C,b,n=pt["C"],pt["b"],pt["n"]
    def fit(K):
        roots=clocks(K,b,n)
        rows=[observables(K,b,x) for x in roots]
        assert all(abs(row["dm31_eV2"]/FREQ["juno_NO_dm31"][0]-1)<1e-8 for row in rows)
        return rows
    row=pt["row"]
    Q=np.linalg.svd(n[:,None],full_matrices=True)[0][:,1:]
    W=pt["D"]@pt["T"]
    def all_gains(k,psi,sgn):
        res=NU["shape"](n,Q,W,pt["A"],b,k,psi)
        def raw(interval):
            return sorted(sgn*e*res["rho_per_x"] for e in interval)
        return {"allowed_by_clocks_only":[raw(i) for i in res["frequency_intervals"]],
                "all_clock_and_angle_intervals":[raw(i["interval"]) for i in res["accepted"]],
                "when_angles_pass_FC68_90":res["FC_if_angles_pass"],
                "no_angle_and_clock_solution":res["full_rejection_certificate"],
                "reactor_only_rejection":res["reactor_only_certificate"]}
    selected=all_gains((4,2),np.pi/3,-1)
    return {"charged":{"gamma_MZ_deg":row["gamma_MZ_deg"],"yb_width":pt["dyb"],"tau_width":pt["tau_width"],
                       "tan_beta":row["positive_zero_stop_tan_beta"],"moment_GeV2":row["moment_GeV2"],
                       "lift_residual":pt["lift"],"strong_angle":pt["strong"],"up_10PeV":pt["A"],
                       "down_10PeV":row["high_down_yukawas"],"null_probabilities":abs(n)**2},
            "two_historical_raw_gains_solar_only":[observables(C,b,rho) for rho in FROZEN],
            "two_clock_roots":fit(C),"original_up_phases_two_clock_roots":fit(pt["original_C"]),
            "all_gains_selected":selected,
            "absolute_mass_outer_box_all_allowed_gains_meV":absolute_window(C,b,selected["all_clock_and_angle_intervals"]),
            "all_gains_original":all_gains((0,0),0,1)}

def main():
    step=CFG["step"]
    cases={name:compute(name,step) for name in CFG["cases"]}
    finer=compute("registered_open",step/2)
    diff=max(abs(x["delta_deg"]-y["delta_deg"]) for x,y in zip(
        cases["registered_open"]["two_clock_roots"],finer["two_clock_roots"]))
    assert diff<1e-5
    return {"old_raw_gains":FROZEN,"solar_calibration_eV2":SOL,"cases":cases,"halfstep_delta_difference_deg":diff}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check=parser.parse_args().check
    result=TRAN["FLOW"]["tidy"](main())
    path=HERE/"results.json"
    if check:
        assert result==json.loads(path.read_bytes())
        print("1731 charged endpoints, physical lifts, two-source fits and holdouts verified")
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
