"""Transport two independent real LNV sources as complete pseudo-Dirac clusters."""
import argparse
import json
import math
from itertools import combinations
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial

HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())

cfg=read("inputs.json")
fl=runpy.run_path(str(HERE/cfg["flavour"]))
zer=runpy.run_path(str(HERE/cfg["zero_lines"]))
cal=read(cfg["calibrating_frequencies"])
src=read("../1698/results.json")["shared_witness"]
base=runpy.run_path(str(HERE/"../1698/check.py"))
D0,C0,H0=base["neutral_graph"](read(cfg["threshold_graph"]))
b,t=src["neutral_b"],src["neutral_t"]
D0,C0,H0=math.sqrt(b)*D0,b*math.sqrt(t)*C0,b*t*H0
U,m,Vh=np.linalg.svd(C0)
U,m,V=U[:,::-1],m[::-1],Vh.conj().T[:,::-1]
vev=read("../1492/inputs.json")["v246_GeV"]/math.sqrt(2)
ANG=np.array(read(cfg["oscillation_angles"])["joint_windows"])
freq=np.array([cal["juno_NO_dm21"][0],cal["juno_NO_dm31"][0]])


def graph(k,dz,role,step):
    ref=read(cfg["endpoints"])
    anchor=(ref["down_one_width_lower"]["direct_k"+str(k)] if dz else
            ref["charged_endpoints"]["direct"][str(k)])
    angle=-np.pi/6
    orig=np.exp(1j*angle*np.array([2,0,1]))
    port=np.exp(1j*np.array([0,2*angle,3*angle+np.angle(3+2*np.exp(4j*angle))]))/orig
    third=zer["kernel"](anchor["large_modulus"]) if role=="shared_third" else 1.
    prior_low=None
    common_gain=1.
    target_moment=src["neutral_only_Higgs_mass_squared_GeV2"]
    for _ in range(18):
        d=common_gain*(port*np.array([1,1,third]))[:,None]*D0
        dirac=d@U.conj()
        yl,yh=dirac.conj().T/vev
        if prior_low is None:
            prior_low=yl.copy()
        fl["physics"]["Y0"]=np.column_stack([prior_low,0*yl,yh,0*yh])
        rep=fl["solve"](k,"direct",seed=anchor["gamma_MZ_deg"],step=step,d_z=dz)
        trace,low_raised=background(rep,dz,step)
        next_h=zer["kernel"](rep["large_modulus"]) if role=="shared_third" else 1.
        err=np.linalg.norm(low_raised-yl)/np.linalg.norm(yl)
        moment=(m[0]**2*np.linalg.norm(prior_low)**2+m[1]**2*np.linalg.norm(yh)**2)/(8*np.pi**2)
        gain_update=math.sqrt(target_moment/moment)
        if err<3e-11 and abs(next_h-third)<2e-11 and abs(gain_update-1)<1e-11:
            break
        prior_low=(prior_low+yl-low_raised)*gain_update
        common_gain*=gain_update
        third=next_h
    else:
        raise RuntimeError("common-scale neutral/charged boundary did not converge")
    common={"dimensionless_common_D_gain":common_gain,"pole_Yukawa_Higgs_moment_GeV2":moment,
            "heavy_weak_weight_leading":np.linalg.norm(dirac/m)**2}
    return dirac,rep,trace,third,err,np.linalg.norm(prior_low-yl)/np.linalg.norm(yl),common


def background(rep,dz,step):
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(fl["central"]("up"))
    down=fl["central"]("down")
    down[0]+=dz*fl["low"]["yukawa"]["down"][0][1]
    y[1,:,:3]=fl["low_ckm"](rep["gamma_MZ_deg"],fl["angles"]("direct"))@np.diag(down)
    y[2,:,:3]=zer["unpack"](rep["charged_left_MZ"])@np.diag(fl["central"]("charged_leptons"))
    hi,lo=src["neutral_poles_GeV"][2],src["neutral_poles_GeV"][0]
    lam=fl["physics"]["quartic"]["lambda_paper"][0]/4
    track=[]
    def stage(y,lg,g):
        h=y@y.conj().transpose(0,2,1)
        E,N=h[2:]
        v=np.array([3,3,1,1])
        T=v@np.trace(h,axis1=1,axis2=2).real
        gy,gw=.6*g[0],g[1]
        T2=v@np.trace(h@h,axis1=1,axis2=2).real
        dl=24*lg**2+(4*T-9*gw-3*gy)*lg-2*T2+3/8*(2*gw**2+(gw+gy)**2)
        return dl,(-1.5*E+.5*N,-1.5*E+3.5*N,
                   2*T-1.5*gy-4.5*gw,4*lg+2*T-3*gw)
    for a,susy,dt,gg in fl["segments"](step):
        if a>=hi:
            return track,y[3,:,0].copy()
        if a==lo:
            y[3,:,0]=fl["physics"]["Y0"][:,0]
        for j in range(0,len(gg)-1,2):
            f1=fl["beta"](y,gg[j],False)
            d1,c1=stage(y,lam,gg[j])
            yy=y+dt*f1/2
            f2=fl["beta"](yy,gg[j+1],False)
            d2,c2=stage(yy,lam+dt*d1/2,gg[j+1])
            yy=y+dt*f2/2
            f3=fl["beta"](yy,gg[j+1],False)
            d3,c3=stage(yy,lam+dt*d2/2,gg[j+1])
            yy=y+dt*f3
            f4=fl["beta"](yy,gg[j+2],False)
            d4,c4=stage(yy,lam+dt*d3,gg[j+2])
            y+=dt*(f1+2*f2+2*f3+f4)/6
            lam+=dt*(d1+2*d2+2*d3+d4)/6
            track.append((dt,a>=lo,(c4,c3,c2,c1)))
    raise ValueError("missing upper cluster")


def pair_sources(d):
    # Work to first order in the LNV spurion. The full Schur kappa is exact.
    # f=(0,1) selects both physical poles in the high cluster of M†M.
    # Its upper block has a finite off-diagonal Loewner divided difference.
    z=d/m
    p=m[0]**2/(m[0]**2-m[1]**2)
    split=np.array([[0.,p],[p,1.]])
    total,high=[],[]
    for h in np.eye(2):
        hs=V.T@(H0*np.diag(h))@V
        total.append(z@hs@z.T*1e9)
        high.append(z@(split*hs)@z.T*1e9)
        # Full spectral-projector contour evaluated with paired Schur complements.
        point=m[1]**2
        radius=(point-m[0]**2)/3
        res=0j
        D=np.diag(m)
        inv=np.diag(1/m)
        for a in 2*np.pi*(np.arange(48)+.5)/48:
            hstep=radius*np.exp(1j*a)
            arg=point+hstep
            P=arg*np.eye(2)-D@D
            B=P-hs.conj().T@hs
            E=D@hs
            K=np.linalg.solve(P-E@np.linalg.solve(B,E.conj().T),
                              -inv@hs@inv+E@np.linalg.solve(B,inv))
            res=res-d@K@d.T*(1e9*hstep/48)
        np.testing.assert_allclose(res,high[-1],rtol=5e-13,atol=1e-15)
    return np.array([total,high])


def lower_tensors(pair,trace):
    # Tensor conjugates here transform covariantly with the left charged
    # Yukawa. Curves carry total kappa and the matched high Weinberg term.
    y=pair.conj().copy()
    for dt,between,ss in reversed(trace):
        def rhs(val,coef):
            P,Q,al,ak=coef
            K,W=val
            T=P@K+K@P.T
            if between:
                return np.array([T+al*K+(ak-al)*W,Q@W+W@Q.T+ak*W])
            return np.array([T+ak*K,np.zeros_like(W)])
        f1=rhs(y,ss[0])
        f2=rhs(y-dt*f1/2,ss[1])
        f3=rhs(y-dt*f2/2,ss[2])
        f4=rhs(y-dt*f3,ss[3])
        y-=dt*(f1+2*f2+2*f3+f4)/6
    return y[0].conj()


def ratio_roots(mats):
    # All solutions for positive real H6/H0 reweighting. Minors avoid
    # subtracting two nearly equal large eigenvalues.
    A,B=mats
    T=Polynomial([np.vdot(A,A).real,2*np.vdot(A,B).real,np.vdot(B,B).real])
    P=Polynomial([0.])
    for i,j in combinations(range(3),2):
        for k,l in combinations(range(3),2):
            el=lambda a,b:Polynomial([A[a,b],B[a,b]])
            mm=el(i,k)*el(j,l)-el(i,l)*el(j,k)
            P+=Polynomial(np.conj(mm.coef))*mm
    poly=Polynomial((P-T*T*np.prod(freq)/sum(freq)**2).coef.real)
    roots=sorted(x.real for x in poly.roots() if abs(x.imag)<1e-9 and x.real>0.)
    signs=np.sign(poly.coef)
    assert len(roots)==2 and sum(signs[1:]!=signs[:-1])==2
    for x in roots:
        assert poly(x*(1-1e-5))*poly(x*(1+1e-5))<0
    return roots,poly.coef.tolist()


def observables(matrix,charged):
    uv,masses,_=np.linalg.svd(matrix)
    pmat=charged.conj().T@uv[:,::-1].conj()
    p=abs(pmat)**2
    s12,s23,s13=p[0,1]/(1-p[0,2]),p[1,2]/(1-p[0,2]),p[0,2]
    J=(pmat[0,0]*pmat[1,1]*np.conj(pmat[0,1]*pmat[1,0])).imag
    sin=J/((1-s13)*math.sqrt(s12*(1-s12)*s23*(1-s23)*s13))
    cos=(p[1,0]-s12*(1-s23)-(1-s12)*s23*s13)/(2*math.sqrt(s12*(1-s12)*s23*(1-s23)*s13))
    assert abs(sin*sin+cos*cos-1)<1e-11
    a=np.array([4*s13*(1-s13),s12,s23])
    return {"masses_eV":masses[::-1].tolist(),"reactor_solar_atm":a.tolist(),
            "windows_pass":((ANG[:,0]<=a)&(a<=ANG[:,1])).tolist(),
            "p1":p[:,0].tolist(),"J":J,"sin_delta":sin,"cos_delta":cos,
            "delta_deg_mod360":float(np.degrees(np.arctan2(sin,cos))%360)}


def sign_budget(matrix,frame):
    left=frame.T@matrix@frame
    gram=left.conj().T@left
    triple=np.array([gram[0,1],gram[1,2],gram[2,0]])
    signal=abs(np.prod(triple).imag)
    amp=abs(triple)
    heavy=np.linalg.norm(matrix,2)
    def gap(eps):
        de=eps*(2*heavy+eps)
        return np.prod(amp+de)-np.prod(amp)-signal
    l,h=0.,heavy
    for _ in range(60):
        mid=(h+l)/2
        if gap(mid)<0: l=mid
        else: h=mid
    return l


def accelerator(obs):
    api=runpy.run_path(str(HERE/"../1705/independent_check.py"))
    poly=read("../1705/independent_results.json")["official_coordinates_rad"]
    x=obs["reactor_solar_atm"][2]
    d=np.radians((obs["delta_deg_mod360"]+180)%360-180)
    conf={name:bool(api["interior"](poly[name],x,d)) for name in ["conf68_0","conf90_0","conf90_1"]}
    old=read("../1705/inputs.json")["t2k_2023"]
    end=old["fc2_near_zero_rad"][1]
    w=old["printed_endpoint_half_unit_rad"]
    conf["T2K2023_FC90"]=bool(old["fc90_rad"][0]+w<=d<=old["fc90_rad"][1]-w)
    conf["T2K2023_two_sigma_near_zero_endpoint"]=("in" if d<end-w else "out" if d>end+w else "rounded")
    conf["delta_rad"]=d
    return conf


def solve_case(k,dz,role,step=.20):
    d,rep,track,h,err,port_corr,common=graph(k,dz,role,step)
    raw=pair_sources(d)
    low=lower_tensors(raw,track)
    roots,poly=ratio_roots(low)
    frame=zer["unpack"](rep["charged_left_MZ"])
    sols=[]
    for x in roots:
        mat=low[0]+x*low[1]
        scale=math.sqrt(sum(freq))/np.linalg.norm(mat)
        obs=observables(mat*scale,frame)
        np.testing.assert_allclose(np.array(obs["masses_eV"])[1:]**2,freq,rtol=1e-10)
        obs["H0_H6_multipliers"]=[scale,scale*x]
        obs["cp_sign_stability_residual_eV"]=sign_budget(mat*scale,frame)
        obs["accelerator"]=accelerator(obs)
        HH=H0@np.diag(obs["H0_H6_multipliers"])
        block=np.block([[C0*0,C0],[C0.T,HH]])
        obs["max_shift_in_tree_heavy_poles_GeV"]=max(abs(np.sort(np.linalg.svd(block,compute_uv=False))-
                                                    np.array(src["neutral_poles_GeV"])))
        sols.append(obs)
    key=f"{'lower_down' if dz else 'central'}_k{k}_{role}"
    past=read(cfg["zero_results"])["cases"][key]
    zerodiff=np.max(abs(np.array(sols[0]["p1"])-np.array(past["p_MZ"])))
    # A separate broken-phase finite h/Z calculation on this very graph. It
    # replaces, and is never added to, the sequential logarithmic result.
    on=runpy.run_path(str(HERE/"../1703/independent_check.py"))
    cf=[]
    orig=d@U.T
    y_left=orig.conj()/vev
    beta_C=(y_left.conj().T@y_left)@C0
    vh=read("../1372/inputs.json")["atlas"]["mass_GeV"][0]
    for i in range(2):
        hs=H0@np.diag(np.eye(2)[i])
        cf.append(on["finite"](orig,C0,hs,fl["MZ"],vh,vev*math.sqrt(2),32)*1e9)
    alt=raw[0]+np.array(cf)
    alt_out=[]
    for x in ratio_roots(alt)[0]:
        tmp=alt[0]+x*alt[1]
        scale=math.sqrt(sum(freq))/np.linalg.norm(tmp)
        o=observables(tmp*scale,frame)
        alt_out.append({"H0_H6_multipliers":[scale,scale*x],"angles":o["reactor_solar_atm"],
                        "J":o["J"],"delta_deg_mod360":o["delta_deg_mod360"]})
    # Norm bound on dimension-four background feedback omitted at relative
    # first order in LNV. Keep a generous envelope over the measured corridor.
    maxH=max(np.linalg.norm(raw[1][0]*s["H0_H6_multipliers"][0]+
                            raw[1][1]*s["H0_H6_multipliers"][1],2) for s in sols)
    assert maxH<.3
    c5=.31e-9/vev**2
    dt=math.log(m[1]/m[0])/fl["LOOP"]
    rel=3*c5*(m[0]+1)*dt
    return {"charged_gamma_MZ_deg":rep["gamma_MZ_deg"],"third_modulus":h,"same_Higgs_calibration":common,
            "canonical_upper_texture_zero_derivative_GeV_per_log":
                abs(beta_C[0,1])/fl["LOOP"],
            "low_port_matching_relative_shift":port_corr,"boundary_relative_residual":err,
            "zero_column_displacement_from_1703":zerodiff,
            "relative_Y_LNV_feedback_first_order_envelope":rel,
            "separate_finite_EW_subsector":alt_out,
            "calibration_quartic":poly,"roots":sols}


def main():
    cases={f"{'lower_down' if dz else 'central'}_k{k}_{role}":solve_case(k,dz,role)
           for dz in (0,-1) for k in (1,2) for role in ("unit_third","shared_third")}
    finer=solve_case(1,-1,"shared_third",step=.1)
    standard=cases["lower_down_k1_shared_third"]["roots"][1]
    fine=finer["roots"][1]
    return {"cluster_reference_scales_GeV":m.tolist(),"cases":cases,
            "halfstep_max_angle_change":max(abs(np.array(fine["reactor_solar_atm"])-standard["reactor_solar_atm"])),
            "halfstep_J_change":abs(fine["J"]-standard["J"])}


def tidy(z):
    if isinstance(z,dict): return {k:tidy(v) for k,v in z.items()}
    if isinstance(z,list): return [tidy(v) for v in z]
    if isinstance(z,(np.floating,float)): return float(f"{z:.11g}")
    return z


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    val=tidy(main())
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==val
        print("1706 verified")
    else:
        path.write_text(json.dumps(val,indent=2)+"\n",encoding="utf-8")
