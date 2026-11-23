"""Test the lightest-state electron fraction over the two-dimensional flavour-phase torus."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
cfg=json.loads((HERE/"inputs.json").read_bytes())


def read(path):
    return json.loads((HERE/path).read_bytes())


def unpack(x):
    return np.array(x["real"])+1j*np.array(x["imag"])


def charged_background(step):
    generator=runpy.run_path(str(HERE/cfg["charged_generator"]))
    flow=generator["FL"]
    original=flow["solve"]
    captured={}
    def record(*args,**kwargs):
        rep=original(*args,**kwargs)
        captured["rep"]=rep
        return rep
    flow["solve"]=record
    row=generator["run_point"](**read(cfg["charged_sources"])["points"][cfg["point"]],step=step)
    ref=read(cfg["frozen_result"])["cases"][cfg["point"]]
    assert abs(row["y_d_MZ"]/ref["y_d_MZ"]-1)<2e-9
    rep=captured["rep"]
    assert abs(rep["gamma_MZ_deg"]-ref["quark_neutrino"]["quarks"]["gamma_deg"])<2e-5
    return flow,rep


def congruence(flow,rep,step):
    y=np.zeros((4,3,4),complex)
    E0=unpack(rep["charged_left_MZ"])
    EH=unpack(rep["charged_left_10PeV"])
    y[0,:,:3]=np.diag(flow["central"]("up"))
    y[1,:,:3]=flow["low_ckm"](rep["gamma_MZ_deg"],flow["angles"]("direct"))@np.diag(flow["central"]("down"))
    y[2,:,:3]=E0@np.diag(flow["central"]("charged_leptons"))
    # Kappa* evolves upwards as B kappa*(low) B^T. The desired inverse dagger solves S'=-P S.
    S=np.eye(3,dtype=complex)
    B=np.eye(3,dtype=complex)
    total_log=0.
    def rhs(z,sv,bv,g,susy):
        E,N=(f@f.conj().T for f in z[2:])
        P=(E+N) if susy else (-1.5*E+.5*N)
        return flow["beta"](z,g,susy), -P@sv,P@bv,np.linalg.norm(P,2)
    for start,susy,dt,gg in flow["segments"](step):
        if start==flow["MS"]:
            tan=flow["TAN"]
            y/=np.array([tan,1,1,tan])[:,None,None]/math.hypot(tan,1)
        for j,m in enumerate(flow["fields"]["neutral_poles_GeV"]):
            if start==m:
                y[3,:,j]=flow["physics"]["Y0"][:,j]
        for j in range(0,len(gg)-1,2):
            k1,f1,h1,p1=rhs(y,S,B,gg[j],susy)
            k2,f2,h2,p2=rhs(y+dt*k1/2,S+dt*f1/2,B+dt*h1/2,gg[j+1],susy)
            k3,f3,h3,p3=rhs(y+dt*k2/2,S+dt*f2/2,B+dt*h2/2,gg[j+1],susy)
            k4,f4,h4,p4=rhs(y+dt*k3,S+dt*f3,B+dt*h3,gg[j+2],susy)
            y+=dt*(k1+2*k2+2*k3+k4)/6
            S+=dt*(f1+2*f2+2*f3+f4)/6
            B+=dt*(h1+2*h2+2*h3+h4)/6
            total_log+=dt*(p1+2*p2+2*p3+p4)/6
    check=np.linalg.norm(S.conj().T@B-np.eye(3))
    assert check<1e-11
    active=np.linalg.svd(y[2],full_matrices=False)[0][:,::-1]
    np.testing.assert_allclose(abs(EH.conj().T@active),np.eye(3),atol=3e-8)
    np.testing.assert_allclose(y[3],unpack(rep["neutral_columns_10PeV"]),atol=3e-10)
    down,leptons,gamma,_,up=flow["extraction"](y)
    np.testing.assert_allclose(up,rep["up_texture_diag_10PeV"],rtol=1e-9,atol=1e-14)
    np.testing.assert_allclose(down[:3],rep["Y_down_at_10PeV"],rtol=1e-9,atol=1e-14)
    np.testing.assert_allclose(leptons,rep["Y_charged_at_10PeV"],rtol=1e-9,atol=1e-14)
    np.testing.assert_allclose(down[3:],rep["CKM_angles_at_10PeV_rad"],rtol=1e-9,atol=1e-13)
    assert abs(gamma-rep["gamma_10PeV_deg"])<2e-8
    return S@E0, {"integral_norm_P":total_log,"inverse_dagger_closure":check,
                   "transport_distance_from_identity":np.linalg.norm(S-np.eye(3),2)}


def phase_bound(W,up,nt):
    X=np.linalg.solve(W,np.diag(np.sqrt(up)))
    norm2=np.sum(abs(X)**2,axis=0)
    ee=abs(X[0])**2
    projected=np.linalg.norm(X[1:,:],axis=0)
    sigma_floor=2*max(ee)-sum(ee)
    q_upper=np.sum(np.sqrt(norm2)*projected)
    assert sigma_floor>q_upper
    worst_tilt=q_upper/sigma_floor
    elementary_bound=1-worst_tilt**2
    # Keep the heaviest inverse contribution fixed by common phase. Cover the two remaining arguments.
    free=np.argsort(norm2)[:2]
    fixed=int(np.argmax(norm2))
    centers=(np.arange(nt)+.5)*2*math.pi/nt
    e1,e2=np.meshgrid(np.exp(1j*centers),np.exp(1j*centers),indexing="ij")
    pieces=np.array([np.outer(x,x) for x in X.T])
    matrices=pieces[fixed]+e1[:,:,None,None]*pieces[free[0]]+e2[:,:,None,None]*pieces[free[1]]
    U,s,_=np.linalg.svd(matrices)
    weights=abs(U[:,:,0,0])**2
    # Perturb the Hermitian squared inverse. The largest-eigenvalue displacement is <= eta.
    delta=2*math.sin(math.pi/(2*nt))*sum(norm2[free])
    eta=(2*s[:,:,0]+delta)*delta
    gap=s[:,:,0]**2-s[:,:,1]**2
    assert np.min(gap-eta)>0
    tilt=eta/(gap-eta)
    cos_center=np.sqrt(np.clip(weights,0.,1.))
    lower=np.maximum(0,cos_center*np.sqrt(1-tilt**2)-np.sqrt(np.maximum(0,1-weights))*tilt)**2
    # Verify the conjugations and the order of singular vectors directly at an off-axis phase.
    test=np.array([1.,np.exp(.39j),np.exp(-1.47j)])
    mass=W.T@np.diag(test/up)@W
    inverse=X@np.diag(1/test)@X.T
    residual=np.linalg.norm(mass@inverse-np.eye(3),2)
    assert residual<16*np.finfo(float).eps*np.linalg.cond(mass)
    left,sv,_=np.linalg.svd(mass)
    v,si,_=np.linalg.svd(inverse)
    np.testing.assert_allclose(abs(left[:,-1])**2,abs(v[:,0])**2,atol=5e-13)
    return {"inverse_column_norms_squared":norm2,"electron_lower_bound_on_largest_inverse_singular_value":sigma_floor,
            "transverse_upper_bound":q_upper,"analytic_electron_probability_floor":elementary_bound,
            "actual_center_probability_range":[weights.min(),weights.max()],
            "torus_certified_probability_floor":lower.min()-1e-8,
            "minimum_squared_gap_over_eta":(gap/eta).min(),"cells":nt**2}


def calculate(step):
    flow,rep=charged_background(step)
    D=unpack(rep["down_texture_10PeV"])
    up=np.array(rep["up_texture_diag_10PeV"])
    # D has row 10 and column 5bar. The transported charged-frame columns carry the open 5bar legs.
    charged,transport=congruence(flow,rep,step)
    high=phase_bound(D@unpack(rep["charged_left_10PeV"]),up,cfg["phase_axis_cells"])
    low=phase_bound(D@charged,up,cfg["phase_axis_cells"])
    old_gate=runpy.run_path(str(HERE/"../1703/check.py"))["growth_gate"]
    leptons=np.array(rep["Y_charged_at_10PeV"])
    ray,charged_ray=old_gate(unpack(rep["neutral_columns_10PeV"]),leptons,
            max(leptons[:-1]/leptons[1:]),flow["MS"],flow["MZ"],
            flow["fields"]["neutral_poles_GeV"],flow["TAN"])
    # R1703 bounds the charged-Yukawa integral in 2*ray and the neutrino integral in
    # 2*charged_ray (including a gap factor >=1 and a SM coefficient >=1/2).
    p_integral_cap=2*(ray+charged_ray)
    eps=math.expm1(2*p_integral_cap)
    transverse=math.sqrt(1-high["analytic_electron_probability_floor"])
    robust_floor=1-((transverse+math.sin(charged_ray)+eps)/(1-eps))**2
    assert transport["integral_norm_P"]<p_integral_cap and robust_floor>.85
    transport.update({"R1703_integral_norm_P_cap":p_integral_cap,"charged_electron_angle_cap":charged_ray,
                      "analytic_low_probability_floor_without_RK_integration":robust_floor})
    angles=read(cfg["oscillation"])["joint_windows"]
    zlo=(1-math.sqrt(1-angles[0][0]))/2
    # The small reactor root is the convention used for the measured solar/atmospheric mixing angles.
    cap=(1-zlo)*(1-angles[1][0])
    s0=math.sqrt(1-low["analytic_electron_probability_floor"])
    st=math.sqrt(1-cap)
    # For an extra norm <= d of the inverse mass matrix in the measured charged frame,
    # Q*inverse <= (s0+d)||inverse|| and new sigma_max >= (1-d)||inverse||.
    extra=(st-s0)/(1+st)
    assert low["torus_certified_probability_floor"]>cap
    return {"source_quark_gamma_deg":rep["gamma_MZ_deg"],"transport":transport,"upper_10PeV":high,
            "low_MZ":low,"largest_allowed_electron_weight_over_all_mass_labels":cap,
            "necessary_extra_relative_inverse_matrix_norm":extra,
            "necessary_extra_nonuniversal_congruence_norm":1-1/math.sqrt(1+extra)}


def tidy(v):
    if isinstance(v,dict):return {k:tidy(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):return [tidy(x) for x in v]
    if isinstance(v,(np.floating,float)):return float(f"{v:.11g}")
    return v


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    fine=calculate(cfg["log_steps"][1])
    coarse=calculate(cfg["log_steps"][0])
    delta=abs(fine["low_MZ"]["torus_certified_probability_floor"]-coarse["low_MZ"]["torus_certified_probability_floor"])
    assert delta<2e-8
    result=tidy({**fine,"halfstep_probability_floor_change":delta})
    target=HERE/"results.json"
    if args.check:
        assert result==json.loads(target.read_bytes())
        print("1715 verified, all full-rank phases and labelled flavour transport")
    else:
        target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")

