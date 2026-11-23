"""Continuum upper bounds and common charged/neutral examples on the old field table."""
from pathlib import Path
import argparse
import json
import math
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
load = lambda name: json.loads((HERE/name).read_bytes())
cfg = load("inputs.json")
m = runpy.run_path(str(HERE/"../1784/check.py"))
old = runpy.run_path(str(HERE/"../1715/check.py"))
Q = runpy.run_path(str(HERE/"../1712/check.py"))
F = Q["FL"]


def interval(record):
    return np.array([record[0]-2*math.hypot(*record[1:]),record[0]+2*math.hypot(*record[1:])])


def high_budget():
    """Use scalar differential inequalities prior to any charged-matrix shooting."""
    mz,soft,cut=F["MZ"],F["MS"],F["END"]
    tau,eta=np.log([soft/mz,cut/soft])/(16*math.pi**2)
    dn=math.log(soft/min(F["fields"]["neutral_poles_GeV"]))/(16*math.pi**2)
    tan=F["TAN"]
    su=math.hypot(tan,1)/tan
    sd=math.hypot(tan,1)
    mass=load(cfg["data_low"])["yukawa"]
    up,dnq,lep=([interval(row) for row in mass[key]] for key in ("up","down","charged_leptons"))
    mod={k:interval(v) for k,v in load(cfg["direct"])["direct_semileptonic_moduli"].items()}
    # Every inserted column is normalized by the same positive squared Higgs moment.
    norm0=8*math.pi**2*Q["MODEL"]["spectra"]["neutral_only_Higgs_mass_squared_GeV2"]/min(F["fields"]["neutral_poles_GeV"])**2
    min_g=np.min([np.min(np.sqrt(g),axis=0) for _,_,_,g in F["segments"](.25)],axis=0)
    assert np.all(min_g>=[.46,.63,.89])
    # Each set lists (top, charm, bottom, tau, squared neutral Frobenius norm).
    T,C,B,E,N=1.,.006,.023,.013,.052
    tr=3*(T*T+2*C*C)+9*B*B+3*E*E+N
    assert 1.5*T*T+tr < 8*.89**2 and up[2][1]<T
    assert up[1][1]*math.exp(tau*(tr+1.5*C*C))<C
    assert dnq[2][1]*math.exp(tau*(tr+1.5*B*B))<B
    assert lep[2][1]*math.exp(tau*(tr+1.5*E*E))<E
    assert norm0*math.exp(2*dn*(tr+1.5*N))<N
    T1,C1,B1,E1,N1=1.15,.008,.15,.074,.078
    st=3*(T1*T1+2*C1*C1)+N1
    sb=9*B1*B1+3*E1*E1
    gt=16*.89**2/3+3*.63**2+13*.46**2/15
    assert T*su*math.exp(eta*(st+3*T1*T1+B1*B1-gt))<T1
    assert C*su*math.exp(eta*(st+3*C1*C1+B1*B1-gt))<C1
    assert B*sd*math.exp(eta*(sb+3*B1*B1+T1*T1))<B1
    assert E*sd*math.exp(eta*(sb+3*E1*E1+N1))<E1
    assert N*su**2*math.exp(2*eta*(st+3*N1+E1*E1-3*.63**2-3*.46**2/5))<N1
    # Singular-ratio drift; universal terms and the common Higgs projection cancel.
    ratio_log=1.5*tau*(T*T+B*B)+eta*(3*B1*B1+T1*T1)
    up_log=1.5*tau*(T*T+B*B)+eta*(3*T1*T1+B1*B1)
    allrat=max(max(dnq[i][1]/dnq[i+1][0] for i in (0,1))*math.exp(ratio_log),
               max(up[i][1]/up[i+1][0] for i in (0,1))*math.exp(up_log))
    gap=(1+allrat**2)/(1-allrat**2)
    cw0=math.hypot(mod["Vub"][1],mod["Vcb"][1])
    heavy=cw0*math.exp(gap*(1.5*tau*(T*T+C*C+2*B*B)+eta*(T1*T1+C1*C1+2*B1*B1)))
    sr,SR=dnq[1][0]/dnq[2][1]*math.exp(-ratio_log),dnq[1][1]/dnq[2][0]*math.exp(ratio_log)
    light_log=1.5*tau*(B**2*SR**2+C*C+T*T*heavy**2)+eta*(3*B1**2*SR**2+C1*C1+T1*T1*heavy**2)
    ds,DS=dnq[0][0]/dnq[1][1]*math.exp(-light_log),dnq[0][1]/dnq[1][0]*math.exp(light_log)
    def mixing_drift(t,c,b):
        return gap*(2*(t*t*heavy**2+c*c)+2*b*b*(heavy**2+SR**2))
    drift=1.5*tau*mixing_drift(T,C,B)+eta*mixing_drift(T1,C1,B1)
    cd=mod["Vus"][0]*math.sqrt(1-mod["Vcb"][1]**2/(1-mod["Vub"][1]**2))-mod["Vcb"][1]*mod["Vub"][1]/math.sqrt(1-mod["Vub"][1]**2)-drift
    # Any allowed combination of columns within this slightly enlarged Cabibbo window.
    cab=mod["Vus"][1]+drift
    assert cd>0 and heavy<.06
    # Tighten the actual spectator bound using the two normalized Dirac pencils.
    model=Q["MODEL"]
    lo=(mod["Vus"][0]-drift)
    hi=cab
    def shape(u):
        one=math.sqrt(model["a"]@model["a"])*u*model["b"]
        two=model["s"]**2*(model["x"]/model["rho"]**2*model["a"]+model["g"]/model["m"]*one)
        return (model["phases"][:,None]*np.column_stack((one,two)))@model["state"]["U"].conj()
    z0=shape(0)
    z1=shape(1)-z0
    mid=shape((lo+hi)/2)
    poles=np.array(model["state"]["m"])
    radius=(hi-lo)/2
    gain=8*math.pi**2*model["spectra"]["neutral_only_Higgs_mass_squared_GeV2"]
    n0=gain*((np.linalg.norm(mid)+radius*np.linalg.norm(z1))/(np.linalg.norm(mid*poles)-radius*np.linalg.norm(z1*poles)))**2
    # The lower column lives at the first pole for the positive Higgs moment.
    # Its profile at the upper pole is the common texture used above.
    max_g=np.max([np.max(np.sqrt(gs),axis=0) for _,_,_,gs in F["segments"](.25)],axis=0)
    walk=(tr+1.5*(N+E**2)+9*max_g[0]**2/20+9*max_g[1]**2/4)*math.log(poles[1]/poles[0])/(16*math.pi**2)
    n0*=math.exp(4*walk)
    N=n0*math.exp(2*dn*(tr+1.5*N))
    N1=N*su**2*math.exp(2*eta*(st+3*N1+E1*E1-3*.63**2-3*.46**2/5))
    p=tau*1.5*E*E+dn*.5*N+eta*(E1*E1+N1)
    e_ratio=1.5*tau*E*E+1.5*dn*N+eta*(3*E1*E1+N1)
    em=lep[0][1]/lep[1][0]*math.exp(e_ratio)
    mt=lep[1][1]/lep[2][0]*math.exp(e_ratio)
    mtmin=lep[1][0]/lep[2][1]*math.exp(-e_ratio)
    charged_turn=(1+mt**2)/(1-mt**2)*(1.5*dn*N+eta*N1)/2
    # Uniform nonnegative trace/determinant. H22 <= |Vcb|²+s², H33 <=1.
    h22=heavy**2+SR**2
    L=1+h22
    A2=SR**2*(DS**2+cab*cab)+heavy**2
    # Simpler sharp upper A² follows from |Vub|<=heavy; large w still bounded well.
    K=A2+(DS*SR/cd)**2+9*h22
    alpha=((em*mt**2)*K*math.sqrt(L)/(ds*sr**2))**2
    assert alpha<4/27
    low_w=1/3
    for test in np.linspace(1/3,1.,1001):
        if test*(1-test)**2>alpha: low_w=test
    # For ds/(|Vcd|²+ds²|Vcs|²), use |Vcs|²>=1-|Vus|²-|Vts|².
    denominator=cd**2+ds**2*(1-cab**2-heavy**2)
    bymu=math.sqrt(em*DS/denominator*math.sqrt(L*(1+mt**2+(em*mt)**2)/low_w))
    tau_leg=bymu*mt
    # Projective angular motion of the normalized inverse image, from its ray velocity.
    move=p/2
    betamu=bymu+charged_turn+move
    betatau=tau_leg+charged_turn+move
    return locals()

def continuum():
    b=high_budget()
    g=lambda name:b[name]
    r,R=load(cfg["osc"])["assignments"]["NO"]["middle_over_heaviest_mass_required"]
    data=load(cfg["osc"])["experimental_rectangles"]
    x,X=data["isolated_weight_roots"][0]
    s,S=data["solar_share"]
    t,T=m["read"](m["CFG"]["atmospheric"])["joint_windows"][2]
    mt,em=g("mt"),g("em")
    # Light-up to heavy-down mixing. The common top piece contracts to Vub
    # by unitarity; its s/d mass-dependent remainder is explicitly bounded.
    ratio=g("SR"); heavy=g("heavy"); gap=g("gap")
    def ub_slope(top,bottom):
        extra=2*ratio**2/(1-ratio**2)*top**2*heavy
        up_piece=gap*bottom**2*(heavy+ratio**2)*(1+heavy)
        return top**2,extra+up_piece
    a1,b1=ub_slope(g("T"),g("B"))
    a2,b2=ub_slope(g("T1"),g("B1"))
    ub=(g("mod")["Vub"][1]+1.5*g("tau")*b1+g("eta")*b2)*math.exp(1.5*g("tau")*a1+g("eta")*a2)
    # Exact F²=minor/A². All ratios here use down eigenvalues normalized to y_b.
    fracF=g("denominator")/(g("DS")**2+g("cab")**2+(ub/g("sr"))**2)
    assert 0<fracF<1
    height=math.sqrt(1-fracF)
    z=math.sqrt(g("h22")/(1-heavy**2))
    # Units y_tau=1 for the charged texture entries.
    l1=(mt*em)**2
    l2=g("mtmin")**2
    bleg=mt*g("bymu")
    rem=(1-g("low_w"))*(1+mt**2+l1)
    # ||lambda A||²<= rem, and a stronger ratio from its quark-normalized trace.
    # A²/B²=A² minor/(d²s²), better obtained from the two measured first-row moduli.
    # |Vcd|² dominates minor; retain its own upper bound.
    cdmax=g("cab")+ub*heavy+g("drift")+.0001
    AoverB2=(g("DS")**2+g("cab")**2+(ub/g("sr"))**2)*(cdmax**2+g("DS")**2+g("DS")**2*g("SR")**2)/g("ds")**2
    aleg2=min(rem,AoverB2*bleg**2)
    ce=math.sqrt(rem)
    # Im P_e is the antisymmetric cross-product matrix from the spectral polynomial.
    ke=z*math.hypot(ce*l1,bleg*(aleg2+l1+height))/((l2-l1)*(1-l1))
    kt=z*math.hypot(ce,bleg*(1+aleg2+height))/((1-l1)*(1-mt**2))
    def ray_norm(k):
        assert 0<=k<.5
        return math.sqrt(2)*k/math.sqrt(1+math.sqrt(1-4*k*k))
    # Ray speeds of signed Takagi columns are bounded B*span(P)/2.
    # A change in two charged rows and two central columns is estimated in their pair norm.
    p=g("p")
    mratio=R*math.exp(2*p)
    turn=.5*(1+mratio)/(1-mratio)*p
    extra=g("charged_turn")+math.sqrt(2)*turn
    ei,ti=ray_norm(ke)+extra,ray_norm(kt)+extra
    mu,ta=g("betamu"),g("betatau")
    be=math.sqrt(1-mu*mu-ta*ta)
    tt=R*math.sqrt(S*(1-S)*T)+math.sqrt(X*(1-t))*(1+R*S)
    tail=(mu*math.sqrt(T*(1-x))+ta*math.sqrt((1-t)*(1-x)))**2
    hmin=(be**2*(r*s*(1-X)-X)-2*ta*tt-tail)/(2*mu)
    aa=[r*math.sqrt(s*(1-s)*(1-T)),R*math.sqrt(S*(1-S)*(1-t))]
    bb=[math.sqrt(x*t)*(1+r*s),math.sqrt(X*T)*(1+R*S)]
    assert hmin>0
    ct=max((a*a+b0*b0-hmin*hmin/(1-x))/(2*a*b0) for a in aa for b0 in bb)
    def sine_cap(cos_upper):
        # max of |tau2 tau3|, retaining t(1-t)<=1/4 and the controlled interference.
        tautau=math.sqrt((1-x)*((1-s)*.25+S*(1-t)**2*X+
              2*(1-t)*math.sqrt(S*(1-s)*.25*X)*max(0,cos_upper)))
        ee=math.sqrt(X*S*(1-x))
        cross_e=math.sqrt(S*(1-x)+X)*ei
        cross_t=ti  # the norm of the two complex tau entries is at most one
        J=tautau*cross_e+ee*cross_t
        pref=(1-X)*math.sqrt(x*s*(1-s)*min(t*(1-t),T*(1-T)))
        return J,J/pref
    J,st=sine_cap(ct)
    assert st<1 and math.sqrt(1-st**2)>ct
    J,st=sine_cap(-math.sqrt(1-st**2))
    m2low,m2high=np.sqrt(data["solar_frequency_eV2"])
    m3low,m3high=np.sqrt(load(cfg["osc"])["assignments"]["NO"]["absolute_atmospheric_frequency_eV2"])
    tilt=math.hypot(mu,ta)
    a_min=math.sqrt(m2low*(1-X)*s)
    b_min,b_max=math.sqrt(m3low*x),math.sqrt(m3high*X)
    distance=tilt*math.sqrt(m2high+m3high)
    mbb=[1000*(a_min*a_min-b_max*b_max),1000*distance*(2*b_max+distance)]
    defect=4*math.asin(math.sqrt((distance**2-(a_min-b_max)**2)/(4*a_min*b_min)))
    return {"mixed_quark_leak":heavy,"max_Vub_high":ub,"min_Vcd_high":g("cd"),
            "high_yd_ys":[g("ds"),g("DS")],"Ynu2_original_bound":g("n0"),"positive_span_integral":p,
            "charged_axis_turn":g("charged_turn"),"Takagi_axis_turn":turn,
            "large_branch_w_lower":g("low_w"),"b_to_mumu":g("bymu"),
            "zero_axis_mu_tau_caps":[mu,ta],"zero_axis_total_tilt_cap":math.hypot(mu,ta),
            "quark_F_fraction_lower":fracF,"charged_A_squared_upper":aleg2,
            "high_e_tau_imaginary_projector_norms":[ke,kt],
            "CP_cos_upper_before_reality":ct,"partial_imaginary_norms":[ei,ti],
            "J_abs_upper":J,"delta_distance_pi_upper_deg":math.degrees(math.asin(st)),
            "global_light_exchange_meV":mbb,"global_relative_Majorana_distance_pi_deg":math.degrees(defect),
            "between_spectator_poles_log_norm_budget":g("walk")}

def endpoint(k, step, shifts=None):
    generator=runpy.run_path(str(HERE/"../1712/check.py"))
    flow=generator["FL"]
    solve=flow["solve"]
    rep={}
    def capture(*args,**kwargs):
        rep.update(solve(*args,**kwargs))
        return rep
    flow["solve"]=capture
    source=load(cfg["source_points"])["points"]["interior_all_one_width"]
    source.update(shifts or {})
    trial=generator["run_point"](**source,k=k,step=step)
    W, transported=old["congruence"](flow,rep,step)
    wanted=np.array(cfg["joint_osc_shape"])
    all_points=[]
    for i in np.linspace(-math.pi,math.pi,7)[:-1]:
        for j in np.linspace(-math.pi,math.pi,11)[:-1]:
            u=np.array([i+.13,j+.29,.1695])
            for _ in range(30):
                if u[2]<=0 or u[2]>=1: break
                beam=m["observable"](W,u)[0]
                if np.linalg.norm(beam-wanted)<1e-9: break
                jac=np.column_stack([(m["observable"](W,u+v)[0]-m["observable"](W,u-v)[0])/2e-5
                                     for v in np.eye(3)*1e-5])
                if np.linalg.cond(jac)>1e9: break
                move=np.linalg.solve(jac,beam-wanted)
                u-=move/max(1,np.linalg.norm(move)/.5)
            if 0<u[2]<1:
                y,p,(_,_,_,_,sv)=m["observable"](W,u)
                if max(abs(y-wanted))<1e-8 and .0803<=p[0]<=.0899:
                    if not any(abs(p[1]-v["delta_deg"])<.0001 for v in all_points):
                        all_points.append({"reactor":p[0],"delta_deg":p[1],"J":p[2],
                           "high_spherical":u.tolist(),"d31_eV2":7.255e-5*(sv[0]/sv[1])**2})
    q=trial["quark_neutrino"]["quarks"]
    exp=load(cfg["quark_holdouts"])
    gold=interval(exp["S_d"])
    qpass=exp["time_integrated_gamma_deg"]["95"][0]<=q["gamma_deg"]<=exp["time_integrated_gamma_deg"]["95"][1]
    assert qpass and gold[0]<=q["golden_S"]<=gold[1]
    return {"gamma_deg":q["gamma_deg"],"golden_S":q["golden_S"],"Bs_phi_rad":q["phi_s_rad"],
            "passes_time_integrated_gamma_68":exp["time_integrated_gamma_deg"]["68"][0]<=q["gamma_deg"]<=exp["time_integrated_gamma_deg"]["68"][1],
            "passes_time_integrated_gamma_95_and_golden_two_widths":True,
            "generated_yd_MZ":trial["y_d_MZ"],"selected_two_source_solutions":all_points}

def main():
    cap=continuum()
    assert cap["zero_axis_total_tilt_cap"] < .075
    assert cap["delta_distance_pi_upper_deg"] < 28
    halfstep=0.
    examples={}
    for k in (1,2):
        fine=endpoint(k,.10,{"cb":.5})
        coarse=endpoint(k,.20,{"cb":.5})
        assert fine["passes_time_integrated_gamma_68"]
        assert abs(fine["golden_S"]-.717)<math.hypot(.013,.008)
        points=fine["selected_two_source_solutions"]
        assert points and len(points)==len(coarse["selected_two_source_solutions"])
        for i,row in enumerate(points):
            assert (-1 if k==1 else 1)*row["J"]>0
            halfstep=max(halfstep,abs(row["delta_deg"]-coarse["selected_two_source_solutions"][i]["delta_deg"]))
        examples[str(k)]=fine
    assert halfstep<2e-7
    atlas=m["read"](m["CFG"]["phase_literature"])["t2k_updated"]["fc1_center_minus_plus_rad"]
    gap=atlas[0]-atlas[1]+math.pi-math.radians(cap["delta_distance_pi_upper_deg"])
    assert gap>0 and gap<.01
    return {"continuum":cap,"two_discrete_branches_common_input":examples,"max_halfstep_delta_change_deg_below":2e-7,
            "gap_to_printed_T2K_FC_lower_rad":gap,"worst_two_rounded_inputs_rad":.01}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    out=m["clean"](main())
    dest=HERE/"results.json"
    if ap.parse_args().check:
        assert out==json.loads(dest.read_bytes())
        print("1786 verified: continuum charged inputs, two opposite signs, oscillation and light-LNV bounds")
    else:
        dest.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf-8")





