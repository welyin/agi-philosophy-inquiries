"""Build an independent 5-flavour family and audit tree calibration and observables."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(n, filename="inputs.json"):
    return json.loads((HERE.parent/str(n)/filename).read_text(encoding="utf-8-sig"))


def rims(row):
    return [row[0]-2*sum(row[1:]),row[0]+2*sum(row[1:])]


def check(x,interval):
    assert interval[0] <= x <= interval[1], (x,interval)


def phase(z):
    return math.degrees(np.angle(z))%180


def ckm(a,b,c,delta):
    # R23 U13 R12 by explicit matrices.
    A,B,C=np.sqrt(1-np.array([a,b,c])**2)
    r12=[[A,a,0],[-a,A,0],[0,0,1]]
    u13=[[C,0,c*np.exp(-1j*delta)],[0,1,0],[-c*np.exp(1j*delta),0,C]]
    r23=[[1,0,0],[0,B,b],[0,-b,B]]
    return np.array(r23)@u13@r12


def ordinary_rows(fixed_rows):
    # Embed the coordinate light basis in the orthogonal complement and correct its metric.
    projector=np.eye(5)-fixed_rows.conj().T@fixed_rows
    w,v=np.linalg.eigh(projector[:3,:3])
    return ((v*w**-.5)@v.conj().T)@projector[:3,:]


def make_mass(M,m5,F,V):
    pol=read(1633,"results.json")
    m1,ms,mb=read(1635)["witness"]["light_masses_GeV"]
    rb=np.sqrt(pol["right_doublet_probability"]); rs=-F/rb
    sb=math.sqrt(2*(pol["left_difference"]/pol["k"]-(mb*rb/M)**2))
    def build(r5):
        r4=math.sqrt(1-rb**2-abs(rs)**2-r5*r5)
        m4=math.sqrt((M**2-(ms*abs(rs))**2-(mb*rb)**2-(m5*r5)**2)/r4**2)
        mass=np.array([m1,ms,mb,m4,m5]); r=np.array([0,rs,rb,r4,r5])
        ell=mass*r/M
        eq=np.array([ell[3:],mass[3:]**2*ell[3:]])
        ss=np.r_[[0,0,sb],np.linalg.solve(eq,-sb*ell[2]*np.array([1,mb**2]))]
        return mass,r,ell,ss
    lo,hi=1e-10,.01
    assert np.linalg.norm(build(lo)[3])>1>np.linalg.norm(build(hi)[3])
    for _ in range(64):
        mid=(lo+hi)/2
        if np.linalg.norm(build(mid)[3])>1:
            lo=mid
        else:
            hi=mid
    m,r,l,s=build((lo+hi)/2)
    rf=np.array([m*s/np.linalg.norm(m*s),r]).conjugate()
    lf=np.array([s,l]).conjugate()
    left=np.vstack([V@ordinary_rows(lf),lf])
    right=np.vstack([ordinary_rows(rf),rf])
    matrix=left@np.diag(m)@right.conj().T
    forbidden=max(np.max(abs(matrix[:4,4])),np.max(abs(matrix[3,:3])),abs(matrix[4,3]))
    assert forbidden<2e-9 and abs(matrix[4,4]-M)<2e-9
    u,sm,vh=np.linalg.svd(matrix)
    order=np.argsort(sm)
    u=u[:,order];v=vh.conj().T[:,order]
    assert np.max(abs(sm[order]-m))<1e-8
    return u,v,sm[order],matrix


def audit(M,m5,main_calibration=False):
    s23=.03724; ratio=.083; delta=math.radians(68.5)
    if main_calibration:
        s23=.0372;delta=math.radians(68.65)
    V=ckm(.223,s23,ratio*s23,delta)
    F0=.000680-.000130j
    if main_calibration:
        F0=.000680-.000125j
    u,v,m,mat=make_mass(M,m5,F0,V)
    K=u[:3,:3]
    blend=lambda vec:vec[1].conjugate()*vec[2]
    lt,lc=blend(K[2]),blend(K[1]); F=-blend(v[4]); fL=-blend(u[3])-2*blend(u[4])
    ms,mb=m[1:3]
    assert abs(fL/F-2*ms*mb/M**2)<1e-9
    data=read(1445); old=read(1381)
    us,ud,cb,ub=abs(K[0,1]),abs(K[0,0]),abs(K[1,2]),abs(K[0,2])
    ff_K=data["K_form"][0]
    check(us*ff_K,rims(data["K_product"]));check(cb,rims(data["BelleII_Vcb"]))
    check(ub/cb,rims(data["ub_over_cb_holdout"]))
    # Rephasing-invariant angles for distinct tree D tags; modulo-pi B strong-phase ambiguity.
    gd=phase(-K[0,0]*K[0,2].conjugate()/(K[1,0]*K[1,2].conjugate()))
    gs=phase(K[0,1]*K[0,2].conjugate()/(K[1,1]*K[1,2].conjugate()))
    prod_ratio=K[0,2]*K[1,1].conjugate()/(K[1,2]*K[0,1].conjugate())
    dCF=K[1,1].conjugate()*K[0,0]
    antiD_DCS=K[1,0]*K[0,1].conjugate()
    dKK=K[1,1].conjugate()*K[0,1]; dpi=K[1,0].conjugate()*K[0,0]
    tags={"CF_DCS_Kpi":phase((prod_ratio*antiD_DCS/dCF).conjugate()),
          "tree_KK":phase((prod_ratio*dKK.conjugate()/dKK).conjugate()),
          "tree_pipi":phase((prod_ratio*dpi.conjugate()/dpi).conjugate())}
    assert abs(tags["CF_DCS_Kpi"]-gd)<1e-10 and abs(tags["tree_KK"]-gs)<1e-10
    for val in tags.values():
        check(val,[68.7-5.1,68.7+5.2]);check(val,rims(data["Run3_gamma_deg"]))
    # Absolute calibrations, retaining an explicit shared f_Bs and an inferred RGI bag.
    c=read(1404);h=read(1407);r=read(1464);mix=read(1468);cfg=read(1406)
    z=r["HPQCD_tableI_eq15"]; cc=z["Xt_ref"]/(math.pi*z["sine_squared_ref"]*z["alpha_inverse_ref"])
    L=-cc*lt; P=(.605 if main_calibration else .60)*L+fL-F; N=L+fL+F
    gc=c["const_GeV"];gf,mw,mes=gc["GF"],gc["MW"],gc["M_Bs"]; hb=h["hbar_GeV_ps"]
    norm=gf**2*mw**2/(4*math.pi**2);factor=2*mes/(3*hb); sc=c["short_distance"]
    mt=sc["mtMS"][0]; xt=sc["xt"][0]; av=mix["alpha_s_ref"]; mz=h["running"]["MZ"]
    at=av/(1+av*(23/3)*math.log(mt/mz)/(2*math.pi))
    am=at/(1+at*7*math.log(M/mt)/(2*math.pi)); qcd=(at/av)**(6/23)*(am/at)**(2/7)
    h1=-(xt-7)/(4*(xt-1))-(xt*xt-2*xt+4)*math.log(xt)/(2*(xt-1)**2)
    light=(1-ms**2/M**2)*(1-mb**2/M**2)
    vev=read(1635)["witness"]["Higgs_vev_GeV"]
    hard=(-2*M**2*light*F/vev**2)**2/(64*math.pi**2*M**2)
    a=norm*lt**2*rims(sc["etaB"])[0]*rims(sc["S0_published"])[0]
    if main_calibration:
        a=norm*lt**2*.5510*2.299
    a+=(.547 if main_calibration else .54)*(gf*(F**2+fL**2)/math.sqrt(2)+qcd*hard)
    dC=F/lt*xt*(math.log(M/mw)+h1)*light
    d=-norm*lt**2*dC*(1.18*(.0534-2*.0032)+1.42*(.0493-2*.0037))
    if main_calibration:
        d=-norm*lt**2*dC*(1.18*.0477+1.42*.0428)
    # Solve |fhat^2*a+d|=Omega/factor for the physical root without using any rare rate.
    roots=np.roots([abs(a)**2,2*(a*d.conjugate()).real,abs(d)**2-(mix["omega_ps_inverse"][0]/factor)**2])
    pos=[x.real for x in roots if abs(x.imag)<1e-12 and x.real>0]
    assert len(pos)==1
    fhat=math.sqrt(pos[0]);fb=.229 if main_calibration else .230;bag=(fhat/fb)**2
    check(fhat,rims(h["HPQCD_inputs"]["fBs_sqrtB_GeV"]));check(fb,rims(r["fBs_GeV"]));check(bag,rims(c["bag"]["s"]))
    W=factor*(fhat**2*a+d); measured_phase=np.angle(W/lc**2)+(.013 if main_calibration else .015)
    check(abs(W),rims(mix["omega_ps_inverse"]))
    for lab in ("LHCb_KK_2015_18","CMS_KK_13TeV_final"):
        check(measured_phase,rims(mix["phase_inputs_rad"][lab]))
    belle=cfg["BelleII"];br_n=1.08*abs(N/(z["lambda_ref"]*cc))**2
    check(br_n,rims([belle["mu"],belle["width"]]))
    mu=gc["m_mu"];beta=math.sqrt(1-4*(mu/mes)**2); th=c["tau_ps"]["s_heavy"][0]
    hm=gf**2*mu**2*mes*beta*fb**2*th/(16*math.pi*hb)
    q=(np.angle(P/lt)-np.angle(W/lt**2)/2)
    resp=read(1464,"results.json")
    acc={}
    cms=cfg["CMS"];run=c["LHCb_run3_2024_25"]["s_in_1e9"]
    bounds={"CMS":[cms["branching_in_1e-9"]-2*sum(cms["lower_widths"]),cms["branching_in_1e-9"]+2*sum(cms["upper_widths"])],
            "Run3":[run[0]-2*(run[1]+run[3]),run[0]+2*(run[2]+run[4])]}
    for lab,w in (("CMS",resp["calculated_light_response"]),("Run3",resp["run3_central_light_response"])):
        axial=1e9*hm*abs(P)**2*(1-(1-w)*math.sin(q)**2)
        scalar=1e9*hm*((mes/125.1)**2*abs(F))**2
        check(axial,bounds[lab]); check(axial+scalar,bounds[lab]); acc[lab]=[axial,axial+scalar]
    yukawa=mat.copy();yukawa[3,3]=yukawa[4,4]=0
    ytop=np.linalg.svd(yukawa,compute_uv=False)[0]*math.sqrt(2)/vev
    assert ytop<math.sqrt(4*math.pi)
    return {"M":M,"m5":m5,"weak_basis_mass_re_im":np.stack((mat.real,mat.imag),axis=-1).tolist(),
            "singular_masses_GeV":m.tolist(),"Ymax":ytop,"tree_us_ud_cb_ub":[us,ud,cb,ub],"rho":ub/cb,
            "Kell3_product":us*ff_K,"gamma_d":gd,"gamma_s":gs,"actual_tag_angles":tags,
            "fBs_GeV":fb,"fhat_GeV":fhat,"Bhat":bag,"frequency_ps_inv":abs(W),
            "mix_phase":np.angle(W/lc**2),"tag_phase_with_common_delta":measured_phase,
            "control_shift_rad":.013 if main_calibration else .015,"Belle_signal_strength":br_n,"detector_Bs_mumu_e9":acc,
            "nuclear_Vud_holdout_two_width":rims(old["ud_beta"]),"nuclear_Vud_miss":ud-rims(old["ud_beta"])[1]}


def first_row():
    old=read(1381); d=read(1445);p=read(1633,"results.json")
    usmax=rims(d["K_product"])[1]/rims(d["K_form"])[0]
    udmax=rims(old["ud_beta"])[1];cbmax=rims(d["BelleII_Vcb"])[1]
    rho=rims(d["ub_over_cb_holdout"])[1];umax=cbmax*rho
    rb=math.sqrt(p["right_doublet_probability"]);smax2=2*p["left_difference"]/p["k"]
    ls=.1*.0027/(1600*rb);lb=4.8*rb/1600
    Gmin=1-smax2-ls**2-lb**2
    deficit=(smax2*umax**2+(ls*usmax+lb*umax)**2)/Gmin
    ubmin=math.sqrt(1-udmax**2-usmax**2-deficit)
    return {"us_max":usmax,"ud_max":udmax,"ub_max":umax,"max_first_row_deficit":deficit,
            "ub_min_with_nuclear":ubmin,"cb_min_with_nuclear":ubmin/rho}


def rounded(x):
    if isinstance(x,dict):return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [rounded(v) for v in x]
    if isinstance(x,(float,np.floating)): return 0. if abs(x)<1e-9 else float(f"{x:.10g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--check",action="store_true");args=p.parse_args()
    data=rounded({"witnesses":[audit(M,m5) for M,m5 in ((1600,2100),(2100,2650),(2600,3200))],
                  "first_row":first_row()})
    data["primary_calibration_SVD"]=rounded([audit(M,1.1*M+340,True) for M in (1600,2100,2600)])
    target=HERE/"independent_results.json"
    if args.check:
        assert data==json.loads(target.read_text(encoding="utf-8"))
        print("1640 independent triple mass SVD, common clocks and tree tags: PASS")
    else:
        target.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
        for item in data["witnesses"]:
            print({k:v for k,v in item.items() if k not in ("weak_basis_mass_re_im","singular_masses_GeV")})
        print(data["first_row"])
