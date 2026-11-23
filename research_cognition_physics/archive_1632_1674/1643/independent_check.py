"""Reconstruct three joint tree, kaon and strange-bottom witnesses from actual SVDs."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())


def bound(row):
    return row[0]-2*sum(row[1:]),row[0]+2*sum(row[1:])


def check_inside(value,box):
    assert box[0]<value<box[1], (value,box)


def pair(z):
    return [z.real,z.imag]


# Candidate choices. Two K_l3 endpoints are compared using the same form factor.
TREE=(.22500186439586628,.04235828886614449,.004178645196645154,248.65)
POLES={1600:[.0003248968497301052,-.00039629623648091485],
       2100:[.00029053824221458044,-.000355409039666075],
       2600:[.0002680567794533893,-.0003243443138708945]}
CAL={"fplus":.9640001,"fhat":.25068110668969057,"fBs":.23362112917039116,
     "etaB":.5511,"S0":1.2696330283011683/.5511,"etaZ":.5438062016051479,
     "U4":.047344269601343705,"U5":.04313821790145036,"k":.6933008318530638,
     "control_shift":.0197}


def witness(md,positive=False):
    author=runpy.run_path(str(HERE.parent/"1640"/"check.py"))
    cfg=read("../1640/inputs.json")
    cal=dict(CAL)
    if positive:
        initial=[.00018855203489523558,-.00020278940620255584]
        tree=[.22500143914676862,.04339263172912409,.004190496690697416,
              math.degrees(1.0344684651107778)]
        cal.update(fplus=.96400005,fBs=.23122270679805626,etaB=.55088,S0=2.2971,
                   etaZ=.5709412689861492,U4=.04713966782241626,U5=.04228637969311237,
                   k=.6993102350574323,control_shift=.0195)
    else:
        initial=POLES[md]
        tree=list(TREE)
    cfg["F"]=initial
    cfg["ordinary_s12_s23_s13_delta_deg"]=tree
    mat,_,_,_,_,_,_,_=author["build"](cfg,md)
    ul,sv,vh=np.linalg.svd(mat)
    order=np.argsort(sv); m=sv[order]; ul=ul[:,order];ur=vh.conj().T[:,order]
    np.testing.assert_allclose(ul@np.diag(m)@ur.conj().T,mat,atol=3e-10)
    K=ul[:3,:3]
    lm=np.array([K[i,1].conjugate()*K[i,2] for i in range(3)])
    lu,lc,lt=lm
    F=-ur[4,1].conjugate()*ur[4,2]
    fl=-ul[3,1].conjugate()*ul[3,2]-2*ul[4,1].conjugate()*ul[4,2]
    ms,mb=m[1:3]
    np.testing.assert_allclose(fl,2*ms*mb/md**2*F,atol=2e-13)
    p=read("../1633/results.json")
    np.testing.assert_allclose(abs(ur[4,2])**2,p["right_doublet_probability"],rtol=1e-10)
    np.testing.assert_allclose(abs(ul[3,2])**2/2+abs(ul[4,2])**2,p["left_difference"]/p["k"])
    g=read("../1636/inputs.json");r=read("../1464/inputs.json");t=read("../1445/inputs.json")
    W=read("../1404/inputs.json");lat=read("../1407/inputs.json");c=read("../1468/inputs.json")
    H=read("../1406/inputs.json");rates=read("../1464/results.json");n=read("../1381/inputs.json")
    nuclear=bound(n["ud_beta"]);Kproduct=bound(t["K_product"])
    ud,us,ub=abs(K[0]);cb=abs(K[1,2])
    check_inside(ud,nuclear)
    check_inside(cal["fplus"],bound(t["K_form"]))
    check_inside(us*cal["fplus"],Kproduct)
    check_inside(cb,bound(t["BelleII_Vcb"]))
    check_inside(cb,bound(cfg["inclusive_holdout"]["Vcb"]))
    check_inside(ub/cb,bound(t["ub_over_cb_holdout"]))
    weakq=K[0,1]*K[1,2]*K[0,2].conjugate()*K[1,1].conjugate()
    gs=np.angle(weakq,deg=True)%360
    gd=np.angle(-K[0,0]*K[0,2].conjugate()/(K[1,0]*K[1,2].conjugate()),deg=True)%360
    gamma_input=cfg["gamma_Bpm"]
    gamwin=[gamma_input["deg"]-2*gamma_input["down"],gamma_input["deg"]+2*gamma_input["up"]]
    for v in (gs,gd):check_inside(v%180,gamwin)
    assert ((gs<180) == positive) and abs(gs-gd)<gamma_input["transfer_cap_deg"]
    d=W["const_GeV"];gf,mw,mbs=d["GF"],d["MW"],d["M_Bs"]
    f=cal["fBs"];hat=cal["fhat"]
    check_inside(f,bound(r["fBs_GeV"]))
    check_inside(cal["etaB"],bound(W["short_distance"]["etaB"]))
    check_inside(cal["S0"],bound(W["short_distance"]["S0_published"]))
    check_inside(cal["etaZ"],c["eta_Z_RGI_range"])
    for v,spec in zip([cal["U4"],cal["U5"]],g["LR_F2B4_B5_GeV2"]):check_inside(v,bound(spec))
    check_inside(cal["k"],H["sm_coefficient_ratio"]["range"])
    assert abs(cal["control_shift"])<c["extra_decay_phase_cap_rad"]
    xt=W["short_distance"]["xt"][0];mt=g["mt_GeV"];mz=g["MZ_GeV"]
    pref=2*mbs/(3*lat["hbar_GeV_ps"]); norm=gf**2*mw**2/(4*np.pi**2)
    alpha=c["alpha_s_ref"]
    atol=alpha/(1+alpha*(23/3)/(2*np.pi)*np.log(mt/mz))
    ah=atol/(1+atol*7/(2*np.pi)*np.log(md/mt))
    rg=(atol/alpha)**(6/23)*(ah/atol)**(2/7)
    h1=-(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)/(2*(xt-1)**2)*np.log(xt)
    light=(1-(ms/md)**2)*(1-(mb/md)**2)
    vr=gf/np.sqrt(2)*(F**2+fl**2)+rg*md**2/(16*np.pi**2*g["vev_GeV"]**4)*light**2*F**2
    # hard matching is (1/4)*(2 M^2/v^2)^2/(16 pi^2 M^2).
    LR= -norm*lt*F*xt*(np.log(md/mw)+h1)*light
    ulr=g["LR_B4_B5_weights"][0]*cal["U4"]+g["LR_B4_B5_weights"][1]*cal["U5"]
    axis=lc.conjugate()**2/abs(lc)**2
    aa=pref*axis*(norm*lt**2*cal["etaB"]*cal["S0"]+cal["etaZ"]*vr)
    bb=pref*axis*LR*ulr
    if positive:
        ov=(aa.conjugate()*bb).real
        om=c["omega_ps_inverse"][0]
        hat=math.sqrt((-ov+math.sqrt(ov**2+abs(aa)**2*(om**2-abs(bb)**2)))/abs(aa)**2)
    bhat=(hat/f)**2
    check_inside(hat,bound(lat["HPQCD_inputs"]["fBs_sqrtB_GeV"]))
    check_inside(bhat,bound(W["bag"]["s"]))
    clk=aa*hat**2+bb
    phi=np.angle(clk)
    check_inside(abs(clk),bound(c["omega_ps_inverse"]))
    for tag in ("CMS_KK_13TeV_final","LHCb_KK_2015_18"):
        check_inside(float(phi)+cal["control_shift"],bound(c["phase_inputs_rad"][tag]))
    ref=r["HPQCD_tableI_eq15"]
    L=-ref["Xt_ref"]/(np.pi*ref["sine_squared_ref"]*ref["alpha_inverse_ref"])*lt
    N=L+fl+F;P=cal["k"]*L+fl-F
    nr=1.1199 if positive else 1.
    assert abs(nr-1)<=2*ref["hadronic_rate_width"]
    bel=rates["Hnu"]*abs(N)**2*nr/(H["BelleII"]["reference_SD_in_1e-6"]*1e-6)
    check_inside(bel,bound([H["BelleII"]["mu"],H["BelleII"]["width"]]))
    tau=1.6329 if positive else W["tau_ps"]["s_heavy"][0]
    tr=1.0199 if positive else 1.
    check_inside(tau,bound(W["tau_ps"]["s_heavy"]))
    assert abs(tr-1)<=r["mu_rate_transfer_width"]
    scale=1e9*rates["Hmu_heavy"]*(f/r["fBs_GeV"][0])**2*tr*tau/W["tau_ps"]["s_heavy"][0]
    smax=scale*(mbs/g["higgs_mass_GeV"])**4*abs(F)**2
    ad=np.angle(P/lt)-np.angle(clk/(axis*lt**2))/2
    muons=[]
    ws=(.9999,.99) if positive else (rates["calculated_light_response"],rates["run3_central_light_response"])
    assert .75<=ws[0]<=1 and .78<=ws[1]<=1
    for w,window in zip(ws,([2.53,5.25],[2.45,4.23])):
        val=scale*abs(P)**2*(1-(1-w)*np.sin(ad)**2)
        check_inside(val,window);check_inside(val+smax,window)
        muons.append([val,val+smax])
    sx=author["read"]("../1635/inputs.json")["witness"]["W_mass_GeV"]
    rate=runpy.run_path(str(HERE.parent/"1635"/"check.py"))["width_factor"]
    branch=abs(ur[4,1]/ur[4,2])**2*rate(ms,md,sx)/rate(mb,md,sx)
    return {"M":md,"U_input":tree,"mass_eigenvalues":m.tolist(),"F_in_input_basis":initial,
        "Kud_Kus_Kub_Kcb":[ud,us,ub,cb],"Kus_fplus":us*cal["fplus"],
        "distance_to_nuclear_upper":nuclear[1]-ud,"distance_to_Kprod_upper":Kproduct[1]-us*cal["fplus"],
        "rho":ub/cb,"gamma_s_gamma_d_deg":[gs,gd],"signed_quartet_imag":weakq.imag,
        "z_F_over_L":pair(F/L),"absolute_clock":abs(clk),"mix_phase":phi,
        "reported_control_phase":phi+cal["control_shift"],"Belle":bel,
        "fhat":hat,"Bhat":bhat,"CMS_Run3":muons,"Ws_over_Wb":branch}


def normalize(x):
    if isinstance(x,dict):return {k:normalize(v) for k,v in x.items()}
    if isinstance(x,list):return [normalize(v) for v in x]
    if isinstance(x,float):return float(f"{x:.11g}")
    return x


def main():
    p=argparse.ArgumentParser();p.add_argument("--check",action="store_true")
    a=p.parse_args()
    out=normalize({"shared":CAL,"witnesses":[witness(m) for m in sorted(POLES)],
                   "positive_branch":witness(1621.1392891179632,True)})
    if a.check:assert out==read("independent_results.json")
    else:(HERE/"independent_results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    for row in out["witnesses"]:print(row["M"],row["Kud_Kus_Kub_Kcb"],row["absolute_clock"],row["Belle"],row["CMS_Run3"])

if __name__=="__main__":
    main()
