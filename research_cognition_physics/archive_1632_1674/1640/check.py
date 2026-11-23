"""Actual five-state witnesses with independently calibrated tree-current readouts."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())


def bounds(row,n):
    return row[0]-n*sum(row[1:]),row[0]+n*sum(row[1:])


def inside(val,limits):
    assert limits[0]<val<limits[1], (val,limits)


def enc(v):
    a=np.asarray(v)
    return np.stack((a.real,a.imag),axis=-1).tolist()


def build(cfg,md):
    prev=read(cfg["original_doublet"])
    wp=prev["witness"]
    old=runpy.run_path(str((HERE/"../1635/check.py").resolve()))
    p=read("../1633/results.json")
    m5=cfg["m5_slope_intercept"][0]*md+cfg["m5_slope_intercept"][1]
    m1,ms,mb=wp["light_masses_GeV"]
    rb=math.sqrt(p["right_doublet_probability"])
    F=complex(*cfg["F"])
    rs=-F/rb
    sb=math.sqrt(2*(p["left_difference"]/p["k"]-(mb*rb/md)**2))
    def trial(r5):
        r4=math.sqrt(1-rb**2-abs(rs)**2-r5**2)
        m4=math.sqrt((md**2-ms**2*abs(rs)**2-mb**2*rb**2-m5**2*r5**2)/r4**2)
        masses=np.array([m1,ms,mb,m4,m5])
        r=np.array([0.,rs,rb,r4,r5])
        ell=masses*r/md
        s4=-ell[2].real*sb*(m5**2-mb**2)/((m5**2-m4**2)*ell[3].real)
        s5=ell[2].real*sb*(m4**2-mb**2)/((m5**2-m4**2)*ell[4].real)
        return masses,r,ell,np.array([0.,0.,sb,s4,s5])
    lo,hi=1e-9,.01
    assert np.linalg.norm(trial(lo)[3])>1>np.linalg.norm(trial(hi)[3])
    for _ in range(75):
        mid=(lo+hi)/2
        if np.linalg.norm(trial(mid)[3])>1:
            lo=mid
        else:
            hi=mid
    m,r,ell,s=trial((lo+hi)/2)
    ts=m*s/np.linalg.norm(m*s)
    s12,s23,s13,deg=cfg["ordinary_s12_s23_s13_delta_deg"]
    u=old["ckm"](s12,s23,s13,math.radians(deg))
    vleft=np.column_stack((s,ell))
    vright=np.column_stack((ts,r))
    wl=np.column_stack((old["polar_columns"](vleft)@u.conj().T,vleft))
    wr=np.column_stack((old["polar_columns"](vright),vright))
    M=wl.conj().T@np.diag(m)@wr
    np.testing.assert_allclose(wl.conj().T@wl,np.eye(5),atol=1e-12)
    np.testing.assert_allclose(wr.conj().T@wr,np.eye(5),atol=1e-12)
    np.testing.assert_allclose(np.linalg.svd(M,compute_uv=False)[::-1],m,rtol=1e-10,atol=1e-10)
    assert max(np.max(abs(M[:4,4])),np.max(abs(M[3,:3])),abs(M[4,3]))<1e-9
    np.testing.assert_allclose([M[4,4],M[3,3]],[md,np.linalg.norm(m*s)],atol=1e-9)
    hard=np.diag([0.,0.,0.,M[3,3],M[4,4]])
    yuk=np.linalg.svd(np.sqrt(2)*(M-hard)/wp["Higgs_vev_GeV"],compute_uv=False)[0]
    assert yuk<math.sqrt(4*math.pi)
    K=wl[:,:3].conj().T[:,:3]
    widths=abs(r)**2*np.array([old["width_factor"](mi,md,wp["W_mass_GeV"]) for mi in m])
    return M,m,K,yuk,widths,r,s,ell


def solve_point(cfg,md):
    n=cfg["working_widths"]
    mat,m,K,yuk,widths,r,s,ell=build(cfg,md)
    t=read(cfg["exclusive_tree"])
    nuc=read(cfg["nuclear_holdout"])
    W=read("../1404/inputs.json")
    H=read("../1406/inputs.json")
    G=read(cfg["operator_matching"])
    C=read("../1468/inputs.json")
    lat=read("../1407/inputs.json")
    ff=read("../1464/inputs.json")
    rare=read("../1464/results.json")
    F=complex(*cfg["F"])
    cb,ub,us=abs(K[1,2]),abs(K[0,2]),abs(K[0,1])
    ud=abs(K[0,0])
    rho=ub/cb
    # The Ds and Dd tree phase conventions differ only by the measured Cabibbo reference.
    gamma_s=np.angle(K[1,2]*K[0,1]*K[0,2].conj()*K[1,1].conj(),deg=True)
    gamma_d=np.angle(-K[0,0]*K[0,2].conj()/(K[1,0]*K[1,2].conj()),deg=True)
    gd=cfg["gamma_Bpm"]
    glimits=(gd["deg"]-n*gd["down"],gd["deg"]+n*gd["up"])
    inside(gamma_s,glimits);inside(gamma_d,glimits)
    assert abs(gamma_s-gamma_d)<gd["transfer_cap_deg"]
    inside(rho,bounds(t["ub_over_cb_holdout"],n))
    inside(cb,bounds(t["BelleII_Vcb"],n))
    inside(us*t["K_form"][0],bounds(t["K_product"],n))
    assert cb<bounds(cfg["inclusive_holdout"]["Vcb"],n)[0] and ud>bounds(nuc["ud_beta"],n)[1]
    lt,lc=[K[j,1].conjugate()*K[j,2] for j in (2,1)]
    ms,mb=m[1:3];e=2*ms*mb/md**2
    fL=-(s[1]*s[2].conjugate()+2*ell[1]*ell[2].conjugate())
    np.testing.assert_allclose([fL,-r[1]*r[2].conjugate()],[e*F,F],atol=1e-15)
    rf=ff["HPQCD_tableI_eq15"]
    co=rf["Xt_ref"]/(np.pi*rf["sine_squared_ref"]*rf["alpha_inverse_ref"])
    L=-co*lt
    N=L+fL+F
    P=cfg["k"]*L+fL-F
    weak=W["const_GeV"]
    gf,mw,mbs=weak["GF"],weak["MW"],weak["M_Bs"]
    pre=2*mbs/(3*lat["hbar_GeV_ps"]); norm=gf**2*mw**2/(4*np.pi**2)
    mt=G["mt_GeV"];xt=W["short_distance"]["xt"][0]
    run=runpy.run_path(str((HERE/"../1636/check.py").resolve()))["alpha_step"]
    az=C["alpha_s_ref"];at=run(az,G["MZ_GeV"],mt,5);am=run(at,mt,md,6)
    qcd=(at/az)**(6/23)*(am/at)**(6/21)
    light=(1-ms**2/md**2)*(1-mb**2/md**2)
    h1=-(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)/(2*(xt-1)**2)*np.log(xt)
    prod=-2*(md/G["vev_GeV"])**2*light*F
    hard=G["RR_hard_Qd_color_weak_factor"]*prod**2/((4*np.pi)**2*md**2)
    values=[cfg[k] for k in ("etaB","S0","eta_Z","f2B4","f2B5","fBs_GeV")]
    box=[bounds(W["short_distance"]["etaB"],n),bounds(W["short_distance"]["S0_published"],n),
         C["eta_Z_RGI_range"],*(bounds(a,n) for a in G["LR_F2B4_B5_GeV2"]),bounds(ff["fBs_GeV"],n)]
    for value,limits in zip(values,box): inside(value,limits)
    a4,a5=G["LR_B4_B5_weights"]; baglr=a4*cfg["f2B4"]+a5*cfg["f2B5"]
    axis=lc.conjugate()**2/abs(lc)**2
    c0=pre*axis*(norm*lt**2*cfg["etaB"]*cfg["S0"]+
                 cfg["eta_Z"]*(gf/np.sqrt(2)*(F**2+fL**2)+qcd*hard))
    c1=-pre*axis*norm*lt*F*xt*(np.log(md/mw)+h1)*light*baglr
    omega=C["omega_ps_inverse"][0]
    ab=(c0.conjugate()*c1).real
    y=(-ab+np.sqrt(ab**2+abs(c0)**2*(omega**2-abs(c1)**2)))/abs(c0)**2
    fhat=np.sqrt(y);fb=cfg["fBs_GeV"];bhat=y/fb**2
    inside(fhat,bounds(lat["HPQCD_inputs"]["fBs_sqrtB_GeV"],n))
    inside(bhat,bounds(W["bag"]["s"],n))
    mix=c0*y+c1;phi=np.angle(mix)
    assert abs(abs(mix)-omega)<1e-11
    phidec=phi+cfg["control_decay_phase_rad"]
    assert abs(cfg["control_decay_phase_rad"])<C["extra_decay_phase_cap_rad"]
    for tag in ("CMS_KK_13TeV_final","LHCb_KK_2015_18"):
        inside(phidec,bounds(C["phase_inputs_rad"][tag],n))
    angle=np.angle(P/lt)-np.angle(mix/(axis*lt**2))/2
    scale=1e9*rare["Hmu_heavy"]*(fb/ff["fBs_GeV"][0])**2
    scalarmax=scale*((mbs/G["higgs_mass_GeV"])**2*abs(F))**2
    mu=lambda w: scale*abs(P)**2*(1-(1-w)*math.sin(angle)**2)
    bid=H["BelleII"]
    bel=rare["Hnu"]*abs(N)**2/(bid["reference_SD_in_1e-6"]*1e-6)
    inside(bel,bounds([bid["mu"],bid["width"]],n))
    cms=H["CMS"];r3=W["LHCb_run3_2024_25"]["s_in_1e9"]
    windows=[(cms["branching_in_1e-9"]-n*sum(cms["lower_widths"]),
              cms["branching_in_1e-9"]+n*sum(cms["upper_widths"])),
             (r3[0]-n*(r3[1]+r3[3]),r3[0]+n*(r3[2]+r3[4]))]
    mu_ranges=[]
    for weight,limits in zip([rare["calculated_light_response"],rare["run3_central_light_response"]],windows):
        val=mu(weight)
        inside(val,limits);inside(val+scalarmax,limits)
        mu_ranges.append([val,val+scalarmax])
    assert abs(mat[4,4]-md)<1e-9 and abs(m[-1]-(1.1*md+340))<1e-9
    return {"M_GeV":md,"physical_masses_GeV":m.tolist(),"weak_basis_mass_GeV":enc(mat),
        "K":enc(K),"F_L_R":enc([fL,F]),"abs_Kus_cb_ub_ud":[us,cb,ub,ud],
        "Kus_times_central_form":us*t["K_form"][0],"Kub_over_Kcb":rho,"gamma_s_d_deg":[gamma_s,gamma_d],
        "abs_lambda_t":abs(lt),"Bhat":bhat,"fBs_sqrt_Bhat_GeV":fhat,"mix_2M12_ps_inverse":enc(mix),
        "mix_frequency":abs(mix),"mix_phase":phi,"KK_phase":phidec,"Belle_strength":bel,
        "CMS_Run3_dimuons_e9":mu_ranges,"Ws_over_Wb":widths[1]/widths[2],
        "Wb_fraction":widths[2]/sum(widths),"largest_yukawa":yuk,
        "dist_below_inclusive_2width":bounds(cfg["inclusive_holdout"]["Vcb"],n)[0]-cb,
        "dist_above_nuclear_2width":ud-bounds(nuc["ud_beta"],n)[1]}


def normalize(obj):
    if isinstance(obj,dict):
        return {k:normalize(x) for k,x in obj.items()}
    if isinstance(obj,list):
        return [normalize(x) for x in obj]
    if isinstance(obj,float):
        return 0. if abs(obj)<1e-12 else float(f"{obj:.11g}")
    return obj


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    cfg=read("inputs.json")
    points=[solve_point(cfg,md) for md in cfg["spectra_GeV"]]
    out=normalize({"round":1640,"common_F":cfg["F"],"points":points})
    if args.check:
        assert out==read("results.json")
    else:
        (HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    for b in points:
        print(b["M_GeV"],b["abs_Kus_cb_ub_ud"][1],b["fBs_sqrt_Bhat_GeV"],b["mix_phase"],b["Belle_strength"],b["CMS_Run3_dimuons_e9"])

if __name__=="__main__":
    main()
