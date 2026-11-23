"""Joint external tree windows and signed, phase-uniform Bs-clock enclosure."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def read(s):
    return json.loads((HERE/s).read_text(encoding="utf-8-sig"))
def width(a,n=2):
    return a[0]-n*sum(a[1:]),a[0]+n*sum(a[1:])

def tree(cfg):
    cc,old,tab,pole=[read(cfg[t]) for t in ("tree","tree_N","semileptonic","fixed_bottom")]
    ud=width(old["ud_beta"]); vp=width(old["K_product"]); ff=width(old["K_form"])
    us=[vp[0]/ff[1],vp[1]/ff[0]]; ratio=width(old["ub_over_cb"])
    rb=math.sqrt(pole["right_doublet_probability"]); hmin=math.sqrt(1-2*pole["left_difference"]/pole["k"])
    ms,mb=.1,4.8; M=1600; el=mb*rb/M; es=ms/M
    dh=(2*es*el+es**2)/(math.sqrt(hmin**2-es**2)+hmin)
    u3min=math.sqrt(max(0,1-ud[1]**2-(us[1]+dh)**2))
    ubmin=hmin*u3min-dh
    gamma=cc["gamma_Bpm"]; cap=gamma["transfer_cap_deg"]
    angles=[gamma["deg"]-2*gamma["down"]-cap,gamma["deg"]+2*gamma["up"]+cap]
    fmax=.0027  # 1639 Belle envelope applies under independent |K_cb| <= .046, |K_ub| <= .005.
    eps=ms*mb/M**2*fmax
    out={}
    for channel,cb in (("exclusive",width(tab["BelleII_Vcb"])),("inclusive",width(cc["inclusive_holdout"]["Vcb"]))):
        cb=[max(cb[0],ubmin/ratio[1]),cb[1]]
        assert cb[0]<cb[1]<.046 and ratio[1]*cb[1]<.005
        s13hi=(ratio[1]*cb[1]+dh)/hmin; c13lo=math.sqrt(1-s13hi**2)
        s12hi=(us[1]+dh)/c13lo
        s12lo=max(us[0]-dh,math.sqrt(max(0,1-(ud[1]/c13lo)**2)))
        s23hi=(cb[1]+dh)/(hmin*c13lo)
        cshi=s12hi*s23hi*s13hi+math.sqrt(1-s12lo**2)+dh
        cslo=math.sqrt((1-s12hi**2)*(1-s23hi**2))-s12hi*s23hi*s13hi-dh
        rlo=max(ratio[0],ubmin/cb[1])
        q=[max(us[0],s12lo*c13lo-dh)*rlo/cshi,us[1]*ratio[1]/cslo]
        lc=[cb[0]*cslo,cb[1]*cshi]
        branches={}
        for sign in (0,1):
            gs=[math.radians(i+180*sign) for i in angles]
            # q+cos(gamma), and sin(gamma), do not change sign on each box.
            assert all(math.cos(v)*math.cos(gs[0])>0 and abs(math.cos(v))>q[1] for v in gs)
            corners=[1+r*cmath.exp(-1j*g) for r,g in itertools.product(q,gs)]
            lo=lc[0]*min(map(abs,corners))-eps; hi=lc[1]*max(map(abs,corners))+eps
            pe=2*math.asin(eps/(lo+eps))+1e-10
            phi=[2*min(map(cmath.phase,corners))-pe,2*max(map(cmath.phase,corners))+pe]
            branches[str(sign)]={"abs_lambda_t":[lo,hi],"top_phase_relative_to_charm_rad":phi}
        # Largest tightening of a single upper endpoint before the row becomes impossible.
        ud_needed=math.sqrt(1-(us[1]+dh)**2-s13hi**2)
        us_needed=math.sqrt(1-ud[1]**2-s13hi**2)-dh
        out[channel]={"actual_abs_Kcb":cb,"abs_Kcs":[cslo,cshi],"q_u_over_c":q,"branches":branches,
                      "max_slack_Kud_upper":ud[1]-ud_needed,"max_slack_Kus_upper":us[1]-us_needed}
    one_f=width(old["K_form"],1); one_p=width(old["K_product"],1)
    one_us=one_p[1]/one_f[0]
    ub1=width(tab["BelleII_Vcb"],1)[1]*width(old["ub_over_cb"],1)[1]
    ub2=width(cc["inclusive_holdout"]["Vcb"],1)[1]*width(old["ub_over_cb"],1)[1]
    rowdef=1-width(old["ud_beta"],1)[1]**2-(one_us+dh)**2-((max(ub1,ub2)+dh)/hmin)**2
    assert rowdef>0
    return {"Kud":ud,"Kus":us,"ub_over_cb":ratio,"operator_perturbation":dh,"bare_s13_floor":u3min,
            "physical_Kub_floor":ubmin,"physical_Kcb_floor":ubmin/ratio[1],"gamma_principal_deg":angles,
            "one_width_required_unit_row_deficit":rowdef,"separate_channels":out}

def square_min(lo,hi):
    return np.where(lo*hi<=0,0,np.minimum(lo*lo,hi*hi))

def covers(cfg,trees):
    data=read(cfg["calibration"])["common_witness"]; q=read(cfg["mixing"])
    weak,hp,clk,ff=[read(cfg[k]) for k in ("weak","lattice","control","nu")]
    wr=weak["const_GeV"]; gf,mw,bs=wr["GF"],wr["MW"],wr["M_Bs"]
    pref=2*bs/(3*hp["hbar_GeV_ps"]); norm=gf**2*mw**2/(4*math.pi**2)
    ref=data["reference_top_modulus"]; oldrad=1.01e-8
    # Convert the archived dimensionless N,P bounds into bounds divided by c alone.
    nabs=data["nmin"]*(ref+oldrad-2e-10)-1e-9
    pabs=data["pmin"]*(ref+oldrad-2e-10)-1e-9
    pmax=data["pmax"]*(ref-oldrad+2e-10)+1e-9
    g=ff["HPQCD_tableI_eq15"]; c=g["Xt_ref"]/(math.pi*g["alpha_inverse_ref"]*g["sine_squared_ref"])
    hat=width(hp["HPQCD_inputs"]["fBs_sqrtB_GeV"]); et=width(weak["short_distance"]["etaB"]); s=width(weak["short_distance"]["S0_published"])
    eb=clk["eta_Z_RGI_range"]; mass=q["mt_GeV"]; mz=q["MZ_GeV"]; al=clk["alpha_s_ref"]; xt=weak["short_distance"]["xt"][0]
    at=al/(1+al*23/(6*math.pi)*math.log(mass/mz))
    h1=-(xt-7)/(4*(xt-1))-(xt*xt-2*xt+4)*math.log(xt)/(2*(xt-1)**2)
    U=[sum(co*width(b)[i] for co,b in zip(q["LR_B4_B5_weights"],q["LR_F2B4_B5_GeV2"])) for i in (0,1)]
    def bcoef(M):
        a=al
        am=at/(1+at*7/(2*math.pi)*math.log(M/mass))
        rr=(at/a)**(6/23)*(am/at)**(2/7)
        light=(1-(.1/M)**2)*(1-(4.8/M)**2)
        return rr*(2*M/q["vev_GeV"]**2*light)**2/(64*math.pi**2)
    def dcoef(M):
        return xt*(math.log(M/mw)+h1)*(1-(.1/M)**2)*(1-(4.8/M)**2)
    counts=cfg["bins"]; MM=np.linspace(1600,2600,counts["mass"]+1)
    ks=np.linspace(.60,.70,counts["k"]+1)
    phi_obs=read("../1639/results.json")["common_witness"]["2"]["mix_phase_range_rad"]
    whi=width(clk["omega_ps_inverse"])[1]
    bins_out={}
    for channel,inp in trees["separate_channels"].items():
        out={}
        for sign,branch in inp["branches"].items():
            lrange=branch["abs_lambda_t"]; phi=branch["top_phase_relative_to_charm_rad"]
            imcap=whi*max(abs(math.sin(phi_obs[0]-phi[1])),abs(math.sin(phi_obs[1]-phi[0])))
            min_gap=1e6; feasible=0; total=0; tleft=1.; tright=0.
            le=np.linspace(*lrange,counts["lambda"]+1)
            # k enters only through the two dimuon disks. Split k and all x values in array form.
            xl=np.linspace(-2.,2.5,counts["x"]+1)[:-1][None,:]
            xh=xl+4.5/counts["x"]
            kl,kh=ks[:-1,None],ks[1:,None]
            for tl,th in zip(le,le[1:]):
                nl=nabs/th; pl=pabs/th; pu=pmax/tl
                assert .70+pu+2e-5<2.5 and .60-pu-2e-5>-2.
                # Bound the tiny finite-mass factors by a common outward addition.
                pad=3.e-6; nl-=pad; pl-=pad; pu+=pad
                ylow=np.maximum.reduce(np.broadcast_arrays(0*xl,
                    nl**2-np.maximum((1+xl)**2,(1+xh)**2),
                    pl**2-np.maximum((kl-xh)**2,(kh-xl)**2)))
                yhigh=pu**2-square_min(kl-xh,kh-xl)
                ar=pref*hat[0]**2*norm*tl**2*et[0]*s[0]
                for ml,mh in zip(MM,MM[1:]):
                    br=[pref*hat[i]**2*(gf/math.sqrt(2)*(1+((2*.1*4.8/1600**2)**2 if i else 0))+bcoef(mm))*eb[i]*(c*l)**2
                        for i,(mm,l) in enumerate(((ml,tl),(mh,th)))]
                    dr=[pref*norm*l**2*c*dcoef(mm)*U[i] for i,(mm,l) in enumerate(((ml,tl),(mh,th)))]
                    # Upper bound Re Z and |Im Z| from a circumscribed clock sector.
                    dl=dr[0]+2*np.minimum(br[0]*xl,br[1]*xl)
                    du=dr[1]+2*np.maximum(br[0]*xh,br[1]*xh)
                    im_eps=cfg["one_sided_residual_ps_inverse"]
                    yc=np.minimum(yhigh,(imcap+im_eps)**2/np.maximum(square_min(dl,du),1e-40))
                    diff=ar+np.minimum(dr[0]*xl,dr[1]*xl)-(whi+im_eps)
                    yr=np.maximum(ylow,square_min(xl,xh)+np.minimum(diff/br[0],diff/br[1]))
                    gap=yr-yc
                    min_gap=min(min_gap,float(np.min(gap)))
                    ok=int(np.count_nonzero(gap<=1e-8))
                    feasible+=ok; total+=gap.size
                    if ok: tleft=min(tleft,tl); tright=max(tright,th)
            out[sign]={"minimum_y2_gap":min_gap,"surviving_outer_cells":feasible,"enclosed_cells":total,
                       "possible_abs_lambda_t_outer":[tleft,tright] if feasible else [],"abs_imag_clock_upper":imcap}
        bins_out[channel]=out
    return bins_out

def witnesses(cfg):
    from runpy import run_path
    build=run_path(str(HERE.parent/"1640"/"check.py"))["build"]
    common=cfg["witness_common"]; base=read(cfg["tree"]); weak=read(cfg["weak"])
    lat,cl,nu,rare=[read(cfg[k]) for k in ("lattice","control","nu","rare")]
    had=read("../1464/results.json"); mc=read(cfg["mixing"]); tr=read(cfg["tree_N"]); tab=read(cfg["semileptonic"])
    inside=lambda value,limits: limits[0]<=value<=limits[1]
    d=weak["const_GeV"]; gt,mw,mbs=d["GF"],d["MW"],d["M_Bs"]
    pre=2*mbs/(3*lat["hbar_GeV_ps"]);norm=gt**2*mw**2/(4*math.pi**2)
    ref=nu["HPQCD_tableI_eq15"];c=ref["Xt_ref"]/(math.pi*ref["sine_squared_ref"]*ref["alpha_inverse_ref"])
    az=cl["alpha_s_ref"];mt,mz=mc["mt_GeV"],mc["MZ_GeV"];xt=weak["short_distance"]["xt"][0]
    h1=-(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)*math.log(xt)/(2*(xt-1)**2)
    at=az/(1+az*23/(6*math.pi)*math.log(mt/mz)); ledger={}
    for name,vals in cfg["physical_witnesses"].items():
        p=dict(zip(cfg["witness_keys"],vals)); M=p["M_GeV"]; F=complex(p["ReF"],p["ImF"])
        state=dict(base); state["F"]=[F.real,F.imag]
        state["ordinary_s12_s23_s13_delta_deg"]=[p["s12"],p["s23"],p["s13"],math.degrees(p["delta_rad"])]
        mat,m,K,yuk,_,_,_,_=build(state,M)
        la=K[:,1].conjugate()*K[:,2];lt,lc=la[2],la[1]
        ub,cb,us,ud=abs(K[0,2]),abs(K[1,2]),abs(K[0,1]),abs(K[0,0])
        gamma_s=math.degrees(cmath.phase(lc/la[0]))%180
        gamma_d=math.degrees(cmath.phase(-K[0,0]*K[0,2].conjugate()/(K[1,0]*K[1,2].conjugate())))%180
        g=base["gamma_Bpm"];gb=(g["deg"]-2*g["down"],g["deg"]+2*g["up"])
        assert inside(gamma_s,gb) and inside(gamma_d,gb) and abs(gamma_s-gamma_d)<g["transfer_cap_deg"]
        assert inside(ud,width(tr["ud_beta"])) and inside(common["fplus"],width(tr["K_form"]))
        assert inside(us*common["fplus"],width(tr["K_product"])) and inside(ub/cb,width(tr["ub_over_cb"]))
        assert inside(cb,width(tab["BelleII_Vcb"])) and inside(cb,width(base["inclusive_holdout"]["Vcb"]))
        assert inside(p["k"],rare["sm_coefficient_ratio"]["range"]) and inside(p["etaZ"],cl["eta_Z_RGI_range"])
        assert all(inside(v,width(b)) for v,b in zip([p["U4"],p["U5"]],mc["LR_F2B4_B5_GeV2"]))
        assert all(inside(common[key],width(weak["short_distance"][key if key=="etaB" else "S0_published"])) for key in ("etaB","S0"))
        assert inside(common["tau"],width(weak["tau_ps"]["s_heavy"])) and inside(p["fBs"],width(nu["fBs_GeV"]))
        assert abs(common["mu_transfer"]-1)<=nu["mu_rate_transfer_width"]
        assert abs(common["nu_transfer"]-1)<=2*ref["hadronic_rate_width"]
        assert inside(common["CMS_weight"],rare["CMS"]["accepted_light_to_heavy_weight_range"]) and nu["run3_working_light_response"]<=common["Run3_weight"]<=1
        assert abs(common["control_shift"])<=cl["extra_decay_phase_cap_rad"]
        ms,mb=m[1:3];e=2*ms*mb/M**2;light=(1-(ms/M)**2)*(1-(mb/M)**2)
        assert abs(sum(la)-e*F/2)<1e-13
        am=at/(1+at*7/(2*math.pi)*math.log(M/mt));qr=(at/az)**(6/23)*(am/at)**(2/7)
        hard=qr*(2*M/mc["vev_GeV"]**2*light)**2/(64*math.pi**2)
        axis=lc.conjugate()**2/abs(lc)**2
        a=pre*axis*(norm*lt**2*common["S0"]*common["etaB"]+
                    p["etaZ"]*(gt/math.sqrt(2)*(1+e**2)+hard)*F**2)
        b=-pre*axis*norm*lt*F*xt*(math.log(M/mw)+h1)*light*(1.18*p["U4"]+1.42*p["U5"])
        w=cl["omega_ps_inverse"][0]; ab=(a.conjugate()*b).real
        fhat2=(-ab+math.sqrt(ab**2+abs(a)**2*(w**2-abs(b)**2)))/abs(a)**2
        y=math.sqrt(fhat2); bhat=fhat2/p["fBs"]**2; Z=a*fhat2+b
        assert inside(y,width(lat["HPQCD_inputs"]["fBs_sqrtB_GeV"])) and inside(bhat,width(weak["bag"]["s"]))
        assert abs(abs(Z)-w)<1e-10
        effective=cmath.phase(Z)+common["control_shift"]
        for tag in ("CMS_KK_13TeV_final","LHCb_KK_2015_18"):
            assert inside(effective,width(cl["phase_inputs_rad"][tag]))
        L=-c*lt; N=L+(1+e)*F; P=p["k"]*L-(1-e)*F
        bnu=common["nu_transfer"]*had["Hnu"]*abs(N)**2/(rare["BelleII"]["reference_SD_in_1e-6"]*1e-6)
        assert inside(bnu,width([rare["BelleII"]["mu"],rare["BelleII"]["width"]]))
        scale=1e9*had["Hmu_heavy"]*(p["fBs"]/nu["fBs_GeV"][0])**2*common["tau"]/weak["tau_ps"]["s_heavy"][0]*common["mu_transfer"]
        angle=cmath.phase(P/lt)-cmath.phase(Z/(axis*lt**2))/2
        mu=[scale*abs(P)**2*(1-(1-common[key])*math.sin(angle)**2) for key in ("CMS_weight","Run3_weight")]
        sc=scale*((mbs/mc["higgs_mass_GeV"])**2*abs(F))**2
        cms=rare["CMS"];r3=weak["LHCb_run3_2024_25"]["s_in_1e9"]
        limits=[(cms["branching_in_1e-9"]-2*sum(cms["lower_widths"]),cms["branching_in_1e-9"]+2*sum(cms["upper_widths"])),
                (r3[0]-2*(r3[1]+r3[3]),r3[0]+2*(r3[2]+r3[4]))]
        assert all(inside(x,l) and inside(x+sc,l) for x,l in zip(mu,limits))
        omega=(gamma_s+math.degrees(cmath.phase(Z)))%180
        combined=read("../1641/inputs.json")["Bs_Run1plus2_Omega_deg"]
        c0,clo,cu=combined
        ledger[name]={"M_GeV":M,"physical_spectrum":m.tolist(),"largest_yukawa":float(yuk),
            "signed_J":float((K[0,0]*K[1,1]*K[0,1].conjugate()*K[1,0].conjugate()).imag),
            "Kud_us_ub_cb_ratio":[ud,us,ub,cb,ub/cb],"Kus_times_fplus":us*common["fplus"],
            "abs_lambda_t":abs(lt),"top_to_charm_phase":cmath.phase((lt/lc)**2),
            "gamma_s_d_mod180":[gamma_s,gamma_d],"f_sqrt_B":y,"Bhat":bhat,"mix_frequency":abs(Z),
            "mix_phase":cmath.phase(Z),"effective_phase":effective,"Belle":bnu,"CMS_Run3_e9":mu,
            "scalar_max_e9":sc,"Omega_deg_mod180":omega,
            "Omega_passes_1641_one_two_width":[inside(omega,(c0-i*clo,c0+i*cu)) for i in (1,2)]}
    return ledger

def main():
    cfg=read("inputs.json"); t=tree(cfg)
    return {"tree":t,"rare_mixing":covers(cfg,t),"physical_witnesses":witnesses(cfg)}
if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--check",action="store_true"); parser.add_argument("--write-new",action="store_true")
    args=parser.parse_args(); val=json.loads(json.dumps(main()),parse_float=lambda v:float(f"{float(v):.10g}"))
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as fd: json.dump(val,fd,indent=2); fd.write("\n")
    elif args.check:
        assert val==read("results.json")
    else:
        print(json.dumps(val,indent=2))
    print("1643 PASS: nuclear first row, physical weak phases and full-bin enclosure")
