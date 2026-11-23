"""Independent Bs clocks for the two fixed mass matrices of round 1635."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent

def load(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))

def enc(z):
    return [z.real, z.imag]

def interval(row, n):
    return [row[0] - n*sum(row[1:]), row[0] + n*sum(row[1:])]

def distance(a, b):
    return max(0, a[0]-b[1], b[0]-a[1])

def alpha_step(alpha, start, end, nf):
    return alpha / (1+alpha*(11-2*nf/3)/(2*math.pi)*math.log(end/start))

def solve():
    cfg=load("inputs.json")
    weak,lat,clock,rare,states=[load(cfg[key]) for key in
        ("weak_inputs","lattice_inputs","clock_inputs","rare_response","weak_states")]
    d=weak["const_GeV"]
    gf,mw,mbs=d["GF"],d["MW"],d["M_Bs"]
    pre=2*mbs/(3*lat["hbar_GeV_ps"])
    norm=gf**2*mw**2/(4*math.pi**2)
    eta,s0,f=weak["short_distance"]["etaB"],weak["short_distance"]["S0_published"],lat["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    a4,a5=cfg["LR_B4_B5_weights"]
    b4,b5=cfg["LR_F2B4_B5_GeV2"]
    al=clock["alpha_s_ref"]; ez=al**(6/23)
    mt,mz=cfg["mt_GeV"],cfg["MZ_GeV"]
    at=alpha_step(al,mz,mt,5); xt=weak["short_distance"]["xt"][0]
    h1=-(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)/(2*(xt-1)**2)*math.log(xt)
    res={}
    for key in ("common_witness","high_mass_witness"):
        b=states[key]; M=b["M_Y_GeV"]
        ms,mb=b["physical_masses_GeV"][1:3]
        K=[[complex(*x) for x in row] for row in b["light_charged_K"]]
        lt=K[2][1].conjugate()*K[2][2]; lc=K[1][1].conjugate()*K[1][2]
        fl,F,_,_,P=[complex(*z) for z in b["F_L_R_and_L_N_P"]]
        r=list(map(lambda z:complex(*z),b["right_D"]))
        assert abs(F+r[1]*r[2].conjugate())<1e-12
        assert abs(fl-2*ms*mb/M**2*F)<1e-15
        light=(1-ms**2/M**2)*(1-mb**2/M**2)
        yprod=-2*(M/cfg["vev_GeV"])**2*light*F
        hard=cfg["RR_hard_Qd_color_weak_factor"]*yprod**2/((4*math.pi)**2*M**2)
        # Leading-log QCD ratio for the same VRR operator across m_t.
        am=alpha_step(at,mt,M,6)
        qcd=(at/al)**(6/23)*(am/at)**(6/21)
        delta_lr=F/lt*xt*(math.log(M/mw)+h1)*light
        # Eq. (6.6) already contains the RG transport of LR1 and the LR2 admixture.
        top=norm*lt**2; z=gf/math.sqrt(2)*(F**2+fl**2)
        loop_lr=-norm*lt**2*delta_lr
        assert (top.conjugate()*z).real>0 and (z.conjugate()*loop_lr).real>0
        axis=lc.conjugate()**2/abs(lc)**2
        tiers={}
        for tier,rc,ec in (("top_Z",0,0),("plus_Qd_hard_RR",1,0),("plus_top_LR",1,1)):
            def amp(ff,tt,zz,wlr):
                return pre*axis*(ff**2*(top*tt+zz*(z+rc*qcd*hard))+ec*loop_lr*wlr)
            central=amp(f[0],eta[0]*s0[0],ez,a4*b4[0]+a5*b5[0])
            phase_rel=cmath.phase(central/(axis*top))
            angle=cmath.phase(P/lt)-phase_rel/2
            scalar_max=(mbs/cfg["higgs_mass_GeV"])**2*abs(F)
            rates={}
            for lab,w in (("CMS",rare["calculated_light_response"]),("LHCb_Run3",rare["run3_central_light_response"])):
                scale=1e9*rare["Hmu_heavy"]
                val=scale*abs(P)**2*(1-(1-w)*math.sin(angle)**2)
                rates[lab]=[val,val+scale*scalar_max**2]
            boxes={}
            for n in (1,2):
                widths=[interval(f,n),list(map(math.prod,zip(interval(eta,n),interval(s0,n)))),
                    clock["eta_Z_RGI_range"],[a4*x+a5*y for x,y in zip(interval(b4,n),interval(b5,n))]]
                pts=[amp(*v) for v in itertools.product(*widths)]
                # All coefficient directions lie in an acute sector. Magnitude is monotone.
                # Argument extrema of a positive multilinear combination are at vertices.
                assert all((x.conjugate()*y).real>0 for x,y in itertools.combinations([top,z,hard,loop_lr],2))
                freq=[min(map(abs,pts)),max(map(abs,pts))]
                ph=[min(map(cmath.phase,pts)),max(map(cmath.phase,pts))]
                data={lab:interval(v,n) for lab,v in clock["phase_inputs_rad"].items()}
                p1,p2=data["LHCb_KK_2015_18"],data["CMS_KK_13TeV_final"]
                phi=[max(p1[0],p2[0]),min(p1[1],p2[1])]
                ob=interval(clock["omega_ps_inverse"],n)
                gap=distance(ph,phi)
                ang=max(0,gap-clock["extra_decay_phase_cap_rad"])
                # Lower bound from a superset allowing angular and radial extrema independently.
                u=min(max(math.cos(ang),ob[0]/freq[1]),ob[1]/freq[0])
                joint=math.hypot(1-u*math.cos(ang),u*math.sin(ang))
                boxes[str(n)]={"frequency":freq,"phase":ph,"common_KK":phi,
                     "freq_gap":distance(freq,ob),"decay_phase_min":gap,
                     "relative_mixing_residual_lower_bound":joint,
                     "phase_gap_per_sample":{k:distance(ph,v) for k,v in data.items()}}
            tiers[tier]={"frequency":abs(central),"phase":cmath.phase(central),
                         "updated_muon_br_e9":rates,"boxes":boxes}
        zrr=gf/math.sqrt(2)*F**2
        res[key]={"M_GeV":M,"lambda_t":enc(lt),"lambda_c":enc(lc),"F":enc(F),
            "Qd_hard_to_tree_RR_high_scale_abs":abs(hard/zrr),"VRR_high_to_Z_RG_ratio":qcd,
            "Qd_hard_to_tree_RR_after_RG_abs":abs(qcd*hard/zrr),
            "top_induced_delta_C_LR1_at_mt":enc(delta_lr),
            "tree_LR_to_tree_RR_bare":abs(2*fl/F),
            "scalar_sum_coeffs_to_tree_RR_bare":(ms+mb)**2/(2*cfg["vev_GeV"]**2*cfg["higgs_mass_GeV"]**2)/(gf/math.sqrt(2)),
            "leading_LR_to_top_amplitude_abs":abs(loop_lr*(a4*b4[0]+a5*b5[0])/(f[0]**2*top*eta[0]*s0[0])),
            "tiers":tiers}
    return res

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--write-new",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=json.loads(json.dumps(solve()),parse_float=lambda x:float(f"{float(x):.12g}"))
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert result==load("results.json")
    print("1636 PASS: two fixed matrices, mixed clocks, complex residuals and dimuon feedback")
