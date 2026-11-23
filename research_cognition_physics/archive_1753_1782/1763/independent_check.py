"""Independent Bd-phase closure and modulus-ratio checks on the archived CKM matrices."""
import argparse
import cmath
import json
import math
from pathlib import Path


HERE=Path(__file__).resolve().parent
ALPHA=[(92.6,4.7,4.5),(177.4,4.5,4.7)]
DALITZ=(.91,.22,.09,.07)


def read(path):
    return json.loads((HERE/path).read_bytes())


def interval_distance(angle, segments):
    return min(max(left-angle-turn,angle+turn-right,0.)
               for left,right in segments for turn in [-360,-180,0,180,360])


def degree(z):
    return math.degrees(cmath.phase(z))%180


def sd(beta):
    return math.sin(math.radians(beta))


def phase_boxes(data):
    s,ss= data["S_d"][0],math.hypot(*data["S_d"][1:])
    c,sc= data["C_d"][0],math.hypot(*data["C_d"][1:])
    cmin=max(0,abs(c)-sc)
    cmax=abs(c)+sc
    def from_sc(s,c):
        return .5*math.degrees(math.asin(s/math.sqrt(1-c*c)))
    low,high=from_sc(s-ss,cmin),from_sc(s+ss,cmax)
    betas=[[low,high],[90-high,90-low]]
    alphas=[[x-a,x+b] for x,a,b in ALPHA]
    boxes=[]
    for bl,bh in betas:
        boxes.append([[180-ah-bh,180-al-bl] for al,ah in alphas])
    return betas,alphas,boxes


def sample(row,data,alphas,boxes):
    if "low_moduli" in row:
        us,ub,cb=[row["low_moduli"]["V"+k] for k in ("us","ub","cb")]
    else:
        us,ub,cb=[row[k] for k in ("us","ub","cb")]
    gamma=row["gamma_deg"]
    c13=math.sqrt(1-ub*ub)
    s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12*s12),math.sqrt(1-s23*s23)
    g=math.radians(gamma)
    delta=g+math.asin(c12*s23*ub*math.sin(g)/(s12*c23))
    e=cmath.exp(1j*delta)
    vud,vub=c12*c13,ub/e
    vcd,vcb=-s12*c23-c12*s23*ub*e,cb
    vtd,vtb=s12*s23-c12*c23*ub*e,c23*c13
    up,charm,top=vud*vub.conjugate(),vcd*vcb.conjugate(),vtd*vtb.conjugate()
    assert abs(up+charm+top)<1e-16
    alpha,beta,gcalc=degree(-top/up),degree(-charm/top),degree(-up/charm)
    assert abs(gcalc-gamma)<1e-10
    assert abs(alpha+beta+gamma-180)<1e-10
    observable=abs(up/top)**2
    assert abs(observable-(sd(beta)/sd(gamma))**2)<1e-13
    signed=boxes[0]
    all_boxes=boxes[0]+boxes[1]
    opposite=180-gamma
    # Rotate a common mixing phase to the obtuse golden-phase branch. Its
    # acutely and obtusely reflected partners retain the individual sin phase.
    c=data["C_d"][0]
    be=.5*math.degrees(math.asin(data["S_d"][0]/math.sqrt(1-c*c)))
    mirror_beta_effective=90-be
    phi=mirror_beta_effective+beta
    alpha_effective=(-alpha-phi)%180
    assert abs(sd(2*mirror_beta_effective)-sd(2*be))<2e-15
    return {"us_ub_cb":[us,ub,cb],"alpha_beta_gamma_deg":[alpha,beta,gamma],
            "rho_isospin_distance_if_standard_mixing_deg":interval_distance(alpha,alphas),
            "ratio_VudVub_over_VtdVtb_squared":observable,
            "main_Btau_clock_CKM_factor":observable/vud**2,
            "arbitrary_common_mixing":dict(
                gamma_residual_with_cos_sign_deg=interval_distance(gamma,signed),
                conjugate_residual_without_cos_sign_deg=interval_distance(opposite,all_boxes),
                conjugate_residual_with_cos_sign_deg=interval_distance(opposite,signed)),
            "conjugate_example_at_golden_central":dict(mixing_extra_half_phase_deg_mod180=phi%180,
                rho_effective_alpha_deg=alpha_effective, rho_range_distance_deg=interval_distance(alpha_effective,alphas),
                effective_cos_two_beta=math.cos(math.radians(2*mirror_beta_effective)))}


def main():
    ka=read("../1707/inputs.json")
    betas,alphas,boxes=phase_boxes(ka)
    # Normal SM mixing, measured acute-gamma branch, and oriented positive tree:
    # beta=180-alpha-gamma.  The near-180° alpha band violates beta>0 here.
    g=ka["time_integrated_gamma_deg"]
    gl,gh=g["68"]
    al,ah=alphas[0]
    def predicted_ratio(alpha,gamma):
        return (sd(alpha+gamma)/sd(gamma))**2
    assert 90<al+gl<ah+gh<180 and 0<gl<gh<90
    sinerule=[predicted_ratio(ah,gh),predicted_ratio(al,gl)]
    onestar=read("../1755/results.json")["full_two_threshold_shooting"]
    two=read("../1757/results.json")
    rows={f"1755_{k}_{j}":sample(r,ka,alphas,boxes) for k,v in onestar.items() for j,r in v.items()}
    rows.update({f"1757_{k}":sample(two[k],ka,alphas,boxes) for k in ("central","only_reactor_at_upper_one_width")})
    rows.update({f"1757_vertex_{k}":sample(r,ka,alphas,boxes) for k,r in two["eight_one_width_vertices"].items()})
    cos1_low=DALITZ[0]-math.hypot(*DALITZ[1:])
    cos1_half=.5*math.degrees(math.acos(cos1_low))
    return {"source":{"Belle_II_rho_isospin":"https://arxiv.org/html/2412.19624v2#S7",
                      "golden_psi_KS":"https://arxiv.org/abs/2309.09728",
                      "cos_sign":"https://arxiv.org/abs/1804.06153",
                      "isospin_corrections":"https://arxiv.org/abs/hep-ph/0502139",
                      "time_integrated_gamma":"../1707/inputs.json"},
            "rho_alpha_branches_deg":alphas,"beta_eff_one_marginal_box_deg":betas,
            "gamma_boxes_modulo_180_deg_small_and_large_beta":boxes,
            "cos2beta_BaBar_Belle":DALITZ,"experimental_cos_positive_evidence_sigma":3.7,
            "extra_Dalitz_minus_golden_half_phase_needed_for_obtuse_branch_deg":{
                "sign_only":betas[1][0]-45., "cosine_at_one_marginal_lower":betas[1][0]-cos1_half},
            "sine_rule_phase_only_input_alpha_and_direct_gamma":dict(
                D_K_68_window=g["68"],ratio_central=predicted_ratio(ALPHA[0][0],g["central"]),
                one_reported_range_box=sinerule),
            "model_results":rows,
            "all_18_frozen_cases":dict(total=len(rows),
                standard_alpha_pass=sum(r["rho_isospin_distance_if_standard_mixing_deg"]==0 for r in rows.values()),
                common_mixing_with_cos_sign_pass=sum(r["arbitrary_common_mixing"]["gamma_residual_with_cos_sign_deg"]==0 for r in rows.values()),
                mirror_without_cos_sign_pass=sum(r["arbitrary_common_mixing"]["conjugate_residual_without_cos_sign_deg"]==0 for r in rows.values()),
                mirror_with_cos_sign_min_needed_decay_shift_deg=min(r["arbitrary_common_mixing"]["conjugate_residual_with_cos_sign_deg"] for r in rows.values()),
                ratio_range=[min(r["ratio_VudVub_over_VtdVtb_squared"] for r in rows.values()),max(r["ratio_VudVub_over_VtdVtb_squared"] for r in rows.values())])}


def tidy(value):
    if isinstance(value,dict):
        return {k:tidy(x) for k,x in value.items()}
    if isinstance(value,(list,tuple)):
        return [tidy(x) for x in value]
    return float(f"{value:.11g}") if isinstance(value,float) else value


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=tidy(main())
    path=HERE/"independent_results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1763 independent OK: rho isospin, common-mixing cancellation, cosine sign and modulus-ratio check")
    else:
        path.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf-8")
