"""Common neutral-scalar weak phases and two neutral-meson clocks."""
from pathlib import Path
from itertools import product
import argparse
import cmath
import json
import math

HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))
def box(row,n=2):
    err=n*sum(row[1:])
    return row[0]-err,row[0]+err
def cosine(s):
    return math.sqrt(1-s*s)
def abs_ckm(us,cb,r):
    s13=(cb[0]*r[0],cb[1]*r[1])
    c13=(cosine(s13[1]),cosine(s13[0]))
    s12=(us[0]/c13[1],us[1]/c13[0])
    s23=(cb[0]/c13[1],cb[1]/c13[0])
    c12=(cosine(s12[1]),cosine(s12[0]))
    c23=(cosine(s23[1]),cosine(s23[0]))
    def times(*args):
        return math.prod(a[0] for a in args),math.prod(a[1] for a in args)
    def sides(a,b):
        return max(0,a[0]-b[1],b[0]-a[1]),a[1]+b[1]
    V=[[times(c12,c13),us,s13],
       [sides(times(s12,c23),times(c12,s23,s13)),sides(times(c12,c23),times(s12,s23,s13)),cb],
       [sides(times(s12,s23),times(c12,c23,s13)),sides(times(c12,s23),times(s12,c23,s13)),times(c23,c13)]]
    return V,(s12,c12,s23,c23,s13,c13)
def heldout_phase_floor(trig,gamma_window,ratio):
    s12,c12,s23,c23,s13,c13=trig
    kmax=c12[1]*s23[1]*s13[1]/(s12[0]*c23[0])
    sinmin=math.sin(math.radians(gamma_window[0])-math.asin(kmax))
    assert math.radians(gamma_window[1])+math.asin(kmax)<math.pi/2
    amin=s12[0]*c23[0]*c13[0]*ratio[0]/c12[1]
    amax=s12[1]*c23[1]*c13[1]*ratio[1]/c12[0]
    bmax=s12[1]*s23[1]*s13[1]/(c12[0]*c23[0])
    return 2*math.atan(amin*sinmin/(1+amax)),2*(math.asin(amax)+math.asin(bmax))
def ckm(us,cb,r,gam):
    s13=cb*r
    c13=cosine(s13)
    s12,s23=us/c13,cb/c13
    c12,c23=cosine(s12),cosine(s23)
    delta=gam+math.asin(c12*s23*s13/(s12*c23)*math.sin(gam))
    exp=cmath.exp(1j*delta)
    return [[c12*c13,s12*c13,s13*exp.conjugate()],
            [-s12*c23-c12*s23*s13*exp,c12*c23-s12*s23*s13*exp,cb],
            [s12*s23-c12*c23*s13*exp,-c12*s23-s12*c23*s13*exp,c23*c13]]
def central_examples(p,tree,belle,ys,gam):
    out=[]
    for branch in [0,180]:
        V=ckm(tree["K_product"][0]/tree["K_form"][0],belle[0],tree["ub_over_cb"][0],math.radians(gam+branch))
        shapes={}
        for label,i,j in [("bd",0,2),("bs",1,2),("ds",0,1)]:
            lam=[row[i].conjugate()*row[j] for row in V]
            absolute=[]
            phase=[]
            for su,sc in product([-1,1],repeat=2):
                X=sum(sign*y*l for sign,y,l in zip((su,sc,1),ys,lam))
                absolute.append(abs(X))
                phase.append(math.remainder(2*cmath.phase(X/lam[2]),math.pi))
            shapes[label]={"abs_X_extrema":[min(absolute),max(absolute)],"phase_line_radian":[min(phase),max(phase)]}
        lam=[row[1].conjugate()*row[2] for row in V]
        shapes["Bs_top_phase_rad"]=math.remainder(2*cmath.phase(lam[2]/lam[1]),2*math.pi)
        out.append(shapes)
    return out
def compatible_point(p,tree,belle,ys,gam,top):
    """Use both clocks to infer two bag inputs only for an existence test."""
    V=ckm(tree["K_product"][0]/tree["K_form"][0],belle[0],tree["ub_over_cb"][0],math.radians(gam))
    const=top["const_GeV"]
    mass=100000.0
    eta4=(p["alpha_mb"][0]/p["alpha_MZ"][0])**(24/23)
    hbar=p["hbar_GeV_ps"]
    point={"M_abs_cos2beta_GeV":mass}
    for q,i in [("d",0),("s",1)]:
        m=const["M_B"+q]
        mb=sum(p["hadron_quark_mass_bounds_GeV"]["b"])/2
        mq=sum(p["hadron_quark_mass_bounds_GeV"][q])/2
        q4=.5*m*((m/(mb+mq))**2+1/6)*p["Q4_f2B_GeV2"][q][0]
        lam=[row[i].conjugate()*row[2] for row in V]
        X=sum(a*b for a,b in zip(ys,lam))
        scalar=X*X*eta4*q4/(mass*mass*hbar)
        top_per_bag=(const["GF"]**2*const["MW"]**2/(6*math.pi**2)*m*
                     top["short_distance"]["etaB"][0]*top["short_distance"]["S0_published"][0]*lam[2]**2/hbar)
        scale=p["clock"][q][0]/abs(top_per_bag)
        constraint=scalar/top_per_bag
        fbag=math.sqrt(scale*scale-constraint.imag**2)-constraint.real
        amplitude=top_per_bag*fbag+scalar
        assert box(p["Q1_f2Bhat_GeV2"][q])[0]<=fbag<=box(p["Q1_f2Bhat_GeV2"][q])[1]
        assert abs(abs(amplitude)-p["clock"][q][0])<1e-12
        point[q]={"f2Bhat_inferred_GeV2":fbag,"frequency_per_ps":abs(amplitude),
                  "phase_rad_vs_charm_squared":cmath.phase(amplitude/lam[1]**2)}
    assert box(p["phi_s"])[0]<=point["s"]["phase_rad_vs_charm_squared"]<=box(p["phi_s"])[1]
    return point


def solve():
    p=read(HERE/"inputs.json")
    tree=read(HERE/p["semileptonic"])
    belle=read(HERE/p["Belle_BGL"])["BelleII_Vcb_BGL"]
    gamma=read(HERE/p["gamma"])["chargedB"]
    yrow=read(HERE/p["ordinary_yukawa_MZ"])["yukawa"]["up"]
    y=[box(item) for item in yrow]
    top=read(HERE/p["top_box"])
    const=top["const_GeV"]
    w=p["transfer_examples"]
    kp,kf=box(tree["K_product"]),box(tree["K_form"])
    V,trig=abs_ckm((kp[0]/kf[1],kp[1]/kf[0]),box(belle),box(tree["ub_over_cb"]))
    kr=heldout_phase_floor(trig,[gamma["central_degrees"]-2*gamma["lower_error"]-gamma["additional_degrees"],
                                 gamma["central_degrees"]+2*gamma["upper_error"]+gamma["additional_degrees"]],box(tree["ub_over_cb"]))
    out={}
    ang={}
    for label,i,j,leader in [("bs",1,2,2),("bd",0,2,2),("ds",0,1,1)]:
        lam=[(row[i][0]*row[j][0],row[i][1]*row[j][1]) for row in V]
        Xlo=y[leader][0]*lam[leader][0]-sum(y[k][1]*lam[k][1] for k in range(3) if k!=leader)
        Xhi=sum(a[1]*b[1] for a,b in zip(y,lam))
        assert Xlo>0
        out[label]={"X_bounds_any_sign_any_delta":[Xlo,Xhi]}
        if label=="ds":
            continue
        # Remove common charm amplitude before bounding the phase.
        eps=(y[0][1]+y[1][1])/(y[2][0]-y[1][1])*lam[0][1]/lam[2][0]
        theta=2*math.asin(eps)
        q=label[-1]
        m=const["M_B"+q]
        had=box(p["Q4_f2B_GeV2"][q])[0]
        mb,mq=p["hadron_quark_mass_bounds_GeV"]["b"][1],p["hadron_quark_mass_bounds_GeV"][q][1]
        Q4overM=0.5*m*((m/(mb+mq))**2+1/6)*had
        eta4=(box(p["alpha_mb"])[0]/box(p["alpha_MZ"])[1])**(24/23)
        om=box(p["clock"][q])
        clock_hi=om[1]*(1+w["abs_clock_2M_relative"])
        clock_lo=om[0]*(1-w["abs_clock_2M_relative"])
        sm_hi=(const["GF"]**2*const["MW"]**2/(6*math.pi**2)*m*
               box(top["short_distance"]["etaB"])[1]*box(top["short_distance"]["S0_published"])[1]*
               box(p["Q1_f2Bhat_GeV2"][q])[1]*lam[2][1]**2/p["hbar_GeV_ps"])
        other=w["phase_abs_background_over_measured"]
        a=theta+math.asin(sm_hi/clock_lo*math.sin(theta)+other)
        mass=math.sqrt(w["absolute_scalar_rate_lower"]*Xlo**2*eta4*Q4overM/
                       (p["hbar_GeV_ps"]*((1+other)*clock_hi+sm_hi)))
        out[label].update({"top_relative_scalar_line_angle_rad_upper":theta,
                           "other_permitted_abs_over_measured":other,"top_clock_upper_per_ps":sm_hi,
                           "total_phase_mod_pi_from_top_rad_upper":a,"Q4_over_MB_GeV3_lower":Q4overM,
                           "leading_QCD_transport_lower":eta4,"degenerate_M_abs_cos2beta_GeV_lower":mass})
        if q=="s":
            phase_box=box(p["phi_s"])
            ang={"top_phase_abs_box_rad":kr,"experimental_two_width":phase_box,
                 "decay_phase_tolerance_for_wrong_branch_separation_rad":kr[0]-a-phase_box[1],
                 "wrong_branch_lower_with_budget_rad":kr[0]-a-w["decay_and_extraction_phase_rad"],
                 "right_branch_upper_with_budget_rad":-kr[0]+a+w["decay_and_extraction_phase_rad"]}
    assert ang["wrong_branch_lower_with_budget_rad"]>ang["experimental_two_width"][1]
    witness=compatible_point(p,tree,belle,[row[0] for row in yrow],gamma["central_degrees"],top)
    return {"finite_compatible_example":witness,"all_delta_and_sign_boxes":out,"signed_Bs_phase_test":ang,
            "central_tree_two_gamma_branches":central_examples(p,tree,belle,[row[0] for row in yrow],gamma["central_degrees"])}
if __name__=="__main__":
    # Recompute without changing any stored reference in check mode.
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    target=HERE/"results.json"
    now=solve()
    if args.check:
        # Arrays already serialized as JSON in the stored reference.
        assert json.loads(json.dumps(now))==json.loads(target.read_text(encoding="utf-8-sig"))
        print("1581 checked")
    else:
        target.write_text(json.dumps(now,indent=2)+"\n",encoding="utf-8")
