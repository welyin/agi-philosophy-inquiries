"""Independent 1217 checks: original populations, exact spectral tails,
rational logarithm enclosures, binary bounds, and deletion witnesses.
No author or earlier check program is imported. All signed decisions are rational.
Default is read-only replay; --write creates the result only once.
"""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache
import argparse, hashlib, json

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=STAGE.parent
LOG_GRID=10**20
LOG_TERMS=24

def require(c,label):
    if not c: raise AssertionError(label)

def floor_int(x): return x.numerator//x.denominator
def ceil_int(x): return -((-x.numerator)//x.denominator)
def enclose(lo,hi,grid=LOG_GRID):
    require(lo<=hi,"ordered enclosure")
    return F(floor_int(lo*grid),grid),F(ceil_int(hi*grid),grid)
def enc(x): return [str(x[0]),str(x[1])]

def series_log_on_unit_interval(y):
    require(1<=y<=2,"normalized logarithm argument")
    z=(y-1)/(y+1)
    s=F(0)
    power=z
    for k in range(LOG_TERMS):
        s+=2*power/F(2*k+1)
        power*=z*z
    remainder=2*power/(F(2*LOG_TERMS+1)*(1-z*z))
    return s,s+remainder

@lru_cache(maxsize=None)
def log_bounds(x):
    x=F(x)
    require(x>0,"positive logarithm")
    if x==1:return F(0),F(0)
    k=0;y=x
    while y<1:y*=2;k-=1
    while y>=2:y/=2;k+=1
    lo,hi=series_log_on_unit_interval(y)
    l2,u2=series_log_on_unit_interval(F(2))
    if k>=0:lo,hi=lo+k*l2,hi+k*u2
    else:lo,hi=lo+k*u2,hi+k*l2
    return enclose(lo,hi)

def V_interval(pop):
    require(sum(pop)==1 and min(pop)>0,"faithful probability vector")
    n=len(pop)
    logs=[log_bounds(x) for x in pop]
    ln=log_bounds(F(n))
    return (-ln[1]-sum(v[1] for v in logs)/n,
            -ln[0]-sum(v[0] for v in logs)/n)

def binary_interval(q,p):
    require(0<q<1 and 0<p<1,"binary interior")
    x=log_bounds(q/p);y=log_bounds((1-q)/(1-p))
    return q*x[0]+(1-q)*y[0],q*x[1]+(1-q)*y[1]

def low_root_bracket(a,V,grid=10**10):
    require(0<a<1 and 0<V[0]<=V[1],"root premises")
    left=0;right=floor_int(a*grid)
    require(binary_interval(a,F(right,grid))[1]<V[0],"upper root bracket")
    while right-left>1:
        middle=(left+right)//2
        db=binary_interval(a,F(middle,grid))
        if db[0]>V[1]:left=middle
        elif db[1]<V[0]:right=middle
        else:break
    require(left>0,"positive lower root")
    lo,hi=F(left,grid),F(right,grid)
    require(binary_interval(a,lo)[0]>V[1],"lower root sign")
    require(binary_interval(a,hi)[1]<V[0],"upper root sign")
    return lo,hi

def tail(pop,alpha):
    require(0<=alpha<=len(pop),"tail domain")
    s=sorted(pop)
    m=floor_int(alpha)
    return sum(s[:m])+(alpha-m)*s[m] if m<len(s) else sum(s)

def tail_effect(pop,alpha):
    order=sorted(range(len(pop)),key=lambda i:pop[i])
    e=[F(0)]*len(pop)
    remaining=F(alpha)
    for j in order:
        take=min(F(1),remaining)
        e[j]=take;remaining-=take
    require(remaining==0,"effect trace")
    return e

def sums_marginals(p,d1,d2):
    return ([sum(p[i*d2+j] for j in range(d2)) for i in range(d1)],
            [sum(p[i*d2+j] for i in range(d1)) for j in range(d2)])

def compute(bind_author=True):
    cold=json.loads((STAGE/"1202/results.json").read_text(encoding="utf-8-sig"))
    hot=json.loads((STAGE/"1204/results.json").read_text(encoding="utf-8-sig"))
    clock=json.loads((STAGE/"1203/results.json").read_text(encoding="utf-8-sig"))
    pc=[F(x) for x in cold["one_marginal_rational"]]
    ph=[F(x) for x in hot["new_reference"]["one_marginal_rational"]]
    require(sum(pc)==sum(ph)==1 and min(pc+ph)>0,"native exact sources")
    joint=[x*y for x in pc for y in ph]
    ordered=sorted(enumerate(joint),key=lambda t:t[1])
    base_tail=tail(joint,F(12))
    require(F("0.540139055")<base_tail<F("0.540139057"),"native zero-error tail")
    g1={"name":"native_exact_joint_spectrum","passed":True,
        "dimension":16,"bad_projection_rank":12,
        "cold_population":[str(x) for x in pc],"hot_population":[str(x) for x in ph],
        "joint_population":[str(x) for x in joint],
        "ascending_original_indices":[i for i,v in ordered],
        "ideal_rank_12_tail":str(base_tail),
        "minimum_population":str(min(joint))}

    # Independent atanh implementation: 24 terms, outward 10^-20 grid.
    ln2=log_bounds(F(2))
    require(F("0.693147180559945309")<ln2[0]<=ln2[1]<F("0.693147180559945310"),
            "log2 rational enclosure")
    require(log_bounds(F(1))==(F(0),F(0)),"log1")
    for x in (F(1,3),F(7,5),F(8),F(1,64)):
        lx=log_bounds(x);li=log_bounds(1/x)
        require(lx[0]+li[0]<=0<=lx[1]+li[1],"reciprocal log overlap")
    vc,vh,vj=V_interval(pc),V_interval(ph),V_interval(joint)
    additive=(vc[0]+vh[0],vc[1]+vh[1])
    require(max(vj[0],additive[0])<=min(vj[1],additive[1]),"joint/product V agreement")
    require(F("0.20400413")<vj[0]<=vj[1]<F("0.20400414"),"native joint V range")
    # Use the wider direct joint interval for all event decisions.
    g2={"name":"independent_rational_log_and_native_V0","passed":True,
        "log_terms":LOG_TERMS,"log_outward_grid":LOG_GRID,
        "log2_interval":enc(ln2),"cold_V0_interval":enc(vc),
        "hot_V0_interval":enc(vh),"direct_joint_V0_interval":enc(vj),
        "sum_of_marginal_V0_intervals":enc(additive),
        "floating_values_used_in_strict_certificates":False}

    event_results=[]
    requested={1:("0.539998258","0.43543064","0.30228764"),
               2:("0.540062252","0.43548341","0.30231162")}
    for v in (1,2):
        eps=F(clock["bounds_by_program"][str(v)]["isometry_error_upper"])
        expected=F("0.000109379990") if v==1 else F("0.000059665807")
        require(eps==expected,"unaltered clock+shape error")
        a=F(3,4)-eps
        alpha=16*a
        spectral=tail(joint,alpha)
        old_linear=min(joint)*alpha
        t_claim,b_claim,e_claim=map(F,requested[v])
        require(spectral>t_claim and spectral>old_linear,"native improved spectral tail")
        btest=binary_interval(a,b_claim)
        require(b_claim<a and btest[0]>vj[1],"binary lower-bound certificate")
        root=low_root_bracket(a,vj)
        require(root[0]>b_claim,"independent binary root bracket")
        log_claim=log_bounds(e_claim)
        exponential_gap=-a*log_claim[1]-vj[1]-ln2[1]
        require(exponential_gap>0,"explicit exponential lower bound")
        effect=tail_effect(joint,alpha)
        require(sum(effect)==alpha and all(0<=x<=1 for x in effect),"native optimal effect")
        require(sum(p*x for p,x in zip(joint,effect))==spectral,"native tail saturation")
        event_results.append({"program_v":v,"rectangular_contact_duration":600//v,
            "separately_covered_events":["actual_A_not_ground","actual_K_not_zero"],
            "same_event_or_shotwise_agreement_claimed":False,
            "epsilon_half_diamond_or_isometry":str(eps),
            "q_lower":str(a),"effect_trace_lower":str(alpha),
            "exact_fractional_spectral_tail":str(spectral),
            "strict_spectral_lower":str(t_claim),
            "old_lambda_linear_bound":str(old_linear),
            "binary_root_interval":enc(root),"strict_binary_lower":str(b_claim),
            "binary_sign_margin":str(btest[0]-vj[1]),
            "strict_explicit_lower":str(e_claim),
            "explicit_log_sign_margin":str(exponential_gap)})
    g3={"name":"same_native_actual_A_and_K_event_certificates","passed":True,
        "events":event_results,
        "initial_AB_independent_of_clock_and_blank_K":True,
        "passive_source_reference_retained":True,
        "real_J_tails_in_clock_bound":True,
        "postselection_used":False,
        "actual_1205_post_preparation_channel_claimed":False}

    # A nonuniform rational spectrum, including fractional filling and alpha=0.
    small=[F(1,12),F(2,12),F(3,12),F(6,12)]
    saturations=[]
    d=4;r=2
    for eps in (F(0),F(1,8),F(1,2)):
        alpha=max(F(0),F(r)-d*eps)
        e=tail_effect(small,alpha)
        q=sum(e)/d
        out_uniform=[q/r]*r+[(1-q)/(d-r)]*(d-r)
        defect=sum(abs(x-F(1,d)) for x in out_uniform)/2
        p=sum(x*y for x,y in zip(small,e))
        require(p==tail(small,alpha),"measure-prepare exact optimum")
        require(defect<=eps and sum(out_uniform)==1,"CPTP information-class defect")
        # Every computational input maps to a normalized positive output.
        for ei in e:
            response=[ei/r]*r+[(1-ei)/(d-r)]*(d-r)
            require(sum(response)==1 and min(response)>=0,"preparation channel")
        saturations.append({"epsilon":str(eps),"alpha":str(alpha),
            "effect":[str(x) for x in e],"q":str(q),"p":str(p),
            "normalized_unital_defect":str(defect)})
    # Rotate both a faithful state and its minimizing effect by rational U.
    U=[[F(3,5),F(-4,5),F(0),F(0)],
       [F(4,5),F(3,5),F(0),F(0)],
       [F(0),F(0),F(1),F(0)],[F(0),F(0),F(0),F(1)]]
    eff=tail_effect(small,F(3,2))
    rho=[[sum(U[i][k]*small[k]*U[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    E=[[sum(U[i][k]*eff[k]*U[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    p_rot=sum(rho[i][j]*E[j][i] for i in range(4) for j in range(4))
    require(p_rot==tail(small,F(3,2)),"coherent-basis exact optimum")
    g4={"name":"fractional_spectral_optimum_information_class","passed":True,
        "saturating_measure_prepare_examples":saturations,
        "coherent_basis_saturating_probability":str(p_rot),
        "original_QQ_or_work_role_reachability_claimed":False}

    # Product additivity, a strictly mixing unital map, and correlated deletion.
    mixed=[(pc[0]+pc[1])/2,(pc[0]+pc[1])/2,
           (pc[2]+pc[3])/2,(pc[2]+pc[3])/2]
    vm=V_interval(mixed)
    require(vm[1]<vc[0],"V decreases under this explicit random swap mixture")
    correlated=[]
    previous=F(-1)
    for delta in (F(1,16),F(1,256),F(1,4096)):
        corr=[(1-delta)/2,delta/2,delta/2,(1-delta)/2]
        marg=sums_marginals(corr,2,2)
        require(marg==([F(1,2)]*2,[F(1,2)]*2),"maxmixed correlated marginals")
        vv=V_interval(corr)
        require(vv[0]>previous,"joint V increases along correlated family")
        previous=vv[1]
        require(tail(corr,F(2))==delta,"correlated exact rank-2 tail")
        correlated.append({"delta":str(delta),"joint_V_interval":enc(vv),
                           "both_marginal_V":0,"rank_2_tail":str(delta)})
    g5={"name":"V0_additivity_unital_mixing_and_joint_correlation","passed":True,
        "mixed_cold_V_interval":enc(vm),"initial_cold_V_interval":enc(vc),
        "correlated_examples":correlated}

    # At fixed d and bounded material spectrum, faithful resources approach pure.
    family=[]
    previous=F(-1)
    for eta in (F(1,16),F(1,256),F(1,4096)):
        nu=[1-3*eta,eta,eta,eta]
        vv=V_interval(nu)
        require(vv[0]>previous,"V grows at fixed dimension")
        previous=vv[1]
        active=[F(1,4)*x for i in range(4) for x in nu]
        # Same-h SWAP: target output is nu. Uniform calibration remains 3/4.
        output_target=nu
        pbad=sum(output_target[1:])
        require(pbad==3*eta and pbad>0,"finite faithful SWAP failure")
        vactive=V_interval(active)
        require(max(vv[0],vactive[0])<=min(vv[1],vactive[1]),"maxmixed factor adds zero V")
        entropy_lo=-sum(p*log_bounds(p)[1] for p in nu)
        require(entropy_lo>0,"positive entropy; forward D bounded by ln4")
        family.append({"eta":str(eta),"fixed_active_dimension":16,
            "q":str(F(3,4)),"exact_failure":str(pbad),"V0_interval":enc(vv),
            "forward_relative_entropy_upper":str(log_bounds(F(4))[1]),
            "energy_within_same_fixed_h_spectral_interval":True})
    pure=[F(1),F(0),F(0),F(0)]
    require(sum(pure[1:])==0,"pure stock deletion")
    g6={"name":"fixed_dimension_resource_family_and_quantifier_deletion","passed":True,
        "faithful_finite_examples":family,
        "pure_limit_has_finite_dimension_and_zero_failure":True,
        "pure_limit_V0": "infinite; budget ineligible, not infinite physical cost",
        "Q1_pointwise_no_uniform_resource_growth_required":True,
        "Q2_uniform_bound_not_implied_by_individual_faithfulness":True,
        "Q3_time_budget_is_an_extra_sufficient_interface":True,
        "infinite_active_input_automatically_in_finite_task_domain":False}

    if bind_author:
        author=json.loads((HERE/"results.json").read_text(encoding="utf-8"))
        ag=author["groups"]
        require([F(x) for x in ag["native_objects_and_quantifiers"]["sorted_joint_spectrum"]]
                ==sorted(joint),"author full joint spectrum")
        require(F(ag["actual_A_and_K_spectral_bounds"]["exact_unital_tail"])==base_tail,
                "author ideal spectral tail")
        for key,own in (("V_C",vc),("V_H",vh),("V_CH_joint_check",vj)):
            av=ag["reverse_relative_entropy_budget"][key]
            saved=(F(av["lower"]),F(av["upper"]))
            require(max(own[0],saved[0])<=min(own[1],saved[1]),"author V enclosure "+key)
        ar={row["contact"]:row for row in ag["actual_A_and_K_spectral_bounds"]["contacts"]}
        br={row["contact"]:row for row in ag["binary_and_additive_uniform_bounds"]["native_contacts"]}
        for own in event_results:
            label=str(own["rectangular_contact_duration"])
            require(F(ar[label]["sharp_tail_exact"])==F(own["exact_fractional_spectral_tail"]),
                    "author actual tail "+label)
            require(F(ar[label]["a"])==F(own["q_lower"]),"author actual q "+label)
            saved=br[label]["binary_root_bracket"]
            lo,hi=map(F,own["binary_root_interval"])
            require(F(saved["lower"])<=lo<hi<=F(saved["upper"]),"author binary root enclosure "+label)
            el=F(br[label]["explicit_exponential_rational_lower"])
            require(-F(own["q_lower"])*log_bounds(el)[1]-vj[1]-ln2[1]>0,
                    "author explicit exponential bound "+label)
        g2["author_V_enclosures_independently_crosschecked"]=True
        g3["author_native_exact_tails_and_binary_brackets_crosschecked"]=True

    groups=[g1,g2,g3,g4,g5,g6]
    require(all(g["passed"] for g in groups),"all independent groups")
    sources=["archive_1201_/1202/results.json","archive_1201_/1203/proof.md",
        "archive_1201_/1203/results.json","archive_1201_/1204/proof.md",
        "archive_1201_/1204/results.json","archive_1201_/1214/proof.md",
        "archive_1201_/1215/proof.md","archive_1201_/1216/proof.md"]
    if bind_author:
        sources += ["archive_1201_/1217/"+f for f in ("proof.md","check.py","results.json","NEXT.md")]
        sources += ["archive_1201_/research_note_1217.md"]
    hashes={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}
    return {"round":1217,"status":"INDEPENDENT_SCOPED_PASS_NOT_ACCEPTANCE",
        "all_checks_passed":True,"check_groups":6,"groups":groups,
        "source_path_base":"research_cognition_physics","source_sha256":hashes,
        "author_program_imported":False,"strict_arithmetic":"integer/Fraction only",
        "scope":"same native independent-AB process; A and K separately; no universal resource-time law"}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    target=HERE/"independent_results.json"
    if args.write and target.exists():
        raise SystemExit("Refuse to overwrite independent_results.json")
    out=compute()
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(out,f,ensure_ascii=False,indent=2,sort_keys=True);f.write("\n")
    saved=json.loads(target.read_text(encoding="utf-8"))
    require(out==saved,"saved independent results mismatch; no file overwritten")
    print(json.dumps({"status":"PASS","all_checks_passed":True,
        "results_match":True,"check_groups":6,
        "mode":"write_and_replay" if args.write else "read_only_replay"}))
if __name__=="__main__":main()
