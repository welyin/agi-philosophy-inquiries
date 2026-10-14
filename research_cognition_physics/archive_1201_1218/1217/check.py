"""1217 finite event certificates. No old program imports.
All decision values use integers/Fraction and outward log intervals.
Default replay is read-only. --write creates results only if absent.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent
GRID=10**24
ROOT_GRID=10**8
TERMS=40


def require(ok,msg):
    if not ok: raise AssertionError(msg)


def floor_frac(x):
    return x.numerator//x.denominator


def ceil_frac(x):
    return -((-x.numerator)//x.denominator)


def outward(lo,hi):
    return F(floor_frac(lo*GRID),GRID), F(ceil_frac(hi*GRID),GRID)


def add_int(x,y):
    return x[0]+y[0],x[1]+y[1]


def scale_int(x,c):
    return (x[0]*c,x[1]*c) if c>=0 else (x[1]*c,x[0]*c)


def atanh_log(y):
    require(F(1)<=y<=F(2),"log range")
    z=(y-1)/(y+1)
    power=z
    total=F(0)
    for j in range(TERMS):
        total+=2*power/F(2*j+1)
        power*=z*z
    rem=2*power/(F(2*TERMS+1)*(1-z*z))
    return total,total+rem


LN2=atanh_log(F(2))


def log_interval(x):
    x=F(x)
    require(x>0,"positive log argument")
    y=x;k=0
    while y>=2:y/=2;k+=1
    while y<1:y*=2;k-=1
    return outward(*add_int(atanh_log(y),scale_int(LN2,F(k))))


def vacancy(p):
    require(sum(p)==1 and min(p)>0,"faithful spectrum")
    total=(F(0),F(0))
    for x in p:total=add_int(total,log_interval(len(p)*x))
    return outward(*scale_int(total,F(-1,len(p))))


def binary_divergence(a,p):
    require(0<a<1 and 0<p<1,"binary interior")
    return outward(*add_int(scale_int(log_interval(a/p),a),
                            scale_int(log_interval((1-a)/(1-p)),1-a)))


def binary_root_bracket(a,V):
    lo=0;hi=floor_frac(a*ROOT_GRID)
    require(binary_divergence(a,F(hi,ROOT_GRID))[1]<V[0],"upper root bracket")
    while hi-lo>1:
        mid=(lo+hi)//2
        bounds=binary_divergence(a,F(mid,ROOT_GRID))
        if bounds[0]>V[1]:lo=mid
        elif bounds[1]<V[0]:hi=mid
        else:raise AssertionError("log precision insufficient; no silent rounding")
    require(lo>0,"positive root certificate")
    require(binary_divergence(a,F(lo,ROOT_GRID))[0]>V[1],"lower root strict")
    require(binary_divergence(a,F(hi,ROOT_GRID))[1]<V[0],"upper root strict")
    return F(lo,ROOT_GRID),F(hi,ROOT_GRID)


def exponential_lower(a,V):
    # Certify p0 < exp(-(V_upper+ln2_upper)/a), using only log intervals.
    c=(V[1]+outward(*LN2)[1])/a
    lo=0;hi=ROOT_GRID
    while hi-lo>1:
        mid=(lo+hi)//2
        b=log_interval(F(mid,ROOT_GRID))
        if b[1]<-c:lo=mid
        elif b[0]>-c:hi=mid
        else:raise AssertionError("exp bracket undecided")
    require(lo>0,"positive explicit lower")
    require(log_interval(F(lo,ROOT_GRID))[1]<-c,"exponential lower strict")
    return F(lo,ROOT_GRID),c


def tail(p,alpha):
    p=sorted(p);d=len(p)
    require(0<=alpha<=d,"fractional trace")
    if alpha==d:return F(1)
    m=floor_frac(alpha)
    return sum(p[:m])+(alpha-m)*p[m]


def effect_weights(d,alpha):
    return [min(F(1),max(F(0),alpha-j)) for j in range(d)]


def interval_json(x):
    return {"lower":str(x[0]),"upper":str(x[1])}


def exact_trace(a):
    return sum(a[i,i] for i in range(len(a)))


def compute():
    sources=[
        "archive_429_466/research_note_457.md",
        "archive_1120_1158/_shared/1157_1158/adoption_and_dependencies.md",
        "archive_1120_1158/1157/proof.md",
        "archive_1201_/1202/results.json",
        "archive_1201_/1203/proof.md",
        "archive_1201_/1203/results.json",
        "archive_1201_/1204/results.json",
        "archive_1201_/1211/proof.md",
        "archive_1201_/1214/proof.md",
        "archive_1201_/1215/proof.md",
        "archive_1201_/1216/proof.md",
    ]
    hashes={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}
    hashes["archive_1201_/1217/check.py"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    cold=json.loads((ROOT/"archive_1201_/1202/results.json").read_text(encoding="utf-8"))
    hot=json.loads((ROOT/"archive_1201_/1204/results.json").read_text(encoding="utf-8"))
    C=[F(v) for v in cold["one_marginal_rational"]]
    H=[F(v) for v in hot["new_reference"]["one_marginal_rational"]]
    require(sum(C)==sum(H)==1 and min(C)>0 and min(H)>0,"native marginals")
    spectrum=sorted(x*y for x in C for y in H)
    groups={}
    groups["native_objects_and_quantifiers"]={
        "pass":True,"material":"same original AB, independent cold/hot pools",
        "D":16,"bad_rank":12,"cold":list(map(str,C)),"hot":list(map(str,H)),
        "sorted_joint_spectrum":list(map(str,spectrum)),
        "auxiliary":"same independent normal pure C and terminal blank K",
        "A_bad_and_K_nonzero_are_separate_effects":True,
        "same_trial_equality_claim":False,"uses_1205_correlated_preparation_as_channel":False,
        "new_hardware_or_axioms":False,
    }

    # A noninteger effect-trace witness, exact in Fraction arithmetic.
    p=[F(1,10),F(1,5),F(3,10),F(2,5)]
    d=4;r=3;eps=F(1,20);alpha=F(r)-d*eps
    e=effect_weights(d,alpha)
    M=np.array([[e[j]/r if i<r else (1-e[j])/(d-r)
                 for j in range(d)] for i in range(d)],dtype=object)
    require(all(sum(M[:,j])==1 for j in range(d)) and min(M.flat)>=0,"CPTP classical channel")
    output=M@np.array(p,dtype=object)
    uniform=M@np.array([F(1,d)]*d,dtype=object)
    defect=sum(abs(x-F(1,d)) for x in uniform)/2
    require(defect==eps,"sharp normalized unital defect")
    require(sum(output[:r])==tail(p,alpha)==F(27,50),"sharp spectral tail")
    O=np.eye(4,dtype=object)
    O[:2,:2]=[[F(3,5),F(-4,5)],[F(4,5),F(3,5)]]
    require(np.array_equal(O.T@O,np.eye(4,dtype=object)),"exact basis rotation")
    rho=O@np.diag(np.array(p,dtype=object))@O.T
    E=O@np.diag(np.array(e,dtype=object))@O.T
    require(exact_trace(rho@E)==tail(p,alpha) and exact_trace(E)==alpha,"basis independent effect")
    groups["sharp_fractional_spectral_certificate"]={
        "pass":True,"spectrum":list(map(str,p)),"alpha":str(alpha),"effect":list(map(str,e)),
        "sharp_value":str(tail(p,alpha)),"normalized_unital_defect":str(defect),
        "non_energy_basis_rotation_checked":True,
        "CP_realization":"measure-prepare with two positive effects; no native access claim",
        "probability_matrix":[list(map(str,row)) for row in M],
    }

    epsilons=[F("0.000109379990"),F("0.000059665807")]
    native=[]
    for label,eps,threshold in zip(("600","300"),epsilons,(F(".539998258"),F(".540062252"))):
        a=F(3,4)-eps;alpha=16*a
        lower=tail(spectrum,alpha)
        old=spectrum[0]*alpha
        require(lower>threshold and lower>old,"native strict improvement")
        native.append({"contact":label,"epsilon":str(eps),"a":str(a),"alpha":str(alpha),
                       "sharp_tail_exact":str(lower),"strict_display_lower":str(threshold),
                       "old_1214_lower":str(old),"event_certificates":["A_excited","K_nonzero"]})
    groups["actual_A_and_K_spectral_bounds"]={
        "pass":True,"exact_unital_tail":str(tail(spectrum,F(12))),"contacts":native,
        "all_free_phases_and_actual_F_J_tails_retained":True,
        "maxmixed_event_from_original_all_input_error":True,
        "source_error_requires_separate_full_trace_budget":True,
    }

    VC=vacancy(C);VH=vacancy(H);Vsum=outward(*add_int(VC,VH));Vjoint=vacancy(spectrum)
    require(max(Vsum[0],Vjoint[0])<=min(Vsum[1],Vjoint[1]),"product additivity enclosure")
    # One strictly dephasing, unital, noncommuting input diagnostic.
    Vin=vacancy(p);Vout=vacancy([rho[i,i] for i in range(4)])
    require(Vout[1]<Vin[0],"strict unital reverse-relative-entropy diagnostic")
    # Measurements on this diagonal source, independent of physical Hamiltonians.
    q=sum(e)/4;p_event=sum(x*y for x,y in zip(p,e))
    Dbin=binary_divergence(q,p_event)
    require(Dbin[1]<Vin[0],"event DPI diagnostic")
    groups["reverse_relative_entropy_budget"]={
        "pass":True,"V_C":interval_json(VC),"V_H":interval_json(VH),
        "V_CH_product_sum":interval_json(Vsum),"V_CH_joint_check":interval_json(Vjoint),
        "unital_dephasing_input":interval_json(Vin),"unital_dephasing_output":interval_json(Vout),
        "binary_measurement":interval_json(Dbin),"definition":"D(I_D/D || Omega)",
        "standard_vacancy_attribution":"Wilming-Gallego PRX7 041033; uniform reference specialization",
        "standard_DPI_additivity_not_claimed_original":True,
    }

    budgets=[]
    binary_native=[]
    for row,lowerB,lowerE in zip(native,(F(".43543064"),F(".43548341")),
                                            (F(".30228764"),F(".30231162"))):
        a=F(row["a"]);br=binary_root_bracket(a,Vsum);ex,c=exponential_lower(a,Vsum)
        require(br[0]>=lowerB and ex>=lowerE,"native vacancy display thresholds")
        binary_native.append({"contact":row["contact"],"binary_root_bracket":interval_json(br),
                              "certified_p_strictly_greater_than":str(br[0]),
                              "explicit_exponential_rational_lower":str(ex),
                              "explicit_exponent_upper":str(c)})
    for copies in (1,2,4):
        vb=outward(*scale_int(Vsum,F(copies)))
        a=F(native[0]["a"]);ex,c=exponential_lower(a,vb)
        budgets.append({"independent_CH_blocks":copies,"aggregate_budget":interval_json(vb),
                        "same_event_a":str(a),"explicit_uniform_lower":str(ex)})
    groups["binary_and_additive_uniform_bounds"]={
        "pass":True,"native_contacts":binary_native,"independent_budget_examples":budgets,
        "log_series_terms":TERMS,"outward_log_grid":GRID,"root_grid":ROOT_GRID,
        "all_log_exp_decisions_use_outward_rational_bounds":True,
        "uniform_budget_does_not_require_fixed_D":True,
        "uniform_time_budget_not_derived":True,
    }

    # A fixed native-h SWAP family: all eta>0 are faithful, no common V bound.
    rows=[];lastV=F(-1)
    for k in (1,2,4,8):
        eta=F(1,10**k);stock=[1-3*eta,eta,eta,eta]
        eig=sorted(x/4 for x in stock for _ in range(4))
        vi=vacancy(stock)
        forward=(F(0),F(0))
        for x in stock:forward=add_int(forward,scale_int(log_interval(4*x),x))
        require(vi[0]>lastV,"vacancy increases");lastV=vi[1]
        require(forward[1]<log_interval(F(4))[0],"forward divergence bounded")
        require(tail(eig,F(12))==3*eta,"swap family sharp tail")
        require(sum(stock)==1 and min(stock)>0,"each finite source faithful")
        rows.append({"eta":str(eta),"failure_exact":str(3*eta),"V0":interval_json(vi),
                     "forward_relative_entropy":interval_json(outward(*forward))})
    sw=np.zeros((16,16),dtype=int)
    for i in range(4):
        for j in range(4):sw[4*j+i,4*i+j]=1
    require(np.array_equal(sw.T@sw,np.eye(16,dtype=int)),"same-h swap")
    groups["bounded_mean_and_forward_entropy_do_not_bound_budget"]={
        "pass":True,"family":rows,"fixed_material_dimension":16,"fixed_native_h":True,
        "source_energy_above_ground_upper":"1.094","forward_divergence_upper":"log(4)",
        "all_eta_instances_have_positive_failure":True,
        "infimum_over_source_preparations_is_zero":True,
        "same_h_SWAP_is_mathematical_permission_not_original_clock_exact_gate":True,
    }

    corr=[]
    for delta in (F(1,10),F(1,100),F(1,10000)):
        ps=[(1-delta)/2,delta/2,delta/2,(1-delta)/2]
        require(ps[0]+ps[1]==ps[2]+ps[3]==F(1,2),"maxmixed marginals")
        require(tail(ps,F(2))==delta,"correlated sharp tail")
        v=vacancy(ps)
        require(v[0]>0,"correlated budget positive")
        corr.append({"delta":str(delta),"joint_V0":interval_json(v),
                     "both_marginal_V0_exact":"0","tail_exact":str(delta)})
    require(tail([F(1,4)]*4+[F(0)]*12,F(12))==0,"pure inventory deletes tail")
    groups["correlations_pure_inventory_and_terminal_record_scope"]={
        "pass":True,"correlated_family":corr,
        "pure_inventory_budget":"infinity (domain failure, not infinite physical cost)",
        "pure_C_and_terminal_K_must_not_be_inserted_into_active_Omega_without_changing_role":True,
        "feedback_alternative":"S has q=0, or K is included in active SK and input is nonfaithful",
        "QND_only_record_allowed_and_accounted":True,
        "no_universe_wide_faithfulness_assumption":True,
    }
    groups["three_quantifiers_and_stopping_line"]={
        "pass":True,"Q1":"Each finite active faithful Omega and same event q>0 gives p>0; no uniform n(t).",
        "Q2":"A class with V0<=B finite and q>=a_*>0 has a common positive bound even for varying dimensions.",
        "Q3":"This sufficient fixed-time route uses separately certified B(T),a(T), or a uniform spectral budget; other methods are not excluded.",
        "CC_scope":"Declared finite task domains and finite generation; not automatic infinite protocol closure.",
        "non_Markov_compatible":True,"needs_GKLS_or_CP_divisibility":False,
        "new_axioms":0,"four_laws_complete":False,
    }
    require(len(groups)==8 and all(g["pass"] for g in groups.values()),"eight groups")
    return {"round":1217,"status":"AUTHOR_CHECKS_PASS_NOT_ACCEPTANCE",
            "all_checks_passed":True,"check_groups":8,
            "source_sha256":hashes,"groups":groups}


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--write",action="store_true")
    args=ap.parse_args();dest=HERE/"results.json"
    if args.write and dest.exists():raise SystemExit("Refuse --write: results.json already exists.")
    result=compute()
    if args.write:
        with dest.open("x",encoding="utf-8") as f:
            json.dump(result,f,ensure_ascii=False,indent=2,sort_keys=True);f.write("\n")
    saved=json.loads(dest.read_text(encoding="utf-8"))
    require(saved==result,"read-only replay mismatch; no file overwritten")
    print(json.dumps({"status":"PASS","all_checks_passed":True,"results_match":True,
                      "check_groups":8,"strict_fields":"all exact; no floating diagnostics"}))


if __name__=="__main__":main()

