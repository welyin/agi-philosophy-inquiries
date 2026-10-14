"""Independent 1204: integer shells, native 16-input contact, record inference.
No import of author 1204 or historical code. Default is read-only recomputation;
--write creates independent_results.json once. NumPy is diagnostic, while exact
integer/Fraction budgets and analytical scope are documented in the review.
"""
from __future__ import annotations
import argparse
from bisect import bisect_right
from decimal import Decimal,localcontext,ROUND_CEILING,ROUND_FLOOR
from fractions import Fraction as F
import hashlib,json,math,random
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
N=100
DEN=10**18
W=F(547,500)
ENTRY=F(1,10**10)

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p): return json.loads(p.read_text(encoding="utf-8-sig"))
def dec(x):
    if isinstance(x,F): return Decimal(x.numerator)/Decimal(x.denominator)
    return Decimal(str(x))
def rounded_up(x,places=12):
    with localcontext() as c:
        c.prec=80
        return str(dec(x).quantize(Decimal(1).scaleb(-places),rounding=ROUND_CEILING))
def entropy(rho):
    v=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert v.min()>-1e-12
    v=v[v>1e-14]
    return float(-np.dot(v,np.log(v)))
def exp_h(h,t):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T
def compare(a,b,p="root"):
    if isinstance(a,dict):
        assert set(a)==set(b),p
        for k in a:compare(a[k],b[k],p+"."+k)
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b),p
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,p+f"[{i}]")
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-12),(p,a,b)
    else:assert a==b,(p,a,b)

def polynomial(x):
    return 32000*x**3-31200*x*x*DEN-2095*x*DEN**2-17*DEN**3
def bracket(lo,hi):
    lo*=DEN//10**6;hi*=DEN//10**6
    assert polynomial(lo)*polynomial(hi)<0
    while hi-lo>1:
        mid=(lo+hi)//2
        if polynomial(lo)*polynomial(mid)>0:lo=mid
        else:hi=mid
    return (lo,hi)
def roots():
    boxes=[bracket(-54072,-54071),bracket(-9461,-9460),
           (DEN,DEN),bracket(1038532,1038533)]
    assert all(boxes[i][1]<boxes[i+1][0] for i in range(3))
    # Three disjoint sign-changing intervals exhaust this cubic's roots.
    return boxes

def count_shell(boxes,lower_twice,upper_twice,seed):
    factorial=[math.factorial(k) for k in range(N+1)]
    counts=[];weights=[];starts=[0]
    one=[0]*4;two=[[0]*4 for _ in range(4)]
    margin=10**100;ambiguous=0
    for a in range(N+1):
      for b in range(N-a+1):
       for c in range(N-a-b+1):
        n=(a,b,c,N-a-b-c)
        lo=sum(n[i]*boxes[i][0] for i in range(4))
        hi=sum(n[i]*boxes[i][1] for i in range(4))
        inside=2*lo>=lower_twice*DEN and 2*hi<=upper_twice*DEN
        outside=2*hi<lower_twice*DEN or 2*lo>upper_twice*DEN
        if not (inside or outside):ambiguous+=1
        margin=min(margin,*[abs(2*x-y*DEN) for x in (lo,hi)
                           for y in (lower_twice,upper_twice)])
        if not inside:continue
        weight=factorial[N]//math.prod(factorial[v] for v in n)
        counts.append(n);weights.append(weight);starts.append(starts[-1]+weight)
        for i in range(4):
            one[i]+=weight*n[i]
            for j in range(4):two[i][j]+=weight*n[i]*(n[j]-(i==j))
    assert ambiguous==0
    dimension=starts[-1]
    r=[F(v,N*dimension) for v in one]
    pair=[[F(v,N*(N-1)*dimension) for v in row] for row in two]
    assert sum(r)==1
    assert all(sum(row)==r[i] for i,row in enumerate(pair))
    typemap={n:i for i,n in enumerate(counts)}

    def unrank(index):
        if not 0<=index<dimension:raise ValueError("outside shell")
        block=bisect_right(starts,index)-1
        residual=index-starts[block]
        ns=list(counts[block]);ways=weights[block];remaining=N;word=[]
        while remaining:
            for letter in range(4):
                num=ways*ns[letter]
                assert num%remaining==0
                size=num//remaining
                if residual>=size:
                    residual-=size
                    continue
                word.append(letter);ns[letter]-=1
                ways=size;remaining-=1
                break
            else:raise AssertionError("missing block")
        assert residual==0 and ways==1
        return tuple(word)
    def rank(word):
        ns=[word.count(i) for i in range(4)]
        block=typemap[tuple(ns)]
        value=starts[block];ways=weights[block];remaining=N
        for letter in word:
            for earlier in range(letter):
                assert ways*ns[earlier]%remaining==0
                value+=ways*ns[earlier]//remaining
            assert ways*ns[letter]%remaining==0
            ways=ways*ns[letter]//remaining
            ns[letter]-=1;remaining-=1
        assert ways==1
        return value

    rng=random.Random(seed)
    tests={0,1,dimension-2,dimension-1}
    for i in range(len(counts)):tests.update([starts[i],starts[i+1]-1])
    tests.update(rng.randrange(dimension) for _ in range(64))
    for idx in sorted(tests):assert rank(unrank(idx))==idx
    for bad in (-1,dimension):
        try:unrank(bad)
        except ValueError:pass
        else:raise AssertionError("invalid index accepted")
    bits=(dimension-1).bit_length()
    q=F(dimension,1<<bits)
    fail=(1-q)**10
    # First-valid successful word is exactly uniform, with all labels retained.
    per_word=sum((1-q)**j/F(1<<bits) for j in range(10))
    assert per_word/(1-fail)==F(1,dimension)
    result={
       "window":[str(F(lower_twice,2)),str(F(upper_twice,2))],
       "dimension":str(dimension),"types":len(counts),
       "all_types_examined":math.comb(N+3,3),"ambiguous_types":ambiguous,
       "boundary_margin_exact":str(F(margin,2*DEN)),
       "seed_bits":bits,"allocated_seed_bits":10*bits,
       "per_seed_acceptance_exact":str(q),"ten_seed_failure_exact":str(fail),
       "ten_seed_failure":float(fail),"rank_unrank_checked":len(tests),
       "uniform_given_success_exact":True,
       "one_rational":[str(x) for x in r],
       "two_rational":[[str(x) for x in row] for row in pair],
       "type_weight_sha256":hashlib.sha256(json.dumps(
           list(zip(counts,weights)),separators=(",",":")).encode()).hexdigest()
    }
    return r,pair,fail,result

def thermal(r,boxes,beta):
    with localcontext() as c:
        c.prec=80
        energies=[Decimal(a+b)/(2*DEN) for a,b in boxes]
        w=[(-dec(beta)*e).exp() for e in energies]
        tau=[x/sum(w) for x in w]
        rr=[dec(x) for x in r]
        radius=Decimal("1e-12")
        chi=sum(rr[i]**2/(tau[i]-radius) for i in range(4))-1
        tv=sum(abs(rr[i]-tau[i]) for i in range(4))/2+2*radius
        actual=sum(rr[i]*(rr[i]/tau[i]).ln() for i in range(4))
        # Root intervals are at most 1e-18 wide. A very loose 1e-12
        # population radius covers Gibbs/root and 80-digit evaluation errors.
        assert max(b-a for a,b in boxes)<=1
        assert all(x>Decimal(".1") for x in tau)
        return {"beta":str(beta),"tau_decimal80_diagnostic":[str(x) for x in tau],
                "tau_entry_radius":"1e-12",
                "TV_upper":rounded_up(tv),"chi_square_D_upper":rounded_up(chi),
                "D_numeric":float(actual)}

def clock_bound(v):
    T=F(606);M=F(10**24);sigma=F(1,10**10)
    L=150+F(3,2)*v*F(22,7)
    spread=sigma*(1+T*T/(8*M*M*sigma**4))
    eps=L*T*spread+F(1,10000*v)
    eta=L*sigma+L*(spread+L*T*T/(2*M))
    with localcontext() as ctx:
        ctx.prec=80
        bad=dec(eps)**2
        entropy_bound=-bad*bad.ln()-(1-bad)*(1-bad).ln()+bad*Decimal(1023).ln()
    return eps,eta,entropy_bound

def run():
    boxes=roots()
    rC,pairC,fC,cold=count_shell(boxes,39,41,1204)
    rH,pairH,fH,hot=count_shell(boxes,79,81,1205)
    prev=load(STAGE/"1202"/"independent_results.json")
    assert cold["dimension"]==prev["exact_counting"]["shell_dimension"]
    assert cold["types"]==prev["exact_counting"]["admissible_types"]
    assert np.max(abs(np.array([float(x) for x in rC])-prev["source_marginals"]["one"]))<1e-15
    assert fC<F(1,4_000_000) and fH<F(1,1024)
    assert hot["seed_bits"]==195

    I=np.eye(4)
    h=np.array([[0,-.2,0,-np.sqrt(3)/80],[-.2,1,0,0],[0,0,1,0],
                [-np.sqrt(3)/80,0,0,-.025]])
    Q=np.zeros((4,4));Q[1,2]=Q[2,1]=1
    e,basis=np.linalg.eigh(h)
    H0=np.kron(h,I)+np.kron(I,h);H=H0+.1*np.kron(Q,Q)
    T=np.kron(basis,basis)
    swap=np.zeros((16,16))
    for a in range(4):
        for b in range(4):swap[4*b+a,4*a+b]=1
    assert np.linalg.norm(H@swap-swap@H,2)<1e-15
    G=np.zeros((64,64));copy=np.zeros((64,64))
    for j in range(4):
        Sj=I.copy()
        if j:Sj[[0,j],:]=Sj[[j,0],:]
        Pj=np.outer(basis[:,j],basis[:,j])
        G+=np.kron(np.kron(Pj,I),math.pi/2*(I-Sj))
        copy+=np.kron(np.kron(Pj,I),Sj)
    swap64=np.kron(swap,I)
    broken=float(np.linalg.norm(G@swap64-swap64@G,2))
    assert broken>1
    assert np.linalg.norm(G@np.kron(H0,I)-np.kron(H0,I)@G,2)<1e-12
    assert np.linalg.norm(exp_h(G,1)-copy,2)<1e-12

    contacts={};Ue_by_time={}
    sw=[4*(j%4)+j//4 for j in range(16)]
    for t in (600,300):
        Ue=T.T@exp_h(H,t)@T
        P=abs(Ue)**2
        KA=P.reshape(4,4,16).sum(axis=1)
        KB=P.reshape(4,4,16).sum(axis=0)
        ea=np.repeat(e,4);eb=np.tile(e,4)
        da=ea@P-ea;db=eb@P-eb
        assert np.max(abs(P.sum(axis=0)-1))<1e-12
        assert np.max(abs(P-P[np.ix_(sw,sw)]))<2e-12
        assert np.max(abs(da[sw]-db))<2e-12
        Ue_by_time[t]=Ue
        contacts[str(t)]={"full_transition_P":P.tolist(),"A_kernel":KA.tolist(),
           "B_kernel":KB.tolist(),"delta_A_coefficients":da.tolist(),
           "delta_B_coefficients":db.tolist(),"work_coefficients":(da+db).tolist(),
           "swap_transition_error":float(np.max(abs(P-P[np.ix_(sw,sw)]))),
           "unitarity_error":float(np.linalg.norm(Ue.conj().T@Ue-np.eye(16),2))}

    thermalC=thermal(rC,boxes,F(6,5));thermalH=thermal(rH,boxes,F(1,3))
    therm_defect=F(thermalC["chi_square_D_upper"])+F(thermalH["chi_square_D_upper"])
    contexts={};flows={};max_fail=F(0)
    for lam in (F(0),F(1,2),F(1)):
        successB=(1-lam)*(1-fC)+lam*(1-fH)
        eff=lam*(1-fH)/successB
        rb=[(1-eff)*c+eff*x for c,x in zip(rC,rH)]
        fail=1-(1-fC)*successB
        max_fail=max(max_fail,fail)
        assert fail<F(1,512)
        contexts[str(lam)]={"conditional_lambda_exact":str(eff),
            "conditional_lambda":float(eff),"bias_exact":str(eff-lam),
            "bias":float(eff-lam),"whole_pair_failure_exact":str(fail),
            "whole_pair_failure":float(fail),"B_one_rational":[str(x) for x in rb]}
        rf=np.array([float(x) for x in rC]);bf=np.array([float(x) for x in rb])
        for t in (600,300):
            v=600//t
            co=contacts[str(t)]
            KA=np.array(co["A_kernel"])
            coeff_a=np.array(co["delta_A_coefficients"])
            coeff_b=np.array(co["delta_B_coefficients"])
            pop=np.kron(rf,bf);reverse_pop=np.kron(bf,rf)
            da=float(coeff_a@pop);db=float(coeff_b@pop);work=da+db;qd=(da-db)/2
            reverse_da=float(coeff_a@reverse_pop)
            assert abs(db-reverse_da)<2e-12
            eps,eta,sc=clock_bound(v)
            a=W*eps+ENTRY
            measurement=2*W*F(5,1000)
            proxy_error=a+measurement
            actual_inference=2*a+measurement
            clock_error=2*W*eps+eta+ENTRY
            key=f"t{t}_lambda{lam}"
            flows[key]={"lambda_effective":float(eff),"delta_A":da,"delta_B":db,
                "Q_dir_ideal":qd,"bare_work_ideal":work,
                "forward_A_pointer":(KA@pop).tolist(),
                "reversed_A_pointer":(KA@reverse_pop).tolist(),
                "actual_single_energy_radius_exact":str(a),
                "actual_single_energy_radius":float(a),
                "actual_Q_interval":[qd-float(a),qd+float(a)],
                "actual_A_interval":[da-float(a),da+float(a)],
                "actual_B_interval":[db-float(a),db+float(a)],
                "clock_change_interval":[-work-float(clock_error),-work+float(clock_error)],
                "clock_loss_certified":work>float(clock_error),
                "read_statistics_vs_ideal_radius":float(proxy_error),
                "inference_radius_for_same_process_actual_Q":float(actual_inference),
                "model_read_statistic_lower":qd-float(proxy_error),
                "ideal_center_minus_actual_inference_radius":qd-float(actual_inference),
                "worst_allowed_observation_actual_confidence_lower":
                    qd-float(proxy_error+actual_inference),
                "all_allowed_observations_force_positive_actual_confidence_lower":
                    qd>float(proxy_error+actual_inference),
                "actual_A_gain_certified":da>float(a),
                "actual_B_loss_certified":db<-float(a),
                "actual_B_gain_certified":db>float(a),
                "clock_rank_upper":1024,
                "clock_entropy_upper":rounded_up(sc+Decimal("1e-60")),
                "Y_CK_rank_upper":256,
                "nominal_lambda_energy_bias_bound":float(W*abs(eff-lam))}
            if lam:
                assert da>float(a) and db<-float(a)
                assert qd>float(proxy_error)
                assert qd>float(actual_inference)
            else:
                assert abs(qd)<2e-12
    assert not flows["t300_lambda1/2"]["all_allowed_observations_force_positive_actual_confidence_lower"]

    # Independent entropy identities on cold-hot contact and subsequent copy.
    diagnostic={}
    rc=np.array([float(x) for x in rC]);rh=np.array([float(x) for x in rH])
    initial=np.kron(rc,rh)
    tauC=np.exp(-1.2*e);tauC/=tauC.sum()
    tauH=np.exp(-e/3);tauH/=tauH.sum()
    Din=float(rc@np.log(rc/tauC)+rh@np.log(rh/tauH))
    # Energy-basis copy is independently explicit, not the author's isometry.
    embed=np.zeros((64,16))
    for a in range(4):
        for b in range(4):embed[(4*a+b)*4+a,4*a+b]=1
    for t,Ue in Ue_by_time.items():
        rho=(Ue*initial)@Ue.conj().T
        rs=rho.reshape(4,4,4,4)
        rA=np.einsum("abcb->ac",rs);rB=np.einsum("abad->bd",rs)
        dA=float(np.trace(np.diag(e)@rA).real-e@rc)
        dB=float(np.trace(np.diag(e)@rB).real-e@rh)
        S0=entropy(np.diag(initial));SA=entropy(rA);SB=entropy(rB)
        Dout=-SA-float(np.diag(rA).real@np.log(tauC))-SB-float(np.diag(rB).real@np.log(tauH))
        lhs=1.2*dA+dB/3
        bare_residual=lhs-(SA+SB-S0+Dout-Din)
        rec=embed@rho@embed.T
        rr=rec.reshape(4,4,4,4,4,4)
        rAr=np.einsum("abkcbk->ac",rr)
        rBr=np.einsum("abkadk->bd",rr)
        rKr=np.einsum("abkabl->kl",rr)
        sar,sbr,sk=entropy(rAr),entropy(rBr),entropy(rKr)
        drec=-sar-float(np.diag(rAr).real@np.log(tauC))-sbr-float(np.diag(rBr).real@np.log(tauH))
        corr=sar+sbr+sk-S0
        read_residual=lhs+sk-(corr+drec-Din)
        assert abs(bare_residual)<2e-12 and abs(read_residual)<2e-12
        diagnostic[str(t)]={"D_initial":Din,
            "contact_identity_residual":bare_residual,
            "post_copy_identity_residual":read_residual,
            "S_K_ideal":sk,"weighted_energy_delta":lhs,
            "post_copy_total_correlation":corr,
            "B_energy_coherence_norm":float(np.linalg.norm(rB-np.diag(np.diag(rB)))),
            "no_claim_of_actual_continuous_clock_simulation":True}

    # Ten flow patterns: 2 times * (3 forward + 2 nonduplicate reversed).
    patterns=2*(3+2);n=250_000;cap=500_000
    assert patterns==10 and cap==2*n
    statistical_failure=2*4*patterns*math.exp(-2*n*.004**2)
    # Per-pattern batch fails only if >=n raw attempts fail. Union bound:
    # 2^(2n)*(1/512)^n=(1/128)^n; over 10 patterns.
    batch_log10=math.log10(patterns)-n*math.log10(128)
    statistics={
       "flow_patterns":patterns,"per_pattern_accepted_trials":n,
       "per_pattern_raw_attempt_cap":cap,"max_raw_attempts":patterns*cap,
       "per_probability_statistical_error":.004,"separate_systematic_error":.001,
       "forty_probabilities_joint_statistical_failure_upper":statistical_failure,
       "joint_success_lower_conditional_batch_available":1-statistical_failure,
       "all_batches_insufficient_success_upper_exact":"10*(1/128)^250000",
       "all_batches_insufficient_success_log10_upper":batch_log10,
       "overall_success_lower_formula":"1-80*exp(-8)-10*(1/128)^250000",
       "maximum_pair_factory_failure":float(max_fail),
       "first_successes_iid_given_success_and_batch_available":True,
       "proof_requires_fresh_independent_attempts_and_no_record_feedback":True,
       "temperature_probe_patterns_not_included_in_ten_flow_patterns":True,
       "maximum_target_plus_stock_materials_over_raw_cap":patterns*cap*1000,
       "seed_bits_upper_over_raw_cap":patterns*cap*3681,
       "readout_systematic_bound_assumed_not_generated":True
    }
    # Existing c thermometer (A energy level index2) remains a separate menu.
    Kc=np.array(contacts["600"]["A_kernel"])[:,8:12]
    temperature_menu={"cold":(Kc@rc).tolist(),"hot":(Kc@rh).tolist(),
       "max_event_separation":float(max(abs(Kc@rc-Kc@rh))),
       "separate_from_ten_flow_mode_confidence_claim":True}

    ledger={}
    for name,r in (("cold",rc),("hot",rh)):
        mean=float(e@r)
        target0=N*e[0]
        inventory0=N*e.sum()
        target1=N*mean
        ledger[name]={"mean_per_target":mean,"initial_target_energy":float(target0),
          "accepted_target_energy":target1,"initial_inventory_energy":float(inventory0),
          "final_inventory_energy":float(inventory0+target0-target1),
          "all_target_stock_total_energy":float(target0+inventory0)}

    return {
      "status":"independent scoped 1204 audit; finite NumPy diagnostics plus exact integer budgets",
      "all_checks_passed":True,
      "root_denominator":DEN,"independent_root_boxes":boxes,
      "cold_shell":cold,"hot_shell":hot,
      "temperature_comparisons":{"cold":thermalC,"hot":thermalH,
          "cold_hot_D_initial_numeric":Din,
          "sum_of_outward_chi_square_upper":str(therm_defect),
          "no_exact_Gibbs_assumption":True},
      "labelled_success_preparations":contexts,"energy_preparation_ledger":ledger,
      "native_material":{"h":h.tolist(),"Q":Q.tolist(),"energies":e.tolist()},
      "symmetry":{"rectangular_AB_swap_commutator_norm":float(np.linalg.norm(H@swap-swap@H,2)),
          "actual_read_G_AB_swap_commutator_norm":broken,
          "complete_parent_AB_symmetry_claimed":False},
      "contacts":contacts,"flows":flows,"entropy_diagnostics":diagnostic,
      "sampling_and_confidence":statistics,"separate_temperature_menu":temperature_menu,
      "scope":{"author_1204_code_imported":False,"historical_code_imported":False,
          "source_counts_reconstructed_with_integers":True,
          "all_16_contact_inputs_rebuilt":True,
          "R_and_failure_records_discarded":False,
          "B_energy_pointer_physically_built":False,
          "same_process_B_inference_uses_two_comparison_errors":True,
          "finite_clock_actual_process_numerically_solved":False,
          "pure_A_rank_bound_from_1203_reused":False,
          "universal_hot_to_cold_theorem_claimed":False,
          "temperature_mixture_assigned_exact_beta":False,
          "all_allowed_data_guaranteed_to_certify_actual_flow":False,
          "no_experiment_claimed":True,"full_source_factory_autonomous":False,
          "geometry_region_identity_completed":False},
      "historical_sha256":{str(p.relative_to(STAGE)).replace("\\","/"):sha(p) for p in
          [STAGE/"1202/independent_results.json",STAGE/"1203/independent_results.json",
           STAGE/"1203/proof.md"]}
    }

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=run();path=HERE/"independent_results.json"
    if args.write:
        if path.exists():raise FileExistsError("--write refuses to overwrite")
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:compare(result,load(path))
    print(json.dumps(result,ensure_ascii=False,indent=2))
if __name__=="__main__":main()
