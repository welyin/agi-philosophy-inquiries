"""1204: independently prepared finite references, native flow and clock ledger.

Default recomputes and compares the author's saved result.  --write creates a
new result only.  It never invokes the write entry points of historical code.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction as F
from bisect import bisect_right
import importlib.util
import argparse
import hashlib
import json
import math
import random
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ERR=F('1e-10')
WIDTH=F('1.094')

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

M1=load_module('native1201',STAGE/'1201/check.py')
M2=load_module('prepare1202',STAGE/'1202/check.py')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dec(x):return D(x.numerator)/D(x.denominator)
def strs(x):return [[str(v) for v in row] for row in x] if x and isinstance(x[0],list) else [str(v) for v in x]
def fs(x):return [F(v) for v in x]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def outer(a,b):return [x*y for x in a for y in b]
def mv(a,b):return [dot(row,b) for row in a]
def up(x,places=12):
    with localcontext() as ctx:
        ctx.prec=80
        return str(dec(x).quantize(D(10)**(-places),rounding=ROUND_CEILING))
def show_interval(x,r):return [str(x-r),str(x+r)]

def new_pool(old):
    """Only the new [79/2,81/2] window is enumerated; old counts are loaded."""
    boxes=[tuple(x) for x in old['preparation']['root_boxes']]
    DEN=old['preparation']['root_denominator'];N=100
    for j,(a,b) in enumerate(boxes):
        if j!=2:assert M2.poly(a)*M2.poly(b)<0 and b-a==1
        else:assert a==b==DEN
    fac=[math.factorial(k) for k in range(N+1)]
    rows=[];total=0;one=[0]*4;two=[[0]*4 for _ in range(4)];margin=10**50
    for a in range(N+1):
        for b in range(N-a+1):
            for c in range(N-a-b+1):
                n=(a,b,c,N-a-b-c)
                lo=sum(n[j]*boxes[j][0] for j in range(4));hi=sum(n[j]*boxes[j][1] for j in range(4))
                inside=2*lo>=79*DEN and 2*hi<=81*DEN
                outside=2*hi<79*DEN or 2*lo>81*DEN
                assert inside or outside
                margin=min(margin,*[abs(2*x-k*DEN) for x in (lo,hi) for k in (79,81)])
                if not inside:continue
                w=fac[N]//math.prod(fac[x] for x in n)
                rows.append((n,w,total));total+=w
                for i in range(4):
                    one[i]+=w*n[i]
                    for j in range(4):two[i][j]+=w*n[i]*(n[j]-(i==j))
    r=[F(x,N*total) for x in one]
    pair=[[F(x,N*(N-1)*total) for x in row] for row in two]
    assert len(rows)==2391 and sum(r)==1
    assert all(sum(pair[i])==r[i] for i in range(4))
    starts=[row[2] for row in rows];lookup={row[0]:row for row in rows}
    rng=random.Random(1204)
    tests=sorted(set([0,total-1,*starts,*[s+w-1 for _,w,s in rows],*[rng.randrange(total) for _ in range(32)]]))
    for idx in tests:
        n,w,start=rows[bisect_right(starts,idx)-1]
        word=M2.unrank_permutation(n,idx-start)
        assert tuple(word.count(i) for i in range(4))==n
        assert lookup[n][2]+M2.rank_permutation(word)==idx
    bits=total.bit_length();q=F(total,2**bits);failure=(1-q)**10
    assert bits==195 and failure<F(1,1024)
    assert sum((1-q)**k/F(2**bits) for k in range(10))==(1-failure)/total
    return boxes,total,r,pair,failure,{
        'energy_window':['79/2','81/2'],'N':N,'selected_types':len(rows),'dimension':str(total),
        'window_boundary_margin_lower_exact':str(F(margin,2*DEN)),
        'seed_bits_each':bits,'seed_registers':10,'seed_bits_total':10*bits,
        'seed_acceptance_exact':str(q),'all_rejected_exact':str(failure),
        'all_rejected_upper':up(failure),'rank_unrank_checked_indices':len(tests),
        'rank_program_reused_without_old_window_reenumeration':True,
        'one_marginal_rational':strs(r),'two_marginal_rational':strs(pair),
        'all_inventory_seed_failure_and_pool_records_retained':True}

def full_contact(t,h,q,e,R):
    """All 16 input energy pairs; old 80-digit certified exp, no new tuning."""
    H0=M1.plus(M1.kron(h,M1.eye(4)),M1.kron(M1.eye(4),h))
    V=M1.kron(q,q);H=M1.plus(H0,M1.sm(V,D('.1')))
    U,certificate=M1.exp_native(H,t);T=M1.kron(R,R)
    Ue=M1.mm(M1.mm(M1.dag(T),U),T)
    P=[[sum(x*x for x in Ue[j][k]) for k in range(16)] for j in range(16)]
    K=[[sum(P[4*j+d][k] for d in range(4)) for k in range(16)] for j in range(4)]
    coeffA=[];coeffB=[];coeffW=[]
    for k in range(16):
        a,b=divmod(k,4);assert abs(sum(P[j][k] for j in range(16))-1)<D('1e-35')
        da=sum(e[j//4]*P[j][k] for j in range(16))-e[a]
        db=sum(e[j%4]*P[j][k] for j in range(16))-e[b]
        coeffA.append(da);coeffB.append(db);coeffW.append(da+db)
    # Swap symmetry, on every full input/output basis state, is a check only;
    # its exact status follows from [H,SWAP_AB]=0 in proof.md.
    sw=lambda k:4*(k%4)+k//4
    assert max(abs(P[j][k]-P[sw(j)][sw(k)]) for j in range(16) for k in range(16))<D('1e-35')
    out={'time':t,'full_transition_P':[[M1.rounded(x) for x in row] for row in P],
         'pointer_K':[[M1.rounded(x) for x in row] for row in K],
         'delta_EA_coefficients':[M1.rounded(x) for x in coeffA],
         'delta_EB_coefficients':[M1.rounded(x) for x in coeffB],
         'bare_work_coefficients':[M1.rounded(x) for x in coeffW],
         'each_saved_entry_error_upper':str(ERR),'exponential_certificate':certificate}
    return out,Ue

def thermal_bound(r,e,beta):
    w=[(-beta*x).exp() for x in e];tau=[F(x/sum(w)) for x in w]
    rad=F('1e-12');assert min(tau)>F('.1')
    tv=sum(abs(a-b) for a,b in zip(r,tau))/2+2*rad
    chi=sum(r[i]**2/(tau[i]-rad) for i in range(4))-1
    return {'beta':str(beta),'gibbs_population_rational':strs(tau),
            'tau_entry_error_upper':str(rad),'local_trace_distance_upper':up(tv),
            'relative_entropy_chi_square_upper':up(chi)},tau,chi

def matrix_entropy(a):
    val=np.linalg.eigvalsh((a+a.conj().T)/2);val=val[val>1e-14]
    return float(-sum(val*np.log(val)))

def diagnostic_state(Ue,rA,rB,e,betaA,betaB):
    U=np.array([[complex(float(x[0]),float(x[1])) for x in row] for row in Ue])
    rin=np.array([float(x) for x in outer(rA,rB)])
    rho=(U*rin)@U.conj().T
    four=rho.reshape(4,4,4,4);a=np.einsum('ibjb->ij',four);b=np.einsum('aiaj->ij',four)
    ein=np.array([float(x) for x in e]);a0=np.array([float(x) for x in rA]);b0=np.array([float(x) for x in rB])
    da=float(np.trace(a@np.diag(ein)).real-ein@a0);db=float(np.trace(b@np.diag(ein)).real-ein@b0)
    s0=matrix_entropy(np.diag(rin));sa=matrix_entropy(a);sb=matrix_entropy(b)
    ta=np.exp(-betaA*ein);ta/=sum(ta);tb=np.exp(-betaB*ein);tb/=sum(tb)
    d0=float(sum(a0*np.log(a0/ta))+sum(b0*np.log(b0/tb)))
    d1=-sa-float(np.diag(a).real@np.log(ta))-sb-float(np.diag(b).real@np.log(tb))
    lhs=betaA*da+betaB*db;rhs=sa+sb-s0+d1-d0
    assert abs(lhs-rhs)<1e-12
    # Ideal energy copy Vrec: |a,b> -> |a,b,k=a>; all cross terms retained.
    W=np.zeros((64,16));
    for j in range(16):W[4*j+j//4,j]=1
    recorded=W@rho@W.T
    rr=recorded.reshape(4,4,4,4,4,4)
    ar=np.einsum('ibkjbk->ij',rr);br=np.einsum('aikajk->ij',rr);kr=np.einsum('abiabj->ij',rr)
    sy=matrix_entropy(kr);sar=matrix_entropy(ar);sbr=matrix_entropy(br)
    dout=-sar-float(np.diag(ar).real@np.log(ta))-sbr-float(np.diag(br).real@np.log(tb))
    totalcorr=sar+sbr+sy-s0
    assert abs(lhs+sy-(totalcorr+dout-d0))<1e-12
    return {'diagnostic_only':True,'S_AB_before_and_after_contact':s0,'S_A_contact':sa,'S_B_contact':sb,
            'I_A_B_contact':sa+sb-s0,'D_initial':d0,'D_final_marginals':d1,
            'weighted_energy_lhs':lhs,'bare_contact_identity_error':abs(lhs-rhs),
            'S_K_after_copy':sy,'I_A_B_K_after_copy':totalcorr,
            'recorded_identity_error':abs(lhs+sy-(totalcorr+dout-d0)),
            'B_energy_after_contact':float(np.trace(b@np.diag(ein)).real),
            'B_energy_basis_coherence_frobenius':float(np.linalg.norm(b-np.diag(np.diag(b))))}

def build():
    old=json.loads((STAGE/'1202/results.json').read_text(encoding='utf-8'))
    native_old=json.loads((STAGE/'1201/results.json').read_text(encoding='utf-8'))
    clock=json.loads((STAGE/'1203/results.json').read_text(encoding='utf-8'))
    boxes,Z,rH,pairH,fH,prepH=new_pool(old)
    rC=fs(old['one_marginal_rational'])
    ZC=int(old['preparation']['dimension']);fC=(1-F(ZC,2**173))**10
    assert fC<F('.000000211')
    active={};diagnostics={};contexts={}
    with localcontext() as ctx:
        ctx.prec=80
        h,q,e,R,basis_certificate=M1.native()
        thermH,tauH,chiH=thermal_bound(rH,e,D(1)/3)
        thermH['beta']='1/3'
        thermC,tauC,chiC=thermal_bound(rC,e,D('1.2'))
        assert F(thermH['local_trace_distance_upper'])<F('.003951')
        assert chiC<F('.000405')
        # Actual late record does not change bare material energy: [G,h_A]=0.
        contacts={};Ues={}
        for t in (600,300):contacts[str(t)],Ues[t]=full_contact(t,h,q,e,R)
        for lam in [F(0),F(1,2),F(1)]:
            successB=(1-lam)*(1-fC)+lam*(1-fH)
            lam_eff=lam*(1-fH)/successB
            rB=[(1-lam_eff)*x+lam_eff*y for x,y in zip(rC,rH)]
            # Success of A is independent of branch choice; retain all failures.
            failure_total=1-(1-fC)*successB
            assert failure_total<F(1,512)
            contexts[str(lam)]={'label_probability_before_selection':str(lam),
                'label_probability_conditional_success':str(lam_eff),
                'label_bias_exact':str(lam_eff-lam),'label_bias_upper':up(abs(lam_eff-lam)),
                'total_preparation_failure_exact':str(failure_total),
                'total_preparation_failure_upper':up(failure_total),
                'B_conditional_marginal':strs(rB),'A_B_active_marginal_product':True,
                'each_active_source_keeps_own_pool_stock_and_program_reference':True}
            for t in (600,300):
                co=contacts[str(t)];eps=F(clock['bounds_by_program'][str(600//t)]['isometry_error_exact_upper'])
                eta=F(clock['bounds_by_program'][str(600//t)]['initial_plus_final_interaction_energy_exact_upper'])
                ae=fs(co['delta_EA_coefficients']);be=fs(co['delta_EB_coefficients']);we=fs(co['bare_work_coefficients'])
                pAB=outer(rC,rB);da=dot(ae,pAB);db=dot(be,pAB);w=dot(we,pAB);qd=(da-db)/2
                kr=[fs(row) for row in co['pointer_K']]
                p=mv(kr,pAB);prev=mv(kr,outer(rB,rC))
                reverse_da=dot(ae,outer(rB,rC))
                assert abs(reverse_da-db)<2*ERR and abs(da+db-w)<2*ERR
                state_energy_error=WIDTH*eps+ERR
                clock_energy_error=2*WIDTH*eps+eta+ERR
                # Each observed four-outcome probability has total .005 error.
                # With normalized p and p_hat, centering energy yields 2*width*.005.
                observed_energy_error=state_energy_error+2*WIDTH*F('.005')
                raw_vs_declared_lam_error=WIDTH*abs(lam_eff-lam)
                key=f't{t}_lambda{lam}'
                active[key]={'t':t,'lambda':str(lam),'lambda_effective':str(lam_eff),
                    'A_pointer_probabilities_center':strs(p),'A_pointer_reversed_probabilities_center':strs(prev),
                    'pointer_probability_error_upper':up(eps+ERR),
                    'delta_EA_center':str(da),'delta_EB_center':str(db),
                    'Q_dir_center':str(qd),'bare_work_center':str(w),
                    'single_energy_error_upper':up(state_energy_error),
                    'Q_dir_actual_interval':show_interval(qd,state_energy_error),
                    'bare_work_actual_interval':show_interval(w,2*WIDTH*eps+ERR),
                    'clock_energy_change_interval':show_interval(-w,clock_energy_error),
                    'clock_energy_change_sign_certified_negative':w>clock_energy_error,
                    'Q_from_two_reversed_A_reads_error_upper':up(observed_energy_error),
                    'Q_observed_lower_after_fixed_probability_budget':str(qd-observed_energy_error),
                    'Q_actual_inference_radius_around_future_observed_proxy':up(2*state_energy_error+2*WIDTH*F('.005')),
                    'inference_radius_is_around_future_measurement_not_a_measured_lower_bound':True,
                    'label_bias_energy_allowance_upper':up(raw_vs_declared_lam_error)}
                if lam:
                    assert da-state_energy_error>0 and db+state_energy_error<0
                    assert qd-observed_energy_error>0
                else:assert abs(qd)<ERR
                # Clock finite rank: input AB rank <=16, ABK dim64, hence 1024.
                pbad=eps*eps;x=dec(pbad)
                sc=-x*x.ln()-(1-x)*(1-x).ln()+x*D(1023).ln()
                active[key]['clock_entropy_rank_bound']=1024
                active[key]['clock_entropy_upper']=str(sc.quantize(D('1e-12'),rounding=ROUND_CEILING))
        for t in (600,300):
            diagnostics[str(t)]=diagnostic_state(Ues[t],rC,rH,e,1.2,1/3)
        energies={}
        for name,r in [('cold',rC),('hot',rH)]:
            mean=sum(dec(r[i])*e[i] for i in range(4))
            energies[name]={'bare_mean':str(mean),'target_pool_excitation':str(100*(mean-e[0])),
                'inventory_initial_excitation':str(100*sum(x-e[0] for x in e)),
                'inventory_final_excitation':str(100*sum(x-e[0] for x in e)-100*(mean-e[0]))}
        mh=mv([[F(x) for x in row] for row in native_old['contacts']['600']['K']],rH)
        mc=mv([[F(x) for x in row] for row in native_old['contacts']['600']['K']],rC)
        temp_readout={'hot_pointer':strs(mh),'cold_pointer':strs(mc),
            'maximum_single_event_separation':str(max(abs(a-b) for a,b in zip(mh,mc))),
            'fixed_h_and_prepared_populations_not_defined_by_temperature':True}
    # Lambda-linearity follows from CPTP linearity, not from a temperature fit.
    for t in (600,300):
        zero=active[f't{t}_lambda0'];one=active[f't{t}_lambda1'];half=active[f't{t}_lambda1/2']
        le=F(half['lambda_effective'])
        assert abs(F(half['Q_dir_center'])-le*F(one['Q_dir_center']))<2*ERR
        assert abs(F(half['bare_work_center'])-((1-le)*F(zero['bare_work_center'])+le*F(one['bare_work_center'])))<2*ERR
    sources=[STAGE/'1201/check.py',STAGE/'1201/results.json',STAGE/'1202/check.py',STAGE/'1202/results.json',
             STAGE/'1202/acceptance.json',STAGE/'1203/proof.md',STAGE/'1203/results.json',STAGE/'1203/acceptance.json']
    return {'schema':'two_reference_operational_flow_1204_v1','status':'AUTHOR_DRAFT_CHECKS_PASS_NOT_ACCEPTANCE',
        'new_reference':prepH,'cold_reference_loaded_from_1202':True,
        'temperature_comparisons':{'cold':thermC,'hot':thermH,'not_unique_exact_temperatures':True},
        'independent_temperature_menu':temp_readout,
        'preparation_energy_ledger':energies,'labelled_preparations':contexts,
        'reference_preparation_resources_per_trial':{'active_target_materials':200,'known_level_inventory_materials':800,
            'cold_factory_seed_bits':1730,'mixed_factory_seed_bits_fixed':1950,'half_mixture_choice_bits':1,
            'mixed_factory_cold_branch_ignores_but_keeps_220_extra_seed_bits':True,
            'factory_resources_follow_prepared_pools_not_A_B_ports':True,
            'forward_ports':'A=cold_factory, B=mixed_factory',
            'reversed_ports':'A=mixed_factory, B=cold_factory; independently prepared pools are exchanged',
            'source_rank_program_inventory_and_all_failures_kept':True,
            'autonomous_preparation_program_claimed':False},
        'full_native_contacts':contacts,'basis_certificate':basis_certificate,
        'actual_flow_predictions':active,
        'multi_reference_thermodynamics':{'cold_hot_initial_relative_entropy_defect_upper':up(chiC+chiH),
            'exact_identity':'beta_A dEA+beta_B dEB+dS_Y=I(A:B:Y)_out+sum_i[D(rho_i_out||tau_i)-D(rho_i_in||tau_i)], Y=C K, initially A,B,Y independent',
            'Y_rank_upper':256,'clock_rank_upper':1024,'actual_clock_not_traced_out_physically':True,
            'direction_proved_from_finite_native_coefficients_not_entropy_bound':True,
            'record_G_commutes_with_bare_hA':True,
            'midpoint_mixture_not_assigned_an_exact_temperature':True},
        'finite_quantum_diagnostics':diagnostics,
        'measurement_contract':{'four_outcome_per_probability_total_error':'.005','statistical_part':'.004',
            'terminal_readout_systematic_part':'.001',
            'systematic_scope':'normalized terminal pointer distribution: readout, classification and calibration only',
            'source_clock_pointer_preparation_errors_included_in_readout_budget':False,
            'preparation_deviations_require_separate_full_state_trace_or_prepared_isometry_contract':True,
            'energy_and_entropy_claims_use_declared_source_clock_pointer_preparation':True,
            'independent_trials_each_pattern':250000,
            'distinct_patterns':10,'all_pattern_success_probability_lower':1-80*math.exp(-8),
            'raw_attempts_reserved_each_pattern':500000,'max_raw_attempts_total':5000000,
            'batch_insufficient_success_probability_upper':'10*(1/128)^250000',
            'first_250000_successes_used_all_other_records_retained':True,
            'fresh_pool_and_clock_resources_each_trial':True,
            'preparation_failure_flagged_not_postselected_silently':True,
            'tiny_work_or_clock_energy_directly_measured':False},
        'remaining_inputs':['1202 finite reversible energy-window preparation and known level inventories',
            '1203 positive-energy clock, fixed gates, pure pointer and matched velocity/program preparation',
            'terminal A-label interrogation with normalized pointer-distribution readout/classification/calibration budget',
            'additional source/clock/pointer preparation deviations need a separate full-state trace or prepared-isometry contract',
            'actual 3D region, signal and metric identity'],
        'not_claimed':['universal hot-to-cold theorem under all driven contacts','exact Gibbs or attracting equilibrium',
            'free repeated use or clock reset','record entropy as heat','area entropy, Unruh or GR'],
        'source_sha256':{str(p.relative_to(HERE.parents[2])):sha(p) for p in sources}}

def stable(result):
    # Floating diagnostics are checked by assertions, not bitwise certificates.
    value=dict(result);value.pop('finite_quantum_diagnostics');return value

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=build();path=HERE/'results.json'
    if args.write:
        if path.exists():raise SystemExit('Refusing to overwrite an existing result')
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(path.read_text(encoding='utf-8'));assert stable(result)==stable(saved),'saved result mismatch'
    print(json.dumps({'status':result['status'],'new_pool_failure_upper':result['new_reference']['all_rejected_upper'],
        'cold_hot_D_upper':result['multi_reference_thermodynamics']['cold_hot_initial_relative_entropy_defect_upper'],
        'hot_temperature_TV_upper':result['temperature_comparisons']['hot']['local_trace_distance_upper'],
        'minimal_observed_half_flow_lower':float(F(result['actual_flow_predictions']['t300_lambda1/2']['Q_observed_lower_after_fixed_probability_budget']))},ensure_ascii=False))
