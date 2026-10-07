"""Exact finite screening of relationship-selection mechanisms.

Classical records form a commuting quantum subtheory. This is not a model
of emergent spacetime or a derivation of a physical Hamiltonian.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import product, combinations
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'coalition_mechanism_screen_results.json'
ETA=F(1,10)
COST=F(1,8)
TRAIN_REPEATS=3
VALIDATE=32
MAX_ERRORS=5
SERVICE=256


def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def number(x):return dict(exact=str(x),decimal=float(x))
def bit(f,x):return (f>>x)&1
def pop(mask):return int(mask).bit_count()


def predictor(table,mask):
    out=[]
    for x in range(4):
        values=[bit(table,y) for y in range(4) if (y&mask)==(x&mask)]
        out.append(int(2*sum(values)>len(values)))
    return out


def risk(truth,pred,independent=False):
    if independent:return F(1,2)
    mismatch=sum(bit(truth,x)!=pred[x] for x in range(4))
    return ETA+(1-2*ETA)*F(mismatch,4)


def policy(learned):
    return min(range(4),key=lambda mask:(risk(learned,predictor(learned,mask))+COST*pop(mask),pop(mask),mask))


def binomial_cdf(n,p,k):
    return sum((F(math.comb(n,j))*p**j*(1-p)**(n-j) for j in range(k+1)),F(0))


def exhaustive_training_error():
    total=F(0)
    for errors in product((0,1),repeat=TRAIN_REPEATS):
        p=ETA**sum(errors)*(1-ETA)**(TRAIN_REPEATS-sum(errors))
        if sum(errors)>TRAIN_REPEATS//2:total+=p
    return total


def evaluate(truth,independent=False):
    # Each truth-table row has its own fresh noise preparation.
    q=F(1,2) if independent else exhaustive_training_error()
    paths=[]
    total=F(0);commit=F(0);pre_validation_commit=F(0)
    expected_loss=F(0);plain_loss=F(0);setup=F(36)*COST
    validation_cost=F(0);best_constant=F(1,2) if independent else risk(truth,predictor(truth,0))
    for learned in range(16):
        errors=pop(learned^truth)
        weight=q**errors*(1-q)**(4-errors)
        mask=policy(learned)
        p=risk(truth,predictor(learned,mask),independent)
        fallback=risk(truth,predictor(learned,0),independent)
        accepted=binomial_cdf(VALIDATE,p,MAX_ERRORS) if mask else F(0)
        loss_if_commit=p+COST*pop(mask)
        outcome_loss=accepted*loss_if_commit+(1-accepted)*fallback
        total+=weight
        commit+=weight*accepted
        pre_validation_commit+=weight*bool(mask)
        plain_loss+=weight*loss_if_commit
        expected_loss+=weight*outcome_loss
        validation_cost+=weight*VALIDATE*(pop(mask)+1)*COST*bool(mask)
        paths.append(dict(learned_table=learned,selected_ports=mask,
            training_probability=number(weight),true_prediction_error=number(p),
            validation_acceptance=number(accepted)))
    assert total==1
    gain=best_constant-expected_loss
    actual_net=SERVICE*gain-setup-validation_cost
    worst_setup=F(36+96)*COST
    conservative_net=SERVICE*gain-worst_setup
    return dict(truth=truth,independent_record=independent,
        normalization=number(total),commit_probability=number(commit),
        pre_validation_commit_probability=number(pre_validation_commit),
        oracle_constant_risk=number(best_constant),
        service_loss_with_port_fees=number(expected_loss),
        no_validation_service_loss=number(plain_loss),
        service_gain=number(gain),expected_setup_cost=number(setup+validation_cost),
        net_256_tasks=number(actual_net),net_using_max_setup=number(conservative_net),
        all_paths=paths)


def greedy(truth):
    current=0
    while True:
        old=risk(truth,predictor(truth,current))+COST*pop(current)
        candidates=[current|b for b in (1,2) if not current&b]
        improving=[s for s in candidates if risk(truth,predictor(truth,s))+COST*pop(s)<old]
        if not improving:return current
        current=min(improving,key=lambda s:(risk(truth,predictor(truth,s))+COST*pop(s),pop(s),s))


def parity_screen():
    rows=[]
    for r in range(2,7):
        # Keep c fixed; this verifies invisibility, not positive gain at every r.
        for size in range(r):
            for indices in combinations(range(r),size):
                counts={}
                for x in product((0,1),repeat=r):
                    key=tuple(x[i] for i in indices)
                    counts.setdefault(key,[0,0])[sum(x)%2]+=1
                assert all(a==b for a,b in counts.values())
        rows.append(dict(arity=r,all_proper_subsets_uninformative=True,
            full_prediction_error=number(ETA),
            gain_with_fixed_cost=number(F(1,2)-ETA-r*COST)))
    return rows


def reversible_query():
    # |f,x,n,z> -> |f,x,n,z XOR f(x) XOR n>; retains truth and noise internally.
    permutation=[]
    for f,x,n,z in product(range(16),range(4),range(2),range(2)):
        permutation.append((((f*4+x)*2+n)*2)+(z^bit(f,x)^n))
    a=np.array(permutation)
    assert np.array_equal(np.sort(a),np.arange(256))
    assert np.array_equal(a[a],np.arange(256))
    return dict(dimension=256,permutation_bijective=True,involution=True,
        preserves_truth_inputs_and_noise=True,
        whole_controller_autonomous_H_constructed=False)


def deterministic_record_witness():
    # A valid zero-noise basis history, not a Monte Carlo success claim.
    truth=6;records=[];observations=[[] for _ in range(4)];credits=900
    for x in range(4):
        for repeat in range(TRAIN_REPEATS):
            prior=observations[x]
            prediction=int(2*sum(prior)>len(prior)) if prior else 0
            # Prediction is appended before this history's target is exposed.
            row=dict(phase='training',ports=3,queried_values=x,prediction=prediction,
                     outcome=None,credits_before=credits)
            records.append(row)
            y=bit(truth,x);credits-=3;row['outcome']=y
            row['credits_after']=credits;observations[x].append(y)
    learned=sum(int(2*sum(v)>len(v))<<x for x,v in enumerate(observations))
    mask=policy(learned);assert mask==3
    failures=0
    for j in range(VALIDATE):
        x=(3*j+1)%4
        forecast=predictor(learned,mask)[x]
        row=dict(phase='validation',ports=mask,queried_values=x&mask,
                 prediction=forecast,outcome=None,credits_before=credits)
        records.append(row)
        y=bit(truth,x);failures+=forecast!=y;credits-=pop(mask)+1
        row['outcome']=y;row['credits_after']=credits
    assert failures==0 and credits==768 and len(records)==44
    return dict(truth_table=truth,learned_table=learned,selected_ports=mask,
        validation_errors=failures,remaining_credits=credits,
        additional_service_tasks_possible=SERVICE,
        retained_records=records,records_or_noise_erased=False,
        this_basis_history_probability_not_used_as_success_estimate=True)


def run():
    q=exhaustive_training_error()
    assert q==F(7,250)
    parity=6
    menu=[]
    for f in range(16):
        costs=[risk(f,predictor(f,s))+COST*pop(s) for s in range(4)]
        menu.append(dict(truth_table=f,oracle_greedy_ports=greedy(f),
            oracle_joint_ports=policy(f),costs=[number(v) for v in costs]))
    assert greedy(parity)==0 and policy(parity)==3
    result=[evaluate(f) for f in range(16)]
    independent=evaluate(0,True)
    xor=result[parity]
    assert F(xor['service_gain']['exact'])>F(1,10)
    assert F(xor['net_using_max_setup']['exact'])>9
    assert F(independent['commit_probability']['exact'])<F(1,40000)
    assert F(independent['service_gain']['exact'])<0
    assert F(independent['no_validation_service_loss']['exact'])==F(9,16)
    h_eta=-float(ETA)*math.log2(float(ETA))-(1-float(ETA))*math.log2(1-float(ETA))
    files=[Path(__file__),HERE/'drafts/mechanism_priority_entry.md',
        STAGE/'research_note_961.md',STAGE/'966/drafts/mechanism_map_v0_2.md',
        STAGE/'968/drafts/internal_organization_increment.md',
        STAGE/'972/drafts/thermal_record_mechanism.md',
        STAGE/'957/drafts/unified_operation_hypotheses_v0_2.md',
        STAGE/'987/drafts/common_overlap_audit_checks.json']
    return dict(round=987,date='2026-10-07',all_scientific_checks_passed=True,
        method='exact rational enumeration, analytic screening, finite permutation check',
        contract=dict(noise=number(ETA),port_cost=number(COST),training_repeats=3,
            training_rows=12,validation_rows=VALIDATE,maximum_validation_errors=MAX_ERRORS,
            service_tasks=SERVICE,maximum_port_uses=900,
            maximum_setup_cost=number(F(33,2)),cost_is_not_heat=True,
            fresh_input_noise_and_blank_records_are_preparation_inputs=True),
        analytic_screen=dict(single_port_information_bits=0,
            joint_information_bits=1-h_eta,oracle_greedy_stalls=True,
            exact_oracle_joint_gain=number(F(3,20)),general_parity=parity_screen()),
        training=dict(majority_error=number(q),
            all_16_functions_used=True,policies=menu),
        finite_trial=dict(all_functions=result,xor_summary={k:v for k,v in xor.items() if k!='all_paths'},
            independent_summary={k:v for k,v in independent.items() if k!='all_paths'},
            independent_false_acceptance_per_proposed_mask=number(binomial_cdf(32,F(1,2),5))),
        query=reversible_query(),witness=deterministic_record_witness(),
        mechanism_decision=dict(reject_strict_edgewise_improvement_as_universal=True,
            retain_bounded_coalition_trial_with_validation=True,
            unconditional_exploration_guarantees_benefit=False,
            fixed_max_coalition_size_is_universal=False,
            learned_functional_relation_is_spacetime=False,
            whole_goal_completed=False),
        references=['https://arxiv.org/abs/1004.2515','https://arxiv.org/abs/1207.1394'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-30),(a,b)
    else:assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    output=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(output,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(output,read(TARGET))
    print(json.dumps(dict(all_scientific_checks_passed=True,
        xor=output['finite_trial']['xor_summary'],
        independent=output['finite_trial']['independent_summary'],
        mechanism_decision=output['mechanism_decision']),ensure_ascii=False,indent=2))
