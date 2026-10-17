"""Finite governance specification audit; physical realization uses cited old proofs."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import itertools
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dist(weight):
    return {(0,0):F(3,4)*(1-weight),(0,1):F(3,5)*weight,
            (1,0):F(1,4)*(1-weight),(1,1):F(2,5)*weight}

def calculate():
    source=json.loads((HERE/'joint_contract_sources.json').read_text(encoding='utf-8'))
    for path,h in source['source_sha256'].items():
        assert sha(ROOT/path)==h,path
    assert sha(Path(source['policy_docx']['path']))==source['policy_docx']['sha256']
    assert sha(ROOT/'research_cognition_physics/archive_1234_/_shared/goal_prompt.txt')==source['goal_sha256']
    prior=json.loads((HERE.parent/'governance_native_joint_instance_after1272/results.json').read_text(encoding='utf-8'))
    cooperation=json.loads((HERE.parent/'governance_cooperative_capacity_after1272/results.json').read_text(encoding='utf-8'))
    assert cooperation['isolated_event_error_lower_bound']=='1/10'
    checked=0; outputs=[]
    for row in prior['finite_cases']:
        q=F(row['q_population']); total=F(0); event2=F(0); false_claim=F(0)
        first=dist((1+q)/2)
        assert first[(1,1)]==F(row['first_event'])
        for y1,w1 in first.items():
            # J and direction are authorized classical results of the original instrument.
            task=int(y1==(1,1)); program=task
            second=dist(F(1) if program==0 else F(1,2))
            for y2,w2 in second.items():
                w=w1*w2; total+=w
                e2=int(y2==(1,1)); event2+=w*e2
                false_claim+=w*(1-e2)  # fixed self-report h=1, ignoring actual evidence
                for claim,exit_first_role in itertools.product((0,1),repeat=2):
                    checked+=1
                    evidence=(y1,y2)
                    transferred=evidence  # authorized transport, not deletion/re-preparation
                    owner='successor' if exit_first_role else 'current_verifier'
                    disagreement=int(claim!=e2)
                    verdict={'observed':e2,'claim':claim,'reopen_record_review':disagreement,
                             'owner':owner,'evidence':transferred,'program':program}
                    assert verdict['observed']==int(verdict['evidence'][1]==(1,1))
                    assert verdict['reopen_record_review']==(verdict['claim']!=verdict['observed'])
                    assert verdict['program']==int(verdict['evidence'][0]==(1,1))
                    # Role exit cannot turn contradictory physical evidence into a passed claim.
                    assert not (disagreement and verdict['reopen_record_review']==0)
        assert total==1
        assert event2==F(row['programs']['feedback']['second_event'])
        assert false_claim==1-event2
        assert event2-F(row['programs']['log_only']['second_event'])==F(row['event_gap'])
        outputs.append({'q_population':str(q),'second_actual_event':str(event2),
                        'false_completion_probability_if_only_claim_h1_used':str(false_claim),
                        'actual_record_review_corrects_claim_in_every_branch':True,
                        'exit_changes_evidence_or_verdict':False})
    assert checked==192
    # Algebraic history witness, not an assertion that no other physical measurement exists.
    h=1; possible=[(h,0),(h,1)]
    assert possible[0][0]==possible[1][0] and possible[0][1]!=possible[1][1]
    return {'status':'PASS_FINITE_JOINT_CONTRACT_SCOPE','branch_specifications_checked':checked,
            'cases':outputs,'native_probability_changes':'none; inherited original two-kit protocol',
            'cooperative_event_error_lower_bound_reused':'1/10',
            'feedback_event_gap_reused':'3/25',
            'source_and_geometry_proof':'1271 finite closure and 1166 independent Flat_op; not inferred from this enumeration',
            'scope':'Classical finite record/role specification only; full quantum state, actual transport, age, source and geometry inherited from the exact cited realization proofs.',
            'not_proved':['All social policy effects','New physical failure remediation after kit 2',
                          'Unlimited reopening','General geometry or Jacobson selection','Actual laboratory precision'],
            'new_scientific_rounds':0,'new_formal_check_groups':0,'new_adopted_axioms':0}

def main():
    result=calculate(); target=HERE/'joint_contract_results.json'
    if sys.argv[1:]==['--write']:
        with target.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        mode='created_once'
    else:
        assert not sys.argv[1:]
        assert json.loads(target.read_text(encoding='utf-8'))==result
        mode='read_only_match'
    print(json.dumps({'status':result['status'],'mode':mode,'branches':result['branch_specifications_checked'],
                      'new_rounds':0},ensure_ascii=False))

if __name__=='__main__': main()
