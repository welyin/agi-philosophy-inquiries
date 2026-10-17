"""Read-only source and arithmetic check for an old finite-closure corollary."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def calculate():
    b=json.loads((HERE/'source_baseline.json').read_text(encoding='utf-8'))
    for path, h in b['source_sha256'].items():
        assert digest(ROOT/path)==h,path
    assert digest(Path(b['policy_docx']['path']))==b['policy_docx']['sha256']
    assert digest(ROOT/'research_cognition_physics/archive_1234_/_shared/goal_prompt.txt')==b['goal_sha256']
    old=json.loads((HERE.parent/'governance_native_joint_instance_after1272/results.json').read_text(encoding='utf-8'))
    old={F(row['q_population']):row for row in old['finite_cases']}
    p0,p1=F(old[F(0)]['first_event']),F(old[F(1)]['first_event'])
    lower=abs(p1-p0)/2
    assert (p0,p1,lower)==(F(1,5),F(2,5),F(1,10))
    midpoint=(p0+p1)/2
    assert abs(midpoint-p0)==abs(midpoint-p1)==lower
    return {'status':'PASS_OLD_COOPERATIVE_NATIVE_COROLLARY',
            'original_native_probabilities':[str(p0),str(p1)],
            'isolated_event_error_lower_bound':str(lower),
            'arithmetic_midpoint':str(midpoint),
            'joint_exact_model_event_error':'0',
            'basis':'Equal local input states imply equal local outputs; same old physical closure permits complete decoding and native task.',
            'not_certified':['New physical primitive','Optimal cost or deadline','Full strong governance conjecture','New geometric selection','Empirical implementation'],
            'new_scientific_rounds':0,'new_formal_check_groups':0,'new_adopted_axioms':0}

def main():
    result=calculate(); output=HERE/'results.json'
    if sys.argv[1:]==['--write']:
        with output.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2); stream.write('\n')
        mode='created_once'
    else:
        assert not sys.argv[1:]
        assert result==json.loads(output.read_text(encoding='utf-8'))
        mode='read_only_match'
    print(json.dumps({'status':result['status'],'mode':mode,'lower_bound':result['isolated_event_error_lower_bound']}))

if __name__=='__main__': main()
