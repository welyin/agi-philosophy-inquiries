"""1032 local verification. External independent-agent review is incomplete."""
from pathlib import Path
import argparse
import ast
import hashlib
import json
import re
from urllib.parse import unquote
import electric_completion_selection as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
ROOT=BASE.parent
OUT=HERE/'research_round_1032_checks.json'
NOTE=HERE.parent/'research_note_1032.md'
OWN=['electric_completion_selection.py','electric_completion_selection_results.json',
     'drafts/electric_completion_derivation.md','selection_audit.md',
     'input_dependency_update_v0_21.md','NEXT.md','review.md','verify_round1032.py']
LINK=re.compile(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)')


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def local_links(text):
    text=re.sub(r'\$\$.*?\$\$|^```.*?^```\s*$','',text,flags=re.S|re.M)
    for m in LINK.finditer(text):
        v=m[1].strip();v=v[1:v.index('>')] if v.startswith('<') else v
        local=unquote(v.split('#',1)[0])
        if local and not re.match(r'^[a-zA-Z]+:',local):yield v,local


def verify(prospective=False):
    fresh=science.run()
    science.compare(fresh,json.loads(science.OUT.read_text('utf8')))
    assert fresh['code_sha256']==sha(HERE/'electric_completion_selection.py')
    assert fresh['external_independent_agent_review_completed'] is False
    assert fresh['new_adopted_cognitive_axioms']==0 and not fresh['goal_completed']
    assert fresh['cumulative_research_groups']==3809
    totals=[]
    for row in fresh['words']['cases']:
        p=row['p'];count=sm=0
        for a in range(5):
            for b in range(5):
                for n in range(5):
                    q0=-2*a+2*b+3*n
                    for q in range(-24,25):
                        if (q-q0)%p==0:
                            count+=1
                            if (q-q0)%6==0:sm+=1
        assert count==row['allowed']==row['completed_generated']
        assert sm==row['SM_generated']==1027
        totals.append(count)
    assert totals==[6125,3075,2045,1027] and sum(totals)==12272
    for row in fresh['center']['cases']:
        p=row['p'];gamma=row['quotient_subgroup_exponents']
        assert gamma==[k for k in range(6) if (k*p)%6==0]
        assert row['completed_kernel_exponents']==gamma
        assert row['sm_word_exists']==(p==6)
        assert row['residual_electric_center_order']==6//p
        assert row['anomaly_cubic']==row['anomaly_gravitational']==row['invariant_mass_charge']==0
    for row in fresh['closure']['cases']:
        assert row['allowed']==72//row['p']
        assert row['reached']==(72//row['p'] if row['extended'] else 12)
        assert row['generated_all']==(row['extended'] or row['p']==6)
    hw=fresh['highest_weights']
    assert len(hw['color_cases'])==13 and len(hw['weak_cases'])==7
    for row in hw['color_cases']:
        assert row['tensor_dimension']==3**(row['a']+row['b'])
        assert row['residual']<1e-12
    for row in hw['weak_cases']:
        assert row['tensor_dimension']==2**row['n'] and row['residual']<1e-12
    for name,digest in fresh['historical_source_sha256'].items():assert sha(BASE/name)==digest
    assets=[NOTE]+[HERE/name for name in OWN]
    links=0
    for p in assets:
        assert p.is_file(),p
        if p.suffix=='.py':ast.parse(p.read_text('utf8'))
        if p.suffix=='.md':
            text=p.read_text('utf8');assert text.count('$$')%2==0,p
            for target,local in local_links(text):
                q=(p.parent/local).resolve()
                assert q.exists() or (prospective and q==OUT),(p,target)
                links+=1
    assert all(f'## {i}.' in NOTE.read_text('utf8') for i in range(1,11))
    assert '独立代理签审' in (HERE/'review.md').read_text('utf8')
    return dict(round=1032,date='2026-10-08',local_delivery_checks_passed=True,
        scientific_result_reproduced=True,external_independent_agent_review_completed=False,
        general_theorem_supported_by_analytic_proof_not_grid=True,
        new_calibration_groups=1,cumulative_research_groups=3809,new_adopted_cognitive_axioms=0,
        allowed_grid_cases=12272,color_highest_weight_cases=13,weak_highest_weight_cases=7,
        finite_character_closures=8,closed_SM_menu_complete_iff_p_6=True,
        expanded_menu_complete_for_all_four_p=True,
        actual_endpoint_preparation_certified=False,uniform_resource_bound_certified=False,
        magnetic_completeness_certified=False,actual_low_energy_error_certified=False,
        full_quantum_gravity_certified=False,goal_completed=False,visual_checks_performed=False,
        frozen_current_files=9,historical_input_files=len(science.HISTORY),local_links_checked=links,
        live_navigation_frozen=False,
        historical_source_sha256=fresh['historical_source_sha256'],
        source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in assets})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-receipt',action='store_true')
    args=parser.parse_args();result=verify(prospective=args.write_receipt)
    if args.write_receipt:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:assert result==json.loads(OUT.read_text('utf8'))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith('sha256')},ensure_ascii=False))
