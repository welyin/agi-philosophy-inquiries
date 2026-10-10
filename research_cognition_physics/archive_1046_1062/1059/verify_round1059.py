"""1059 author receipt: explicit exclusive freeze, otherwise read-only."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import qed_recoil_check as science

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RECEIPT=HERE/'research_round_1059_checks.json'
OWN=['../research_note_1059.md','analysis.md','public_inputs.json','qed_recoil_check.py',
     'results.json','sources.md','dependency_update.md','NEXT.md','verify_round1059.py']
HISTORY=['archive_956_989/981/drafts/common_parent_contract_v1.md',
         'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
         'archive_1046_/research_note_1049.md','archive_1046_/research_note_1052.md',
         'archive_1046_/_admission/after1058_qed_prediction/selection.md']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run(freeze=False):
    actual=science.calculate()
    assert actual==science.calculate(90)
    assert actual==json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    assert actual['all_checks_passed']
    local_links=0
    for name in OWN:
        path=HERE/name
        assert path.is_file(), name
        if path.suffix!='.md':
            continue
        body=path.read_text(encoding='utf-8')
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',body):
            target=target.strip().strip('<>')
            if target.startswith(('http://','https://','#')):
                continue
            assert (path.parent/target.split('#',1)[0]).is_file(),(name,target)
            local_links+=1
    record={
        'round':1059,'date':'2026-10-09','author_checks_passed':True,
        'owned_sha256':{name:digest(HERE/name) for name in OWN},
        'historical_sha256':{name:digest(ROOT/name) for name in HISTORY},
        'local_links_checked':local_links,
        'independent_review_complete_at_author_freeze':False,
        'project_empirical_calibration_increment':1,'new_cognitive_axioms':0,
        'new_QED_theorem_claimed':False,'whole_roadmap_complete':False,
        'scope':'fixed-version retrospective shared-alpha calibration and free-electron anomaly test',
        'main_rejection_uses_only_two_marginal_coverage_models':True,
        'full_QED_and_experimental_budget_rebuilt':False,
        'two_ae_residuals_combined_as_independent':False
    }
    if freeze:
        with RECEIPT.open('x',encoding='utf-8') as out:
            json.dump(record,out,ensure_ascii=False,indent=2);out.write('\n')
    else:
        assert record==json.loads(RECEIPT.read_text(encoding='utf-8'))
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true');args=p.parse_args()
    r=run(args.freeze)
    print(json.dumps({'round':1059,'author_checks_passed':True,'owned_assets':len(OWN),
                      'historical_inputs':len(HISTORY),'local_links':r['local_links_checked'],
                      'mode':'exclusive_freeze' if args.freeze else 'read_only',
                      'independent_review':'pending_at_author_freeze'}))
