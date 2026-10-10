"""Exclusive author freeze; default verifies existing bytes and recalculates."""
from pathlib import Path
from hashlib import sha256
import argparse
import json
import re
import ckm_closure_check as science

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RECEIPT = HERE / 'research_round_1061_checks.json'
OWN = ['../research_note_1061.md','analysis.md','sources.md','public_inputs.json',
       'ckm_closure_check.py','results.json','dependency_update.md','verify_round1061.py']
HISTORY = ['archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
           'archive_956_989/981/drafts/common_parent_contract_v1.md',
           'archive_1046_/research_note_1060.md',
           'archive_1046_/1059/statistical_review.md',
           'archive_1046_/_admission/after1060_ckm_closure/selection.md',
           'archive_1046_/_admission/after1060_ckm_closure/sources.md']


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def run(freeze=False):
    result = science.calculate()
    assert result == science.calculate(100)
    assert result == json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    links = 0
    for name in OWN:
        path = HERE/name
        assert path.is_file(), name
        if path.suffix != '.md':
            continue
        for target in re.findall(r'\[[^\]\n]*\]\(([^)\n]+)\)',path.read_text(encoding='utf-8')):
            target = target.strip().strip('<>')
            if target.startswith(('http://','https://','#')):
                continue
            assert (path.parent/target.split('#',1)[0]).is_file(), (name,target)
            links += 1
    record = {
        'round':1061,'date':'2026-10-09','author_checks_passed':True,
        'owned_sha256':{name:digest(HERE/name) for name in OWN},
        'historical_sha256':{name:digest(ROOT/name) for name in HISTORY},
        'local_links_checked':links,
        'project_empirical_calibration_increment':1,'new_cognitive_axioms':0,
        'independent_review_complete_at_author_freeze':False,
        'main_four_summary_outer_region_nonempty':True,
        'full_shared_nuisance_matching_passed':False,
        'new_physics_discovery_or_CKM_generation_claimed':False,
        'whole_roadmap_complete':False,
    }
    if freeze:
        with RECEIPT.open('x',encoding='utf-8') as out:
            json.dump(record,out,ensure_ascii=False,indent=2)
            out.write('\n')
    else:
        assert record == json.loads(RECEIPT.read_text(encoding='utf-8'))
    return record


if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--freeze',action='store_true')
    args=p.parse_args();r=run(args.freeze)
    print(json.dumps({'round':1061,'author_checks_passed':True,'assets':len(OWN),
        'history':len(HISTORY),'links':r['local_links_checked'],
        'mode':'exclusive_freeze' if args.freeze else 'read_only'}))
