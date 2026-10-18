"""Root acceptance of the public conditional goal, after independent reviews."""
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
STAGE = HERE.parent
OUT = HERE / 'conditional_sm_stage_acceptance.json'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--integration-receipt', required=True)
    args = parser.parse_args()
    if OUT.exists():
        raise SystemExit('Acceptance already exists; no overwrite.')
    checks = json.loads((HERE/'conditional_sm_joint_checks.json').read_text(encoding='utf-8'))
    assert checks['status'] == 'PASS_WITH_SCOPE'
    assert checks['accepted_rounds'] == 18
    for rel, expected in checks['files_sha256'].items():
        assert sha(ROOT/rel) == expected, rel
    review_bindings = {}
    for name in ('conditional_sm_paper_math_review.md', 'conditional_sm_paper_scope_review.md'):
        path = HERE/name
        text = path.read_text(encoding='utf-8')
        assert 'PASS' in text, name
        count = 0
        for line in text.splitlines():
            if not line.startswith('|'):
                continue
            cells = [v.strip().strip('`') for v in line.strip('|').split('|')]
            if len(cells) != 2 or not re.fullmatch('[0-9a-f]{64}', cells[1]):
                continue
            link = re.fullmatch(r'\[[^\]]+\]\(([^)]+)\)', cells[0])
            target = (HERE/(link.group(1) if link else cells[0])).resolve()
            assert target.is_file() and sha(target) == cells[1], (name, cells[0])
            count += 1
        assert count >= 3, (name, 'No sufficient version bindings')
        review_bindings[name] = count
    files = dict(checks['files_sha256'])
    for name in ('conditional_sm_joint_check.py', 'conditional_sm_joint_checks.json',
                 'accept_conditional_sm_stage.py', 'finalize_conditional_sm_navigation.py'):
        path = HERE/name
        files[path.relative_to(ROOT).as_posix()] = sha(path)
    result = dict(
        schema='conditional_sm_stage_acceptance_v1', date='2026-10-18',
        status='COMPLETE_WITH_SCOPE',
        goal='Standard Model conditional reconstruction and common consistency; goal_after1306.md',
        scope='SM major structure with explicit inputs; SM+3N protected Dirac common branch. Full parent structure and its declared background quantum process are distinct levels with a checked restriction map.',
        criteria={
            'A': {'status':'PASS', 'evidence':['conditional_sm_premise_ledger.md', 'goal_after1306_scope_review.md']},
            'B': {'status':'PASS', 'evidence':['../../13_标准模型结构的条件重建与共同来源_阶段论文.md', '../1307/acceptance.json']},
            'C': {'status':'PASS', 'evidence':['../1307/proof.md', '../1305/acceptance.json', 'conditional_sm_joint_checks.json']},
            'D': {'status':'PASS', 'evidence':['conditional_sm_paper_math_review.md', 'conditional_sm_paper_scope_review.md', 'conditional_sm_stage_acceptance.json']}
        },
        scientific_rounds=list(range(1290,1308)), author_check_groups=18, independent_check_groups=18,
        new_scientific_rounds_this_goal=[1307], synthesis_scientific_round_increment=0,
        new_axioms=0,
        conditions_retained=[
            'Standard G6 group/global quotient, standard representations and three generations',
            'One Higgs, healthy local field class and dimension<=4 leading action',
            'Exact discrete r protection and stable nonzero Higgs vacuum branch',
            'Coupling, potential and Yukawa parameters; full-column-rank Dirac Ynu for common neutral source weight',
            'Declared leading source class and equivalence; total normalization, constants, improvements and higher-order choices remain',
            'Specified backgrounds, reference preparations, finite access and separate error/resource contracts'
        ],
        conclusions=[
            'Gauge and matter vertices obey the stated consistency relations',
            'One Higgs/Yukawa structure supplies tree masses and mixing relations',
            'Declared fermionic anomaly classes vanish on the stated joint backgrounds',
            'Ynu(R-kI)=0 is necessary and sufficient in the leading neutral source class; kernel dimension d leaves d^2 real directions',
            'The same parent has the declared finite background quantum record/source/thermal-geometry certificate'
        ],
        limits=[
            'No unique derivation of group, representations, generations or numerical parameters',
            'The common branch extends the strictly minimal Standard Model by three right-handed neutrinos',
            'No full interacting gauge/Higgs quantum dynamics or controlled approximation to all such dynamics',
            'No unrestricted UV completion, empirical fit, full hadron spectrum or theory-of-everything claim'
        ],
        review_version_bindings_verified=review_bindings,
        frozen_scientific_bindings_verified=checks['frozen_scientific_bindings_verified'],
        historical_support_bindings_verified=checks['historical_support_bindings_verified'],
        actual_root_execution_receipts={'author_1307':'0e955f','independent_1307':'367195',
                                       'round_1307_acceptance':'85a694','integration':args.integration_receipt},
        application_goal_status_at_acceptance='active; complete only after final navigation checks',
        next_scope='No further scientific round opened. Stronger input generation, empirical parameters or full dynamic matching require a separately scoped goal.',
        files_sha256=files)
    with OUT.open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(result,stream,ensure_ascii=False,indent=2)
        stream.write('\n')
    print(json.dumps({'status':result['status'],'reviews':review_bindings,'sha256':sha(OUT)},ensure_ascii=False))

if __name__ == '__main__':
    main()
