"""Reproducible integration checks, not a new experiment or research round.

Checks the common SM(+3N) normalization, branch restrictions, source certificate,
immutable round evidence and the reviewed paper/condition ledger. No UV or
instrument construction is inferred from these finite algebraic checks.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import re
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
BASE = STAGE.parent
ROOT = BASE.parent
OUT = HERE / 'conditional_sm_joint_checks.json'
PAPER = BASE / '13_标准模型结构的条件重建与共同来源_阶段论文.md'

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def compute():
    # Field data are explicit structural inputs, not inferred by this script.
    y = dict(Q=F(1,6), uc=F(-2,3), dc=F(1,3), L=F(-1,2), ec=F(1), N=F(0), H=F(1,2))
    allowed = [y['Q']+y['H']+y['uc'], y['Q']-y['H']+y['dc'],
               y['L']-y['H']+y['ec'], y['L']+y['H']+y['N']]
    assert allowed == [0]*4
    # All left Weyl fields have Z4 charge 1, either H or H* has charge 2.
    assert (1+1+2)%4 == 0 and (1+1)%4 == 2 and (1+2+1+2)%4 == 2
    # Independent direct Hessian of |D H|^2 at H=(0,1), hence v=sqrt(2).
    p1 = np.array([[0,1],[1,0]], complex)
    p2 = np.array([[0,-1j],[1j,0]], complex)
    p3 = np.diag([1.,-1.]).astype(complex)
    g1, g2 = F(1,2), F(2,3)
    h = np.array([0.,1.], complex)
    gens = [float(g2)*p/2 for p in (p1,p2,p3)] + [float(g1)*np.eye(2)/2]
    mass2 = np.array([[2*np.vdot(a@h,b@h).real for b in gens] for a in gens])
    expected = [0., float(F(2,9)), float(F(2,9)), float(F(25,72))]
    residual = float(np.max(np.abs(np.linalg.eigvalsh(mass2)-expected)))
    assert residual < 1e-14
    photon = np.array([0.,0.,float(g1),float(g2)])
    assert np.max(np.abs(mass2@photon)) < 1e-14
    # Bounds on the inherited family; no numerical joint (s,kappa) is claimed.
    mu = F(1,4096)
    lower = min(1+F(1,4)*(j*mu-1) for j in (1,2,3))
    assert lower > F(3,4)
    source = json.loads((STAGE/'1307/results.json').read_text(encoding='utf-8'))
    independent = json.loads((STAGE/'1307/independent_results.json').read_text(encoding='utf-8'))
    assert source['status'] == 'PASS'
    assert independent['status'] == 'PASS_INDEPENDENT_LEADING_SOURCE_CLASS'
    assert [x['all_dynamic_weight_dimension'] for x in source['exact_source_kernels']] == [10,5,2,1]
    assert [x['Hermitian_extra_source_dimension'] for x in independent['exact_kernel_dimension_table']] == [9,4,1,0]
    assert independent['exact_complex_dictionary_M_AR_equals_conjugate_YA']
    assert all(independent['Hilbert_momentum_commutator_identity_each_spatial_axis'])
    frozen_count = 0
    for n in range(1290,1308):
        ap = STAGE/str(n)/'acceptance.json'
        a = json.loads(ap.read_text(encoding='utf-8'))
        assert a['status'].startswith('ACCEPTED'), (n,a['status'])
        for name, expected_hash in a['files_sha256'].items():
            p = ap.parent/name
            assert digest(p)==expected_hash, str(p)
            frozen_count += 1
    support = json.loads((BASE/'_migration/structure_selection_active_1288/layer.json').read_text(encoding='utf-8'))
    for name, expected_hash in support['unchanged_support_sha256'].items():
        assert digest(ROOT/name)==expected_hash, name
    files = [PAPER, HERE/'goal_after1306.md', HERE/'goal_after1306_scope_review.md',
             HERE/'conditional_sm_premise_ledger.md', HERE/'conditional_sm_paper_math_review.md',
             HERE/'conditional_sm_paper_scope_review.md', STAGE/'research_note_1307.md',
             STAGE/'1307/acceptance.json', STAGE/'1307/proof.md', STAGE/'1305/acceptance.json']
    local_links = 0
    for p in [PAPER,HERE/'conditional_sm_premise_ledger.md']:
        text = p.read_text(encoding='utf-8')
        assert text.count('$$')%2 == 0, 'Unbalanced display math: '+str(p)
        assert '\\[' not in text and '\\]' not in text, 'Old math delimiters'
        for link in re.findall(r'\]\(([^)]+)\)', text):
            if '://' in link or link.startswith('#'):
                continue
            target = (p.parent/link.split('#')[0]).resolve()
            # Acceptance is written immediately after this first exclusive check.
            if target == HERE/'conditional_sm_stage_acceptance.json' and not target.exists():
                continue
            assert target.exists(), (str(p),link)
            local_links += 1
    assert digest(HERE/'goal.md') == '4bd9fa800ac8385b0850e13526f4f5356b62a287464bb13b1cacc1a22370072e'
    return dict(schema='conditional_sm_integration_checks_v1', status='PASS_WITH_SCOPE',
                scientific_round_increment=0, new_axioms=0, accepted_rounds=18,
                frozen_scientific_bindings_verified=frozen_count,
                historical_support_bindings_verified=len(support['unchanged_support_sha256']),
                paper_and_ledger_local_links=local_links,
                input_branch='G6 standard representations, 3 generations, single Higgs, SM+3N, exact r, Dirac dimension<=4',
                yukawa_hypercharge_sums=[str(v) for v in allowed],
                z4=dict(yukawa=0,majorana=2,weinberg=2),
                gauge_hessian_eigenvalues_squared=[str(x) for x in [F(0),F(2,9),F(2,9),F(25,72)]],
                gauge_hessian_tolerance=1e-14,
                inherited_mass_lower_bound=str(lower),
                source_rank_dimensions=[10,5,2,1],
                meaning='Integration audit only. Inputs stay inputs; background certificate is not a full interacting QFT error estimate.',
                files_sha256={str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in files})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=compute()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
    else:
        old=json.loads(OUT.read_text(encoding='utf-8'))
        # Acceptance creation adds two existing links after the first write.
        result.pop('paper_and_ledger_local_links',None)
        old.pop('paper_and_ledger_local_links',None)
        assert result==old, 'Saved integration certificate mismatch'
    print(json.dumps(dict(status='PASS_WITH_SCOPE',accepted_rounds=18, new_scientific_groups=0,
                         frozen_bindings=result['frozen_scientific_bindings_verified']),ensure_ascii=False))

if __name__=='__main__':
    main()
