"""Append-only mainline provenance/dependency verification for round 1042."""
from pathlib import Path
import argparse
import json
import math
import runpy

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
BASE=STAGE.parent
OUT=HERE/'mainline_acceptance.json'
PRIOR=STAGE/'1041/mainline_acceptance.json'
OLD=runpy.run_path(str(PRIOR.parent/'verify_mainline_acceptance.py'))
sha,read,relative,check_map=[OLD[n] for n in ('sha','read','relative','check_map')]
REVIEWS=[HERE/'independent_mathematical_review.md',HERE/'independent_scope_review.md']
ADMISSIONS=[STAGE/'1041'/n for n in ('relay_information_admission.md',
    'species_capacity_admission.md','thermodynamic_selection_admission.md','weyl_matching_admission.md')]
NAV=[BASE/'README.md',BASE/'research_direction.md',BASE/'RESEARCH_STATE.md',
     STAGE/'README.md',STAGE/'文件索引.md']


def evidence():
    assert sha(PRIOR)=='589badc27ac57a0ae6a87c837f303728ccefc2d83085a31f0759bae8c591c32d'
    old=read(PRIOR); old.pop('navigation_links_at_acceptance')
    assert old==OLD['evidence'](), 'Prior evidence changed'
    author_path=HERE/'research_round_1042_checks.json'
    author=read(author_path)
    assert author['passed'] and author['author_science_replay_passed']
    assert author['scientific_groups_proposed']==1
    science=check_map(HERE,author['owned_sha256'])
    history=check_map(BASE,author['historical_sha256'])
    delta_path=HERE/'dependency_delta_1042.json'
    delta=read(delta_path)
    prior_delta_path=BASE/delta['base_delta']
    assert sha(prior_delta_path)==delta['base_delta_sha256']
    d1=read(prior_delta_path)
    d0=read(BASE/d1['base_delta'])
    rows=read(BASE/d0['base_ledger'])['rows']+d0['rows_append']+d1['rows_append']
    assert len(rows)==delta['base_resolved_rows']==33
    rows+=delta['rows_append']
    assert len(rows)==delta['resolved_row_count']==34
    assert len({r['id'] for r in rows})==34
    for row in rows:
        for source in row['source_paths']: assert (BASE/source).exists(),source
    for row in delta['rows_append']:
        assert row['requires'] and row['not_proved']
        for flag in ('cognitive_origin_closed','full_parent_process_certified','physical_implementation_certified'):
            assert row[flag] is False
    assert not delta['goal_completed'] and delta['new_adopted_cognitive_axioms']==0
    for p in REVIEWS:
        body=p.read_text('utf8')
        assert '通过' in body and '1042' in body
        assert author['owned_sha256']['../research_note_1042.md'] in body
    data=read(HERE/'complementary_entropy_results.json')
    assert data['families'][0]['p']==[0.5,0.25]
    assert data['families'][2]['p']==[0.125,0.375]
    quotient_gap=max(data['families'][0]['central_entropy'])-max(data['families'][2]['central_entropy'])
    assert 3**3*5**5>4**8
    assert abs(quotient_gap-math.log((3**3*5**5)/(4**8))/8)<1e-13
    supplements=[Path(__file__),delta_path,HERE/'mainline_integration.md',
                 HERE/'goal_scope_review.md',HERE/'representation_quotient_clarification.md',
                 *REVIEWS,*ADMISSIONS]
    return dict(schema='append_only_mainline_acceptance_v1',date='2026-10-08',
                rounds_accepted=[1042],cumulative_before=3818,cumulative_after=3819,
                new_scientific_calibration_groups=1,integration_calibration_groups=0,
                new_adopted_cognitive_axioms=0,prior_acceptance=relative(PRIOR),
                prior_acceptance_sha256=sha(PRIOR),author_receipt=relative(author_path),
                author_receipt_sha256=sha(author_path),
                independent_agent_reviews=[relative(p) for p in REVIEWS],human_peer_review=False,
                historical_author_receipts_preserved=True,dependency_delta=relative(delta_path),
                dependency_base_rows=33,dependency_added_rows=1,dependency_resolved_rows=34,
                full_common_parent_model_established=False,physical_implementation_certified=False,
                goal_completed=False,additional_frozen_science_sha256=science,
                quotient_freedom_witness='original_first_and_third_encoding_families',
                largest_central_entropy_eigenvalue_gap=quotient_gap,
                additional_historical_sha256=history,
                supplemental_sha256={relative(p):sha(p) for p in supplements})


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=evidence()
    checker=OLD['OLD']['check_links']
    checker.__globals__['OUT']=OUT
    links=sum(checker(p,args.write) for p in NAV)
    for p in [HERE/'mainline_integration.md',HERE/'goal_scope_review.md',
              HERE/'representation_quotient_clarification.md',*REVIEWS,*ADMISSIONS]:
        checker(p,args.write)
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(dict(result,navigation_links_at_acceptance=links),f,
                      ensure_ascii=False,indent=2,allow_nan=False); f.write('\n')
    else:
        old=read(OUT); old.pop('navigation_links_at_acceptance')
        assert old==result,'Frozen mainline evidence changed'
    print(json.dumps(dict(passed=True,round=1042,cumulative=3819,dependency_rows=34,
                         author_assets=len(result['additional_frozen_science_sha256']),
                         agent_reviews=2,live_navigation_links=links)))
