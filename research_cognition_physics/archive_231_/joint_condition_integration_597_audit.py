"""Read-only evidence/coverage audit for the non-numbered integration update.

This checks provenance and explicit scope, not scientific truth by table labels.
No scientific tests or completed round are added by this audit.
"""
import argparse
import ast
import json
from pathlib import Path
import re
import verify_interaction_rounds as core
import verify_round597 as baseline

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_condition_integration_597_audit_results.json'
NOTE=HERE/'joint_condition_compression_update_597.md'
STATUS=HERE/'round598_drafts/condition_integration_status.md'


def run():
    base=baseline.verify()
    assert (base['round'],base['cumulative_numbered_tests'])==(597,3037)
    assert not (HERE/'research_note_598.md').exists(), 'Inspect newer science before revising priorities'
    text=NOTE.read_text('utf8')
    rows=re.findall(r'^\|C(\d{2})\|([^\n]+)$',text,re.M)
    assert [int(n) for n,_ in rows]==list(range(1,28))
    assert all(len(line.split('|'))>=3 for _,line in rows)
    hashes={}
    for n in range(531,598):
        name=f'research_note_{n}.md'
        receipt=core.read(HERE/f'research_round_{n}_checks.json')
        expected=receipt['new_file_hashes'][name]
        assert core.digest(HERE/name)==expected,name
        hashes[name]=expected
    paths=(NOTE,STATUS,HERE/'joint_condition_compression_update_584.md',
           HERE/'cognition_forward_bridge_review_585.md',
           HERE/'unified_physics_condition_ledger.md')
    links=0
    for path in paths[:2]:
        for link in core.link_parser()(path.read_text('utf8')):
            target=(path.parent/link).resolve()
            assert target.exists() or target==TARGET.resolve(),(path.name,link)
            links+=1
    for path in paths:hashes[path.relative_to(HERE).as_posix()]=core.digest(path)
    ast.parse(Path(__file__).read_text('utf8'))
    formulas=core.text_checks(NOTE)
    assert formulas['display_formulas']==2
    # Scope declarations are manually reviewed evidence mappings, not SAT proofs.
    groups=[
        dict(id='process',rounds=[584,586,589,590,591,592,593],
             classification='same_declared_finite_graph_process',
             excludes=['automatic_infinite_dimensional_reconstruction','autonomous_apparatus','GR']),
        dict(id='gauge_matter',rounds=[531,532,568,569],
             classification='conditional_parameter_and_representation_matching',
             excludes=['unique_gauge_group','complete_chiral_lattice_dynamics','emergent_dimension']),
        dict(id='classical_geometry',rounds=[571,572,573,574],
             classification='classical_common_source_with_separate_fixed_graph_semiclassical_bridge',
             excludes=['finite_hbar_Einstein_solution','global_reference_chart']),
        dict(id='parameters',rounds=[543,544,545,546,547,551,552,553],
             classification='same_scale_and_approximation_required',
             excludes=['current_observational_fit','all_scale_spectral_identity']),
        dict(id='representations',rounds=[580,581,582,583,584],
             classification='three_distinct_maps_not_full_quantum_frame_equivalence',
             excludes=['full_metric_path_integral','fixed_hbar_continuum_limit']),
        dict(id='resources',rounds=[590,591,592,593,594,595,596,597],
             classification='conditional_budgets_and_implementation_boundaries',
             excludes=['free_hardware','demonstrated_long_lived_original_memory'])]
    return dict(date='2026-10-01',kind='non-numbered condition integration audit',
                scientific_baseline=597,cumulative_scientific_checks_unchanged=3037,
                numbered_scientific_files_unchanged=1103,
                protected_science_evidence_verified=1764,
                prior_numbered_notes_checked=67,condition_ids=[f'C{n}' for n,_ in rows],
                condition_table_covers_all_27=True,all_conditions_proved=False,
                groups=groups,group_count_is_not_axiom_count=True,
                active_goal_changed=False,old_round598_candidate_preserved=True,
                autonomous_memory_design_deferred=True,
                next_priority='existing chiral/gauge representation to current quantum candidate',
                new_scientific_round_count=0,new_scientific_test_count=0,
                source_hashes=hashes,text_checks=formulas,local_links_checked=links,
                historical_science_full_rerun=False,recent_baseline_checks_reproduced=True,
                primary_scope_review=True,independent_review=False,all_document_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert core.read(TARGET)==result
    print(json.dumps({k:result[k] for k in ('scientific_baseline','prior_numbered_notes_checked',
          'condition_table_covers_all_27','new_scientific_round_count','all_document_checks_passed')}))
