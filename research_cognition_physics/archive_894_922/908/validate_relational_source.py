"""908 fast reproduction plus saved full refinement and physical-family evidence."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import relational_source_checks as source
import physical_family_source_checks as family
from material_background import Background
HERE=Path(__file__).resolve().parent;TARGET=HERE/'relational_source_validation_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(rerun=False):
    raw=json.loads(source.TARGET.read_text('utf-8'));physical=json.loads(family.TARGET.read_text('utf-8'))
    for name,digest in raw['source_file_sha256'].items():assert sha(HERE/name)==digest,name
    if rerun:assert source.run()==raw and family.run()==physical
    fresh,q=source.row(Background(16),2,8)
    saved=raw['rows'][0]
    assert fresh==saved
    jets=source.jet_checks(q)
    assert len(raw['rows'])==6 and len(physical['rows'])==2
    assert raw['refinement_comparisons']['path16_to32']<raw['refinement_comparisons']['path8_to16']
    assert all(abs(r['full_source']['total'][2])>1e-8 for r in physical['rows'])
    names=('magnetic_reference_source.py','material_background.py','relational_loop_source.py','relational_source_checks.py',
        'relational_source_checks_results.json','physical_family_source_checks.py','physical_family_source_checks_results.json',
        'drafts/initial_anchor_failure.json','drafts/source_implementation_working.md')
    return dict(round=908,date='2026-10-06',all_checks_passed=True,cumulative_numbered_groups=3693,
        argument_scope=raw['argument_scope'],one_actual_background_and_loop_case_reproduced=True,
        independent_offshell_source_checks_reproduced=jets,
        full_refinement_and_family_runs_previously_executed=True,all_large_runs_rerun_by_default=False,
        full_check_command='Set OPENBLAS_NUM_THREADS=1; python -B -X utf8 validate_relational_source.py --rerun-all',
        refinement_diagnostics=raw['refinement_comparisons'],physical_family_smallest_epsilon=physical['rows'][-1],
        source_file_sha256={name:sha(HERE/name) for name in names},
        original_first_jet_weak_source_implemented=True,strong_spacetime_source_array_computed=False,
        full_support_and_dynamical_error_certified=False,original872_forced_feedback_computed=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');p.add_argument('--rerun-all',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run(a.rerun_all)
    if a.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_file_sha256','physical_family_smallest_epsilon')},ensure_ascii=False,indent=2))
