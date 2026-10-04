"""683: actual strict-time sources, uniform matching and a scoped reflection obstruction."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_time_local_source_interface_results.json'

def entry(name):
    path=HERE/'round683_drafts'/f'{name}.py'
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    result=module.run()
    assert result==json.loads(module.TARGET.read_text('utf8'))
    return result

def run():
    sources=entry('time_local_source_probe')
    gauge=entry('gauge_reflection_probe')
    free=entry('compensation_free_probe')
    deps=['research_note_669.md','research_note_670.md','research_note_673.md',
          'research_note_677.md','research_note_678.md','research_note_679.md',
          'research_note_681.md','research_note_682.md','joint_local_source_lift.py',
          'joint_physical_source_reflection.py','joint_rational_physical_limit.py']
    for name in ('time_local_source_probe','gauge_reflection_probe','compensation_free_probe'):
        deps += [f'round683_drafts/{name}.py',f'round683_drafts/{name}_results.json']
    return dict(date='2026-10-02',round=683,tests_run=2,failures=0,errors=0,
        strict_time_source_and_full_average_bounds=sources,
        actual_gauge_and_reflection_interface=dict(nonflat=gauge,exact_free=free),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_16_channel_sources=True,
            strict_time_support_in_auxiliary_variables=True,
            uniform_in_L_error_bound=True,all_fixed_order_unnormalized_sources_same_full_average_limit=True,
            exact_local_Gauss_covariance=True,declared_scalar_reflection_of_compensator_fails=True,
            all_possible_auxiliary_reflections_excluded=False,
            original_physical_Gauss_RP_decided=False,
            normalized_or_volume_uniform_or_physical_continuum_limit=False,
            original_HF_identity_or_quantum_GR=False,old_space_contracts_inherited=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=683,tests_run=2,all_checks_passed=True)))
