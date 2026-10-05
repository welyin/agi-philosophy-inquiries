"""Aggregate round 802 coefficient checks and previously uncounted diagnostics."""
from pathlib import Path
import argparse
import json
import original_bff_vertex as vertex
import reference_vertex_probe as reference
import dual_record_response_probe as joint

HERE=Path(__file__).resolve().parent
RESULT=HERE/'physical_bff_response_results.json'

def run():
    actual=vertex.run()
    assert actual==json.loads((HERE/'original_bff_vertex_results.json').read_text('utf-8'))
    ref=reference.run();dual=joint.run()
    assert ref==json.loads((HERE/'reference_vertex_probe_results.json').read_text('utf-8'))
    assert dual==json.loads((HERE/'dual_record_response_probe_results.json').read_text('utf-8'))
    return dict(round=802,all_checks_passed=True,fresh_test_groups=2,
        original_coefficient=actual,reference_slice=ref,joint_output=dual,
        full_local_BFF_coefficient_mapped_to_original_slice=True,
        on_shell_fermion_cross_density_vanishes=True,
        nonzero_original_continuum_record_response_proven=False,
        graph_continuum_or_autonomous_apparatus_proven=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not RESULT.exists(),'Do not overwrite evidence.'
        RESULT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(RESULT.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
