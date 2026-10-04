"""748: exact Duhamel reduction; numerical total coefficient is not certified."""
import argparse,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round748_drafts'))
import scalar_edge_response as entry
TARGET=HERE/'joint_remote_total_coefficient_results.json'
def run():
    result=entry.run()
    assert result==json.loads(entry.TARGET.read_text('utf8'))
    rows=result['quadrature_rows']
    assert abs(rows[-1]['total_sixth_for_common_epsilon_psi_hbar_1']
               -rows[-2]['total_sixth_for_common_epsilon_psi_hbar_1'])<1e-9
    deps=('research_note_577.md','research_note_746.md','research_note_747.md',
          'round745_drafts/exact_history_density_results.json',
          'round747_drafts/remote_hop_sign_certificate_results.json',
          'round748_drafts/scalar_edge_response.py')
    return dict(round=748,tests_run=2,failures=0,errors=0,
        checks=['native_geodesic_force_and_exact_Duhamel_time_factor',
                'joint_original_channel_quadrature_convergence_diagnostic'],
        result=result,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Analytic reduction of the original scalar-edge sixth coefficient to the inherited local fourth density and the original edge force, then addition to747 original hopping contribution. Numerical evidence at one declared common geometry gives same-sign channels. Convergence is not a rigorous integration error certificate; total remote distinguishability remains unproven. No coupling, probe, record reset or geometric reconstruction has been added. Return to the unified condition ledger before further specialized readout optimization.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=748,tests_run=2,total_signal_certified=False)))
