"""749: specific geometry/source erasure map fails the inherited nonlinear constraint."""
import argparse,hashlib,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round749_drafts'))
import joint_coarse_source_entry as first
import positive_constraint_defect as positive
TARGET=HERE/'joint_record_geometry_averaging_results.json'
def run():
    a=first.run();assert a==json.loads(first.TARGET.read_text('utf8'))
    b=positive.run();assert b==json.loads(positive.TARGET.read_text('utf8'))
    deps=('research_note_634.md','research_note_706.md','research_note_708.md',
          'research_note_724.md','research_note_725.md','research_note_731.md',
          'research_note_741.md','research_note_748.md','joint_source_constraint_response.py',
          'round749_drafts/joint_coarse_source_entry.py','round749_drafts/positive_constraint_defect.py')
    return dict(round=749,tests_run=2,failures=0,errors=0,
        checks=['original_all_source_family_and_mean_constraint_diagnostic',
                'positive_analytic_energy_profile_original_constraint_check'],
        original_family=a,positive_family=b,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Within731 original prescribed smooth source initial-data family, averaging conformal psi and canonical matter fields after erasing an equal-weight branch label is not closed under the Hamiltonian constraint. A positive analytic second-order coefficient proves the limited failure; numerical original-background checks calibrate it. This is not a derivation of Einstein dynamics, not a physical collapse claim, not a proof that actual quantum records realize these sources, and not a no-go for all coarse-graining maps. Old634 record/source covariance is reused.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=749,tests_run=2,specific_map_counterexample=True)))
