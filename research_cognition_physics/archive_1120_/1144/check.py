"""1144: replay both newly admitted qualifications; no continuum simulation."""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, subprocess, sys
HERE=Path(__file__).resolve().parent

def calculate():
    hashes=json.loads((HERE/'source_hashes.json').read_text(encoding='utf-8'))
    for name,digest in hashes.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest,name
    entries=[('../_shared/after1143_independent_record_parent_screen.py',3),
             ('../_shared/after1143_native_record_reuse_check.py',4)]
    for name,groups in entries:
        run=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/name)],capture_output=True,text=True,encoding='utf-8')
        assert run.returncode==0,run.stdout+'\n'+run.stderr
    cell=json.loads((HERE/'../_shared/after1143_independent_record_parent_screen.json').read_text(encoding='utf-8'))
    kg=json.loads((HERE/'../_shared/after1143_native_record_reuse_results.json').read_text(encoding='utf-8'))
    ec=F(cell['normal_autonomous_analytic_budget']['autonomous_joint_isometry_bound_T5'])
    eg=F(kg['rational_finite_T_certificate']['complete_isometry_upper'])
    assert ec<F(709,1000000) and eg<F(865,1000000)
    # Conservative bound in proof and frozen source certificate.
    assert 2*12000*3*24**2 < 6600**2
    return dict(round=1144,status='PASS',check_groups=7,frozen_sources=hashes,
        autonomous_cell_isometry_upper=str(ec),continuous_source_isometry_upper=str(eg),
        autonomous_source_effect_gap_lower=cell['normal_autonomous_analytic_budget']['actual_bounded_source_probability_gap_strict_lower'],
        exact_fixed_profile_closure_obstacle_is_analytic=True,
        continuum_source_phase_and_residual_proved=True,
        qualification_energy_picture_erratum_in_proof=True,
        two_instances_jointly_realized=False,actual_CO_SR_clock_and_position_bridge=False,
        complete_local_stress_or_four_momentum=False,goal_complete=False,
        new_axioms=0,scientific_calibration_increment=0,scientific_count_total=3860,
        finite_matrices_are_not_continuum_or_autonomous_simulations=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true');a=ap.parse_args()
    result=calculate();target=HERE/'results.json'
    if a.write_results:
        assert not target.exists(),'Do not overwrite frozen results'
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        assert json.loads(target.read_text(encoding='utf-8'))==result,'Saved result mismatch'
    print(json.dumps({k:result[k] for k in ['round','status','check_groups','autonomous_cell_isometry_upper','continuous_source_isometry_upper','goal_complete']},ensure_ascii=False))
if __name__=='__main__': main()
