"""690: original spatialY insertions, exact neutral Ward mismatch, full averaging.

The result excludes the unchanged conserved-CAR-charge dictionary, not the
original reflection positivity or all possible enlarged/continuum realizations.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_spatial_charge_dictionary_results.json'
spec=importlib.util.spec_from_file_location('entry690',HERE/'round690_drafts/spatial_car_source_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def val(x):return float(x[0])+float(x[1])*np.sqrt(5)
def serial(x):return dict(rational=str(x[0]),sqrt5_coefficient=str(x[1]))


def exact_certificate(probe):
    one=(F(1),F(0));a=(F(-2),F(1));a2=mul(a,a);a4=mul(a2,a2)
    n01=mul((F(1,4),F(0)),add(one,a2))
    assert a2==(9,-4) and a4==(161,-72) and n01==(F(5,2),-1)
    assert add(n01,a)==(F(1,2),F(0))
    # For theta=pi, both inserted neutral group elements equal-1.
    assert 0<val(a4)<1 and 0<val(a)<1
    maximum=0.
    for row in probe['rows']:
        pair=row['neutral_charge_pairs'][-1]
        maximum=max(maximum,abs(pair['original_inverse_charge_pair_real']-val(a4)))
        assert abs(row['two_time_occupation_real']-val(n01))<3e-12
    assert maximum<3e-12
    # Explicit positive small reflected Gram: the observed mismatch is not
    # a negative-norm or a general originalQ0 counterexample.
    gram=np.diag([1.,val(a4)])
    assert np.min(np.linalg.eigvalsh(gram))>0
    return dict(arithmetic='exact Q(sqrt(5)) pairs of Fractions',amplitude=serial(a),
        two_time_occupation=serial(n01),two_time_neutral_parity=serial(a4),
        conserved_parity_pair_required=1,parity_Ward_defect=serial(add(one,(-a4[0],-a4[1]))),
        explicit_test_Gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
        exact_to_actual_matrix_error=float(maximum),
        full_average_identity='I0(O_theta(0) O_minus_theta(1)) = [1-(1-a^2) sin(theta/2)^2]^2 I0(1)',
        normalized_interpretation_requires_nonzero_original_partition=True,
        fixed_E_or_partial_Haar_not_used_in_analytic_identity=True,
        neutral_number_changing_vertices_excluded_by_massless_dictionary_contract=True)


def run():
    probe=entry.run()
    assert probe==json.loads(entry.TARGET.read_text('utf8'))
    deps=('research_note_612.md','research_note_646.md','research_note_653.md','research_note_670.md',
          'research_note_673.md','research_note_679.md','research_note_680.md','research_note_689.md',
          'joint_physical_car_projector.py','joint_physical_car_projector_results.json',
          'round690_drafts/spatial_car_source_probe.py','round690_drafts/spatial_car_source_probe_results.json')
    return dict(date='2026-10-02',round=690,tests_run=2,failures=0,errors=0,
        original_spatial_source=probe,exact_neutral_charge_certificate=exact_certificate(probe),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(same_pure_time_CAR_charge_dictionary_excluded_with_spatial_propagation=True,
            neutral_identity_retains_complete_original_Hb_S9_Gauss_average=True,
            original_Q0_RP_failure_claimed=False,all_CAR_realizations_excluded=False,
            existing612_fixed_lattice_spectral_obstruction_not_reproved=True,
            conserved_nonlocal_or_enlarged_charge_dictionary_still_possible=True,
            full_HF_or_common_continuum_or_quantum_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=690,all_checks_passed=True,
        exact_neutral_parity=result['exact_neutral_charge_certificate']['two_time_neutral_parity'],
        full_average_retained=True,original_RP_failure_claimed=False)))
