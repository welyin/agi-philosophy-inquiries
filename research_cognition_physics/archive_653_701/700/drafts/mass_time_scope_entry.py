"""700 entry: apply inherited mass/source contracts to the new699 certificate.
No new numbered science round or repeated experiment is claimed.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_offdiagonal_haar_certificate as prior
TARGET=HERE/'mass_time_scope_entry_results.json'


def run():
    data=json.loads(prior.TARGET.read_text('utf8'))
    assert prior.run()==data
    integral=data['exact_integral']
    a,b=map(F,integral['exact_B10_Qsqrt2'])
    lo=F(1414213562373095,10**15);assert lo*lo<2 and b>0
    lower=a+b*lo;alpha=F(3,10000)
    bound=F(integral['inherited_B00'])-2*alpha*lower+alpha*alpha*F(integral['inherited_B11_upper'])
    assert lower==F(integral['exact_B10_bounds'][0])
    assert bound==F(integral['quadratic_upper']) and bound<0
    deps=('research_note_643.md','research_note_655.md','research_note_667.md',
        'research_note_673.md','research_note_675.md','research_note_680.md','research_note_684.md',
        'research_note_699.md','joint_offdiagonal_haar_certificate.py',
        'joint_offdiagonal_haar_certificate_results.json','research_round_699_checks.json')
    return dict(date='2026-10-02',entry_round=700,latest_formal_round=699,new_formal_round=False,
        inherited_negative_bound=str(bound),prior_results_reproduced=True,
        derived_contract=dict(fixed_Hb_and_Wilson_and_auxiliary_baseline=True,
            mass_path_alone_cannot_restore_all_algebra_RP_at_bad_tau=True,
            inverse_half_mass_insertion_preserves_original_physical_Gauss_algebra=True,
            no_exact_positive_dilation_with_same_physical_reflection_functional=True,
            no_fixed_bad_tau_pointwise_limit_of_positive_functionals=True,
            changing_Hb_or_Wilson_or_measure_is_a_different_joint_path=True,
            shared_continuum_scaling_NOT_refuted=True,
            no_numeric_tau_threshold=True,unification_goal_unchanged=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry=700,formal=699,certificate_reproduced=True)))
