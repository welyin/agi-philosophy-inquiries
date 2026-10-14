"""Source-realization screen, not a new GR or dust quantum model.

Checks the necessary initial momentum-compression coefficient of a proposed
dust completion of the preceding static-profile ADM preflight. Never treats
the formal internal-M substitution as a complete quantum dust implementation.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
PRE = HERE.parent / 'after1199_self_source_preflight'


def run():
    prior = json.loads((PRE/'results.json').read_text(encoding='utf-8-sig'))
    eps, m0, g = F(1, 1000), F(526, 5), F(2, 5)
    m2 = m0*m0+g*g
    xvals = [v[0] for v in prior['numerical_diagnostics']['X_moments_plus_minus']]
    force_moment = [-float(eps*m0/2)*v-float(eps*m2/48) for v in xvals]
    delta_formula = float(eps*m0)*3*np.sqrt(3)/(16*np.pi)
    assert abs(force_moment[0]-force_moment[1]-delta_formula) < 1e-14
    # Exact coefficients from integration by parts with Dirichlet boundary.
    # ∫(z-1/2)(GA)'=-∫GA=-X/2; ∫(z-1/2)(1-2z)=-1/6.
    self_integral = F(-2, 3)+F(2, 2)-F(1, 2)
    assert self_integral == F(-1, 6)
    density_self_coefficient = self_integral/8
    assert density_self_coefficient == F(-1, 48)
    delta_lower = eps*m0*F(357, 1760)/2
    target = F(1, 100)
    assert delta_lower > target
    # These are conditional operator/form bounds in the retained source symbol
    # model. Positive X and M>=m_min imply a nonzero compressive necessary term.
    m_min = m0-g
    compression_magnitude_lower = eps*m_min*m_min/48
    assert compression_magnitude_lower > F(22, 100)
    # Independent polynomial quadrature verifies the self-density moment and
    # vanishing total self-force, not the evolution of a quantum dust body.
    zz, ww = np.polynomial.legendre.leggauss(32)
    zz, ww = (zz+1)/2, ww/2
    total = float(np.dot(ww, (1-2*zz)/8))
    moment = float(np.dot(ww, (zz-.5)*(1-2*zz)/8))
    assert abs(total) < 1e-14 and abs(moment+1/48) < 1e-14
    return {
        'schema': 'after1199_static_profile_source_realization_screen_v1',
        'status': 'PASS_NECESSARY_INITIAL_FORCE_SCREEN',
        'numbered_round': None,
        'previous_goal_turn_classification': 'progress: mixed source and original-clock jet preflight with independent review',
        'prior_results_sha256': hashlib.sha256((PRE/'results.json').read_bytes()).hexdigest(),
        'exact_coefficients': {'M_X': '-epsilon/2', 'M_squared': str(density_self_coefficient)+' * epsilon'},
        'formal_retained_operator_bounds': {
            'compression_magnitude_lower': str(compression_magnitude_lower),
            'plus_minus_compression_difference_strict_lower': str(delta_lower),
            'difference_target': str(target), 'difference_margin': str(delta_lower-target)},
        'numerical_diagnostics': {'formal_compression_plus_minus': force_moment,
            'formal_compression_difference': delta_formula,
            'integrated_self_force': total, 'self_force_moment': moment},
        'decision': 'Do not extend fixed-profile oscillator integration into a new numbered GR result; prioritize independently defined operational source and actual conversion.',
        'not_certified': ['quantum_dust_preparation_or_complete_hamiltonian',
            'invariant_volume_acceleration', 'finite_time_error_of_the_retained_clock',
            'GR_generation_or_a_general_no_go'],
    }


if __name__ == '__main__':
    print(json.dumps(run(), ensure_ascii=False, indent=2))
