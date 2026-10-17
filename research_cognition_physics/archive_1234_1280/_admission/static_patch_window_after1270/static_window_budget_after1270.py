"""Exact rational sufficient budgets; no quantum shell implementation claim.

Default is a read-only comparison with the adjacent saved JSON.
--write creates that result once and refuses replacement.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import math


def compute():
    G, L, R = F(1, 10**10), F(1), F(1, 10)
    N, xmax_energy, x0, a2 = F(643, 2), F(3, 2), F(9, 10), F(99, 100)
    # exp(20) exceeds this positive Taylor sum exactly.
    exp20_lower = sum((F(20)**j / math.factorial(j) for j in range(41)), F(0))
    gamma_tail_upper = F(4663, 3) / exp20_lower
    entropy_upper = (gamma_tail_upper + F(9, 4096)) / (1 - gamma_tail_upper)
    # pi < 22/7; a^2 = .99; r_c/L >= x0 follows from mu <= G N/L.
    geometry_upper = F(22, 7) * G * N**2 * (L/R + 4*a2/(x0*x0*(1+x0)))
    # beta=2*pi*sqrt(.99)>6. e>2, N-Emax=320.
    central_tail_upper = F(1, 2**1920)
    normalization_upper = central_tail_upper / (1-central_tail_upper)
    # Conditional fixed-radius Hamiltonian comparison at boundary time t=1.
    phase_upper = G*N*N/(2*R)
    # Same one-EOS source window; classical radial stability only.
    z, bmax = F(1,100), G*N/R
    a_lower, b_redshift_lower = F(9949,10000), F(9948,10000)
    w_upper = z / (2*a_lower*b_redshift_lower)
    stability_lower = 4*a2-2*z-2*z*z/a2-16*a2*F(51,10000)**2-6*bmax
    assert a_lower*a_lower < a2 and a_lower-bmax>b_redshift_lower
    assert bmax/4 < z/2  # w<0 because a<1
    assert w_upper==F(125000,24743163) and w_upper<F(51,10000)
    assert stability_lower==F(1949995105937,495000000000)>F(39,10)
    margin = F(3, 1000) - entropy_upper - geometry_upper - normalization_upper
    assert G*N/L < (x0-x0**3)/2
    assert (G*N/R)**2 < a2 and G*N/R < F(1, 4)
    assert F(99, 100)**2 < a2 and 2*F(314,100)*F(99,100)>6
    assert R/L < x0
    assert margin > 0
    def item(v):
        return {'numerator': str(v.numerator), 'denominator': str(v.denominator), 'diagnostic_float': float(v)}
    return {
        'scope': 'Conditional finite-source dictionary only; exact rational sufficient entropy/area budgets, no apparatus or quantum shell certification',
        'parameters': {'G':'1/10000000000','L':'1','R':'1/10','N':'643/2','Emax':'3/2','k':4,'theta':16,'B':320,'beta':'2*pi*sqrt(99/100)'},
        'gamma_tail_upper': item(gamma_tail_upper),
        'entropy_remainder_upper': item(entropy_upper),
        'area_dictionary_error_upper': item(geometry_upper),
        'corner_normalization_error_upper': item(normalization_upper),
        'fixed_radius_phase_error_upper_t1': item(phase_upper),
        'classical_source_w_absolute_upper': item(w_upper),
        'classical_source_R2Vpp_lower': item(stability_lower),
        'strict_margin_below_3_over_1000': item(margin),
        'reference_window_and_work_window_same': False,
        'classical_surface_source_defined': True,
        'same_quantum_static_source_implemented': False,
        'uncertified_inputs': ['same quantum embedding and normal radial preparation of the stated classical shell action','physical reference and clock identification','independent actual area readout'],
        'proved': ['positive shell branch on stated energy support','cosmological root above x0','combined conditional entropy-area bound strictly below 3/1000'],
    }


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    path=Path(__file__).with_suffix('.json'); result=compute()
    if args.write:
        with path.open('x',encoding='utf-8') as f: json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        action='CREATED'
    else:
        assert json.loads(path.read_text('utf-8'))==result
        action='READ_ONLY_MATCH'
    print(json.dumps({'action':action,'conditional_bound':'<0.003','entropy_upper':result['entropy_remainder_upper']['diagnostic_float'],'area_upper':result['area_dictionary_error_upper']['diagnostic_float'],'phase_upper':result['fixed_radius_phase_error_upper_t1']['diagnostic_float'],'same_quantum_static_source_implemented':False},ensure_ascii=False))
