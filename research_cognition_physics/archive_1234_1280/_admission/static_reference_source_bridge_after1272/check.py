"""Finite arithmetic check of a mature modular/source dictionary.

The normal-state and source identities are analytic (proof.md). Quadrature
convergence is a numerical diagnostic, not a rigorous total physical error.
Default execution does not modify the saved result. No detector is simulated.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
H, XI, R0 = 1.0, 1.0/6.0, 0.5

def profile(r):
    r = np.asarray(r, dtype=float)
    h = np.zeros_like(r); hp = np.zeros_like(r); hpp = np.zeros_like(r)
    inside = np.abs(r) < R0
    x = r[inside]
    d = 1.0 - (x/R0)**2
    z = np.exp(1.0 - 1.0/d)
    lp = -2.0*x/(R0**2*d**2)
    lpp = -2.0/(R0**2*d**2) - 8.0*x*x/(R0**4*d**3)
    h[inside] = z
    hp[inside] = z*lp
    hpp[inside] = z*(lp*lp+lpp)
    return h, hp, hpp

def boundary(b):
    h, hp, _ = (float(v) for v in profile(np.array(b)))
    f = 1.0 - H*H*b*b
    return -4.0*math.pi*XI*b*b*(2.0*f*h*hp+H*H*b*h*h)

def charges(b, n):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    r = b*(nodes+1.0)/2.0
    weights = b*weights/2.0
    h, hp, hpp = profile(r)
    f = 1.0-H*H*r*r
    u, up, upp = h*h, 2.0*h*hp, 2.0*(hp*hp+h*hpp)
    # Q_eff = integral T_tt / f * r^2 dr dOmega, with dot(phi)=0.
    eff_density = 2.0*math.pi*r*r*(f*hp*hp+2.0*H*H*u)
    # Direct tt component of the improved Hilbert tensor difference.
    delta_density = 4.0*math.pi*XI*(-3.0*H*H*r*r*u
        -(2.0*r*f-H*H*r**3)*up-r*r*f*upp)
    eff = float(weights@eff_density)
    delta = float(weights@delta_density)
    return {'nodes': n, 'radius': b, 'effective_charge': eff,
            'physical_charge': eff+delta, 'difference_volume': delta,
            'difference_boundary': boundary(b),
            'identity_residual': abs(delta-boundary(b))}

def calculate():
    rows = [charges(b,n) for b in (0.25,0.5) for n in (64,128,256)]
    # At b=1/4, u=e^(-2/3), u'/u=-64/9. These are exact rationals.
    b = Fraction(1,4); xi = Fraction(1,6)
    bracket = (1-b*b)*Fraction(-64,9)+b
    coefficient = -4*xi*b*b*bracket
    assert coefficient == Fraction(77,288)
    exact_cut_difference = float(coefficient)*math.pi*math.exp(-2.0/3.0)
    full = rows[-1]
    cut = rows[2]
    energy_lower = math.pi*math.exp(-2.0/3.0)/48.0
    entropy_lower = 2*math.pi*energy_lower
    assert full['effective_charge'] > energy_lower > 0
    assert abs(cut['difference_volume']-exact_cut_difference) < 1e-11
    assert max(r['identity_residual'] for r in rows) < 1e-10
    refinement = [abs(rows[2]['effective_charge']-rows[1]['effective_charge']),
                  abs(rows[5]['effective_charge']-rows[4]['effective_charge'])]
    assert max(refinement) < 1e-10
    # Finite time-slice representation can be supported away from the horizon.
    maximum_radius = math.tanh(math.atanh(R0)+1.0/8.0)
    assert maximum_radius < 0.6 < 1.0
    return {
        'status': 'PASS_DICTIONARY_DIAGNOSTIC_WITH_SCOPE',
        'parameters': {'H':H,'xi':XI,'support_radius':R0,'dot_phi_initial':0},
        'charges': rows,
        'beta_central_proper_time':2*math.pi,
        'full_patch_relative_entropy':2*math.pi*full['effective_charge'],
        'analytic_energy_lower_bound':energy_lower,
        'analytic_relative_entropy_lower_bound':entropy_lower,
        'finite_cut_exact_coefficient_of_pi_exp_minus_two_thirds':str(coefficient),
        'finite_cut_source_difference':exact_cut_difference,
        'quadrature_refinement_differences':refinement,
        'timeslice_support_radius_upper_at_delta_one_eighth':maximum_radius,
        'new_numbered_rounds':0,
        'scope': [
            'Mature theorem plus an explicitly matched normal coherent-state source.',
            'Quadrature is a diagnostic; analytic proof carries the exact identities.',
            'Finite subball does not inherit the full-patch modular Hamiltonian.',
            'No actual preparation/controller source or finite detector certificate.',
            'No new gravitational response, coefficient selection or adopted axiom.']}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write:
        assert not path.exists(), 'Preserve the previous result.'
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        status = 'SAVED'
    else:
        expected = json.loads(path.read_text(encoding='utf-8'))
        assert result == expected, 'Recomputation differs from saved result.'
        status = 'READ_ONLY_MATCH'
    print(json.dumps({'status':status,'relative_entropy':result['full_patch_relative_entropy'],
        'finite_cut_source_difference':result['finite_cut_source_difference'],
        'max_identity_residual':max(x['identity_residual'] for x in result['charges'])}))

if __name__ == '__main__': main()
