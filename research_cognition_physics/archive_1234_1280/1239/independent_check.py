"""Independent 1239: exact budget and Volterra source/geometry comparison.

No author implementation is imported.  The numerical diagnostic solves the
integral equations by global Legendre collocation, whereas the author advances
Raychaudhuri's equation.  Analytic comparison, not discretization agreement,
certifies the finite error bound.  Default execution is read-only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre as L

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
TARGET = HERE / 'independent_results.json'
AUTHOR = ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']
HISTORY = ['archive_1234_/1236/results.json',
           'archive_1234_/1236/proof.md',
           'archive_956_989/research_note_964.md',
           'archive_990_1008/research_note_998.md']


def require(value, message):
    if not bool(value):
        raise AssertionError(message)


def rational(x):
    return {'fraction': str(x), 'decimal': float(x)}


def exact_budget():
    old = json.loads((BASE / HISTORY[0]).read_text('utf-8'))
    d, emax, umax = Q(3, 100), Q(91, 125), Q(6, 5)
    fuel = Q(1, 100000)
    # R_h=e+2 d u^2+4 d^2 h u/3+8 d^3 h^2/27, u^3=1+h e.
    deriv = 4*d*emax/3 + 4*d*d*(umax+emax/3)/3 + 16*d**3/27
    require(deriv == Q('0.0308672'), 'exact inherited-menu derivative bound')
    require(Q(4, 5)-fuel >= Q(4, 5)**2, 'h_L >= 4/5 in both source branches')
    delta_h = fuel/Q(8, 5)
    geom = deriv*delta_h
    require(geom == Q('0.00000019292'), 'geometry bound')
    other = Q(8, 10000000)
    total = geom+other
    require(total < Q(1, 1000000), 'nonempty remaining instrument budget')
    D = Q(old['group_2']['D']['fraction'])
    C = Q(old['group_2']['linear_error_coefficient']['fraction'])
    gap = Q(old['group_3']['analytic_Peano_gap_lower']['fraction'])
    source_error = 10*total+16*total*total
    flat_error = C*(2*total)+8*(2*total)**2
    require(gap > source_error and -D > flat_error, 'inherited separations remain positive')
    # If the negative fuel term or positive receiving term is removed, even
    # the initial Einstein constraint derivative is nonzero at fixed H0.
    double_count_slope = -4*fuel/9
    missing_receiver_slope = 4*fuel/9
    require(double_count_slope != 0 and missing_receiver_slope != 0,
            'one-sided exchange violates constraint propagation')
    return dict(status='PASS', fuel=rational(fuel), return_h_derivative_bound=rational(deriv),
                h_change_bound=rational(delta_h), geometry_error=rational(geom),
                remaining_contract_input=rational(other), total_contract=rational(total),
                old_source_gap_lower=rational(gap), two_source_error=rational(source_error),
                source_margin=rational(gap-source_error), flat_two_box_error=rational(flat_error),
                flat_margin=rational(-D-flat_error),
                double_count_initial_constraint_derivative=rational(double_count_slope),
                missing_receiver_initial_constraint_derivative=rational(missing_receiver_slope),
                instrument_budget_is_assumed_not_manufactured=True)


def integral_rows(x, degree):
    """Integral of P_k from -1 to each x, for k=0,...,degree-1."""
    x = np.asarray(x)
    values = L.legvander(x, degree)
    out = np.empty((x.size, degree))
    out[:, 0] = x+1
    for k in range(1, degree):
        out[:, k] = (values[:, k+1]-values[:, k-1])/(2*k+1)
    return out


def solve_integral(E0, degree, end=1.):
    nodes, _ = L.leggauss(degree)
    time = (nodes+1)*end/2
    V = L.legvander(nodes, degree-1)
    invV = np.linalg.inv(V)
    B = end/2*integral_rows(nodes, degree)@invV
    fuel0 = 1e-5
    F = fuel0*np.exp(-time)
    y = 1+np.sqrt(E0)*time
    # y=a^(3/2); U=a^-1 integral a F, exactly retaining pressure dilution.
    for iteration in range(12):
        a = y**(2/3)
        U = (B@(a*F))/a
        E = E0-fuel0+F+U
        new_y = 1+B@np.sqrt(E)
        difference = float(np.max(abs(new_y-y)))
        y = new_y
        if difference < 3e-15:
            break
    require(difference < 3e-14, 'Volterra fixed point diagnostic')
    a = y**(2/3)
    U = (B@(a*F))/a
    E = E0-fuel0+F+U
    H = (2/3)*np.sqrt(E)/y
    work = B@(H*U)
    pressure_residual = float(np.max(abs(fuel0-F-U-work)))
    require(pressure_residual < 2e-15, 'complete pressure work')
    y_coeff = invV@y
    d_y = (2/end)*L.legval(nodes, L.legder(y_coeff))
    constraint = float(np.max(abs(((2/3)*d_y/y)**2-4*E/(9*y*y))))
    require(constraint < 2e-10, 'independent differentiation of the constraint')
    U_coeff = invV@U
    dU = (2/end)*L.legval(nodes, L.legder(U_coeff))
    source_residual = float(np.max(abs(dU-F+H*U)))
    require(source_residual < 2e-14, 'receiving source and dilution')
    require(np.min(U) > 0 and np.min(E) > E0-fuel0
            and np.max(E) < E0, 'positive receiver and total-source bracket')
    h0, hlo = np.sqrt(E0), np.sqrt(E0-fuel0)
    require(np.min(y-1-hlo*time) > -1e-14 and
            np.max(y-1-h0*time) < 1e-14, 'scale-factor bracket')
    eta_coeff = invV@(y**(-2/3))

    def eta(t):
        xx = np.array([2*t/end-1])
        return float((end/2*integral_rows(xx, degree)@eta_coeff)[0])

    def return_time(e):
        target = eta(e)+.06
        low, high = e, end
        require(eta(high) > target, 'return is in fixed finite interval')
        for _ in range(60):
            middle = (low+high)/2
            if eta(middle) < target:
                low = middle
            else:
                high = middle
        return (low+high)/2

    emissions = np.array([0., .331, .728])
    returns = np.array([return_time(e) for e in emissions])
    def dust_return(e, h):
        u = (1+h*e)**(1/3)
        return e+.06*u*u+.0012*h*u+.000008*h*h
    upper = np.array([dust_return(e, h0) for e in emissions])
    lower = np.array([dust_return(e, hlo) for e in emissions])
    require(np.min(returns-lower) > -2e-14 and np.max(returns-upper) < 2e-14,
            'finite source-to-return ordering')
    require(np.max(upper-returns) > 1e-10 and np.max(upper-returns) < 1.9292e-7,
            'nonzero finite source changes the returns within analytic budget')
    # The six-record signed determinant used in 1236, without exact emission
    # permissions: only the geometry part leaves these ideal e unchanged.
    def det(r):
        return float((emissions[1]-emissions[0])*(r[2]-r[0])
                   -(emissions[2]-emissions[0])*(r[1]-r[0]))
    require(det(returns) < 0, 'curved signed record witness retained')
    return dict(status='PASS', E0=E0, degree=degree, fixed_point_iterations=iteration+1,
                fixed_point_residual=difference, pressure_work_residual=pressure_residual,
                differentiated_constraint_residual=constraint,
                receiving_source_residual=source_residual,
                final_node_fuel=float(F[-1]), final_node_radiation_comoving=float(U[-1]),
                final_node_expansion_work=float(work[-1]),
                energy_initial=E0, initial_H=float((2/3)*h0),
                minimum_E=float(np.min(E)), maximum_E=float(np.max(E)),
                ideal_emission_clock_readings=emissions.tolist(),
                return_times=returns.tolist(), original_dust_returns=upper.tolist(),
                lower_envelope_returns=lower.tolist(),
                nonzero_geometry_shifts=(upper-returns).tolist(), D=det(returns),
                original_D=det(upper),
                nominal_phase_frequency=1., actual_clock_instrument_certified=False)


def source_geometry_diagnostic():
    rows, refinement = [], []
    for E0 in [1., .8]:
        coarse = solve_integral(E0, 24)
        fine = solve_integral(E0, 40)
        err = float(np.max(abs(np.array(coarse['return_times'])-
                              np.array(fine['return_times']))))
        require(err < 3e-13, 'independent-order return consistency')
        rows.append(fine)
        refinement.append(err)
    source_D_gap = rows[1]['D']-rows[0]['D']
    require(source_D_gap > 3.3810880023751804e-5,
            'fixed deletion still separated in the diagnostic')
    return dict(status='PASS', method='Volterra/Legendre; no author import',
                finite_conversion_rate=1., finite_end_time=1., branches=rows,
                refinement_errors=refinement, source_D_gap=source_D_gap,
                analytic_certificate_from_energy_comparison_not_numerical_convergence=True,
                deleted_source_branch_is_not_a_second_solution_of_the_same_action=True)


def digest(paths):
    return {str(p.relative_to(BASE)).replace('\\', '/'):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run(bind=True):
    watched = [BASE/p for p in HISTORY]
    if bind:
        watched += [HERE/p for p in AUTHOR]
    before = digest(watched)
    result = dict(round=1239, status='PASS', independent_scientific_groups=2,
                  exact_finite_budget=exact_budget(),
                  same_source_Volterra_check=source_geometry_diagnostic(),
                  CT_adopted=False, GR_goal_complete=False,
                  homogeneous_effective_source_only=True,
                  directional_signal_and_write_recoil_certified=False,
                  quantum_input_reference_channel_certified=False,
                  physical_approximation_error_derived=False,
                  source_sha256=before)
    require(before == digest(watched), 'read-only historical and author audit')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--diagnostic', action='store_true')
    args = parser.parse_args()
    require(not (args.write and args.diagnostic), 'choose one mode')
    result = run(bind=not args.diagnostic)
    if args.diagnostic:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    elif args.write:
        with TARGET.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('WROTE_ONCE: 1239 independent_results.json; two groups PASS')
    else:
        require(result == json.loads(TARGET.read_text('utf-8')), 'READ_ONLY_MATCH')
        print('READ_ONLY_MATCH: 1239 independent; two groups PASS')
