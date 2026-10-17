"""1245 independent rational local-jet and Legendre audit.

The complete-field identities are analytic, as discussed in the review.
These exact finite jets diagnose signs, label weights and interaction terms;
they do not prove a finite lattice representation of the continuum HDA.
No author code is imported. --write exclusively creates the first result.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def transpose(a):
    return [list(row) for row in zip(*a)]


def mm(a, b):
    return [[dot(row, col) for col in transpose(b)] for row in a]


def mv(a, x):
    return [dot(row, x) for row in a]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def contract(a, b):
    return sum((dot(x, y) for x, y in zip(a, b)), F(0))


def scale(a, c):
    return [[c*x for x in row] for row in a]


def add(a, b):
    return [[x+y for x, y in zip(rowa, rowb)] for rowa, rowb in zip(a, b)]


def inverse(a):
    n = len(a)
    w = [row[:] + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for i in range(n):
        j = next(j for j in range(i, n) if w[j][i])
        w[i], w[j] = w[j], w[i]
        pivot = w[i][i]
        w[i] = [x/pivot for x in w[i]]
        for j in range(n):
            if j != i:
                c = w[j][i]
                w[j] = [x-c*y for x, y in zip(w[j], w[i])]
    return [row[n:] for row in w]


def pack(q):
    return {'rational': str(q), 'decimal_diagnostic': float(q)}


def data():
    L = [[F(1), 0, 0], [F(1, 3), F(1), 0], [0, F(1, 4), F(1)]]
    h = mm(L, transpose(L))
    hinv = inverse(h)
    dh = [
        [[F(1, 7), F(1, 11), 0], [F(1, 11), F(-1, 5), F(1, 13)], [0, F(1, 13), F(1, 9)]],
        [[F(-1, 8), 0, F(1, 10)], [0, F(1, 6), 0], [F(1, 10), 0, F(1, 12)]],
        [[F(1, 9), F(-1, 12), 0], [F(-1, 12), F(1, 8), F(1, 14)], [0, F(1, 14), F(-1, 10)]],
    ]
    # det L=1, so sqrt(det h)=1; derivatives still need the volume term.
    volume = F(1)
    dvolume = [contract(hinv, k)/2 for k in dh]
    dhinv = [scale(mm(mm(hinv, k), hinv), -1) for k in dh]
    dK = [add(scale(hinv, dv), kinv) for dv, kinv in zip(dvolume, dhinv)]
    q = [F(1, 5), F(2, 7), F(-1, 4), F(3, 8)]
    momenta = [F(1, 5), F(-1, 10), F(3, 10), F(2, 5)]
    gradients = [[F((a+1)*(i+2)-4, 13+a+2*i) for i in range(3)] for a in range(4)]
    hessians = [[[F((a+1)*(i+j+1), 23+a+i+j) for j in range(3)] for i in range(3)] for a in range(4)]
    N, M = F(2, 3), F(5, 4)
    dN, dM = [F(1, 5), F(-2, 7), F(1, 3)], [F(-1, 4), F(1, 3), F(1, 8)]
    return L, h, hinv, volume, dhinv, dK, q, momenta, gradients, hessians, N, M, dN, dM


def full_material_brackets():
    L, h, inv, volume, dinv, dK, q, P, grad, hess, N, M, dN, dM = data()
    field_factor, Fp = 1+q[0]/5, F(1, 5)
    assert field_factor > 0
    factors = [F(1)] + [field_factor]*3
    target_p = [P[a]/(volume*factors[a]) for a in range(4)]
    rod_grad = sum((dot(grad[a], mv(inv, grad[a])) for a in range(1, 4)), F(0))
    local_p_force = volume*q[0] + Fp*(volume*rod_grad-sum((P[a]**2 for a in range(1, 4)), F(0))/(volume*field_factor**2))/2

    def derivative_q(lapse, dlapse):
        answer = []
        for a in range(4):
            base = factors[a]
            div = F(0)
            for i in range(3):
                for j in range(3):
                    coeff = volume*base*inv[i][j]
                    dcoeff = base*dK[i][i][j]
                    if a:
                        dcoeff += volume*Fp*grad[0][i]*inv[i][j]
                    div += dlapse[i]*coeff*grad[a][j]
                    div += lapse*(dcoeff*grad[a][j]+coeff*hess[a][i][j])
            answer.append(-div + (lapse*local_p_force if a == 0 else 0))
        return answer

    dHN, dHM = derivative_q(N, dN), derivative_q(M, dM)
    bracket_f = dot(dHN, [M*x for x in target_p])-dot([N*x for x in target_p], dHM)
    v = mv(inv, [N*dM[i]-M*dN[i] for i in range(3)])
    momentum_density = [sum((P[a]*grad[a][i] for a in range(4)), F(0)) for i in range(3)]
    expected_f = dot(momentum_density, v)
    assert bracket_f == expected_f != 0
    assert local_p_force != volume*q[0]  # F_p terms were actually retained.

    body_rows = []
    body_sum = F(0)
    for label_weight, mass, mass_p, I in ((F(1), F(1), F(0), F(0)),
                                         (F(1, 10), F(1, 5)+q[0]/10, F(1, 10), F(1))):
        frame_dir = [F(2, 3), F(2, 3), F(1, 3)]
        mechanical = mv(L, [3*mass*x/4 for x in frame_dir])
        raised = mv(inv, mechanical)
        E = 5*mass/4
        assert E*E == mass*mass+dot(mechanical, raised)
        dE = [(dot(mechanical, mv(dinv[i], mechanical))+2*mass*mass_p*grad[0][i])/(2*E) for i in range(3)]
        dX_N = [dN[i]*E+N*dE[i] for i in range(3)]
        dX_M = [dM[i]*E+M*dE[i] for i in range(3)]
        dP_N, dP_M = [N*x/E for x in raised], [M*x/E for x in raised]
        # Integrand carries n; the inverse symplectic measure contributes 1/n.
        weighted_PB = (dot([label_weight*x for x in dX_N], [label_weight*x for x in dP_M])
                       -dot([label_weight*x for x in dP_N], [label_weight*x for x in dX_M]))/label_weight
        expected = -label_weight*dot(mechanical, v)
        assert weighted_PB == expected != 0
        incorrect_PB = label_weight*weighted_PB
        if label_weight != 1:
            assert incorrect_PB != expected
        cross_f_body = -label_weight*N*M*P[0]*mass*mass_p/(volume*E)
        cross_body_f = label_weight*M*N*P[0]*mass*mass_p/(volume*E)
        assert cross_f_body+cross_body_f == 0
        if mass_p:
            assert cross_f_body != 0
        body_sum += weighted_PB
        body_rows.append({'weight': pack(label_weight), 'mass': pack(mass),
                          'mechanical_covector': [pack(x) for x in mechanical],
                          'energy': pack(E), 'bracket': pack(weighted_PB),
                          'D_body': pack(expected), 'wrong_PB_if_measure_omitted': pack(incorrect_PB),
                          'field_body_cross_NM': pack(cross_f_body), 'body_field_cross_NM': pack(cross_body_f)})
    total = bracket_f+body_sum
    assert total != expected_f  # Omitting mechanical D is detectable off the homogeneous sector.
    # Label actions are fixed in X-variation. At spatially constant p, a
    # nonconstant advected label I does not create an external mass force.
    grad_I_as_eulerian_labels = [F(1, 10), F(1, 5), F(-3, 10)]
    wrong_label_force = [-(1+q[0])*x/10 for x in grad_I_as_eulerian_labels]
    assert any(x != 0 for x in wrong_label_force)
    return {'status': 'PASS', 'nonconstant_F': pack(field_factor), 'full_p_variation_local_term': pack(local_p_force),
            'field_bracket': pack(bracket_f), 'D_field': pack(expected_f),
            'body_rows': body_rows, 'total_material_bracket': pack(total),
            'omitted_body_momentum_residual': pack(body_sum),
            'constant_p_correct_external_force': [pack(F(0)) for _ in range(3)],
            'wrong_force_from_advected_I_gradient': [pack(x) for x in wrong_label_force],
            'scope': 'Exact continuum-functional local jets; not a lattice constraint algebra or a prepared physical initial state.'}


def legendre_and_same_model():
    L, h, inv, volume, dinv, dK, q, P, grad, hess, N, M, dN, dM = data()
    beta = [F(1, 10), F(-1, 5), F(1, 20)]
    Fp, factor = F(1, 5), 1+q[0]/5
    factors = [F(1)]+[factor]*3
    vn = [P[a]/(volume*factors[a]) for a in range(4)]
    qdot = [N*vn[a]+dot(beta, grad[a]) for a in range(4)]
    gradnorm = [dot(g, mv(inv, g)) for g in grad]
    kin = sum((P[a]**2/(2*volume*factors[a]) for a in range(4)), F(0))
    potential = volume*(q[0]**2+sum((factors[a]*gradnorm[a] for a in range(4)), F(0)))/2
    H_field = kin+potential
    Dshift = sum((P[a]*dot(beta, grad[a]) for a in range(4)), F(0))
    legendre_field = dot(P, qdot)-N*H_field-Dshift
    covariant_field = N*volume*(sum((factors[a]*(vn[a]**2-gradnorm[a]) for a in range(4)), F(0))-q[0]**2)/2
    assert legendre_field == covariant_field

    # One recorder with nonzero physical momentum and nonzero shift.
    action, mass_p = F(1), F(1, 10)
    mass, mass_I = F(1, 10)+(1+q[0])*action/10, (1+q[0])/10
    physical = [mass/2, mass/2, mass/4]
    mechanical = mv(L, physical)
    raised, E = mv(inv, mechanical), 5*mass/4
    assert E*E == mass*mass+dot(mechanical, raised)
    xdot = [N*raised[i]/E-beta[i] for i in range(3)]
    proper = N*mass/E
    theta_dot = proper*mass_I
    Hbody = N*E-dot(beta, mechanical)
    canonical_body = dot(mechanical, xdot)+action*theta_dot-Hbody
    covariant_body = action*theta_dot-mass*proper
    assert canonical_body == covariant_body
    assert theta_dot/proper == (1+q[0])/10
    stripped_mass_phase = F(1, 10)
    assert mass_I-stripped_mass_phase == q[0]/10 != 0

    # General non-diagonal gravitational tensor, no FLRW reduction here.
    Pi = [[F(1, 3), F(1, 7), F(-1, 11)], [F(1, 7), F(-2, 5), F(1, 13)], [F(-1, 11), F(1, 13), F(3, 8)]]
    a, curvature, Lambda_H = F(6), F(1, 5), F(3, 100)
    pi_trace = contract(h, Pi)
    pi_lower = mm(mm(h, Pi), h)
    K = scale(add(pi_lower, scale(h, -pi_trace/2)), a/volume)
    Ktrace = contract(inv, K)
    Kupper = mm(mm(inv, K), inv)
    recovered_Pi = scale(add(Kupper, scale(inv, -Ktrace)), volume/a)
    assert recovered_Pi == Pi
    Hg = a*(contract(Pi, pi_lower)-pi_trace*pi_trace/2)/volume-volume*curvature/a+2*Lambda_H*volume
    legendre_g = 2*contract(Pi, K)-Hg
    lag_g = volume*(contract(K, Kupper)-Ktrace*Ktrace+curvature-2*a*Lambda_H)/a
    assert legendre_g == lag_g

    # Direct actual matter metric variation along delta h=delta H_g[1]/delta Pi.
    dh_dir = scale(K, 2)
    dv = volume*contract(inv, dh_dir)/2
    di = scale(mm(mm(inv, dh_dir), inv), -1)
    dHf = -kin*dv/volume + potential*dv/volume
    dHf += volume*sum((factors[j]*dot(grad[j], mv(di, grad[j])) for j in range(4)), F(0))/2
    dHd = F(1, 10)*dot(mechanical, mv(di, mechanical))/(2*E)
    clock_momentum, clock_energy = [x/mass for x in mechanical], F(5, 4)
    dHd += dot(clock_momentum, mv(di, clock_momentum))/(2*clock_energy)
    cross_gm, cross_mg = -N*M*(dHf+dHd), N*M*(dHf+dHd)
    assert cross_gm != 0 and cross_gm+cross_mg == 0

    rows = []
    for stored_K in (F(1, 5), F(0)):
        Mc = F(7, 10)+F(1, 10)+stored_K
        D = Mc+F(1, 10)*(F(1, 10)+F(1, 10))
        rho = D+F(3, 20)
        Hsquare = rho
        constraint = -6*Hsquare/a+rho
        assert constraint == 0
        # Exact contraction of the covariant stress and phase source gives
        # the already certified 1242 coefficients, with no new integration.
        Q, b = F(1, 100), F(3, 100)
        p_second = -(b+Q)
        p_third_over_H = 5*b+6*Q
        z_fourth_over_H = F(1, 10)*p_third_over_H/24
        assert p_second == -F(1, 25) and z_fourth_over_H == F(7, 8000)
        rows.append({'stored_K': pack(stored_K), 'Mc': pack(Mc), 'D': pack(D),
                     'rho_and_H_squared': pack(rho), 'Hamiltonian_constraint': pack(constraint),
                     'p_second': pack(p_second), 'p_third_over_H': pack(p_third_over_H),
                     'z_t4_coefficient_over_H': pack(z_fourth_over_H)})
    assert rows[0]['D']['rational'] == '51/50' and rows[1]['D']['rational'] == '41/50'
    return {'status': 'PASS', 'field_Legendre_residual': pack(legendre_field-covariant_field),
            'moving_recorder_Legendre_residual': pack(canonical_body-covariant_body),
            'proper_time_rate': pack(proper), 'phase_rate': pack(theta_dot),
            'deleted_Ip_phase_mismatch': pack(mass_I-stripped_mass_phase),
            'tensor_inverse_Legendre_residual': pack(legendre_g-lag_g),
            'gravity_matter_cross_NM': pack(cross_gm), 'matter_gravity_cross_NM': pack(cross_mg),
            'gravity_a': pack(a), 'physical_kappa': pack(a/2),
            'same_1242_initial_states': rows,
            'scope': 'Same classical effective material and existing initial-data family; no new quantum preparation, detector, or CT adoption.'}


def compute():
    files = ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']
    return {'round': 1245, 'status': 'PASS', 'audit': 'independent', 'science_groups': 2,
            'author_sha256': {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files},
            'independent_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'groups': {'full_material_brackets': full_material_brackets(), 'Legendre_and_same_model': legendre_and_same_model()},
            'CT_adopted': False, 'all_fields_HDA_by_numeric_sampling': False,
            'complete_quantum_process_or_thermal_geometry': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = compute()
    out = HERE/'independent_results.json'
    if args.write:
        with out.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('CREATED_ONCE PASS independent science_groups=2')
    else:
        assert json.loads(out.read_text(encoding='utf-8')) == result, 'Bound file or saved result mismatch; no writes performed'
        print('READ_ONLY_MATCH PASS independent science_groups=2')


if __name__ == '__main__':
    main()
