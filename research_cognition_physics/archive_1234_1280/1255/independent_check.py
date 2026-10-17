"""Independent round 1255 audit; standard library only; default read-only replay.

Checks exact constrained-source identities and reconstructs the finite initial
geometry independently of the author's implementation.  Analytic scope and
symplectic pullback are audited in independent_review.md, not proved by samples.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def s(x):
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), F(0))


def mv(a, v):
    return [dot(row, v) for row in a]


def vv(a, b):
    return [[x * y for y in b] for x in a]


def matrix_strings(a):
    return [[s(v) for v in row] for row in a]


def canonical_source():
    # A non-diagonal inverse metric of determinant one; density n is NOT one.
    hi = [[F(1), F(-1), F(0)], [F(-1), F(2), F(0)], [F(0), F(0), F(1)]]
    n, volume, C = F(3, 2), F(1), F(-3)
    D = [F(1, 3), F(-1, 5), F(1, 7)]
    raised = mv(hi, D)
    spatial = dot(D, raised)
    root2 = C * C - spatial
    assert root2 > 0 and C < 0
    # Solve original constraints, keeping E and M distinct.
    P = [d / n for d in D]
    E = -C / n
    mass2 = root2 / (n * n)
    assert dot(P, mv(hi, P)) + mass2 == E * E
    assert C + n * E == 0
    assert all(d - n * p == 0 for d, p in zip(D, P))
    lapse2 = C * C / root2
    beta_norm2 = spatial / root2
    assert -lapse2 + beta_norm2 == -1
    assert lapse2 * mass2 / (E * E) == 1  # positive theta-dot is +1
    assert all((-C) * p / E == d for p, d in zip(mv(hi, P), raised))
    # Source from the original material stress, not from Einstein's equation.
    rest_density = root2 / (volume * (-C))
    Tnn = rest_density * lapse2
    Tni = [-d / volume for d in D]
    Sij = [[di * dj / (volume * (-C)) for dj in D] for di in D]
    assert Tnn == n * E / volume
    assert all(Tni[i] == -n * P[i] / volume for i in range(3))
    assert all(Sij[i][j] == n * P[i] * P[j] / (volume * E)
               for i in range(3) for j in range(3))

    # Independent directional finite difference of n U - sqrt(C^2 - hDD).
    U, dU, dC = F(7, 8), F(2, 9), F(1, 7)
    dD = [F(1, 11), F(-1, 13), F(1, 17)]
    dhi = [[F(1, 13), F(1, 17), F(0)],
           [F(1, 17), F(-1, 19), F(1, 23)],
           [F(0), F(1, 23), F(1, 29)]]
    root = math.sqrt(float(root2))
    explicit_metric = float(dot(D, mv(dhi, D))) / (2 * root)
    derivative = float(n * dU) + float(-C * dC + dot(raised, dD)) / root + explicit_metric

    def hamiltonian(t):
        ct = float(C) + t * float(dC)
        dt = [float(D[i]) + t * float(dD[i]) for i in range(3)]
        hit = [[float(hi[i][j]) + t * float(dhi[i][j]) for j in range(3)] for i in range(3)]
        return float(n) * (float(U) + t * float(dU)) - math.sqrt(ct * ct - dot(dt, mv(hit, dt)))

    steps = [0.01, 0.005]
    errors = [abs((hamiltonian(t) - hamiltonian(-t)) / (2 * t) - derivative) for t in steps]
    assert errors[1] < 1e-6 and 3.5 < errors[0] / errors[1] < 4.5
    assert abs(explicit_metric) > 1e-5

    # Full internal mass exchange.  q,r and Y,K carry the ORIGINAL n-weighted PB.
    q, osc_r, B, p, X, Y, K = map(F, ['2/5', '-1/7', '3/2', '1/9', '-2/11', '1/13', '-1/17'])
    w, wp, f, fp, vp, vX = map(F, ['3/8', '-2/9', '4/7', '1/6', '2/15', '-1/19'])
    alpha, eps, g = B - 1, F(1, 4), F(1, 10)
    Jdot = -alpha * eps * wp * q * p - g * fp * X * Y
    qdot, rdot, Ydot, Kdot = osc_r, -q - alpha * eps * w * p, K, -g * f * X
    mdot = (Jdot + q * qdot + osc_r * rdot + K * Kdot
            + alpha * eps * (wp * q * p + w * qdot * p + w * q * vp)
            + g * (fp * X * Y + f * vX * Y + f * X * Ydot))
    desired = alpha * eps * w * q * vp + g * f * Y * vX
    assert mdot == desired
    # Reduced Hamilton equations have no fictitious lapse factor on internal pairs.
    assert (n * K) / n == K and (-n * g * f * X) / n == -g * f * X
    return {
        'status': 'PASS', 'n': s(n), 'C0': s(C), 'D0': [s(d) for d in D],
        'inverse_metric': matrix_strings(hi), 'root_squared': s(root2),
        'normal_energy': s(Tnn), 'normal_spatial_stress': [s(x) for x in Tni],
        'spatial_stress_covariant': matrix_strings(Sij), 'g00': '-1',
        'gauge_matrix_det_squared': s(root2 / (C * C)),
        'variation_error_h': format(errors[0], '.12g'),
        'variation_error_h_over_2': format(errors[1], '.12g'),
        'omitted_explicit_metric_term': format(explicit_metric, '.12g'),
        'mass_exchange_residual': '0',
        'scope': 'finite algebra and diagnostic differences; analytic pullback needs the proof'
    }


def geometric_witness():
    k2, kappa, H2 = F(1, 10), F(3), F(111, 100)
    X, q = F(1, 10), F(1, 2)
    n_clock, n_record, record_mass = F(1), F(1, 10), F(1, 5)
    U = F(7, 10) + F(1, 100) + q * q / 2
    assert U == F(167, 200)
    out = {}
    for d in [F(0), F(1, 10)]:
        # Rebuild every non-clock matter contribution from the initial material.
        X_potential = X * X / 2
        rod_gradient = 3 * k2 / 2
        rod_time_kinetic = d * d / (2 * k2)  # Pi_1=d/k; F(p)=1, sqrt(h)=1
        record_energy = n_record * record_mass  # P_record = 0
        Cg = -3 * H2 / kappa
        C0 = Cg + X_potential + rod_gradient + rod_time_kinetic + record_energy
        A = -C0
        root2 = A * A - d * d
        assert A > 0 and root2 > U * U > 0
        # Positive full clock mass is sqrt(root2); J is mass minus U, not mass itself.
        assert C0 + A == 0 and d - n_clock * d == 0
        rho = X_potential + rod_gradient + rod_time_kinetic + record_energy + A
        assert rho == H2 and 3 * H2 == kappa * rho
        # Scalar/rod spatial stress derived directly, plus moving clock and resting record.
        iso = -X_potential - k2 / 2 + rod_time_kinetic
        S = [iso + d * d / A, iso, iso]
        trace = sum(S)
        stf = [v - trace / 3 for v in S]
        E = [-kappa * v / 2 for v in stf]
        magnetic_squared = F(0)  # flat spatial metric and constant isotropic extrinsic K
        weyl_squared = 8 * (sum(v * v for v in E) - magnetic_squared)
        assert (-A * A + d * d) / root2 == -1
        out[s(d)] = {
            'clock_independent_U': s(U), 'C0': s(C0), 'D01': s(d),
            'clock_E': s(A), 'clock_mass_squared': s(root2),
            'positive_J_certificate_M2_minus_U2': s(root2 - U * U),
            'rho': s(rho), 'rod_time_kinetic': s(rod_time_kinetic),
            'clock_normal_spatial_component': s(-d),
            'rod_normal_spatial_component': s(d), 'record_normal_spatial_component': '0',
            'total_spatial_stress_diagonal': [s(v) for v in S],
            'electric_Weyl_diagonal': [s(v) for v in E],
            'magnetic_Weyl_squared': '0', 'Weyl_squared': s(weyl_squared),
            'lapse_squared': s(A * A / root2),
            'shift1_squared': s(d * d / root2),
            'record_coordinate_velocity': '-beta (P_record=0; not comoving with clock)',
            'record_phase_rate': 'N/10 at p=0; not the old homogeneous 1/10 law'
        }
    assert out['0']['Weyl_squared'] == '0'
    assert out['1/10']['electric_Weyl_diagonal'] == ['-2/177', '1/177', '1/177']
    assert out['1/10']['Weyl_squared'] == '16/10443'
    assert out['1/10']['clock_mass_squared'] == '30929/40000'
    return {'status': 'PASS', 'branches': out,
            'scope': 'initial common-source geometry; not a propagated quantum/readout result'}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write', action='store_true', help='exclusive first save only')
    args = ap.parse_args()
    inputs = [HERE / 'admission.md', HERE / 'proof.md', HERE / 'NEXT.md', HERE / 'check.py', HERE / 'results.json',
              HERE.parent / 'research_note_1255.md', Path(__file__).resolve()]
    missing = [str(p) for p in inputs if not p.is_file()]
    if missing:
        raise SystemExit('Await author files before final binding: ' + ', '.join(missing))
    result = {'round': 1255, 'status': 'PASS', 'independent_groups': 2,
              'canonical_source': canonical_source(), 'finite_geometric_witness': geometric_witness(),
              'input_sha256': {str(p.relative_to(HERE.parent)): sha(p) for p in inputs}}
    dst = HERE / 'independent_results.json'
    if args.write:
        with dst.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write('\n')
    else:
        saved = json.loads(dst.read_text(encoding='utf-8'))
        if result != saved:
            raise SystemExit('FAIL: saved independent results or bound author versions differ')
    print(('SAVED' if args.write else 'READ_ONLY_MATCH') + ': 1255 independent canonical/source and finite geometric witness; two groups; PASS')


if __name__ == '__main__':
    main()