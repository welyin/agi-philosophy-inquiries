"""1129: exact budgets and independent finite diagnostics.

No finite matrix certifies the continuum scalar vacuum, CCR, or ultraviolet limit.
Default is read-only comparison. Use --write-results to register results explicitly.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import itertools
import json
from math import factorial
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TOL = 2e-10


def norm(a):
    return float(np.linalg.norm(a))


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def evolution(h, t=1.0):
    e, q = np.linalg.eigh(h)
    return (q*np.exp(-1j*t*e)) @ q.conj().T


def annihilators(n):
    ans = []
    for j in range(n):
        a = np.zeros((2**n, 2**n), complex)
        for bits in range(2**n):
            if (bits >> j) & 1:
                a[bits ^ (1 << j), bits] = (-1)**((bits & ((1 << j)-1)).bit_count())
        ans.append(a)
    return ans


def creation(v, ann):
    return sum(v[j]*ann[j].conj().T for j in range(len(ann)))


def quadratic(h, ann):
    return sum(h[i, j]*ann[i].conj().T @ ann[j]
               for i in range(len(h)) for j in range(len(h)))


def lift(u):
    n = len(u)
    out = np.zeros((2**n, 2**n), complex)
    out[0, 0] = 1
    for k in range(1, n+1):
        subs = list(itertools.combinations(range(n), k))
        for rows in subs:
            for cols in subs:
                out[sum(1 << i for i in rows), sum(1 << j for j in cols)] = np.linalg.det(u[np.ix_(rows, cols)])
    return out


def exact_budget():
    epsilon, amplitude, lam = F(1, 10**5), F(10**34), F(11, 7)
    q = 8/epsilon**3
    c = 8*q*524
    vacuum_std = 3/epsilon
    material = 80*lam*epsilon
    source = lam*c*512*vacuum_std/amplitude
    back = lam*512*c*q/(2*amplitude**2)
    endpoint = lam*c*512*(vacuum_std/amplitude+back)
    tail = lam**10/F(factorial(10))/(1-lam/11)
    causal = 4*tail+2*source
    handoff = material+source
    vnorm = 16*lam*c
    input_h0 = 268+1073*amplitude**2
    output_h0 = input_h0+2*vnorm
    assert handoff < F(1, 100)
    assert causal < F(1, 1000)
    assert endpoint < F(84, 10**8)
    assert F(9)**2*2*21**2 < 268**2
    items = dict(epsilon=epsilon, amplitude=amplitude, rho_l2_sq_upper=q,
                 C_upper=c, material_error=material, source_error=source,
                 handoff_error=handoff, dyson_tail_n10=tail, causal_signal=causal,
                 normalized_backreaction=back, endpoint_interaction=endpoint,
                 interaction_norm_upper=vnorm, input_free_energy_upper=input_h0,
                 output_free_energy_upper=output_h0)
    return dict(scope="Exact rational upper bounds; pi < 22/7, no floating-point certification.",
                fractions={k: str(v) for k, v in items.items()},
                decimal_views={k: float(v) for k, v in items.items()},
                handoff_target="1/100", finite_menu_causal_target="1/1000")


def kernel_spectrum():
    nodes, weights = np.polynomial.legendre.leggauss(240)
    q = np.exp(-1/(1-nodes**2))
    dq = -2*nodes*q/(1-nodes**2)**2
    q2 = float(weights @ (q*q))
    dq2 = float(weights @ (dq*dq))
    gradient_sq = 3*dq2/q2
    radial = (nodes+1)/2
    rw = weights/2
    bump = np.exp(-1/(1-radial**2))
    z = float(4*np.pi*(rw @ (radial**2*bump)))
    rho_sq = float(4*np.pi*(rw @ (radial**2*bump**2))/z**2)
    mu1 = float(4*np.pi*(rw @ (radial**3*bump))/z)
    s = (nodes+1)/2
    gate_b = np.exp(-1/(s*(1-s)))
    ib = float(rw @ gate_b)
    gate_sup = float(2*np.exp(-4)/ib)
    gate_derivative = 4*gate_b*(1-2*s)/(s*s*(1-s)**2)/ib
    assert q2 > 1/16 and dq2 < 8 and gradient_sq < 400
    assert z > 1/8 and rho_sq < 8 and mu1 < 1
    assert gate_sup < 16 and max(abs(gate_derivative)) < 512

    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.])
    alpha = [np.kron(sx, u) for u in (sx, sy, sz)]
    beta = np.kron(sz, np.eye(2))
    es, vs = np.linalg.eigh(alpha[2])
    xi = vs[:, np.argmax(es)]
    minus_weight, total_weight = 0., 0.
    gap, spectrum_error = float("inf"), 0.
    negative_levels = positive_levels = 0
    for p in itertools.product((-1.2, -0.4, 0.4, 1.2), repeat=3):
        p = np.array(p)
        h = beta+sum(p[j]*alpha[j] for j in range(3))
        e, v = np.linalg.eigh(h)
        w = float(np.exp(-np.dot(p, p)))
        sm = v[:, e < 0] @ v[:, e < 0].conj().T
        minus_weight += w*float(np.vdot(xi, sm @ xi).real)
        total_weight += w
        gap = min(gap, min(abs(e)))
        spectrum_error = max(spectrum_error, norm(h @ h-(1+np.dot(p, p))*np.eye(4)))
        negative_levels += int(sum(e < 0))
        positive_levels += int(sum(e > 0))
    occupation = minus_weight/total_weight
    assert abs(occupation-.5) < TOL and spectrum_error < TOL

    # A finite range diagnostic: every power < 10 has zero matrix entry 0 -> 10.
    hopping = (np.diag(np.ones(11), 1)+np.diag(np.ones(11), -1))/2
    powers = [np.linalg.matrix_power(hopping, n)[10, 0] for n in range(10)]
    tail = (11/7)**10/factorial(10)/(1-(11/7)/11)
    actual = abs(evolution(hopping, np.pi/2)[10, 0])
    longlink = hopping.copy()
    longlink[10, 0] = longlink[0, 10] = .5
    wrong_range = abs(evolution(longlink, np.pi/2)[10, 0])
    assert max(abs(np.array(powers))) == 0 and actual <= tail
    assert wrong_range > 100*tail
    return dict(scope="Quadrature and paired Dirac momenta are diagnostics, not compact continuum Fourier certification.",
                q_l2_sq=q2, q_derivative_l2_sq=dq2, normalized_gradient_sq=gradient_sq,
                radial_normalizer=z, rho_l2_sq=rho_sq, rho_first_moment=mu1,
                gate_sup=gate_sup, sampled_gate_derivative_sup=float(max(abs(gate_derivative))),
                paired_momentum_negative_weight=occupation, positive_levels=positive_levels,
                negative_levels=negative_levels, sampled_dirac_gap=float(gap),
                dirac_square_identity_error=spectrum_error,
                range_diagnostic_amplitude=float(actual), range_tail_bound=tail,
                forbidden_longlink_amplitude=float(wrong_range))


def car_code():
    ann = annihilators(6)
    basis = np.eye(6, dtype=complex)
    vacuum = np.zeros(64, complex)
    vacuum[0] = 1
    for v in ((basis[:, 0]+basis[:, 1])/np.sqrt(2),
              (basis[:, 2]+basis[:, 3])/np.sqrt(2)):
        vacuum = creation(v, ann) @ vacuum
    number = [a.conj().T @ a for a in ann]
    raw = (np.eye(64)-number[0]) @ (np.eye(64)-number[2]) @ vacuum
    probability = norm(raw)**2
    base = raw/np.sqrt(probability)
    code = np.column_stack([ann[j].conj().T @ base for j in (0, 2)])
    ideal, actual = np.eye(6, dtype=complex), np.eye(6, dtype=complex)
    mode_distances, rank2_errors = [], []
    for (source, target), angle, phase in (((0, 4), .07, .04), ((2, 5), .09, -.03)):
        pair = (source, target)
        sx = np.array([[0., 1.], [1., 0.]])
        block_ideal = evolution(sx, np.pi/2)
        block_actual = evolution(sx, np.pi/2-angle) @ np.diag([np.exp(1j*phase), np.exp(-1j*phase)])
        ideal[np.ix_(pair, pair)] = block_ideal
        actual[np.ix_(pair, pair)] = block_actual
        rr = np.eye(6, dtype=complex)
        rr[np.ix_(pair, pair)] = block_ideal @ block_actual.conj().T
        distance = norm(block_actual[:, 0]-block_ideal[:, 0])
        mode_distances.append(distance)
        rank2_errors.append(abs(opnorm(lift(rr)-np.eye(64))-distance))
    rr = lift(ideal @ actual.conj().T)
    ua = lift(actual)
    out = rr @ ua @ code
    remainder = rr @ ua @ base
    unit_error = 0.
    for i, ti in enumerate((4, 5)):
        for j, tj in enumerate((4, 5)):
            expected = np.zeros((2, 2), complex)
            expected[i, j] = 1
            unit_error = max(unit_error, norm(out.conj().T @ ann[ti].conj().T @ ann[tj] @ out-expected))
    empty_error = max(norm(ann[i] @ remainder) for i in (4, 5))
    factor_error = norm(out-np.column_stack([-1j*ann[i].conj().T @ remainder for i in (4, 5)]))
    # Coherent reference witness, plus all matrix units above (the latter are complete).
    ref_out = (np.kron(out[:, 0], [1, 0])+np.kron(out[:, 1], [0, 1]))/np.sqrt(2)
    ref_expected = (np.kron(-1j*ann[4].conj().T @ remainder, [1, 0])+
                    np.kron(-1j*ann[5].conj().T @ remainder, [0, 1]))/np.sqrt(2)
    reference_error = norm(ref_out-ref_expected)
    pure = np.eye(6, dtype=complex)
    pure[0, 0], pure[1, 1] = np.exp(.21j), np.exp(-.21j)
    pure_error = abs(opnorm(lift(pure)-np.eye(64))-abs(np.exp(.21j)-1))
    full_distance = opnorm(rr-np.eye(64))
    assert abs(probability-.25) < TOL
    assert max(rank2_errors+[unit_error, empty_error, factor_error, reference_error, pure_error]) < TOL
    assert full_distance <= sum(mode_distances)+TOL
    return dict(scope="Full six-mode CAR, correlated negative-level vacuum, complete logical matrix units.",
                fock_dimension=64, empty_probability=probability,
                code_isometry_error=norm(code.conj().T @ code-np.eye(2)),
                rank_two_lift_norm_error=max(rank2_errors), pure_phase_pair_error=pure_error,
                two_rail_full_fock_distance=full_distance, two_rail_distance_upper=sum(mode_distances),
                all_matrix_units_error=unit_error, target_empty_remainder_error=empty_error,
                complete_output_factorization_error=factor_error, entangled_reference_error=reference_error)


def bounded_source():
    ann = annihilators(4)
    sx = np.array([[0., 1.], [1., 0.]])
    h = np.kron(np.diag([.8, -.8]), np.eye(2))
    hf = quadratic(h, ann)+1.6*np.eye(16)
    kernel = np.array([[.8, .25], [.25, .5]])
    j = quadratic(np.kron(kernel, sx), ann)
    n = 6
    a = np.diag(np.sqrt(np.arange(1, n)), 1)
    x = (a+a.T)/np.sqrt(2)
    p = 1j*(a.T-a)/np.sqrt(2)
    hs = np.diag(np.arange(n, dtype=float))
    ex, vx = np.linalg.eigh(x)
    gate = (vx*np.exp(-(ex-.4)**2/2)) @ vx.T
    source0 = np.eye(n)[:, 0]
    mean_gate = float(np.vdot(source0, gate @ source0).real)
    variance = norm((gate-mean_gate*np.eye(n)) @ source0)
    coupling, t = .7, .6
    field = np.kron(hf, np.eye(n))
    source = np.kron(np.eye(16), hs)
    interaction = coupling*np.kron(j, gate)
    total = field+source+interaction
    reference = np.kron(hf+coupling*mean_gate*j, np.eye(n))+source
    ua, ur = evolution(total, t), evolution(reference, t)
    embedding = np.kron(np.eye(16), source0[:, None])
    joint_error = opnorm((ua-ur) @ embedding)
    bound = t*coupling*opnorm(j)*variance
    field0 = np.zeros(16, complex)
    field0[13], field0[14] = 1/np.sqrt(2), np.exp(.4j)/np.sqrt(2)
    initial = np.kron(field0, source0)
    final = ua @ initial
    mean = lambda op, state: float(np.vdot(state, op @ state).real)
    changes = [mean(op, final)-mean(op, initial) for op in (field, source, interaction)]
    matrix = final.reshape(16, n)
    reduced_source = matrix.T @ matrix.conj()
    source_change = float(sum(abs(np.linalg.eigvalsh(reduced_source-np.outer(source0, source0))))/2)
    purity = float(np.trace(reduced_source @ reduced_source).real)
    pm = np.kron(np.eye(16), p)
    force = 1j*(total @ pm-pm @ total)
    free_force = 1j*(source @ pm-pm @ source)
    interaction_force = coupling*np.kron(j, 1j*(gate @ p-p @ gate))
    force_error = norm(force-free_force-interaction_force)
    ccr_defect = norm(x @ p-p @ x-1j*np.eye(n))
    assert joint_error <= bound+TOL
    assert abs(sum(changes)) < TOL and force_error < TOL
    assert source_change > 1e-4 and purity < .99999
    assert abs(changes[2]) > 1e-5 and ccr_defect > 1
    return dict(scope="Full four-mode CAR plus six-dimensional bounded source; not a scalar field or canonical CCR cutoff.",
                total_dimension=96, joint_isometry_error=joint_error, stationary_source_duhamel_bound=bound,
                field_energy_change=changes[0], source_energy_change=changes[1], interaction_energy_change=changes[2],
                total_energy_error=abs(sum(changes)), deleting_interaction_balance_error=abs(changes[0]+changes[1]),
                reduced_source_trace_distance=source_change, reduced_source_purity=purity,
                source_force_identity_error=force_error, source_interaction_force_norm=opnorm(interaction_force),
                finite_source_ccr_defect=ccr_defect)


def compare(actual, expected, path="results"):
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys(), path
        for k in actual:
            compare(actual[k], expected[k], path+"."+k)
    elif isinstance(actual, list):
        assert len(actual) == len(expected), path
        for i, (a, b) in enumerate(zip(actual, expected)):
            compare(a, b, path+f"[{i}]")
    elif isinstance(actual, float):
        assert np.isfinite(actual) and np.isclose(actual, expected, atol=1e-10, rtol=1e-10), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = dict(round=1129, schema=1, status="PASS", exact_budget=exact_budget(),
                  kernel_and_spectrum=kernel_spectrum(), full_car=car_code(), bounded_source=bounded_source())
    path = HERE/"results.json"
    if args.write_results:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
        print("PASS: four groups; results.json registered explicitly.")
    else:
        compare(result, json.loads(path.read_text(encoding="utf-8")))
        print("PASS: four groups; recorded results match; no files written.")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
