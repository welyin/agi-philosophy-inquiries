"""555: common energy-moment bounds for the interacting finite reference menu.

Polynomial wave functions and independent phase-space quadrature check the
identities. No discretized-time evolution, continuum limit, or source
preparation is claimed. The mathematical all-time statement is proved in note.
"""
import argparse
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path
import numpy as np
import joint_matter_reference_quantum_readout as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_matter_energy_moment_control_results.json'


def const(c, d):
    return {(0,) * d: complex(c)} if c else {}


def add(*polys):
    out = {}
    for p in polys:
        for k, c in p.items():
            out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c != 0}


def scale(p, c):
    return {k: c * x for k, x in p.items() if c * x != 0}


def mul(p, q):
    out = {}
    for k, c in p.items():
        for l, e in q.items():
            key = tuple(x + y for x, y in zip(k, l))
            out[key] = out.get(key, 0) + c * e
    return {k: c for k, c in out.items() if c != 0}


def deriv(p, i):
    out = {}
    for k, c in p.items():
        if k[i]:
            key = list(k); key[i] -= 1
            out[tuple(key)] = c * k[i]
    return out


def var(i, d):
    key = [0] * d; key[i] = 1
    return {tuple(key): 1.}


def conj(p):
    return {k: np.conj(c) for k, c in p.items()}


@lru_cache(None)
def normal_moment(n, variance):
    if n % 2:
        return 0.
    return math.prod(range(1, n, 2)) * variance ** (n // 2)


def expect(p, variance):
    return sum(c * math.prod(normal_moment(n, variance) for n in key)
               for key, c in p.items())


def real(x):
    assert abs(complex(x).imag) < 1e-8 * max(1., abs(x))
    return float(complex(x).real)


def norm2(p, variance):
    return real(expect(mul(conj(p), p), variance))


def quadratic(q, matrix):
    return add(*(scale(mul(q[i], q[j]), matrix[i, j])
                 for i in range(len(q)) for j in range(len(q))
                 if matrix[i, j] != 0))


def matter():
    row = json.loads((HERE / 'joint_singlet_common_mass_rg_results.json')
                     .read_text('utf8'))['examples'][2]['state']
    L = np.array([[row['lambda_H'], row['p']], [row['p'], row['lambda_s']]])
    C = .25 * np.array([row['x'], row['y']])
    return L, C, np.linalg.solve(L, C)


def constants(N, v, K, L, C, u, E1, E2, omega, r=1.):
    ell = float(np.linalg.eigvalsh(L)[0])
    alpha = max(0., 3 * L[0, 0] + L[0, 1], 3 * L[1, 1] + L[0, 1])
    A = 2 * alpha / (ell * r)
    B = v * 2 * np.trace(K) + v * alpha * (N * sum(u) + N * r) + v * N * sum(abs(C))
    a = 1 / (2 * v); V = N * v
    graph = E2 + a * A * E1 + a * B
    S = N * sum(u) + N * r + 2 * E1 / (v * ell * r)
    tq, tp = 1 / (2 * omega), omega / 2
    noise = (8 * E1 / N**2 * (tq * np.linalg.norm(K, 2) / v + tp / v**3)
             + 2 / N**2 * (tq**2 * np.trace(K @ K) + tp**2 * N / v**4)
             - np.trace(K) / (N**2 * v**2))
    records = [S / N + tq / N] * 2 + [8 * graph / V**2 + noise] * 2
    return dict(A=float(A), B=float(B), graph=float(graph),
                position_sum=float(S), intrinsic_square=float(8 * graph / V**2),
                record_second_moment_bounds=[float(x) for x in records])


def hermite_wave(n, x, omega, d):
    h0 = const(1, d)
    if n == 0:
        return h0
    h1 = scale(x, 2 * math.sqrt(omega))
    for k in range(1, n):
        h0, h1 = h1, add(scale(mul(x, h1), 2 * math.sqrt(omega)), scale(h0, -2 * k))
    return scale(h1, 1 / math.sqrt(2**n * math.factorial(n)))


def witness(N, v, K, occupations, omega=1.3, mu=None, momentum=None):
    L, C, u = matter(); d = 2 * N; variance = 1 / (2 * omega)
    mu = np.zeros(d) if mu is None else np.asarray(mu)
    momentum = np.zeros(d) if momentum is None else np.asarray(momentum)
    xs = [var(i, d) for i in range(d)]
    qs = [add(xs[i], const(mu[i], d)) for i in range(d)]
    P = const(1, d)
    for i, n in enumerate(occupations):
        P = mul(P, hermite_wave(n, xs[i], omega, d))
    assert abs(norm2(P, variance) - 1) < 1e-11
    gradients, seconds = [], []
    for i in range(d):
        g = add(scale(xs[i], -omega), const(1j * momentum[i], d))
        dp = deriv(P, i)
        gradients.append(add(dp, mul(g, P)))
        seconds.append(add(deriv(dp, i), scale(mul(g, dp), 2),
                           mul(add(mul(g, g), const(-omega, d)), P)))
    Wgrad = add(*(scale(quadratic(qs[f*N:(f+1)*N], K), v / 2) for f in (0, 1)))
    Wpot = {}
    for j in range(N):
        delta = [add(mul(qs[f*N+j], qs[f*N+j]), const(-u[f], d)) for f in (0, 1)]
        Wpot = add(Wpot, scale(quadratic(delta, L), v / 4))
    W = add(Wgrad, Wpot)
    Tpsi = scale(add(*seconds), -1 / (2 * v)); Wpsi = mul(W, P)
    Hpsi = add(Tpsi, Wpsi)
    E1 = real(expect(mul(conj(P), Hpsi), variance)); E2 = norm2(Hpsi, variance)
    T2, W2 = norm2(Tpsi, variance), norm2(Wpsi, variance)
    lapW = add(*(deriv(deriv(W, i), i) for i in range(d)))
    grad_weight = sum(real(expect(mul(W, mul(conj(g), g)), variance)) for g in gradients)
    lap_mean = real(expect(mul(mul(conj(P), P), lapW), variance))
    rhs = T2 + W2 + grad_weight / v - lap_mean / (2 * v)
    residual = abs(E2 - rhs) / max(1., E2)
    assert residual < 2e-10
    bounds = constants(N, v, K, L, C, u, E1, E2, omega)
    assert T2 + W2 <= bounds['graph'] + 1e-8 * max(1., E2)
    means, source2, records2 = [], [], []
    tq, tp = 1 / (2 * omega), omega / 2
    for f in (0, 1):
        Q = scale(add(*qs[f*N:(f+1)*N]), 1 / N)
        Qpsi = mul(Q, P)
        means.append(real(expect(mul(conj(P), Qpsi), variance)))
        source2.append(norm2(Qpsi, variance))
        records2.append(source2[-1] + tq / N)
    for f in (0, 1):
        part = qs[f*N:(f+1)*N]
        G = scale(quadratic(part, K), 1 / N)
        Fpsi = add(mul(G, P), scale(add(*seconds[f*N:(f+1)*N]), 1 / (N * v*v)))
        mean = real(expect(mul(conj(P), Fpsi), variance)); F2 = norm2(Fpsi, variance)
        qK2q = real(expect(mul(mul(conj(P), P), quadratic(part, K @ K)), variance))
        p2 = sum(norm2(g, variance) for g in gradients[f*N:(f+1)*N])
        noise = (4 * tq * qK2q / N**2 + 4 * tp * p2 / (N**2 * v**4)
                 + 2 * (tq*tq * np.trace(K @ K) + tp*tp * N / v**4) / N**2
                 - np.trace(K) / (N**2 * v*v))
        means.append(mean); source2.append(F2); records2.append(float(F2 + noise))
        assert F2 <= bounds['intrinsic_square'] + 1e-8 * max(1., E2)
    assert all(x <= b + 1e-8 * max(1., b)
               for x, b in zip(records2, bounds['record_second_moment_bounds']))
    return dict(N=N, v=v, omega=omega, occupations=list(occupations),
                mean_position=mu.tolist(), mean_momentum=momentum.tolist(),
                E1=E1, E2=E2, T2_plus_W2=T2+W2, identity_relative_residual=residual,
                menu_means=means, intrinsic_second_moments=source2,
                actual_record_second_moments=records2, bounds=bounds)


def husimi_quadrature_moment(qpower, ppower, n, omega):
    # Actual displaced-vacuum POVM density: e^-x x^n/n! dx dtheta/(2pi).
    # This integrates its positive output density, not a positive-Wigner ansatz.
    if (qpower + ppower) % 2 or qpower % 2 or ppower % 2:
        return 0.
    x, w = np.polynomial.laguerre.laggauss(14)
    theta = 2 * np.pi * np.arange(24) / 24
    angular = np.mean(np.cos(theta)**qpower * np.sin(theta)**ppower)
    radial = np.sum(w * x**(n+(qpower+ppower)//2)) / math.factorial(n)
    return float((2/omega)**(qpower/2) * (2*omega)**(ppower/2) * angular * radial)


def independent_fock_records(N, v, K, occupations, omega):
    d = 4*N; z = [var(i, d) for i in range(d)]
    tq, tp = 1/(2*omega), omega/2
    menu = [scale(add(*z[:N]), 1/N), scale(add(*z[N:2*N]), 1/N)]
    for f in (0, 1):
        fpoly = add(scale(quadratic(z[f*N:(f+1)*N], K), 1/N),
                    scale(add(*(mul(z[i], z[i]) for i in range((f+2)*N, (f+3)*N))), -1/(N*v*v)),
                    const(-tq*np.trace(K)/N + tp/v**2, d))
        menu.append(fpoly)
    cache = {}
    def integrate(poly):
        answer = 0j
        for key, c in poly.items():
            factor = c
            for i, n in enumerate(occupations):
                args = (key[i], key[i+2*N], n)
                if args not in cache:
                    cache[args] = husimi_quadrature_moment(*args, omega)
                factor *= cache[args]
            answer += factor
        return real(answer)
    return [integrate(p) for p in menu], [integrate(mul(p, p)) for p in menu]


def run():
    tests = []; L, C, u = matter()
    assert np.linalg.eigvalsh(L)[0] > 0 and min(u) > 0
    rng = np.random.default_rng(555)
    majorant_gap = -math.inf
    for N, v, K in ((1, 1., np.zeros((1, 1))), (2, .7, np.array([[1., -1.], [-1., 1.]]))):
        cb = constants(N, v, K, L, C, u, 0., 0., 1.3)
        for q in rng.normal(size=(60, 2, N)) * 4:
            delta = (q*q).T-u
            W = v*sum(x@L@x for x in delta)/4 + v*sum(x@K@x for x in q)/2
            alpha = np.array([3*L[0, 0]+L[0, 1], 3*L[1, 1]+L[0, 1]])
            lap = v*2*np.trace(K) + v*sum(alpha[f]*sum(q[f]**2)-N*C[f] for f in (0, 1))
            gap = lap - cb['A']*W - cb['B']
            majorant_gap = max(majorant_gap, float(gap))
            assert gap < 1e-12
    tests.append('same_positive_interacting_potential_and_explicit_laplacian_majorant')
    rows = []
    K2 = np.array([[1., -1.], [-1., 1.]])
    for N, v, K, occ in ((1, 1., np.zeros((1, 1)), [1, 2]),
                         (2, .7, K2, [1, 0, 2, 1]),
                         (2, 1., K2, [0, 0, 0, 0])):
        rows.append(witness(N, v, K, occ))
    xyz, D = old.geometry(); K8 = sum(d.T@d for d in D)
    mu = np.r_[np.full(8, math.sqrt(u[0])), math.sqrt(u[1])+xyz[:, 0]+xyz[:, 0]*xyz[:, 2]]
    mom = np.r_[1+xyz[:, 1], np.zeros(8)]
    rows.append(witness(8, 1., K8, [0]*16, omega=2., mu=mu, momentum=mom))
    tests.append('graph_norm_identity_and_full_menu_operator_bounds_on_nonGaussian_and_eight_cell_sources')
    independent_errors = []
    for row, K in zip(rows[:3], (np.zeros((1, 1)), K2, K2)):
        m, s = independent_fock_records(row['N'], row['v'], K, row['occupations'], row['omega'])
        err = max(np.max(abs(np.array(m)-row['menu_means'])),
                  np.max(abs(np.array(s)-row['actual_record_second_moments'])))
        independent_errors.append(float(err))
        assert err < 1e-8
    # Cross-check the displaced eight-cell case against independently inherited
    # Gaussian phase-space moments, including the meter subtraction only.
    N=8; omega=2.; _, ls, As=old.menu(N, 1., D)
    G=np.diag([1/(2*omega)]*16+[omega/2]*16)
    m, cov=old.gaussian_moments(np.r_[mu,mom],2*G,ls,As)
    m-=np.array([np.trace(A@G) for A in As])
    err=float(np.max(abs(np.diag(cov)+m*m-rows[-1]['actual_record_second_moments'])))
    assert err<1e-8; independent_errors.append(err)
    tests.append('independent_positive_Husimi_density_quadrature_and_Moyal_noise_correction')
    certificates=[]
    for row in rows:
        budgets=np.array(row['bounds']['record_second_moment_bounds'])
        actual=np.array(row['actual_record_second_moments'])
        radius=math.sqrt(sum(budgets)/.05)
        assert sum(actual)/radius**2<=.05+1e-12
        certificates.append(dict(N=row['N'],raw_menu_radius=radius,
                                 certified_single_read_failure_at_most=.05,
                                 witnessed_Markov_bound=float(sum(actual)/radius**2)))
    tests.append('one_common_energy_budget_bounds_all_actual_records_and_finite_confidence')
    base=witness(2,1.,K2,[0]*4)
    tail=[]
    for M in (2.,4.,16.,64.):
        p=np.array([M/math.sqrt(2),-M/math.sqrt(2),0.,0.])
        displaced=witness(2,1.,K2,[0]*4,momentum=p)
        weight=1/(M*M)
        E1=(1-weight)*base['E1']+weight*displaced['E1']
        E2=(1-weight)*base['E2']+weight*displaced['E2']
        rec=(1-weight)*base['actual_record_second_moments'][2]+weight*displaced['actual_record_second_moments'][2]
        assert abs(E1-base['E1']-.5)<1e-10
        tail.append(dict(M=M,E1=E1,E2=E2,actual_A_record_second_moment=rec))
    assert all(tail[i+1]['E2']>tail[i]['E2'] and
               tail[i+1]['actual_A_record_second_moment']>tail[i]['actual_A_record_second_moment']
               for i in range(3))
    tests.append('fixed_mean_energy_tail_failure_is_excluded_by_second_energy_moment_budget')
    # With W=0 and N=1, position translations keep every kinetic-energy moment.
    omega=1.3
    free=[dict(position_shift=M,E1=omega/2,E2=omega**2/2,
               actual_Q_record_second_moment=M*M+1/omega) for M in (0.,4.,16.,64.)]
    assert len({x['E1'] for x in free})==1 and len({x['E2'] for x in free})==1
    assert free[-1]['actual_Q_record_second_moment']>4000
    tests.append('without_confinement_energy_moments_do_not_control_field_value_records')
    dependencies=('joint_singlet_common_mass_rg_results.json','joint_matter_reference_quantum_readout.py',
                  'joint_matter_reference_quantum_readout_results.json')
    return dict(round=555,tests_run=len(tests),failures=0,errors=0,checks=tests,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
                L=L.tolist(),C=C.tolist(),u=u.tolist(),maximum_sampled_majorant_gap=majorant_gap,
                witness_sources=rows,independent_record_moment_errors=independent_errors,
                finite_record_certificates=certificates,tail_counterfamily=tail,
                nonconfining_counterfamily=free,
                scope=dict(finite_cells_only=True,hbar_in_numerics=1,
                           arbitrary_source_theorem_proved_in_note=True,
                           no_Gaussian_preservation_assumed=True,
                           all_time_bound_is_for_unmeasured_static_H_evolution=True,
                           instrument_and_preparation_are_inputs=True,
                           no_continuum_or_volume_uniform_limit_claimed=True,
                           no_spacetime_coordinate_or_gravity_generation_claimed=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=result['round'],tests=result['tests_run'],
                         max_record_residual=max(result['independent_record_moment_errors']),
                         max_graph_identity_residual=max(r['identity_relative_residual']
                                                         for r in result['witness_sources']))))
