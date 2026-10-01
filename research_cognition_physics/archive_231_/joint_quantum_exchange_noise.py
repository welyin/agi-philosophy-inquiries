"""620: joint scalar/energy noise of the original 598 CAR mass sector.

Conditional prescribed backgrounds, not a dynamical Gauss/GR solution.
Preserves original complex Dirac/Majorana parameters and Nambu half weight.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quantum_response_matching as old

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_quantum_exchange_noise_results.json'
IDS = [24, 25, 30, 31]
BETA = 2.0
TIME = 2.3


def size(a):
    return float(np.max(np.abs(a)))


def occupation(b):
    e, v = np.linalg.eigh(b)
    return (v * (1 / (1 + np.exp(BETA * e)))) @ v.conj().T


def mean(f, a):
    return float(.5 * np.trace(f @ a).real)


def covariance(f, a, b):
    # Wick covariance for physical self-dual quadratic BdG matrices.
    return float(.5 * np.trace(f @ a @ (np.eye(len(f)) - f) @ b).real)


def cov_matrix(f, obs):
    return np.array([[covariance(f, a, b) for b in obs] for a in obs])


def frame_sources():
    rows = []
    for s in (old.S0, .9, 1.55):
        (b, j, _), _ = old.matrices(s)
        phi = np.array([0., old.HIGGS, 0., 0., s])
        fval = float(old.matter.original.F(phi))
        root = np.sqrt(fval)
        dr = -s / (6 * root)
        # Same Hamiltonian and same initial state on both sides.
        f = occupation(root * b)
        einstein = [b, root * j]
        jac = np.array([[root, 0.], [dr, 1.]])
        jordan = [root * b, dr * b + root * j]
        eps = 1e-5
        def numerator(z):
            phi1 = np.array([0., old.HIGGS, 0., 0., z])
            return np.sqrt(old.matter.original.F(phi1)) * old.matrices(z)[0][0]
        direct = (numerator(s + eps) - numerator(s - eps)) / (2 * eps)
        ne = cov_matrix(f, einstein)
        nj = cov_matrix(f, jordan)
        err = size(nj - jac @ ne @ jac.T)
        derivative_error = size(direct - jordan[1])
        assert max(err, derivative_error) < 3e-11
        # Drop only the cross term, retaining both marginal contributions.
        independent = dr**2 * ne[0, 0] + ne[1, 1]
        cross = 2 * dr * ne[0, 1]
        assert abs(independent - nj[1, 1]) > 1e-5
        # At fixed Higgs the Jordan source is purely original nu pairing.
        neutral_bdg = IDS + [i + 32 for i in IDS]
        outside = np.ones(64, dtype=bool); outside[neutral_bdg] = False
        assert size(jordan[1][outside]) < 1e-14
        rows.append(dict(s=float(s), F=fval, same_state_BdG_dimension=64,
            Einstein_source_covariance=ne.tolist(), Jordan_source_covariance=nj.tolist(),
            source_jacobian=jac.tolist(), covariance_pullback_error=err,
            direct_mass_derivative_error=derivative_error,
            omitted_mixed_contribution=float(cross),
            independent_source_prediction=float(independent),
            correct_Jordan_scalar_variance=float(nj[1, 1]),
            spurious_charged_source_norm=size(jordan[1][outside])))
    return dict(rows=rows, all_32_original_modes_retained=True,
                state_held_identical_under_parameter_pullback=True)


def path(t):
    return old.S0 + .28 * np.sin(.7 * t), .196 * np.cos(.7 * t)


def exponential(h, dt):
    e, v = np.linalg.eigh(h)
    return (v * np.exp(-1j * dt * e)) @ v.conj().T


def evolve(steps, ids=None, fock=False):
    def matrices(s):
        mats, pairs = old.matrices(s, ids)
        return [old.fock(p) for p in pairs[:2]] if fock else mats[:2]
    b0, _ = matrices(path(0)[0])
    u = np.eye(len(b0), dtype=complex)
    work = np.zeros_like(b0)
    dt = TIME / steps
    for k in range(steps):
        t = (k + .5) * dt
        s, ds = path(t)
        b, j = matrices(s)
        half = exponential(b, dt / 2)
        middle = half @ u
        work += dt * ds * (middle.conj().T @ j @ middle)
        u = half @ middle
    end, jend = matrices(path(TIME)[0])
    end_heis = u.conj().T @ end @ u
    delta = end_heis - b0
    je = u.conj().T @ jend @ u
    return b0, delta, work, je, u


def exchange_check():
    rows = []
    for steps in (64, 128, 256):
        b0, delta, work, je, u = evolve(steps)
        f = occupation(b0)
        n = cov_matrix(f, [delta, work])
        error = size(delta - work)
        residual = covariance(f, delta - work, delta - work)
        independent = n[0, 0] + n[1, 1]
        comm = size(b0 @ old.matrices(old.S0)[0][1] - old.matrices(old.S0)[0][1] @ b0)
        assert comm > 1e-3
        assert size(u.conj().T @ u - np.eye(64)) < 2e-12
        assert np.linalg.eigvalsh(n).min() > -1e-12
        assert n[0, 0] > 1e-4 and independent > 1e-3
        assert residual >= 0 and residual < 1e-8
        obs = [b0, delta, work, je]
        assert np.linalg.eigvalsh(cov_matrix(f, obs)).min() > -1e-12
        rows.append(dict(steps=steps, operator_balance_error=error,
            energy_change=mean(f, delta), integrated_scalar_work=mean(f, work),
            exchange_noise=n.tolist(), joint_constraint_variance=residual,
            independent_noise_constraint_variance=float(independent),
            noncommuting_mass_source_norm=comm,
            full_joint_covariance_min_eigenvalue=float(np.linalg.eigvalsh(cov_matrix(f, obs)).min())))
    ratios = [rows[i]['operator_balance_error'] / rows[i+1]['operator_balance_error'] for i in range(2)]
    assert min(ratios) > 3.8 and max(ratios) < 4.2
    assert rows[-1]['operator_balance_error'] < 3e-6
    return dict(path=dict(s_initial=float(old.S0), amplitude=.28, frequency=.7, duration=TIME),
        initial_beta=BETA, rows=rows, second_order_convergence_ratios=ratios,
        numerical_errors_are_not_rigorous_interval_bounds=True,
        integrated_work_is_a_Heisenberg_operator_not_a_TPM_distribution=True)


def independent_fock_check():
    steps = 128
    b0, delta, work, je, u = evolve(steps, IDS)
    hf, df, wf, jf, uf = evolve(steps, IDS, fock=True)
    f = occupation(b0)
    _, _, _, rho = old.thermal(hf, BETA)
    obs_b = [b0, delta, work, je]
    obs_f = [hf, df, wf, jf]
    mb = np.array([mean(f, a) for a in obs_b])
    mf = np.array([np.trace(rho @ a).real for a in obs_f])
    nf = np.array([[np.trace(rho @ a @ b).real - ma * mbb
                    for b, mbb in zip(obs_f, mf)] for a, ma in zip(obs_f, mf)])
    nb = cov_matrix(f, obs_b)
    err_mean = size(mb - mf); err_cov = size(nb - nf)
    assert max(err_mean, err_cov) < 2e-12
    assert size(df - wf) < 1e-5
    variance = float(np.trace(rho @ (df - wf) @ (df - wf)).real - np.trace(rho @ (df - wf)).real**2)
    return dict(CAR_indices=IDS, Fock_dimension=16, steps=steps,
        mean_error=err_mean, symmetrized_covariance_error=err_cov,
        exact_Fock_covariance=nf.tolist(), work_identity_matrix_error=size(df-wf),
        residual_constraint_variance=variance,
        full_original_Gauss_state_not_used=True)


def run():
    a = frame_sources(); b = exchange_check(); c = independent_fock_check()
    dependencies = ('joint_quantum_response_matching.py','joint_fermion_gauss_completion.py',
                    'research_note_590.md','research_note_600.md','research_note_602.md',
                    'research_note_603.md','research_note_618.md','research_note_619.md')
    return dict(round=620, tests_run=3, failures=0, errors=0,
        frame_sources=a, driven_exchange=b, independent_Fock=c,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in dependencies},
        scope=dict(original_quadratic_mass_sector=True, same_initial_state_and_source_retained=True,
            prescribed_background_not_a_closed_universe_solution=True,
            no_Gauss_continuum_matching_or_GR_completion=True,
            covariance_not_a_complete_quantum_process=True,
            covariant_Ward_extension_requires_renormalized_anomaly_free_source_identity=True,
            established_Ward_and_Wick_methods_not_claimed_as_new_fundamental_theorems=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args(); result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(result, ensure_ascii=False))
