"""1137: five mathematical verification groups, not scientific trials.

Default execution reads and compares results.json; only --write-results writes it.
Frozen preflight computations are reused after SHA256 checks and loaded without
creating bytecode caches. Finite diagnostics do not replace the continuous proof.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import types
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FROZEN = {
    '1133/proof.md': '476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
    '1134/proof.md': '4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666',
    '1135/proof.md': 'f3577777e8983a02d02bf66ef5f869b207ac4f94df1610202c1ce6e995e17068',
    '1136/proof.md': '360c450f85f40798c6fb9524b00124364b604434d814c770503e0e97e2196462',
    '_shared/receiver1137_preflight_author.md': '24c1903ddaedd651ca2ad137158c7c5511b92069c20a5e2bb448ed14b30b28a7',
    '_shared/receiver1137_preflight_author.py': 'b5c750404c2cfe793356811b244d9403859f7501ad3bb8a49bce895b5effa90d',
    '_shared/receiver1137_preflight_author.json': '71606653ed4dd237a76c538725d9a525ccb6880bc759e82edfefecc7b4d255d5',
    '_shared/receiver1137_preflight_review.md': '0beca63835b664196ac008aa719f8d7c6359448be2e123b6f2841852efa8e3ac',
    '_shared/receiver1137_preflight_review.py': 'ba0ee1df1fa3e3700aae4c413f5a907c2a51e6e5de80249c655a0a396e896a30',
    '_shared/receiver1137_preflight_review.json': '5b123bc19c1bcccd0405f3e7ea4a1584b2784442dbff65cec2071b62cc00ee19',
}


def load_frozen(name):
    path = ROOT / '_shared' / ('receiver1137_preflight_' + name + '.py')
    mod = types.ModuleType('frozen_' + name)
    mod.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), 'exec'), mod.__dict__)
    return mod


def sqrt_hermitian(a):
    values, vectors = np.linalg.eigh(a)
    assert values.min() > 0
    return (vectors * np.sqrt(values)) @ vectors.conj().T


def noncommuting_mass_and_terminal():
    """Noncommuting internal mass, its pointer conjugate, and full Kraus output.

    Direct spectral evaluation is independent of the analytic semigroup proof.
    The terminal finite model retains a two-dimensional passive reference.
    """
    rng = np.random.default_rng(7371137)
    d = 10
    z = np.kron(np.eye(5), np.diag([1., -1.]))
    raw = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    v = (raw + raw.conj().T)/2
    v = 3*v/np.linalg.norm(v, 2)
    a = np.diag(np.linspace(20., 34., d)) + v
    b = z @ a @ z
    c = a-b
    delta = float(np.linalg.norm(c, 2))
    mass = float(min(np.linalg.eigvalsh(a).min(), np.linalg.eigvalsh(b).min()))
    noncomm = float(np.linalg.norm(a@b-b@a, 2))
    assert noncomm > 1
    rows = []
    for k in [.125, 2., 11., 80.]:
        ha = sqrt_hermitian(a@a+k*k*np.eye(d))
        hb = sqrt_hermitian(b@b+k*k*np.eye(d))
        diff = ha-hb
        ga, gb = ha-a, hb-b
        rhs = -(ga@c+c@gb)
        wa, ua = np.linalg.eigh(ha)
        wb, ub = np.linalg.eigh(hb)
        # Solve the two-sided equation in two independent spectral bases.
        integral = ua @ ((ua.conj().T@rhs@ub)/(wa[:, None]+wb[None, :])) @ ub.conj().T
        sylv = float(np.linalg.norm(diff-c-integral, 2))
        derivative = -k*np.linalg.solve(ha, diff)@np.linalg.inv(hb)
        step = 1e-4
        plus = sqrt_hermitian(a@a+(k+step)**2*np.eye(d))-sqrt_hermitian(b@b+(k+step)**2*np.eye(d))
        minus = sqrt_hermitian(a@a+(k-step)**2*np.eye(d))-sqrt_hermitian(b@b+(k-step)**2*np.eye(d))
        fd = float(np.linalg.norm((plus-minus)/(2*step)-derivative, 2))
        conjugacy = float(np.linalg.norm(hb-z@ha@z, 2))
        assert sylv < 1e-10 and fd < 1e-8 and conjugacy < 1e-10
        assert np.linalg.norm(diff, 2) <= 2*delta+1e-10
        assert np.linalg.norm(derivative, 2) <= delta/mass+1e-10
        rows.append(dict(momentum=k,negative_Sylvester_residual=sylv,
                         derivative_finite_difference_residual=fd,
                         conjugacy_residual=conjugacy,
                         mass_difference_norm=float(np.linalg.norm(diff, 2)),
                         derivative_norm=float(np.linalg.norm(derivative, 2))))
    # Two position outcomes x two bare pointer outcomes, with the complete
    # internal state and passive reference retained. This checks the square-root
    # tail coefficient, not an unproved conversion of tail probability to norm.
    positions = np.array([-.3, -.12, -.07, 0., .02, .11, .25])
    ell = .1
    amplitudes = np.exp(-positions**2/(4*.04**2))
    amplitudes /= np.linalg.norm(amplitudes)
    internal = rng.normal(size=(d, 2))+1j*rng.normal(size=(d, 2))
    internal /= np.linalg.norm(internal)
    state = amplitudes[:, None, None]*internal[None, :, :]
    u = np.abs(positions)/ell-1
    ramp = np.zeros_like(u)
    ramp[u >= 1] = 1
    middle = (u > 0) & (u < 1)
    left = np.exp(-1/u[middle]); right = np.exp(-1/(1-u[middle]))
    ramp[middle] = left/(left+right)
    lin, lout = np.cos(math.pi*ramp/2), np.sin(math.pi*ramp/2)
    projectors = [(np.eye(d)+z)/2, (np.eye(d)-z)/2]
    out = np.zeros((2, 2, len(positions), d, 2), complex)
    target = np.zeros_like(out)
    for j, p in enumerate(projectors):
        branch = np.einsum('ab,ybr->yar', p, state)
        out[0, j] = lin[:, None, None]*branch
        out[1, j] = lout[:, None, None]*branch
        target[0, j] = branch
    tail_norm = float(np.linalg.norm(lout[:, None, None]*state))
    distance = float(np.linalg.norm(out-target))
    rms = float(np.sqrt(np.dot(amplitudes**2, positions**2)))
    gram = float(abs(np.linalg.norm(out)**2-1))
    assert gram < 1e-13
    assert distance <= math.sqrt(2)*tail_norm+1e-14
    assert tail_norm <= rms/ell
    return dict(internal_dimension=d,passive_reference_dimension=2,
                mass_lower=mass,mass_conjugacy_delta=delta,
                noncommuting_mass_commutator_norm=noncomm,rows=rows,
                full_four_Kraus=dict(output_norm_residual=gram,
                     full_poststate_distance=distance,outcome_tail_norm=tail_norm,
                     squared_tail_probability=tail_norm**2,position_RMS=rms,
                     square_root_tail_bound=math.sqrt(2)*tail_norm),
                continuous_boost_and_field_proof='proof.md, sections 4 and 8; diagnostic is finite only')


def calculate():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert hashes == FROZEN
    author, review = load_frozen('author'), load_frozen('review')
    first = author.calculate()
    budgets = review.exact_budget()
    assert budgets['fractions']['full_task_error'] == first['task_rows']['code']['fractions']['complete_task_error']
    assert F(budgets['fractions']['full_task_error']) < F(1, 100)
    assert F(budgets['fractions']['direction_probability_gap_lower']) > F(2, 5)
    assert F(budgets['fractions']['complete_boost_joint_instrument_error']) < F(1, 100)
    assert F(first['task_rows']['far_local_menu']['fractions']['finite_causal_error']) < F(1, 10**6)
    # Independent selected-axis test includes a counterexample to a fixed x axis.
    axis, car = review.spinor_axis(), review.full_car()
    assert abs(car['logical_effect_commutator_norm']-.5) < 1e-12
    return dict(round=1137,status='PASS',new_axioms=0,science_experiments=0,
        mathematical_check_groups=5,stage='active',groups={
          '1_frozen_identity':dict(count=len(hashes),sha256=hashes),
          '2_exact_budgets':dict(common=first['common'],task_rows=first['task_rows'],
                               joint_direction_and_boost=budgets),
          '3_normal_Haar_and_derivative_domains':dict(angular=first['angular'],
                               independent_quadrature=first['Haar_quadrature']),
          '4_complete_normal_CAR_and_instruments':dict(selected_axis=axis,full_CAR=car),
          '5_noncommuting_mass_and_terminal_poststate':noncommuting_mass_and_terminal()},
        scope=dict(material_contract='M1137',source_selector_b=0,T=4,
          old_D_dimension=288,unknown_Q_and_passive_R='uniform full output',
          source_Z_O_Y_and_all_field_remainders='retained',
          settings=['N_Z','N_X'],continuation='exp(-i pi Y_g / 4)',
          rapidity='log(2)',boost_event=[5,0,0,3],position_resolution=.001,
          old_Mr_at_T4_claimed=False,b1_extended_to_T4=False,
          arbitrary_remote_source_controls=False,local_stress_or_shell_flux_claimed=False,
          terminal_detector_autonomous_energy_closed=False))


def compare(a, b, path=()):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a: compare(a[key], b[key], path+(key,))
    elif isinstance(a, list):
        assert len(a) == len(b)
        for index, (x, y) in enumerate(zip(a, b)): compare(x, y, path+(str(index),))
    elif isinstance(a, float):
        # Tiny physical numbers must not disappear under a residual tolerance.
        # Exact-budget fractions are also compared as exact strings above.
        physical = ('2_exact_budgets' in path or '3_normal_Haar_and_derivative_domains' in path)
        rel_tol, abs_tol = (1e-8, 0.) if physical else (1e-7, 1e-10)
        assert math.isclose(a, b, rel_tol=rel_tol, abs_tol=abs_tol), (path, a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    data = calculate()
    result_path = HERE/'results.json'
    if args.write_results:
        result_path.write_text(json.dumps(data, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        compare(data, json.loads(result_path.read_text(encoding='utf-8')))
    print('1137 PASS: 5 mathematical groups; 10 frozen sources; normal CAR full records, continuation, rest/direction/boost/causal budgets. New axioms 0; science trials 0; stage active.')
