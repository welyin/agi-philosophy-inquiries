"""Independent 1169 review: affine Hamilton flow and polynomial FK projectors.

Only writes review_checks.json with --write-once.  Default is read-only.
The author check.py is not imported or executed by this independent algorithm.
The separate author default run was observed to exit 0 before this review.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import re
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
EXPECTED = {
    'archive_1161_/research_note_1169.md': '49ef3dbe394569a3db113f1d358156753a3ad81274313af8c77e31192b4b008a',
    'archive_1161_/1169/admission.md': '648a3cbf3ca6c76811d6b05cddc24fa6d852f3110446343d5294036ccefe09d0',
    'archive_1161_/1169/proof.md': '33f4384f5a851a5b8e35f816e577fc4101eefcf8365d9024bab0a70e739d5871',
    'archive_1161_/1169/check.py': '67d8f3f2fb886de162e5e906035239a93a4c6fe94fc409bc08adc5e5e61f7c16',
    'archive_1161_/1169/results.json': 'ffb6d6ae5199a31cee25d54e2eb89f062a31da18de1281354f8c2690c0dcd6f5',
    'archive_1161_/1169/NEXT.md': 'adb6c24d98afb449bd176560a4ec09f3079512760b6b644fd01f8b3e5c729284',
}
M, MU, K, A, SIGMA = 1e8, 1e12, 10., 100., .01
T = 32*math.pi+math.pi/2
WIDTH, THRESHOLD = .001, .032


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def coeff(r):
    s = r*r+A*A
    return s**-.5, -r*s**-1.5, [(2*r*r-A*A)*s**-2.5, -s**-1.5, -s**-1.5]


def affine_flow(radius, time, mass=M, coupling=K):
    """Exponentiate the 5x5 affine Hamilton generator directly in long double.

    Coordinates are deviations from the initial means, [X,Y,PX,PY,1].
    This uses neither the author's COM-relative conversion nor f/g functions.
    """
    _, grad, hs = coeff(radius)
    initial = np.diag(np.array([SIGMA**2, SIGMA**2,
                               1/(4*SIGMA**2), 1/(4*SIGMA**2)], dtype=np.longdouble))
    means, covs, last = [], [], []
    for axis, h in enumerate(hs):
        generator = np.zeros((5, 5), dtype=np.longdouble)
        generator[0, 2], generator[1, 3] = 1/mass, 1/MU
        v = np.longdouble(coupling*mass*h)
        generator[2:4, 0:2] = [[v, -v], [-v, v]]
        generator[2, 4] = -coupling*mass*grad if axis == 0 else 0
        generator[3, 4] = -generator[2, 4]
        z = generator*np.longdouble(time)
        term = np.eye(5, dtype=np.longdouble)
        exp = term.copy()
        for n in range(1, 37):
            term = term@z/n
            exp += term
        last.append(float(np.max(np.abs(term))))
        means.append(np.array(exp[:4, 4], dtype=float))
        covs.append(np.array(exp[:4, :4]@initial@exp[:4, :4].T, dtype=float))
    assert max(last) < 1e-80
    return np.array(means), np.array(covs), max(last)


def independent_phase(radius):
    """Composite Simpson on direct affine trajectories, not sinh integrals."""
    u, grad, h = coeff(radius)
    intervals = 512
    rates = []
    for n in range(intervals+1):
        mean, _, _ = affine_flow(radius, T*n/intervals)
        rel = mean[:, 1]-mean[:, 0]
        ux = u+grad*rel[0]+.5*np.dot(h, rel*rel)
        rates.append(1-K*ux-.5*np.sum(mean[:, 2]**2)/M**2)
    vals = np.array(rates)
    integral = T/(3*intervals)*(vals[0]+vals[-1]+4*vals[1:-1:2].sum()+2*vals[2:-1:2].sum())
    return float(integral)


def independent_bound(radius):
    d, s, p, sp = (.052 if radius == 1000 else .014), .0143, (105000. if radius == 1000 else 27000.), 52.
    u, b, h = coeff(radius)
    hs = max(map(abs, h))
    soft = 4*T*K*(M+3)/(radius-1)**4*(d+105**(1/6)*s)**3
    fluctuations = K*((abs(b)+hs*d)*s+hs*math.sqrt(15)*s*s/2)
    fluctuations += (2*p*sp+math.sqrt(15)*sp*sp)/(2*M*M)+3*p*p/(2*M**3)
    clock = 3*T*fluctuations
    log_tail_p = math.log(6)-(1-d)**2/(6*s*s)
    tail_polynomial = 1/A+u+abs(b)*(d+15**.25*s)+hs*(d+945**.125*s)**2/2
    log_tail_norm = math.log(T*K*(M+3)*tail_polynomial)+log_tail_p/4
    assert log_tail_norm < math.log(1e-40)
    return soft, clock, soft+clock+1e-40, log_tail_norm


def polynomial_clock(ages):
    """Exact integer-spectrum Lagrange projectors instead of eigh."""
    spec = importlib.util.spec_from_file_location('review_1169_old_1112', BASE/'archive_1106_1119/1112/check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    material = module.Material(1)
    h, ident = material.hfk, np.eye(288, dtype=complex)
    projectors = []
    for j in range(4):
        p = ident.copy()
        for k in range(4):
            if j != k:
                p = p@(h-k*ident)/(j-k)
        projectors.append(p)
    algebra_error = max(float(np.linalg.norm(sum(projectors)-ident)),
                        *(float(np.linalg.norm(p@p-p)) for p in projectors))
    ready = material.encode_joint(np.eye(3, dtype=complex)[:, None, :], 0).reshape(288, 3)
    mask = np.zeros((4, 8, 3, 3))
    mask[:, [1, 3, 5, 7], :, :] = 1
    mask = mask.reshape(288)
    effects, grams = [], []
    for age in ages:
        evolved = sum(np.exp(-1j*j*age)*projectors[j]@ready for j in range(4))
        effect = evolved.conj().T@(mask[:, None]*evolved)
        effects.append(float(np.linalg.norm(effect-math.sin(age/2)**6*np.eye(3))))
        grams.append(float(np.linalg.norm(evolved.conj().T@evolved-np.eye(3))))
    assert max(effects+grams+[algebra_error]) < 2e-11
    return dict(projector_algebra_error=algebra_error, complete_qutrit_effect_errors=effects,
                complete_qutrit_isometry_errors=grams, method='integer-spectrum polynomial projectors; no eigendecomposition')


def local_records():
    """Local independent noise convolution and nonselective derivative cost."""
    nodes, weights = np.polynomial.hermite.hermgauss(64)
    weights = weights/math.sqrt(math.pi)
    sd = WIDTH/math.sqrt(2)
    max_error = 0.
    for offset in [-3., -1.1, 0., 1.3, 3.]:
        d = THRESHOLD+offset*WIDTH
        # Condition on A's independent Gaussian noise and integrate B's CDF.
        vals = [.5*math.erfc(-(d-THRESHOLD+math.sqrt(2)*sd*x)/(sd*math.sqrt(2))) for x in nodes]
        conditional = float(weights@vals)
        direct = .5*math.erfc(-offset/math.sqrt(2))
        max_error = max(max_error, abs(conditional-direct))
    derivative_norm = float(weights@((math.sqrt(2)*sd*nodes)**2/(4*sd**4)))
    energy = derivative_norm*(1/M+1/MU)/2
    expected = (1/M+1/MU)/(4*WIDTH**2)
    assert max_error < 1e-12 and abs(energy-expected) < 1e-15
    return dict(two_local_noise_convolution_max_error=max_error,
                local_gaussian_derivative_norm=derivative_norm,
                ideal_local_records_kinetic_energy_bound=energy,
                full_local_records_are_not_relative_lueders=True,
                physical_instrument_contract_is_assumed_not_constructed=True)


def text_audit():
    link_count = 0
    md = [BASE/n for n in EXPECTED if n.endswith('.md')]
    eq_line = None
    for p in md:
        data = p.read_bytes()
        assert not any(b < 32 and b not in (9, 10, 13) for b in data), p
        txt = data.decode('utf-8')
        assert txt.count('$$') % 2 == 0, p
        assert '\\[' not in txt and '\\]' not in txt, p
        for target in re.findall(r'\]\(([^)]+)\)', txt):
            if '://' in target or target.startswith('#'):
                continue
            target = target.split('#', 1)[0]
            assert (p.parent/target).exists(), (p, target)
            link_count += 1
        if p.name == 'proof.md':
            line = next(x for x in data.splitlines() if x.startswith(b'u(r)&='))
            eq_line = len(line)-len(line.rstrip(b'\\'))
    assert eq_line == 2
    return dict(local_link_occurrences=link_count, eq3_terminal_backslashes=eq_line,
                block_delimiters_paired=True, control_character_scan_passed=True)


def run():
    before = {n: sha(BASE/n) for n in EXPECTED}
    assert before == EXPECTED
    saved = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    sources = saved['sources_sha256']
    for n, h in sources.items():
        assert sha(BASE/n) == h
    max_mean = max_sd = max_momentum = max_phase = max_joint = max_bound = 0.
    outputs, ages = [], [math.pi/4]
    for radius, stored, predicted, cert in zip([1000, 2000], saved['group2_moving_two_body_reference'],
            saved['group4_conditional_actual_joint_prediction']['rows'], saved['group3_all_unknown_input_error_certificates']):
        for j, row in enumerate(stored['complete_mass_fibers']):
            means, covs, last = affine_flow(radius, T, M+j)
            shortening = means[0, 0]-means[0, 1]
            relative_var = covs[0, 0, 0]+covs[0, 1, 1]-2*covs[0, 0, 1]
            max_mean = max(max_mean, abs(shortening-row['mean_shortening']),
                           abs(means[0, 0]-row['mean_x']), abs(means[0, 1]-row['source_shift']))
            max_sd = max(max_sd, abs(math.sqrt(relative_var)-row['relative_sd']))
            max_momentum = max(max_momentum, abs(means[0, 2]-row['px']), abs(means[0, 3]-row['py']))
        mean, cov, _ = affine_flow(radius, T)
        phase = independent_phase(radius)
        ages.append(math.pi/4+phase)
        max_phase = max(max_phase, abs(phase-stored['phase']['tau']))
        relative_var = cov[0, 0, 0]+cov[0, 1, 1]-2*cov[0, 0, 1]
        z = (mean[0, 0]-mean[0, 1]-THRESHOLD)/math.sqrt(relative_var+WIDTH**2)
        joint = math.sin(ages[-1]/2)**6*.5*math.erfc(-z/math.sqrt(2))
        max_joint = max(max_joint, abs(joint-predicted['reference_joint_probability']))
        soft, clock, bound, log_tail = independent_bound(radius)
        max_bound = max(max_bound, abs(bound-cert['bounds']['full_prepared_isometry_error']))
        total = bound+.003+.0005+.000002
        assert total < .004
        outputs.append(dict(radius=radius, independently_integrated_clock_age=ages[-1],
                            joint_reference=joint, soft_norm_bound=soft, clock_norm_bound=clock,
                            total_error=total, log_tail_norm=log_tail))
    assert max_mean < 1e-12 and max_sd < 1e-12 and max_momentum < 3e-10
    assert max_phase < 3e-11 and max_joint < 3e-12 and max_bound < 1e-14
    lower = outputs[0]['joint_reference']-outputs[1]['joint_reference']-sum(x['total_error'] for x in outputs)
    assert lower > .02
    free_mean, free_cov, _ = affine_flow(1000, T, coupling=0.)
    free_var = free_cov[0, 0, 0]+free_cov[0, 1, 1]-2*free_cov[0, 0, 1]
    free_joint = math.sin((math.pi/4+T)/2)**6*.5*math.erfc(THRESHOLD/math.sqrt(2*(free_var+WIDTH**2)))
    assert abs(free_joint-saved['group4_conditional_actual_joint_prediction']['zero_coupling_control_joint_probability_each']) < 1e-14
    stats = math.ceil(math.log(400)/(2*.0005**2))
    assert stats == saved['group4_conditional_actual_joint_prediction']['independent_repeats_per_preparation']
    extra = polynomial_clock(ages)
    records = local_records()
    text = text_audit()
    assert before == {n: sha(BASE/n) for n in EXPECTED}
    for n, h in sources.items():
        assert sha(BASE/n) == h
    return dict(schema='round1169_independent_review_v1', round=1169, date='2026-10-13',
        verdict='PASS_WITH_SCOPE', blocker_count=0,
        independence='Reviewer provided pre-admission feedback, including relative-readout risks and residual decomposition; not blind review. Independent implementation owns no author core.',
        core_sha256=before, source_sha256=sources, reviewer_code_sha256=sha(Path(__file__)),
        author_default_run_observed=dict(command='python -B -X utf8 archive_1161_/1169/check.py', exit_code=0,
            message='1169 READ_ONLY PASS: 4 groups; saved results and all 7 sources unchanged'),
        independent_methods=['affine 5x5 Hamilton matrix exponential in long double',
            'Simpson integration on independent trajectories', '288D integer-spectrum polynomial projectors',
            'two local Gaussian records convolution and derivative norm'],
        independent_numerics=dict(max_mean_error=max_mean,max_relative_sd_error=max_sd,
            max_momentum_error=max_momentum,max_phase_error=max_phase,max_joint_error=max_joint,
            max_dynamic_bound_error=max_bound,rows=outputs,certified_contrast_lower=lower,
            zero_coupling_joint_probability_each=free_joint,independent_repeats_per_preparation=stats),
        complete_clock=extra, local_records=records, text_and_links=text,
        analytic_review=dict(direct_duhamel_on_common_prepared_reference=True,
            inverse_mass_residual_includes_momentum_fluctuations=True,
            reference_uniform_prepared_isometry_bound_not_all_position_inputs=True,
            probability_error_le_isometry_norm_for_purified_inputs=True,
            gaussian_tails_are_bounded_not_deleted=True,
            actual_poststates_not_identified_with_relative_lueders=True,
            contract_point003_is_additional_assumption=True,
            no_model_error_claim_to_full_sr_pn_or_gr=True,
            no_arbitrary_remote_intervention_causality_claim=True,
            no_instrument_hardware_or_new_axiom_claim=True),
        checks_complete=True, author_and_source_files_unchanged=True)


def same(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a: same(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): same(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=2e-9, abs_tol=3e-12), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    # Normalize NumPy scalar subclasses before opening the write-once receipt.
    # In particular Windows NumPy may preserve long-double scalar arithmetic.
    result = json.loads(json.dumps(run(), ensure_ascii=False, allow_nan=False,
                                  default=lambda x: float(x) if isinstance(x, np.floating) else int(x)))
    path = HERE/'review_checks.json'
    if args.write_once:
        with path.open('x', encoding='utf-8', newline='\n') as fp:
            json.dump(result, fp, ensure_ascii=False, indent=2, allow_nan=False)
            fp.write('\n')
    else:
        before = sha(path)
        same(result, json.loads(path.read_text(encoding='utf-8')))
        assert sha(path) == before
    print('1169 INDEPENDENT REVIEW PASS_WITH_SCOPE')
    print(json.dumps(result['independent_numerics'], ensure_ascii=False, indent=2))
