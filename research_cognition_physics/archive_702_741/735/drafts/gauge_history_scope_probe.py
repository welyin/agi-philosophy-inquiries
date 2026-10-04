"""735 working probe: full original local symbols under time-dependent gauge changes.

This is a scope diagnostic, NOT a continuum anomaly or Ward-normalization proof.
Frozen research counts are not advanced by this file.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_dynamic_continuum_reference as ref
import joint_relative_source_development as prior
sys.path.insert(0, str(ARCHIVE / 'round732_drafts'))
import source_feedback_entry as record

TARGET = HERE / 'gauge_history_scope_probe_results.json'


def maxabs(x):
    return float(np.max(np.abs(x)))


def bump(t, lo, hi):
    s = (t - lo) / (hi - lo)
    if not 0 < s < 1:
        return 0., 0.
    value = np.exp(4 - 1 / (s * (1 - s)))
    derivative = value * (1 - 2 * s) / (s*s*(1-s)**2*(hi-lo))
    return float(value), float(derivative)


def frame(t):
    x, dx = bump(t, -.16, -.02)
    y, dy = bump(t, -.12, -.005)
    x *= .8; dx *= .8; y *= .55; dy *= .55
    w1 = np.cos(x/2)*np.eye(2) + 1j*np.sin(x/2)*ref.SIG[0]
    w2 = np.cos(y/2)*np.eye(2) + 1j*np.sin(y/2)*ref.SIG[1]
    w = w1 @ w2
    dw = .5j*dx*ref.SIG[0]@w + w1@(.5j*dy*ref.SIG[1])@w2
    temporal = 1j * dw @ w.conj().T
    coefficients = np.array([np.trace(s @ temporal).real for s in ref.SIG])
    r = ref.old.matter.representation(np.eye(3), w, 1.+0j)
    v64 = ref.old.block([r, r.conj()])
    v = prior.assemble([v64, v64])
    physical_temporal = ref.gauge_h(np.array([coefficients, np.zeros(3), np.zeros(3)]), np.zeros(3))[0]
    temporal64 = ref.old.bdg(physical_temporal, np.zeros((32, 32)))
    return w, v, prior.assemble([temporal64, temporal64])


def transformed_background(t):
    g, phi, a, a0, _ = ref.collar(t, .08)
    w, v, temporal = frame(t)
    x = w @ (phi[:2] + 1j*phi[2:4])
    rotated_phi = np.r_[x.real, x.imag, phi[4]]
    rotated_a = np.zeros_like(a)
    for i in range(3):
        component = w @ sum(a[i, j]*ref.SIG[j]/2 for j in range(3)) @ w.conj().T
        rotated_a[i] = [np.trace(s @ component).real for s in ref.SIG]
    return (g, rotated_phi, rotated_a, a0), v, temporal


def operator(bg):
    return prior.matrix_and_source(*bg, *[np.zeros_like(z) for z in bg])[0]


def local_dictionary_check():
    error = 0.; derivative_error = 0.; nonzero_temporal = 0.
    for t in (-.135, -.1, -.075, -.03):
        bg, v, temporal = transformed_background(t)
        raw = ref.collar(t, .08)[:4]
        b = operator(raw)
        independent = operator(bg) + temporal
        error = max(error, maxabs(independent - v@b@v.conj().T - temporal))
        h = 1.e-6
        dv = (frame(t+h)[1] - frame(t-h)[1]) / (2*h)
        derivative_error = max(derivative_error, maxabs(temporal - 1j*dv@v.conj().T))
        nonzero_temporal = max(nonzero_temporal, float(np.linalg.norm(temporal)))
    assert error < 3.e-13 and derivative_error < 2.e-6 and nonzero_temporal > 1.
    return dict(full_nambu_dimension=128, independent_mass_and_spatial_connection_error=error,
                temporal_connection_vs_frame_derivative_error=derivative_error,
                temporal_connection_norm_witness=nonzero_temporal,
                temporal_connection_is_distinct_from_original_spatial_U1_a0=True)


def unitary(b, dt):
    e, v = np.linalg.eigh(b)
    return (v*np.exp(-1j*dt*e)) @ v.conj().T


@lru_cache(maxsize=2)
def transport(steps):
    start, end = -.24, 0.
    dt = (end-start)/steps
    p = ref.projector(operator(ref.collar(start, .08)[:4]))
    exact = p.copy(); actual = p.copy(); missing = p.copy()
    for j in range(steps):
        left = start + j*dt; t = left + dt/2
        b = operator(ref.collar(t, .08)[:4])
        rotated, _, temporal = transformed_background(t)
        rb = operator(rotated)
        u = unitary(b, dt)
        uv = frame(left+dt)[1] @ u @ frame(left)[1].conj().T
        up = unitary(rb+temporal, dt)
        um = unitary(rb, dt)
        p = u@p@u.conj().T
        exact = uv@exact@uv.conj().T
        actual = up@actual@up.conj().T
        missing = um@missing@um.conj().T
    return p, exact, actual, missing


def same_past_record_check():
    rows = []
    b = operator(ref.collar(0., .08)[:4])
    phi = ref.collar(0., .08)[1]
    mass64 = ref.old.bdg(*ref.old.matter.mass_matrices(phi))
    geometry_source = -(b - prior.assemble([mass64, mass64]))
    for n in (192, 384):
        base, transported, direct, wrong = transport(n)
        posts = [p + record.record_change(p) for p in (base, transported, direct, wrong)]
        values = [prior.source(p, geometry_source) for p in posts]
        err = float(np.linalg.norm(direct-base, 'fro'))
        mismatch = float(np.linalg.norm(wrong-base, 'fro'))
        exact_error = max(maxabs(transported-base), abs(values[1]-values[0]))
        assert exact_error < 4.e-10 and mismatch > .001
        rows.append(dict(steps=n, exact_transport_covariance_and_record_source_error=exact_error,
                         independent_midpoint_covariance_error=err,
                         omitted_temporal_connection_covariance_difference=mismatch,
                         original_post_record_geometric_source=values[0],
                         independent_midpoint_source_difference=values[2]-values[0],
                         omitted_temporal_connection_source_difference=values[3]-values[0]))
    ratio = rows[0]['independent_midpoint_covariance_error']/rows[1]['independent_midpoint_covariance_error']
    assert 3.5 < ratio < 4.5 and rows[1]['independent_midpoint_covariance_error'] < 1.e-4
    return dict(rows=rows, integration_error_ratio=ratio,
                finite_symbol_check_not_continuum_gauge_anomaly_calculation=True)


def phase_blindness_check():
    # A deliberately gauge-breaking c-number Hamiltonian addition.
    # It is not an allowed final renormalization or a physical new coupling.
    # Its four-dimensional local analogue is kappa*sqrt(-g)*phi[0]*phi[1].
    t = -.07
    bg = ref.collar(t, .08)[:4]
    phi = bg[1]
    x = phi[:2]+1j*phi[2:4]
    tangent = .5j*ref.SIG[1]@x
    delta = np.r_[tangent.real, tangent.imag, 0.]
    exact_variation = float(delta[0]*phi[1]+phi[0]*delta[1])
    h = 1.e-5
    def r(alpha):
        w = np.cos(alpha/2)*np.eye(2)+1j*np.sin(alpha/2)*ref.SIG[1]
        xx = w@x
        return float(xx[0].real*xx[1].real)
    difference = (r(h)-r(-h))/(2*h)
    assert abs(difference-exact_variation) < 1.e-10 and abs(exact_variation) > .01
    n = 384; dt = .24/n
    phase_integral = 0.
    for j in range(n):
        tt = -.24+(j+.5)*dt
        original_phi = ref.collar(tt, .08)[1]
        transformed_phi = transformed_background(tt)[0][1]
        phase_integral += dt*(transformed_phi[0]*transformed_phi[1] - original_phi[0]*original_phi[1])
    p = transport(384)[0]
    # Overall many-body scalar phases cancel in density matrices; no scalar
    # identity is added to a BdG matrix (that would violate the CAR dictionary).
    phase = np.exp(-1j*phase_integral)
    equality = maxabs(phase*p*phase.conjugate()-p)
    assert equality < 2.e-15 and abs(phase-1) > 1.e-4
    return dict(purpose='Counterexample to certifying absolute source normalization by state/record covariance alone',
                gauge_variation_of_forbidden_local_term=exact_variation,
                independent_finite_difference=difference,
                pure_gauge_history_scalar_phase_integral=float(phase_integral),
                nontrivial_scalar_phase_distance_from_one=float(abs(phase-1)),
                unchanged_covariance_error=equality,
                c_number_is_on_many_body_Fock_space_not_a_BdG_identity_shift=True,
                this_does_not_establish_a_nonzero_anomaly_in_the_original_model=True)


def run():
    checks = ('local_dictionary_check', 'same_past_record_check', 'phase_blindness_check')
    result = {name: globals()[name]() for name in checks}
    dependencies = ('joint_dynamic_continuum_reference.py', 'joint_relative_source_development.py',
                    'round732_drafts/source_feedback_entry.py', 'research_note_623.md',
                    'research_note_628.md', 'research_note_734.md', 'round735_drafts/STATUS.md',
                    'round735_drafts/counterterm_scope_entry.md', 'round735_drafts/entry_checks.json')
    return dict(round_in_progress=735, latest_completed_round=734, completed_new_round=False,
                probe_checks=len(checks), failures=0, results=result,
                dependencies={name: hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in dependencies},
                scope='Original complete finite local-symbol gauge-history and record calibration; phase-blindness limitation. No continuum Wick normalization or self-consistent gravity existence is claimed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as handle:
            handle.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result, ensure_ascii=False))
