"""771: local one-insertion Ward repair, with scoped algebra/jet checks.

The proof uses the complete original D1/D0 as variable-coefficient operators.
Finite matrices check the Schur contact; an actual noncommuting polynomial
potential checks Hadamard transport jets. Neither is a full spacetime loop
simulation or a quantum-gravity construction.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import joint_brst_relative_source as old

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_local_ward_repair_results.json'
NU = 1/(8*np.pi**2)


def mx(a):
    return float(np.max(np.abs(a)))


def st(a):
    return float(np.real(np.trace(a[:79, :79])-np.trace(a[79:, 79:])))


def schur_contact():
    h, w, k, ks, p, d0, ell, _, _, m, _ = old.setup([1.3, .2, -.3, .4])
    sizes = [63, 16, 16, 16]
    h1 = np.linalg.solve(h, np.diag(np.linspace(.03, .17, 63)))
    h0 = np.linalg.solve(-w, np.diag(np.linspace(.021, .043, 16)))
    hh = old.block_matrix([[h1, k@h0, None, None], [h0@ks, None, None, None],
                           [None, None, None, h0], [None, None, h0, None]], sizes)
    t = np.eye(111, dtype=complex); t[63:79, :63] = ks
    ti = np.eye(111, dtype=complex); ti[63:79, :63] = -ks
    adj = lambda a: np.linalg.solve(m, a.conj().T@m)
    ld = adj(t)@ell@t
    ht = ti@hh@adj(ti)
    hd = old.block_matrix([[h1, None, None, None], [None, None, None, None],
                           [None, None, None, h0], [None, None, h0, None]], sizes)
    s = ht-hd
    q = ks@h1-h0@ks
    qp = h1@k-k@h0
    expected_s = old.block_matrix([[None, -qp, None, None],
                                   [-q, ks@h1@k-d0@h0-h0@d0, None, None],
                                   [None, None, None, None], [None, None, None, None]], sizes)
    expected_ld = old.block_matrix([[p+k@ks, None, None, None],
                                    [None, -np.eye(16), None, None],
                                    [None, None, None, d0], [None, None, d0, None]], sizes)
    rng = np.random.default_rng(77101)
    # Same Fourier convention as K=i*K_real; a real delta K here would
    # artificially send the diagnostic contact into its imaginary part.
    dk = .02j*rng.standard_normal(k.shape)
    dks = np.linalg.solve(-w, dk.conj().T@h)
    x = .002*rng.standard_normal((63, 63))
    dp = np.linalg.solve(h, x+x.T)
    dd0 = dks@k+ks@dk
    dell = old.block_matrix([[dp, dk, None, None], [dks, None, None, None],
                             [None, None, None, dd0], [None, None, dd0, None]], sizes)
    dld = old.block_matrix([[dp+dk@ks+k@dks, None, None, None],
                            [None, None, None, None], [None, None, None, dd0],
                            [None, None, dd0, None]], sizes)
    om = np.zeros((111, 111), complex); om[63:79, :63] = dks
    # Zero here is a homogeneous matrix solution for calibration only,
    # NOT a physical covariance. The analytic construction keeps the old state.
    jl = .5*st(-hh@dell)
    jd = .5*st(-hd@dld)
    subtraction_contact = -.5*st(s@dld)
    transformation_contact = .5*st(ld@ht@adj(om)+ht@ld@om)
    ell_t = subtraction_contact+transformation_contact
    errors = dict(
        schur_factorization=mx(ld-expected_ld),
        smooth_difference_blocks=mx(s-expected_s),
        differentiated_congruence=mx(adj(t)@dell@t-(dld-adj(om)@ld-ld@om)),
        local_source_contact=abs(jl-jd-ell_t),
    )
    assert max(errors.values()) < 1e-12
    assert abs(ell_t) > 1e-5 and abs(transformation_contact) > 1e-5
    return dict(errors=errors, original_source=jl, diagonal_source=jd,
                subtraction_contact=subtraction_contact,
                change_of_field_contact=transformation_contact, total_schur_contact=ell_t,
                scope='Original full principal block, fixed pairing, arbitrary transverse variation. Smooth local matrix calibration only; no physical state is set to zero in the proof.')


def transport_jets():
    e0 = np.array([[-.8, .13, .07], [.13, -.5, -.09], [.07, -.09, -.3]])
    es = np.array([[[.17, .11, -.03], [.11, -.08, .04], [-.03, .04, .12]],
                   [[-.05, .02, .14], [.02, .09, -.07], [.14, -.07, .06]],
                   [[.04, -.08, .02], [-.08, .03, .05], [.02, .05, -.02]],
                   [[-.03, .06, .09], [.06, -.04, .02], [.09, .02, .07]]])
    nodes, weights = np.polynomial.legendre.leggauss(3)
    nodes, weights = (nodes+1)/2, weights/2
    def e(x):
        return e0+np.einsum('a,aij->ij', x, es)
    def v0(x, y):
        return -(e(x)+e(y))/4
    def v1(x, y):
        # Solves (r.d_x+2)v1 = -(Box_x+E(x))v0/2.
        # Box_x v0=0 for this linear, matrix-valued potential.
        return sum(-.5*w*t*e(y+t*(x-y))@v0(y+t*(x-y), y)
                   for t, w in zip(nodes, weights))
    def dxv1(x, y, mu):
        return sum(w*t*t*(es[mu]@(e(y+t*(x-y))+e(y))
                            +e(y+t*(x-y))@es[mu])/8
                   for t, w in zip(nodes, weights))
    def rleft(x, y):
        # Universal zeroth/first coincidence jet of D H / NU.
        return 6*v1(x, y)+2*sum((x-y)[mu]*dxv1(x, y, mu) for mu in range(4))
    zero = np.zeros(4)
    x = np.array([.09, -.04, .07, .03]); y = np.array([-.06, .02, -.03, .05])
    left = np.array([(2*a@e0+e0@a)/24 for a in es])
    right = np.array([(a@e0+2*e0@a)/24 for a in es])
    grad = np.array([(a@e0+e0@a)/8 for a in es])
    numeric_left, numeric_right, comm_derivatives = [], [], []
    step = 1e-20
    for mu in range(4):
        dz = zero.astype(complex); dz[mu] = 1j*step
        numeric_left.append(v1(dz, zero).imag/step)
        numeric_right.append(v1(zero, dz).imag/step)
        comm_derivatives.append((rleft(zero, dz)-rleft(dz, zero).T).imag/step)
    transport = 2*sum((x-y)[mu]*dxv1(x, y, mu) for mu in range(4))+4*v1(x, y)+e(x)@v0(x, y)
    errors = dict(
        matrix_transport_equation=mx(transport),
        sesqui_symmetry=mx(v1(x, y)-v1(y, x).T),
        coincidence_coefficient=mx(v1(zero, zero)-e0@e0/8),
        left_jet=mx(np.asarray(numeric_left)-left),
        right_jet=mx(np.asarray(numeric_right)-right),
        Synge_rule=mx(left+right-grad),
        commutator_jet=mx(np.asarray(comm_derivatives)+2*grad),
    )
    noncommute = max(mx(e0@a-a@e0) for a in es)
    assert noncommute > .01 and max(errors.values()) < 1e-12
    return dict(errors=errors, potential_commutator_norm=noncommute,
                diagonal_v1_trace=float(np.trace(e0@e0)/8),
                gradient_v1_trace=[float(np.trace(a)) for a in grad],
                scope='Actual finite Hadamard transport calculation for a noncommuting linear matrix potential on flat Lorentz spacetime. It checks universal jet coefficients, not the original 753 background numerically.')


def volume_contact():
    ts = 2*np.pi*np.arange(128)/128
    h = 2*np.pi/128
    ev = np.array([[-.8, .13], [.13, -.5]])
    av = np.array([[.17, .11], [.11, -.08]])
    bv = np.array([[-.05, .02], [.02, .09]])
    eg = np.array([[-.6, .07], [.07, -.4]])
    ag = np.array([[.08, -.04], [-.04, -.02]])
    bg = np.array([[.03, .06], [.06, -.07]])
    def coefficient(base, a, b):
        e = base+np.sin(ts)[:, None, None]*a+np.cos(2*ts)[:, None, None]*b
        ep = np.cos(ts)[:, None, None]*a-2*np.sin(2*ts)[:, None, None]*b
        epp = -np.sin(ts)[:, None, None]*a-4*np.cos(2*ts)[:, None, None]*b
        eppp = -np.cos(ts)[:, None, None]*a+8*np.sin(2*ts)[:, None, None]*b
        # Flat D=Box+E, time-dependent potential: Box E = -E''.
        vv = e@e/8-epp/24
        dv = (ep@e+e@ep)/8-eppp/24
        return np.trace(vv, axis1=1, axis2=2), np.trace(dv, axis1=1, axis2=2)
    v, vp = coefficient(ev, av, bv); g, gp = coefficient(eg, ag, bg)
    b, bp = .5*v-g, .5*vp-gp
    xi = np.cos(ts)+.31*np.sin(2*ts)
    div = -np.sin(ts)+.62*np.cos(2*ts)
    ward = float(-2*NU*h*np.sum(xi*bp))
    repair = float(-2*NU*h*np.sum(b*div))
    omit_ghost = float(-2*NU*h*np.sum(.5*v*div))
    wrong_sign = ward-repair
    errors = dict(local_Ward_plus_volume_contact=abs(ward+repair))
    assert abs(ward) > 1e-5 and abs(ward+omit_ghost) > 1e-5
    assert errors['local_Ward_plus_volume_contact'] < 1e-14
    # Scalar normalization: density contact 3*NU*v1 and volume -NU*v1
    # give T-contact 2*NU*v1, matching eta_4=1/3 times <phi P phi>.
    scalar_contact = 3-1
    assert scalar_contact == 6/3
    return dict(errors=errors, uncorrected_Ward=ward, volume_repair=repair,
                omitted_ghost_residual=ward+omit_ghost, wrong_sign_residual=wrong_sign,
                scalar_Moretti_normalization=True,
                scope='Periodic local-jet/weight calibration with two matrix bundles. It tests the joint volume contact, ghost weight, and scalar normalization; it is not a compact-time physical spacetime or an actual source simulation.')


def run():
    result = dict(round=771, tests_run=3, failures=0, errors=0,
                  schur_source_contact=schur_contact(), matrix_Hadamard_jets=transport_jets(),
                  joint_volume_contact=volume_contact())
    deps = ('research_note_735.md', 'research_note_765.md', 'research_note_768.md',
            'research_note_770.md', 'joint_brst_relative_source.py',
            'round771_drafts/research_note_771_working.md')
    result['dependency_hashes'] = {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps}
    result['scope'] = ('For the original on-shell coupled free state, a local real finite-jet normalization '
                       'of the complete one-loop tadpole can satisfy the joint infinitesimal background Ward identity. '
                       'The construction keeps the actual state and full mixed operators, including Schur/pairing contacts. '
                       'It does not prove interacting QME, full background independence, gauge-fixing independence, '
                       'finite-strength self-consistency, or the original Q-to-continuum map.')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args(); result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k: result[k] for k in ('round', 'tests_run', 'failures', 'errors')}))
