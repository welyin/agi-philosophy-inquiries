"""770: BRST-adapted subtraction blocks and same-action source identities.

Checks finite algebra only. Actual curved-space Hadamard coefficients and a
conservation-restoring local counterterm are NOT computed here.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_brst_relative_source as old
import joint_covariant_gauge_complex as cov
import joint_frame_hessian_matching as frame

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_local_brst_subtraction_results.json'


def mx(a):
    return float(np.max(np.abs(a)))


def subtraction_blocks():
    h, w, k, ks, p, d0, ell, gamma, ga, pairing, k2 = old.setup([1.3, .2, -.3, .4])
    # Arbitrary fiber-self-adjoint local matrix coefficients, not propagators.
    # No intertwining or field equation is imposed on these calibration data.
    x1 = np.linalg.solve(h, np.diag(np.linspace(.03, .17, 63)))
    x0 = np.linalg.solve(-w, np.diag(np.linspace(.021, .043, 16)))
    sizes = [63, 16, 16, 16]
    hx = old.block_matrix([[x1, k@x0, None, None], [x0@ks, None, None, None],
                           [None, None, None, x0], [None, None, x0, None]], sizes)
    d1 = p+k@ks
    r1, r0 = d1@x1, d0@x0
    q = ks@x1-x0@ks
    predicted = old.block_matrix([[r1-k@q, None, None, None], [q, r0, None, None],
                                   [None, None, r0, None], [None, None, None, r0]], sizes)
    errors = dict(
        BRST_exact=mx(gamma@hx-hx@ga),
        fiber_adjoint=mx(pairing@hx-hx.conj().T@pairing),
        equation_remainder_blocks=mx(ell@hx-predicted),
        free_gauge_fixing_insertion=mx(hx[:63, 63:79]-k@hx[79:95, 95:]),
    )
    assert max(errors.values()) < 2e-13
    assert mx(q) > .01 and mx(ell@hx) > .01
    return dict(errors=errors, equation_remainder_norm=mx(ell@hx),
                intertwining_remainder_norm=mx(q),
                scope='Full 111-field matrix identities for arbitrary local coefficients. Exact linear BRST subtraction does not force its equation remainder to vanish; no actual parametrix coefficient is simulated.')


def scalar_source_and_split():
    # Original common scalar vacuum, potential V/F^2, and original SU(2) action.
    phi, f0, metric, hessian, _ = frame.vacuum()
    generators = cov.scalar_generators()[8:11]
    r = np.stack([t@phi for t in generators], axis=1)
    # A declared transverse linear gauge with a small physical-direction mixing.
    # It is a finite algebra calibration, not the complete 753 gauge fixing.
    mix = np.array([[.07, -.12, .09, .03, .11],
                    [-.04, .08, .13, -.06, .05],
                    [.06, .09, -.05, .12, -.08]])
    g = r.T@metric+.25*mix
    m = g@r
    a = hessian+g.T@g
    c = np.linalg.inv(a)
    mi = np.linalg.inv(m)
    assert np.linalg.eigvalsh(a)[0] > 0 and abs(np.linalg.det(m)) > 1e-5
    chol = np.linalg.cholesky(c)
    jb = sum(old.potential_gradient_jet(phi, chol[:, j])[2] for j in range(5))
    jghost = np.array([-np.trace(mi@g@np.stack([t[:, i] for t in generators], axis=1))
                       for i in range(5)])
    wb, wg = r.T@jb, r.T@jghost
    # Change only the auxiliary gauge functional. The two determinant terms
    # cancel separately from the genuine background potential variation.
    dg = np.array([[.13, -.04, .02, -.03, .08],
                   [.07, .11, -.09, .04, -.06],
                   [-.05, .01, .06, .10, .02]])
    split_b = float(np.trace(dg@c@g.T))
    split_ghost = float(-np.trace(mi@dg@r))
    # Independent differences of the finite one-loop expression.
    def determinant_part(offset):
        gt = g+offset*dg
        aa, mm = hessian+gt.T@gt, gt@r
        sa, la = np.linalg.slogdet(aa); sm, lm = np.linalg.slogdet(mm)
        assert sa > 0 and sm != 0
        return .5*la-lm
    step = 1e-4
    derivative = (determinant_part(step)-determinant_part(-step))/(2*step)
    errors = dict(
        vacuum_equation=mx(old.potential_gradient_jet(phi, np.zeros(5))[0]),
        on_shell_gauge_identity=mx(hessian@r),
        inverse_gauge_identity=mx(c@g.T-r@mi),
        boson_plus_ghost_Ward=mx(wb+wg),
        background_gauge_fixing_cancellation=abs(split_b+split_ghost),
        independent_logdet_derivative=abs(derivative),
    )
    assert max(v for k, v in errors.items() if k != 'independent_logdet_derivative') < 3e-13
    assert abs(derivative) < 5e-9
    assert mx(wb) > .01 and abs(split_b) > .01
    return dict(errors=errors, bosonic_Ward_terms=wb.tolist(),
                ghost_Ward_terms=wg.tolist(),
                gauge_fixing_boson_term=split_b, gauge_fixing_ghost_term=split_ghost,
                ghost_matrix_determinant=float(np.linalg.det(m)),
                minimum_bosonic_eigenvalue=float(np.linalg.eigvalsh(a)[0]),
                scope='Original five-scalar potential and SU(2) generators in a zero-dimensional Gaussian calibration, with a declared auxiliary gauge matrix. This is neither the physical Hadamard covariance nor a continuum loop computation.')


def run():
    result = dict(round=770, tests_run=2, failures=0, errors=0,
                  exact_subtraction_algebra=subtraction_blocks(),
                  same_potential_source_and_gauge_fixing=scalar_source_and_split())
    deps = ('research_note_734.md', 'research_note_735.md', 'research_note_765.md',
            'research_note_768.md', 'research_note_769.md',
            'joint_brst_relative_source.py', 'joint_covariant_gauge_complex.py',
            'joint_frame_hessian_matching.py', 'round770_drafts/research_note_770_working.md')
    result['dependency_hashes'] = {p: hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps}
    result['scope'] = ('The same linear BRST extension admits a locally constructed subtraction with exact free BRST identities. '
                       'The finite same-state all-sector tadpole candidate and background quadratic insertion agree at this level. '
                       'The remaining on-shell Ward defect is a local equation-remainder trace, not proven to vanish or to have a permitted local primitive. '
                       'Interacting Ward identities, absolute conserved backreaction and actual instruments remain open.')
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
