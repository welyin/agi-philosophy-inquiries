"""1238: the original relational read vertex is not metric-velocity free.

Two diagnostics: exact rational old-jet calibration; full original 64-component
Dirac density, including spin, internal connection, mass and Nambu half weight.
Default is read-only replay.  --write-once creates results.json exclusively.
No original file, process, preparation, or calibration is changed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
TARGET = HERE / 'results.json'
JET = BASE / 'archive_894_922/907/material_reference_jets_results.json'
INITIAL = BASE / 'archive_894_922/900/current_reference_initial_budget_results.json'
SOURCES = [
    'archive_554_584/research_note_573.md',
    'archive_797_818/research_note_802.md',
    'archive_797_818/802/original_bff_vertex.py',
    'archive_819_853/research_note_847.md',
    'archive_819_853/research_note_852.md',
    'archive_854_872/research_note_860.md',
    'archive_854_872/research_note_861.md',
    'archive_894_922/research_note_920.md',
    str(JET.relative_to(BASE)).replace('\\', '/'),
    str(INITIAL.relative_to(BASE)).replace('\\', '/'),
]


def need(value, label):
    if not bool(value):
        raise AssertionError(label)


def inverse(a):
    n = len(a)
    m = [list(row) + [F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for k in range(n):
        p = next(i for i in range(k, n) if m[i][k])
        m[k], m[p] = m[p], m[k]
        z = m[k][k]
        m[k] = [x / z for x in m[k]]
        for i in range(n):
            if i != k:
                z = m[i][k]
                m[i] = [x - z * y for x, y in zip(m[i], m[k])]
    return [row[n:] for row in m]


def fixture():
    jet = json.loads(JET.read_text('utf-8'), parse_float=F)['rows'][0]
    ini = next(r for r in json.loads(INITIAL.read_text('utf-8'), parse_float=F)['rows']
               if r['N'] == 24)
    j = [[F(x) for x in row] for row in jet['data']['old']['Jacobian']]
    psi = ini['psi']
    e = inverse(j)
    # b_z is the original 847 orthonormal z metric test. kappa*v_i=1
    # only removes a fixed positive normalization from this coefficient check.
    c = 2 * psi**4
    w = e[3]
    velocity = [F(0), F(0), -psi**-8 * j[0][2]**2,
                -psi**-8 * j[1][2]**2]
    a = sum(w[i] * velocity[i] for i in range(4))
    k2 = c * a * a
    return jet, psi, j, e, c, w, velocity, a, k2


def exact_old_calibration():
    jet, psi, j, e, c, w, velocity, a, k2 = fixture()
    row = [sum(w[i] * j[i][k] for i in range(4)) for k in range(4)]
    need(row == [0, 0, 0, 1], 'old fixed tensor calibration')
    need(k2 > 0, 'nonzero metric-velocity quadratic coefficient')
    b = [[c * w[i] * w[k] for k in range(4)] for i in range(4)]
    rows = []
    for t in (F(-1, 50), F(-1, 100), F(0), F(1, 100), F(1, 50)):
        q = [[j[i][k] + (t * velocity[i] if k == 0 else 0)
              for k in range(4)] for i in range(4)]
        tensor = [[sum(q[i][mu] * b[i][k] * q[k][nu]
                       for i in range(4) for k in range(4))
                   for nu in range(4)] for mu in range(4)]
        z = [t * a, F(0), F(0), F(1)]
        expected = [[c * z[mu] * z[nu] for nu in range(4)] for mu in range(4)]
        need(tensor == expected, 'all components of the old tensor')
        need(tensor[0][0] == k2 * t * t, 'quadratic coefficient')
        rows.append({'velocity_label': str(t), 'C00': float(tensor[0][0]),
                     'C0z': float(tensor[0][3]), 'exact_identity': True})
    return dict(
        status='PASS', old_grid_N=24, old_time=0,
        point=[float(x) for x in jet['point']],
        jacobian_rows='h,s,(dh)^2,(ds)^2', jacobian_columns='t,x,y,z',
        psi=float(psi), inverse_z_X4=float(e[3][3]), B44=float(b[3][3]),
        d_dtX_d_dotgamma_yy=[float(x) for x in velocity],
        coefficient_C00_velocity_squared=float(k2),
        coefficient_is_exact_nonzero_rational=True,
        rational_check_rows=rows, fixed_B_not_recalibrated=True,
        all_saved_spatial_dh_roundoff_retained=True,
        continuum_background_error_certified=False,
    )


def old_matrices():
    sys.path.insert(0, str(ROOT / 'scripts'))
    sys.path.insert(0, str(BASE / 'archive_797_818/802'))
    from run_research_gr_response_active import GRResponseLayout, GRResponseRuntime
    with GRResponseRuntime(GRResponseLayout()).installed():
        import original_bff_vertex
    return original_bff_vertex


def full_dirac_contraction():
    old = old_matrices()
    jet, ps, j, e, c, w, dv, a, k2 = fixture()
    psi, c, a, k2 = map(float, (ps, c, a, k2))
    e0 = np.diag([1., psi**2, psi**2, psi**2])
    g0 = e0.T @ old.ETA @ e0
    gi0 = np.linalg.inv(g0)
    density = float(np.linalg.det(e0))
    phi = np.array([0., float(jet['data']['old']['value'][0]), 0., 0.,
                    float(jet['data']['old']['value'][1])])
    mass = old.mass(phi)
    rng = np.random.default_rng(1238002)
    internal = old.gauge_connection(rng)
    dg0 = rng.normal(size=(4, 4, 4)) * .025
    dg0 += dg0.transpose(0, 2, 1)
    metric_velocity = np.zeros_like(dg0)
    metric_velocity[0, 2, 2] = 1.
    # A nonzero direction derivative is retained. Its spin contraction vanishes
    # by 802, for every symmetric direction jet, not by freezing the connection.
    db = rng.normal(size=(4, 4, 4)) * .03
    db += db.transpose(0, 2, 1)
    _, vectors = np.linalg.eigh(old.ALPHA[1])
    u = vectors[:, -1]
    need(abs((u.conj() @ old.ALPHA[1] @ u).real - 1) < 2e-13,
         'original Clifford eigenvector')
    dpsi0 = (rng.normal(size=(4, 64)) + 1j*rng.normal(size=(4, 64))) * .02
    coupling, q = F(2, 7), F(3, 11)
    factor = float(coupling*q)

    def data(b, dg, t):
        g = g0 + t*b
        metric_jet = dg + t*db
        root = old.root_near_one(gi0@g)
        frame = e0@root
        df = np.array([e0@old.sylvester(root, gi0@metric_jet[mu]) for mu in range(4)])
        inv = np.linalg.inv(frame)
        christoffel = old.christoffel(g, metric_jet)
        omega = np.array([frame@christoffel[mu]@inv-df[mu]@inv for mu in range(4)])
        C = np.array([sum((inv[mu, aa]*old.ALPHA[aa] for aa in range(4)),
                         np.zeros((64, 64), complex)) for mu in range(4)])
        spin = np.array([old.spin_lift(x) for x in omega])
        return float(abs(np.linalg.det(frame))), C, spin+internal

    def action_value(b, dg, k, t, keep_density=True):
        vol, C, connection = data(b, dg, t)
        derivatives = dpsi0.copy()
        derivatives[1] += 1j*k*u
        val = .5*(vol if keep_density else density)*old.bilinear(
            C, connection, mass, u, derivatives, u, derivatives)
        need(abs(val.imag) < 2e-11, 'Hermitian full density')
        return float(val.real)

    taus = (-.02, 0., .02)
    def response(tau, k, step, keep_density=True):
        z = np.array([tau*a, 0., 0., 1.])
        b = c*np.outer(z, z)
        dg = dg0 + tau*metric_velocity
        # The spin-variation contraction is zero in the symmetric section,
        # so -lambda*q*delta_b S_F is this Hilbert read contraction locally.
        return -factor*(action_value(b, dg, k, step, keep_density)
                        -action_value(b, dg, k, -step, keep_density))/(2*step)

    expected = -factor*density*k2/(4*psi**2)
    rows = []
    for step in (4e-5, 2e-5):
        slopes = [response(tau, 2., step)-response(tau, 1., step) for tau in taus]
        measured = (slopes[2]+slopes[0]-2*slopes[1])/(2*taus[2]**2)
        need(abs(measured-expected) < 3e-7, 'full Dirac mixed coefficient')
        rows.append(dict(metric_variation_step=step, full_response_k_slopes=slopes,
                         velocity_squared_k_coefficient=measured,
                         absolute_error=abs(measured-expected)))
    wrong = [response(tau, 2., 2e-5, False)-response(tau, 1., 2e-5, False)
             for tau in taus]
    wrong_coeff = (wrong[2]+wrong[0]-2*wrong[1])/(2*taus[2]**2)
    need(abs(wrong_coeff) < 3e-7 and abs(expected) > 1.,
         'deleting density loses this actual coefficient')
    nu, dc, ds = old.frame_tangent(e0, c*np.diag([0., 0., 0., 1.]), db)
    _, C0, _ = data(np.zeros((4, 4)), dg0, 0.)
    spin_contraction = sum((C0[mu]@ds[mu]-ds[mu].conj().T@C0[mu]
                            for mu in range(4)), np.zeros((64, 64), complex))
    spin_error = float(np.max(np.abs(spin_contraction)))
    need(spin_error < 1e-12, 'retained differentiated spin connection contraction')
    need(np.linalg.norm(mass) > 0 and np.linalg.norm(internal) > 0,
         'original full mass and nonzero internal connection retained')
    return dict(status='PASS', original_selfdual_dimension=64,
                independent_CAR_modes=32, nambu_half_weight=.5,
                all_original_mass_matrices_retained=True,
                full_density_tetrad_LC_internal_connection_retained=True,
                nonzero_background_connection_is_local_diagnostic=True,
                actual_background_LC_jet_claimed=False,
                k_is_off_shell_derivative_label_not_preparation=True,
                lambda_label=str(coupling), q_label=str(q),
                expected_velocity_squared_k_coefficient=expected,
                finite_difference_rows=rows,
                spin_connection_variation_contraction_error=spin_error,
                density_deleted_coefficient=wrong_coeff,
                nonzero_record_statistic_claimed=False)


def run():
    before = {p: hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in SOURCES}
    first, second = exact_old_calibration(), full_dirac_contraction()
    after = {p: hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in SOURCES}
    need(before == after, 'historical inputs untouched')
    return dict(round=1238, status='PASS', new_scientific_check_groups=2,
                exact_saved_material_jet=first, full_original_Dirac_vertex=second,
                historical_source_sha256=before,
                analytic_scope='Fixed old off-shell Lagrangian representative; direct CT Legendre inheritance fails.',
                canonical_or_integration_by_parts_no_go=False,
                finite_EFT_reorganization_excluded=False,
                CT_adopted=False, Einstein_equations_derived=False,
                actual_device_or_total_quantum_state_certified=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_once:
        with TARGET.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('WROTE_ONCE: 1238 results.json; two scientific groups PASS')
    else:
        need(result == json.loads(TARGET.read_text('utf-8')), 'READ_ONLY_MATCH')
        print('READ_ONLY_MATCH: 1238; two scientific groups PASS')
