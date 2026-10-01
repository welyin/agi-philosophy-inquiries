"""591: original finite-graph vacuum response, records, and ADM comparison.

The full spectrum and ground wave function are NOT computed. Spectral
existence/convergence are analytic results in research_note_591.md. Checks
use original geometry coefficients and actual instruments/packet, plus
exact scalar rational identities; there is no surrogate bath Hamiltonian.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_full_spatial_metric as model
import joint_geometry_work_noise as previous
import joint_record_source_compression as records

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_ground_geometry_response_results.json'
spec = importlib.util.spec_from_file_location('entry591', HERE / 'round591_drafts/gapped_source_entry.py')
entry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(entry)
HBAR, W = previous.HBAR, previous.W


def inherited_confinement_check():
    result = entry.run()
    assert result == json.loads(entry.TARGET.read_text('utf8'))
    return dict(entry_exactly_reproduced=True, inherited_entry=result,
                full_gap_is_analytic_not_numerically_estimated=True)


def original_contact_check():
    q = model.original.shared_source(8)
    x, y, z = np.moveaxis(q['grid'], -1, 0)
    psi = 1.1 + .05*np.cos(x) + .02*np.sin(y)
    shape = model.shape_field(q['grid'])
    base = model.energy(model.graph_data(q, psi), shape)
    powers = dict(scalar_kinetic=-6, onsite=6, gradient=2, electric=-2, magnetic=-2)
    contact = sum(powers[k]**2 * base[k] for k in powers)
    rows = []
    for h in (.02, .01, .005):
        plus = model.energy(model.graph_data(q, psi*np.exp(h)), shape)['total']
        minus = model.energy(model.graph_data(q, psi*np.exp(-h)), shape)['total']
        observed = (plus - 2*base['total'] + minus)/h**2
        rows.append(dict(step=h, relative_contact_error=abs(observed-contact)/contact))
    assert rows[-1]['relative_contact_error'] < rows[0]['relative_contact_error']/12
    assert 0 < contact <= 36*base['total']
    return dict(original_full_classical_Gauss_source=True, contact_coefficient=contact,
                finite_difference_rows=rows, ground_E0_hessian_not_computed=True,
                contact_is_not_ground_energy_hessian=True)


def resolvent_remainder_check():
    # Exact rational identity for x=hbar*omega/Delta. These x values are
    # dimensionless kernel arguments, NOT invented original eigenvalues.
    eta = Fraction(3, 5)
    rows = []
    for x in (Fraction(1, 10), Fraction(1, 5), Fraction(2, 5), eta):
        remainder = 1/(1-x*x) - 1 - x*x
        exact = x**4/(1-x*x)
        bound = x**4/(1-eta*eta)
        assert remainder == exact and 0 <= remainder <= bound
        rows.append(dict(x=str(x), exact_remainder=str(exact), uniform_bound=str(bound)))
    # Inverse-gap weights after the form estimate: no spectrum is sampled.
    c, delta = Fraction(12, 5), Fraction(3, 10)
    weights = []
    for p in (1, 2, 3, 5):
        bound = delta**(1-p) + c*delta**(-p)
        values = [(d+c)/d**p for d in (delta, 2*delta, 10*delta, 1000*delta)]
        assert max(values) <= bound
        weights.append(dict(power=p, exact_bound=str(bound)))
    return dict(exact_rational_kernel_rows=rows, inverse_gap_scalar_bounds=weights,
                low_frequency_error_is_linear_response_only=True,
                numerical_ground_spectrum_used=False)


def record_escape_check():
    rows = []
    for nodes in (48, 80, 128):
        h, s, p, kss = previous.packet_data(nodes)
        lf, df, lc, dc = records.instruments(s)
        result = {}
        for name, L in (('fine', lf), ('direct_parity', lc)):
            means = np.sum(p[..., None]*L, axis=(0, 1))
            second = np.sum(p[..., None]*L*L, axis=(0, 1))
            escape = 1-float(means@means)
            variance = float(np.sum(second-means*means))
            assert escape > 0 and abs(escape-variance) < 8e-15
            result[name] = dict(survival=float(means@means), escape=escape,
                                sum_kraus_variances=variance)
        rows.append(dict(nodes=nodes, packet_results=result))
    assert abs(rows[-1]['packet_results']['fine']['escape']-rows[-2]['packet_results']['fine']['escape']) < 2e-14
    return dict(rows=rows, wavefunction_is_original_compact_Gauss_packet=True,
                packet_not_claimed_to_be_ground_state=True,
                strict_ground_escape_is_proved_using_ground_positivity=True)


def actual_record_force_check():
    h, s, p, kss = previous.packet_data(100)
    lf, df, lc, dc = records.instruments(s)
    rows = []
    for name, derivative in (('fine', df), ('direct_parity', dc)):
        A = np.sum(derivative**2, axis=-1)
        Q = HBAR**2/(2*W)*float(np.sum(p*kss*A))
        source_jump = -6*Q
        errors = []
        for step in (.01, .005, .0025):
            plus = HBAR**2/(2*W*np.exp(6*step))*float(np.sum(p*kss*A))
            minus = HBAR**2/(2*W*np.exp(-6*step))*float(np.sum(p*kss*A))
            errors.append(abs((plus-minus)/(2*step)-source_jump))
        assert Q > 0 and errors[-1] < errors[0]/12
        rows.append(dict(instrument=name, original_injection=Q,
                         conformal_source_jump=source_jump, derivative_errors=errors))
    assert rows[0]['original_injection'] > rows[1]['original_injection']
    return dict(rows=rows, original_packet_not_ground=True,
                exact_identity_holds_on_any_finite_form_energy_state=True,
                postmeasurement_source_not_replaced_by_vacuum_source=True)


def adm_signature_check():
    basis = []
    for i in range(3):
        a = np.zeros((3, 3)); a[i, i] = 1; basis.append(a)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        a = np.zeros((3, 3)); a[i, j] = a[j, i] = 1/np.sqrt(2); basis.append(a)
    gamma = model.shape_exp(np.array([[.17, .04, -.02], [.04, -.12, .03], [-.02, .03, .08]]))
    inv = np.linalg.inv(gamma)
    def bilinear(a, b):
        return float(np.trace(inv@a@inv@b) - np.trace(inv@a)*np.trace(inv@b))
    D = np.array([[bilinear(a, b) for b in basis] for a in basis])
    eigen = np.linalg.eigvalsh(D)
    conformal = bilinear(gamma, gamma)
    assert sum(eigen < 0) == 1 and sum(eigen > 0) == 5
    assert abs(conformal+6) < 2e-14
    # Pullbacks by the actual six-direction metric parameterization preserve
    # inertia; this is not an arbitrary positive surrogate mass tensor.
    Q0 = np.array([[.17, .04, -.02], [.04, -.12, .03], [-.02, .03, .08]])
    columns = []
    for a in basis:
        v = (model.shape_exp(Q0+1e-5*a)-model.shape_exp(Q0-1e-5*a))/(2e-5)
        columns.append(np.array([np.sum(v*b) for b in basis]))
    J = np.stack(columns, axis=1)
    pulled = np.linalg.eigvalsh(J.T@D@J)
    assert sum(pulled < 0) == 1 and sum(pulled > 0) == 5
    return dict(ADM_metric_velocity_form_eigenvalues=eigen.tolist(),
                conformal_value=conformal, actual_coordinate_pullback_eigenvalues=pulled.tolist(),
                ground_induced_inertia_PSD_is_analytic=True,
                excludes_only_unconstrained_six_component_direct_identification=True,
                constrained_GR_not_excluded=True)


def run():
    checks = (inherited_confinement_check, original_contact_check, resolvent_remainder_check,
              record_escape_check, actual_record_force_check, adm_signature_check)
    evidence = {f.__name__: f() for f in checks}
    names = ('research_note_579.md', 'research_note_588.md', 'research_note_589.md', 'research_note_590.md',
             'joint_full_spatial_metric.py', 'joint_geometry_work_noise.py',
             'joint_record_source_compression.py', 'round591_drafts/gapped_source_entry.py',
             'round591_drafts/gapped_source_entry_results.json')
    return dict(round=591, tests_run=len(checks), failures=0, errors=0, evidence=evidence,
                dependency_hashes={n: hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original full finite-graph model. Analytic compact resolvent, unique positive Gauss ground, and convergent subgap linear geometric response. Actual read instruments excite the ground. Positive induced parameter inertia is not unreduced ADM inertia. No numerical full spectrum, nonlinear adiabatic error, fixed-hbar continuum, or Einstein derivation.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--write-results', action='store_true'); args = p.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            f.write(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8')) == result
    print(json.dumps(result, ensure_ascii=False, indent=2))
