"""1165: full finite material spectra versus common isospectral transport.

Default recomputes and compares the saved result read-only. --write-once
creates results. No new geometry instrument or physical W is implemented.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCES = {
    '../1161/check.py': '52edf7a1dfc530996daa88ee3715efa936f7c5845eefe067124a544688621119',
    '../1164/proof.md': 'c20ace2c8163d278933e2f3232d2df791fe60d9046a9d14fac5b1d6dae8e68cc',
    '../1164/check.py': 'dd11c2df91bff7ceb3049673409b4d4de560ee6eb92f04d26367cb1a111742dd',
    '../1164/results.json': '24ef020b4b33538dc26bab116de2e13651b01aba0dc3264dc45e7eb0da55366e',
}
TOL = 2e-9


def small(value, label):
    error = float(np.linalg.norm(value))
    assert error < TOL, (label, error)
    return error


def energies(q):
    ex = math.sqrt(9 + 16*math.exp(-2*q))
    ey = math.sqrt(9 + 16*math.exp(2*q))
    delta = 32*math.sinh(2*q)/(ex+ey)
    emax = math.sqrt(9 + 16*math.exp(2*abs(q)))
    distance = 32*math.expm1(2*abs(q))/(emax+5)
    return ex, ey, delta, distance


def compare(actual, expected):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys()
        for key in expected:
            compare(actual[key], expected[key])
    elif isinstance(expected, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            compare(a, b)
    elif isinstance(expected, float):
        assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=3e-10), (actual, expected)
    else:
        assert actual == expected, (actual, expected)


def calculate():
    for name, expected in SOURCES.items():
        assert hashlib.sha256((HERE/name).read_bytes()).hexdigest() == expected, name
    spec = importlib.util.spec_from_file_location('round1161_source_for_1165', HERE/'../1161/check.py')
    source = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    old, data = source.reuse()  # Runs and compares the original frozen 1160.
    hd = data['m']['hout']
    baseline = data['m']['baseline']
    ev = np.linalg.eigvalsh(hd)
    assert hd.shape == (288, 288)
    d_spectrum_error = small(ev-baseline-np.repeat(np.arange(4), 72), 'full D spectrum')
    spectrum0 = np.sort(np.repeat(ev+13, 4*12))
    assert len(spectrum0) == 13824
    old_multiplicities = [int(np.count_nonzero(abs(spectrum0-(baseline+13+n)) < TOL)) for n in range(4)]
    assert old_multiplicities == [3456]*4
    source_gram = max(small(state.reshape(-1, 2).conj().T @ state.reshape(-1, 2)-np.eye(2),
                            'same complete unknown input') for state in data['outputs'].values())

    # Full eigenvalue lists come from the original physical D spectrum and all
    # AB/M/E/L factors. No preparation-code compression or dense mother needed.
    q_values = [-.01, -.003, 0., .003, .01]
    orbit_rows = []
    scale_rows = []
    split_rows = []
    max_orbit_error = max_scale_error = 0.
    for q in q_values:
        ex, ey, delta, distance = energies(q)
        ab = np.array([2*ex, ex+ey, ex+ey, 2*ey])
        spectrumq = np.sort((ev[:, None, None]+ab[None, :, None]+3*np.ones((1, 1, 12))).ravel())
        assert len(spectrumq) == 13824
        sorted_lower_bound = float(np.max(abs(spectrumq-spectrum0)))
        identity_upper_bound = float(np.max(abs(ab-10)))
        upper_endpoint_bound = float(spectrumq[-1]-spectrum0[-1])
        for value in (sorted_lower_bound, identity_upper_bound, upper_endpoint_bound):
            max_orbit_error = max(max_orbit_error, small(value-distance, 'exact full orbit distance'))
        orbit_rows.append(dict(q=q, analytic_distance=distance, sorted_spectral_bound=sorted_lower_bound,
                               identity_upper_bound=identity_upper_bound,
                               top_spectral_lower_bound=upper_endpoint_bound))
        if q == 0:
            continue
        emin, emax = min(ex, ey), max(ex, ey)
        counts = []
        for n in range(4):
            band = spectrumq[n*3456:(n+1)*3456]
            small(band[0]-(baseline+n+3+2*emin), 'band lower edge')
            small(band[-1]-(baseline+n+3+2*emax), 'band upper edge')
            values = [baseline+n+3+2*emin, baseline+n+3+ex+ey, baseline+n+3+2*emax]
            count = [int(np.count_nonzero(abs(band-value) < TOL)) for value in values]
            assert count == [864, 1728, 864]
            counts.append(count)
            if n < 3:
                assert spectrumq[(n+1)*3456]-band[-1] > 0
        split_rows.append(dict(q=q, multiplicities_by_D_level=counts,
                               interband_gap=1-2*abs(delta)))
        for scale in [.9, 1., 1.1]:
            differences = spectrumq-scale*spectrum0
            shift = float((differences.max()+differences.min())/2)
            best_shifted_error = float(np.max(abs(differences-shift)))
            # Noninterlacing gives this stronger finite calibration. The proof
            # only needs every old repeated block to span both new endpoints.
            expected = abs(delta)+1.5*abs(1-scale)
            max_scale_error = max(max_scale_error, small(best_shifted_error-expected, 'affine full spectra'))
            if scale == 1:
                max_scale_error = max(max_scale_error, small(shift-(ex+ey-10), 'common optimal offset'))
                max_scale_error = max(max_scale_error, small(best_shifted_error-abs(delta), 'affine optimum'))
            scale_rows.append(dict(q=q, scale=scale, optimum_offset_for_scale=shift,
                                   sorted_spectral_error=best_shifted_error,
                                   lower_bound_all_positive_scales=abs(delta)))

    interval_bound = Fraction(64, 3)*Fraction(1, 100)/(1-Fraction(1, 50))
    assert interval_bound == Fraction(32, 147) < 1

    # The 16 spectral representatives have the physical multiplicity 864 each.
    # This is a calibration of the analytic projector argument, not a new model.
    hrep = np.diag(np.repeat(baseline+13+np.arange(4), 4))
    s_ab = np.array([-32/5, 0., 0., 32/5])
    srep = np.diag(np.tile(s_ab, 4))
    rng = np.random.default_rng(1165)
    z = rng.normal(size=(16, 16))+1j*rng.normal(size=(16, 16))
    generator = (z+z.conj().T)/2
    generator /= np.linalg.norm(generator, 2)
    tangent = -1j*(generator@hrep-hrep@generator)
    tangent_block_error = 0.
    affine_block_radii = []
    for n in range(4):
        block = slice(4*n, 4*n+4)
        tangent_block_error = max(tangent_block_error, small(tangent[block, block], 'zero degenerate tangent'))
        difference = srep-tangent-hrep/7+2*np.eye(16)/9
        block_eigenvalues = np.linalg.eigvalsh(difference[block, block])
        radius = float((block_eigenvalues[-1]-block_eigenvalues[0])/2)
        small(radius-32/5, 'unremoved traceless block radius')
        affine_block_radii.append(radius)
    h = 1e-5
    ex1, ey1, _, _ = energies(h)
    ex2, ey2, _, _ = energies(-h)
    derivative = (np.array([2*ex1, ex1+ey1, ex1+ey1, 2*ey1])
                  -np.array([2*ex2, ex2+ey2, ex2+ey2, 2*ey2]))/(2*h)
    derivative_error = float(np.max(abs(derivative-s_ab)))
    assert derivative_error < 2e-8

    # Bounded waits of the proposed boundary operators: all old material and
    # original output columns remain. This does not grant a new physical wait.
    approximation_rows = []
    for q in [.01, .0001, .000001]:
        ex, ey, delta, distance = energies(q)
        shifts = np.array([2*ex-10, ex+ey-10, ex+ey-10, 2*ey-10])
        for duration in [.5, 2.]:
            phase_difference = np.exp(-1j*duration*shifts)-1
            norms = []
            for program in range(3):
                state = data['outputs'][program]
                difference = state*phase_difference[None, :, None, None, None, None]
                value = float(np.linalg.norm(difference.reshape(-1, 2), 2))
                assert value <= duration*distance+2e-12
                norms.append(value)
            approximation_rows.append(dict(q=q, duration=duration, energy_operator_distance=distance,
                                           affine_energy_distance=abs(delta),
                                           duhamel_isometry_upper_bound=duration*distance,
                                           original_output_column_norms=norms))
    slope_error = abs(energies(1e-6)[3]/1e-6-32/5)
    assert slope_error < 5e-6
    accepted = json.loads((HERE/'../1164/results.json').read_text(encoding='utf-8-sig'))
    for row in accepted['full_process']['samples']:
        small(energies(row['q'])[2]-row['delta'], 'same old shear dictionary')

    return dict(schema='round1165_isospectral_transport_distance_v1', status='PASS', check_groups=4,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), source_sha256=SOURCES,
        original_process_replay='1160 calculate matched frozen results through 1161 reuse',
        full_material=dict(dimension=13824, D_spectrum_error=d_spectrum_error,
                           old_level_multiplicities=old_multiplicities,
                           actual_unknown_input_gram_error=source_gram,
                           deformed_level_multiplicities=split_rows),
        unitary_orbit=dict(formula='2*(sqrt(9+16*exp(2*abs(q)))-5)',
                           maximum_identity_error=max_orbit_error, cases=orbit_rows),
        common_scale_and_tangent=dict(affine_optimum='abs(Ey-Ex) for abs(q)<=1/100',
            noninterlacing_two_delta_bound=str(interval_bound), affine_error=max_scale_error,
            cases=scale_rows, zero_block_commutator_error=tangent_block_error,
            unremoved_block_half_widths=affine_block_radii, source_derivative_error=derivative_error,
            source_slope=32/5),
        finite_precision=dict(cases=approximation_rows, distance_slope_error=slope_error,
                              nonzero_operator_lower_bound_is_not_measurement_probability=True),
        new_adopted_axioms=0, new_empirical_calibration=0,
        scope='Frozen finite boundary-material spectrum; no new geometry instrument, physical W, or GR completion.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    path = HERE/'results.json'
    if args.write_once:
        with path.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
    else:
        compare(result, json.loads(path.read_text(encoding='utf-8-sig')))
    print(json.dumps(dict(status='PASS', mode='WRITE_ONCE' if args.write_once else 'READ_ONLY_MATCH',
                         groups=4, full_material_dimension=13824,
                         orbit_error=result['unitary_orbit']['maximum_identity_error'],
                         affine_error=result['common_scale_and_tangent']['affine_error']), ensure_ascii=False))


if __name__ == '__main__':
    main()
