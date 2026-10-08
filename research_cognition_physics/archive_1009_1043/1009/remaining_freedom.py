"""1009: audit parameter fibers using the existing 956/981 finite material.

Analytic inference criteria are in the note. Numerical checks do not derive
the material Hamiltonian, the lapse coupling, or a physical apparatus.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
OUT = HERE / 'remaining_freedom_results.json'
MATERIAL = BASE / 'archive_956_989/956/native_material_interface.py'


def alpha(delta, d, w):
    return 2*d*delta/(delta*delta-w*w)


def beta(delta, d, w):
    return 2*d*delta*(delta*delta+w*w)/(delta*delta-w*w)**2


def parameters(delta, d):
    assert delta > 0 and 0 < d < .5
    return delta*(1-2*d)/(1-d), delta/2*math.sqrt(d/(1-d))


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a: compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b): compare(x, y)
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=2e-8, abs_tol=2e-12), (a, b)
    else:
        assert a == b, (a, b)


def run():
    spec = importlib.util.spec_from_file_location('frozen_material_956', MATERIAL)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    baseline = old.material(1., .1)
    lower = 1 + baseline['J']
    d0 = baseline['occupancy']
    static = 2*d0/lower
    upper = 2*lower
    frequency = .25*lower
    assert static*upper < 1
    cs = [old.annihilate(i) for i in range(4)]
    ns = [c.T@c for c in cs]
    sector = np.eye(16)[:, [i for i in range(16) if i.bit_count() == 2]]
    charge = sector.T@(ns[0]+ns[1]-ns[2]-ns[3])@sector/2
    rows = []
    spectral_error = material_error = static_error = source_error = 0.
    for delta in np.linspace(lower, upper, 17):
        delta = float(delta)
        d = static*delta/2
        U, hopping = parameters(delta, d)
        m = old.material(U, hopping)
        material_error = max(material_error, abs(m['occupancy']-d), abs(U+m['J']-delta))
        values, vectors = np.linalg.eigh(m['H'])
        ground = vectors[:, 0]
        gaps = values-values[0]
        weights = abs(vectors.conj().T@charge@ground)**2
        good = gaps > 1e-9
        for w in (0., frequency/2, frequency):
            spectral = float(np.sum(2*gaps[good]*weights[good]/(gaps[good]**2-w*w)))
            spectral_error = max(spectral_error, abs(spectral-alpha(delta, d, w)))
        static_error = max(static_error, abs(alpha(delta, d, 0)-static), abs(beta(delta, d, 0)-static))
        # Independent derivative of the inherited, fixed-coordinate-field kernel.
        step = 1e-4
        def kernel(N):
            return 2*d*N*delta/((N*delta)**2-frequency**2)
        derivative = (kernel(1-2*step)-8*kernel(1-step)+8*kernel(1+step)-kernel(1+2*step))/(12*step)
        source_error = max(source_error, abs(-derivative-beta(delta, d, frequency)))
        a = alpha(delta, d, frequency)
        r = a/static
        inferred_delta = math.sqrt(frequency**2*r/(r-1))
        inferred_d = static*inferred_delta/2
        assert abs(inferred_delta-delta) < 1e-12 and abs(inferred_d-d) < 1e-13
        # A genuine finite-time matrix comparison; this is a prescribed task,
        # not a constructed autonomous instrument or full parent realization.
        field, duration = .02, math.pi/lower
        output = old.evolution(m['H']-field*charge, duration)@ground
        excitation = float(1-abs(np.vdot(ground, output))**2)
        assert abs(np.linalg.norm(output)-1) < 1e-13
        rows.append(dict(delta=delta, oscillator_weight=d, U=U, hopping=hopping,
                         static_alpha=alpha(delta,d,0), dynamic_alpha=a,
                         dynamic_lapse_beta=beta(delta,d,frequency),
                         finite_pulse_excitation=excitation))
    assert spectral_error < 1e-12 and material_error < 1e-13
    assert static_error < 1e-14 and source_error < 1e-10
    a_min, a_max = rows[-1]['dynamic_alpha'], rows[0]['dynamic_alpha']
    b_min, b_max = rows[-1]['dynamic_lapse_beta'], rows[0]['dynamic_lapse_beta']
    # Endpoints are exact by monotonicity in Delta; grid only checks calibration.
    for row in rows:
        assert a_min-1e-14 <= row['dynamic_alpha'] <= a_max+1e-14
        assert b_min-1e-14 <= row['dynamic_lapse_beta'] <= b_max+1e-14
    r = rows[0]['dynamic_alpha']/static
    ratio_error = .002
    rlo, rhi = r-ratio_error, r+ratio_error
    identified = [math.sqrt(frequency**2*rhi/(rhi-1)),
                  math.sqrt(frequency**2*rlo/(rlo-1))]
    assert identified[0] < lower < identified[1]
    tiny_ratio = alpha(lower,d0,.01*lower)/static
    # If the ratio interval intersects 1, the data alone gives no finite
    # upper gap bound. A predeclared material prior still gives Delta < 1/A.
    assert tiny_ratio-ratio_error <= 1 < tiny_ratio+ratio_error
    # Independent finite differences check information ranks, not theoremhood.
    def observables(x, dynamic):
        D, d = x
        return np.array([alpha(D,d,0), beta(D,d,0)] +
                        ([alpha(D,d,frequency)] if dynamic else []))
    ranks=[]
    for dynamic in (False, True):
        x=np.array([lower,d0]); eps=1e-6
        jac=np.column_stack([(observables(x+eps*np.eye(2)[k],dynamic)-
                              observables(x-eps*np.eye(2)[k],dynamic))/(2*eps) for k in range(2)])
        ranks.append(int(np.linalg.matrix_rank(jac,tol=1e-8)))
    assert ranks == [1,2]
    pulse_difference = abs(rows[0]['finite_pulse_excitation']-rows[-1]['finite_pulse_excitation'])
    assert pulse_difference > 1e-5
    evidence = [Path(__file__), MATERIAL,
        BASE/'archive_956_989/981/native_response_source.py',
        BASE/'archive_956_989/research_note_981.md',
        BASE/'archive_702_741/research_note_708.md',
        BASE/'archive_923_934/research_note_930.md',
        BASE/'archive_935_955/research_note_954.md',
        BASE/'archive_990_1008/1008/overall_completion_results.json']
    return dict(round=1009, date='2026-10-08', all_scientific_checks_passed=True,
        fresh_model_calibration_groups=1, cumulative_test_groups=3787,
        input_classification='conditional dependency and task-relative identifiability',
        static_summary=static, gap_interval=[lower,upper], fixed_frequency=frequency,
        static_summary_rank=ranks[0], extended_summary_rank=ranks[1], rows=rows,
        max_errors=dict(material_inverse=material_error, spectral_response=spectral_error,
                        common_static_summary=static_error, inherited_source_derivative=source_error),
        exact_fiber_prediction_intervals=dict(alpha=[a_min,a_max], lapse_beta=[b_min,b_max]),
        minimax_absolute_errors=dict(alpha=(a_max-a_min)/2, lapse_beta=(b_max-b_min)/2),
        minimax_errors_over_static=dict(alpha=(a_max-a_min)/(2*static), lapse_beta=(b_max-b_min)/(2*static)),
        ratio_measurement_example=dict(ratio=r, absolute_error=ratio_error,
            inferred_gap_interval=identified, low_frequency_ratio=tiny_ratio,
            low_frequency_data_alone_gives_finite_upper_bound=False,
            existing_material_prior_upper_gap=1/static),
        finite_pulse=dict(field=.02, duration=math.pi/lower,
                          endpoint_probability_difference=pulse_difference,
                          evidence_kind='floating_point_finite_matrix_calibration'),
        limitations=dict(new_cognitive_axioms=False, material_parameters_derived_from_cognition=False,
            entire_cognitive_axiom_countermodel=False, physical_dimension_derived=False,
            full_source_recovered=False, autonomous_instrument_constructed=False,
            empirical_parameter_measurement=False, exact_finite_precision_identification=False,
            generative_goal_completed=False),
        source_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in evidence})


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2); stream.write('\n')
    else:
        compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps({k:result[k] for k in ('round','all_scientific_checks_passed',
        'static_summary_rank','extended_summary_rank','minimax_errors_over_static','finite_pulse')},ensure_ascii=False,indent=2))
