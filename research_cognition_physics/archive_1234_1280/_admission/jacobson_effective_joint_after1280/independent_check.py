"""Independent exact audit; no author imports, simulation, or field recomputation."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir())
AUTHOR_FILES = ('fixed_spec.md', 'proof.md', 'check.py', 'results.json')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def encode(x):
    return {'exact': str(x), 'decimal': float(x)}

def calculate():
    low = Q(4) * Q(499, 500)**2
    high = Q(72)
    d = Q(591628623, 124500500000)
    calibration, department = Q(1, 1000), Q(1, 5000)
    standard_deviation, noise_event = Q(1, 100000), Q(1, 1000)
    grid = Q(1, 10000)
    source_budget = 5 * department
    e = calibration + source_budget + noise_event + grid / 2
    failure = 3 * (standard_deviation / noise_event)**2
    assert e == Q(61, 20000) and failure == Q(3, 10000)
    assert Q(0) < d < Q(1, 8)

    # The stipulated 3-by-5 source projection matrix has rank one.
    # Check all independent department sign corners, not three separate sources.
    source_corners = []
    for signs in product((-1, 1), repeat=5):
        residual = sum((s * department for s in signs), Q(0))
        vector = (residual, residual, residual)
        assert max(abs(v) for v in vector) <= source_budget
        assert vector[0] - vector[1] == vector[1] - vector[2] == 0
        source_corners.append(vector)
    witness = Q(5, 10000)
    assert Q(0) < witness < source_budget

    # Derive the stronger joint-latent-parameter interval separation.
    joint_g = low / 4 - Q(7, 2) * e
    joint_t = (Q(1, 4) - 2*d) * low - (Q(7, 2) + 2*d) * e
    published_g = low / 4 - Q(15, 4) * e
    published_t = (Q(1, 4) - 2*d) * low - Q(15, 4) * e
    assert joint_g - published_g == e / 4 > 0
    assert joint_t - published_t == (Q(1, 4) - 2*d) * e > 0
    assert published_g > Q(94, 100) and published_t > Q(94, 100)
    pair_relative = d + 2*e/low
    precision_margin = (Q(1, 100) - d)*low - 2*e
    assert pair_relative < Q(1, 100) and precision_margin > 0
    assert (Q(1, 100) - d)*low - 2 < 0

    def feasible(data, fraction, use_geometry=True, use_thermal=True):
        xb, xy, xm = data
        left, right = max(low, xb-e), min(high, xb+e)
        if use_geometry:
            left, right = max(left, (xy-e)/fraction), min(right, (xy+e)/fraction)
        if use_thermal:
            left = max(left, (xm-e)/(fraction+d))
            right = min(right, (xm+e)/(fraction-d))
        return left <= right

    # Finite exact diagnostic over a conservative error polytope.  The universal
    # proof is the preceding affine interval inequality, not this enumeration.
    corners3 = tuple(product((-1, 1), repeat=3))
    accepted_cases = 0
    for b, msign, zsign, cs, ns, rs in product(
            (low, high), (-1, 1), (-1, 1), corners3, corners3, corners3):
        targets = (b, b, b*(1+msign*d))
        data = tuple(targets[j] + cs[j]*calibration + zsign*source_budget
                     + ns[j]*noise_event + rs[j]*grid/2 for j in range(3))
        assert feasible(data, Q(1))
        assert not feasible(data, Q(3, 4), use_thermal=False)
        assert not feasible(data, Q(3, 4), use_geometry=False)
        accepted_cases += 1
    assert accepted_cases == 4096

    # A finite push-forward of the normalized positive Gaussian is normalized.
    # Positive covariance pivots, range margin, and exact output cardinality.
    covariance_pivots = [standard_deviation**2] * 3
    assert all(v > 0 for v in covariance_pivots)
    good_max = high*(1+d)+e
    assert good_max < 100
    grid_count = int(Q(200)/grid)+1
    codes = grid_count+1  # one explicit overflow outcome
    bits = (codes-1).bit_length()
    assert 2**(bits-1) < codes <= 2**bits and bits == 21

    # Incidence-matrix conservation; no inference from energy to stress.
    names = ('clock_support', 'control', 'store', 'pointer', 'record', 'outgoing')
    start = tuple(map(Q, (2, 2, 4, 1, 1, 0)))
    edges = ((2, 1, Q(3, 10)), (1, 3, Q(1, 10)),
             (1, 4, Q(1, 5)), (1, 5, Q(1, 10)))
    changes = [Q(0)]*6
    for sender, receiver, transfer in edges:
        column = [Q(0)]*6
        column[sender], column[receiver] = -transfer, transfer
        assert sum(column) == 0
        changes = [a+b for a,b in zip(changes, column)]
    end = tuple(a+b for a,b in zip(start, changes))
    assert sum(start) == sum(end) == 10 and min(end) >= 0
    assert end == (Q(2), Q(19,10), Q(37,10), Q(11,10), Q(6,5), Q(1,10))

    # Cross-check only saved numbers after deriving them independently.
    author = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    assert author['fixed_spec_sha256'] == sha(HERE/'fixed_spec.md')
    expected = {'B_min':low, 'B_max_strict':high,
                'relative_reference_bound':d, 'channel_radius':e,
                'joint_failure_probability_upper':failure,
                'one_percent_margin':precision_margin,
                'geometry_omission_gap':published_g,
                'thermal_omission_gap':published_t}
    for key, value in expected.items():
        assert Q(author['bounds'][key]) == value, key
    assert author['finite_record'] == {'codes_per_channel':codes,
                                      'bits_per_channel':bits, 'total_bits':3*bits}
    for j, name in enumerate(names):
        assert Q(author['energy_ledger']['initial'][name]) == start[j]
        assert Q(author['energy_ledger']['final'][name]) == end[j]

    files = [HERE/name for name in AUTHOR_FILES]+[Path(__file__).resolve()]
    return {
        'status':'PASS_WITH_SCOPE',
        'independent_diagnostic_groups':1,
        'formal_scientific_group_increment':0,
        'files_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files},
        'exact_bounds':{key:encode(value) for key,value in expected.items()},
        'joint_latent_parameter':{
            'geometry_gap_stronger':encode(joint_g),
            'thermal_gap_stronger':encode(joint_t),
            'geometry_extra_margin':encode(joint_g-published_g),
            'thermal_extra_margin':encode(joint_t-published_t),
            'corner_diagnostics':accepted_cases,
            'full_acceptance_and_both_omission_rejections':True},
        'shared_source':{'projection_rank':1, 'department_count':5,
                         'exact_corner_checks':len(source_corners),
                         'per_channel_budget':str(source_budget),
                         'nonzero_common_residual_witness':str(witness)},
        'finite_output':{'positive_covariance_pivots':[str(v) for v in covariance_pivots],
                         'good_event_maximum':encode(good_max),
                         'codes_per_channel':codes,'bits_per_channel':bits,'total_bits':3*bits,
                         'success_probability_lower':str(1-failure)},
        'pairing':{'worst_relative_error':encode(pair_relative),
                   'coarse_radius_one_does_not_certify_one_percent':True},
        'resource_ledger':{'initial':dict(zip(names,map(str,start))),
                           'final':dict(zip(names,map(str,end))),
                           'closed_total':'10','registered_outflow':'1/10'},
        'scope':{'public_effective_observation_contract':True,
                 'microscopic_povm_constructed':False,
                 'resource_energy_implies_stress_budget':False,
                 'source_and_geometric_transport_are_registered_inputs':True,
                 'gaussian_normalization_proves_apparatus_qualification':False,
                 'finite_output_implies_global_finite_field_resource':False,
                 'omission_is_an_alternative_physical_solution':False,
                 'whole_goal_automatically_completed':False}}

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args=parser.parse_args()
    data=calculate()
    destination=HERE/'independent_results.json'
    if args.write:
        with destination.open('x', encoding='utf-8', newline='\n') as handle:
            json.dump(data,handle,ensure_ascii=False,indent=2)
            handle.write('\n')
        mode='EXCLUSIVE_FIRST_WRITE'
    else:
        saved=json.loads(destination.read_text(encoding='utf-8'))
        assert saved==data, 'independent replay or source binding mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps({'status':data['status'],'mode':mode,
                      'corner_diagnostics':data['joint_latent_parameter']['corner_diagnostics'],
                      'geometry_gap':str(Q(data['exact_bounds']['geometry_omission_gap']['exact'])),
                      'thermal_gap':str(Q(data['exact_bounds']['thermal_omission_gap']['exact']))},
                     ensure_ascii=False))
