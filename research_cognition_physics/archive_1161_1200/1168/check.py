"""1168: real old clock transport, flat-cylinder ambiguity, direction holdout.

The new common piecewise-flat transport permission is explicit input.
Default is read-only comparison; --write-once creates results.json once.
No dynamic source, free turnaround, clock reset, or unknown-rail measurement.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
RADIUS, SPEED, GAMMA = .25, .6, 1.25
THETA, CAPTURE = math.pi/3, math.pi/4
TOL = 3e-9
SOURCES = {
 'archive_1106_1119/1106/proof.md': '3cf66b629a49e00a6503b68e516f67bdd5b9fecd258447e90156b15cf87fb8cb',
 'archive_1106_1119/1112/check.py': '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f',
 'archive_1106_1119/1119/proof.md': 'b4a97154a8302a49c4f4f6788f615909b4849820bfa42e3bf11574c16808008d',
 'archive_1120_1158/1158/proof.md': '8d05617c800014afd8965eaf51874fa3235e16cf2e49a2918ad00ebe5a8e7608',
 'archive_1161_/1166/proof.md': '02335c08c46865b438b075a0d8475e143fc113e01eb711f725b96e4290bb0812',
 'archive_1161_/1166/results.json': 'd1bf0dc11287685e6eb54b0b4c6a9af36c11b7bca82a51af86d8fd57ec3e7cf3',
 'archive_1161_/1167/acceptance.json': 'cab43b4c5421791c800cebfe012c8d42ebf8358da8050752167f340b51b6b74a',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def small(x, label):
    error = float(np.linalg.norm(x))
    assert error < TOL, (label, error)
    return error


def rotation(angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.array([[c, -s], [s, c]])


def polar(radius, angle):
    return radius*np.array([math.cos(angle), math.sin(angle)])


def trace_cone(start, target, winding):
    """Real straight-line face intersections; no target length/time argument."""
    assert (start+winding) % 4 == target
    point = polar(RADIUS, THETA/2)
    target_image = polar(RADIUS, THETA/2+winding*THETA)
    direction = target_image-point
    direction /= np.linalg.norm(direction)
    initial = direction.copy()
    cell, segments, elapsed, proper = start, [], 0., 0.
    for _ in range(8):
        candidates = []
        for side, angle in [('lower', 0.), ('upper', THETA)]:
            ray = polar(1., angle)
            matrix = np.column_stack([direction, -ray])
            if abs(np.linalg.det(matrix)) < 1e-12:
                continue
            step, radial = np.linalg.solve(matrix, -point)
            if step > 1e-11 and radial > 1e-11:
                candidates.append((float(step), side))
        if cell == target:
            difference = polar(RADIUS, THETA/2)-point
            step = float(difference@direction)
            if step > 1e-11 and np.linalg.norm(difference-step*direction) < 1e-10:
                candidates.append((step, 'actual_target'))
        assert candidates, (start, target, cell, point, direction)
        length, side = min(candidates)
        finish = point+length*direction
        coordinate = length/SPEED
        sigma = length/(GAMMA*SPEED)
        dz = np.array([coordinate, *(length*direction), 0.])
        u = np.array([GAMMA, *(GAMMA*SPEED*direction), 0.])
        small(dz-sigma*u, 'one local actual SR free leg')
        small(dz[0]**2-dz[1:]@dz[1:]-sigma**2, 'same proper interval')
        nearest = point+np.clip(-float(point@direction), 0., length)*direction
        segments.append(dict(cell=cell, entry=point.tolist(), exit=finish.tolist(),
            direction=direction.tolist(), length=length, coordinate_duration=coordinate,
            proper_duration=sigma, starts_at=elapsed, finishes_at=elapsed+coordinate,
            termination=side, smallest_radius=float(np.linalg.norm(nearest))))
        elapsed += coordinate
        proper += sigma
        if side == 'actual_target':
            small(finish-polar(RADIUS, THETA/2), 'actual fixed receiving station')
            break
        step_cell = 1 if side == 'upper' else -1
        transform = rotation(-step_cell*THETA)
        moved_point, moved_direction = transform@finish, transform@direction
        # Each seam acts on actual endpoints and the very same tangent.
        small(np.linalg.norm(moved_point)-np.linalg.norm(finish), 'shared face radial match')
        small(np.linalg.norm(moved_direction)-1., 'same speed at face')
        before_u = np.r_[GAMMA, GAMMA*SPEED*direction, 0.]
        after_u = np.r_[GAMMA, GAMMA*SPEED*moved_direction, 0.]
        lorentz = np.eye(4)
        lorentz[1:3, 1:3] = transform
        small(after_u-lorentz@before_u, 'same complete four-P coefficient at seam')
        eta = np.diag([1., -1., -1., -1.])
        small(lorentz.T@eta@lorentz-eta, 'local face isometry')
        segments[-1]['next_cell'] = (cell+step_cell) % 4
        segments[-1]['face_rotation'] = transform.tolist()
        point, direction = moved_point, moved_direction
        cell = (cell+step_cell) % 4
    else:
        raise AssertionError('finite route did not reach fixed endpoint')
    return dict(start=start, target=target, signed_cell_steps=winding,
        initial_direction=initial.tolist(), final_direction=direction.tolist(),
        segments=segments, length=sum(seg['length'] for seg in segments),
        coordinate_time=elapsed, proper_time=proper)


def cylinder_route(start, target, winding):
    circumference = 2*math.sqrt(3)*RADIUS
    station = lambda j: np.array([j*circumference/4, (-1)**j*RADIUS/4, 0.])
    source = station(start)
    destination = station(target)
    destination[0] = source[0]+winding*circumference/4
    displacement = destination-source
    length = float(np.linalg.norm(displacement))
    return dict(start=start, target=target, signed_cell_steps=winding,
        cover_source=source.tolist(), cover_target=destination.tolist(),
        direction=(displacement/length).tolist(), length=length,
        coordinate_time=length/SPEED, proper_time=length/(GAMMA*SPEED))


def group_geometry():
    cone, cylinder = [], []
    largest_match = 0.
    for start in range(4):
        for target in range(start+1, 4):
            steps = (target-start) % 4
            choices = [1] if steps == 1 else ([-1] if steps == 3 else [2, -2])
            for step in choices:
                a, b = trace_cone(start, target, step), cylinder_route(start, target, step)
                # These expressions are verification only, after ray tracing.
                expected = RADIUS if abs(step) == 1 else math.sqrt(3)*RADIUS
                small(a['length']-expected, 'closed-form check after actual intersections')
                largest_match = max(largest_match, small(a['length']-b['length'], 'cone-cylinder same distance'),
                    small(a['proper_time']-b['proper_time'], 'same full clock duration'))
                cone.append(a); cylinder.append(b)
    min_radius = min(seg['smallest_radius'] for route in cone for seg in route['segments'])
    small(min_radius-RADIUS/2, 'all actual paths stay away from excluded axis')
    assert len(cone) == 8
    # Reverse experiments are independent predeployments, not a physical turn.
    reverse_error = max(small(trace_cone(a['target'], a['start'], -a['signed_cell_steps'])['length']-a['length'],
                              'independently prepared reverse route') for a in cone)
    return dict(radius=RADIUS, wedge_angle=THETA, cone_total_angle=4*THETA,
        cylinder_circumference=2*math.sqrt(3)*RADIUS, speed=SPEED, gamma=GAMMA,
        cone_routes=cone, cylinder_routes=cylinder, distinct_station_pairs=6,
        route_count_including_diagonal_alternatives=8, smallest_radius=min_radius,
        complete_time_match_error=largest_match, independently_prepared_reverse_error=reverse_error,
        times_computed_from_face_intersections=True), cone, cylinder


def group_clock(cone, cylinder):
    path = BASE/'archive_1106_1119/1112/check.py'
    spec = importlib.util.spec_from_file_location('old1112_for1168', path)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    material = old.Material(1)
    hd = material.baseline*np.eye(288)+material.hfk
    hin = material.baseline*np.eye(288)+np.kron(material.hc, np.eye(72))
    values, vectors = np.linalg.eigh(hd)
    incoming = material.encode_joint(np.eye(3, dtype=complex)[:, None, :], 0.).reshape(288, 3)
    capture = material.gamma@old.exp_h(hin, CAPTURE)
    captured = capture@incoming
    capture_energy_error = small(hd@capture-capture@hin, 'original complete capture energy')
    transformed = vectors.conj().T@captured

    def output(s):
        return vectors@(np.exp(-1j*values*s)[:, None]*transformed)

    mask = np.zeros((4, 8, 3, 3), bool)
    mask[:, [1, 3, 5, 7], :, :] = True
    mask = mask.reshape(-1)

    def effect(v):
        return v[mask].conj().T@v[mask]

    rows, states = [], []
    max_effect = max_gram = max_match = max_phase = 0.
    for a, b in zip(cone, cylinder):
        v, w = output(a['proper_time']), output(b['proper_time'])
        e = effect(v)
        p = math.sin((CAPTURE+a['proper_time'])/2)**6
        station = output(a['coordinate_time'])
        ps = math.sin((CAPTURE+a['coordinate_time'])/2)**6
        max_effect = max(max_effect, small(e-p*np.eye(3), 'all input original clock effect'),
            small(effect(station)-ps*np.eye(3), 'same original stationary clock effect'))
        max_gram = max(max_gram, small(v.conj().T@v-np.eye(3), 'complete unknown Q/R isometry'))
        max_match = max(max_match, small(v-w, 'same complete internal output, matched roles'))
        inferred_s = 2*math.asin(float(e[0, 0].real)**(1/6))-CAPTURE
        inferred_t = 2*math.asin(float(effect(station)[0, 0].real)**(1/6))-CAPTURE
        small(inferred_s-a['proper_time'], 'proper time from actual Mr')
        small(inferred_t-a['coordinate_time'], 'station time from actual Mr')
        small(inferred_s/inferred_t-.8, 'common matter-clock rate')
        phase_coeff = sum(GAMMA*seg['coordinate_duration']-
                          GAMMA*SPEED*seg['length'] for seg in a['segments'])
        max_phase = max(max_phase, small(values*(phase_coeff-a['proper_time']),
                                         'complete E dt minus P dx equals HD proper time'))
        ages = [CAPTURE+a['proper_time'], CAPTURE+a['coordinate_time']]
        assert all(0 < age < math.pi for age in ages)
        rows.append(dict(start=a['start'], target=a['target'], signed_steps=a['signed_cell_steps'],
            moving_probability=p, stationary_probability=ps, inferred_proper_time=inferred_s,
            inferred_station_time=inferred_t, inferred_distance=GAMMA*SPEED*inferred_s,
            total_clock_ages=ages))
        states.append(v)
    return dict(material_dimension=288, baseline_mass=material.baseline,
        positive_complete_mass_spectrum=[float(values.min()), float(values.max())],
        capture_duration=CAPTURE, capture_energy_intertwiner_error=capture_energy_error,
        all_qutrit_effect_error=max_effect, complete_isometry_error=max_gram,
        cone_cylinder_complete_internal_matrix_difference=max_match,
        complete_four_P_phase_error=max_phase, task_readouts=rows,
        largest_two_flight_plus_capture_station_age=CAPTURE+2*max(a['coordinate_time'] for a in cone),
        actual_roundtrip_or_clock_reset_claimed=False), states, output, effect


def group_rail(cone, states):
    indices = [i for i, r in enumerate(cone) if r['start'] == 0 and r['target'] == 2]
    assert len(indices) == 2
    vp, vm = [states[i] for i in indices]
    small(vp-vm, 'equal diagonal proper time, complete old matrices')
    columns = [np.kron(vp[:, q], vm[:, 0]) for q in (1, 2)]
    columns += [np.kron(vp[:, 0], vm[:, q]) for q in (1, 2)]
    f = np.column_stack(columns)
    gram = small(f.conj().T@f-np.eye(4), 'all unknown path, qubit and arbitrary reference')
    # Normalized maximally entangled input/reference has full reference I4/4.
    reference = small((f.conj().T@f).T/4-np.eye(4)/4, 'retained arbitrary-reference identity')
    coherent = np.array([1., 0., 0., 1.])/math.sqrt(2)
    density = np.outer(coherent, coherent)
    dephased = density.copy()
    dephased[:2, 2:] = 0; dephased[2:, :2] = 0
    trace_difference = float(np.sum(abs(np.linalg.eigvalsh(density-dephased)))/2)
    small(trace_difference-.5, 'different channel when path coherence is deleted')
    off_diagonal_norm = float(np.linalg.norm(f[:, 0])*np.linalg.norm(f[:, 3]))
    small(off_diagonal_norm-1., 'full cross-path matrix unit retained')
    return dict(complete_two_material_dimension=288**2, open_rail_qubit_dimension=4,
        output_shape=list(f.shape), all_input_gram_error=gram,
        reference_reduced_state_error=reference, cross_path_matrix_unit_HS_norm=off_diagonal_norm,
        full_output_trace_distance_after_actual_path_dephasing=trace_difference,
        both_vac_carriers_and_receiver_velocity_types_retained=True,
        equal_internal_maps_do_not_mean_all_geometric_types_equal=True,
        no_new_cross_mode_interference_or_free_turning_permission=True)


def medium_proper(route, points):
    nodes, weights = np.polynomial.legendre.leggauss(points)
    total = 0.
    for seg in route['segments']:
        length = seg['length']
        x = np.array(seg['entry'])[None, :]+((nodes+1)*length/2)[:, None]*np.array(seg['direction'])
        radii = np.linalg.norm(x, axis=1)
        total += float(length/(2*SPEED)*np.dot(weights, np.sqrt(1-(2/3)*SPEED**2*radii)))
    return total


def group_identifiability(cone, cylinder, output, effect):
    diagonal = [i for i, row in enumerate(cone) if row['start'] == 0 and row['target'] == 2]
    cone_n = [rotation(-THETA/2)@np.array(cone[i]['initial_direction']) for i in diagonal]
    cylinder_n = [np.array(cylinder[i]['direction']) for i in diagonal]
    cone_delta = float(np.linalg.norm(sum(cone_n)))
    cylinder_delta = float(np.linalg.norm(sum(cylinder_n)))
    small(cone_delta-math.sqrt(3), 'actual cone initial directions')
    small(cylinder_delta, 'actual cylinder diagonal paths are antipodal')
    a, eps_w = .05, .001
    leg = SPEED*a
    assert leg < RADIUS*math.sin(THETA/2)
    eps_direction = 2*eps_w/leg
    eps_delta = 2*eps_direction
    intervals = [[-eps_delta, eps_delta], [cone_delta-eps_delta, cone_delta+eps_delta]]
    assert intervals[0][1] < intervals[1][0]
    # Literal measured displacements in the same small local SR chart.
    trial = [leg*n+eps_w*np.array([1., 0.]) for n in cone_n]
    measured = float(np.linalg.norm(sum(w/np.linalg.norm(w) for w in trial)))
    assert abs(measured-cone_delta) <= eps_delta

    fit = RADIUS*np.array([[-math.sqrt(3)/2, 0.], [0., .5], [math.sqrt(3)/2, 0.], [0., -.5]])
    distances = np.linalg.norm(fit[:, None, :]-fit[None, :, :], axis=2)
    sides = [(0, 1), (1, 2), (2, 3), (0, 3)]
    small(np.array([distances[i, j] for i, j in sides])-RADIUS, 'five-task Euclidean sides fit')
    small(distances[0, 2]-math.sqrt(3)*RADIUS, 'five-task first diagonal fit')
    small(distances[1, 3]-RADIUS, 'Euclidean holdout bound attained')
    violation = 2*(math.sqrt(3)*RADIUS)**2-4*RADIUS**2
    small(violation-2*RADIUS**2, 'strict global Euclidean incompatibility')

    lo, hi = .01, .15
    eps_p, eps_capture, eps_k, eps_other = 1e-5, 1e-5, 1e-4, 2e-4
    lip = 1/(3*lo**(5/6)*math.sqrt(1-hi**(1/3)))
    smax = max(route['proper_time'] for route in cone)
    eps_s = lip*eps_p+eps_capture
    eps_distance = (GAMMA*SPEED+eps_k)*eps_s+smax*eps_k+eps_other
    robust = 2*(math.sqrt(3)*RADIUS-eps_distance)**2-4*(RADIUS+eps_distance)**2
    assert robust > 0
    for route in cone:
        for duration in (route['proper_time'], route['coordinate_time']):
            assert lo < math.sin((CAPTURE+duration)/2)**6 < hi
    optical = []
    for index in (0, diagonal[0]):
        route = cone[index]
        medium = medium_proper(route, 96)
        small(medium-medium_proper(route, 48), 'finite smooth optical integral cross-check')
        p_flat = float(effect(output(medium))[0, 0].real)
        p_cone = float(effect(output(route['proper_time']))[0, 0].real)
        bound = route['coordinate_time']*(math.sqrt(1-(2/3)*SPEED**2*RADIUS)-.8)
        assert medium-route['proper_time'] >= bound-1e-12
        assert p_flat-p_cone > 2*eps_p
        optical.append(dict(start=route['start'], target=route['target'], same_coordinate_time=route['coordinate_time'],
            cone_moving_proper_time=route['proper_time'], flat_medium_moving_proper_time=medium,
            proper_time_gap_lower_bound=bound, old_Mr_probability_cone=p_cone,
            old_Mr_probability_flat_medium=p_flat, probability_gap=p_flat-p_cone))
    return dict(direction_holdout=dict(local_coordinate_duration=a, local_free_length=leg,
        distance_to_nearest_cone_face=RADIUS*math.sin(THETA/2), cone_initial_radial_tangent_directions=[n.tolist() for n in cone_n],
        cylinder_initial_directions=[n.tolist() for n in cylinder_n], cone_antipodal_defect=cone_delta,
        cylinder_antipodal_defect=cylinder_delta, displacement_error_contract=eps_w,
        direction_error_bound=eps_direction, defect_error_bound=eps_delta,
        separated_prediction_intervals=intervals, example_perturbed_defect=measured),
        global_Euclidean_screen=dict(five_task_fit_points=fit.tolist(), fitted_holdout_distance=float(distances[1, 3]),
            cone_and_cylinder_holdout_distance=math.sqrt(3)*RADIUS, squared_sum_violation=violation),
        finite_accuracy=dict(probability_window=[lo, hi], total_probability_error_contract=eps_p,
            capture_age_error=eps_capture, gamma_v_error=eps_k, other_endpoint_seam_unit_length_error=eps_other,
            inverse_clock_Lipschitz_bound=lip, total_length_error_bound=eps_distance,
            positive_robust_Euclidean_separation_margin=robust,
            these_are_declared_independent_implementation_errors_not_empirical_estimates=True),
        single_flat_optical_medium_countercheck=optical,
        optical_massive_path_and_proper_clock_rule_are_extra_conditional_inputs=True,
        no_claim_to_exclude_all_flat_topologies_or_forced_paths=True)


def compare(actual, saved):
    if isinstance(saved, dict):
        assert actual.keys() == saved.keys()
        for key in saved:
            compare(actual[key], saved[key])
    elif isinstance(saved, list):
        assert len(actual) == len(saved)
        for a, b in zip(actual, saved):
            compare(a, b)
    elif isinstance(saved, float):
        assert math.isclose(actual, saved, rel_tol=2e-9, abs_tol=3e-10), (actual, saved)
    else:
        assert actual == saved, (actual, saved)


def calculate():
    for path, expected in SOURCES.items():
        assert digest(BASE/path) == expected, path
    geometry, cone, cylinder = group_geometry()
    clock, states, output, effect = group_clock(cone, cylinder)
    rail = group_rail(cone, states)
    diagnostics = group_identifiability(cone, cylinder, output, effect)
    for path, expected in SOURCES.items():
        assert digest(BASE/path) == expected, ('source changed', path)
    return dict(schema='round1168_complete_clock_flat_topology_direction_v1', status='PASS',
        check_groups=4, code_sha256=digest(Path(__file__)), source_sha256=SOURCES,
        group1_shared_faces_and_flat_cylinder=geometry,
        group2_same_complete_old_clock=clock,
        group3_retained_unknown_path_and_reference=rail,
        group4_direction_holdout_and_false_positives=diagnostics,
        scope=dict(new_adopted_axioms=0, new_empirical_calibration=0, GR_hard_deliverables_completed=0,
            additional_candidate_permission='K1168 common static piecewise-flat actual transport',
            inherited_global_SR_A_not_claimed_on_quotient_domain=True,
            six_timing_and_complete_internal_channels_do_not_identify_geometry=True,
            local_direction_distinguishes_only_the_stated_cone_and_cylinder_models=True,
            no_dynamic_source_geometry_update_or_Einstein_equation_claim=True,
            no_nonzero_local_Riemann_claim_in_propagation_domain=True))


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
        mode = 'CREATED_ONCE'
    else:
        compare(result, json.loads(path.read_text(encoding='utf-8')))
        mode = 'READ_ONLY_MATCH'
    direction = result['group4_direction_holdout_and_false_positives']['direction_holdout']
    print(json.dumps(dict(status='PASS', mode=mode, groups=4,
        geometry_match=result['group1_shared_faces_and_flat_cylinder']['complete_time_match_error'],
        complete_clock_match=result['group2_same_complete_old_clock']['cone_cylinder_complete_internal_matrix_difference'],
        cone_direction_defect=direction['cone_antipodal_defect'],
        cylinder_direction_defect=direction['cylinder_antipodal_defect'],
        new_adopted_axioms=0, GR_hard_deliverables_completed=0), ensure_ascii=False))


if __name__ == '__main__':
    main()
