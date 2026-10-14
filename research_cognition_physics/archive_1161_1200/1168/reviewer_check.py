"""Read-only independent audit of 1168; not a new scientific round/group.

Independent oracles: universal-cover closest lifts; unfolded chord crossings;
spin-3/2 closed amplitudes and Lagrange spectral polynomials, without eigh.
Only this review script and the review receipts were authored by the reviewer.
"""
from __future__ import annotations
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import re
from urllib.parse import unquote
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
CORE = {
 '../research_note_1168.md': '9859deb6fb58547dd107b9c14624b941107de6f1ac41b059c3701345e883c441',
 'admission.md': '70b0331643ef7f2168010fcd1a13e2e27960eef083b257e83a4702b871e6abc1',
 'proof.md': 'c5674849ab457dd3b6a4bba2d0ce81b03263f94c1b3881303aa76821e46c7711',
 'check.py': 'fbd033b53b3df2484fd359867dee43e3bea7fd9c6e06983847e19ac261e74b6d',
 'results.json': '64822d73d9dadef64d53742c1468a380aa0ac58ba2bdd12e0f7385dc9adc0777',
 'NEXT.md': '32b722438b79555f31349a3370f06e24a628c2078b3ab86625ab2f50b838083e',
}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def small(value, label, tol=3e-10):
    x = float(np.linalg.norm(value))
    assert x < tol, (label, x)
    return x

def projector_polynomials(h):
    eye = np.eye(len(h), dtype=complex)
    ps = []
    for j in range(4):
        p = eye.copy()
        for k in range(4):
            if k != j:
                p = p @ ((h-k*eye)/(j-k))
        ps.append(p)
    return ps

def main():
    for name, expected in CORE.items():
        assert digest(HERE/name) == expected, name
    author = load('author1168_for_review', HERE/'check.py')
    old = load('old1112_for_independent1168', BASE/'archive_1106_1119/1112/check.py')
    saved = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    replay = author.calculate()
    stats = {'leaves': 0, 'floats': 0, 'maximum_absolute_float_difference': 0.}
    def compare(a, b):
        if isinstance(b, dict):
            assert a.keys() == b.keys()
            for key in b:
                compare(a[key], b[key])
        elif isinstance(b, list):
            assert len(a) == len(b)
            for x, y in zip(a, b):
                compare(x, y)
        else:
            stats['leaves'] += 1
            if isinstance(b, float):
                stats['floats'] += 1
                stats['maximum_absolute_float_difference'] = max(stats['maximum_absolute_float_difference'], abs(a-b))
                assert math.isclose(a, b, rel_tol=2e-9, abs_tol=3e-10)
            else:
                assert a == b
    compare(replay, saved)
    geometry = saved['group1_shared_faces_and_flat_cylinder']
    cone, cylinder = geometry['cone_routes'], geometry['cylinder_routes']
    r, v, gamma, theta, delta = .25, .6, 1.25, math.pi/3, math.pi/4
    alpha, circumference = 4*theta, 2*math.sqrt(3)*r
    station = lambda j: np.array([j*circumference/4, (-1)**j*r/4, 0.])
    max_lift = max_chord = max_crossing = max_direction = 0.
    lift_counts = []
    for start in range(4):
        for target in range(start+1, 4):
            # Outside this finite range the longitudinal distance is greater
            # than L; closest lifts are floor/ceil of -(s_j-s_i)/L.
            images = [station(target)-station(start)+np.array([m*circumference, 0., 0.]) for m in range(-2, 3)]
            norms = [float(np.linalg.norm(x)) for x in images]
            best = min(norms)
            nearest = [x/best for x, norm in zip(images, norms) if abs(norm-best) < 1e-12]
            rows = [x for x in cone if x['start'] == start and x['target'] == target]
            lift_counts.append(len(nearest))
            assert len(nearest) == len(rows)
            for row in rows:
                max_lift = max(max_lift, small(row['length']-best, 'nearest cylindrical lift'))
            # A non-apex geodesic on the punctured cone has angular sweep
            # strictly less than pi. Other windings cannot be shorter chords.
            sweeps = [(target-start)*theta+m*alpha for m in range(-2, 3)]
            lengths = [2*r*abs(math.sin(b/2)) for b in sweeps if abs(b) < math.pi-1e-12]
            max_chord = max(max_chord, small(min(lengths)-best, 'independent cone cover minimum'))
    for row in cone:
        w = row['signed_cell_steps']
        z0 = r*np.exp(1j*theta/2)
        z1 = r*np.exp(1j*(theta/2+w*theta))
        dz = z1-z0
        direction = dz/abs(dz)
        max_direction = max(max_direction, small(np.array([direction.real, direction.imag])-row['initial_direction'], 'unfolded initial direction'))
        # Cross-product intersection, rather than the author's 2x2 solve.
        cuts = []
        for m in range(-4, 5):
            ray = np.exp(1j*m*theta)
            den = (np.conj(ray)*dz).imag
            if abs(den) < 1e-14:
                continue
            t = -(np.conj(ray)*z0).imag/den
            z = z0+t*dz
            if 1e-12 < t < 1-1e-12 and (np.conj(ray)*z).real > 0:
                cuts.append(t)
        cuts = sorted(set(round(x, 14) for x in cuts))
        pieces = np.diff([0., *cuts, 1.])*abs(dz)
        max_crossing = max(max_crossing, small(pieces-np.array([s['length'] for s in row['segments']]), 'independent unfolded face intervals'))
        max_chord = max(max_chord, small(abs(dz)-row['length'], 'unfolded total chord'))
    assert lift_counts == [1, 2, 1, 1, 2, 1]
    cone_dirs = [np.array([-math.sqrt(3)/2, .5, 0.]), np.array([-math.sqrt(3)/2, -.5, 0.])]
    cylinder_dirs = [np.array([1., 0., 0.]), np.array([-1., 0., 0.])]
    angular_gap = [math.acos(float(cone_dirs[0]@cone_dirs[1])), math.acos(float(cylinder_dirs[0]@cylinder_dirs[1]))]
    small(angular_gap-np.array([math.pi/3, math.pi]), 'independent angular holdout')

    material = old.Material(1)
    h, hc, baseline = material.hfk, material.hc, material.baseline
    ps, pcs = projector_polynomials(h), projector_polynomials(hc)
    proj_error = max(small(sum(ps)-np.eye(288), 'spectral completeness'),
        max(small(h@p-j*p, 'full FK spectral projector') for j, p in enumerate(ps)),
        max(small(p@p-p, 'projector idempotence') for p in ps))
    def polynomial_u(t):
        return np.exp(-1j*baseline*t)*sum(np.exp(-1j*j*t)*p for j, p in enumerate(ps))
    eye3 = np.eye(3, dtype=complex)
    aux = np.kron(np.kron(np.array([1, 1j])/math.sqrt(2), [1, 0]), [1, 0])
    data = np.column_stack([np.kron(aux, np.kron(eye3[:, q], [1, 0, 0])) for q in range(3)])
    initial = np.kron(np.eye(4)[:, 0:1], data)
    uc = sum(np.exp(-1j*j*delta)*p for j, p in enumerate(pcs))
    captured = np.exp(-1j*baseline*delta)*material.gamma@np.kron(uc, np.eye(72))@initial
    def closed_output(s):
        age = delta+s
        phi = np.array([math.sqrt(math.comb(3, j))*math.cos(age/2)**(3-j)*(-1j*math.sin(age/2))**j for j in range(4)])
        phi *= np.exp(-1j*(baseline+1.5)*age)
        return material.gamma@np.kron(phi[:, None], data)
    _, author_states, _, _ = author.group_clock(cone, cylinder)
    mask = np.zeros((4, 8, 3, 3), bool)
    mask[:, [1, 3, 5, 7], :, :] = True
    mask = mask.reshape(288)
    max_state = max_segment = max_effect = max_station = 0.
    for row, author_state in zip(cone, author_states):
        state = closed_output(row['proper_time'])
        max_state = max(max_state, small(state-author_state, 'closed-spin amplitudes versus original full output'))
        cumulative = captured.copy()
        for seg in row['segments']:
            cumulative = polynomial_u(seg['proper_duration'])@cumulative
        max_segment = max(max_segment, small(cumulative-state, 'one capture then every actual free leg, no resets'))
        p = math.sin((delta+row['proper_time'])/2)**6
        max_effect = max(max_effect, small(state[mask].conj().T@state[mask]-p*eye3, 'entire qutrit code clock effect'))
        sta = closed_output(row['coordinate_time'])
        psta = math.sin((delta+row['coordinate_time'])/2)**6
        max_station = max(max_station, small(sta[mask].conj().T@sta[mask]-psta*eye3, 'independent same-material station clock'))
    diagonal = [i for i, row in enumerate(cone) if row['start'] == 0 and row['target'] == 2]
    vp, vm = [closed_output(cone[i]['proper_time']) for i in diagonal]
    f = np.column_stack([np.kron(vp[:, q], vm[:, 0]) for q in (1, 2)]+[np.kron(vp[:, 0], vm[:, q]) for q in (1, 2)])
    ap, am = [author_states[i] for i in diagonal]
    af = np.column_stack([np.kron(ap[:, q], am[:, 0]) for q in (1, 2)]+[np.kron(ap[:, 0], am[:, q]) for q in (1, 2)])
    rail_error = small(f-af, 'all coherent rail columns with both complete D bodies')
    rail_gram = small(f.conj().T@f-np.eye(4), 'all rail/path/payload/reference inputs')
    full_h = baseline*np.eye(288)+h
    ff = f.reshape(288, 288, 4)
    hf = (np.einsum('ab,bcn->acn', full_h, ff)+np.einsum('ab,cbn->can', full_h, ff)).reshape(288**2, 4)
    mass_error = small(f.conj().T@hf-(2*baseline+3)*np.eye(4), 'vac and occupied full mass contributions')

    accuracy = saved['group4_direction_holdout_and_false_positives']['finite_accuracy']
    max_observed = 0.
    lo, hi = accuracy['probability_window']
    for row in cone:
        s = row['proper_time']
        p = math.sin((delta+s)/2)**6
        for a, b, c, d in itertools.product((-1, 1), repeat=4):
            phat = min(hi, max(lo, p+a*accuracy['total_probability_error_contract']))
            shat = 2*math.asin(phat**(1/6))-(delta+b*accuracy['capture_age_error'])
            dhat = (gamma*v+c*accuracy['gamma_v_error'])*shat+d*accuracy['other_endpoint_seam_unit_length_error']
            max_observed = max(max_observed, abs(dhat-row['length']))
    assert max_observed <= accuracy['total_length_error_bound']
    optical_error = 0.
    for entry in saved['group4_direction_holdout_and_false_positives']['single_flat_optical_medium_countercheck']:
        row = next(a for a in cone if a['start'] == entry['start'] and a['target'] == entry['target'])
        z0 = r*np.exp(1j*theta/2)
        z1 = r*np.exp(1j*(theta/2+row['signed_cell_steps']*theta))
        t = np.linspace(0., 1., 4097)
        values = np.sqrt(1-(2/3)*v*v*np.abs(z0+t*(z1-z0)))
        integral = (values[0]+values[-1]+4*sum(values[1:-1:2])+2*sum(values[2:-1:2]))/(3*4096)
        estimate = abs(z1-z0)*integral/v
        optical_error = max(optical_error, small(estimate-entry['flat_medium_moving_proper_time'], 'independent whole-chord Simpson integral'))

    links, issues = [], []
    for name in ('../research_note_1168.md', 'admission.md', 'proof.md', 'NEXT.md'):
        path = HERE/name
        text = path.read_text(encoding='utf-8')
        assert sum(line.strip() == '$$' for line in text.splitlines()) % 2 == 0, name
        assert not any(ord(c) < 32 and c not in '\n\r\t' for c in text), name
        for match in re.finditer(r'\[[^\]]*\]\(([^)]+)\)', text):
            target = match.group(1).strip('<>')
            if re.match(r'^[a-z]+://', target) or target.startswith('#'):
                continue
            dest = (path.parent/unquote(target.split('#')[0])).resolve()
            links.append({'source': name, 'target': target, 'exists': dest.exists()})
            if not dest.exists():
                issues.append([name, target])
    assert not issues, issues
    for name, expected in CORE.items():
        assert digest(HERE/name) == expected, ('core changed during review', name)
    for name, expected in author.SOURCES.items():
        assert digest(BASE/name) == expected, ('source changed during review', name)
    print(json.dumps(dict(status='PASS', kind='independent_review_not_new_scientific_group',
        core_sha256=CORE, source_sha256=author.SOURCES, replay_comparison=stats,
        closest_cylinder_lift_error=max_lift, independently_selected_lift_counts=lift_counts,
        independent_cone_chord_error=max_chord, unfolded_face_interval_error=max_crossing,
        unfolded_initial_direction_error=max_direction, diagonal_angles_radians=angular_gap,
        full_FK_projector_residual=proj_error, closed_spin_full_D_state_error=max_state,
        one_capture_then_segmented_propagation_error=max_segment,
        entire_qutrit_clock_effect_error=max_effect, same_station_clock_effect_error=max_station,
        independent_full_rail_output_error=rail_error, full_rail_reference_gram_error=rail_gram,
        retained_two_positive_body_minimum_mass=2*baseline,
        full_rail_both_body_energy_effect_error=mass_error,
        finite_error_corner_maximum=max_observed, declared_length_error_bound=accuracy['total_length_error_bound'],
        independent_optical_integral_error=optical_error,
        metadata={'local_links_checked': len(links), 'links': links, 'issues': issues},
        frozen_core_and_seven_sources_unchanged=True), ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
