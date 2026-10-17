"""1240 independent check: one continuous event history and a priori bounds.

Uses Gauss--Legendre free-flight integration and direct normal-frame vertices;
does not import the author's Simpson integrator.  Finite numerical histories
diagnose the analytic comparison, not real apparatus or source backreaction.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
TARGET = HERE / 'independent_results.json'
AUTHOR = ['admission.md', 'proof.md', 'check.py', 'results.json', 'NEXT.md']
HISTORY = ['archive_1234_/1236/results.json', 'archive_1234_/1239/results.json',
           'archive_1161_1200/1162/proof.md',
           'archive_1161_1200/1169/proof.md',
           'archive_956_989/research_note_971.md',
           'archive_1219_1233/1231/proof.md']


def need(value, message):
    if not bool(value):
        raise AssertionError(message)


def rat(value):
    return {'fraction': str(value), 'decimal': float(value)}


def exact_bounds_and_vertices():
    r, T, A = Q(1, 10**8), Q(41, 50), Q(3, 2)
    V = 9*r/(1-6*r)
    motion = 4*A*T*V+(A+1)*T*V*V
    need(motion < Q(4429, 10**10), 'published uniform motion upper bound')
    lag, drift = T*V*V, T*V
    old = json.loads((BASE/HISTORY[0]).read_text('utf-8'))
    source = json.loads((BASE/HISTORY[1]).read_text('utf-8'))
    gap = Q(old['group_3']['analytic_Peano_gap_lower']['fraction'])
    D = Q(old['group_2']['D']['fraction'])
    C = Q(old['group_2']['linear_error_coefficient']['fraction'])
    geometry = Q(source['group_1']['per_return_geometry_budget']['fraction'])
    eta = Q(3, 10**7)
    total = geometry+motion+eta
    need(total < Q(1, 10**6), 'conditional combined record budget')
    source_error = 10*total+16*total*total
    flat_error = C*2*total+8*(2*total)**2
    need(source_error < gap and flat_error < -D, 'conditional source/flat separations')
    e = [Q(0), Q(331, 1000), Q(91, 125)]
    reference = [Q(7651, 125000), Q(6327, 15625), Q(101981, 125000)]
    # First-exit contradiction: even using the coarser dust upper geometry,
    # the available conformal time by T exceeds the largest path perturbation.
    conformal_slack = (T-reference[-1])/A
    path_and_trigger = 4*T*V+T*V*V
    need(conformal_slack > path_and_trigger, 'all final receptions precede T')
    gaps = [e[1]-reference[0]-motion, e[2]-reference[1]-motion]
    need(min(gaps) > 0, 'return-before-next-proper-trigger induction')
    # Independent exact local mass-shell checks with nonzero incoming motion.
    checks = []
    for boost in [Q(1), Q(4, 5), Q(9, 10)]:
        mass, q = Q(1), r
        energy = mass*(boost+1/boost)/2
        momentum = mass*(boost-1/boost)/2
        emitted = q*boost
        after_energy, after_p = energy-emitted, momentum-emitted
        need(after_energy**2-after_p**2 == mass**2-2*mass*q, 'emission shell')
        need(after_p < 0 and emitted <= q, 'emitter sign and pulse budget')
        received = q/2
        energy2, p2 = energy+received, momentum-received
        need(energy2**2-p2**2 == mass**2+2*received*(energy+momentum),
             'absorption shell and positive stored energy')
        checks.append({'type': 'emission_absorption', 'boost': str(boost)})
    for boost in [Q(1), Q(5, 4), Q(11, 10)]:
        mass, incident = Q(1), r
        energy = mass*(boost+1/boost)/2
        p = mass*(boost-1/boost)/2
        reflected = incident*(energy-p)/(energy+p+2*incident)
        out_energy, out_p = energy+incident-reflected, p+incident+reflected
        need(out_energy**2-out_p**2 == mass**2, 'finite-mass reflection')
        need(0 < reflected <= incident and out_p > p, 'mirror recedes and pulse loses energy')
        checks.append({'type': 'reflection', 'boost': str(boost)})
    phase_kick = Q(1, 500)
    phase_D_change = e[1]*phase_kick
    need(D+phase_D_change > 0, 'four-momentum does not bound a clock phase kick')
    return dict(status='PASS', finite_ratio=rat(r), velocity_bound=rat(V),
                coordinate_drift_bound=rat(drift), proper_lag_bound=rat(lag),
                motion_bound=rat(motion), inherited_uniform_geometry_bound=rat(geometry),
                unimplemented_remaining_contract=rat(eta), conditional_total=rat(total),
                source_gap_margin=rat(gap-source_error), flat_gap_margin=rat(-D-flat_error),
                first_exit_conformal_slack=rat(conformal_slack),
                first_exit_path_and_trigger_budget=rat(path_and_trigger),
                inter_round_proper_time_gaps=[rat(x) for x in gaps],
                exact_nonzero_velocity_vertex_cases=checks,
                last_vertex_phase_kick=rat(phase_kick), phase_kick_D_change=rat(phase_D_change),
                continuous_phase_is_an_extra_vertex_contract=True,
                conditional_sum_is_not_a_self_consistent_geometric_solution=True)


def determinant(e, rr):
    return float((e[1]-e[0])*(rr[2]-rr[0])-(e[2]-e[0])*(rr[1]-rr[0]))


def history(h, order):
    q, end, distance = 1e-8, .82, .03
    emissions = [0., .331, .728]
    nodes, weights = leggauss(order)
    state = {'t': 0., 'A': {'m': 1., 'p': 0., 'x': 0., 'lag': 0.},
             'B': {'m': 1., 'p': 0., 'x': distance, 'lag': 0.}}
    events, returns, emission_times, reflect_times, receive_times = [], [], [], [], []
    losses, stored = [], []
    max_energy_error, max_momentum_error, max_p, max_comoving_p = 0., 0., 0., 0.
    active_K, expansion_work, complete_rest_change = 0., 0., 0.

    def scale(t):
        return (1+h*np.asarray(t))**(2/3)

    def conformal_interval(t0, t1):
        # Independent expression as the integral of 1/a by Gauss quadrature.
        z = (t0+t1)/2+(t1-t0)*nodes/2
        return float((t1-t0)/2*np.dot(weights, 1/scale(z)))

    def free(body, t0, t1):
        duration = t1-t0
        z = (t0+t1)/2+duration*nodes/2
        a = scale(z)
        conserved = float(scale(t0))*body['p']
        p = conserved/a
        energy = np.sqrt(body['m']**2+p*p)
        v = p/energy
        displacement = float(duration/2*np.dot(weights, v/a))
        # Stable integration of t-tau, retaining effects below an ulp of t.
        delay = float(duration/2*np.dot(weights, v*v/(1+body['m']/energy)))
        return {'m': body['m'], 'p': conserved/float(scale(t1)),
                'x': body['x']+displacement, 'lag': body['lag']+delay}

    def advance(t):
        nonlocal expansion_work
        duration = t-state['t']
        z = (state['t']+t)/2+duration*nodes/2
        a = scale(z)
        energy_pressure = active_K/a
        for label in ['A', 'B']:
            b = state[label]
            p = float(scale(state['t']))*b['p']/a
            energy_pressure = energy_pressure+p*p/np.sqrt(b['m']**2+p*p)
        H = (2/3)*h/(1+h*z)
        expansion_work += float(duration/2*np.dot(weights, H*energy_pressure))
        for label in ['A', 'B']:
            state[label] = free(state[label], state['t'], t)
        state['t'] = t

    def root(function, low, high):
        need(function(low) <= 2e-15 and function(high) >= -2e-15, 'finite event root bracket')
        for _ in range(56):
            middle = (low+high)/2
            if function(middle) < 0:
                low = middle
            else:
                high = middle
        return (low+high)/2

    def vertex(label, kind, photon_energy):
        nonlocal max_energy_error, max_momentum_error, max_p, max_comoving_p, complete_rest_change
        b = state[label]
        m, p = b['m'], b['p']
        energy = np.sqrt(m*m+p*p)
        if kind == 'emit':
            local_energy = photon_energy*(energy+p)/m
            new_m = np.sqrt(m*m-2*m*photon_energy)
            dm = -2*m*photon_energy/(new_m+m)
            dp, de = -local_energy, -local_energy
            outgoing = local_energy
            losses.append(-dm)
        elif kind == 'reflect':
            outgoing = photon_energy*(energy-p)/(energy+p+2*photon_energy)
            dm = 0.
            dp, de = photon_energy+outgoing, photon_energy-outgoing
        else:
            increment_m2 = 2*photon_energy*(energy+p)
            dm = increment_m2/(np.sqrt(m*m+increment_m2)+m)
            dp, de = -photon_energy, photon_energy
            outgoing = 0.
            stored.append(dm)
        new_energy = np.sqrt((m+dm)**2+(p+dp)**2)
        # Stable difference of the complete on-shell energy, retaining recoil.
        actual_de = (dm*(2*m+dm)+dp*(2*p+dp))/(new_energy+energy)
        max_energy_error = max(max_energy_error, abs(actual_de-de))
        max_momentum_error = max(max_momentum_error, abs((p+dp)-p-dp))
        complete_rest_change += dm
        b['m'], b['p'] = m+dm, p+dp
        max_p = max(max_p, abs(b['p']))
        max_comoving_p = max(max_comoving_p, float(scale(state['t']))*abs(b['p']))
        need(b['m'] > 1-6*q-1e-15, 'same finite source not reset or depleted')
        need((state['A']['p'] <= 0 and state['B']['p'] >= 0), 'endpoint kick signs')
        need(outgoing <= q*(1+1e-14), 'finite normal photon budget')
        events.append({'kind': kind, 'endpoint': label, 't': state['t'],
                       'A_tau': state['t']-state['A']['lag'],
                       'A': deepcopy(state['A']), 'B': deepcopy(state['B']),
                       'outgoing_photon_energy': outgoing})
        return outgoing

    for i, trigger in enumerate(emissions):
        current = state['t']
        need(current-state['A']['lag'] <= trigger+2e-15,
             'one history returns before its next proper trigger')
        def trigger_equation(t):
            aa = free(state['A'], current, t)
            return t-aa['lag']-trigger
        te = current if i == 0 else root(trigger_equation, current, end)
        advance(te)
        emission_times.append(te)
        outgoing = vertex('A', 'emit', q)
        photon_x, photon_t = state['A']['x'], te
        photon_comoving_energy = outgoing*float(scale(te))
        active_K = photon_comoving_energy
        current = state['t']
        def mirror_equation(t):
            bb = free(state['B'], current, t)
            return photon_x+conformal_interval(photon_t, t)-bb['x']
        tm = root(mirror_equation, current, end)
        advance(tm)
        reflect_times.append(tm)
        outgoing = vertex('B', 'reflect', photon_comoving_energy/float(scale(tm)))
        photon_x, photon_t = state['B']['x'], tm
        photon_comoving_energy = outgoing*float(scale(tm))
        active_K = photon_comoving_energy
        current = state['t']
        def reception_equation(t):
            aa = free(state['A'], current, t)
            return aa['x']-(photon_x-conformal_interval(photon_t, t))
        tr = root(reception_equation, current, end)
        advance(tr)
        receive_times.append(tr)
        vertex('A', 'absorb', photon_comoving_energy/float(scale(tr)))
        active_K = 0.
        returns.append(tr-state['A']['lag'])

    e = np.array(emissions)
    u = (1+h*e)**(1/3)
    nominal = e+.06*u*u+.0012*h*u+.000008*h*h
    delta = np.array(returns)-nominal
    V = 9e-8/(1-6e-8)
    bound = 4.92*V+2.05*V*V
    need(np.max(abs(delta)) < bound and np.max(abs(delta)) > 1e-10,
         'one finite pulse history has nonzero shifts inside analytic budget')
    need(max_energy_error < 1e-22 and max_momentum_error < 1e-22,
         'normal-frame complete four-momentum at every vertex')
    need(max_comoving_p <= 9*q and max_p <= 6*q, 'full cumulative momentum')
    need(sum(losses) <= 6*q and min(stored) > 0, 'fuel spent and absorption stored')
    need(max(receive_times) < end and returns[0] < emissions[1]
         and returns[1] < emissions[2], 'strict nonoverlap and finite domain')
    need(state['A']['lag'] <= end*V*V, 'continuous proper-time lag')
    complete_energy_change = complete_rest_change
    for label in ['A', 'B']:
        b = state[label]
        complete_energy_change += b['p']**2/(np.sqrt(b['m']**2+b['p']**2)+b['m'])
    pressure_residual = complete_energy_change+expansion_work
    need(expansion_work > 0 and abs(pressure_residual) < 2e-21,
         'same branch total energy and redshift work')
    signed = determinant(e, returns)
    kicked = list(returns)
    kicked[-1] += .002
    kicked_D = determinant(e, kicked)
    need(abs((kicked_D-signed)-.331*.002) < 2e-16 and kicked_D > 0,
         'last-event phase deletion preserves mechanical history but defeats bound')
    return dict(status='PASS', h=float(h), quadrature_order=order,
                initial_mass_A=1., initial_mass_B=1., same_rest_pulse_energy=q,
                same_proper_trigger_menu=emissions, no_resets=True,
                emission_coordinate_times=emission_times, mirror_coordinate_times=reflect_times,
                reception_coordinate_times=receive_times, reception_proper_readings=returns,
                original_readings=nominal.tolist(), motion_shifts=delta.tolist(),
                source_fuel_spent=float(sum(losses)), received_rest_energy=float(sum(stored)),
                complete_branch_energy_change=float(complete_energy_change),
                cosmological_redshift_work=float(expansion_work),
                full_energy_work_residual=float(pressure_residual),
                total_comoving_momentum_bound=float(max_comoving_p),
                normal_vertex_energy_residual=float(max_energy_error),
                normal_vertex_momentum_residual=float(max_momentum_error),
                D=signed, phase_kicked_D=kicked_D, final_state=deepcopy(state), events=events,
                actual_clock_reading_vertex_permission=False,
                geometry_backreaction_of_endpoints_and_pulses_certified=False)


def independent_histories():
    rows, refinement = [], []
    for h in [1., np.sqrt(.8)]:
        coarse = history(h, 12)
        fine = history(h, 20)
        error = float(np.max(abs(np.array(coarse['reception_proper_readings'])-
                                 np.array(fine['reception_proper_readings']))))
        need(error < 2e-13, 'different quadrature order agrees')
        rows.append(fine)
        refinement.append(error)
    need(rows[0]['same_proper_trigger_menu'] == rows[1]['same_proper_trigger_menu']
         and rows[0]['same_rest_pulse_energy'] == rows[1]['same_rest_pulse_energy'],
         'source comparison really shares one menu and resource rule')
    return dict(status='PASS', method='Gauss--Legendre free legs; direct normal-frame vertices',
                cases=rows, quadrature_differences=refinement,
                source_D_gap=rows[1]['D']-rows[0]['D'],
                finite_numerics_do_not_prove_instrument_or_backreaction_contract=True)


def digest(paths):
    return {str(p.relative_to(BASE)).replace('\\', '/'):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def run(bind=True):
    watched = [BASE/p for p in HISTORY]
    if bind:
        watched += [HERE/p for p in AUTHOR]
    before = digest(watched)
    result = dict(round=1240, status='PASS', independent_scientific_groups=2,
                  exact_budget_and_local_vertices=exact_bounds_and_vertices(),
                  continuous_history=independent_histories(),
                  CT='candidate conjecture, not adopted',
                  full_GR_delivery_complete=False, actual_vertex_implementation=False,
                  complete_apparatus_covariant_source=False,
                  unknown_quantum_clock_reference_channel=False,
                  source_sha256=before)
    need(before == digest(watched), 'author and historical files remain read-only')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--diagnostic', action='store_true')
    args = parser.parse_args()
    need(not (args.write and args.diagnostic), 'choose one mode')
    result = run(bind=not args.diagnostic)
    if args.diagnostic:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    elif args.write:
        with TARGET.open('x', encoding='utf-8') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write('\n')
        print('WROTE_ONCE: 1240 independent_results.json; two groups PASS')
    else:
        need(result == json.loads(TARGET.read_text('utf-8')), 'READ_ONLY_MATCH')
        print('READ_ONLY_MATCH: 1240 independent; two groups PASS')
