"""1219 candidate: thermal rest energy, actual recoil, records and the old clock.

Finite effective K1219 permissions are inputs, not microscopic generation.
Default recomputes read-only; --write-once creates results.json exactly once.
Old scientific suites are not rerun. Only round1111.model() is reused.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEPS = {
    'archive_1106_1119/1111/check.py': '704ae4ca016c77e0df83fc9c2d39e90d79e9fbf3d6469c8c5b9bf1dad8a5bd44',
    'archive_1106_1119/1118/proof.md': 'cbf928ac8a1dc8997869eddc50dae4774c58e9baecdfad904ad813a03ae3cd60',
    'archive_1106_1119/1119/proof.md': 'b4a97154a8302a49c4f4f6788f615909b4849820bfa42e3bf11574c16808008d',
    'archive_1106_1119/_shared/1116_1119/adoption_and_dependencies.md': '1c1395c16cf2392be0c6e7275d660f454113290bd0e23d6892e08b3aa41b4a0e',
    'archive_1159_1160/1159/proof.md': '20abbc0de0feeafc4ff51163e646beeaa91ead8596e6264d4a274177d04aaa19',
    'archive_1159_1160/1160/check.py': '8af83d78311b77c804693e055f2ece2c9f195e1d217b5b67e41fa6b3b6753170',
    'archive_1201_1218/1205/proof.md': '192e0a6256b618fccfbc1c1754056c1017481415bbd9eafce532f34f9f681872',
}
TOL = 3e-10


def small(x, label):
    value = float(np.linalg.norm(x))
    assert value < TOL, (label, value)
    return value


def load_old():
    for rel, digest in DEPS.items():
        assert hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() == digest, rel
    spec = importlib.util.spec_from_file_location('round1111_model_only_1219', ROOT / 'archive_1106_1119/1111/check.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def carrier():
    # Fixed, distinguishable carrier slots; local order n=0,1 then r=+,-.
    p = np.array([[5., 4., 0., 0.], [5., -4., 0., 0.],
                  [20/3, 16/3, 0., 0.], [20/3, -16/3, 0., 0.]])
    h = np.array([0., 0., 1., 1.])
    mass = 3 + h
    return p, h, mass


def core():
    p, h, mass = carrier()
    swap = np.zeros((16, 16), complex)
    for a, b in np.ndindex(4, 4):
        swap[4*b+a, 4*a+b] = 1
    half = (np.eye(16)-1j*swap)/math.sqrt(2)
    # Active core order AB, M, E, J. All records have positive mass1,
    # constant across values. Program L is a separate retained mass1 body.
    record = np.zeros((128, 128), complex)
    pointer = np.zeros_like(record)
    for ab, rm, re, j in np.ndindex(16, 2, 2, 2):
        col = 8*ab+4*rm+2*re+j
        a = ab//4
        direction = a % 2
        record[8*ab+4*(rm ^ direction)+2*(re ^ direction)+j, col] = 1
        pointer[8*ab+4*rm+2*re+(j ^ (a//2)), col] = 1
    collision = np.kron(half, np.eye(8))
    charges = np.array([p[a]+p[b] for a, b in np.ndindex(4, 4)])
    pa = np.array([p[a] for a, b in np.ndindex(4, 4)])
    pb = np.array([p[b] for a, b in np.ndindex(4, 4)])
    ha = np.array([h[a] for a, b in np.ndindex(4, 4)])
    hb = np.array([h[b] for a, b in np.ndindex(4, 4)])
    return locals()


def initial_core(na, nb):
    v = np.zeros(128, complex)
    a, b = 2*na, 2*nb+1  # A incoming +, B incoming -.
    v[8*(4*a+b)] = 1
    return v


def execute(c, v, q, program, leak=False, terminal=True):
    # Conservative gate core only. Every physical window also includes
    # exp(-i delta P_all^0). All gates commute with that full free factor.
    # On the certified retained-source input the added energy is 92/3,
    # hence its omission here is one global phase. On arbitrary mother
    # inputs the full factor MUST be retained; diagonal event probabilities
    # below are invariant, but arbitrary full state amplitudes need not be.
    first = c['collision'] @ v if q == 1 else v.copy()
    rec = c['record'] @ first if program in (1, 2) else first.copy()
    # Leak means a distinct extra retained environment, not deleting a factor
    # and pretending to have the same physical experiment.
    if leak:
        parts = []
        for ev in (0, 1):
            branch = rec.reshape(16, 2, 2, 2).copy()
            branch[:, :, 1-ev, :] = 0
            branch = branch.reshape(128)
            if program == 2:
                branch = c['record'] @ branch
            branch = c['collision'] @ branch
            if terminal:
                branch = c['pointer'] @ branch
            parts.append(branch)
        return np.stack(parts, axis=1), first, rec
    if program == 2:
        rec = c['record'] @ rec
    out = c['collision'] @ rec
    if terminal:
        out = c['pointer'] @ out
    return out[:, None], first, rec


def event_statistics(v):
    # v may have many orthogonal kept-source/environment columns.
    weights = np.sum(np.abs(v)**2, axis=tuple(range(1, v.ndim))).reshape(16, 2, 2, 2)
    joint = np.zeros((2, 2))  # pointer J then true A outgoing direction.
    for ab, rm, re, j in np.ndindex(16, 2, 2, 2):
        joint[j, (ab//4) % 2] += weights[ab, rm, re, j]
    ab = weights.sum(axis=(1, 2, 3))
    return joint, ab


def thermal_outputs(c, qa, qb, q, program, leak=False):
    out = []
    for na, nb in np.ndindex(2, 2):
        weight = (qa if na else 1-qa)*(qb if nb else 1-qb)
        result, _, _ = execute(c, initial_core(na, nb), q, program, leak)
        out.append(math.sqrt(weight)*result)
    return np.concatenate(out, axis=1)


def prep_certificate(nseed, selected, direction):
    # Whole A/stock 4x4 mother, seed nseed. Act only if the two incoming
    # velocity labels agree, so the mother extension conserves fourP too.
    p, _, _ = carrier()
    d = 16*nseed
    permutation = np.zeros((d, d), complex)
    charges = np.zeros((d, 4))
    for a, stock, seed in np.ndindex(4, 4, nseed):
        col = (4*a+stock)*nseed+seed
        oa, os = a, stock
        if seed < selected and a % 2 == stock % 2:
            oa, os = 2*(stock//2)+a % 2, 2*(a//2)+stock % 2
        row = (4*oa+os)*nseed+seed
        permutation[row, col] = 1
        # Seed is a same-velocity mass1 material.
        useed = np.array([5/3, (1-2*direction)*4/3, 0, 0])
        charges[col] = p[a]+p[stock]+useed
    unit = small(permutation.conj().T@permutation-np.eye(d), 'prep full mother unitary')
    cons = max(small((charges[:, mu, None]-charges[None, :, mu])*permutation, 'prep fourP') for mu in range(4))
    before = np.zeros(d, complex)
    for seed in range(nseed):
        before[(4*direction+(2+direction))*nseed+seed] = 1/math.sqrt(nseed)
    after = permutation@before
    state = after.reshape(4, 4, nseed)
    rhoa = state.reshape(4, -1)@state.reshape(4, -1).conj().T
    target = np.zeros(4)
    target[direction], target[2+direction] = 1-selected/nseed, selected/nseed
    residual = small(rhoa-np.diag(target), 'exact finite-count Gibbs marginal')
    # Store+seed not discarded; source remains a pure full vector.
    fullcharge = np.sum(np.abs(after[:, None])**2*charges, axis=0)
    expected = np.array([40/3, (1-2*direction)*32/3, 0, 0])
    small(fullcharge-expected, 'complete finite source fourP')
    return dict(seed_dimension=nseed, selected_labels=selected, target_excited=selected/nseed,
                whole_source_energy=fullcharge[0], whole_source_Px=fullcharge[1],
                seed_purity_resource='known pure uniform superposition, retained',
                target_stock_initial_occupations=[0, 1],
                mother_unitarity=unit, mother_fourP=cons, Gibbs_marginal_residual=residual)


def calculate():
    old = load_old()
    m = old.model()
    c = core()
    source_a = prep_certificate(4, 1, 0)
    source_b = prep_certificate(5, 2, 1)
    reverse_a = prep_certificate(5, 2, 0)
    reverse_b = prep_certificate(4, 1, 1)
    p, h, mass = carrier()
    id128 = np.eye(128)
    unit = max(small(g.conj().T@g-id128, 'primitive full unitary')
               for g in (c['collision'], c['record'], c['pointer']))
    cons = 0.0
    present = 0.0
    for mu in range(4):
        dtotal = np.repeat(c['charges'][:, mu], 8)
        da = np.repeat(c['pa'][:, mu], 8)
        for gate in (c['collision'], c['record'], c['pointer']):
            cons = max(cons, small((dtotal[:, None]-dtotal[None, :])*gate, 'fourP full mother'))
        for gate in (c['record'], c['pointer']):
            present = max(present, small((da[:, None]-da[None, :])*gate, 'record current A fourP'))
    int_total = np.repeat(c['ha']+c['hb'], 8)
    intcons = small((int_total[:, None]-int_total[None, :])*c['collision'], 'rest internal exchange')
    beta_a, beta_b = math.log(3), math.log(1.5)
    small(1/(1+math.exp(beta_a))-0.25, 'betaA')
    small(1/(1+math.exp(beta_b))-0.4, 'betaB')

    rows = []
    outputs = {}
    boost = old.boost(np.array([0.6, 0., 0.]))
    for label, qa, qb in [('cold_A_hot_B', .25, .4), ('hot_A_cold_B', .4, .25)]:
        for q in (0, 1):
            for pr in (0, 1, 2):
                v = thermal_outputs(c, qa, qb, q, pr)
                outputs[label, q, pr] = v
                joint, probs = event_statistics(v)
                lam = .5 if q == 0 or pr == 1 else 1.
                expected_joint = np.array([[(1-lam)*(1-qa), lam*(1-qb)],
                                           [(1-lam)*qa, lam*qb]])
                small(joint-expected_joint, 'heat pointer and original direction joint')
                delta_pa = probs@c['pa'] - (1-qa)*p[0]-qa*p[2]
                delta_pb = probs@c['pb'] - (1-qb)*p[1]-qb*p[3]
                small(delta_pa+delta_pb, 'actual carrier recoil')
                dq = float(probs@c['ha']-qa)
                small(dq-lam*(qb-qa), 'thermal rest exchange')
                small(delta_pa-lam*np.array([(5/3)*(qb-qa), -(4/3)*(6+qa+qb), 0, 0]), 'joint thermal mechanical exchange')
                small(joint[1].sum()-(qa+dq), 'terminal actual energy pointer')
                rows.append(dict(preparation=label, unknown_control_basis=q, program=pr,
                                 exchange_weight=lam, joint_J_direction=joint.tolist(),
                                 J_excited=float(joint[1].sum()), A_reflected=float(joint[:, 1].sum()),
                                 A_delta_rest_internal=dq, A_delta_fourP=delta_pa.tolist(),
                                 B_delta_fourP=delta_pb.tolist(), A_boost_delta_fourP=(boost@delta_pa).tolist()))

    # Sources stay as separate orthogonal seed+stock records. No postselection.
    source_outputs = {}
    norm_res = gram_res = sampled_ref_res = 0.0
    rng = np.random.default_rng(1219)
    for pr in (0, 1, 2):
        phi = []
        for q in (0, 1):
            columns = []
            for sa, sb in np.ndindex(4, 5):
                na, nb = int(sa < 1), int(sb < 2)
                col, _, _ = execute(c, initial_core(na, nb), q, pr)
                columns.append(col[:, 0]/math.sqrt(20))
            phi.append(np.column_stack(columns))
        source_outputs[pr] = phi
        norm_res = max(norm_res, *(small(np.vdot(f, f)-1, 'retained full source vector') for f in phi))
        for q in (0, 1):
            j, _ = event_statistics(phi[q])
            jr, _ = event_statistics(outputs['cold_A_hot_B', q, pr])
            small(j-jr, 'actual retained-source probabilities equal marginal calculation')
        # The retained source Q output is orthogonal; the old Gamma is an
        # isometry on that full source, so this is the exact logical Gram.
        gram = np.diag([np.vdot(f, f) for f in phi])
        gram_res = max(gram_res, small(gram-np.eye(2), 'all unknown Q and arbitrary R'))
        psi = rng.normal(size=(2, 3))+1j*rng.normal(size=(2, 3))
        psi /= np.linalg.norm(psi)
        out = np.stack([phi[q][:, :, None]*psi[q][None, None, :] for q in (0, 1)])
        rout = np.einsum('qijr,qijs->rs', out.conj(), out)
        rin = psi.conj().T@psi
        sampled_ref_res = max(sampled_ref_res, small(rout-rin, 'entangled R reduced state unchanged'))

    recovery = max(small(source_outputs[0][q]-source_outputs[2][q], 'complete record inverse') for q in (0, 1))
    leaked = thermal_outputs(c, .25, .4, 1, 2, leak=True)
    joint_leak, _ = event_statistics(leaked)
    joint_kept, _ = event_statistics(outputs['cold_A_hot_B', 1, 2])
    small(joint_leak-np.array([[.375, .3], [.125, .2]]), 'kept inaccessible history changes actual future')
    future_gap = float(joint_kept[1, 1]-joint_leak[1, 1])
    assert future_gap > .19999999

    # Nonthermal coherent A has identical h populations. This is a separate
    # mathematical input diagnostic, NOT a claimed stock preparation/readout.
    coherent = []
    for nb in (0, 1):
        v = math.sqrt(.75)*initial_core(0, nb)+.5*initial_core(1, nb)
        vv, _, _ = execute(c, v, 1, 0)
        coherent.append(math.sqrt(.4 if nb else .6)*vv)
    coherent = np.concatenate(coherent, axis=1)
    coherent_joint, _ = event_statistics(coherent)
    blind = small(coherent_joint-joint_kept, 'population readout does not establish Gibbs qualification')
    sigma_a = np.array([[.75, math.sqrt(3)/4], [math.sqrt(3)/4, .25]])
    nonthermal_stationarity = float(np.linalg.norm(np.diag([0, 1])@sigma_a-sigma_a@np.diag([0, 1])))
    assert nonthermal_stationarity > .6

    # Old D is used unchanged, with 7 positive boundary windows ending at pi:
    # two finite source preparations, four collision/record slots, final h copy.
    hd, gamma, hc = m['hout'], m['gamma'], m['hc']
    intertwiner = small(hd@gamma-gamma@(m['baseline']*np.eye(288)+np.kron(hc, np.eye(72))), 'same old full H')
    t0, total, delta = math.pi/2, math.pi, math.pi/14
    ev, vec = np.linalg.eigh(hc)
    def head_wait(t):
        return (vec*np.exp(-1j*ev*t))@vec.conj().T
    head_word = np.linalg.matrix_power(head_wait(delta), 7)@head_wait(t0)
    head_res = small(head_word-head_wait(total), 'old clock full seven-window word')
    def old_state(t):
        decoded = m['ready'].reshape(4, 72, 2)
        evolved = np.einsum('ab,bkq->akq', head_wait(t), decoded).reshape(288, 2)
        return gamma@evolved*np.exp(-1j*m['baseline']*t)*m['phase']
    start, final = old_state(t0), old_state(total)
    # Independent old time function called only once, not calculate().
    direct_old = small(final-m['w'](total), 'unchanged old clock exact function')
    oldgram = small(final.conj().T@final-np.eye(2), 'physical full D isometry')
    def mr_prob(v):
        state = v.reshape(4, 8, 3, 3, 2)
        return np.sum(abs(state[:, [1, 3, 5, 7], ...])**2, axis=(0, 1, 2, 3))
    old_mr_start, old_mr_final = mr_prob(start), mr_prob(final)
    small(old_mr_start-np.ones(2)/8, 'original Mr at pi/2')
    small(old_mr_final-np.ones(2), 'original Mr at pi')
    # Physical dressed gates commute with H_D because their decoded control
    # acts only on source Qc. Check its full old D representative once.
    cc = np.kron(np.eye(32), np.kron(np.diag([0, 0, 1]), np.eye(3)))
    dressed_c = gamma@cc@gamma.conj().T
    control_comm = small(hd@dressed_c-dressed_c@hd, 'physical full D control not bare receiver')

    # Actual endpoints, not fitted coordinates: each whole carrier has
    # sigma=6/5, u=P/M, including both rest-energy values.
    finish = np.array([m['arrival']+total, 1., 0., 0.])
    begin = finish-np.array([7*delta, 0., 0., 0.])
    endpoints = np.array([finish+(6/5)*p[i]/mass[i] for i in range(4)])
    cov = shell = phase_res = 0.0
    for i in range(4):
        pp = boost@p[i]
        dz = endpoints[i]-finish
        cov = max(cov, small(boost@dz-(6/5)*pp/mass[i], 'same full endpoint transport'))
        shell = max(shell, small(old.mink(pp, pp)-mass[i]**2, 'same internal-energy mass shell'))
        phase_res = max(phase_res, small(old.mink(p[i], dz)-(6/5)*mass[i], 'proper-time full energy phase'))
    endpoint_gap = float(np.linalg.norm(endpoints[0, 1:]-endpoints[1, 1:]))
    small(endpoint_gap-16/5, 'real endpoint separation')
    # Incoming positions depend on the common velocity, not the unknown h.
    incoming_a = begin[1:]-np.array([.8, 0., 0.])*begin[0]
    incoming_b = begin[1:]+np.array([.8, 0., 0.])*begin[0]
    small(incoming_a+.8*begin[0]*np.array([1, 0., 0.])-begin[1:], 'same input A arrival')
    small(incoming_b-.8*begin[0]*np.array([1, 0., 0.])-begin[1:], 'same input B arrival')

    # Delete energy from mass while keeping the actual h clock: scalar-only
    # P then predicts no relative phase for a sigma=pi energy-superposition.
    wrong_phase_gap = abs(np.exp(-1j*math.pi)-1)/2
    wrong_shell = p[2, 0]**2-p[2, 1]**2-9  # Actual excited mass is 4, not3.
    assert abs(wrong_shell-7) < TOL and wrong_phase_gap > .99999999
    # A heat-only swap with unchanged velocity has a full-mother momentum
    # defect; for initial |1,+;0,-> the full-swap defect is -8/3.
    heatonly = np.zeros((16, 16), complex)
    for a, b in np.ndindex(4, 4):
        oa, ob = 2*(b//2)+a % 2, 2*(a//2)+b % 2
        heatonly[4*oa+ob, 4*a+b] = 1
    mom = c['charges'][:, 1]
    heatonly_norm = float(np.linalg.norm((mom[:, None]-mom[None, :])*heatonly, 2))
    assert abs(heatonly_norm-8/3) < TOL
    first_col = 4*2+1
    after_probs = abs(heatonly[:, first_col])**2
    witness_defect = float(after_probs@mom-mom[first_col])
    small(witness_defect+8/3, 'heat-only wrong recoil witness')
    heatonly_thermal_defect = 0.0
    for na, nb in np.ndindex(2, 2):
        col = 4*(2*na)+(2*nb+1)
        weight = (.25 if na else .75)*(.4 if nb else .6)
        heatonly_thermal_defect += weight*(abs(heatonly[:, col])**2@mom-mom[col])
    small(heatonly_thermal_defect-.4, 'same thermal input heat-only net momentum failure')

    qa, qb = F(1, 4), F(2, 5)
    dq = qb-qa
    de, dpx = F(5, 3)*dq, -F(4, 3)*(6+qa+qb)
    deboost = F(5, 4)*(de-F(3, 5)*dpx)
    exact = dict(qA=str(qa), qB=str(qb), delta_hA=str(dq), delta_EA=str(de),
                 delta_PxA=str(dpx), delta_boost_EA=str(deboost),
                 recorded_delta_hA=str(dq/2), recorded_delta_EA=str(de/2),
                 joint_J1_reflected_full=str(qb), joint_J1_reflected_recorded=str(qb/2),
                 joint_J1_forward_recorded=str(qa/2), J1_recorded=str((qa+qb)/2),
                 endpoint_gap=str(F(16, 5)), heatonly_momentum_defect=str(-F(8, 3)),
                 added_all_resource_energy=str(F(92, 3)))
    return dict(schema='round1219_thermal_recoil_candidate_v1', status='PASS', check_groups=5,
        adopted_axioms=0, empirical_calibrations=0, source_sha256=DEPS,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        exact_witnesses=exact,
        group1_same_thermal_preparation=dict(A=source_a, B=source_b, reverse_A=reverse_a, reverse_B=reverse_b,
            betaA=beta_a, betaB=beta_b,
            initial_global_source_energy=80/3, positive_stationary_record_program_mass=4,
            old_complete_D_energy_upper=m['baseline']+3,
            full_initial_support_energy_upper=m['baseline']+3+92/3,
            full_mother_dimension=288*4**4*4*5*2**3*3,
            actual_retained_source_labels=20, no_postselection=True,
            local_Gibbs_only=True, initial_source_correlations_retained=True),
        group2_complete_mechanical_gates=dict(unitarity=unit, fourP=cons, rest_internal_sum=intcons,
            records_preserve_current_A_fourP=present, maximum_active_core_dimension=128),
        group3_joint_predictions_and_future=dict(rows=rows, coherent_record_recovery=recovery,
            leaked_record_joint=joint_leak.tolist(), actual_joint_future_gap=future_gap,
            leakage_contrast_added_positive_record_mass=1,
            leakage_contrast_requires_separate_same_window_copy_permission=True,
            nonthermal_population_blindness=blind, nonthermal_stationarity_defect=nonthermal_stationarity,
            nonthermal_is_separate_mathematical_diagnostic_not_claimed_seed_preparation=True),
        group4_original_clock_and_unknown_reference=dict(original_FK_intertwiner=intertwiner,
            dressed_source_control_commutator=control_comm, seven_window_clock_word=head_res,
            independent_old_final_state=direct_old, old_complete_isometry=oldgram,
            old_Mr_start=old_mr_start.tolist(), old_Mr_finish=old_mr_final.tolist(),
            retained_source_normalization=norm_res, all_input_Gram=gram_res,
            entangled_reference_sample=sampled_ref_res, windows=7, positive_window=delta,
            original_elapsed_start=t0, original_elapsed_finish=total),
        group5_actual_ports_and_deletions=dict(begin_event=begin.tolist(), end_event=finish.tolist(),
            A_initial_position=incoming_a.tolist(), B_initial_position=incoming_b.tolist(),
            outgoing_endpoints_by_level_direction=endpoints.tolist(), endpoint_gap=endpoint_gap,
            transported_endpoint_residual=cov, mass_shell_residual=shell, proper_phase_residual=phase_res,
            deleting_internal_mass_contribution_shell_defect=float(wrong_shell),
            deleting_internal_mass_clock_halftrace_witness=float(wrong_phase_gap),
            heatonly_swap_momentum_defect_norm=heatonly_norm,
            heatonly_swap_specific_momentum_defect=witness_defect,
            heatonly_swap_same_thermal_input_mean_momentum_defect=float(heatonly_thermal_defect)),
        limitations=['K1219 finite positive-window mechanical permissions are new model inputs',
            'The old clock does not autonomously trigger these effective collisions',
            'Known degenerate seeds, excited stocks and blank records are input resources',
            'Rest-internal Gibbs is not global moving-carrier or whole-process Gibbs',
            'Access to a retained stock/seed changes reference and correlation qualifications',
            'Actual finite incoming preparations are narrower than the full algebraic mother',
            'No accelerated detector, Unruh response, local stress, area entropy or GR is proved'])


def compare(a, b, path='root'):
    if isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            compare(a[k], b[k], path+'.'+k)
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            compare(x, y, path+f'[{i}]')
    elif isinstance(a, float):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-10), (path, a, b)
    else:
        assert a == b, (path, a, b)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-once', action='store_true')
    args = parser.parse_args()
    result = calculate()
    output = HERE/'results.json'
    if args.write_once:
        with output.open('x', encoding='utf-8') as fp:
            json.dump(result, fp, ensure_ascii=False, indent=2)
            fp.write('\n')
        print('PASS: 5 groups; first candidate results written')
    else:
        compare(result, json.loads(output.read_text(encoding='utf-8')))
        print('PASS: 5 groups; READ_ONLY_MATCH')
    print('Finite K1219 candidate; no new axiom or geometry/GR conclusion')


if __name__ == '__main__':
    main()
