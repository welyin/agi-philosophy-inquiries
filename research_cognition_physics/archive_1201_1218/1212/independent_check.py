"""Independent 1212 audit: frozen-source intervals and complete work balance.

No import or execution of the author's calculation.  Finite matrix tests are
diagnostics; the entropy identity and interval inequalities carry the proof.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def read(p): return json.loads(p.read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def entropy(r):
    vals = np.linalg.eigvalsh((r+r.conj().T)/2)
    assert vals.min() > -1e-12
    vals = vals[vals > 1e-15]
    return float(-np.dot(vals, np.log(vals)))


def marginal(r, dims, keep):
    remaining = list(range(len(dims)))
    out = r.reshape(tuple(dims)*2)
    for i in reversed(range(len(dims))):
        if i not in keep:
            j = remaining.index(i)
            out = np.trace(out, axis1=j, axis2=j+len(remaining))
            remaining.pop(j)
    size = int(np.prod([dims[i] for i in remaining]))
    return out.reshape(size, size)


def tensor(items):
    out = np.ones((1, 1))
    for a in items: out = np.kron(out, a)
    return out


def embed(a, i, dims):
    return tensor([a if j == i else np.eye(d) for j, d in enumerate(dims)])


def lower_square(lo, hi):
    assert lo <= hi
    return F(0) if lo <= 0 <= hi else min(lo*lo, hi*hi)


def evaluate():
    a2, a4, a5 = [read(STAGE/str(n)/'results.json') for n in (1202, 1204, 1205)]
    den = a2['preparation']['root_denominator']
    roots = [(F(a, den), F(b, den)) for a, b in a2['preparation']['root_boxes']]
    width = F('1.094')
    assert roots[-1][1]-roots[0][0] < width
    populations = {'cold': list(map(F, a2['one_marginal_rational'])),
                   'hot': list(map(F, a4['new_reference']['one_marginal_rational']))}
    offsets = {}
    for name in populations:
        tc = a4['temperature_comparisons'][name]
        tau = list(map(F, tc['gibbs_population_rational']))
        err = F(tc['tau_entry_error_upper'])
        lo, hi = F(0), F(0)
        for p, t, (elo, ehi) in zip(populations[name], tau, roots):
            corners = [(p-u)*e for u in (t-err, t+err) for e in (elo, ehi)]
            lo += min(corners); hi += max(corners)
        offsets[name] = (lo, hi)
    assert offsets['hot'][1] < 0
    din = sum(F(a4['temperature_comparisons'][s]['relative_entropy_chi_square_upper'])
              for s in populations)
    assert din == F('.000502765845')
    eta = F(9, 10**15)
    initial_D_extra = F(2, 10**12)
    assert F(a5['normal_clock']['conditional_joint_trace_error_upper']) < eta
    bounds = {}
    for t in (300, 600):
        key = f't{t}_lambda1'
        old = a4['actual_flow_predictions'][key]
        current = a5['native_flow_transport'][key]
        radius = F(current['new_actual_energy_radius_upper'])
        qc = F(old['delta_EA_center']); qh = -F(old['delta_EB_center'])
        hot = lower_square(offsets['hot'][0]-qh-radius-width*eta,
                           offsets['hot'][1]-qh+radius+width*eta)*2/width**2
        cold = lower_square(offsets['cold'][0]+qc-radius-width*eta,
                            offsets['cold'][1]+qc+radius+width*eta)*2/width**2
        # These endpoint bounds use certified energy transport, not a claim
        # that the correlated preparation has the initial free-resource din.
        hot_gap = (hot-din-initial_D_extra)/F('1.2')
        assert hot_gap > F('.0008') if t == 300 else hot_gap > F('.007')
        p0 = F(old['A_pointer_probabilities_center'][0])
        pointer_error = F(current['pointer_probability_error_upper'])
        cond_not_blank = 1-p0-pointer_error
        all_attempt_not_blank = F(a5['fixed_program']['worst_joint_success_probability_lower'])*cond_not_blank
        assert all_attempt_not_blank > (F('.625') if t == 300 else F('.653'))
        p_success = F(a5['fixed_program']['worst_joint_success_probability_lower'])
        energy_radius = 2*F(1094)*eta + 2*width*F(current['rectangular_shape_error'])
        energy_radius += F(1, 10**10)+2*F(a5['normal_clock']['endpoint_total_interaction_upper'])/p_success
        clock_center = -F(old['bare_work_center'])
        s_bound = F(a5['normal_clock']['conditional_clock_entropy_upper'])
        clock_F_interval = (clock_center-energy_radius-s_bound/F('1.2'),
                            clock_center+energy_radius+s_bound/F('1.2'))
        if t == 300: assert clock_F_interval[1] < F('-.000033')
        else: assert clock_F_interval[0] < 0 < clock_F_interval[1]
        bounds[str(t)] = dict(
            q_hot_outer=[str(qh-radius), str(qh+radius)],
            q_cold_outer=[str(qc-radius), str(qc+radius)],
            hot_D_lower_exact=str(hot), two_reference_D_lower_exact=str(hot+cold),
            independent_reference_hypothetical_work_deficit_lower_exact=str(hot_gap),
            actual_factory_initial_correlations_still_required=True,
            actual_conditional_nonblank_lower_exact=str(cond_not_blank),
            actual_all_attempt_nonblank_lower_exact=str(all_attempt_not_blank),
            full_material_clock_work_value_interval=[str(x) for x in clock_F_interval])

    claimed = read(HERE/'results.json')
    for row in claimed['finite_flow_bounds']:
        if row['family'] != '1205_conserved_S_contact': continue
        actual = bounds[str(row['t'])]
        assert F(row['hot_relative_entropy_lower']) <= F(actual['hot_D_lower_exact'])
        assert F(row['hot_only_work_value_gap_lower']) <= F(actual['independent_reference_hypothetical_work_deficit_lower_exact'])
    for row in claimed['actual_factory_cycle_obstruction']:
        actual = bounds[str(row['t'])]
        assert F(row['raw_K_nonblank_lower']) == F(actual['actual_all_attempt_nonblank_lower_exact'])
        assert list(map(F, row['contact_clock_work_value_change_interval'])) == list(map(F, actual['full_material_clock_work_value_interval']))

    # Reconstruct original h, QQ and their finite spectral evolution afresh.
    h = np.array([[0., -.2, 0., -np.sqrt(3)/80], [-.2, 1., 0., 0.],
                  [0., 0., 1., 0.], [-np.sqrt(3)/80, 0., 0., -.025]])
    q = np.zeros((4, 4)); q[1, 2] = q[2, 1] = 1
    e, eigenbasis = np.linalg.eigh(h)
    H = np.kron(h, np.eye(4))+np.kron(np.eye(4), h)+.1*np.kron(q, q)
    energies, vectors = np.linalg.eigh(H)
    pc, ph = [np.array([float(x) for x in populations[s]]) for s in ('cold', 'hot')]
    rc = (eigenbasis*pc) @ eigenbasis.conj().T
    rh = (eigenbasis*ph) @ eigenbasis.conj().T
    rho = np.kron(rc, rh)
    direct = {}
    for t in (300, 600):
        U = (vectors*np.exp(-1j*energies*t)) @ vectors.conj().T
        out = U@rho@U.conj().T
        mc, mh = [marginal(out, [4, 4], [i]) for i in (0, 1)]
        delta = [float(np.trace(h@(r-r0)).real) for r, r0 in ((mc, rc), (mh, rh))]
        old = a4['actual_flow_predictions'][f't{t}_lambda1']
        max_error = max(abs(delta[i]-float(F(old[k]))) for i, k in enumerate(('delta_EA_center','delta_EB_center')))
        assert max_error < 1e-10
        deficits = []
        for name, beta, out_r in (('cold', 1.2, mc), ('hot', 1/3, mh)):
            z = np.exp(-beta*e).sum()
            relative = -entropy(out_r)+beta*np.trace(h@out_r).real+np.log(z)
            mean_tau = np.dot(e, np.exp(-beta*e)/z)
            pin = 2*(np.trace(h@out_r).real-mean_tau)**2/float(width)**2
            assert relative >= pin-1e-12
            deficits.append(float(relative))
        direct[str(t)] = dict(delta_energy=delta, source_center_error=max_error,
                              endpoint_relative_entropies=deficits)

    # Independently test the *complete* balance with correlated initial data,
    # nonzero interactions and entropy-changing machinery/work storage.
    rng = np.random.default_rng(12120043)
    dims = [2, 3, 2, 2]  # hot, cold, clock/work carrier, machine
    d = int(np.prod(dims)); bc, bh = 1.2, 1/3
    local = [np.diag([0., 1.3]), np.diag([0., .4, 1.7]),
             np.diag([0., .8]), np.diag([0., .5])]
    bare = sum(embed(a, i, dims) for i, a in enumerate(local))
    errors, deleted_terms = [], []
    for k in range(24):
        x = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
        rin = x@x.conj().T; rin /= np.trace(rin)
        z = rng.normal(size=(d, d))+1j*rng.normal(size=(d, d))
        interaction = .07*(z+z.conj().T)
        full = bare+interaction
        es, vs = np.linalg.eigh(full)
        U = (vs*np.exp(-1j*es*(.1+.08*k)))@vs.conj().T
        rout = U@rin@U.conj().T
        ends = []
        for state in (rin, rout):
            ms = [marginal(state, dims, [j]) for j in range(4)]
            ss = [entropy(r) for r in ms]
            ee = [float(np.trace(a@r).real) for a, r in zip(local, ms)]
            dd = sum(-ss[j]+beta*ee[j]+np.log(np.exp(-beta*np.diag(local[j])).sum())
                     for j, beta in ((0, bh), (1, bc)))
            ends.append(dict(E=ee, S=ss, D=dd, I=sum(ss)-entropy(state)))
        x0, x1 = ends
        de = np.array(x1['E'])-x0['E']; ds = np.array(x1['S'])-x0['S']
        dv = float(np.trace(interaction@(rout-rin)).real)
        lhs = bc*(de[2]+de[3]+dv)-ds[2]-ds[3]
        rhs = (bc-bh)*(-de[0])+x0['D']-x1['D']+x0['I']-x1['I']
        assert abs(np.sum(de)+dv) < 1e-12
        errors.append(abs(lhs-rhs))
        deleted_terms.append(dict(initial_correlation=x0['I'],
                                  machine_and_interaction=bc*(de[3]+dv)-ds[3],
                                  clock_entropy=ds[2]))
    assert max(errors) < 1e-11
    assert max(abs(x['machine_and_interaction']) for x in deleted_terms) > .01
    assert min(x['initial_correlation'] for x in deleted_terms) > .1
    assert max(abs(x['clock_entropy']) for x in deleted_terms) > .001
    sources = [STAGE/str(n)/'results.json' for n in (1202, 1204, 1205, 1212)]
    return dict(status='PASS', all_checks_passed=True,
                source_sha256={p.relative_to(STAGE).as_posix(): sha(p) for p in sources},
                independent_code_sha256=sha(Path(__file__)),
                initial_energy_offsets={s:[str(x) for x in v] for s,v in offsets.items()},
                frozen_source_bounds=bounds, native_QQ_reconstruction=direct,
                generic_complete_balance_cases=len(errors), max_balance_error=max(errors),
                deleted_terms_material=True,
                qualifications=dict(matrix_checks_prove_general_identity=False,
                    actual_original_program_is_cycle=False, actual_work_engine_certified=False,
                    relative_free_energy_change_equals_extractable_work=False,
                    finite_factory_initial_debt_ignored=False,
                    finite_bath_penalty_is_preexisting_method=True))


if __name__ == '__main__':
    value = evaluate()
    out = HERE/'independent_results.json'
    if '--write' in sys.argv:
        assert not out.exists(), 'Refuse to overwrite an existing independent result'
        out.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert out.exists() and read(out) == value, 'Saved independent result differs'
    print(json.dumps(dict(status='PASS', all_checks_passed=True, results_match=True,
                         round=1212, actual_windows=2, generic_balance_cases=24)))
