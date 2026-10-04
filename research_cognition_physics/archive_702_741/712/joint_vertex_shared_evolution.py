"""712: original QQQL insertion, forced channels and quantum record jets.

Sparse CAR checks use original 32/64 modes. Pointwise bosonic checks do not
replace the full Hamiltonian or its Gibbs state; the note proves those links.
"""
import argparse
import hashlib
from itertools import product
import json
from pathlib import Path
import numpy as np
import joint_charge_changing_vertex as vertex

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_vertex_shared_evolution_results.json'
old = vertex.old
matter = old.matter
geom = matter.original
P = {ids: a/np.sqrt(vertex.NORM2) for ids, a in vertex.COEFF.items()}


def creator(poly, state):
    out = {}
    for ids, coefficient in poly.items():
        term = state
        for index in reversed(ids):
            term = old.create(term, index)
        out = old.add(out, term, coefficient)
    return out


def lift(poly, h):
    """[dGamma(h), pure creation polynomial], preserving canonical CAR signs."""
    result = {}
    for ids, coefficient in poly.items():
        for slot, index in enumerate(ids):
            for target in np.flatnonzero(abs(h[:, index]) > 1e-15):
                new = list(ids)
                new[slot] = int(target)
                if len(set(new)) < len(new):
                    continue
                key = tuple(sorted(new))
                result[key] = result.get(key, 0)+coefficient*h[target, index]*vertex.parity(new)
    return {k: v for k, v in result.items() if abs(v) > 1e-14}


def quadratic(state, h, delta):
    out = {}
    for i, j in zip(*np.nonzero(abs(h) > 1e-15)):
        out = old.add(out, old.create(old.annihilate(state, int(j)), int(i)), h[i, j])
    for i, j in zip(*np.nonzero(np.triu(abs(delta) > 1e-15, 1))):
        out = old.add(out, old.create(old.create(state, int(j)), int(i)), delta[i, j])
        out = old.add(out, old.annihilate(old.annihilate(state, int(i)), int(j)), delta[i, j].conjugate())
    return out


def difference(a, b):
    return vertex.distance(a, b)


def mass_norm(phi):
    r2 = np.dot(phi[:4], phi[:4])
    y = matter.Y
    return float(r2/geom.F(phi)*(1.5*(abs(y['u'])**2+abs(y['d'])**2)
                                   +.5*(abs(y['e'])**2+abs(y['nu'])**2)))


def channel_check():
    phi = np.array([.4, -.3, .2, .1, .35])
    other = np.array([-.2, .15, .1, -.25, .3])
    h, delta = matter.mass_matrices(phi)
    c = lift(P, h)
    got = sum(abs(a)**2 for a in c.values())
    expected = mass_norm(phi)
    assert abs(got-expected) < 1e-13
    module_norms = {}
    for name in ('u', 'd', 'e', 'nu'):
        sl = matter.SLICES[name]
        module_norms[name] = float(sum(abs(a)**2 for ids, a in c.items() if any(sl.start <= i < sl.stop for i in ids)))
        factor = 1.5 if name in ('u', 'd') else .5
        target = factor*abs(matter.Y[name])**2*sum(phi[:4]**2)/geom.F(phi)
        assert abs(module_norms[name]-target) < 1e-13
    states = [{0: 1.}, {1 << 30: 1.}, {(1 << 4)|(1 << 29): 1.}]
    residuals = []
    for s in states:
        lhs = old.add(quadratic(creator(P, s), h, delta), creator(P, quadratic(s, h, delta)), -1)
        residuals.append(difference(lhs, creator(c, s)))
    # Full 64-mode two-node coefficients, including all six original modules.
    h2, d2 = matter.mass_matrices(other)
    H = np.zeros((64, 64), complex)
    D = np.zeros_like(H)
    H[:32, :32], H[32:, 32:] = h, h2
    D[:32, :32], D[32:, 32:] = delta, d2
    C = matter.gauge.group_exp(np.array([.1, -.3, .07, .2, -.1, .04, .02, .13]), 3)
    W = matter.gauge.group_exp(np.array([.2, -.15, .3]), 2)
    R = matter.representation(C, W, np.exp(.17j))
    t = .23-.11j
    H[32:, :32] = t*R
    H[:32, 32:] = t.conjugate()*R.conj().T
    full = lift(P, H)
    hopping_norm = sum(abs(a)**2 for ids, a in full.items() if ids[-1] >= 32)
    assert abs(hopping_norm-4*abs(t)**2) < 1e-13
    full_norm = float(sum(abs(a)**2 for a in full.values()))
    assert abs(full_norm-expected-4*abs(t)**2) < 1e-13
    for s in ({0: 1.}, {(1 << 30)|(1 << 60): 1.}):
        lhs = old.add(quadratic(creator(P, s), H, D), creator(P, quadratic(s, H, D)), -1)
        residuals.append(difference(lhs, creator(full, s)))
    assert max(residuals) < 1e-12
    return dict(original_modes_per_node=32,full_two_node_modes=64,
        diagnostic_phi=phi.tolist(),mass_channel_norms=module_norms,
        exact_mass_norm_squared=expected,computed_mass_norm_squared=float(got),
        hopping_norm_squared=float(hopping_norm),expected_hopping_norm_squared=4*abs(t)**2,
        full_vacuum_commutator_norm_squared=full_norm,commutator_errors=residuals,
        original_vacuum_not_assumed_stationary=True,full_Gibbs_not_replaced_by_vacuum=True)


def pairing_replacement(poly, delta, state):
    """Leibniz action [pairing, C1]; no normal-ordering approximation."""
    out = {}
    for ids, coefficient in poly.items():
        for slot, index in enumerate(ids):
            for j in np.flatnonzero(abs(delta[index]) > 1e-15):
                term = state
                for k in reversed(range(len(ids))):
                    term = (old.annihilate(term, int(j)) if k == slot else old.create(term, ids[k]))
                out = old.add(out, term, coefficient*delta[index, j].conjugate())
    return out


def majorana_check():
    phi = np.array([.4, -.3, .2, .1, .35])
    h, d = matter.mass_matrices(phi)
    zero = np.zeros_like(h)
    c = lift(P, h)
    rows = []
    for index in (30, 31):
        s = {1 << index: 1.}
        first = old.add(quadratic(creator(P, s), zero, d), creator(P, quadratic(s, zero, d)), -1)
        second = old.add(quadratic(creator(c, s), zero, d), creator(c, quadratic(s, zero, d)), -1)
        expected = pairing_replacement(c, d, s)
        assert difference(first, {}) < 1e-13 and difference(second, expected) < 1e-13
        norm2 = float(old.inner(second, second).real)
        assert norm2 > 1e-7
        # Mixed creators/annihilator preserve quark-charge +3, not total CAR number +4.
        assert all(vertex.nq(mask)-vertex.nq(next(iter(s))) == 3 for mask in second)
        assert all(mask.bit_count()-1 == 2 for mask in second)
        rows.append(dict(input_neutrino_mode=index,first_pair_commutator_norm=difference(first, {}),
                         second_pair_commutator_norm_squared=norm2,replacement_error=difference(second, expected)))
    return dict(original_pair_coefficient_real=float(d[30, 31].real),
        original_pair_coefficient_imag=float(d[30, 31].imag),checks=rows,
        no_new_Majorana_parameter=True,quark_charge_change=3,total_number_change_in_this_channel=2)


def derivative_matrices(phi):
    h, _ = matter.mass_matrices(phi)
    f = geom.F(phi)
    result = []
    for j in range(5):
        basis = np.eye(5)[j]
        numerator, _ = matter.mass_matrices(basis)
        numerator *= np.sqrt(geom.F(basis))
        result.append(numerator/np.sqrt(f)+phi[j]/(6*f)*h)
    return np.array(result)


def geometry_check():
    rng = np.random.default_rng(7123)
    rows = []
    w = .73**3*np.exp(6*.12)
    for _ in range(5):
        phi = rng.normal(size=5)*.3
        h, _ = matter.mass_matrices(phi)
        dh = derivative_matrices(phi)
        inv = geom.inverse(phi)/w
        ds = inv[4]
        df = np.r_[2*phi[:4], 0.]@inv
        singlet = np.einsum('i,ijk->jk', ds, dh)
        higgs = np.einsum('i,ijk->jk', df, dh)
        expected = 2*geom.F(phi)/w*h
        err_s = float(np.linalg.norm(singlet))
        err_h = float(np.linalg.norm(higgs-expected))
        eps = 1e-5
        plus = matter.mass_matrices(phi+eps*df)[0]
        minus = matter.mass_matrices(phi-eps*df)[0]
        fd_error = float(np.linalg.norm((plus-minus)/(2*eps)-expected)/max(np.linalg.norm(expected), 1))
        assert max(err_s, err_h) < 1e-12 and fd_error < 2e-8
        # A coordinate derivative ds h is usually nonzero; metric contraction cancels it.
        naive = float(np.linalg.norm(dh[4]))
        assert naive > 1e-4
        c = lift(P, h)
        record = lift(P, higgs)
        scale = 2*geom.F(phi)/w
        error = difference(creator(record, {0: 1.}), {m: scale*a for m, a in creator(c, {0: 1.}).items()})
        assert error < 1e-12
        rows.append(dict(phi=phi.tolist(),singlet_metric_contraction_error=err_s,
            nonzero_naive_partial_s_norm=naive,higgs_metric_identity_error=err_h,
            directional_finite_difference_relative_error=fd_error,
            original_CAR_record_jet_error=error,higgs_jet_multiplier=float(scale)))
    return dict(original_H5_inverse_metric=True,cell_weight=w,checks=rows,
        original_s_record_second_jet_zero=True,
        higgs_radius_is_comparison_observable_not_installed_instrument=True,
        only_second_order_time_jets_claimed=True)


def run():
    results = [channel_check(), majorana_check(), geometry_check()]
    deps = ('research_note_574.md', 'research_note_598.md', 'research_note_614.md',
            'research_note_623.md', 'research_note_629.md', 'research_note_710.md',
            'research_note_711.md', 'joint_charge_changing_vertex.py',
            'joint_fermion_gauss_completion.py', 'joint_curved_quantum_source.py',
            'round712_drafts/minimal_left_vertex_entry.md',
            'round712_drafts/minimal_left_vertex_entry_results.json')
    return dict(round=712,tests_run=3,failures=0,errors=0,
        forced_channels=results[0],majorana_second_step=results[1],actual_geometry_record=results[2],
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original finite full H0; exact insertion commutators and local-core record jets. Thermal sum rule is analytic. No instanton coefficient, chiral continuum, full Gibbs simulation or all-time fixed-background closure.')


if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as output:json.dump(result,output,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
