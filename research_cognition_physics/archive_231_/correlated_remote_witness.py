"""Unnumbered remote-witness interface derived from the round 503 receipt.

No new numbered physical result: explicit additional local coupling, known
memory preparation, and a history-dependent encoding. No autonomous controller,
event IDs at the remote party, or space dimension is generated.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import io
import itertools
import json
from pathlib import Path
import unittest

import numpy as np
import uniform_identity_receipt as old
import all_scale_monitored_graph_source as family_model

HERE = Path(__file__).resolve().parent
TARGET = HERE/'correlated_remote_witness_results.json'
OBS = {}


def colour_bit(colour):
    return 0 if colour == 0 else 1 << (colour-1)


def action(trees, flips, n, v, graph, colour, memory):
    """Physical single-packet terms, independent of any encoded-space matrix."""
    out = defaultdict(int)
    for gp, weight in enumerate(flips[:, graph]):
        if weight:
            out[(v, gp, colour, memory)] += int(weight)
    degree = 0
    for a, b in trees[graph]:
        if v in (a, b):
            dest = b if v == a else a
            out[(dest, graph, colour, memory)] += 1
            degree += 1
    out[(v, graph, colour, memory)] -= degree
    if colour in (0, v+1):
        dest = v+1 if colour == 0 else 0
        out[(v, graph, dest, memory ^ (1 << v))] += 1
    return {row:w for row,w in out.items() if w}


def encoded_rows(n, ng, sector):
    return np.asarray([((v*ng+g)*(n+1)+c)*(1 << n)
                       + (sector ^ colour_bit(c))
                       for v in range(n) for g in range(ng) for c in range(n+1)])


def unitary(h, t):
    values, vectors = np.linalg.eigh(h)
    return (vectors*np.exp(-1j*t*values))@vectors.conj().T


def encode(x, n, ng, sector):
    answer = np.zeros((n*ng*(n+1)*(1 << n), x.shape[1]), complex)
    answer[encoded_rows(n, ng, sector)] = x
    return answer


class Audit(unittest.TestCase):
    def test_01_all_sectors_exact_local_generator_and_charges(self):
        rows = []
        for internal in (1, 2):
            n, trees, transitions = family_model.family(internal)
            ng = len(trees)
            flips = np.zeros((ng, ng), dtype=int)
            for g, transitions_g in enumerate(transitions):
                for gp, weight in transitions_g.items():
                    flips[gp, g] = weight
            h = old.model(trees, flips, n)
            checked = 0
            for sector in range(1 << n):
                for column in range(len(h)):
                    v, rem = divmod(column, ng*(n+1))
                    g, c = divmod(rem, n+1)
                    memory = sector ^ colour_bit(c)
                    actual = action(trees, flips, n, v, g, c, memory)
                    expected = {}
                    for row in np.flatnonzero(h[:, column]):
                        vp, remp = divmod(int(row), ng*(n+1))
                        gp, cp = divmod(remp, n+1)
                        expected[(vp, gp, cp, sector ^ colour_bit(cp))] = int(h[row, column])
                    self.assertEqual(actual, expected)
                    # Exact conserved signs, on every physical basis state.
                    for vp, gp, cp, mp in actual:
                        self.assertEqual(mp ^ colour_bit(cp), memory ^ colour_bit(c))
                    checked += 1
            indices = np.concatenate([encoded_rows(n, ng, m) for m in range(1 << n)])
            np.testing.assert_array_equal(np.sort(indices), np.arange(len(h)*(1 << n)))
            rows.append(dict(vertices=n, graphs=ng, sectors=1 << n,
                             physical_columns_checked=checked,
                             physical_active_dimension=len(h)*(1 << n)))
        OBS['integer_full_sector_certificates'] = rows

    def test_02_physical_three_step_history_with_local_resets(self):
        # Three graph basis states and an unknown reference; this small diagnostic
        # uses static graphs. Actual NNI generator equivalence is checked above.
        n = 3
        trees = sorted({old.graphs.prufer_tree(n, (v,)) for v in range(n)},
                       key=lambda g:sorted(g))
        ng, nc, nm = len(trees), n+1, 1 << n
        flips = np.zeros((ng, ng), dtype=int)
        h = old.model(trees, flips, n)
        full = np.zeros((len(h)*nm, len(h)*nm), dtype=int)
        for v, g, c, memory in itertools.product(range(n), range(ng), range(nc), range(nm)):
            col = ((v*ng+g)*nc+c)*nm+memory
            for (vp, gp, cp, mp), weight in action(trees, flips, n, v, g, c, memory).items():
                row = ((vp*ng+gp)*nc+cp)*nm+mp
                full[row, col] += weight
        np.testing.assert_array_equal(full, full.T)
        u, uf = unitary(h, 0.2), unitary(full, 0.2)
        rng = np.random.default_rng(20260930)
        gr = rng.normal(size=(ng, 2))+1j*rng.normal(size=(ng, 2))
        gr /= np.linalg.norm(gr)
        x = np.zeros((len(h), 2), complex)
        x[np.arange(ng)*nc] = gr
        reference = gr.T@gr.conj()
        branches = [(x, encode(x, n, ng, 0), 0)]
        max_error = 0.0
        statistics = []
        for depth in range(1, 4):
            fresh = []
            for a, physical, sector in branches:
                a, physical = u@a, uf@physical
                for colour in range(-1, nc):
                    mask = np.zeros(len(h), dtype=bool)
                    if colour == -1:
                        mask[ng*nc:] = True
                    else:
                        mask[np.arange(ng)*nc+colour] = True
                    aa = a*mask[:, None]
                    pp = physical*np.repeat(mask, nm)[:, None]
                    nxt = sector
                    if colour >= 0:
                        # Observed-colour-controlled local swap c <-> q at root.
                        order = np.arange(len(h))
                        for g in range(ng):
                            q, r = g*nc, g*nc+colour
                            order[q], order[r] = order[r], order[q]
                        physical_order = (order[:, None]*nm+np.arange(nm)).reshape(-1)
                        aa, pp = aa[order], pp[physical_order]
                        nxt ^= colour_bit(colour)
                    error = float(np.linalg.norm(pp-encode(aa, n, ng, nxt)))
                    max_error = max(max_error, error)
                    self.assertLess(error, 2e-12)
                    fresh.append((aa, pp, nxt))
            branches = fresh
            mass = sum(float(np.vdot(p, p).real) for _,p,_ in branches)
            rr = sum((p.T@p.conj() for _,p,_ in branches), np.zeros((2, 2), complex))
            self.assertAlmostEqual(mass, 1.0, places=11)
            self.assertLess(float(np.linalg.norm(rr-reference)), 2e-11)
            statistics.append(dict(depth=depth, all_histories=len(branches),
                                   total_mass=mass, reference_residual=float(np.linalg.norm(rr-reference))))
        OBS['complete_physical_history'] = dict(
            static_graph_diagnostic=True, no_postselection=True, histories=statistics,
            full_encoded_state_max_error=max_error,
            arbitrary_reference_checked=True,
            not_a_claim_of_unchanged_unencoded_activity_marginal=True)

    def test_03_remote_value_window_and_inherited_source_budget(self):
        trees, _, _, flips = old.old500.six_model()
        n, ng = 6, len(trees)
        h = old.model(trees, flips, n)
        u = unitary(h, 1/8)
        beta = F(3, 128)
        max_loss = 0.0
        for root in range(n):
            for b in range(n):
                if root == b:
                    continue
                cols = (root*ng+np.arange(ng))*(n+1)+b+1
                rows = np.asarray([k for k in range(len(h)) if k % (n+1) != b+1])
                block = u[np.ix_(rows, cols)]
                loss = float(np.linalg.eigvalsh(block.conj().T@block)[-1])
                max_loss = max(max_loss, loss)
                self.assertLessEqual(loss, float(beta*beta)+1e-12)
        degree, f, v = 3, 24, F(7)
        c = F(1, 6)
        remainder = v*v*48/4+v**4/12
        t = 1/(16*(remainder/c+f))
        epsilon = 4*(remainder/c+f)**2*t*t
        self.assertLessEqual(v*t, 1)
        self.assertLessEqual(t, c/(2*remainder))
        self.assertEqual(epsilon, F(1, 64))
        OBS['inherited_bounds'] = dict(
            source_step=str(t), all_size_source_lower=str(c*c*t**6/4),
            conditional_current_nonedge_upper=str(epsilon),
            remote_value_window='1/8', value_loss_upper=str(beta*beta),
            six_graph_all_roots_targets=30, numerical_max_value_loss=max_loss,
            no_new_source_rate_or_preservation_theorem_claim=True)

    def test_04_unknown_initial_flags_and_no_free_mutual_acknowledgment(self):
        n, ng, root, b = 2, 1, 0, 1
        d = n*ng*(n+1)
        old_reply = np.zeros((d, 1))
        old_reply[root*(n+1)+b+1, 0] = 1
        w0 = encode(old_reply, n, ng, 0)
        wb = encode(old_reply, n, ng, 1 << b)
        self.assertEqual(int(np.argmax(w0)) % (1 << n), 1 << b)
        self.assertEqual(int(np.argmax(wb)) % (1 << n), 0)
        # Same colour readout but opposite remote flag: known initial sector matters.
        reset = np.zeros_like(w0)
        reset[root*(n+1)*(1 << n)+(1 << b), 0] = 1
        self.assertEqual(float(np.linalg.norm(reset[encoded_rows(n, ng, 0)])), 0.0)
        query = np.zeros_like(old_reply); query[root*(n+1), 0] = 1
        np.testing.assert_array_equal(reset, encode(query, n, ng, 1 << b))
        # The same M_b=1 occurs while the packet is at b, before a root return.
        remote_packet = np.zeros_like(old_reply)
        remote_packet[b*(n+1)+b+1, 0] = 1
        self.assertEqual(int(np.argmax(encode(remote_packet,n,ng,0))) % (1 << n), 1 << b)
        self.assertEqual(float(np.vdot(old_reply, remote_packet).real), 0.0)
        OBS['scope_witnesses'] = dict(
            same_reply_opposite_flag_for_unknown_sector=True,
            local_reset_leaves_zero_sector_but_enters_known_xor_sector=True,
            remote_flag_alone_does_not_announce_root_receipt=True,
            unmonitored_remote_reader_not_free=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    dependencies = ['uniform_identity_receipt.py', 'uniform_identity_receipt_results.json',
        'research_note_503.md', 'research_note_502.md', 'research_note_514.md',
        'internal_implementation_scope_review.md', 'joint_event_conjecture_review.md']
    return dict(date='2026-09-30', scientific_baseline_round=520,
        numbered_round_created=False, numbered_scientific_test_increment=0,
        diagnostic_tests=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        dependency_sha256={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in dependencies},
        scope=dict(explicit_added_local_coupling=True, original_510_H_unchanged_claim=False,
            all_finite_size_encoding_proof=True, root_history_tracks_memory_sector=True,
            remote_event_identifier_or_autonomous_read_trigger_generated=False,
            spatial_dimension_or_GR_derived=False),
        observations=OBS.copy())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    answer = json.loads(json.dumps(run()))
    if args.check:
        assert answer == json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            json.dump(answer, f, ensure_ascii=False, indent=2); f.write('\n')
    print(json.dumps(answer, ensure_ascii=False, indent=2))
