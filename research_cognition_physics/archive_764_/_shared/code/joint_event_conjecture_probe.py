"""Unnumbered audit: shared classical records do not imply the proposed physics.

These small exact matrix identities support the review, not an autonomous
universe, a new quantum reconstruction, or a countermodel to all FUCP premises.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_event_conjecture_probe_results.json'
OBS = {}


def gate(control, target, bits=3):
    out = np.zeros((2**bits, 2**bits), dtype=int)
    for value in range(2**bits):
        destination = value ^ ((1 << (bits-1-target))
                               if value & (1 << (bits-1-control)) else 0)
        out[destination, value] = 1
    return out


def swap_readers():
    out = np.zeros((8, 8), dtype=int)
    for value in range(8):
        bits = [(value >> shift) & 1 for shift in (2, 1, 0)]
        dest = 4*bits[0] + 2*bits[2] + bits[1]
        out[dest, value] = 1
    return out


def partial_transpose(x, bit, count=3):
    axes = list(range(2*count))
    axes[bit], axes[count+bit] = axes[count+bit], axes[bit]
    return x.reshape([2]*(2*count)).transpose(axes).reshape(x.shape)


class Audit(unittest.TestCase):
    def test_01_faithful_reversible_shared_record_without_entanglement(self):
        # Event E and two initially blank records A,B. Integer weights keep
        # all preparation, copying, symmetry and readout identities exact.
        u = gate(0, 2) @ gate(0, 1)
        initial = np.zeros((8, 8), dtype=int)
        initial[0, 0], initial[4, 4] = 1, 2
        final = u @ initial @ u.T
        expected = np.zeros_like(final)
        expected[0, 0], expected[7, 7] = 1, 2
        np.testing.assert_array_equal(final, expected)
        np.testing.assert_array_equal(u.T @ final @ u, initial)
        s = swap_readers()
        np.testing.assert_array_equal(s @ u @ s.T, u)
        np.testing.assert_array_equal(s @ final @ s.T, final)
        for bit in range(3):
            np.testing.assert_array_equal(partial_transpose(final, bit), final)
        # Separability follows from the explicit product-state mixture, NOT
        # from the partial-transpose diagnostic alone.
        for event in (0, 1):
            self.assertEqual(int(np.argmax(u[:, 4*event])), 7*event)
        rho = final / 3
        entropy = -sum(p*np.log2(p) for p in (1/3, 2/3))
        self.assertAlmostEqual(float(np.trace(rho @ rho)), 5/9)
        OBS['shared_record'] = dict(
            exact_integer_identities=True, faithful_for_both_events=True,
            reader_exchange_covariant=True, separable_product_mixture=True,
            protocol_reversible=True, global_entropy_before_after_bits=float(entropy),
            no_unknown_nonorthogonal_state_cloning_claim=True)

    def test_02_ghz_redundancy_does_not_detect_all_local_tampering(self):
        # Unnormalised vectors: integer calculations avoid a tolerance proof.
        plus = np.zeros(8, dtype=int); plus[0] = plus[7] = 1
        minus = plus.copy(); minus[7] = -1
        z_a = np.diag([1 if not (v & 2) else -1 for v in range(8)])
        np.testing.assert_array_equal(z_a @ plus, minus)
        self.assertEqual(int(plus @ minus), 0)
        np.testing.assert_array_equal(plus**2, minus**2)
        # A classical bit flip IS detected by comparing the three records.
        flip_a = np.eye(8, dtype=int)[[v ^ 2 for v in range(8)]]
        corrupted = flip_a @ (plus**2)
        self.assertEqual(int(corrupted[0] + corrupted[7]), 0)
        OBS['tampering_scope'] = dict(
            local_phase_flip_changes_to_orthogonal_GHZ=True,
            complete_Z_record_distribution_unchanged=True,
            classical_value_flip_detected_by_record_comparison=True,
            threat_model_and_verification_menu_required=True)

    def test_03_shared_facts_allow_private_quantum_information(self):
        # The event register itself is internal. R can purify its classical
        # mixture. The diagonal shared record is not an assertion that the
        # entire closed universe has lost all quantum coherence.
        phi = np.zeros((8, 2), dtype=int)
        phi[0, 0] = phi[7, 1] = 1
        marginal = phi @ phi.T
        expected = np.zeros((8, 8), dtype=int)
        expected[0, 0] = expected[7, 7] = 1
        np.testing.assert_array_equal(marginal, expected)
        record = expected / 2
        private0 = np.diag([1, 0]); private1 = np.diag([0, 1])
        joint0 = np.kron(record, private0)
        joint1 = np.kron(record, private1)
        for joint in (joint0, joint1):
            marginal = np.einsum('aibi->ab', joint.reshape(8, 2, 8, 2))
            np.testing.assert_array_equal(marginal, record)
        distance = np.sum(np.abs(np.linalg.eigvalsh(joint0-joint1))) / 2
        self.assertEqual(float(distance), 1.0)
        OBS['private_information'] = dict(
            identical_complete_shared_record_law=True,
            orthogonal_private_inputs_trace_distance=float(distance),
            internal_purification_explicit=True,
            task_agreement_does_not_require_global_state_learning=True)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    dependencies = [HERE.parents[1]/'猜想'/'2029-9-29.md',
                    HERE/'research_note_502.md', HERE/'research_note_516.md',
                    HERE/'research_note_520.md',
                    HERE.parent/'可组合认知结构与有限维量子理论_阶段论文.md']
    return dict(date='2026-09-30', scientific_baseline_round=520,
        numbered_round_created=False, numbered_scientific_test_increment=0,
        diagnostic_tests=result.testsRun, failures=len(result.failures),
        errors=len(result.errors),
        dependency_sha256={str(p.relative_to(HERE.parents[1])).replace('\\','/'):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in dependencies},
        observations=OBS.copy(), scope=dict(
            minimal_shared_event_axiom_does_not_force_entangled_records=True,
            full_FUCP_reconstruction_refuted=False,
            autonomous_internal_implementation_generated=False,
            spacetime_dimension_or_GR_or_standard_model_derived=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    answer = json.loads(json.dumps(run()))
    if args.check:
        assert answer == json.loads(TARGET.read_text(encoding='utf8'))
    else:
        with TARGET.open('x', encoding='utf8', newline='\n') as f:
            json.dump(answer, f, ensure_ascii=False, indent=2)
            f.write('\n')
    print(json.dumps(answer, ensure_ascii=False, indent=2))
