"""Round 455: a three-exchange conversion of an unknown relational subject.

Baseline 453. Existing three-spin subsystem and four-spin singlet encodings
are reused. The supplied singlet is an explicit internal preparation resource;
the fixed Hamiltonian gives an endpoint and a finite window, not a permanent
handoff or an autonomous preparation factory.
"""
import argparse
from functools import lru_cache
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest

import numpy as np

import exchange_relation_audit as old
import encoded_exchange_response_audit as four

TARGET = Path(__file__).with_name("relational_interface_conversion_audit_results.json")
OBS = {}


def short(value):
    return float(f"{float(value):.12g}")


def local(n, position, operator):
    out = np.ones((1, 1), dtype=complex)
    for i in range(n):
        out = np.kron(out, operator if i == position else np.eye(2))
    return out


@lru_cache(None)
def model():
    v3 = old.encoding()  # columns G,L
    singlet = np.array([0, 1, -1, 0], complex)/math.sqrt(2)
    vin = np.kron(v3, singlet[:, None])
    # First four physical spins are old A=(0,1,2) and new helper 3.
    w4 = np.column_stack([
        (np.kron(v3[:, ell], [0, 1])-np.kron(v3[:, 2+ell], [1, 0]))/math.sqrt(2)
        for ell in range(2)])
    target = np.column_stack([
        np.kron(w4[:, ell], np.eye(2)[:, gauge])
        for gauge, ell in itertools.product(range(2), repeat=2)])
    embedding = np.kron(v3, np.eye(4))  # order G,L,helper3,carrier4
    h = sum(old.swap(5, a, 4) for a in range(3))
    swap_g4 = old.swap(4, 0, 3)
    pout = np.kron(w4@w4.conj().T, np.eye(2))
    pin = np.kron(v3@v3.conj().T, np.eye(4))
    return dict(v3=v3, s=singlet, vin=vin, w4=w4, target=target,
                embedding=embedding, h=h, swap_g4=swap_g4, pout=pout, pin=pin)


class Audit(unittest.TestCase):
    def close(self, left, right, tolerance=4e-12):
        self.assertLess(float(np.linalg.norm(left-right)), tolerance)

    def test_01_same_old_codes_and_collective_swap_identity(self):
        m = model()
        self.close(m["w4"], four.code_block())
        self.close(m["vin"].conj().T@m["vin"], np.eye(4))
        self.close(m["target"].conj().T@m["target"], np.eye(4))
        e = m["embedding"]
        self.close(m["h"]@e, e@(np.eye(16)+m["swap_g4"]))
        self.close(m["h"]@m["pin"], m["pin"]@m["h"])
        for mu in old.PAULI:
            jtotal = sum(local(5, i, mu) for i in range(5))
            self.close(m["h"]@jtotal, jtotal@m["h"])
            jfour = sum(local(4, i, mu) for i in range(4))
            self.close(jfour@m["w4"], np.zeros((16, 2)))
        OBS["exact_collective_identity"] = dict(raw_dimension=32,
            original_code_plus_aux_dimension=16, input_logical_gauge_dimension=4,
            physical_edges=[[0, 4], [1, 4], [2, 4]], all_weights=1,
            original_code_preserved_all_times=True,
            compressed_generator="I + SWAP_(G,carrier4)",
            four_spin_basis_identical_to_round_434=True,
            common_SU2_invariance=True)

    def test_02_all_unknown_inputs_and_reference_exact_conversion(self):
        m = model()
        u = old.evolve(m["h"], math.pi/2)
        self.close(u@m["vin"], m["target"])
        # Matrix-column identity proves every reference dimension; this mixed
        # GLR instance is an independent check, not the universal argument.
        rng = np.random.default_rng(455)
        rho = old.density(12, rng)
        left = np.kron(u@m["vin"], np.eye(3))
        right = np.kron(m["target"], np.eye(3))
        actual, ideal = left@rho@left.conj().T, right@rho@right.conj().T
        self.close(actual, ideal)
        self.close(old.partial(actual, [16, 2, 3], [1, 2]),
                   old.partial(rho, [2, 2, 3], [0, 2]))
        reduced_four = old.partial(actual, [16, 2, 3], [0, 2])
        logic_ref = old.partial(rho, [2, 2, 3], [1, 2])
        wref = np.kron(m["w4"], np.eye(3))
        self.close(reduced_four, wref@logic_ref@wref.conj().T)
        # The entire G-L-reference state is recovered by the output isometry.
        self.close(right.conj().T@actual@right, rho)
        OBS["endpoint"] = dict(time="pi/2", isometry="U_T V_in = V_out",
            original_gauge_retained_in_physical_spin=4,
            unknown_G_L_reference_correlations_preserved=True,
            mixed_reference_check_dimension=3, postselection=False,
            unknown_input_tomography=False, intermediate_pulses=False,
            output_joint_error=short(np.linalg.norm(actual-ideal)))

    def test_03_readable_logic_survives_the_conversion_all_times(self):
        m = model()
        _, _, _, _, x3, y3, z3 = old.operators()
        effects = [(np.eye(8)+q)/2 for q in (x3, y3, z3)]
        # Z can also be the ordinary physical singlet comparison (I-S01)/2.
        singlet_effect = (np.eye(8)-old.swap(3, 0, 1))/2
        errors = []
        for effect, pauli in zip(effects, old.PAULI):
            self.assertGreater(np.linalg.eigvalsh(effect).min(), -1e-14)
            self.assertLess(np.linalg.eigvalsh(effect).max(), 1+1e-14)
            raw = np.kron(effect, np.eye(4))
            self.close(raw@m["h"], m["h"]@raw)
            expected = np.kron(np.eye(2), (np.eye(2)+pauli)/2)
            self.close(m["vin"].conj().T@raw@m["vin"], expected)
            self.close(m["target"].conj().T@raw@m["target"], expected)
            errors.append(float(np.linalg.norm(raw@m["h"]-m["h"]@raw)))
        self.close(m["v3"].conj().T@singlet_effect@m["v3"],
                   np.kron(np.eye(2), (np.eye(2)+old.PAULI[2])/2))
        self.close(m["w4"].conj().T@np.kron(singlet_effect, np.eye(2))@m["w4"],
                   (np.eye(2)+old.PAULI[2])/2)
        OBS["logic_effects"] = dict(unchanged_physical_effects_for_all_three_components=True,
            logic_channel_exact_identity_for_every_time=True,
            max_commutator_residual=short(max(errors)),
            effect_permissions_input=True,
            internal_pointer_or_full_tomography_implemented=False)

    def test_04_finite_window_full_channel_error_and_nonpermanent_code(self):
        m = model()
        values = [-.37, -.125, 0., .125, .52]
        errors = []
        for delta in values:
            u = old.evolve(m["h"], math.pi/2+delta)@m["vin"]
            overlap = m["target"].conj().T@u
            scalar = np.exp(-1j*delta)*(math.cos(delta)-.5j*math.sin(delta))
            self.close(overlap, scalar*np.eye(4))
            success = u.conj().T@m["pout"]@u
            self.close(success, (1-.75*math.sin(delta)**2)*np.eye(4))
            reference_vector = np.eye(4).reshape(-1)/2
            va = np.kron(u, np.eye(4))@reference_vector
            vb = np.kron(m["target"], np.eye(4))@reference_vector
            distance = old.distance(np.outer(va, va.conj()), np.outer(vb, vb.conj()))
            exact = math.sqrt(3)/2*abs(math.sin(delta))
            self.assertAlmostEqual(distance, exact, places=12)
            errors.append(abs(distance-exact))
        # Exact code occupancy has period pi and returns to 1/4 at t=pi.
        later = old.evolve(m["h"], math.pi)@m["vin"]
        self.close(later.conj().T@m["pout"]@later, np.eye(4)/4)
        OBS["finite_window"] = dict(center="pi/2", half_width="1/8",
            code_probability="1 - (3/4)*cos(t)^2",
            code_probability_window_lower="253/256",
            half_diamond_error="sqrt(3)/2 * abs(sin(delta))",
            half_diamond_window_upper="sqrt(3)/16",
            maximum_numeric_distance_error=short(max(errors)),
            exact_permanent_handoff=False,
            code_probability_at_pi=.25)

    def test_05_minimal_internal_carriers_and_pure_resource_boundary(self):
        m = model()
        ps = np.outer(m["s"], m["s"].conj())
        u = old.evolve(m["h"], math.pi/2)
        rho_old = m["v3"]@m["v3"].conj().T/4
        cases = []
        for p in (0., .1, .25, .7, 1.):
            tau = p*ps+(1-p)*(np.eye(4)-ps)/3
            initial = np.kron(rho_old, tau)
            final = u@initial@u.conj().T
            success = np.trace(m["pout"]@final).real
            self.assertAlmostEqual(success, p, places=12)
            ideal = m["target"]@np.eye(4)@m["target"].conj().T/4
            self.assertAlmostEqual(old.distance(final, ideal), 1-p, places=12)
            eigenvalues = np.linalg.eigvalsh(initial)
            max_rank4_weight = sum(eigenvalues[-4:])
            bound = max(p, (1-p)/3)
            self.assertAlmostEqual(max_rank4_weight, bound, places=12)
            input_rank = int(np.linalg.matrix_rank(initial, tol=1e-10))
            if p < 1:
                self.assertGreater(input_rank, 4)
                self.assertLess(bound, 1)
            cases.append(dict(singlet_fraction=p, input_rank=input_rank,
                fixed_conversion_success=short(success),
                any_unitary_rank4_capture_upper=short(bound)))
        # The full 32-dimensional maximally mixed state cannot be cooled by U.
        self.close(u@(np.eye(32)/32)@u.conj().T, np.eye(32)/32)
        self.assertAlmostEqual(np.trace(m["pout"]).real/32, 1/8)
        OBS["resources"] = dict(minimum_total_physical_spins=5,
            minimum_added_spins_for_closed_full_G_L_retention=2,
            minimum_complementary_carrier_dimension=2,
            invariant_two_spin_pure_resource_unique="singlet",
            deterministic_five_spin_conversion_requires_pure_aux=True,
            no_additional_sink_measurement_postselection_in_bound=True,
            maximally_mixed_full_state_code_probability="1/8",
            mixed_invariant_resource_cases=cases,
            source_of_supplied_singlet_derived=False,
            noisy_resource_endpoint_half_diamond_error="1-p",
            thermodynamic_work_minimum_claimed=False)

    def test_06_pure_resource_correlation_rerouting_and_energy_ledger(self):
        m = model()
        # G entangled with a two-dimensional internal R, L=(|0>+i|1>)/sqrt(2).
        bell = np.array([1, 0, 0, 1], complex)/math.sqrt(2)
        logic = np.array([1, 1j])/math.sqrt(2)
        source = np.array([bell[2*g+r]*logic[ell]
                           for g, ell, r in itertools.product(range(2), repeat=3)])
        output = np.kron(m["target"], np.eye(2))@source
        rho = np.outer(output, output.conj())
        self.close(old.partial(rho, [16, 2, 2], [1, 2]),
                   np.outer(bell, bell.conj()))
        # Old helper pair no longer contains the original free singlet.
        aux_pair = old.partial(rho, [2]*6, [3, 4])
        self.close(aux_pair, np.eye(4)/4)
        self.close(m["vin"].conj().T@m["h"]@m["vin"], 1.5*np.eye(4))
        self.close(m["vin"].conj().T@m["h"]@m["h"]@m["vin"], 3*np.eye(4))
        OBS["internal_ledger"] = dict(raw_spins=5, supplied_singlet_pairs=1,
            old_gauge_reference_Bell_pair_moved_exactly=True,
            former_aux_pair_singlet_probability_after_this_input=.25,
            singlet_resource_returned_unchanged=False,
            global_information_erased=False, initial_and_final_mean_H="3/2",
            variance_H="3/4", relative_energy_unit="unit coefficient of each SWAP",
            equal_coupling_layout_and_access_are_inputs=True,
            automatic_handoff_to_round_434_implemented=False)


def run():
    OBS.clear()
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise AssertionError(stream.getvalue())
    return dict(round=455, baseline_round=453, date="2026-09-24",
        runtime=dict(python=platform.python_version(), numpy=np.__version__),
        tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
        observations=OBS.copy(),
        scope=dict(existing_three_spin_unknown_G_L_input=True,
            invariant_singlet_aux_resource_explicit=True,
            exact_three_constant_exchange_conversion=True,
            every_internal_reference_preserved=True,
            finite_window_full_channel_error_proved=True,
            minimum_carrier_and_closed_five_spin_purity_conditions_proved=True,
            no_new_abstract_encoding_universality_claimed=True,
            depends_on_round_454=False,
            singlet_preparation_source_derived=False,
            permanent_autonomous_handoff=False,
            macro_subject_entangling_network_implemented=False,
            full_soca_or_spatial_dimension_or_GR_generated=False,
            full_GR_goal_completed=False, phase_closure_triggered=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    report = run()
    if not args.dry_run:
        with TARGET.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))
