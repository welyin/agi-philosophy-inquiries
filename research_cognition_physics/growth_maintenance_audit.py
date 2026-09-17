"""Round 206: sensitivity of finite growth to declared maintenance costs.

This is an exact rational budget simulation, not quantum dynamics. A and B
are abstract implementations, not real/complex labels. Their task equivalence
and cost functions are inputs. A switch does not establish emergence of J.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

from reference_maintenance_budget import maximum_age, refresh_resources
from reusable_network_reference import resources


@dataclass(frozen=True)
class CostLaw:
    coefficient: F
    exponent: int

    def __post_init__(self):
        if self.coefficient <= 0 or not isinstance(self.exponent, int) or self.exponent < 0:
            raise ValueError("Positive coefficient and nonnegative integer exponent required.")

    def value(self, size):
        if not isinstance(size, int) or size < 1:
            raise ValueError("Positive integer active size required.")
        return F(self.coefficient)*size**self.exponent


def fraction(value):
    value = F(value)
    return {"exact": str(value), "diagnostic": float(value)}


def simulate(laws, initial_mode="A", total_credits=1000, capacity=16,
             switch_cost=6, growth_cost=2, max_steps=64):
    """Greedy single-epoch switch-and-grow controller; no optimality claim.

    Active size starts at one, all dormant slots are already inside the whole.
    Quota 4*n is processing capacity, not newly created consumable stock.
    Growth adds one declared new capability and keeps old labels. Whether
    quantum implementations satisfy that contract is outside this model.
    """
    if initial_mode not in laws or not isinstance(capacity, int) or capacity < 1:
        raise ValueError("Known initial mode and positive capacity required.")
    if not isinstance(max_steps, int) or max_steps < 0:
        raise ValueError("Nonnegative integer horizon required.")
    stock, switch_cost, growth_cost = map(F, (total_credits, switch_cost, growth_cost))
    if stock < 0 or switch_cost < 0 or growth_cost <= 0:
        raise ValueError("Nonnegative stock/switch cost and positive growth cost required.")
    initial_stock = stock
    size, mode, events = 1, initial_mode, []
    status = "finite_horizon_reached"
    for step in range(max_steps):
        if size == capacity:
            status = "internal_candidate_slots_exhausted"
            break
        quota = F(4*size)
        options = []
        for index, (candidate, law) in enumerate(laws.items()):
            setup = F(0) if candidate == mode else switch_cost
            cost = law.value(size)+setup+growth_cost
            if cost <= quota and cost <= stock:
                options.append((cost, candidate != mode, index, candidate, setup))
        before = size
        if options:
            cost, changed, _, selected, setup = min(options)
            maintenance = laws[selected].value(size)
            old_mode, mode = mode, selected
            size += 1
            grew = True
        else:
            maintenance = laws[mode].value(size)
            if maintenance > quota or maintenance > stock:
                status = "maintenance_not_affordable"
                break
            cost, setup, changed, old_mode, grew = maintenance, F(0), False, mode, False
            status = "growth_blocked_under_declared_policy"
        stock -= cost
        events.append({"step": step+1, "size_before": before, "size_after": size,
                       "mode_before": old_mode, "mode_after": mode, "switched": changed,
                       "grew": grew, "quota": str(quota), "maintenance": str(maintenance),
                       "switch_debit": str(setup), "growth_debit": str(growth_cost if grew else 0),
                       "total_debit": str(cost), "remaining_credits": str(stock)})
        if not grew:
            break
    if size == capacity:
        status = "internal_candidate_slots_exhausted"
    return {"initial_active_size": 1, "final_active_size": size,
            "initial_credits": str(initial_stock), "remaining_credits": str(stock),
            "spent_credits": str(initial_stock-stock), "internal_slot_limit": capacity,
            "final_mode": mode, "status": status,
            "switch_sizes": [row["size_before"] for row in events if row["switched"]],
            "events": events}


def legacy_reference_rows():
    """Old counts vs repeated task number T, not number of agents n."""
    span = maximum_age(F(1, 10**6))+1
    rows = []
    for trials in (1, 14, 28, 1000):
        reused = resources(trials)
        refreshed = refresh_resources(trials, span)
        rows.append({"task_trials_T": trials,
                     "noiseless_alignment_to_task_read_ratio": fraction(F(reused["alignment_reads_once"], reused["target_reads"])),
                     "noisy_alignment_to_task_read_ratio": fraction(F(5*refreshed["alignment_blocks"], 20*trials)),
                     "alignment_blocks": refreshed["alignment_blocks"]})
    return {"axis_is_fixed_group_task_repetitions_not_system_size": True,
            "ratio_counts_readouts_only_not_total_cost": True,
            "old_declared_noise": "1/1000000", "certified_tasks_per_block": span,
            "asymptotic_alignment_to_task_read_ratio": fraction(F(1, 4*span)),
            "rows": rows}


def report():
    fixed, linear, quadratic = CostLaw(F(1), 0), CostLaw(F(1), 1), CostLaw(F(1), 2)
    scenarios = {
        "fixed_maintenance": simulate({"A": fixed}),
        "linear_maintenance": simulate({"A": linear}),
        "quadratic_maintenance": simulate({"A": quadratic}),
        "quadratic_with_linear_alternative": simulate({"A": quadratic, "B": linear}),
        "identical_linear_alternatives": simulate({"A": linear, "B": linear}),
        "quadratic_with_expensive_switch": simulate({"A": quadratic, "B": linear}, switch_cost=100),
        "finite_stock_exhaustion": simulate({"A": linear}, total_credits=8),
    }
    return {"round": 206,
            "model_kind": "Finite rational accounting sensitivity model; not quantum evolution or an empirical experiment",
            "declared_inputs": {"per_epoch_processing_quota": "4*n", "initial_consumable_work_credits": 1000,
                "total_internal_capability_slots": 16, "growth_cost": 2, "switch_cost": 6,
                "old_capability_preservation_and_new_capability_are_abstract_contracts": True,
                "physical_realizability_of_all_architecture_contracts_proved": False,
                "costs_derived_from_quantum_theory_or_Marx": False,
                "controller_clock_and_log_physical_costs_derived": False,
                "unit_is_joules_or_landauer_bits": False,
                "policy": "grow if one-epoch growth is affordable; minimize current debit, prefer current mode on ties; no staged switching"},
            "maintenance_share_formula": "M/B=(a/b)*n**(alpha-beta)",
            "maintenance_share_increases_iff_for_positive_power_laws": "alpha>beta",
            "ratio_table": [{"active_size_n": n,
                             "fixed_over_quota": fraction(F(1, 4*n)),
                             "linear_over_quota": fraction(F(1, 4)),
                             "quadratic_over_quota": fraction(F(n, 4))} for n in (1, 2, 4, 8, 16)],
            "scenarios": scenarios, "legacy_reference_read_account": legacy_reference_rows(),
            "conclusions": {"cost_pressure_alone_identifies_complex_structure": False,
                "spontaneous_J_or_Y_emergence_simulated": False,
                "generic_real_theory_instability_proved": False,
                "matched_models_must_charge_complex_resources_too": True,
                "budget_threshold_is_a_proved_physical_phase_transition": False,
                "finite_stock_can_limit_any_mode": True}}


class GrowthMaintenanceAuditTests(unittest.TestCase):
    def test_three_scalings_have_different_maintenance_share_trends(self):
        shares = [[law.value(n)/F(4*n) for n in (1, 2, 4, 8)]
                  for law in (CostLaw(F(1), 0), CostLaw(F(1), 1), CostLaw(F(1), 2))]
        self.assertTrue(all(a>b for a,b in zip(shares[0], shares[0][1:])))
        self.assertEqual(set(shares[1]), {F(1, 4)})
        self.assertTrue(all(a<b for a,b in zip(shares[2], shares[2][1:])))

    def test_quadratic_budget_boundary_blocks_growth_at_four(self):
        result = simulate({"A": CostLaw(F(1), 2)})
        self.assertEqual(result["final_active_size"], 4)
        self.assertEqual(result["status"], "growth_blocked_under_declared_policy")
        self.assertEqual(result["events"][-1]["maintenance"], "16")
        self.assertEqual(result["events"][-1]["quota"], "16")

    def test_linear_and_fixed_models_keep_growing_to_declared_inventory(self):
        for exponent, debit in ((0, 45), (1, 150)):
            result = simulate({"A": CostLaw(F(1), exponent)})
            self.assertEqual(result["final_active_size"], 16)
            self.assertEqual(F(result["spent_credits"]), debit)
            self.assertEqual(result["status"], "internal_candidate_slots_exhausted")

    def test_paid_switch_relieves_a_bottleneck_without_a_theory_label(self):
        result = simulate({"A": CostLaw(F(1), 2), "B": CostLaw(F(1), 1)})
        self.assertEqual(result["switch_sizes"], [4])
        self.assertEqual(result["final_active_size"], 16)
        self.assertEqual(F(result["spent_credits"]), 164)

    def test_switch_cost_is_not_omitted(self):
        result = simulate({"A": CostLaw(F(1), 2), "B": CostLaw(F(1), 1)}, switch_cost=100)
        self.assertEqual(result["switch_sizes"], [])
        self.assertEqual(result["final_active_size"], 4)

    def test_equal_cost_labels_cannot_create_a_transition(self):
        law = CostLaw(F(1), 1)
        first = simulate({"A": law, "B": law})
        second = simulate({"renamed": law, "alternative": law}, initial_mode="renamed")
        self.assertEqual(first["switch_sizes"], [])
        self.assertEqual(first["spent_credits"], second["spent_credits"])
        self.assertEqual(first["final_active_size"], second["final_active_size"])

    def test_finite_stock_limits_even_linear_growth(self):
        result = simulate({"A": CostLaw(F(1), 1)}, total_credits=8)
        self.assertEqual(result["final_active_size"], 3)
        self.assertEqual(result["status"], "maintenance_not_affordable")
        self.assertEqual(result["remaining_credits"], "1")

    def test_conservation_at_every_recorded_step(self):
        for exponent in (0, 1, 2):
            result = simulate({"A": CostLaw(F(1), exponent), "B": CostLaw(F(1), 1)})
            spent = F(0)
            for row in result["events"]:
                cost = sum(F(row[key]) for key in ("maintenance", "switch_debit", "growth_debit"))
                self.assertEqual(cost, F(row["total_debit"]))
                self.assertLessEqual(cost, F(row["quota"]))
                spent += cost
                self.assertEqual(spent+F(row["remaining_credits"]), F(result["initial_credits"]))
                self.assertEqual(row["size_after"]-row["size_before"], int(row["grew"]))

    def test_initially_full_inventory_does_not_charge_fake_activity(self):
        result = simulate({"A": CostLaw(F(1), 1)}, capacity=1)
        self.assertEqual(result["events"], [])
        self.assertEqual(result["remaining_credits"], result["initial_credits"])

    def test_old_refresh_budget_ratio_reuses_certified_fourteen_task_span(self):
        result = legacy_reference_rows()
        self.assertEqual(result["certified_tasks_per_block"], 14)
        self.assertEqual(result["asymptotic_alignment_to_task_read_ratio"]["exact"], "1/56")
        row = result["rows"][-1]
        self.assertEqual(row["alignment_blocks"], 72)
        self.assertEqual(row["noiseless_alignment_to_task_read_ratio"]["exact"], "1/4000")
        self.assertEqual(row["noisy_alignment_to_task_read_ratio"]["exact"], "9/500")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    checks = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(GrowthMaintenanceAuditTests))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result = report()
    result["automated_checks"] = {"run": checks.testsRun, "failures": 0, "errors": 0}
    if args.write_results:
        Path(__file__).with_name("growth_maintenance_audit_results.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"round": 206, "scenarios": {name: {k: v for k, v in r.items() if k != "events"}
                     for name, r in result["scenarios"].items()}}, indent=2))


if __name__ == "__main__":
    main()
