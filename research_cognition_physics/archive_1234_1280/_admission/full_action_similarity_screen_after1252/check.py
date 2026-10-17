"""Exact algebra replay of the non-numbered full-action similarity screen.

Only rational Laurent-polynomial identities are tested.  No ODE, quantum
instrument, field-theory existence claim or additional scientific round.
Default mode checks saved output without writing; --write is exclusive.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
SCREEN = HERE.parent / "full_action_similarity_screen_after1252.md"
NAMES = ("s", "J", "q", "r", "B", "p", "X", "Y", "K", "I",
         "W", "f", "Wp", "fp", "u", "v", "V", "H", "kappa")
INDEX = {name: i for i, name in enumerate(NAMES)}
ZERO = (0,) * len(NAMES)


class Poly:
    """Small sparse Laurent algebra; negative exponents used only for s."""

    def __init__(self, value=0):
        if isinstance(value, Poly):
            self.terms = dict(value.terms)
        elif isinstance(value, dict):
            self.terms = {m: F(c) for m, c in value.items() if c}
        else:
            self.terms = {} if not value else {ZERO: F(value)}

    @classmethod
    def var(cls, name):
        m = list(ZERO)
        m[INDEX[name]] = 1
        return cls({tuple(m): F(1)})

    def __add__(self, other):
        out = dict(self.terms)
        for m, c in Poly(other).terms.items():
            out[m] = out.get(m, F(0)) + c
        return Poly(out)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + (-Poly(other))

    def __rsub__(self, other):
        return Poly(other) - self

    def __mul__(self, other):
        out = {}
        for a, ca in self.terms.items():
            for b, cb in Poly(other).terms.items():
                m = tuple(x + y for x, y in zip(a, b))
                out[m] = out.get(m, F(0)) + ca * cb
        return Poly(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * (Poly(other) ** -1)

    def __rtruediv__(self, other):
        return Poly(other) * (self ** -1)

    def __pow__(self, power):
        assert isinstance(power, int)
        if power < 0:
            assert len(self.terms) == 1, "only monomial inverses are used"
            m, c = next(iter(self.terms.items()))
            return Poly({tuple(power * i for i in m): c ** power})
        out = Poly(1)
        for _ in range(power):
            out = out * self
        return out

    def diff(self, name):
        index = INDEX[name]
        out = {}
        for m, c in self.terms.items():
            if m[index]:
                new = list(m)
                new[index] -= 1
                out[tuple(new)] = c * m[index]
        return Poly(out)

    def scaled(self, weights):
        """Substitute x_s=s**weight*x, retaining the explicit physical s."""
        out = {}
        for m, c in self.terms.items():
            new = list(m)
            new[INDEX["s"]] += sum(m[INDEX[k]] * w for k, w in weights.items())
            new = tuple(new)
            out[new] = out.get(new, F(0)) + c
        return Poly(out)

    def at(self, values):
        total = F(0)
        for m, c in self.terms.items():
            term = c
            for name, power in zip(NAMES, m):
                if power:
                    term *= F(values[name]) ** power
            total += term
        return total

    def __eq__(self, other):
        return self.terms == Poly(other).terms


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    sy = {name: Poly.var(name) for name in NAMES}
    s, J, q, r, B, p, X, Y, K, I = [sy[n] for n in NAMES[:10]]
    W, f, Wp, fp, u, v, V, H, kappa = [sy[n] for n in NAMES[10:]]
    identities = []

    def same(name, left, right):
        assert Poly(left) == Poly(right), name
        identities.append(name)

    # The exponent ledger checks every canonical pair, not just mass density.
    pair_weights = {
        "clock_J_dtheta": [1, 1],
        "oscillator_r_dq": [2, 0],
        "program_B_dpsi": [2, 0],
        "memory_K_dY": [2, 0],
        "recorder_I_dvartheta": [1, 1],
        "moving_worldline_P_dx": [1, 1],
    }
    assert all(sum(pair) == 2 for pair in pair_weights.values())
    action_weights = {
        "EH_fixed_kappa_Lambda_zero": 4 - 2,
        "p_kinetic": 4 + 2 * -1,
        "X_kinetic": 4 + 2 * -1,
        "p_mass_potential": 4 + 2 * -1,
        "X_mass_potential": 4 + 2 * -1,
        "all_three_rods_fixed_target_period": 4 + 2 * -1,
        "clock_mass_fixed_label_measure": 0 + 1 + 1,
        "recorder_mass_fixed_label_measure": 0 + 1 + 1,
    }
    assert set(action_weights.values()) == {2}

    Mc = F(71, 100) + J + (q*q + r*r + K*K)/2 \
        + (B-1)*W*q*p/4 + f*X*Y/10
    Mr = F(1, 10) + (1+p)*I/10
    # These are functions of the physical scaled canonical coordinates.
    Mcs = 71*s/100 + J + s*q*q/2 + (r*r + K*K)/(2*s**3) \
        + (B-s*s)*W*q*p/(4*s) + s*f*X*Y/10
    Mrs = s/10 + (1+p)*I/10
    weights = {"J": 1, "r": 2, "B": 2, "K": 2, "I": 1}
    same("full_clock_mass", Mcs.scaled(weights), s*Mc)
    same("full_recorder_mass", Mrs.scaled(weights), s*Mr)
    same("clock_p_force", Mcs.diff("p").scaled(weights), s*Mc.diff("p"))
    same("clock_X_force", Mcs.diff("X").scaled(weights), s*Mc.diff("X"))

    def internal_flow(mass, gate_scale):
        return {
            "J": -(mass.diff("W")*Wp + mass.diff("f")*fp)/gate_scale,
            "q": mass.diff("r"), "r": -mass.diff("q"),
            "B": Poly(0), "Y": mass.diff("K"), "K": -mass.diff("Y"),
        }

    old_flow = internal_flow(Mc, Poly(1))
    new_flow = internal_flow(Mcs, s)
    variable_weights = {"J": 1, "q": 0, "r": 2, "B": 2, "Y": 0, "K": 2}
    for name, old_rhs in old_flow.items():
        same("canonical_flow_" + name, new_flow[name].scaled(weights),
             s**(variable_weights[name]-1)*old_rhs)
    same("clock_proper_speed", Mcs.diff("J").scaled(weights), Mc.diff("J"))
    same("program_phase_no_extra_B_minus_one", Mcs.diff("B").scaled(weights),
         Mc.diff("B")/s)
    same("recorder_phase", Mrs.diff("I").scaled(weights), Mr.diff("I"))

    def chain_derivative(mass, flow, gate_scale):
        return sum((mass.diff(name)*rhs for name, rhs in flow.items()), Poly(0)) \
            + (mass.diff("W")*Wp + mass.diff("f")*fp)/gate_scale \
            + mass.diff("p")*v/gate_scale + mass.diff("X")*V/gate_scale

    transferred_power = Mc.diff("p")*v + Mc.diff("X")*V
    same("old_complete_internal_power_cancellation",
         chain_derivative(Mc, old_flow, Poly(1)), transferred_power)
    same("scaled_complete_internal_power_cancellation",
         chain_derivative(Mcs, new_flow, s).scaled(weights), transferred_power)
    same("recorder_exchange", (Mrs.diff("p")*v/s).scaled(weights),
         Mr.diff("p")*v)
    # Omitting the clock control work leaves a nonzero symbolic source.
    without_clock_work = dict(new_flow)
    without_clock_work["J"] = Poly(0)
    deleted_power = chain_derivative(Mcs, without_clock_work, s).scaled(weights) \
        - transferred_power
    same("omitted_clock_work_defect", deleted_power,
         Mc.diff("W")*Wp + Mc.diff("f")*fp)
    assert deleted_power != Poly(0)

    Fp = 1 + p/5
    rho = (v*v + p*p + V*V + X*X)/2 + 3*Fp*u*u/20 \
        + (Mc+Mr/10)*u**3
    pressure = (v*v-p*p+V*V-X*X)/2 - Fp*u*u/20
    rho_s = (v*v + p*p + V*V + X*X)/(2*s*s) + 3*Fp*u*u/(20*s*s) \
        + (Mcs.scaled(weights)+Mrs.scaled(weights)/10)*u**3/s**3
    pressure_s = (v*v-p*p+V*V-X*X)/(2*s*s) - Fp*u*u/(20*s*s)
    same("full_source_density", rho_s, rho/s**2)
    same("full_source_pressure", pressure_s, pressure/s**2)
    same("Friedmann_constraint_same_kappa", H*H/s**2-kappa*rho_s/3,
         (H*H-kappa*rho/3)/s**2)
    same("null_source_same_kappa", kappa*(rho_s+pressure_s),
         kappa*(rho+pressure)/s**2)
    source_weights = {
        "scalar_kinetic": -2,
        "scalar_mass": -2,
        "rods": -2,
        "clock_density_times_mass": -3+1,
        "recorder_density_times_mass": -3+1,
        "divergence_and_exchange": -2-1,
        "total_material_inventory": -2+3,
        "action_over_fixed_hbar": 2,
    }
    assert source_weights["clock_density_times_mass"] == -2
    assert source_weights["recorder_density_times_mass"] == -2
    assert source_weights["divergence_and_exchange"] == -3
    assert source_weights["total_material_inventory"] == 1

    # Initial-value decomposition is read from the stated mass, not fitted.
    initial = {name: F(0) for name in NAMES}
    initial.update(s=F(1), J=F(1, 10), q=F(1, 2), B=F(1),
                   X=F(1, 10), I=F(1), u=F(1))
    mc0, mr0 = Mc.at(initial), Mr.at(initial)
    dust0 = mc0+mr0/10
    field_rod0 = (rho-(Mc+Mr/10)*u**3).at(initial)
    assert (mc0, mr0, dust0, field_rod0) == \
        (F(187, 200), F(1, 5), F(191, 200), F(31, 200))
    assert rho.at(initial) == F(111, 100)
    same("initial_full_inventory", 1000*s**3*(F(111, 100)/s**2), 1110*s)
    bad_rho = field_rod0/s**2+dust0/s**3
    bad_constraint = kappa*F(111, 100)/(3*s**2)-kappa*bad_rho/3
    same("deleted_mass_scaling_constraint", bad_constraint,
         kappa*dust0*(s**-2-s**-3)/3)
    deletion = bad_constraint.at({"s": F(2), "kappa": F(3)})
    assert deletion == F(191, 1600) > 0

    return {
        "status": "PASS",
        "scope": "NON_NUMBERED_EXACT_ALGEBRA_SCREEN",
        "formal_scientific_groups": 0,
        "arithmetic": "fractions.Fraction; sparse Laurent polynomials; no ODE",
        "bindings": {"screen_sha256": sha(SCREEN), "check_sha256": sha(Path(__file__))},
        "symplectic_pair_weights": pair_weights,
        "action_scaling_exponents": action_weights,
        "source_scaling_exponents": source_weights,
        "exact_identities": identities,
        "initial_data": {
            "clock_mass": str(mc0), "recorder_mass": str(mr0),
            "field_and_rods_density": str(field_rod0), "all_dust_density": str(dust0),
            "full_density": "111/100", "volume": "1000*s^3",
            "total_material_inventory": "1110*s",
        },
        "deletion": {
            "rule": "keep mapped geometry and density, omit each label mass scaling",
            "formula": "kappa/3 * 191/200 * (s^-2-s^-3)",
            "s": "2", "kappa": "3", "initial_constraint_defect": str(deletion),
        },
        "limitations": [
            "Exact algebra checks the proposed classical dictionary, not new dynamics.",
            "Coordinate and label domains, Hilbert variation and finite phase lifts remain in the screen.",
            "No local quantum model, actual terminal precision or finite-body realization is certified.",
            "No numerical integration, parameter scan or new formal scientific group.",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    target = HERE / "results.json"
    if args.write:
        with target.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
        print("INITIAL_SAVE")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert result == saved, "saved algebra certificate or bindings changed"
        print("READ_ONLY_MATCH")
    print("PASS: exact classical similarity algebra; formal_scientific_groups=0")


if __name__ == "__main__":
    main()
