"""Exact finite-coupling certificate for the unnumbered quadrupole parent.

This is not a field simulation or an actual detector certificate. Fractions
certify smearing/field norm bounds and the vacuum prepared-state Dyson budget.
Default checks the saved record without writing; --write is exclusive.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_integral(power, denominator_power):
    return sum((F((-1)**j*comb(denominator_power, j), power+2*j+1)
                for j in range(denominator_power+1)), F())


def qr(value):
    return {"exact": str(value), "display_float": float(value)}


def smearing_bounds():
    radius, k0 = F(1, 4), F(4)
    radial_norm = polynomial_integral(2, 4)
    i0_pi = polynomial_integral(2, 8)/(4*radius**3*radial_norm**2)
    i2_pi = 64*polynomial_integral(4, 6)/(4*radius**5*radial_norm**2)
    assert i0_pi == F(776160, 4199)
    assert i2_pi == F(10644480, 221)
    # pi>3 supplies rational upper bounds; no numerical quadrature.
    i0_upper, i2_upper = i0_pi/3, i2_pi/3
    im2_upper = k0/18+i0_upper/k0**2
    field_norm_squared = 2*(k0**2/72+i0_upper/(2*k0))
    l_squared_upper = 72**2*field_norm_squared
    s_squared_upper = 72**2*im2_upper
    # Full total-source self + cross + mass terms, uniform |xi|<=2.
    w_upper = 2*(49*im2_upper+194*i0_upper+i2_upper)
    assert l_squared_upper < 300**2
    assert s_squared_upper < 150**2
    assert w_upper < 60000
    return {
        "radius": str(radius), "split_k": str(k0),
        "pi_times_I0": str(i0_pi), "pi_times_I2": str(i2_pi),
        "I_minus2_upper": qr(im2_upper),
        "L_squared_upper": qr(l_squared_upper),
        "S_squared_upper": qr(s_squared_upper),
        "W_norm_upper": qr(w_upper),
        "declared_rational_bounds": {"L": 300, "S": 150, "w": 60000},
        "mass_monopole_and_self_terms_included": True,
        "momentum_or_fock_cutoff": False,
    }


def coefficient():
    # sin^6(t) is the middle head population up to 5/16.
    # Integral sin^(2n)(t) from 0 to pi is pi*choose(2n,n)/4^n.
    c = F(5, 16)
    ipi = c*c*F(comb(12, 6), 4**6)
    dpi = 36*c*c*(F(comb(10, 5), 4**5)-F(comb(12, 6), 4**6))
    assert ipi == F(5775, 262144)
    assert dpi == F(4725, 65536)
    R, r, m2 = F(4), F(1, 4), F(3, 104)
    ic = (6*ipi-(3*R**2+m2)*dpi)/(4*R**5)
    assert ic == -F(22700475, 27917287424)
    assert F(22, 7) < R-2*r
    # Matrix-free Pauli check: [Z,(I+Y)/2]=-iX, hence i[Z,F]=X.
    # A +/-Z and B +X yield the factor -2 in H_inst=-g^2 xi C_AB.
    return {
        "population_prefactor": str(c),
        "integral_a_squared_div_pi": str(ipi),
        "integral_adot_squared_div_pi": str(dpi),
        "I_C": qr(ic),
        "setting_difference_leading_coefficient": qr(-2*ic),
        "sign_convention": "Delta p=p(A+)-p(A-); coefficient multiplies g^2*(xi-1)",
        "readout": "all-head all-field nonselective B Y effect",
        "whole_operator_replaced_by_mean": False,
    }, ic


def dyson_budget(ic):
    g, pi_upper = F(1, 10**15), F(22, 7)
    L, S, w = F(300), F(150), F(60000)
    u, v = 2*g*L*pi_upper, g*g*w*pi_upper
    assert 0 < u < 1 and 0 < v < 1
    d = u*u+v
    r3 = u**3/(1-u)+u*v/((1-u)*(1-v))+v*v/(1-v)
    probability_tail = 2*r3*(1+u+d)+r3*r3+2*u*d+d*d
    pair_tail = 2*probability_tail
    # Uniform outside |xi-1|<1/2, within the fixed xi in [0,2] menu.
    leading = g*g*abs(ic)
    lower = leading-pair_tail
    assert lower > 0
    simple_lower = F(7, 10**34)
    assert lower > simple_lower
    positivity_margin = 2-g*g*(2*S*S+w)
    assert positivity_margin > 1
    # A required future matching budget, NOT a verified implementation.
    hypothetical_each_match = F(1, 10**35)
    conditional_lower = lower-2*hypothetical_each_match
    assert conditional_lower > simple_lower
    return {
        "g": str(g), "T": "pi", "T_upper": str(pi_upper),
        "xi_range": [0, 2], "deviation_test": "abs(xi-1)>=1/2",
        "u": qr(u), "v": qr(v), "second_vector_bound_d": qr(d),
        "vector_tail_r3": qr(r3),
        "single_probability_tail": qr(probability_tail),
        "two_setting_probability_tail": qr(pair_tail),
        "minimal_leading_difference": qr(leading),
        "normal_parent_difference_lower": qr(lower),
        "simple_strict_lower": str(simple_lower),
        "hamiltonian_positivity_margin": qr(positivity_margin),
        "hypothetical_actual_probability_error_each": str(hypothetical_each_match),
        "conditional_actual_difference_lower": qr(conditional_lower),
        "actual_matching_certified": False,
        "uniform_unknown_matter_and_reference": True,
        "field_preparation": "normal vacuum; full Fock output retained",
        "zero_source_head_boundary_scope": "prepared head image only, not every head state",
    }


def calculate():
    coeff, ic = coefficient()
    sources = {
        "archive_935_955/research_note_940.md":
            sha(BASE/"archive_935_955/research_note_940.md"),
        "archive_1161_/1171/proof.md": sha(BASE/"archive_1161_/1171/proof.md"),
        "archive_1161_/_admission/after1185_quadrupole_screen/proof.md":
            sha(HERE.parent/"after1185_quadrupole_screen/proof.md"),
    }
    return {
        "schema": "unnumbered_dynamic_quadrupole_normal_parent_v1",
        "status": "PASS_MATHEMATICAL_PARENT_WITH_SCOPE",
        "numbered_round": None, "new_scientific_groups": 0,
        "check_groups": 3,
        "smearing_and_complete_norm": smearing_bounds(),
        "coefficient": coeff,
        "finite_coupling_budget": dyson_budget(ic),
        "source_sha256": sources,
        "proof_sha256": sha(HERE/"response_proof.md"),
        "code_sha256": sha(Path(__file__)),
        "not_certified": [
            "Dynamic COM, body frame, support and full recoil matching.",
            "Finite fresh setting operation and actual terminal instrument accuracy.",
            "All-order Ward identity, nonlinear Einstein dynamics or physical xi selection.",
            "Feasibility of detecting a probability difference of order 1e-34.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = calculate()
    target = HERE/"response_results.json"
    if args.write:
        with target.open("x", encoding="utf-8", newline="\n") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        mode = "WRITE_ONCE_PASS"
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        if saved != result:
            raise AssertionError("Saved response certificate or bound sources changed")
        mode = "READ_ONLY_MATCH"
    print(mode+"; unnumbered; groups=3; actual_matching=False")
    print("parent difference lower =", result["finite_coupling_budget"]
          ["normal_parent_difference_lower"]["display_float"])


if __name__ == "__main__":
    main()
