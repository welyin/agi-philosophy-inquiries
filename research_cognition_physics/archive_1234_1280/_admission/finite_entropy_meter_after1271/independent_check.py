"""Independent finite entropy-meter audit; no author-code imports.

One independent mathematical group. Fraction bounds prove the continuous
guardrails; split Gauss quadrature is a numerical diagnostic, not a certified
physical implementation or a uniform quadrature-error theorem.
Default: read-only result comparison. --write: exclusive first creation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents
            if (p / "research_cognition_physics").is_dir()
            and (p / "scripts").is_dir())
OUT = HERE / "independent_results.json"


def frac(x):
    return {"numerator": x.numerator, "denominator": x.denominator}


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return out


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def derivative(a):
    return [F(i)*a[i] for i in range(1, len(a))]


def evaluate(a, x):
    y = F(0)
    for coefficient in reversed(a):
        y = y*x + coefficient
    return y


def rational_checks():
    f = list(map(F, [0, 0, 0, 10, -15, 6]))
    fp, fpp = derivative(f), derivative(derivative(f))
    fp_factor = mul([F(0), F(0), F(30)], [F(1), F(-2), F(1)])
    assert fp == fp_factor
    # fp = 30 [s(1-s)]^2, and s(1-s) <= 1/4 on [0,1].
    assert evaluate(fp, F(1, 2)) == F(15, 8)
    assert [evaluate(f, x) for x in [F(0), F(1)]] == [F(0), F(1)]
    assert all(evaluate(poly, x) == 0
               for poly in (fp, fpp) for x in [F(0), F(1)])
    # fpp = 60s(1-s)(1-2s), hence |q| has one peak.
    assert fpp == mul([F(0), F(60)], mul([F(1), F(-1)],
                                                     [F(1), F(-2)]))

    d, tau = F(1, 10), F(1000)
    mlo, mhi = F(10), F(101, 10)
    nu = F(15, 8)*d/tau
    tlo, thi = F(49, 500), F(51, 500)
    assert nu == F(3, 16000)
    assert 1/mhi - nu > tlo and 1/mlo + nu < thi
    zi_f, zi_r = F(1000)/mlo**3, F(1000)/mhi**3
    z_exc = F(3, 4)*d/tlo
    zlo, zhi = min(zi_f, zi_r)-z_exc, max(zi_f, zi_r)+z_exc
    assert zlo > F(1, 5) and zhi < F(9, 5)
    # Actual off-equilibrium C_gamma and the equilibrium test at same M,X
    # are distinct. This second bound addresses the latter.
    c_gamma_actual_upper = 4*F(9, 5)
    c_gamma_equilibrium_upper = 4*F(9, 5)*((1/mlo)/tlo)**3
    assert c_gamma_actual_upper < 8
    assert c_gamma_equilibrium_upper < 8 < mlo**2

    entropy_target = (mhi*mhi-mlo*mlo)/2
    irreversibility = nu*d/tlo**2
    width = 2*irreversibility
    tv_t = d/(mlo*mhi)+2*nu
    tv_work = d/4+F(9, 5)*tv_t
    assert entropy_target == F(201, 200)
    assert irreversibility == F(75, 38416)
    assert width == F(75, 19208) < F(1, 250)
    assert tv_work < F(7, 250) < 1

    # Independent exact derivation at equilibrium, fixed total E and X:
    # d(Tgas - 1/M)/dM = -T/[4(E-M)] + T/M.
    # Multiplication by positive 4M(E-M)/T leaves 4E - 5M.
    # Check as polynomial, not by sample fitting.
    thermal_polynomial = add([F(0), F(-1)],
                             [F(4), F(-4)])  # E=1, normalized M=x
    assert thermal_polynomial == [F(4), F(-5)]
    branch_examples = {}
    for name, total, mass in [("left_unstable", F(5), F(3)),
                              ("right_stable", F(5), F(9, 2))]:
        e_rad, temp = total-mass, 1/mass
        response = -temp/(4*e_rad)+temp/mass
        expected = temp*(4*total-5*mass)/(4*mass*(total-mass))
        assert response == expected
        heat_capacity_gamma, heat_capacity_abs_bh = 4*e_rad/temp, mass**2
        stable = response < 0
        assert stable == (mass > 4*total/5)
        assert stable == (heat_capacity_gamma < heat_capacity_abs_bh)
        assert stable == (name == "right_stable")
        branch_examples[name] = {
            "total_energy": frac(total), "mass": frac(mass),
            "gas_energy": frac(e_rad), "response_derivative": frac(response),
            "heat_capacity_gamma": frac(heat_capacity_gamma),
            "heat_capacity_abs_bh": frac(heat_capacity_abs_bh),
            "stable": stable,
        }

    # Independent coefficient check of energy/work and entropy identities.
    # Coefficients of (q, Z*Tdot): E_gamma_dot=(-3/4,1),
    # p*Vdot=(-1/4,-1), Mdot=(1,0), hence complete sum zero.
    e_dot, work_dot, m_dot = (F(-3, 4), F(1)), (F(-1, 4), F(-1)), (F(1), F(0))
    assert all(e_dot[i]+work_dot[i]+m_dot[i] == 0 for i in range(2))
    assert F(4, 3)*F(-3, 4) == -1
    return {
        "quintic_endpoint_orders": {"q_zero": True, "qdot_zero": True,
                                    "maximum_fprime": frac(F(15, 8))},
        "nu": frac(nu), "temperature_lower": frac(tlo),
        "temperature_upper": frac(thi),
        "proved_z_lower": frac(zlo), "proved_z_upper": frac(zhi),
        "declared_z_interval": [frac(F(1, 5)), frac(F(9, 5))],
        "equilibrium_Cgamma_upper": frac(c_gamma_equilibrium_upper),
        "Cgamma_current_upper": frac(c_gamma_actual_upper),
        "abs_Cbh_lower": frac(mlo**2),
        "delta_wald": frac(entropy_target),
        "one_leg_entropy_error_upper": frac(irreversibility),
        "two_leg_interval_width_upper": frac(width),
        "temperature_total_variation_upper": frac(tv_t),
        "controller_work_total_variation_upper": frac(tv_work),
        "initial_controller_balance": frac(F(1)),
        "controller_balance_lower": frac(1-tv_work),
        "energy_coefficient_cancellation": True,
        "radiation_entropy_identity": True,
        "fixed_energy_branch_examples": branch_examples,
        "right_branch_threshold_mass_fraction": frac(F(4, 5)),
    }


def profile(s, initial, final):
    tau = 1000.0
    delta = final-initial
    f = s**3*(10+s*(-15+6*s))
    fp = 30*s*s*(1-s)**2
    fpp = 60*s*(1-s)*(1-2*s)
    mass = initial+delta*f
    q = delta*fp/tau
    qdot = delta*fpp/(tau*tau)
    th = 1/mass
    temp = th+q
    tdot = -q/mass**2+qdot
    return mass, q, qdot, th, temp, tdot


def path(initial, final, panels):
    # 12-point Gauss on every panel. Prefix integrals at the Gauss nodes
    # are independently computed on [panel_left,node].
    z, w = np.polynomial.legendre.leggauss(12)
    unit, weight = (z+1)/2, w/2
    left = np.arange(panels, dtype=float)/panels
    width = 1.0/panels
    grid = left[:, None]+width*unit[None, :]
    mass, q, qdot, th, temp, tdot = profile(grid, initial, final)
    integrand = 1000*q/temp
    panel_int = width*(integrand@weight)
    prefix = np.concatenate(([0.0], np.cumsum(panel_int)[:-1]))
    # Partial integral for each evaluation node.
    inner = left[:, None, None]+width*unit[None, :, None]*unit[None, None, :]
    _, iq, _, _, it, _ = profile(inner, initial, final)
    partial = width*unit[None, :]*np.sum((1000*iq/it)*weight, axis=2)
    cumulative = prefix[:, None]+partial
    z0 = 1000/initial**3
    zz = z0-.75*cumulative
    zdot = -.75*q/temp
    xx = zz/temp**3
    xdot = zdot/temp**3-3*zz*tdot/temp**4

    def integrate(values):
        return float(1000*width*np.sum(values*weight[None, :]))

    # Compute from E=X*T^4 and p*Vdot=T^4 Xdot/3, not from the
    # already simplified source coefficients.
    edot = xdot*temp**4+4*xx*temp**3*tdot
    piston_power = temp**4*xdot/3
    controller_dot = -q/4-zz*tdot
    entropy_dot = F(4, 3).__float__()*(xdot*temp**3+3*xx*temp**2*tdot)
    entropy_prod = q*q/(temp*th)
    target_signed = (final**2-initial**2)/2
    integral = float(np.sum(panel_int))
    zfinal = z0-.75*integral
    egamma_initial, egamma_final = z0/initial, zfinal/final
    control_change = integrate(controller_dot)
    heat_residual = q+edot+piston_power
    total_residual = q+edot+controller_dot
    entropy_residual = entropy_dot+q/temp
    endpoint_total = final-initial+egamma_final-egamma_initial+control_change
    result = {
        "initial_mass": initial, "final_mass": final, "panels": panels,
        "q_over_T_integral": integral,
        "irreversibility_integral": integrate(entropy_prod),
        "wald_minus_meter_signed": target_signed-integral,
        "radiation_entropy_change": F(4, 3).__float__()*(zfinal-z0),
        "Z_final": zfinal, "X_final": zfinal*final**3,
        "gas_energy_change": egamma_final-egamma_initial,
        "controller_energy_change": control_change,
        "controller_work_total_variation": integrate(np.abs(controller_dot)),
        "temperature_total_variation": integrate(np.abs(tdot)),
        "temperature_min_sampled": float(np.min(temp)),
        "temperature_max_sampled": float(np.max(temp)),
        "Z_min_sampled": float(min(z0, zfinal, np.min(zz))),
        "Z_max_sampled": float(max(z0, zfinal, np.max(zz))),
        "equilibrium_capacity_ratio_max_sampled":
            float(np.max(4*zz*(th/temp)**3/mass**2)),
        "endpoint_total_energy_residual": endpoint_total,
        "max_energy_balance_residual": float(np.max(np.abs(total_residual))),
        "max_radiation_heat_residual": float(np.max(np.abs(heat_residual))),
        "max_control_work_residual": float(np.max(np.abs(piston_power-controller_dot))),
        "max_entropy_balance_residual": float(np.max(np.abs(entropy_residual))),
    }
    assert result["temperature_min_sampled"] > .098
    assert result["temperature_max_sampled"] < .102
    assert result["Z_min_sampled"] > .2 and result["Z_max_sampled"] < 1.8
    assert result["equilibrium_capacity_ratio_max_sampled"] < 1
    assert result["controller_work_total_variation"] < .028
    assert abs(result["irreversibility_integral"] -
               result["wald_minus_meter_signed"]) < 5e-13
    assert abs(endpoint_total) < 5e-13
    for name in ["max_energy_balance_residual", "max_radiation_heat_residual",
                 "max_control_work_residual", "max_entropy_balance_residual"]:
        assert result[name] < 2e-14
    # Endpoint profile q=qdot=0 implies Tdot=Xdot=0 independently of Z.
    for endpoint in [np.array(0.0), np.array(1.0)]:
        _, qe, qde, _, _, tde = profile(endpoint, initial, final)
        assert qe == qde == tde == 0
    return result


def main_payload():
    exact = rational_checks()
    branches = {}
    for name, start, finish in [("forward", 10., 10.1), ("reverse", 10.1, 10.)]:
        coarse, fine = path(start, finish, 128), path(start, finish, 256)
        differences = {}
        for key in ["q_over_T_integral", "irreversibility_integral",
                    "controller_energy_change", "controller_work_total_variation"]:
            differences[key] = abs(coarse[key]-fine[key])
            assert differences[key] < 2e-10
        branches[name] = {"diagnostic": fine, "refinement_differences": differences}
    jp = branches["forward"]["diagnostic"]["q_over_T_integral"]
    jm = -branches["reverse"]["diagnostic"]["q_over_T_integral"]
    target = 201/200
    bound = float(F(75, 19208))
    assert jp < target < jm and jm-jp < bound
    sources = {
        "independent_code": Path(__file__).resolve(),
        "author_proof_only": HERE / "proof.md",
        "old317": ROOT / "research_cognition_physics/archive_301_341/research_note_317.md",
        "old1250": ROOT / "research_cognition_physics/archive_1234_/1250/proof.md",
        "old1271": ROOT / "research_cognition_physics/archive_1234_/1271/proof.md",
    }
    return {
        "schema": "independent_finite_entropy_meter_v1",
        "status": "PASS_CONDITIONAL_SCALAR_ENTROPY_METER_ONLY",
        "independent_groups": 1,
        "author_code_imported": False,
        "exact_fraction_guardrails": exact,
        "numerical_diagnostics": branches,
        "numerical_finite_interval": {"J_plus": jp, "J_minus": jm,
                                      "delta_wald": target,
                                      "width": jm-jp,
                                      "width_analytic_upper": bound},
        "scope": [
            "Exact arithmetic certifies the stated scalar effective equations only.",
            "Gauss refinement is a numerical diagnostic, not a rigorous quadrature tail bound.",
            "No actual controller, finite POVM, wall stress, or covariant complete source is certified.",
            "Forward and reverse use distinct gas preparations; no free resetting.",
            "No claim that the off-equilibrium state itself is equilibrium.",
            "Right stable branch is derived at fixed total energy and volume.",
            "Operational target entropy squeeze is conditional on actual adiabatic admissibility and additive complete accounting.",
            "No Einstein generation or rejection of the 1271 flat model.",
        ],
        "source_bindings": {
            label: {"path": str(path.relative_to(ROOT)).replace("\\", "/"),
                    "sha256": digest(path)}
            for label, path in sources.items()
        },
    }


def compare(saved, actual, path="root"):
    if isinstance(actual, float):
        if not isinstance(saved, (int, float)) or not math.isclose(
                saved, actual, rel_tol=2e-10, abs_tol=2e-12):
            raise AssertionError((path, saved, actual))
    elif isinstance(actual, dict):
        assert saved.keys() == actual.keys(), path
        for key in actual:
            compare(saved[key], actual[key], path+"."+key)
    elif isinstance(actual, list):
        assert len(saved) == len(actual), path
        for i, item in enumerate(actual):
            compare(saved[i], item, path+"["+str(i)+"]")
    else:
        assert saved == actual, (path, saved, actual)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true",
                        help="create independent_results.json exclusively")
    args = parser.parse_args()
    actual = main_payload()
    if args.write:
        with OUT.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(actual, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        print("PASS_EXCLUSIVE_WRITE", OUT)
    else:
        compare(json.loads(OUT.read_text(encoding="utf-8")), actual)
        print("PASS_READ_ONLY_MATCH", OUT)
    print(json.dumps(actual["numerical_finite_interval"], sort_keys=True))
