"""Independent round-1222 certificate audit; no author modules are imported.

Default is read-only.  --write-results creates, but never overwrites, the JSON.
Analytic inequalities belong to the proof; Decimal and quadrature are numerical
cross-checks, not interval arithmetic or an implemented detector probability.
"""
from __future__ import annotations

import argparse
from decimal import Decimal as D, localcontext
import hashlib
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PI_TEXT = (
    "3.141592653589793238462643383279502884197169399375105820974944592307"
    "8164062862089986280348253421170679821480865132823066470938446"
)


def s(x):
    if isinstance(x, D):
        return format(x, ".25g")
    if isinstance(x, dict):
        return {k: s(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [s(v) for v in x]
    return x


def calculate(precision):
    with localcontext() as ctx:
        ctx.prec = precision
        pi = +D(PI_TEXT)
        beta, eps, eps_mat = 2*pi, D(".1"), D(".1")
        m, mass_source, coupling = D(100), D(500), D(".02")
        gaps, times, budgets = [D(".1"), D(".2")], [D(20), D(30)], [D("1e5"), D("1e8")]

        def f(e):
            return 1/beta if e == 0 else e/((beta*e).exp()-1)

        def bias(e, delta):
            return ((1+delta/f(e))/(1+delta/f(-e))).ln()

        def acosh(x):
            assert x >= 1
            return (x+(x*x-1).sqrt()).ln()

        # Exact spectral identities tested on both signs, and even discrete
        # spectral measures with a finite second moment.  This tests constants;
        # it is not a replacement for the all-H1 analytic argument.
        spectral = []
        weights = [D(".1"), D(".3"), D(".4"), D(".2")]
        omega = [D(0), D(".1"), D(1), D(4)]
        r = sum(w*x*x for w, x in zip(weights, omega))
        for e in [D("-.7"), D(0), *gaps, D(2)]:
            assert abs(f(-e)-f(e)-e) < D("1e-55")
            avg = sum(w*(f(e+x)+f(e-x))/2 for w, x in zip(weights, omega))
            avg_minus = sum(w*(f(-e+x)+f(-e-x))/2 for w, x in zip(weights, omega))
            delta = avg-f(e)
            assert -D("1e-55") <= delta <= beta*r/12
            assert abs(avg_minus-avg-e) < D("1e-55")
            spectral.append({"E": e, "delta": delta, "upper": beta*r/12})
        curvature = []
        for y in [D(0), D(".001"), D(".1"), D(1), D(5), D(20)]:
            ey = y.exp()
            value = beta/6 if y == 0 else beta*ey*((y-2)*ey+y+2)/(ey-1)**3
            assert 0 <= value <= beta/6
            curvature.append({"beta_E": y, "f_second": value})

        thresholds, entries = [], []
        rstars = []
        for gap in gaps:
            fp, fm = f(gap), f(-gap)
            assert 0 < eps < beta*gap
            denominator = fm-eps.exp()*fp
            assert denominator > 0
            kstar = fp*fm*(eps.exp()-1)/denominator
            rstar = 12*kstar/beta
            assert abs(bias(gap, kstar)-eps) < D("1e-55")
            assert bias(gap, kstar/2) < eps < bias(gap, 2*kstar)
            rstars.append(rstar)
            thresholds.append({
                "gap": gap, "f_plus": fp, "f_minus": fm,
                "K_star": kstar, "r_star": rstar,
                "minimum_T_H1_certificate": pi/rstar.sqrt(),
                "minimum_T_cos_squared_certificate": (4*pi*pi/(3*rstar)).sqrt(),
            })
            for t in times:
                norm = 3*t/8
                r_cos2 = 4*pi*pi/(3*t*t)
                k = beta*r_cos2/12
                qmin = coupling*coupling*norm/(2*pi)*fp
                material_allowance = ((eps_mat.exp()-1)/(eps_mat.exp()+1))*qmin
                worst_log = ((qmin+material_allowance)/(qmin-material_allowance)).ln()
                assert abs(worst_log-eps_mat) < D("1e-55")
                entries.append({
                    "gap": gap, "T": t, "window_norm_squared": norm,
                    "r_cos_squared": r_cos2, "K": k,
                    "log_bias_upper": bias(gap, k),
                    "thermal_certificate_pass": k <= kstar,
                    "leading_reference_q_lower": qmin,
                    "conditional_absolute_probability_error_allowance": material_allowance,
                    "worst_material_log_error": worst_log,
                })

        resource = []
        for budget in budgets:
            assert budget >= m+mass_source
            gamma = (budget*budget+m*m-mass_source*mass_source)/(2*budget*m)
            tmax = 2*acosh(gamma)
            pmax = m*(gamma*gamma-1).sqrt()
            recovered_b = (m*m+pmax*pmax).sqrt()+(mass_source*mass_source+pmax*pmax).sqrt()
            assert abs(recovered_b-budget) < D("1e-45")
            by_time = []
            for t in times:
                p = m*((t/2).exp()-(-t/2).exp())/2
                bmin = (m*m+p*p).sqrt()+(mass_source*mass_source+p*p).sqrt()
                assert (t <= tmax) == (bmin <= budget)
                by_time.append({"T": t, "B_required": bmin, "energy_pass": bmin <= budget})
            tmin_h1 = pi/min(rstars).sqrt()
            tmin_cos2 = (4*pi*pi/(3*min(rstars))).sqrt()
            resource.append({
                "B": budget, "T_max": tmax, "p_max": pmax,
                "fixed_times": by_time,
                "continuous_H1_certificate_interval": [tmin_h1, tmax] if tmin_h1 <= tmax else None,
                "continuous_cos_squared_certificate_interval": [tmin_cos2, tmax] if tmin_cos2 <= tmax else None,
            })

        matrix = []
        for item in resource:
            for test in item["fixed_times"]:
                same_t = [x for x in entries if x["T"] == test["T"]]
                thermal = all(x["thermal_certificate_pass"] for x in same_t)
                matrix.append({
                    "B": item["B"], "T": test["T"],
                    "both_gap_thermal_certificate_pass": thermal,
                    "two_body_energy_pass": test["energy_pass"],
                    "combined_fixed_plan_pass": thermal and test["energy_pass"],
                })
        assert not any(row["combined_fixed_plan_pass"] for row in matrix)
        assert resource[0]["continuous_cos_squared_certificate_interval"] is None
        assert resource[1]["continuous_cos_squared_certificate_interval"] is not None
        return s({
            "parameters": {"a": D(1), "m": m, "M": mass_source, "gaps": gaps,
                "epsilon_switch_log": eps, "epsilon_material_log": eps_mat,
                "coupling": coupling, "times": times, "budgets": budgets},
            "spectral_identity_diagnostics": spectral,
            "curvature_diagnostics": curvature,
            "thermal_thresholds": thresholds,
            "fixed_window_certificates": entries,
            "resource_certificates": resource,
            "fixed_plan_admission": matrix,
            "uniform_actual_error_contract_required": min(x["conditional_absolute_probability_error_allowance"] for x in entries),
        })


def window_quadrature():
    # Independent direct integrals of the two real windows, not Fourier tails.
    nodes, weights = np.polynomial.legendre.leggauss(128)
    out = []
    for t in (20.0, 30.0):
        x = nodes*t/2
        measure = weights*t/2
        for power in (1, 2):
            chi = np.cos(np.pi*x/t)**power
            derivative = -power*np.pi/t*np.cos(np.pi*x/t)**(power-1)*np.sin(np.pi*x/t)
            norm = float(measure @ (chi*chi))
            grad = float(measure @ (derivative*derivative))
            norm_exact = t/2 if power == 1 else 3*t/8
            grad_exact = np.pi**2/(2*t)
            r_exact = (1 if power == 1 else 4/3)*np.pi**2/t**2
            assert abs(norm-norm_exact) < 1e-11
            assert abs(grad-grad_exact) < 1e-12
            assert abs(grad/norm-r_exact) < 1e-12
            out.append({"T": t, "cosine_power": power, "norm_integral": norm,
                        "gradient_integral": grad, "r": grad/norm})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-results", action="store_true")
    args = ap.parse_args()
    a, b = calculate(70), calculate(95)
    assert a == b, "70/95-digit results disagree at 25 displayed significant digits"
    payload = {
        "status": "PASS", "check_groups": 6,
        "check_group_names": ["spectral_identities_and_constants", "H1_window_integrals",
            "threshold_inversion_and_domain", "closed_COM_energy", "fixed_plan_admission",
            "conditional_material_log_error_transport"],
        "scope": "Conditional finite-source/certificate interface; no actual Unruh detector or measured probability.",
        "numerical_status": "Decimal 70/95 precision cross-check plus float64 quadrature; not directed interval arithmetic.",
        "resource_domain": "Constant collinear proper acceleration with optimally centered rapidity interval; endpoint additive complete four-momentum decomposes into the on-shell detector of mass at least m and a future-causal remainder of invariant mass at least M, with opposite spatial momenta in the complete COM frame. If additional binding is allowed, the same decomposition must be proved separately; an interaction energy lower bound alone is insufficient.",
        "failure_interpretation": "Fixed-plan failure is failure of these sufficient thermal bounds plus this necessary resource bound, not proof of absent physical thermality.",
        "window_optimality_scope": "Wirtinger minimizes this second-moment certificate in H1_0; it does not prove a shortest actual thermalization time.",
        "readout_contract_status": "Required absolute error budget only; no actual source/clock/field/readout implementation is supplied.",
        "source_sha256": {"archive_1219_/1222/independent_check.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "calculations": a,
        "window_integral_diagnostics": window_quadrature(),
    }
    target = HERE / "independent_results.json"
    if args.write_results:
        with target.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write("\n")
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert saved == payload, "Independent result mismatch; refusing to rewrite"
    print(json.dumps({"status": "PASS", "check_groups": 6,
        "fixed_pairs_passed": sum(x["combined_fixed_plan_pass"] for x in a["fixed_plan_admission"]),
        "independent_results_sha256": hashlib.sha256(target.read_bytes()).hexdigest()}, ensure_ascii=False))


if __name__ == "__main__":
    main()
