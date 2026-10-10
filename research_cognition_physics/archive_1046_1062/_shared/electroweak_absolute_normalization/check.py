"""Small read-only matching/units check; no modern SM input fit or CMS residual.

Awramik hep-ph/0311148v3, Eqs. (6), (7), (9). All printed coefficients
are in GeV. --write exclusively creates the first results.json.
"""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import argparse
import json


C = tuple(map(D, (
    "80.3779", "0.05263", "0.010239", "0.000954", "-0.000054",
    "1.077", "0.5252", "0.0700", "0.004102", "0.000111",
    "0.0774", "115.0",
)))
REF = {
    "MH": D("100"), "mt": D("174.3"), "MZ": D("91.1875"),
    "Delta_alpha": D("0.05907"), "alpha_s": D("0.119"),
}


def prediction(p):
    H = (p["MH"]/100).ln()
    h = (p["MH"]/100)**2
    t = (p["mt"]/D("174.3"))**2-1
    z = p["MZ"]/D("91.1875")-1
    a = p["Delta_alpha"]/D("0.05907")-1
    s = p["alpha_s"]/D("0.119")-1
    (m0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = C
    return (m0-c1*H-c2*H**2+c3*H**4+c4*(h-1)-c5*a
            +c6*t-c7*t*t-c8*H*t+c9*h*t-c10*s+c11*z)


def derivatives(p):
    H = (p["MH"]/100).ln()
    h = (p["MH"]/100)**2
    t = (p["mt"]/D("174.3"))**2-1
    (_, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11) = C
    return {
        "MH": (-c1-2*c2*H+4*c3*H**3-c8*t)/p["MH"]
              +2*p["MH"]/10000*(c4+c9*t),
        "mt": 2*p["mt"]/D("174.3")**2*(c6-2*c7*t-c8*H+c9*h),
        "MZ": c11/D("91.1875"),
        "Delta_alpha": -c5/D("0.05907"),
        "alpha_s": -c10/D("0.119"),
    }


def show(x):
    return format(x, ".28g")


def compute():
    with localcontext() as ctx:
        ctx.prec = 65
        assert prediction(REF) == C[0]
        illustrative = dict(REF, MH=D("125"))
        m = prediction(illustrative)
        der = derivatives(illustrative)
        hstep = D("1e-7")
        finite_difference_error = {}
        for key, val in der.items():
            plus, minus = dict(illustrative), dict(illustrative)
            plus[key] += hstep
            minus[key] -= hstep
            numeric = (prediction(plus)-prediction(minus))/(2*hstep)
            err = abs(numeric-val)
            assert err < D("1e-15")
            finite_difference_error[key] = show(err)
        assert der["MH"] < 0 and der["mt"] > 0 and der["MZ"] > 0
        assert der["Delta_alpha"] < 0 and der["alpha_s"] < 0

        # Exact rational real/imaginary parts of the running-width pole.
        mz, gz = F("91.1875"), F("2.4952")
        gamma = gz/mz
        real, minus_imag = mz*mz/(1+gamma*gamma), mz*gz/(1+gamma*gamma)
        assert real*(1+gamma*gamma) == mz*mz
        assert minus_imag == real*gamma
        pole_mass = (D(real.numerator)/D(real.denominator)).sqrt()
        pole_width = D("2.4952") / (1+(D("2.4952")/D("91.1875"))**2).sqrt()
        assert abs(pole_mass*pole_width
                   - D(minus_imag.numerator)/D(minus_imag.denominator)) < D("1e-55")

        # Natural-energy dimensions: G_mu^2 m_mu^5 is a rate (E^1);
        # alpha/G_mu and the mass-matching LHS both have dimension E^2.
        dim_G, dim_mass = -2, 1
        assert 2*dim_G+5*dim_mass == 1
        assert -dim_G == 2*dim_mass

        # Algebraic fixed-Delta-r sensitivity only. This does not update Eq.(9).
        mw, mz_d = C[0], REF["MZ"]
        ratio = (mw/mz_d)**2
        sensitivity = mw*(1-ratio)/(2*(2*ratio-1))
        g0, g_mulan = D("1.166379e-5"), D("1.1663788e-5")
        dg_log = (g_mulan/g0).ln()
        dM_MeV = sensitivity*dg_log*1000
        assert sensitivity > 0 and abs(dM_MeV) < D("0.004")
        # The small MuLan/reference mismatch is diagnostic, not a refit.

        return {
            "scope": "published_parametrization_and_algebra_only_no_new_MW_test",
            "theory_version": "hep-ph/0311148v3, Eq.(6),(7),(9)",
            "coefficients_GeV": [str(v) for v in C],
            "reference_inputs": {k: str(v) for k, v in REF.items()},
            "reference_MW_GeV": show(prediction(REF)),
            "illustrative_MH125_other_inputs_reference_MW_GeV": show(m),
            "derivatives_at_illustrative_point_GeV_per_input_unit":
                {k: show(v) for k, v in der.items()},
            "derivative_finite_difference_max_errors":
                finite_difference_error,
            "mass_scheme": {
                "running_Z_M_GeV": str(mz),
                "running_Z_width_GeV": str(gz),
                "pole_Z_M_GeV": show(pole_mass),
                "pole_Z_width_GeV": show(pole_width),
                "complex_pole_identity_exact": True,
            },
            "fixed_Delta_r_diagnostic_only": {
                "dMW_dlogG_GeV": show(sensitivity),
                "log_MuLan_G_over_formula_G": show(dg_log),
                "linear_shift_MeV_not_applied_to_prediction": show(dM_MeV),
            },
            "MuLan_transcription": {
                "tau_ps": "2196980.3", "tau_uncertainty_ps": "2.2",
                "G_GeVminus2": "1.1663788e-5", "G_uncertainty_GeVminus2": "7e-12",
            },
            "CMS_v2_transcription_only": {
                "MW_MeV": "80360.2", "uncertainty_MeV": "9.9",
                "mass_scheme": "running_width",
                "no_residual_or_new_significance_computed": True,
            },
            "natural_units_consistent": True,
            "new_scientific_groups": 0,
            "all_assertions_passed": True,
        }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = compute()
    path = Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")
        print("First-save algebra and parametrization checks passed; no new fit.")
    else:
        assert json.loads(path.read_text(encoding="utf-8")) == result
        print("Read-only checks passed; CMS likelihood and SM inputs not refitted.")
