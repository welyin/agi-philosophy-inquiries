"""Count-zero algebra checks for the acoustic-transport mature adoption.

No Boltzmann solver, cosmological integration or data likelihood is run.
Default is read-only; --write exclusively creates the first result.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
REPORT = HERE.parent / "acoustic_transport_adoption.md"
RESULT = HERE / "results.json"
ROOT = HERE.parents[2]
HISTORY = (
    ROOT / "archive_1046_/_admission/acoustic_transport_after1055/selection.md",
    ROOT / "archive_990_1008/research_note_1005.md",
    ROOT / "archive_1046_/_shared/neutrino_decoupling_adoption.md",
)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def number(x):
    return {"exact": str(x), "decimal": float(x)}

def check(allow_pending_result=False):
    # One exact diagnostic point; numbers are not cosmological estimates.
    rho_gamma, rho_b = F(3, 4), F(3, 5)
    inertia_gamma = F(4, 3) * rho_gamma
    R = rho_b / inertia_gamma
    Gamma, k, Hc = F(1000), F(1), F(1, 10)
    vg, vb = F(3, 100), F(-1, 50)
    gdot, bdot = Gamma * (vb-vg), Gamma/R * (vg-vb)
    assert inertia_gamma*gdot + rho_b*bdot == 0

    # Entire 2x2 local collision operator, not just one velocity pair.
    L = [[-Gamma, Gamma], [Gamma/R, -Gamma/R]]
    weights = [inertia_gamma, rho_b]
    left_null = [sum(weights[i]*L[i][j] for i in range(2)) for j in range(2)]
    right_null = [sum(row) for row in L]
    trace = L[0][0]+L[1][1]
    determinant = L[0][0]*L[1][1]-L[0][1]*L[1][0]
    assert left_null == right_null == [0, 0]
    assert determinant == 0 and trace == -Gamma*(1+1/R)
    # False independent baryon rate creates a retained-order momentum source.
    false_momentum = inertia_gamma*gdot + rho_b*(2*bdot)
    assert false_momentum != 0

    # Leading tight-coupling polarization algebra in MB F_l/G_l normalization.
    g0, g2 = F(5, 4), F(1, 4)  # G_0/F_2 and G_2/F_2
    total = 1 + g0 + g2
    assert g0 == total/2 and g2 == total/10
    f2 = F(9, 10)-(g0+g2)/10
    visc = 4/(5*f2)
    assert f2 == F(3, 4) and visc == F(16, 15)

    cs2 = 1/(3*(1+R))
    numerator = R*R+visc*(1+R)
    diffusion = numerator/(6*Gamma*(1+R)**2)
    assert F(0) < cs2 < F(1, 3) and diffusion > 0
    slip_term = R*R/(6*Gamma*(1+R)**2)
    shear_term = visc/(6*Gamma*(1+R))
    assert diffusion == slip_term+shear_term
    # Algebraic opacity rescaling leaves the leading sound speed unchanged.
    assert numerator/(6*(2*Gamma)*(1+R)**2) == diffusion/2
    sound_unchanged = cs2 == 1/(3*(1+R))
    assert sound_unchanged

    # Local fully-ionized scaling Gamma'=-2 Hc Gamma. No recombination history.
    Gamma_prime = -2*Hc*Gamma
    eps = {
        "k_over_Gamma": k/Gamma,
        "Hc_over_Gamma": Hc/Gamma,
        "abs_Gamma_prime_over_Gamma_squared": abs(Gamma_prime)/Gamma**2,
    }
    assert max(eps.values()) == F(1, 1000)
    assert max(eps.values()) < F(1, 100)
    # This is a scale diagnostic, not a theorem certifying an absolute error.

    text = REPORT.read_text(encoding="utf-8")
    local = [x for x in re.findall(r"\]\(([^)]+)\)", text)
             if not x.startswith(("http://", "https://", "#"))]
    targets = [(REPORT.parent / x.split("#")[0]).resolve() for x in local]
    assert all(p.exists() or (allow_pending_result and p == RESULT.resolve()) for p in targets)

    return {
        "date": "2026-10-09",
        "counting": {"new_science_groups": 0, "new_cognitive_axioms": 0},
        "scope": "Exact mature-relation calibration at one dimensionless diagnostic point; not a cosmological evolution, full error certificate or empirical likelihood.",
        "all_checks_passed": True,
        "momentum": {
            "R_b": number(R),
            "photon_collision_acceleration": number(gdot),
            "baryon_collision_acceleration": number(bdot),
            "weighted_collision_sum": number(F(0)),
            "left_null": [str(x) for x in left_null],
            "right_null": [str(x) for x in right_null],
            "fast_slip_eigenvalue": number(trace),
            "false_double_baryon_rate_momentum_residual": number(false_momentum),
        },
        "polarization": {"G0_over_F2": number(g0), "G2_over_F2": number(g2),
                         "effective_f2": number(f2), "shear_coefficient": number(visc)},
        "acoustic": {"sound_speed_squared": number(cs2),
                     "diffusion_integrand": number(diffusion),
                     "slip_contribution": number(slip_term),
                     "shear_contribution": number(shear_term),
                     "opacity_double_diffusion_half": True,
                     "leading_sound_speed_unchanged": sound_unchanged},
        "tight_coupling_local_diagnostics": {a: number(b) for a, b in eps.items()},
        "history_sha256": {str(p.relative_to(ROOT)): sha(p) for p in HISTORY},
        "report_sha256": sha(REPORT),
        "local_links_checked": len(local),
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = check(allow_pending_result=args.write)
    if args.write:
        with RESULT.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
    else:
        assert result == json.loads(RESULT.read_text(encoding="utf-8"))
    print(json.dumps(result, ensure_ascii=False, indent=2))

if __name__ == "__main__":
    main()
