"""Round 1262: finite one-particle source insertion; default replay is read-only.

Two author groups. Finite samples diagnose algebra only; proof.md gives
the continuous uniform estimates. No ODE, QFT evolution, or postselection.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
SIGN = [-1, 1, 1, 1]
SOURCE_HASHES = {
    "../1253/proof.md": "422c736a1c3d65671fcef1795cc889e44514516344510048be7c8240acc29161",
    "../1261/parent.md": "1fd650bd75b2c345d599c602d6a13cf35f3263bfe5b14a6c6c75b6d2cc73776d",
    "../1261/proof.md": "cc8f311b405f637cdb1a29a6408b86ccf2385d61ce95736a1af0a833c34dd4ef",
    "../../archive_935_955/research_note_942.md": "0f51b118cb170d547f9ec5edbdef1c78bca5c35ba8379b56e60904b98356374e",
}

def rational_shell(m, v):
    vv = sum(x*x for x in v)
    assert vv < 1
    return [m*(1+vv)/(1-vv)] + [2*m*x/(1-vv) for x in v]

def numerator(m, p, q):
    dot = sum(F(SIGN[i])*p[i]*q[i] for i in range(4))
    return [
        [p[i]*q[j]+q[i]*p[j] -
         (F(SIGN[i]) if i == j else 0)*(dot+m*m) for j in range(4)]
        for i in range(4)
    ]

def diagnostic_pair(m, p, q):
    Ep = np.sqrt(m*m + np.dot(p,p))
    Eq = np.sqrt(m*m + np.dot(q,q))
    P = (p+q)/2
    E0 = np.sqrt(m*m + np.dot(P,P))
    d2 = float(np.dot(p-q,p-q))
    W = np.zeros((4,4))
    W[0,0] = E0
    W[0,1:] = P
    W[1:,0] = P
    if E0:
        W[1:,1:] = np.outer(P,P)/E0
    if Ep*Eq:
        pv = np.r_[Ep,p]
        qv = np.r_[Eq,q]
        eta = np.diag(SIGN)
        K = (np.outer(pv,qv)+np.outer(qv,pv)
             -eta*(pv@eta@qv+m*m))/(2*np.sqrt(Ep*Eq))
    else:
        # Continuous massless kernel extension at zero energy.
        assert m == 0
        K = np.zeros((4,4))
    error = float(np.max(np.abs(K-W)))
    if d2 == 0:
        assert error < 1e-10
        return 0., 0.
    fine_ratio = error*min(Ep,Eq)/d2
    coarse_ratio = error/np.sqrt(d2)
    t = d2/4
    a = E0*E0
    scale = max(a,t,1.)
    assert a-t-1e-8*scale <= Ep*Eq <= a+t+1e-8*scale
    if min(Ep,Eq) > 0:
        assert fine_ratio <= .75+1e-10
    assert coarse_ratio <= 1.75+1e-10
    return fine_ratio, coarse_ratio

def calculate():
    sources = {}
    for rel, expected in SOURCE_HASHES.items():
        got = hashlib.sha256((HERE/rel).read_bytes()).hexdigest()
        assert got == expected, ("Frozen source changed", rel, got)
        sources[rel] = got

    # Group 1: exact mass-shell identities plus bounded diagnostic samples.
    count = 0
    for m in [F(26),F(31),F(12),F(10**20)]:
        p = rational_shell(m,[F(1,5),F(2,7),F(-1,6)])
        q = rational_shell(m,[F(-1,4),F(1,5),F(2,9)])
        N = numerator(m,p,q)
        for j in range(4):
            assert sum(SIGN[i]*(p[i]-q[i])*N[i][j]
                       for i in range(4)) == 0
        D = numerator(m,p,p)
        for i in range(4):
            for j in range(4):
                assert D[i][j]/(2*p[0]) == p[i]*p[j]/p[0]
        count += 1
    for p,q in [
        ([F(1),F(1),F(0),F(0)],[F(5),F(3),F(4),F(0)]),
        ([F(1),F(1),F(0),F(0)],[F(1),F(-1),F(0),F(0)])
    ]:
        N = numerator(F(0),p,q)
        for j in range(4):
            assert sum(SIGN[i]*(p[i]-q[i])*N[i][j]
                       for i in range(4)) == 0
        count += 1
    # m_* = 17+2sqrt(19): original six masses lie strictly in (26,32).
    assert F(4)**2 < 19 < F(9,2)**2
    rng = np.random.default_rng(126101)
    pairs = []
    for m in [0.,1.,26.,31.,12.,1e5]:
        for _ in range(120):
            p = rng.normal(size=3)*float(rng.uniform(.1,50))
            q = rng.normal(size=3)*float(rng.uniform(.1,50))
            pairs.append(diagnostic_pair(m,p,q))
    pairs.append(diagnostic_pair(0.,np.array([1.,0,0]),np.array([-1.,0,0])))
    zero_pairs = [
        diagnostic_pair(0.,np.zeros(3),np.array([3.,4.,0.])),
        diagnostic_pair(0.,np.zeros(3),np.zeros(3))
    ]

    # Group 1 continued: constants of the proven Schur estimates.
    # Fine kernel = 3|d|^2/(4 e_min), Gaussian second moment = 3/L^2.
    assert F(3,4)*3 == F(9,4)
    # Coarse case split: sqrt(2)<3/2, (2+sqrt(2))/2<7/4.
    assert F(2) < F(3,2)**2
    assert (F(2)+F(3,2))/2 == F(7,4)
    # Gaussian first moment = 2sqrt(2/pi)/L < 2/L (pi>2).
    assert F(7,4)*2 == F(7,2)

    # Group 2: finite normal preparation and finite source-only budgets.
    b = F(10**20)
    sigma = F(10**13)
    k = F(10**8)
    mu = b
    assert sigma**20 == b**13
    assert k**5 == b**2
    gap = k-4/sigma
    assert gap > k/2
    fine_two = F(9,4)/b**2*(F(1,26)+2/k)
    fine_three = fine_two + F(9,4)/(mu*b**2)
    assert fine_two < F(9,10**42)
    assert fine_three < F(9,10**42)
    no_band_two = F(9,4)/(26*b**2)+F(7,2)/b
    no_band_three = no_band_two + F(9,4)/(mu*b**2)
    assert no_band_two < F(36,10**21)
    assert no_band_three < F(36,10**21)

    return {
        "round": 1262,
        "status": "PASS",
        "author_scientific_groups": 2,
        "group_1_kernel": {
            "exact_Ward_and_shell_cases": count,
            "positive_energy_numerical_diagnostics": len(pairs),
            "zero_energy_extension_diagnostics": len(zero_pairs),
            "max_fine_ratio_bound_three_quarters": max(x[0] for x in pairs),
            "max_coarse_ratio_bound_seven_quarters":
                max(x[1] for x in pairs+zero_pairs),
            "numeric_samples_prove_uniform_bound": False,
            "original_material_dimension": 864,
            "original_mass_sectors": 6,
            "original_mass_bounds": "26<M<32",
        },
        "group_1_smear": {
            "fine_kernel_coefficient": "3/(4 e_min)",
            "Gaussian_second_moment": "3/L^2",
            "fine_source_bound": "9/(4 e_min L^2)",
            "coarse_global_kernel_coefficient": "7/4",
            "Gaussian_first_moment": "2 sqrt(2/pi)/L < 2/L",
            "global_no_band_source_bound": "7/(2L)",
            "Gaussian_compact_support": False,
            "individual_source_bounded": False,
            "source_difference_bounded": True,
        },
        "group_2_finite_task": {
            "b_L": str(b), "sigma": str(sigma), "k": str(k),
            "initial_signal_gap": str(gap),
            "fine_two_body_bound_exact": str(fine_two),
            "fine_two_body_bound_float": float(fine_two),
            "optional_finite_B_bound_exact": str(fine_three),
            "no_band_two_body_bound_exact": str(no_band_two),
            "no_band_two_body_bound_float": float(no_band_two),
            "optional_finite_B_no_band_exact": str(no_band_three),
            "error_units": "energy of dimensionless spatially weighted source",
            "finite_scattering_time_budget_inherited": False,
        },
        "scope": {
            "minimal_KG_xi": 0,
            "one_particle_relative_vacuum_source_only": True,
            "all_mass_degeneracies_retained": True,
            "cross_mass_coherence_retained": True,
            "uniform_unknown_internal_reference": True,
            "positive_band_required_for_L_minus_2_massless_bound": True,
            "positive_band_required_for_any_finite_source_bound": False,
            "AW_source_matching_proved": False,
            "actual_interacting_time_band_proved": False,
            "finite_scattering_unitary_proved": False,
            "actual_Mr_source_instrument_proved": False,
            "quantum_field_energy_source_certified": False,
        },
        "source_sha256": sources,
    }

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results","--write",action="store_true",
                        help="Create results.json exclusively; never overwrite.")
    args = parser.parse_args()
    result = calculate()
    target = HERE/"results.json"
    if args.write_results:
        with target.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,sort_keys=True)
            stream.write("\n")
        print("WRITE_ONCE_PASS_1262: 2 author groups; exclusive results.json created")
    else:
        saved = json.loads(target.read_text(encoding="utf-8-sig"))
        assert saved == result, "Saved result mismatch; no files changed."
        print("READ_ONLY_MATCH_1262: 2 author groups; no files changed")
    print("fine two-body source bound:", result["group_2_finite_task"]["fine_two_body_bound_exact"])
    print("global no-band bound:", result["group_2_finite_task"]["no_band_two_body_bound_exact"])

if __name__ == "__main__":
    main()
