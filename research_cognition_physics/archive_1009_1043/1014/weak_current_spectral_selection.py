"""Round 1014: finite matrix calibration of a conditional weak-current test.

This is not an experimental event-rate fit. All strengths are coefficients of
the specified canonical left-handed current after a complete final-state sum.
The continuum-in-r theorem is proved in the report, not by these samples.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
OUT = HERE / "weak_current_spectral_selection_results.json"
OLD_RESULT = BASE / "archive_1009_/1012/lepton_mass_operator_classification_results.json"
N = 9
J = np.zeros((N, N), dtype=complex)
J[3:6, :3] = np.eye(3)
PN = np.diag([1., 1., 1., 0., 0., 0., 0., 0., 0.])


def norm(a):
    return float(np.linalg.norm(a, 2))


def serialize_complex(a):
    return [[[float(z.real), float(z.imag)] for z in row] for row in a]


def decode(a):
    a = np.asarray(a, dtype=float)
    return a[..., 0] + 1j*a[..., 1] if a.ndim == 3 else a.astype(complex)


def takagi(m):
    """Full-rank complex symmetric Takagi factorization via a real eigenproblem.

    M u = s conjugate(u); its positive real eigenvectors give an orthonormal
    complex basis even when positive singular values repeat. No SVD columns
    are silently assumed to satisfy the Takagi phase/congruence condition.
    """
    n = m.shape[0]
    assert norm(m-m.T) < 1e-11
    real_problem = np.block([[m.real, -m.imag], [-m.imag, -m.real]])
    vals, vecs = np.linalg.eigh(real_problem)
    positive = np.flatnonzero(vals > 1e-11)
    assert len(positive) == n, "Calibration samples are explicitly full rank."
    masses = vals[positive]
    u = vecs[:n, positive] + 1j*vecs[n:, positive]
    # Only the allowed Takagi sign is fixed, not an arbitrary complex phase.
    for col in range(n):
        pivot = u[np.argmax(np.abs(u[:, col])), col]
        if pivot.real < -1e-13 or (abs(pivot.real) <= 1e-13 and pivot.imag < 0):
            u[:, col] *= -1
    unitary_residual = norm(u.conj().T@u-np.eye(n))
    takagi_residual = norm(u.T@m@u-np.diag(masses))
    assert unitary_residual < 2e-10
    assert takagi_residual < 2e-10*max(1., norm(m))
    return masses, u, unitary_residual, takagi_residual


def clusters(masses):
    """A numerical grouping for finite calibration, not an experimental test."""
    tol = 2e-8*max(1., float(np.max(masses)))
    out = []
    for idx, mass in enumerate(masses):
        if out and abs(mass-masses[out[-1][0]]) <= tol:
            out[-1].append(idx)
        else:
            out.append([idx])
    return out


def analyze(m, q, current=J):
    masses, u, ur, tr = takagi(m)
    sets = clusters(masses)
    simple = [group[0] for group in sets if len(group) == 1 and masses[group[0]] > 1e-10]
    s = u[:, simple]
    columns = current@s
    strengths = np.sum(np.abs(columns)**2, axis=0)
    gram = columns.conj().T@columns
    transported = u.conj().T@current@u
    by_final_mass_states = np.sum(np.abs(transported[:, simple])**2, axis=0)
    total = float(np.sum(np.abs(transported)**2))
    charge_leakage = norm(q@s) if simple else 0.
    ward = norm(q.T@m+m@q)
    assert ward < 1e-10*max(1., norm(m)*norm(q))
    assert charge_leakage < 2e-9
    assert np.max(np.abs(strengths-by_final_mass_states), initial=0.) < 2e-10
    assert abs(total-3.) < 2e-10
    return dict(masses=masses, u=u, simple=simple, groups=sets,
                strengths=strengths, gram=gram, total=total,
                F=float(np.sum(strengths)), R=float(np.sum(strengths))/total,
                charge_leakage=charge_leakage, ward_residual=ward,
                unitary_residual=ur, takagi_residual=tr,
                final_basis_transport_residual=float(np.max(np.abs(strengths-by_final_mass_states), initial=0.)))


def public_analysis(row):
    return dict(masses=row["masses"].tolist(),
                masses_and_multiplicities=[[float(row["masses"][g[0]]), len(g)] for g in row["groups"]],
                positive_simple_indices=row["simple"],
                positive_simple_strengths=row["strengths"].tolist(),
                simple_current_gram=serialize_complex(row["gram"]),
                simple_current_gram_eigenvalues=np.linalg.eigvalsh(row["gram"]).tolist(),
                simple_current_gram_rank=int(np.linalg.matrix_rank(row["gram"], tol=1e-9)),
                F=row["F"], R=row["R"], F_all=row["total"],
                zero_charge_residual=row["charge_leakage"],
                ward_residual=row["ward_residual"],
                takagi_residual=row["takagi_residual"],
                unitarity_residual=row["unitary_residual"],
                final_basis_transport_residual=row["final_basis_transport_residual"])


def random_unitary(rng, n):
    z = rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q, r = np.linalg.qr(z)
    return q@np.diag(np.diag(r)/np.abs(np.diag(r)))


def compare(fresh, saved):
    if isinstance(fresh, dict):
        assert fresh.keys() == saved.keys(), "Dictionary keys differ."
        for key in fresh:
            compare(fresh[key], saved[key])
    elif isinstance(fresh, list):
        assert len(fresh) == len(saved)
        for a, b in zip(fresh, saved):
            compare(a, b)
    elif isinstance(fresh, float):
        assert math.isclose(fresh, saved, rel_tol=4e-8, abs_tol=4e-10), (fresh, saved)
    else:
        assert fresh == saved, (fresh, saved)


def run():
    saved = json.loads(OLD_RESULT.read_text(encoding="utf8"))
    old = saved["isospectral_counterexample"]
    assert norm(J.conj().T@J-PN) == 0
    standard = decode(old["standard_mass_matrix"])
    nonstandard = decode(old["nonstandard_mass_matrix"])
    q0 = np.diag([float(Fraction(x)) for x in old["standard_electric_charges"]])
    q2 = np.diag([float(Fraction(x)) for x in old["nonstandard_electric_charges"]])
    a0, a2 = analyze(standard, q0), analyze(nonstandard, q2)
    assert np.max(np.abs(a0["masses"]-a2["masses"])) < 2e-13
    assert len(a0["simple"]) == len(a2["simple"]) == 3
    assert np.max(np.abs(a0["strengths"]-np.ones(3))) < 2e-12
    assert np.max(np.abs(a2["strengths"]-[.25,.5,.25])) < 2e-12
    assert norm(a0["gram"]-np.eye(3)) < 2e-12
    assert np.max(np.abs(np.linalg.eigvalsh(a2["gram"])-[0,0,1])) < 2e-12
    assert abs(a0["F"]-3) < 2e-12 and abs(a2["F"]-1) < 2e-12

    # Whole degenerate cluster, never an arbitrarily chosen vector from it.
    cluster11 = [i for i, mass in enumerate(a2["masses"]) if abs(mass-.11)<1e-10]
    single04 = [i for i, mass in enumerate(a2["masses"]) if abs(mass-.04)<1e-10]
    assert len(cluster11)==2 and len(single04)==1
    u2 = a2["u"]
    p11 = u2[:,cluster11]@u2[:,cluster11].conj().T
    p04 = u2[:,single04]@u2[:,single04].conj().T
    f11 = float(np.trace(PN@p11).real)
    f04 = float(np.trace(PN@p04).real)
    assert abs(f11-1.) < 1e-12 and abs(f04-.5) < 1e-12
    assert abs(float(np.trace(PN@(p11+p04)).real)-1.5) < 1e-12

    rng = np.random.default_rng(1014)
    covariance = []
    for name, m, q, prior in [("standard",standard,q0,a0),("nonstandard",nonstandard,q2,a2)]:
        w = random_unitary(rng, N)
        transformed = analyze(w.T@m@w, w.conj().T@q@w, w.conj().T@J@w)
        residual = max(abs(transformed["F"]-prior["F"]), abs(transformed["R"]-prior["R"]))
        assert residual < 2e-10
        covariance.append(dict(candidate=name, shared_field_basis_change_residual=residual,
                               transformed_takagi_residual=transformed["takagi_residual"]))

    # Independent exact charge-pair mass basis; its dimension must agree with
    # the operator-image dimension proved and saved in 1012. This does not
    # import/reexecute that operator implementation.
    branch_rows = []
    maxima = dict(ward=0., zero_charge=0., unitary=0., takagi=0., final_basis=0.)
    for branch in saved["parameter_branches"]:
        r = Fraction(branch["r"])
        h = Fraction(1,2)
        x = [h*r, -h*r, Fraction(0)]
        charges = x+[v-2*h for v in x]+[2*h-v for v in x]
        pairs = [(i,j) for i in range(N) for j in range(i,N) if charges[i]+charges[j]==0]
        assert len(pairs) == branch["independent_mass_matrix_dimension"]
        q = np.diag([float(v) for v in charges])
        pn_zero_dim = sum(v==0 for v in x)
        assert pn_zero_dim == (3 if r==0 else 1)
        samples=[]
        for sample in range(4):
            m = np.zeros((N,N), complex)
            for i,j in pairs:
                value = rng.normal()+1j*rng.normal()
                m[i,j] = m[j,i] = value
            row = analyze(m,q)
            assert row["F"] <= pn_zero_dim+2e-9
            if r==0:
                assert abs(row["F"]-len(row["simple"]))<2e-9
            samples.append(dict(sample=sample, positive_simple_count=len(row["simple"]),
                                F=row["F"], R=row["R"],
                                simple_current_gram_rank=int(np.linalg.matrix_rank(row["gram"],tol=1e-8))))
            for key, value in [("ward",row["ward_residual"]),("zero_charge",row["charge_leakage"]),
                               ("unitary",row["unitary_residual"]),("takagi",row["takagi_residual"]),
                               ("final_basis",row["final_basis_transport_residual"])]:
                maxima[key]=max(maxima[key],value)
        branch_rows.append(dict(r=str(r), symmetric_mass_dimension=len(pairs),
                                zero_charge_dimension=sum(v==0 for v in charges),
                                upper_zero_charge_dimension=pn_zero_dim,
                                simple_strength_bound=pn_zero_dim, samples=samples))

    common_g=[]
    for coupling in [.2, .7, 2.5]:
        for name, row in [("standard",a0),("nonstandard",a2)]:
            # The same g multiplies the complete transported current.
            transformed=coupling*(row["u"].conj().T@J@row["u"])
            numerator=float(np.sum(np.abs(transformed[:,row["simple"]])**2))
            denominator=float(np.sum(np.abs(transformed)**2))
            ratio=numerator/denominator
            assert abs(ratio-row["R"])<2e-12
            common_g.append(dict(g=coupling,candidate=name,weighted_target=numerator,
                                 weighted_all=denominator,R=ratio))

    # Incompletely collected, but absolutely normalized, final-state effects
    # provide lower bounds. They cannot evade the nonstandard bound by losing
    # channels; completeness is only needed for equality with PN and F_all=3.
    partial_final=[]
    observed_unitary=random_unitary(rng,N)
    random_projector=observed_unitary[:,:6]@observed_unitary[:,:6].conj().T
    effects=[("omit_last_final_mass_mode",None),("random_rank6_projector",random_projector),
             ("uniform_known_efficiency",.7*np.eye(N))]
    for name,row in [("standard",a0),("nonstandard",a2)]:
        target=row["u"][:,row["simple"]]
        for label,effect in effects:
            if effect is None:
                last=row["u"][:,-1]
                effect=np.eye(N)-np.outer(last,last.conj())
            eigenvalues=np.linalg.eigvalsh(effect)
            assert min(eigenvalues)>-1e-12 and max(eigenvalues)<1+1e-12
            f_observed=float(np.trace(target.conj().T@J.conj().T@effect@J@target).real)
            assert -1e-12<=f_observed<=row["F"]+1e-12
            if name=="nonstandard":
                assert f_observed<=1+1e-12
            partial_final.append(dict(candidate=name,effect=label,F_observed=f_observed,
                                      F_inclusive=row["F"],absolute_lower_bound=True,
                                      success_conditioned=False))

    # A general instrument example: A=aV is an isometric success filter.
    # The missing branch B=sqrt(1-a^2)I completes the instrument. Conditioning
    # only on success erases a. This is not claimed to be a second SM model.
    omega=np.exp(2j*np.pi/3)
    v=np.array([[1,1],[1,omega],[1,omega**2]],complex)/math.sqrt(3)
    assert norm(v.conj().T@v-np.eye(2)) < 1e-14
    rho=np.array([[.6,.15+.1j],[.15-.1j,.4]],complex)
    assert min(np.linalg.eigvalsh(rho))>0
    conditional=[]
    for amplitude in [.5,1.]:
        a=amplitude*v
        b=math.sqrt(1-amplitude**2)*np.eye(2)
        assert norm(a.conj().T@a+b.conj().T@b-np.eye(2))<1e-14
        success=a@rho@a.conj().T
        probability=float(np.trace(success).real)
        normalized=success/probability
        assert norm(normalized-v@rho@v.conj().T)<1e-14
        conditional.append(dict(amplitude=amplitude,success_probability=probability,
                                absolute_two_input_strength=float(np.trace(a.conj().T@a).real),
                                normalized_success_state=serialize_complex(normalized)))
    assert abs(conditional[0]["absolute_two_input_strength"]-.5)<1e-12
    assert abs(conditional[1]["absolute_two_input_strength"]-2)<1e-12

    finite_error=[]
    for label,estimate,error,bound,expected in [
            ("F_clear",1.1,.05,1.,True), ("F_inconclusive",1.05,.06,1.,False),
            ("F_boundary",1.05,.05,1.,False),
            ("R_clear",.4,.02,1/3,True), ("R_inconclusive",.35,.03,1/3,False)]:
        certified = estimate-error > bound+1e-14
        assert certified==expected
        finite_error.append(dict(label=label,estimate=estimate,absolute_error_bound=error,
                                 null_bound=bound,lower_endpoint=estimate-error,
                                 excludes_nonstandard_if_error_contract_is_valid=certified))

    sources=[BASE/"archive_990_1008/research_note_1007.md",
             BASE/"archive_990_1008/1007/neutrino_observation_adoption_v1.md",
             BASE/"archive_531_553/research_note_545.md",
             BASE/"archive_531_553/research_note_547.md",
             BASE/"archive_1009_/research_note_1012.md",OLD_RESULT,
             BASE/"archive_1009_/research_note_1013.md"]
    return dict(round=1014,new_calibration_groups=1,cumulative_test_groups=3792,
        new_cognitive_axioms=0,all_scientific_calibrations_passed=True,
        scope="canonical nine left-handed Weyl fields; constant symmetric mass; exact residual Ward identity; specified chiral current; complete final-state coefficient sum",
        full_realistic_experiment_claimed=False,ordinary_event_rate_identity_claimed=False,
        mass_simplicity_inferred_from_numerical_tolerance_claimed=False,
        continuum_parameter_theorem_proved_by_sampling_claimed=False,
        charged_current=serialize_complex(J),current_gram_equals_PN_exact=True,
        complete_current_strength=3,
        reused_isospectral_candidates=dict(standard=public_analysis(a0),nonstandard=public_analysis(a2),
            equal_mass_residual=float(np.max(np.abs(a0["masses"]-a2["masses"]))),
            same_current_strengths=False),
        shared_field_basis_covariance=covariance,
        full_spectrum_simplicity_deletion_counterexample=dict(r="2",cluster_mass=.11,cluster_rank=2,
            cluster_strength=f11,positive_simple_mass=.04,single_strength=f04,
            combined_rank=3,combined_strength=f11+f04,
            complete_degenerate_cluster_used=True,refutes_theorem_with_simple_condition_removed=True),
        random_complex_ward_space_calibration=dict(seed=1014,samples_per_branch=4,
            branches=branch_rows,maximum_residuals=maxima),
        common_coupling_cancellation=common_g,
        partial_final_state_lower_bounds=partial_final,
        success_conditioning_counterexample=dict(type="general filter instrument, not a physical hypercharge candidate",
            cases=conditional,success_conditioned_output_is_identical=True),
        finite_error_decisions=finite_error,
        analytic_obligations="positive simple Takagi columns have zero residual charge; orthogonality bounds their total PN weight; nonstandard N intersect kerQ has dimension one",
        retained_inputs=["field content and weak representation","canonical quadratic description",
            "exact unbroken charge and original anomaly/Yukawa branch","complete-spectrum positive simple states",
            "specified current and calibrated absolute lower bounds or complete-total relative strengths",
            "common source, kinematic and detection matching if compared with data"],
        historical_source_sha256={str(p.relative_to(BASE)).replace("\\","/"):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
    else:
        compare(result,json.loads(OUT.read_text(encoding="utf8")))
    print(json.dumps(dict(round=1014,passed=result["all_scientific_calibrations_passed"],
        old_equal_mass_residual=result["reused_isospectral_candidates"]["equal_mass_residual"],
        standard_F=result["reused_isospectral_candidates"]["standard"]["F"],
        nonstandard_F=result["reused_isospectral_candidates"]["nonstandard"]["F"],
        nonsimple_counterexample_F=result["full_spectrum_simplicity_deletion_counterexample"]["combined_strength"],
        random_samples=44,maximum_residuals=result["random_complex_ward_space_calibration"]["maximum_residuals"]),indent=2))
