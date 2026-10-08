"""Reproduce the limited 1021 assets; no active navigation is frozen."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import math
import re
from urllib.parse import unquote

import quotient_sector_reference as science

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
ROOT = BASE.parent
NOTE = HERE.parent/"research_note_1021.md"
OUT = HERE/"research_round_1021_checks.json"
OWN = ["quotient_sector_reference.py", "quotient_sector_reference_results.json",
       "review.md", "selection_audit.md", "input_dependency_update_v0_10.md",
       "NEXT.md", "verify_round1021.py"]
HISTORICAL = {
    "archive_1009_/research_note_1020.md",
    "archive_1009_/1020/passivity_reference_selection_results.json",
    "archive_531_553/research_note_531.md",
    "archive_585_628/research_note_603.md",
    "archive_585_628/research_note_617.md",
    "archive_629_652/research_note_637.md",
    "archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md",
    "archive_1009_/1009/input_dependency_ledger_v0_1.md",
}
LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)|^\[[^\]\n]+\]:\s*(\S+)\s*$", re.M)
EXCLUDED = re.compile(r"\\\[.*?\\\]|\$\$.*?\$\$|^```.*?^```\s*$", re.S | re.M)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(a, b):
    assert math.isclose(a, b, rel_tol=3e-10, abs_tol=3e-12), (a, b)


def local_links(text):
    masked = EXCLUDED.sub(lambda match: " "*len(match[0]), text)
    for match in LINK.finditer(masked):
        value = match[1] if match[1] is not None else match[2]
        stripped = value.strip()
        if stripped.startswith("<") and ">" in stripped:
            target = stripped[1:stripped.index(">")]
        else:
            target = re.split(r"\s+[\"']", stripped, 1)[0]
        local = unquote(target.split("#", 1)[0])
        if not local or re.match(r"^[a-zA-Z]+:", local) or local.startswith("//"):
            continue
        yield target, local


def verify(prospective=False):
    fresh = science.run()
    science.compare(fresh, json.loads((HERE/"quotient_sector_reference_results.json").read_text("utf8")))
    assert fresh["round"] == 1021 and fresh["all_scientific_calibrations_passed"]
    assert fresh["new_calibration_groups"] == 1 and fresh["cumulative_test_groups"] == 3799
    assert fresh["new_cognitive_axioms"] == 0
    for key in ("quotient_group_uniquely_selected_by_matter_kernel",
                "quotient_reference_equal_to_full_cover_Gibbs", "all_SM_full_H_spectrum_computed",
                "continuum_or_nonperturbative_field_theory_proved",
                "general_global_topology_classification_proved",
                "actual_thermalization_or_measurement_device_built"):
        assert fresh[key] is False, key
    assert fresh["same_sector_representation_equivalence_preserved"]

    gauss = fresh["Gauss_and_cycle_count"]
    assert gauss["cover_enumerated_electric_triples"] == 1728
    assert gauss["cover_Gauss_basis"] == [[q,q,q] for q in range(12)]
    assert gauss["quotient_Gauss_basis"] == [[0,0,0],[1,1,1]]
    assert gauss["kernel_ring_cycle_flows"] == [[r,r,r] for r in range(6)]
    assert gauss["kernel_tree_cycle_flows"] == [[0,0]]
    assert (gauss["cover_dimension"], gauss["quotient_dimension"]) == (12,2)
    assert (gauss["ring_cycle_sector_count"], gauss["tree_cycle_sector_count"]) == (6,1)
    assert not gauss["full_configuration_matrix_constructed"]

    processes = fresh["process_intertwining"]
    assert [row["lambda_value"] for row in processes] == [0.,.2]
    for row in processes:
        assert len(row["histories"]) == 4
        assert max(row["errors"].values()) < 3e-13
        assert row["Hamiltonian_source_commutator_norm"] > 1
        assert row["unknown_inputs_and_untouched_reference_retained"]
        assert row["analytic_all_time_all_finite_descending_histories"]
        assert not row["claim_for_arbitrary_cover_sector_input"]
        close(row["history_total_probability"], 1)
        close(sum(h["branch_probability_on_maximally_entangled_input"] for h in row["histories"]), 1)
        for history in row["histories"]:
            assert history["reference_Choi_error"] < 3e-13
            assert history["full_Kraus_intertwining_error"] < 3e-13

    tw = fresh["coherence_twirl"]
    assert tw["phase_average_mask_error"] < 3e-13
    assert not tw["twirl_is_sector_postselection"]
    assert tw["coherences_within_each_sector_retained"]
    assert tw["normalized_selective_branch_can_reweight_sectors"]
    assert not tw["sector_transfer_in_any_branch"]
    close(tw["initial_purity"], 1)
    assert 0 < tw["twirled_purity"] < 1
    close(sum(tw["initial_sector_weights"]), 1)
    close(sum(tw["selected_plus_weights"]), 1)
    for initial, averaged, nonselective in zip(tw["initial_sector_weights"],
            tw["final_sector_weights"], tw["nonselective_instrument_weights"]):
        close(initial, averaged); close(initial, nonselective)
    assert max(abs(a-b) for a,b in zip(tw["initial_sector_weights"],tw["selected_plus_weights"])) > 1e-3
    close(tw["distance_to_sector0_conditioned_state"], 1-tw["initial_sector_weights"][0])

    cert = fresh["finite_probability_certificate"]
    K, v = F(27,10), F(2,5)
    radius_squared = 4*K*K+v*v
    assert F(cert["R0_squared"]) == radius_squared == F(733,25)
    assert F(27,5)**2 < radius_squared < F(11,2)**2
    lower = sum((F(108,25)**j/F(math.factorial(j)) for j in range(13)), F(0))
    x=F(11,5)
    upper = sum((x**j/F(math.factorial(j)) for j in range(9)), F(0))+x**9/F(math.factorial(9))/(1-x/10)
    assert F(cert["exp_108_over25_lower"]) == lower > F(197,3)
    assert F(cert["exp_11_over5_geometric_tail_upper"]) == upper < F(19,2)
    assert 1+x > 2 and x/10 < 1
    signal = F(27,55)*F(97,100)/30
    assert F(cert["probability_gap_strict_lower_exact"]) == signal == F(873,55000)
    close(cert["probability_gap_strict_lower"],float(signal))
    assert cert["all_bounds_exact_rational"] and cert["finite_nonzero_signal_certified"]
    assert not cert["physical_thermal_preparation_and_instrument_implementation_proved"]

    cases=fresh["thermal_cases"]
    assert [(r["lambda_value"],r["beta"]) for r in cases] == [(lam,beta) for lam in (0.,.2) for beta in (.4,1.2,3.)]
    maximum_identity_error=0.; maximum_derivative_crosscheck_error=0.
    for row in cases:
        beta=row["beta"]; K=2.7*math.exp(-2*row["lambda_value"]); v=.4
        close(row["total_K"],K)
        weights=[]; z_values=[]; e_values=[]; s_values=[]; second_values=[]
        for r, block in enumerate(row["sectors"]):
            x=math.cos(math.pi*r/6)**2
            R=math.sqrt(4*K*K*x+v*v); t=math.tanh(beta*R)
            Z=2*math.exp(-beta*(2*K+v))*math.cosh(beta*R)
            effect=.5-K*x*t/R
            source=-4*K+8*K*K*x*t/R
            # Independently differentiate the scalar closed free energy twice.
            second=8*K-32*K*K*x*t/R+64*K**4*x*x*t/R**3-64*beta*K**4*x*x/(R*R*math.cosh(beta*R)**2)
            close(block["Z"],Z); close(block["closed_form_Z"],Z)
            close(block["energy_effect"],effect); close(block["source"],source)
            close(block["F2"],second); close(block["spectral_radius"],R)
            weights.append(block["weight"]); z_values.append(Z); e_values.append(effect)
            s_values.append(source); second_values.append(second)
        Ztotal=sum(z_values)
        for w,z in zip(weights,z_values):
            close(w,z/Ztotal); assert w>0
        close(row["full_Z"],Ztotal); close(row["quotient_Z"],z_values[0])
        close(sum(weights),1); w0=weights[0]; assert 0<w0<1
        close(row["sector0_weight"],w0)
        close(row["trace_distance"],1-w0)
        close(row["relative_entropy_supported"],-math.log(w0))
        full_effect=sum(w*e for w,e in zip(weights,e_values))
        close(row["full_energy_effect_probability"],full_effect)
        close(row["quotient_energy_effect_probability"],e_values[0])
        gap=full_effect-e_values[0]; close(row["probability_gap"],gap)
        assert all(e>e_values[0] for e in e_values[1:])
        R0=math.sqrt(4*K*K+v*v)
        bound=K*math.tanh(beta*R0)/(6*R0*math.cosh(beta*R0))
        close(row["analytic_single_sector_gap_lower_bound"],bound)
        assert gap>bound>0
        source=sum(w*s for w,s in zip(weights,s_values))
        variance=sum(w*s*s for w,s in zip(weights,s_values))-source*source
        second=sum(w*s for w,s in zip(weights,second_values))-beta*variance
        close(row["full_mean_source"],source)
        close(row["between_sector_source_variance"],variance)
        close(row["full_F_second"],second)
        close(row["sector_mixture_F_second"],second)
        close(row["full_F_second"],-2*source-beta*row["full_KM_covariance"])
        assert row["full_ordinary_variance"]-row["full_KM_covariance"]>1e-4
        assert abs(row["wrong_ordinary_variance_F_second"]-row["full_F_second"])>1e-4
        assert not row["support_regularization_used"]
        for key,value in row["errors"].items():
            if key.startswith("finite_difference"):
                maximum_derivative_crosscheck_error=max(maximum_derivative_crosscheck_error,value)
                assert value<3e-6
            else:
                maximum_identity_error=max(maximum_identity_error,value)
                assert value<2e-11
    assert cases[0]["probability_gap"]>float(signal)

    assert set(fresh["historical_source_sha256"])==HISTORICAL
    for name,digest in fresh["historical_source_sha256"].items():
        assert sha(BASE/name)==digest,name
    assets=[NOTE]+[HERE/name for name in OWN]
    assert len(assets)==8
    assert all(path.name not in {"README.md","RESEARCH_STATE.md","research_direction.md"} for path in assets)
    count=0
    for path in assets:
        assert path.is_file(),path
        if path.suffix==".md":
            text=path.read_text("utf-8-sig")
            assert text.count("$$")%2==0,path
            for target,local in local_links(text):
                resolved=(path.parent/local.replace("\\","/")).resolve()
                assert resolved.exists() or (prospective and resolved==OUT),(path,target)
                count+=1
    assert all(f"## {i}." in NOTE.read_text("utf8") for i in range(1,11))
    assert "独立科学签审通过" in (HERE/"review.md").read_text("utf8")
    return dict(round=1021,date="2026-10-08",all_delivery_checks_passed=True,
        scientific_result_reproduced=True,new_calibration_groups=1,cumulative_research_groups=3799,
        local_links_checked=count,frozen_current_files=len(assets),historical_input_files=len(HISTORICAL),
        live_navigation_frozen=False,neighboring_round_frozen=False,
        cover_electric_configurations_enumerated=1728,Gauss_cover_dimension=12,Gauss_quotient_dimension=2,
        kernel_cycle_sectors=6,tree_cycle_sectors=1,thermal_cases=6,
        reference_Choi_histories=8,exact_probability_gap_lower=str(signal),
        quotient_supported_process_intertwines=True,unknown_quantum_inputs_retained=True,
        center_twirl_keeps_sector_probabilities=True,selective_reweighting_distinguished=True,
        Kubo_covariance_not_ordinary_variance=True,
        maximum_thermal_identity_error=maximum_identity_error,
        maximum_derivative_crosscheck_error=maximum_derivative_crosscheck_error,
        derivatives_are_equilibrium_family_not_fixed_weight_realtime=True,
        physical_thermalization_or_instrument_built=False,all_SM_spectrum_computed=False,
        all_continuous_bundle_sectors_proved=False,global_group_uniquely_generated=False,
        new_cognitive_axiom=False,goal_completed=False,visual_checks_performed=False,
        analytic_scope=fresh["scope"],historical_source_sha256=fresh["historical_source_sha256"],
        source_sha256={str(path.relative_to(ROOT)).replace("\\","/"):sha(path) for path in assets})


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true")
    args=parser.parse_args(); result=verify(prospective=args.write)
    if args.write:
        serialized=json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n"
        with OUT.open("x",encoding="utf8") as stream:
            stream.write(serialized)
    else:
        assert result==json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")},ensure_ascii=False,indent=2))
