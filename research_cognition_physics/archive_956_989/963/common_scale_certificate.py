"""963: finite common clock-scale consistency, exact rational certificates.

The records below are synthetic declared contracts, not experimental data.
Difference constraints are established mathematics, not a new algorithm.
No clock hardware, spacetime, causal propagation, or gravity is generated.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
ROOT = HERE.parents[2]
RESEARCH = STAGE.parent
TARGET = HERE/"common_scale_certificate_results.json"
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
BASE = [F(0), F(3,100), F(5,100), F(2,100)]


def differences(s):
    return [s[v]-s[u] for u, v in EDGES]


def make_records(groups, eps):
    return [dict(clock=clock, u=u, v=v, b=b, eps=F(eps))
            for clock, values in groups.items()
            for (u, v), b in zip(EDGES, values)]


def constraints(records):
    arcs = []
    for k, r in enumerate(records):
        arcs.append(dict(u=r["u"], v=r["v"], weight=r["b"]+r["eps"],
                         record=k, kind="upper", clock=r["clock"]))
        arcs.append(dict(u=r["v"], v=r["u"], weight=-r["b"]+r["eps"],
                         record=k, kind="reverse_lower", clock=r["clock"]))
    return arcs


def certify(records, nodes=4):
    """Bellman-Ford from an implicit zero-weight source to every node."""
    arcs = constraints(records)
    dist = [F(0) for _ in range(nodes)]
    pred = [None]*nodes
    for _ in range(nodes):
        changed = None
        for index, arc in enumerate(arcs):
            candidate = dist[arc["u"]] + arc["weight"]
            if dist[arc["v"]] > candidate:
                dist[arc["v"]] = candidate
                pred[arc["v"]] = index
                changed = arc["v"]
        if changed is None:
            potential = [x-dist[0] for x in dist]
            residuals = [potential[r["v"]]-potential[r["u"]]-r["b"] for r in records]
            assert all(abs(x) <= r["eps"] for x, r in zip(residuals, records))
            return dict(feasible=True, potential=potential, residuals=residuals)
    # Follow predecessors through n steps to enter a repeated directed cycle.
    vertex = changed
    for _ in range(nodes):
        assert pred[vertex] is not None
        vertex = arcs[pred[vertex]]["u"]
    start = vertex
    cycle = []
    while True:
        arc = arcs[pred[vertex]]
        cycle.append(arc)
        vertex = arc["u"]
        if vertex == start:
            break
        assert len(cycle) <= nodes
    cycle.reverse()
    assert all(cycle[i]["v"] == cycle[(i+1) % len(cycle)]["u"]
               for i in range(len(cycle)))
    weight = sum((a["weight"] for a in cycle), F(0))
    assert weight < 0
    return dict(feasible=False, negative_cycle=cycle, total_cycle_weight=weight)


def serial(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k:serial(v) for k,v in x.items()}
    if isinstance(x, list):
        return [serial(v) for v in x]
    return x


def run():
    base = differences(BASE)
    changed = BASE.copy()
    changed[2] += F(1,50)
    species = {"A":base, "B":differences(changed)}
    twisted = base.copy()
    twisted[3] += F(1,25)
    common_twist = {"A":twisted, "B":twisted}
    cases = [
        ("common_exact", {"A":base, "B":base}, F(0), True),
        ("each_species_exact_but_no_common_scale", species, F(0), False),
        ("species_tolerance_below_boundary", species, F(9,1000), False),
        ("species_tolerance_at_boundary", species, F(1,100), True),
        ("common_path_dependence", common_twist, F(0), False),
        ("path_tolerance_below_boundary", common_twist, F(1,100), False),
        ("path_tolerance_at_boundary", common_twist, F(1,75), True)]
    reports = []
    for name, groups, eps, expected in cases:
        records = make_records(groups, eps)
        answer = certify(records)
        assert answer["feasible"] is expected, name
        reports.append(dict(name=name, records=records, certificate=answer))
    # The species conflict is not an individually inconsistent apparatus.
    individual = {a:certify(make_records({a:values}, 0)) for a,values in species.items()}
    assert all(r["feasible"] for r in individual.values())
    # Exact lower bounds: a pair of species differs by 1/50, shared intervals
    # need eps >=1/100. The shortcut triangle has mismatch 1/25 and 3 edges.
    assert abs(species["A"][1]-species["B"][1]) == F(1,50)
    cycle0123 = sum(twisted[:4])
    cycle023 = twisted[4]+twisted[2]+twisted[3]
    assert cycle0123 == cycle023 == F(1,25)
    # Same local unit change for every species leaves feasibility invariant.
    chi = [F(7,10), F(-1,5), F(3,10), F(-4,5)]
    gauge = differences(chi)
    for name, groups, eps, expected in cases:
        transformed = {a:[b+d for b,d in zip(bs,gauge)] for a,bs in groups.items()}
        assert certify(make_records(transformed, eps))["feasible"] is expected
    exact_shift = certify(make_records({"A":[b+d for b,d in zip(base,gauge)]}, 0))
    assert exact_shift["potential"] == [BASE[i]+chi[i]-chi[0] for i in range(4)]
    # A check of the dimensionless observations after reconstructing scales.
    witness = certify(make_records({"A":base, "B":base}, 0))["potential"]
    scales = np.exp(np.array([float(x) for x in witness]))
    rate_error = max(abs(math.log(scales[v]/scales[u])-float(b))
                     for (u,v),b in zip(EDGES,base))
    assert rate_error < 3e-15
    # Flat connection on a ring need not be globally exact: all local cells
    # absent/contractible curvature zero, but total scale holonomy nonzero.
    ring_cycle = F(1,10)
    ring_records = [dict(clock="A",u=i,v=(i+1)%4,b=ring_cycle/4,eps=F(0))
                    for i in range(4)]
    ring = certify(ring_records)
    assert not ring["feasible"]
    # Rate equality never requires elapsed proper times to coincide.
    # Twin protocol: two segments at +/-v, total inertial-coordinate time 2.
    # v=3/5, Lorentz gamma=5/4 => elapsed 8/5 versus 2, later rest rates equal.
    v, coordinate_time = F(3,5), F(2)
    lapse_squared = 1-v*v
    lapse = F(math.isqrt(lapse_squared.numerator), math.isqrt(lapse_squared.denominator))
    assert lapse*lapse == lapse_squared
    stationary_elapsed = coordinate_time
    moving_elapsed = coordinate_time*lapse
    twin = dict(v=str(v), elapsed_stationary=str(stationary_elapsed),
                elapsed_traveling=str(moving_elapsed), subsequent_rate_ratio="1",
                is_second_clock_effect=False,
                ideal_instant_turnaround_and_identical_rest_clocks_are_inputs=True)
    assert moving_elapsed == F(8,5) and stationary_elapsed-moving_elapsed == F(2,5)
    out = dict(round=963,date="2026-10-07",all_scientific_checks_passed=True,
        kind="finite conditional calibration certificate, not a spacetime model",
        cases=reports,individual_species=individual,
        sharp_uniform_error_thresholds=dict(species="1/100",loop="1/75",
            proof="cycle lower bounds plus exact feasible potentials at equality"),
        shared_unit_change_invariant=True,scale_reconstruction_log_error=rate_error,
        flat_nonexact_ring=dict(total_log_holonomy=ring_cycle,certificate=ring),
        rate_vs_elapsed_time=twin,
        scope=dict(synthetic_declared_records=True,
            errors_are_measurement_contracts_not_physical_cutoffs=True,
            no_clock_material_dynamics_or_source_derived=True,
            zero_scale_holonomy_does_not_remove_phase_or_spacetime_curvature=True,
            EPS_requires_additional_geometry_and_clock_identification=True,
            not_all_Weyl_matter_theories_ruled_out=True,
            no_continuum_or_dimension_or_Einstein_equation_proved=True,
            full_goal_completed=False,stop_calibration_tool_optimization=True))
    sources = [Path(__file__), HERE/"drafts/calibration_bridge_decision.md",
               RESEARCH/"archive_259_300/research_note_290.md",
               RESEARCH/"archive_301_341/research_note_324.md",
               RESEARCH/"archive_342_369/research_note_352.md",
               STAGE/"research_note_930.md",STAGE/"research_note_933.md",
               STAGE/"research_note_941.md",STAGE/"research_note_962.md"]
    out["source_hashes"] = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                            for p in sources}
    return serial(out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as dest:
            json.dump(result,dest,ensure_ascii=False,indent=2)
            dest.write("\n")
    else:
        saved = json.loads(TARGET.read_text("utf-8"))
        assert abs(saved.pop("scale_reconstruction_log_error") -
                   result["scale_reconstruction_log_error"]) < 1e-14
        current = {k:v for k,v in result.items() if k != "scale_reconstruction_log_error"}
        assert saved == current
    print(json.dumps(dict(round=963,all_scientific_checks_passed=True,
        outcomes={c["name"]:c["certificate"]["feasible"] for c in result["cases"]},
        sharp_uniform_error_thresholds=result["sharp_uniform_error_thresholds"],
        scope=result["scope"]),ensure_ascii=False,indent=2))
