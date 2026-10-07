"""962: reuse 430/438; screen clock protection against relation feedback.

No spatial geometry, physical gravity or internal readout hardware is derived.
Default is read-only recomputation; --write creates the result exclusively.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent
TARGET = HERE / "clock_relation_screen_results.json"
OLD = RESEARCH / "archive_429_466"
spec = importlib.util.spec_from_file_location("round430_reuse", OLD / "430/exchange_relation_audit.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def norm(a):
    return float(np.linalg.norm(a, 2))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_close(a, b, tol=5e-12):
    e = norm(a - b)
    assert e < tol, e
    return e


def run():
    i2, i4, i8 = np.eye(2), np.eye(4), np.eye(8)
    z = np.diag([1., -1.])
    # Same raw three spins as 430, graph sector exactly the N=3, C=1
    # resource-conserving sector of 438. No newly appended clock qubit.
    a, b, c, p, _, _, _ = old.operators()
    w = old.encoding()
    jz = sum(np.kron(np.kron(z if k == 0 else i2,
                            z if k == 1 else i2),
                     z if k == 2 else i2) for k in range(3)) / 2
    n = np.diag([0., 1., 1., 1.])
    f = np.zeros((4, 4))
    f[0, 1:] = 1
    f[1:, 0] = 1
    swaps = (a, b, c)
    logical = (-z, (np.sqrt(3)*old.PAULI[0] + z)/2,
               (-np.sqrt(3)*old.PAULI[0] + z)/2)
    structural = {}
    structural["encoding"] = check_close(w.conj().T @ w, i4)
    structural["clock_on_same_spins"] = check_close(jz @ w, w @ np.kron(z/2, i2))
    for k, (raw, log) in enumerate(zip(swaps, logical)):
        structural[f"swap_{k}"] = check_close(raw @ w, w @ np.kron(i2, log))
    h0_raw = np.kron(i8, f)
    h0 = np.kron(i2, f)
    for e, (raw, log) in enumerate(zip(swaps, logical), 1):
        proj = np.zeros((4, 4))
        proj[e, e] = 1
        h0_raw += np.kron(i8 + raw.real, proj)
        h0 += np.kron(i2 + log.real, proj)
    wc = np.kron(w, i4)
    structural["old438_sector_intertwiner"] = check_close(
        h0_raw @ wc, wc @ np.kron(i2, h0))
    # Exact local resource charges f_i + degree_i = 1 in every legal graph.
    degrees = np.array([[0, 1, 0, 1], [0, 1, 1, 0], [0, 0, 1, 1]])
    resources = 1 - degrees
    assert np.array_equal(resources + degrees, np.ones((3, 4)))

    h_clock = np.kron(z/2, np.eye(8))
    h_conversion = np.kron(i4, f)
    h_rest = np.kron(i2, h0)
    source = np.kron(np.kron(z/2, i2), n)
    results = []
    unitary_errors = []
    for eta in (0., .5):
        raw = h0_raw + np.kron(jz, i4 + eta*n)
        h = h_rest + h_clock + eta*source
        structural[f"full_intertwiner_eta_{eta}"] = check_close(raw @ wc, wc @ h)
        assert norm(raw @ np.kron(p, i4) - np.kron(p, i4) @ raw) < 1e-12
        parts = [h_conversion, h_rest - h_conversion, h_clock, eta*source]
        currents = [1j*(h@part - part@h) for part in parts]
        balance = norm(sum(currents))
        assert balance < 1e-12
        for t in (.005, .5):
            u = old.evolve(h, t)
            uraw = old.evolve(raw, t)
            unitary_errors.append(check_close(uraw @ wc, wc @ u))
            check_close(u.conj().T @ u, np.eye(16))
            check_close(u.conj().T @ h @ u, h)
            e0 = np.zeros(8, complex)
            e0[0] = 1  # logical0, empty graph
            # Strip the common +/- omega/2 phase to inspect relation branches.
            hp = h0 + eta/2*np.kron(i2, n)
            hm = h0 - eta/2*np.kron(i2, n)
            ep = old.evolve(hp, t) @ e0
            em = old.evolve(hm, t) @ e0
            overlap = np.vdot(em, ep)
            visibility = float(abs(overlap))
            rqp = old.partial(np.outer(ep, ep.conj()), (2, 4), (1,))
            rqm = old.partial(np.outer(em, em.conj()), (2, 4), (1,))
            dq = old.distance(rqp, rqm)
            # All complement degrees, not only the visible graph occupation.
            dre = old.distance(np.outer(ep, ep.conj()), np.outer(em, em.conj()))
            occp = float(np.trace(rqp @ n).real)
            occm = float(np.trace(rqm @ n).real)
            assert dq <= dre + 2e-12
            assert abs(dre**2 + visibility**2 - 1) < 1e-11
            # Fixed-basis coherence equals bare phase times branch overlap.
            inp = np.kron(np.array([1., 1.])/np.sqrt(2), e0)
            state = u @ inp
            rc = old.partial(np.outer(state, state.conj()), (2, 2, 4), (0,))
            check_close(rc, np.array([[.5, np.exp(-1j*t)*overlap/2],
                                     [np.exp(1j*t)*overlap.conjugate()/2, .5]]))
            reference_clock = old.evolve(z/2, t)
            factor_error = norm(u - np.kron(reference_clock, old.evolve(h0, t)))
            if eta == 0:
                assert factor_error < 1e-12 and dq < 1e-12
            else:
                assert factor_error > 1e-5 and dre > 1e-6
            results.append(dict(eta=eta, time=t, clock_visibility=visibility,
                graph_trace_distance=dq, full_complement_trace_distance=dre,
                occupied_probability_plus=occp, occupied_probability_minus=occm,
                occupied_difference_plus_minus=occp-occm,
                isolated_clock_factorization_error=factor_error,
                total_current_balance_error=balance,
                interaction_current_norm=norm(currents[-1])))
    # The graph responds to unknown relational data even at eta=0.
    probs = []
    for logical_index in (0, 1):
        inp = np.zeros(8, complex)
        inp[4*logical_index] = 1
        psi = old.evolve(h0, .5) @ inp
        rhoq = old.partial(np.outer(psi, psi.conj()), (2, 4), (1,))
        probs.append(float(rhoq[1, 1].real))
    assert abs(probs[0] - probs[1]) > .005
    # Independent exact finite-series certificate for total occupied
    # probability: p_plus - p_minus = -t^4/4 + O(t^6).
    # ||Hbranch|| <= sqrt(3)+2+1/4 < 4; real H and input imply even series.
    t = Fraction(1, 200)
    leading = t**4/4
    remainder = 2*Fraction(25, 24)*(8*t)**6/720  # exp(1/25) <= 25/24
    lower = leading - remainder
    assert lower > 0
    # Check the fourth derivative without fitting a time series.
    obs = np.kron(i2, n).astype(complex)
    fourth = []
    for sign in (1, -1):
        branch = h0 + sign*.25*np.kron(i2, n)
        heis = obs.copy()
        for _ in range(4):
            heis = 1j*(branch@heis - heis@branch)
        fourth.append(float(heis[0, 0].real))
    assert abs(fourth[0] - fourth[1] + 6) < 1e-12
    out = dict(round=962, date="2026-10-07", all_scientific_checks_passed=True,
        structural_errors=structural, raw_and_encoded_unitary_max_error=max(unitary_errors),
        cases=results,
        relation_data_response_at_eta_zero=dict(edge12_probabilities=probs,
            absolute_contrast=abs(probs[0]-probs[1])),
        analytic_nonzero_clock_source_certificate=dict(time=str(t),
            fourth_derivative_difference=fourth[0]-fourth[1],
            signed_leading_term="-t^4/4", absolute_remainder_upper=str(remainder),
            negative_contrast_magnitude_lower=str(lower), lower_float=float(lower)),
        scope=dict(old430_encoding_and_old438_resource_sector_reused=True,
            old449_compatibility_not_reproved=True, same_raw_spins_carry_clock_and_relation=True,
            uniform_clock_field_and_eta_coupling_are_inputs=True,
            clock_preparation_readout_and_axis_are_not_generated=True,
            arbitrary_unknown_encoded_inputs_and_references_in_intertwiner=True,
            source_visibility_example_uses_a_declared_product_preparation=True,
            isolated_clock_requirement_is_stronger_than_H1_H3=True,
            no_metric_no_expansion_no_universal_gravity_derived=True,
            no_autonomous_recovery_apparatus_certified=True,
            full_goal_completed=False, stop_this_candidate_optimization=True))
    inputs = [Path(__file__), OLD/"430/exchange_relation_audit.py",
              OLD/"research_note_430.md", OLD/"research_note_438.md",
              OLD/"research_note_449.md", OLD/"research_note_459.md",
              STAGE/"research_note_961.md", STAGE/"961/operational_scale_screen_results.json"]
    out["source_hashes"] = {str(q.relative_to(ROOT)): sha(q) for q in inputs}
    return out


def compare(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            compare(a[key], b[key])
    elif isinstance(a, list):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            compare(x, y)
    elif isinstance(a, float):
        assert abs(a-b) <= 2e-11*(1+abs(a)), (a, b)
    else:
        assert a == b, (a, b)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = run()
    if args.write:
        with TARGET.open("x", encoding="utf-8") as dest:
            json.dump(out, dest, ensure_ascii=False, indent=2)
            dest.write("\n")
    else:
        compare(out, json.loads(TARGET.read_text("utf-8")))
    print(json.dumps({k:v for k,v in out.items() if k not in ("source_hashes", "structural_errors")},
                     ensure_ascii=False, indent=2))
