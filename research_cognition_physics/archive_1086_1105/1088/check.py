"""Finite common clock/signal witness; no exact CO or Lorentz counterclaim.

Run without options to compare with the saved result. --write writes results.json.
The classical loop evaluates actual basis-state transitions of the reversible
extension specified in clock_signal_proof.md; it does not allocate that enormous
unitary. The payload part checks the full linear map, including a reference.
"""
from pathlib import Path
from fractions import Fraction as F
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
SPACE = 16384
PERIOD = 4250
STEPS = 1100


def simulate(epsilon):
    velocity = {"A": 0, "B": 3}
    # Same rule for both clocks; velocity is a retained internal register.
    rates = {observer: 25 + epsilon*w*w for observer, w in velocity.items()}
    phase = {k: 0 for k in rates}
    elapsed = {k: 0 for k in rates}
    emitted = {k: 0 for k in rates}
    position = {"A": (0, 0, 0), "B": (340, 0, 0)}
    active = []
    emissions = []
    arrivals = []
    max_coordinate = 340
    for tick in range(STEPS + 1):
        if tick:
            b = position["B"]
            position["B"] = ((b[0] + velocity["B"]) % SPACE, b[1], b[2])
            for observer in rates:
                phase[observer] = (phase[observer] + rates[observer]) % PERIOD
                elapsed[observer] += rates[observer]
            for pulse in active:
                pulse["position"] += pulse["step"]
                # This check certifies that the selected history never wraps.
                assert 0 <= pulse["position"] < SPACE
            for pulse in list(active):
                receiver = pulse["receiver"]
                if pulse["position"] == position[receiver][0]:
                    arrivals.append({
                        "pulse": pulse["id"], "receiver": receiver,
                        "tick": tick, "clock_units": elapsed[receiver],
                        "clock_cycles": str(F(elapsed[receiver], PERIOD)),
                    })
                    active.remove(pulse)
        # Sending is driven by the local clock and the finite source counter.
        # No analytically predicted emission/arrival time is used here.
        for sender in rates:
            if phase[sender] == 0 and emitted[sender] < 3:
                receiver = "B" if sender == "A" else "A"
                identifier = sender + str(emitted[sender])
                active.append({"id": identifier, "receiver": receiver,
                               "position": position[sender][0],
                               "step": 5 if sender == "A" else -5})
                emissions.append({"pulse": identifier, "tick": tick,
                                  "clock_cycles": str(F(elapsed[sender], PERIOD))})
                emitted[sender] += 1
        max_coordinate = max(max_coordinate, position["B"][0],
                             *(p["position"] for p in active), 0)
    assert not active and len(arrivals) == 6
    assert len({a["pulse"] for a in arrivals}) == 6
    for sender, step in [("A", 170), ("B", 170 if epsilon == 0 else 125)]:
        assert [a["tick"] for a in emissions if a["pulse"].startswith(sender)] == [step*k for k in range(3)]
    expected = {"A": [68 + (272 if epsilon == 0 else 200)*k for k in range(3)],
                "B": [170 + 425*k for k in range(3)]}
    intervals = {}
    for receiver in rates:
        subset = [a for a in arrivals if a["receiver"] == receiver]
        assert [a["tick"] for a in subset] == expected[receiver]
        readings = [F(a["clock_cycles"]) for a in subset]
        differences = [readings[i+1]-readings[i] for i in range(2)]
        assert differences[0] == differences[1]
        intervals[receiver] = differences[0]
    assert intervals["A"] == (F(8, 5) if epsilon == 0 else F(20, 17))
    assert intervals["B"] == (F(5, 2) if epsilon == 0 else F(17, 5))
    return {
        "epsilon": epsilon, "emissions": emissions, "arrivals": arrivals,
        "received_intervals_cycles": {k: str(v) for k, v in intervals.items()},
        "reciprocity_gap": str(abs(intervals["A"] - intervals["B"])),
        "largest_coordinate": max_coordinate, "spatial_wraps": 0,
        "pulse_deliveries": 6, "copying_unknown_payload": False,
    }


def payload_check():
    # Source S, empty receiver D, passive reference R, each a qubit.
    inject = np.zeros((8, 4), complex)
    swap = np.zeros((8, 8), complex)
    for s in range(2):
        for r in range(2):
            inject[4*s+r, 2*s+r] = 1
        for d in range(2):
            for r in range(2):
                swap[4*d+2*s+r, 4*s+2*d+r] = 1
    assert np.array_equal(swap.conj().T @ swap, np.eye(8))
    paulis = [np.array([[0, 1], [1, 0]], complex),
              np.array([[0, -1j], [1j, 0]], complex),
              np.diag([1., -1.]).astype(complex)]
    points = [(0., 0., 0.), (.4, 0., 0.), (-.4, 0., 0.),
              (0., .5, 0.), (0., 0., -.5), (.3, -.2, .1)]
    errors = []
    for xyz in points:
        v = np.array(xyz)
        v = np.diag([1., .75, .5]) @ (v/max(1., np.linalg.norm(v)))
        effect = (np.eye(2)+sum(v[j]*paulis[j] for j in range(3)))/2
        roots = []
        for mat in [effect, np.eye(2)-effect]:
            w, u = np.linalg.eigh(mat)
            assert w.min() >= -1e-14
            roots.append((u*np.sqrt(np.maximum(w, 0))) @ u.conj().T)
        errors.append(float(np.max(np.abs(sum(k.conj().T@k for k in roots)-np.eye(2)))))
        # All matrix units, not a random sample of unknown states.
        for a in range(4):
            for b in range(4):
                matrix = np.zeros((4, 4), complex)
                matrix[a, b] = 1
                for k in roots:
                    kr = np.kron(k, np.eye(2))
                    branch = kr @ matrix @ kr.conj().T
                    joint = swap @ inject @ branch @ inject.conj().T @ swap.conj().T
                    output = np.trace(joint.reshape(2, 4, 2, 4), axis1=0, axis2=2)
                    errors.append(float(np.max(np.abs(output-branch))))
    assert max(errors) < 1e-12
    return {"effects": len(points), "full_input_matrix_units_per_effect": 16,
            "branches_per_input": 2, "swap_unitary_exact": True,
            "reference_included": True, "max_residual": max(errors),
            "note": "General arbitrary-reference statement follows from tensor identity, not finite sampling."}


def main():
    result = {
        "schema": "round1088_finite_common_clock_signal_v1",
        "scope": "one finite quantized task; full finite quantum carrier exists; exact CO not certified",
        "spatial_ring_per_axis": SPACE, "clock_modulus": PERIOD,
        "last_tick": STEPS, "tick_in_tau_units": "1/170",
        "spatial_step_unit": "1/850", "branches": [simulate(0), simulate(1)],
        "quantum_payload": payload_check(),
        "exact_FUCP_plus_CO_counterexample": False,
        "continuous_between_tick_causal_speed_certified": False,
        "full_large_unitary_numerically_allocated": False,
        "scientific_calibration_increment": 0,
    }
    path = HERE / "results.json"
    if "--write" in sys.argv:
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf8")
    else:
        saved = json.loads(path.read_text(encoding="utf8"))
        # Floating residual is diagnostic, all rational records match exactly.
        current_error = result["quantum_payload"].pop("max_residual")
        saved_error = saved["quantum_payload"].pop("max_residual")
        assert current_error < 1e-12 and saved_error < 1e-12
        assert result == saved
    print(json.dumps({"passed": True, "histories": 2, "deliveries": 12,
                      "full_payload_basis_checked": True,
                      "exact_joint_goal_proved": False}))


if __name__ == "__main__":
    main()
