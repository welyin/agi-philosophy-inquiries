"""Working calibration, not a completed scientific round or a Hadamard theorem."""
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE/"smooth_positivity_calibration_results.json"


def run():
    rows = []
    for n in (8, 16, 24):
        s = np.exp(-np.sqrt(1+np.arange(1, n+1, dtype=float)**2))
        eye, zero = np.eye(n), np.zeros((n, n))
        S = np.diag(s); noise = np.block([[S, zero], [zero, S]])
        Q = np.block([[eye, zero], [zero, -eye]])
        plus = np.block([[eye-S, zero], [zero, -S]])
        minus = np.block([[-S, zero], [zero, eye-S]])
        swap = np.block([[zero, eye], [eye, zero]])
        error = float(np.max(abs(plus-minus-Q)))
        reality = float(np.max(abs(swap@plus@swap-minus)))
        corrected = [float(np.linalg.eigvalsh(x+noise).min()) for x in (plus, minus)]
        # Correct the first four modes only: another negative direction survives.
        limited = np.zeros_like(noise)
        limited[np.arange(4), np.arange(4)] = s[:4]
        limited[n+np.arange(4), n+np.arange(4)] = s[:4]
        remaining = float(np.linalg.eigvalsh(plus+limited).min())
        assert error < 1e-14 and reality == 0 and min(corrected) >= -1e-14
        assert remaining < -1e-4
        rows.append(dict(modes_per_sign=n, commutator_error=error,
                         conjugacy_error=reality,
                         corrected_min_eigenvalues=corrected,
                         remaining_negative_after_four_modes=remaining))
    return dict(status="working calibration; official completed round remains 765",
                rows=rows, all_checks_passed=True,
                infinite_rank_claim_basis="Kernel intersection proof in working report, not finite truncation.",
                actual_coupled_Hadamard_state_constructed=False)


if __name__ == "__main__":
    result = run()
    if TARGET.exists():
        assert json.loads(TARGET.read_text("utf8")) == result
    else:
        with TARGET.open("x", encoding="utf8") as f:
            json.dump(result, f, ensure_ascii=False, indent=2)
    print(json.dumps(dict(all_checks_passed=True, completed_round=765)))
