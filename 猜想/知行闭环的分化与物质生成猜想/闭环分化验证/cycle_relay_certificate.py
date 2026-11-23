"""Compare cyclic relational readout with certification of an unknown four-dimensional handoff."""
import argparse
import json
from pathlib import Path

import numpy as np

from cycle_memory_readout import CODE, CYCLE, H, POWERS, V, W, code_outputs


def projector(v):
    return np.outer(v, v.conj())


def close(a, b):
    np.testing.assert_allclose(a, b, rtol=0, atol=2e-12)


def mean_channel(operators, matrix):
    return sum(u @ matrix @ u.conj().T for u in operators) / len(operators)


def effect_spectrum(states):
    operator = np.zeros((16, 16))
    for v in states:
        operator += 4 * projector(np.kron(v, v)) / len(states)
    return np.linalg.eigvalsh(operator)


def checks():
    ident = np.eye(2)
    x = np.array([[0, 1], [1, 0]])
    z = np.diag([1, -1])
    p = W.T @ CYCLE @ W
    f = W.T @ np.eye(4)[[0, 2, 1, 3]] @ W
    close(f @ p @ f, p.T)
    common = [np.kron(ident, W.T @ step @ W) for step in POWERS]
    reflection = np.kron(x, f)
    group = common + [u @ reflection for u in common]
    effects = []
    for step in POWERS:
        k = W.T @ step @ W
        effects.append(np.block([[ident, k.T], [k, ident]]) / 3)
    close(sum(effects), np.eye(4))
    for e in effects:
        for u in group:
            close(u.T @ e @ u, e)
    for (a, b), state in code_outputs().items():
        rho = projector(CODE.T @ state)
        for ensemble in (common, group):
            output = mean_channel(ensemble, rho)
            close([np.trace(e @ output) for e in effects],
                  [2/3 if k == (b-a) % 3 else 1/6 for k in range(3)])

    route_plus = projector(H[:, 0])
    for i, j in np.ndindex(2, 2):
        cell = np.eye(2)[:, [i]] @ np.eye(2)[[j], :]
        close(mean_channel(group, np.kron(route_plus, cell)),
              np.kron(route_plus, np.trace(cell)*ident/2))
    omega = np.eye(4).ravel()/2
    fe = {}
    for name, ensemble in (("cyclic", common), ("dihedral", group)):
        output = mean_channel([np.kron(np.eye(4), u) for u in ensemble], projector(omega))
        fe[name] = round(float(omega @ output @ omega), 9)

    routes = [*np.eye(2).T, *H.T]
    spins = [W.T @ step @ V for step in POWERS]
    probes = [np.kron(c, w) for c in routes for w in spins]
    spectrum = effect_spectrum(probes)
    close(spectrum, [0]*7 + [0.25]*4 + [0.5]*4 + [1])
    for u in (np.kron(x, ident), np.kron(z, ident)):
        close(sum(abs(v @ u @ v)**2 for v in probes) / 12, 0.5)
        close(np.trace(u), 0)
    bad_scores = {}
    for name, ensemble in (("cyclic", common), ("dihedral", group)):
        bad_scores[name] = round(float(np.mean(
            [v @ mean_channel(ensemble, projector(v)) @ v for v in probes])), 9)

    # Test leakage out of the four-dimensional code as well as both scalar fields.
    rng = np.random.default_rng(28)
    embed = np.eye(5, 4)
    for complex_entries in (False, True):
        for _ in range(24):
            draw = rng.normal(size=(15, 4))
            if complex_entries:
                draw = draw + 1j * rng.normal(size=(15, 4))
            ops = np.linalg.qr(draw)[0].reshape(3, 5, 4)
            p_test = np.mean([sum(abs(np.vdot(embed @ v, u @ v))**2 for u in ops)
                              for v in probes])
            ent = sum(abs(np.trace(embed.T @ u))**2 for u in ops)/16
            assert ent >= 2*p_test - 1 - 1e-12
            # Mix with the ideal isometry to also test informative, near-unit scores.
            q = 0.003
            score = (1-q) + q*p_test
            witness = min(1, np.sqrt(8*(1-np.sqrt(max(0, 2*score-1)))))
            ket = rng.normal(size=12)
            if complex_entries:
                ket = ket + 1j*rng.normal(size=12)
            ket = ket / np.linalg.norm(ket)
            before = projector(np.kron(np.eye(3), embed) @ ket)
            after = sum(projector(np.kron(np.eye(3), u) @ ket) for u in ops)
            distance = q * np.abs(np.linalg.eigvalsh(after-before)).sum()/2
            assert distance <= witness + 1e-12
    n = 30000
    lower = 0.01**(1/n)
    return {"verified_private_setting_pairs": 9, "product_probe_settings": len(probes),
            "comparison_success_both_relays": "2/3", "entanglement_fidelities": fe,
            "probe_acceptance": bad_scores, "process_spectrum": {"0": 7, "1/4": 4, "1/2": 4, "1": 1},
            "random_channels_with_leakage": 48, "zero_failures_99pct": {
                "trials": n, "p_lower": round(lower, 10),
                "partner_error_upper": round(float(np.sqrt(8*(1-np.sqrt(2*lower-1)))), 10)}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = checks()
    path = Path(__file__).with_name("cycle_relay_certificate_results.json")
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8")) == result
    else:
        with path.open("x", encoding="utf-8") as output:
            output.write(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))
