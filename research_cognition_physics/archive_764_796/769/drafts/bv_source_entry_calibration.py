"""769 working entry: original-potential BV contact and gauge-fixing jets.

These are finite algebra calibrations, NOT quantum anomalies, loop integrals,
or a completed research round. Historical modules are imported read-only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_brst_relative_source as old

TARGET = HERE / "bv_source_entry_calibration_results.json"


def canonical_contact():
    phi = np.array([.13, .64, -.11, .08, .37])
    shift = np.array([.27, -.14, .31, .18, -.09])
    direction = np.array([-.17, .12, .07, -.21, .16])
    t = old.old.scalar_generators()[8]
    e, j, _ = old.potential_gradient_jet(phi, shift)
    _, h_direction, q_direction = old.potential_gradient_jet(phi, direction)
    _, _, q_shift = old.potential_gradient_jet(phi, shift)
    _, _, q_sum = old.potential_gradient_jet(phi, shift + direction)
    # For B' = B + hbar * shift: Gamma_1 = E . shift,
    # J = Hessian * shift, R_1 = t * shift. The cubic mixed derivative
    # is obtained by exact Taylor-jet polarization, not finite differencing.
    r, r1 = t @ phi, t @ shift
    dj = q_sum - q_shift - q_direction
    bare = float(r @ j)
    contact = float(e @ r1)
    derivative_bare = float((t @ direction) @ j + r @ dj)
    derivative_contact = float(h_direction @ r1)
    errors = {
        "source_Ward_with_generator_contact": abs(bare + contact),
        "differentiated_Ward_with_contact": abs(derivative_bare + derivative_contact),
    }
    assert max(errors.values()) < 2e-14
    assert abs(bare) > 1e-6 and abs(derivative_bare) > 1e-6
    return dict(errors=errors, bare_source_contraction=bare,
                generator_contact=contact, bare_response_contraction=derivative_bare,
                response_contact=derivative_contact,
                scope="Off-shell original H5 potential under a local coordinate translation. This is a BV canonical-coordinate calibration, not an allowed physical deformation, actual quantum counterterm, or computed anomaly.")


def scalar_gauge_fixing_jet():
    phi = np.array([.13, .64, -.11, .08, .37])
    direction = np.array([-.17, .12, .07, -.21, .16])
    chi = np.array([.21, -.08, .19, .06, -.14])
    t = old.old.scalar_generators()[8]
    f = 2 - phi @ phi / 6
    dot = phi @ direction
    metric = old.old.scalar_metric(phi)
    dmetric = (np.eye(5) * dot / (3*f*f)
               + (np.outer(direction, phi) + np.outer(phi, direction)) / (6*f*f)
               + np.outer(phi, phi) * dot / (9*f**3))
    direct = float((t @ direction) @ metric @ chi + (t @ phi) @ dmetric @ chi)
    # Antisymmetry gives (t phi)^T K chi = (t phi . chi)/F.
    reduced = float((t @ direction) @ chi / f + (t @ phi) @ chi * dot / (3*f*f))
    value = float((t @ phi) @ metric @ chi)
    error = abs(direct - reduced)
    assert error < 2e-14 and abs(direct) > 1e-5
    return dict(derivative_error=error, scalar_gauge_fixing_coefficient=value,
                background_derivative=direct,
                scope="One actual scalar term of K-star, paired with the auxiliary bundle, on a declared off-shell local jet. Nonzero variation is not a quantum split anomaly; full metric/gauge/derivative terms are not simulated.")


def run():
    dependencies = ("joint_brst_relative_source.py", "joint_covariant_gauge_complex.py",
                    "joint_scalar_propagation_matching.py", "research_note_600.md",
                    "research_note_629.md", "research_note_735.md", "research_note_768.md")
    return dict(working_round=769, completed_round=768, numbered_scientific_tests_added=0,
                entry_calibrations=2, failures=0, canonical_contact=canonical_contact(),
                scalar_gauge_fixing_jet=scalar_gauge_fixing_jet(),
                dependency_hashes={n: hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest()
                                   for n in dependencies},
                all_sector_quantum_Ward_proven=False,
                background_split_anomaly_computed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open("x", encoding="utf8") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2)
    else:
        assert result == json.loads(TARGET.read_text("utf8"))
    print(json.dumps(result, ensure_ascii=False, indent=2))
