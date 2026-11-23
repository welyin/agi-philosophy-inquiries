"""Triangular orthogonality as a CP test of a real potential, independent of the author construction."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent


def potential(p, q):
    return (p[0]-p[1]-p[2])**2+(p[1]-p[2])**2+sum(
        (q[i]-q[j])**2 for i in range(3) for j in range(i+1, 3))


def projectors(h):
    spectrum = np.linalg.eigvalsh(h)
    unit = np.eye(3, dtype=complex)
    modes = []
    for i, x in enumerate(spectrum):
        projector = unit.copy()
        for j, y in enumerate(spectrum):
            if j != i:
                projector = projector @ ((h-y*unit)/(x-y))
        modes.append(projector)
    return modes


def energy(a, b):
    rows, columns = projectors(a), projectors(b)
    overlaps = np.array([[np.trace(e@f).real for f in columns] for e in rows])
    return potential(overlaps[:, 0], overlaps[:, 1]), overlaps


def finite_gap(eta):
    """eta is a bound on sqrt(U/kappa); these separate factor bounds are conservative."""
    p1low, p1high = (1-eta)/2, (1+eta)/2
    other_low, other_high = (1-p1high-eta)/2, (1-p1low+eta)/2
    qlo, qhi = (1-2*eta)/3, (1+2*eta)/3
    low = np.sqrt([p1low*qlo, other_low*qlo, other_low*qlo])
    high = np.sqrt([p1high*qhi, other_high*qhi, other_high*qhi])
    sides = [low.sum()]
    for k in range(3):
        sides.append(sum(low[i] for i in range(3) if i != k)-high[k])
    assert min(sides) > 0
    return float(math.sqrt(math.prod(sides))/2), sides


def audit():
    a = np.diag([1., 2., 8.])
    pvec = np.array([1/math.sqrt(2), .5, .5])
    qvec = np.exp(1j*np.array([0., 3*np.pi/4, -3*np.pi/4]))/math.sqrt(3)
    last = np.cross(pvec.conjugate(), qvec.conjugate())
    v = np.column_stack((pvec, qvec, last))
    assert np.linalg.norm(v.conj().T@v-np.eye(3)) < 1e-14
    d = np.diag([.3, 2.5, 17.])
    b = v@d@v.conj().T
    e, moduli = energy(a, b)
    conjugate_energy = energy(a, b.conj())[0]
    assert max(e, conjugate_energy) < 1e-25
    j = (v[0, 0]*v[1, 1]*v[0, 1].conjugate()*v[1, 0].conjugate()).imag
    weights = moduli[:, 0]*moduli[:, 1]
    heron = math.sqrt((sum(weights)**2-2*sum(weights**2))/4)
    assert math.isclose(j, 1/12) and math.isclose(heron, abs(j))
    # Independent nonzero basis-invariant of this A, B pair.
    c = a@b-b@a
    cube = np.trace(c@c@c).imag
    assert abs(cube) > 1 and np.isclose(np.trace(c.conj()@c.conj()@c.conj()).imag, -cube)
    rnd = np.random.default_rng(1501)
    for _ in range(16):
        q, _ = np.linalg.qr(rnd.normal(size=(3, 3))+1j*rnd.normal(size=(3, 3)))
        rotated = energy(q@a@q.conj().T, q@b@q.conj().T)
        assert np.allclose(moduli, rotated[1], atol=1e-13)
    threshold, factors = finite_gap(.04)
    data = json.loads((BASE.parent/"1381"/"inputs.json").read_text(encoding="utf-8-sig"))
    ud_min = data["ud_beta"][0] - 2*sum(data["ud_beta"][1:])
    assert np.max(moduli) < ud_min**2
    return {"ground_probabilities": moduli.tolist(), "conjugate_vacuum_J_values": [j, -j],
            "tr_commutator_cubed_imag_example": cube, "U_over_kappa_certified_at_most": .0016,
            "lower_abs_J_entire_sublevel": threshold, "lower_Heron_factors": factors,
            "tree_beta_two_work_widths_Vud_squared_lower": ud_min**2,
            "maximum_squared_modulus_under_any_relabelling_ground": .5}


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    m = p.add_mutually_exclusive_group()
    m.add_argument("--write-new", action="store_true")
    m.add_argument("--check", action="store_true")
    args = p.parse_args()
    path = BASE/"independent_results.json"
    out = json.loads(json.dumps(audit(), default=float))
    if args.write_new:
        with path.open("x", encoding="utf-8") as f:
            f.write(json.dumps(out, indent=2)+"\n")
    elif args.check:
        assert out == json.loads(path.read_text(encoding="utf-8"))
        print("Independent tests passed.")
    else:
        print(json.dumps(out, indent=2))
