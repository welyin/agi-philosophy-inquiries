"""Full chiral Dirac resolvent, static matching and finite-endpoint pole audit."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
SPEC = runpy.run_path(str(ARCH / "1661/independent_check.py"))


def massive(z, m, g):
    return np.array([[m, g*z], [-g*z.conjugate(), m]], complex)


def chiral_inverse(b, p):
    """Euclidean inverse at a momentum along a single axis; spin degeneracy omitted."""
    zero = np.zeros_like(b)
    h = np.block([[zero, b], [b.conj().T, zero]])
    return np.linalg.inv(h+1j*p*np.eye(len(h)))


def transfer(z, m, g, p):
    # light L couples R1, light R couples L2. The effective mass adds a minus sign.
    return -chiral_inverse(massive(z, m, g), p)[2, 1]


def full(z, m, g, t):
    mat = np.zeros((3, 3), complex)
    mat[1:, 1:] = massive(z, m, g)
    mat[0, 1] = mat[2, 0] = t
    return mat


def eigenvalues(zs, m, g, t):
    return [float(np.linalg.svd(full(z, m, g, t), compute_uv=False)[-1]) for z in zs]


def main():
    phi = -math.pi/6
    z = (3+2*np.exp(4j*phi))/5
    q = abs(z)
    legs = [1., np.exp(2j*phi), np.exp(3j*phi)*z]
    g = 1.
    m = math.sqrt(q)
    momenta = np.array([0., .2, .5, 1., 3.])
    response = np.array([[transfer(a, m, g, p) for a in legs] for p in momenta])
    for i, p in enumerate(momenta):
        expected = np.array([g*a/(m*m+g*g*abs(a)**2+p*p) for a in legs])
        np.testing.assert_allclose(response[i], expected, atol=5.e-16)
    ratio = abs(response[:, 2]/response[:, 0])
    np.testing.assert_allclose(ratio, 1-(1-q)*momenta**2/(q*(1+q)+momenta**2))
    np.testing.assert_allclose(response[0]/response[0, 0], np.array(legs)/abs(np.array(legs)))
    heavy = [np.linalg.svd(massive(a, m, g), compute_uv=False) for a in legs]
    np.testing.assert_allclose(heavy[2]/heavy[0], [math.sqrt(q)]*2)
    # Use an independently diagonalised neutral source and charged SVD.
    d=SPEC["read"]("1657/independent_results.json")["extended_box_GJ_test"]
    e, _ = SPEC["OLD"]["charged"](d)
    _, neutral, _ = SPEC["OLD"]["mass_source"](SPEC["roots"](5.74)[0], 1, 7.53e-5)
    angles = []
    mass_ratio = []
    for row in response:
        k = np.diag(row/row[0])
        n = k @ neutral @ k.T
        ev, v = np.linalg.eigh(n.conj().T@n)
        mix = abs(e.conj().T@v)**2
        pe = mix[0, 2]
        angles.append([4*pe*(1-pe), mix[0, 1]/(1-pe), mix[1, 2]/(1-pe)])
        mass_ratio.append(math.sqrt(ev[2]/ev[1]))
    np.testing.assert_allclose(angles[0], [.080731633658, .30499020104, .55282862668], atol=2.e-11)
    np.testing.assert_allclose(mass_ratio[0], 5.74)
    # Identical finite couplings t on both ends for each family, in units g=1.
    t = .4
    poles = eigenvalues(legs, m, g, t)
    can = [t*t*abs(a)/(m*m+abs(a)**2+t*t) for a in legs]
    for a, pole in zip(legs, poles):
        np.testing.assert_allclose(pole*(m*m+abs(a)**2+t*t-pole**2), t*t*abs(a))
    assert abs(poles[2]/poles[0]-1) > .04
    assert t*t < q*(1+q)
    target = t*t/(1+q)
    tuned_m = math.sqrt(q-t*t+target**2)
    retuned = eigenvalues(legs, tuned_m, g, t)
    np.testing.assert_allclose(retuned, target, atol=2.e-16)
    eps = .05
    return dict(q=q, heavy_mass_ratio=math.sqrt(q), momenta_over_g=momenta.tolist(),
                relative_modulus=ratio.tolist(), fixed_seed_mass_ratio=mass_ratio,
                fixed_seed_angles=angles,
                five_percent_limit_p_over_g=math.sqrt(eps*q*(1+q)/(1-q-eps)),
                finite_end_t_over_g=t, static_wilson_mass_over_g=t*t/(1+q),
                canonical_leading_masses_over_g=can, actual_light_poles_over_g=poles,
                actual_pole_ratio=poles[2]/poles[0],
                retuned_M_squared_over_g_squared=tuned_m**2, retuned_light_poles_over_g=retuned)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    result = main()
    target = HERE/"independent_results.json"
    if args.check:
        old = json.loads(target.read_bytes())
        for name, value in result.items():
            np.testing.assert_allclose(value, old[name], atol=1.e-12, rtol=1.e-11)
        print("PASS: complete Dirac inverse, light poles, canonical residues and rank-two mass")
    else:
        with target.open("x", encoding="utf8") as f:
            f.write(json.dumps(result, indent=2)+"\n")
