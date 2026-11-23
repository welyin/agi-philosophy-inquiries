"""Check weak-basis covariants and the determinant phase of the two frozen mass matrices."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE / path).read_bytes())


def unpack(values):
    array = np.array(values, dtype=float)
    return array[..., 0] + 1j * array[..., 1]


def covariants(A, B, C, projector):
    bridge = C @ np.linalg.solve(A, B)
    metric = A.conj().T @ A
    covectors = [C @ np.linalg.matrix_power(metric, j) @ A.conj().T @ B for j in range(3)]
    tagged = C @ A.conj().T @ projector @ B
    return bridge, covectors, tagged


def unitary(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)))
    return q @ np.diag(np.exp(-1j * np.angle(np.diag(r))))


def matrix_test(w, cfg):
    mass = unpack(w["weak_basis_mass_GeV"])
    vev = cfg["v_GeV"] / math.sqrt(2)
    A, B, C = mass[:3, :3]/vev, mass[:3, 3]/vev, mass[4, :3]/vev
    s, d = mass[3, 3], mass[4, 4]
    hard_phase = (s*d)/abs(s*d)
    projector = np.diag([0., 0., 1.])
    bridge, covectors, tagged = covariants(A, B, C, projector)
    metric_eig = np.linalg.eigvalsh(A.conj().T @ A)
    # Cayley-Hamilton reproduces the inverse using three allowed finite polynomials.
    e2 = sum(metric_eig[i]*metric_eig[j] for i,j in [(0,1),(1,2),(0,2)])
    inverse_poly = (covectors[2] - sum(metric_eig)*covectors[1] + e2*covectors[0])
    assert np.isclose(inverse_poly/np.prod(metric_eig), bridge, rtol=3.e-5)
    assert all(abs(z) > 0 for z in covectors)
    corrections = np.array([vev**2*bridge*z.conjugate()/abs(s*d) for z in covectors])
    top = vev**2*bridge*tagged.conjugate()/abs(s*d)
    assert max(abs(corrections.imag)/abs(corrections)) < 2.e-12
    assert top.imag > 1.e-7

    # Arbitrary changes of axes among like gauge representations; regenerate the top projector.
    rng = np.random.default_rng(1638)
    for _ in range(7):
        ul, ur = unitary(rng), unitary(rng)
        uD, uS = np.exp(1j * rng.uniform(-math.pi, math.pi, size=2))
        W, polys, T = covariants(ul @ A @ ur.conj().T, ul @ B * uS.conjugate(),
                                 uD*C @ ur.conj().T, ul @ projector @ ul.conj().T)
        for v1,v2 in zip([bridge,*covectors,tagged],[W,*polys,T]):
            assert np.isclose(v2,uD*v1*uS.conjugate(),rtol=5.e-10,atol=1.e-15)

    test = cfg["example_real_response"]
    reference = np.linalg.slogdet(mass)[0]
    for z, coefficient in [(covectors[0],corrections[0]),(tagged,top)]:
        new = mass.copy()
        new[3,4] += vev*test*hard_phase*z.conjugate()
        ratio = np.linalg.det(new)/np.linalg.det(mass)
        assert abs(ratio-(1+test*coefficient)) < 5.e-14
        assert abs(np.angle(np.linalg.slogdet(new)[0]/reference)-np.angle(1+test*coefficient))<5.e-15
    # This explicit symmetric interval is necessary and sufficient on the connected branch around 0.
    tan = math.tan(cfg["working_angle"])
    radius = tan/(abs(top.imag)+abs(top.real)*tan)
    for sign in (-1,1):
        assert abs(np.angle(1+sign*radius*top)) <= cfg["working_angle"]*(1+1e-11)
    assert max(abs(np.angle(1+sign*1.001*radius*top)) for sign in (-1,1)) > cfg["working_angle"]
    neutron=read("../1579/independent_results.json")["neutron_theta_upper_in_declared_matching"]
    neutron_tan=math.tan(neutron)
    neutron_radius=neutron_tan/(abs(top.imag)+abs(top.real)*neutron_tan)
    return {"weak_bridge_abs":float(abs(bridge)),
            "down_only_lowest_polynomial_abs":float(abs(covectors[0])),
            "tagged_coefficient_real_imag":[float(top.real),float(top.imag)],
            "angle_for_example_real_response_rad":float(np.angle(1+test*top)),
            "beta_radius_for_all_real_signs":radius,
            "beta_radius_using_archived_neutron_transfer":neutron_radius}


def output():
    cfg=read("inputs.json")
    w=read(cfg["matrices"])
    return {key:matrix_test(w[key],cfg) for key in cfg["keys"]}


def rounded(obj):
    if isinstance(obj,dict): return {k:rounded(v) for k,v in obj.items()}
    if isinstance(obj,list): return [rounded(v) for v in obj]
    return float(f"{obj:.10g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=rounded(output())
    if args.check:
        assert result==read("results.json")
    else:
        with (HERE/"results.json").open("x",encoding="utf-8") as out:
            out.write(json.dumps(result,indent=2)+"\n")
    print("1638 PASS: finite covariants, physical flavour reorientation and exact phase")
