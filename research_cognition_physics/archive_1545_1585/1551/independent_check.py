"""Heavy-pole determinant and a non-aligned exact six-state inverse test."""
import argparse
import json
import math
from itertools import permutations
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def source(path):
    return json.loads((HERE.parent / path).read_text(encoding="utf-8-sig"))


def mixing(angles, phase):
    s12, s13, s23 = np.sin(angles)
    c12, c13, c23 = np.cos(angles)
    z = np.exp(1j * phase)
    return np.array([
        [c12*c13, s12*c13, s13/z],
        [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
        [s12*s23-c12*c23*s13*z, -c12*s23-s12*c23*s13*z, c23*c13],
    ])


def mass_volume():
    previous = source("1550/inputs.json")
    yd_ye = np.array(source("1547/inputs.json")["light_scales"]["1000"])
    yu = np.array(previous["up_at_1TeV"])
    nu = source("1370/inputs.json")
    scale = previous["illustrative_common_v_GeV"] / math.sqrt(2)
    pars = previous["ckm_theta12_theta13_theta23_delta"]
    vckm = mixing(pars[:3], pars[3])
    y = np.vstack((yu, yd_ye))
    contracted = {}
    for key, widths in (("central", 0), ("lower", 2)):
        up, down, lepton = (scale * (y[:, 0]-widths*y[:, 1])).reshape(3, 3)
        diag = up + lepton
        vals = [np.prod(diag + down[list(p)]) for p in permutations(range(3))]
        minimal = np.prod(diag+down)
        assert np.isclose(min(vals), minimal)
        with_v = np.diag(diag)+(vckm*down)@vckm.conj().T
        direct = np.linalg.det(with_v).real
        cof_a, cof_d = np.prod(diag)/diag, np.prod(down)/down
        cross = np.outer(cof_a, down) + np.outer(diag, cof_d)
        expanded = np.prod(diag)+np.prod(down)+np.sum(abs(vckm)**2*cross)
        np.testing.assert_allclose(expanded, direct, rtol=2e-13)
        assert minimal < direct < max(vals)
        contracted[key] = {"min_any_ckm_GeV3":float(minimal),
                           "actual_ckm_GeV3":float(direct),
                           "max_any_ckm_GeV3":float(max(vals))}

    volume = contracted["lower"]["min_any_ckm_GeV3"]
    def offsets(order):
        dat = nu[order]
        # Pick simultaneous upper bounds on both offsets from the lightest.
        if order == "NO":
            return np.array([0., dat["solar"][2], dat["atm_abs"][2]])
        return np.array([dat["atm_abs"][2]-dat["solar"][1], dat["atm_abs"][2], 0.])
    rows = {}
    for order in ("NO", "IO"):
        shift = offsets(order)
        def product(x):
            return float(np.prod(np.sqrt(x*x+shift)))
        beta = nu["katrin"]["m_beta_upper_eV"]
        lower_heavy_max = (volume**2 / (product(beta)*1e-27))**(1/3)
        cap = 1e10
        target = volume**2 / cap**3 * 1e27
        lo, hi = 0., 1.
        for _ in range(70):
            mid = (lo+hi)/2
            if product(mid) < target:
                lo = mid
            else:
                hi = mid
        rows[order] = {"heavy_max_min_katrin_GeV":lower_heavy_max,
                       "lightest_min_at_cap_1e10_GeV_in_eV":lo,
                       "zero_mass_allowed": False}
    return {"volume":contracted,"pole_conditions":rows}


def hermitian_basis():
    basis = [np.diag(v).astype(complex) for v in np.eye(3)]
    for a,b in ((0,1),(0,2),(1,2)):
        for val in (1.,1j):
            k = np.zeros((3,3),complex)
            k[a,b], k[b,a] = val, val.conjugate()
            basis.append(k)
    return np.array(basis)


def coordinates(h):
    return np.array([*h.diagonal().real, *(q for i,j in ((0,1),(0,2),(1,2))
                                           for q in (h[i,j].real,h[i,j].imag))])


def exact_small_system():
    # Comparable scales keep all six physical singular values within double precision.
    D = (mixing([.22,.39,.71],.94)*[.4,.8,1.3]) @ mixing([.8,.2,.5],1.43).conj().T
    U = mixing([.62,.31,.78],-1.92) @ np.diag(np.exp(1j*np.array([.2,-.45,.73])))
    mass = np.array([.02,.08,.1])
    m = np.diag(mass)
    T = np.linalg.solve(D,U.conj())
    G = T.conj().T @ T
    h = np.eye(3,dtype=complex)
    basis = hermitian_basis()
    for _ in range(8):
        res = h@h + m@h.conj()@G@h.conj()@m - np.eye(3)
        derivatives = [h@k+k@h + m@(k.conj()@G@h.conj()+h.conj()@G@k.conj())@m
                       for k in basis]
        jac = np.array([coordinates(d) for d in derivatives]).T
        h -= np.einsum("i,ijk->jk",np.linalg.solve(jac,coordinates(res)),basis)
    n = U@h
    s = np.linalg.solve(D,n.conj()@m)
    to_right = s.conj()@m-D.T@n
    right = np.linalg.solve(s.T,to_right.T).T
    full = np.block([[np.zeros((3,3)),D],[D.T,right]])
    light = np.vstack((n,s))
    np.testing.assert_allclose(light.conj().T@light,np.eye(3),atol=2e-13)
    np.testing.assert_allclose(full@light,light.conj()@m,atol=2e-13)
    np.testing.assert_allclose(full,full.T,atol=3e-13)
    left,_,right_h = np.linalg.svd(n)
    np.testing.assert_allclose(left@right_h,U,atol=2e-13)
    singular = np.sort(np.linalg.svd(full,compute_uv=False))
    np.testing.assert_allclose(singular[:3],mass,atol=2e-13)
    np.testing.assert_allclose(np.prod(singular),abs(np.linalg.det(D))**2,rtol=2e-12)
    return {"D_singular":[.4,.8,1.3],"light":mass.tolist(),"heavy":singular[3:].tolist(),
            "least_eigenvalue_h":float(np.linalg.eigvalsh(h)[0]),
            "active_loss_op":float(np.linalg.norm(s.conj().T@s,2)),
            "determinant_ratio":float(np.prod(singular)/abs(np.linalg.det(D))**2)}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.9g}")
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = rounded({"charged_volume":mass_volume(),"finite_exact_test":exact_small_system()})
    path = HERE / "independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2)
            stream.write("\n")
    elif args.check:
        assert result == json.loads(path.read_text(encoding="utf-8"))
        print("1551 independent: physical mass product and exact six-state inverse verified")
    else:
        print(json.dumps(result,indent=2))
