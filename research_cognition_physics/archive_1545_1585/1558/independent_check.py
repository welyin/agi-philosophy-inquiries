"""Check common Higgs matching with basis invariants and a seesaw-preserving rescaling."""
import argparse
import json
from math import pi
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
VEV = 174.0
HIGGS_POLE_GEV = 125.11


def frob2(matrix):
    return float(np.vdot(matrix, matrix).real)


def positive_threshold(y, right):
    # Right Takagi basis V satisfies right* right = V diag(M^2) V*, with * adjoint.
    return frob2(y @ right.conj().T)/(8*pi**2)


def light_eV(y, right):
    return -(1e9*VEV**2)*y @ np.linalg.solve(right, y.T)


def independently_rebuild(t):
    f = np.zeros((3, 3), dtype=complex)
    i, j = np.triu_indices(3)
    f[i,j] = np.asarray(t["F_real_upper"]) + 1j*np.asarray(t["F_imag_upper"])
    f[j,i] = f[i,j]
    y = t["r"] * (np.diag(t["H_diag"]) - 3*complex(t["s_real"],t["s_imag"])*f)
    right = -(VEV*t["r"])**2*f/(1e-9*t["t_eV"])
    return y, right


def check_point(point):
    y, right = independently_rebuild(point["tensors"])
    _, descending, vh = np.linalg.svd(right)
    axes = vh.conj().T
    by_pole = ((abs(y @ axes)**2).sum(axis=0)*descending**2)/(8*pi**2)
    threshold = positive_threshold(y, right)
    np.testing.assert_allclose(sum(by_pole), threshold, rtol=4e-13)
    light = light_eV(y, right)
    light_masses = np.linalg.svd(light, compute_uv=False)
    floor = float(np.dot(descending[::-1]**3, light_masses*1e-9) / (8*pi**2*VEV**2))
    assert threshold >= floor*(1-1e-10)
    c = (HIGGS_POLE_GEV**2 / threshold)**(1/6)
    alternative = light_eV(c*y, c*c*right)
    np.testing.assert_allclose(alternative, light, rtol=2e-11, atol=2e-13)
    np.testing.assert_allclose(positive_threshold(c*y, c*c*right), HIGGS_POLE_GEV**2)
    mix_before = frob2(VEV*np.linalg.solve(right.T,y.T).T)
    mix_after = frob2(VEV*np.linalg.solve(c*c*right.T,c*y.T).T)
    np.testing.assert_allclose(mix_after*c*c, mix_before, rtol=2e-11,atol=0)
    # Masses fixed, reassign the light directions to the heavy poles in reverse order.
    # Check the entire mass matrix, including phases, under a complex angle.
    poles = descending[::-1]
    _, _, light_vh = np.linalg.svd(light)
    u = light_vh.conj().T
    u = u*np.exp(-0.5j*np.angle(np.diag(u.T@light@u)))
    def column_set(angle):
        o = np.eye(3,dtype=complex)
        a, b = np.cosh(angle), 1j*np.sinh(angle)
        o[:2,:2] = [[a,b],[-b,a]]
        return 1j*(u.conj()*np.sqrt(light_masses*1e-9))@o@np.diag(np.sqrt(poles))/VEV
    for x in [0.,0.25,1.]:
        yy = column_set(x)
        np.testing.assert_allclose(light_eV(yy,np.diag(poles)),light,rtol=1e-9,atol=1e-12)
        extra = (sum(poles[:2]**3)*sum(light_masses[:2])*1e-9*np.sinh(x)**2)/(8*pi**2*VEV**2)
        np.testing.assert_allclose(positive_threshold(yy,np.diag(poles)),floor+extra,rtol=2e-11)
    return {"common_tensor_positive_threshold_GeV2":threshold,"all_yukawas_floor_GeV2":floor,
            "actual_to_floor":threshold/floor,"largest_threshold_fraction":float(max(by_pole)/threshold),
            "coupling_rescale_for_125_GeV":c,"new_heavy_spectrum_GeV":(c*c*poles).tolist(),
            "weak_mixing_trace_before":mix_before,"weak_mixing_trace_after":mix_after}


def compare(a, b):
    if isinstance(a,dict):
        assert a.keys() == b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-8,atol=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    check = parser.parse_args().check
    raw = json.loads((HERE.parent/"1556"/"results.json").read_text(encoding="utf-8"))
    result = {key:check_point(v) for key,v in raw["candidates"].items()}
    author = json.loads((HERE/"results.json").read_text(encoding="utf-8"))
    conversion = (author["contract"]["u_GeV"]/VEV)**4
    for key, check_point_values in result.items():
        saved = author["candidates"][key]
        pair = [check_point_values["common_tensor_positive_threshold_GeV2"],
                check_point_values["all_yukawas_floor_GeV2"]]
        np.testing.assert_allclose(np.asarray(pair)*conversion,
                                   [saved["sum_GeV2"], saved["freely_reoriented_min_GeV2"]],
                                   rtol=2e-10, atol=0)
    path = HERE/"independent_results.json"
    if check:
        compare(result,json.loads(path.read_text(encoding="utf-8")))
        print("1558 six trace comparisons, six sharp minima and seesaw-preserving mass rescalings checked")
    else:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
