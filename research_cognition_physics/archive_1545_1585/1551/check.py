"""Finite exact type-I completion at fixed Dirac block and specified light polar unitary."""
import argparse
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def rotation(s12, s13, s23, phase):
    c12, c13, c23 = np.sqrt(1-np.square([s12,s13,s23]))
    z = np.exp(1j*phase)
    return np.array([
        [c12*c13, s12*c13, s13*z.conjugate()],
        [-s12*c23-c12*s23*s13*z, c12*c23-s12*s23*s13*z, s23*c13],
        [s12*s23-c12*c23*s13*z, -c12*s23-s12*c23*s13*z, c23*c13],
    ])


def norm(matrix):
    return np.linalg.norm(matrix,ord=2)


def inverse_light(d, u, light):
    eps = max(light)/np.linalg.svd(d,compute_uv=False)[-1]
    assert eps < 0.5
    h = np.eye(3)
    m = np.diag(light)
    for step in range(1,65):
        s = np.linalg.solve(d,(u@h).conj()@m)
        deficit = s.conj().T@s
        x, vec = np.linalg.eigh(deficit)
        # Evaluate sqrt(1-x)-1 without cancellation.
        delta = -x/(1+np.sqrt(1-x))
        next_h = np.eye(3)+(vec*delta)@vec.conj().T
        change = norm(h-next_h)
        h = next_h
        if change < 3e-15:
            break
    else:
        raise AssertionError("fixed point did not converge")
    n = u@h
    s = np.linalg.solve(d,n.conj()@m)
    rhs = s.conj()@m-d.T@n
    right = np.linalg.solve(s.T,rhs.T).T
    assert norm(right-right.T) < 2e-10*norm(right)
    # The exact proof gives symmetry. Remove only numeric skew before checking residuals.
    right = (right+right.T)/2
    columns = np.vstack((n,s))
    np.testing.assert_allclose(columns.conj().T@columns,np.eye(3),atol=5e-14)
    top_residual = norm(d@s-n.conj()@m)/max(light)
    bottom_residual = norm(d.T@n+right@s-s.conj()@m)/norm(d)
    assert top_residual < 1e-10 and bottom_residual < 3e-9
    return right,columns,step,float(max(top_residual,bottom_residual))


def evaluate():
    cfg = load(HERE/"inputs.json")
    high = load(HERE/cfg["charged_and_CKM"])
    down_leptons = load(HERE/cfg["down_and_charged"])["light_scales"]["1000"]
    osc = load(HERE/cfg["oscillation"])["NO"]
    scale = cfg["illustrative_vev_GeV"]/np.sqrt(2)
    rows = np.array(high["up_at_1TeV"]+down_leptons)*scale
    up, down, leptons = rows[:,0].reshape(3,3)
    v = rotation(*high["ckm_theta12_theta13_theta23_delta"])
    md = (v*down)@v.conj().T
    mu, me = np.diag(up), -np.diag(leptons)
    adj = (md-me)/4
    sing_up, sing_down = mu+adj, md-adj
    d = sing_up+3*adj
    np.testing.assert_allclose(mu,sing_up-adj)
    np.testing.assert_allclose(me,sing_down-3*adj,atol=1e-14)
    np.testing.assert_allclose(md,sing_down+adj)
    u = rotation(np.sqrt(osc["s12"][0]),np.sqrt(osc["s13"][0]),
                 np.sqrt(cfg["sine_squared_23"]),cfg["dirac_phase"])
    u = u@np.diag(np.exp(1j*np.array(cfg["majorana_column_phases"])))
    m0 = cfg["lightest_eV"]
    light = 1e-9*np.sqrt([m0*m0,m0*m0+osc["solar"][0],m0*m0+osc["atm_abs"][0]])
    _, columns, _, scaled_error = inverse_light(d,u,light)
    lower = float(sum(rows[i,0]-2*rows[i,1] for i in (0,3,6)))
    epsilon = cfg["comparison_light_upper_eV"]*1e-9/lower
    heavy_floor = lower*(np.sqrt(1-epsilon**2)-epsilon**2)/epsilon
    # Well-conditioned, finite-mixing full-spectrum experiment, including a nonsymmetric D.
    dtest = np.array([[1.2,0.2j,0.07],[-0.1,0.8+0.12j,0.09j],[0.15j,0.05,0.53]])
    mtest = np.array([.06,.11,.19])
    right, sub, iterations, stress_err = inverse_light(dtest,u,mtest)
    full = np.block([[np.zeros((3,3)),dtest],[dtest.T,right]])
    singular = np.sort(np.linalg.svd(full,compute_uv=False))
    np.testing.assert_allclose(singular[:3],mtest,atol=1e-13)
    ul, _, vh = np.linalg.svd(sub[:3])
    np.testing.assert_allclose(ul@vh,u,atol=2e-14)
    eps_test = max(mtest)/np.linalg.svd(dtest,compute_uv=False)[-1]
    assert np.linalg.norm(sub[3:],ord=2) > 0.05
    return {
        "two_width_uniform_D_min_GeV":lower,
        "any_CKM_any_polar_U_active_loss_upper":float(epsilon**2),
        "uniform_heavy_floor_for_positive_witness_GeV":float(heavy_floor),
        "benchmark_light_masses_eV":list(light*1e9),
        "actual_D_min_GeV":float(np.linalg.svd(d,compute_uv=False)[-1]),
        "benchmark_active_loss_from_sterile_rows":float(norm(columns[3:])**2),
        "benchmark_block_scaled_residual":float(scaled_error),
        "stress_full_singular_values":singular.tolist(),
        "stress_eps":float(eps_test),"stress_iterations":iterations,
        "stress_relative_residual":stress_err,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args = parser.parse_args()
    result = evaluate()
    path = HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        expected = load(path)
        assert result.keys()==expected.keys()
        for key in expected:
            if key in {"benchmark_block_scaled_residual", "stress_relative_residual"}:
                assert max(result[key], expected[key]) < 3e-9
            else:
                np.testing.assert_allclose(result[key],expected[key],atol=0,rtol=1e-9)
        print("1551: open charged-parameter construction, certified light subspace and physical stress spectrum passed")
