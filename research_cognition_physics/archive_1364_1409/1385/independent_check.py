"""Independent projected weak-charge check for aligned singlets, at tree level."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def mix_with_singlet(direction, probability):
    rotation = np.eye(4)
    cosine = sqrt(1 - probability)
    sine = sqrt(probability)
    rotation[:3, :3] += (cosine - 1) * np.outer(direction, direction)
    rotation[:3, 3] = -sine * direction
    rotation[3, :3] = sine * direction
    rotation[3, 3] = cosine
    return rotation


def calculate():
    old = json.loads((HERE.parent / "1381" / "results.json").read_text(encoding="utf-8"))
    complete = (np.array(old["explicit_four_unitary_real"])
                + 1j * np.array(old["explicit_four_unitary_imag"]))
    assert np.allclose(complete @ complete.conj().T, np.eye(4), atol=1e-13, rtol=0)
    # Archive orders the columns d,b,s,x. Choose d,s,b for both mass bases.
    row = np.abs(complete[0, [0, 2, 1]])
    loss = 1 - row @ row
    assert abs(loss - abs(complete[0, 3])**2) < 2e-15
    v = row / np.linalg.norm(row)
    r = v[0]**2
    e = np.array([1., 0., 0.])
    n = (e - v) / np.linalg.norm(e - v)
    charged = np.zeros((4, 4))
    charged[:3, :3] = np.eye(3) - 2 * np.outer(n, n)
    p = np.diag([1., 1., 1., 0.])
    null_x = 2*r*loss / (2+r+sqrt((2+r)**2-8*r*loss))
    samples = {}
    for label, x in (("up_singlet", loss), ("proton_null", null_x), ("down_singlet", 0.)):
        y = (loss-x) / (1-x)
        u = mix_with_singlet(e, x)
        d = mix_with_singlet(v, y)
        current = u.T @ charged @ d
        up = u.T @ p @ u
        down = d.T @ p @ d
        assert np.allclose(current @ current.T, up, atol=1e-13, rtol=0)
        assert np.allclose(current.T @ current, down, atol=1e-13, rtol=0)
        assert np.allclose(current[0, :3], row, atol=1e-13, rtol=0)
        shift_u, shift_d = up[0, 0]-1, 1-down[0, 0]
        proton, neutron = 2*shift_u+shift_d, shift_u+2*shift_d
        assert abs(proton - (r*y - 2*x)) < 1e-13
        assert abs(abs(down[0, 1]) - y*abs(v[0]*v[1])) < 1e-13
        samples[label] = dict(X=float(x), Y=float(y), proton=float(proton), neutron=float(neutron),
                              Cs133=float(55*proton + 78*neutron),
                              ds_projector=float(abs(down[0, 1])))
    assert abs(samples["proton_null"]["proton"]) < 1e-13
    assert abs(samples["proton_null"]["neutron"] - 3*null_x) < 1e-13
    assert samples["proton_null"]["ds_projector"] > 0
    assert samples["up_singlet"]["proton"] < 0 < samples["down_singlet"]["proton"]
    # A different alignment in the author's model couples only d to the new singlet.
    # At minimax the two opposite target residuals have equal weighted magnitude.
    width_p = np.linalg.norm([.0045, .0003])
    width_cs = np.linalg.norm([.29, .32, .02])
    k = (2*width_cs+188*width_p)/(width_cs+211*width_p)
    t = row[0]**2
    scale = 1 + k*(t+loss)
    x_best = 2*loss/(scale + sqrt(scale**2-4*k*loss))
    y_best = k*x_best
    assert abs(loss - x_best - t*y_best/(1-y_best)) < 1e-13
    score = (2*x_best-y_best)/width_p
    assert abs(score - (211*y_best-188*x_best)/width_cs) < 1e-13
    # Direct scan on the full admissible curve checks the minimax branch.
    grid = np.linspace(0., loss, 10001)
    y_grid = (loss-grid)/(t+loss-grid)
    scan = np.min(np.maximum(abs(y_grid-2*grid)/width_p, abs(211*y_grid-188*grid)/width_cs))
    assert score-1e-12 < scan < score+.001
    return dict(row_dsb=row.tolist(), loss=loss, branches=samples, d_aligned_minimax=score,
                full_doublet_tree_shifts=dict(proton=0, neutron=0, ds_projector=0))


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--check", action="store_true")
    options = cli.parse_args()
    findings = calculate()
    target = HERE / "independent_results.json"
    if options.check:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert findings.keys() == saved.keys()
        np.testing.assert_allclose(findings["row_dsb"], saved["row_dsb"], rtol=0, atol=1e-13)
        np.testing.assert_allclose(findings["loss"], saved["loss"], rtol=0, atol=1e-13)
        np.testing.assert_allclose(findings["d_aligned_minimax"], saved["d_aligned_minimax"],
                                   rtol=0, atol=1e-13)
        assert findings["full_doublet_tree_shifts"] == saved["full_doublet_tree_shifts"]
        assert findings["branches"].keys() == saved["branches"].keys()
        for label, branch in findings["branches"].items():
            assert branch.keys() == saved["branches"][label].keys()
            for key, value in branch.items():
                assert abs(value - saved["branches"][label][key]) < 1e-13
        print("1385 independent neutral and flavor-changing current checks passed.")
    else:
        print(json.dumps(findings, ensure_ascii=False, indent=2))
