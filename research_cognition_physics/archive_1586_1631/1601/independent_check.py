"""A protected massless neutrino read through the shared right-handed plane."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / "independent_results.json"


def previous(round_id):
    return json.loads((HERE.parent / str(round_id) / "inputs.json").read_text(encoding="utf8"))


def frame(line):
    return np.linalg.svd(line.conj()[None, :], full_matrices=True)[2].conj().T[:, 1:]


def eigframe(y, line):
    q = frame(line)
    values, rotation = np.linalg.eigh((y[:, None]*q).conj().T @ (y[:, None]*q))
    return values, q @ rotation


def plane_shape(y, weights):
    h = (y / y[-1])**2
    tr = np.dot(h, 1-weights)
    det = np.dot(weights, [h[1], h[0], h[0]*h[1]])
    return det / tr**2


def unitary_ckm(angles):
    s12, s13, s23 = np.sin(angles[:3])
    c12, c13, c23 = np.cos(angles[:3])
    phase = np.exp(1j * angles[3])
    return np.array([[c12*c13, s12*c13, s13/phase],
                     [-s12*c23-c12*s23*s13*phase, c12*c23-s12*s23*s13*phase, s23*c13],
                     [s12*s23-c12*c23*s13*phase, -c12*s23-s12*c23*s13*phase, c23*c13]])


def witness(down, charged, vrow, reactor, atmospheric):
    tag = vrow.conj()/down
    tag = tag/np.linalg.norm(tag)
    target = plane_shape(down, abs(tag)**2)

    def lepton_column(theta):
        return np.sqrt([reactor, (1-reactor)*theta, (1-reactor)*(1-theta)])

    lo, hi = atmospheric
    assert plane_shape(charged, lepton_column(lo)**2) > target > plane_shape(charged, lepton_column(hi)**2)
    for _ in range(70):
        mid = (lo + hi)/2
        if plane_shape(charged, lepton_column(mid)**2) > target:
            lo = mid
        else:
            hi = mid
    neutral = lepton_column((lo + hi)/2)
    sd, qd = eigframe(down, tag)
    sl, ql = eigframe(charged, neutral)
    k = math.sqrt(sl.sum()/sd.sum())
    # First match the common family sphere on the input; then the two
    # output columns by a left unitary. The correction is free on the tag.
    rotation = ql @ qd.conj().T + np.outer(neutral, tag.conj())
    initial = np.diag(charged) @ rotation
    cd = down[:, None]*qd / np.sqrt(sd)
    cl = charged[:, None]*ql / np.sqrt(sl)
    rd = np.linalg.svd(cd.conj().T, full_matrices=True)[2].conj().T[:, 2]
    rl = np.linalg.svd(cl.conj().T, full_matrices=True)[2].conj().T[:, 2]
    y = (cd @ cl.conj().T + np.outer(rd, rl.conj())) @ initial
    correction = y-k*np.diag(down)
    a = correction @ tag
    u, s, vh = np.linalg.svd(y)
    np.testing.assert_allclose(s[::-1], charged, rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(correction, np.outer(a, tag.conj()), rtol=2e-11, atol=2e-15)
    np.testing.assert_allclose(abs(vh[::-1] @ tag)**2, neutral**2, atol=3e-12)
    lepton_mass2 = charged**2
    weights = [np.prod(x-k*k*sd)/np.prod(x-np.delete(lepton_mass2, i))
               for i, x in enumerate(lepton_mass2)]
    np.testing.assert_allclose(weights, neutral**2, rtol=3e-11, atol=1e-13)
    return {"sin2_theta23": float((lo+hi)/2), "sin2_theta13": reactor,
            "charged_singular_values": charged.tolist(), "down": down.tolist(),
            "k": k, "rank1_correction_singular_value": float(np.linalg.norm(a)),
            "p_e_mu_tau": (neutral**2).tolist(), "scale_free_plane_shape": float(target)}


def evaluate():
    charged_rows = np.array(previous(1402)["scales"]["1e16"])
    angles = np.array(previous(1397)["high_scale"]["angles_theta12_theta13_theta23_delta"])
    centers, errs = charged_rows.T
    low = centers-2*errs
    high = centers+2*errs
    angle_low, angle_high = angles[:, 0]-2*angles[:, 1], angles[:, 0]+2*angles[:, 1]
    l12, l13, l23 = np.sin(angle_low[:3])
    u12, u13, u23 = np.sin(angle_high[:3])
    c12, c13, c23 = np.cos(angle_high[:3])
    rows = [[(c12*c13, 1.), (l12*c13, u12), (l13, u13)],
            [(l12*c23-u23*u13, u12+u23*u13), (c12*c23-u12*u23*u13, 1.), (l23*c13, u23)],
            [(l12*l23-u13, u12*u23+u13), (c12*l23-u12*u13, u23+u12*u13), (c23*c13, 1.)]]
    osc = previous(1370)
    atm = previous(1580)["theta23_sin2"]["IO"][1:]
    pe_no_low = (1-osc["NO"]["s12"][2]) * (1-osc["NO"]["s13"][2])
    re_hi, rm_hi = (high[3]/low[5])**2, (high[4]/low[5])**2
    rm_low = (low[4]/high[5])**2
    quark_trace_max = 1 + (high[1]/low[2])**2 + (high[0]/low[2])**2
    v = unitary_ckm(angles[:, 0])
    result = {}
    for j, row in enumerate(rows):
        vmin = np.array(row)[:, 0]
        vmax = np.array(row)[:, 1]
        def prob_min(i):
            return 1/(1+sum((vmax[k]*high[i]/(vmin[i]*low[k]))**2 for k in range(3) if k != i))
        def prob_max(i):
            return 1/(1+sum((vmin[k]*low[i]/(vmax[i]*high[k]))**2 for k in range(3) if k != i))
        pd_min, pb_max = prob_min(0), prob_max(2)
        qhi = (high[1]/low[2])**2/(1-pb_max)
        qlo = pd_min*(low[1]/high[2])**2 / quark_trace_max**2
        pe_cap = qhi*(1+re_hi+rm_hi)**2/rm_low
        gap = math.acos(math.sqrt(pe_cap)) - math.acos(math.sqrt(pe_no_low))
        io_pe_lo, io_pe_hi = osc["IO"]["s13"][1:]
        trace_cap = math.sqrt((io_pe_hi*rm_hi + re_hi*(1+rm_hi))/qlo)
        theta23_cap = (trace_cap-io_pe_lo)/(1-io_pe_lo)
        # A shared charged-spectrum point on the declared marginal edges.
        down = np.array([centers[0], low[1], high[2]])
        charged = np.array([centers[3], high[4], low[5]])
        fit = witness(down, charged, v[j], io_pe_hi, atm)
        result["uct"[j]] = {"quark_shape_bounds": [float(qlo),float(qhi)],
            "protected_electron_weight_upper": float(pe_cap),
            "normal_angle_min_deg": math.degrees(gap), "IO_sin2theta23_upper": float(theta23_cap),
            "IO_joint_witness": fit}
        assert gap > math.radians(27) and theta23_cap < 0.584
    return {"NO_observed_electron_weight_lower":pe_no_low, "rows": result}


def check(a, b):
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for k in a:
            check(a[k], b[k])
    else:
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=1e-13)


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    g = cli.add_mutually_exclusive_group(required=True)
    g.add_argument("--write", action="store_true")
    g.add_argument("--check", action="store_true")
    args = cli.parse_args()
    data = evaluate()
    if args.write:
        with RESULT.open("x", encoding="utf8") as output:
            json.dump(data, output, indent=2)
            output.write("\n")
    else:
        check(data, json.loads(RESULT.read_text(encoding="utf8")))
        print("1601 independent plane compression, NO limit and all three IO completions passed")
