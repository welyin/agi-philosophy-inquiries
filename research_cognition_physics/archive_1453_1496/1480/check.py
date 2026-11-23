"""Verify matrix anomalous dimensions, heavy thresholds, and light-mass holdout."""
import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def close(x, y):
    np.testing.assert_allclose(x, y, atol=2e-12, rtol=2e-12)


def kernel_frame(projection):
    """Return an orthonormal frame whose final vector is the specified light state."""
    v = np.r_[projection, np.sqrt(1 - np.dot(projection, projection))]
    basis = np.linalg.svd(v[None, :], full_matrices=True)[2][1:].T
    return np.column_stack((basis, v))


def family_test(v, data):
    g = data["g"]
    u, d = np.sqrt(8/5)*g, np.sqrt(6/5)*g
    unit = np.eye(3)
    up = [u * np.outer(e, e) for e in unit]
    down = [d * np.outer(v @ e, e) for e in unit]
    def gram(matrices):
        return np.array([[np.vdot(a,b) for b in matrices] for a in matrices])
    # The matter legs, both Higgs Gram blocks, and the pure adjoint equation.
    close(sum(3*a@a.conj().T + 2*b@b.conj().T for a,b in zip(up,down)), 36*g**2/5*unit)
    close(4*sum(b.conj().T@b for b in down), 24*g**2/5*unit)
    close(3*gram(up), 24*g**2/5*unit)
    close(4*gram(down), 24*g**2/5*unit)
    # Fourth Higgs has only the adjoint interaction, f44=g.
    close(24*abs(g)**2/5, 24*g**2/5)
    close(g**2 + 21/5*(15/7*g**2), 10*g**2)
    close(2*3+4+5-15, 0)
    h = np.array(data["up_projection"])
    k = np.array(data["down_projection"])
    yu = sum(z*a for z,a in zip(h,up))
    yd = sum(z*b for z,b in zip(k,down))
    close(yu, u*np.diag(h))
    close(yd, d*v@np.diag(k))
    lepton = yd.T
    close(np.linalg.svd(lepton, compute_uv=False), np.linalg.svd(yd, compute_uv=False))
    left, right = kernel_frame(h), kernel_frame(k)
    heavy = np.array(data["doublet_masses"])
    doublets = left @ np.diag(np.r_[heavy,0]) @ right.T
    triplets = doublets.copy()
    triplets[3,3] += data["adjoint_split"]
    effective = abs(np.linalg.det(triplets))/np.prod(heavy)
    expected = data["adjoint_split"]*np.sqrt((1-(np.linalg.norm(yu)/u)**2)*(1-(np.linalg.norm(yd)/d)**2))
    close(effective, expected)
    # Same masses and threshold, independently variable CP in the down orientation.
    j = np.imag(v[0,0]*v[1,1]*v[0,1].conj()*v[1,0].conj())
    return {"jarlskog":round(float(j),12), "heavy_mass_ratio":round(float(effective),12)}


def solve():
    data = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    q = data["lattice"]
    md, ed = q["md"]
    ms, es = q["ms"]
    mu, em = data["muonium"]["pole_ratio"]
    k = data["empirical_widths"]
    low = (md-k*ed)*(mu-k*em)/(ms+k*es)
    high = (md+k*ed)*(mu+k*em)/(ms-k*es)
    f = np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
    out = [family_test(v,data["example"]) for v in (np.eye(3),f,f.conj())]
    close(out[1]["jarlskog"], 1/(6*np.sqrt(3)))
    close(out[2]["jarlskog"], -out[1]["jarlskog"])
    assert 10 < low < md*mu/ms < high < 11
    return {"light_double_ratio":md*mu/ms,"two_width_box":[low,high],
            "minimal_four_factor_relative_radius":(np.sqrt(low)-1)/(np.sqrt(low)+1),
            "minimum_flavour_factor_spread":np.sqrt(low),"examples":out}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write-new", action="store_true")
    args = parser.parse_args()
    result = solve()
    dest = HERE/"results.json"
    if args.write_new:
        with dest.open("x", encoding="utf-8") as stream:
            json.dump(result, stream, indent=2)
            stream.write("\n")
    else:
        assert result == json.loads(dest.read_text(encoding="utf-8"))
        print("1480: matrix finiteness, arbitrary mixing witnesses, common threshold and mass-ratio test verified.")
