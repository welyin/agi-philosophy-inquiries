"""OPAL two-bin test of an arbitrary neutral transverse J=1 pole amplitude."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
HEL = (1, 0, -1)
PAIRS = [(a, b) for a in HEL for b in HEL if abs(a - b) <= 1]
CUT = math.cos(math.radians(CFG["cosine_cut_degrees"]))


def wigner(row, col, z):
    if row == col == 0:
        return z
    if row == 0:
        return col * np.sqrt(np.maximum(0., (1 - z*z) / 2))
    if col == 0:
        return -row * np.sqrt(np.maximum(0., (1 - z*z) / 2))
    return (1 + row*col*z) / 2


def quad_data(n):
    x, w = leggauss(n)
    k = np.concatenate([.5*(1-CUT)*x + s*.5*(1+CUT) for s in (-1, 1)])
    phi = np.arange(4*n) * (2*math.pi/(4*n))
    scale = np.sqrt(1-k*k)[:, None]
    return (scale*np.cos(phi), scale*np.sin(phi), k[:, None],
            np.tile(w, 2)[:, None] * ((1-CUT)/2) * (2*math.pi/(4*n)))


def accepted_spin(z, beta, grid):
    # Integrate two lab veto caps, using the inverse boost into the rest frame.
    lx, ly, lz, weights = grid
    sine = math.sqrt(1-z*z)
    longitudinal = sine*lx + z*lz
    gamma = 1/math.sqrt(1-beta**2)
    den = 1-beta*longitudinal
    t = np.clip((longitudinal-beta)/den, -1., 1.)
    phase = np.exp(1j*np.arctan2(ly, z*lx-sine*lz))
    weighted = (3/(4*math.pi) * weights/(gamma*den)**2).ravel()
    matrices = []
    for helicity in (-1, 1):
        a = np.array([(phase**h*wigner(h, helicity, t)).ravel() for h in HEL])
        matrices.append(np.eye(3) - ((a*weighted)@a.conj().T).real)
    return matrices


def rate_matrices(bounds, beta, order):
    grid = quad_data(order)
    nodes, weights = leggauss(order)
    a, b = bounds
    out = np.zeros((2,7,7))
    for xi, wi in zip(nodes, weights):
        z = (b+a)/2 + (b-a)*xi/2
        minus, plus = accepted_spin(z, beta, grid)
        decay = np.zeros((7,7))
        for i,(l,r) in enumerate(PAIRS):
            for j,(k,t) in enumerate(PAIRS):
                if r == t:
                    decay[i,j] += minus[1-l,1-k]/2
                if l == k:
                    decay[i,j] += plus[1-r,1-t]/2
        for index,m in enumerate((-1,1)):
            d = np.array([wigner(m,l-r,z) for l,r in PAIRS])
            out[index] += (wi*(b-a)/2) * decay * np.outer(d,d)
    return out


def fixed_energy(energy, order):
    beta = math.sqrt(1 - (2*CFG["pole_mass_GeV_for_shape_examples"]/energy)**2)
    middle, forward = [rate_matrices(CFG[key], beta, order)
                       for key in ("z_middle","z_forward")]
    maxima = []
    for m, f in zip(middle, forward):
        chol = np.linalg.cholesky(m)
        op = np.linalg.solve(chol, np.linalg.solve(chol,f).T).T
        maxima.append(float(np.linalg.eigvalsh(op)[-1]))
    return max(maxima)


def inclusive_max():
    middle, front = (CFG[key] for key in ("z_middle","z_forward"))
    def integ(bounds):
        a,b=bounds
        return np.array([((1+b)**3-(1+a)**3)/3,
                         ((1-a)**3-(1-b)**3)/3, b-a-(b**3-a**3)/3])
    return float(max(integ(front)/integ(middle)))


def analytic_bound():
    # Interval upper squares: extrema of every elementary factor on its cell.
    beta = np.linspace(*CFG["common_beta_range_for_certificate"], 513)
    z = np.linspace(*CFG["z_middle"], 257)
    low = beta[:-1,None]*z[None,:-1]
    high = beta[1:,None]*z[None,1:]
    d2 = (1-beta[:-1,None]**2)*(1-CUT*CUT)
    upper = np.zeros_like(low)
    for v in (CUT-high, CUT+low):
        cap = (1-v/np.sqrt(v*v+d2))/2
        # Spin-1 coherent-state effect on the cap. Highest eigenvalue.
        upper += 1-(1-cap)**3
    acceptance = 1-float(np.max(upper))-1e-10
    return acceptance, inclusive_max()/acceptance+1e-10


def compute():
    acceptance, cap_bound = analytic_bound()
    out={"inclusive_sharp":inclusive_max(), "uniform_middle_acceptance_lower":acceptance,
         "uniform_fiducial_ratio_upper":cap_bound,"energies":{}}
    for label,(m,f) in CFG["OPAL_middle_forward"].items():
        energy=float(label)
        sharp=fixed_energy(energy,12)
        refined=fixed_energy(energy,18)
        assert abs(refined-sharp)<3e-7
        row={"fiducial_fixed_energy_sharp":refined,"observed_central":f[0]/m[0],
             "rectangles":{}}
        for w in (1,2):
            lo=(f[0]-w*sum(f[1:]))/(m[0]+w*sum(m[1:]))
            row["rectangles"][str(w)]={"ratio_lower":lo,
               "uniform_extra_transfer_limit":max(0.,(lo-cap_bound)/(lo+cap_bound)),
               "fixed_energy_extra_transfer_limit":max(0.,(lo-refined)/(lo+refined))}
        out["energies"][label]=row
    assert out["energies"]["198.4"]["rectangles"]["2"]["uniform_extra_transfer_limit"]>0
    return out


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for key in a:
            compare(a[key],b[key])
    else:
        assert math.isclose(a,b,rel_tol=5e-9,abs_tol=5e-10),(a,b)


if __name__ == "__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write",action="store_true")
    args=cli.parse_args()
    result=compute()
    path=HERE/"results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2)
            handle.write("\n")
    else:
        compare(result,json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(result,ensure_ascii=False))
