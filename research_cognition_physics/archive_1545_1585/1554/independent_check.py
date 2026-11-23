"""Independently audit the flavor rank at published hierarchical matrices."""
import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# Dueck and Rodejohann, arXiv:1306.4468v3, eqs. (15) and (16).
# Order of F is 11,12,13,22,23,33. Each printed decimal is held exact.
POINTS = {
    "MN_no_RGE_eq15": [
        "1.22387e-6 5.92428e-5 6.29473e-3",
        "-2.95102e-6 1.27484e-5 1.07772e-4 -1.538e-4 -2.67281e-4 -7.38503e-4",
        "-3.48291e-6 -7.53714e-8 6.02931e-5 6.75236e-5 2.48978e-4 -1.44559e-3",
        "68.9624 .370726 .063044 3.014e-16"],
    "MN_RGE_eq16_UV_only": [
        "1.15249e-6 6.71983e-5 6.70159e-3",
        "-2.25817e-6 1.22057e-5 -1.49653e-4 -2.06635e-4 3.76355e-4 -7.01333e-4",
        "7.40559e-6 -1.39062e-5 8.30809e-5 -1.34686e-5 2.15049e-4 -7.53673e-4",
        "63.9043 .409807 -.0420522 3.39236e-16"],
}


def fractions(lines):
    x = [Fraction(a) for line in lines for a in line.split()]
    # An independently chosen mass unit, numerically close to the historical v=174 GeV.
    x[-1] *= x[15]**2 * 174**2 * 10**9
    return x


def form(x, iunit):
    h = np.diag(x[:3])
    trip = np.triu_indices(3)
    f = np.zeros((3, 3), dtype=np.result_type(x, iunit))
    f[trip] = x[3:9] + iunit*x[9:15]
    f[trip[1], trip[0]] = f[trip]
    fs = f * (x[16] + iunit*x[17])
    return (x[15]*(h+fs), h+f, h-3*f, h-3*fs, f)


def adj(m, p):
    return np.array([np.cross(m[(i+1)%3], m[(i+2)%3]) for i in range(3)]).T % p


def numerators(x, p, j):
    groups, dets = [], []
    for sign in (1, -1):
        u, d, e, w, f = [m % p for m in form(x, sign*j)]
        a = adj(f, p)
        dets.append(int(f[0] @ a[:, 0]) % p)
        kn = ((w @ a % p) @ w.T % p) * x[18] % p
        groups.append((u, d, e, kn))
    den = dets[0]*dets[1] % p
    # Convert *all* observables to polynomial numerators over den^3.
    den3 = pow(den, 3, p)
    hu, hd, he, hn = [(minus.T @ plus)%p for plus, minus in zip(*groups)]
    b = (hn - np.eye(3, dtype=int)*(int(hn.trace())*pow(3, -1, p)%p))%p
    def tr(*ms):
        t = ms[0]
        for m in ms[1:]:
            t = t @ m % p
        return int(t.trace())%p
    charged = []
    for h in (hu, hd, he):
        charged.extend([tr(h), tr(h,h), tr(h,h,h)])
    charged.extend([tr(hu,hd),tr(hu,hd,hd),tr(hu,hu,hd),tr(hu,hu,hd,hd)])
    values = [q*den3 for q in charged]
    values += [tr(b,b)*den,tr(b,b,b)]
    values += [tr(he,b)*den**2,tr(he,b,b)*den,tr(he,he,b)*den**2,tr(he,he,b,b)*den]
    return np.array([v%p for v in values]+[den3], dtype=np.int64)


def interpolate_jacobian(x0, p):
    x = np.array([(v.numerator*pow(v.denominator,-1,p))%p for v in x0],dtype=np.int64)
    j = next(k for k in range(2,p) if k*k%p == p-1)
    n = 42
    # Exact derivative at zero in the Lagrange basis of nodes 0,...,n.
    weights = [(-sum(pow(k,-1,p) for k in range(1,n+1)))%p]
    weights += [((-1)**(k+1)*math.comb(n,k)*pow(k,-1,p))%p for k in range(1,n+1)]
    base = numerators(x,p,j)
    inv = pow(int(base[-1]),-1,p)
    columns = []
    for e in np.eye(19,dtype=np.int64):
        diff = sum(w*numerators((x+k*e)%p,p,j) for k,w in enumerate(weights))%p
        columns.append(((diff[:-1]*base[-1]-base[:-1]*diff[-1])*inv**2)%p)
    rows = np.array(columns).T.copy()
    determinant = 1
    for i in range(19):
        pivot = next((j for j in range(i,19) if rows[j,i]),None)
        if pivot is None:
            return 0
        if pivot != i:
            rows[[i,pivot]] = rows[[pivot,i]]
            determinant *= -1
        elt = int(rows[i,i])
        determinant = determinant*elt%p
        for k in range(i+1,19):
            t = int(rows[k,i])*pow(elt,-1,p)%p
            rows[k]=(rows[k]-t*rows[i])%p
    return int(determinant)


def physical(x):
    u,d,e,w,f = form(np.array([float(a) for a in x]), 1j)
    n = float(x[-1])*w@np.linalg.solve(f,w.T)
    evec = np.linalg.eigh(e@e.conj().T)[1][:,0]
    m2,axes = np.linalg.eigh(n@n.conj().T)
    probs = np.abs(evec.conj()@axes)**2
    return {"lightest_meV":float(1000*np.sqrt(m2[0])),
            "m_beta_meV":float(1000*np.sqrt(probs@m2)),
            "m_bb_meV":float(1000*abs(evec.conj()@n@evec.conj())),
            "solar_to_span":float((m2[1]-m2[0])/(m2[2]-m2[0])),
            "s13_squared":float(probs[2])}


def evaluate():
    data = {}
    for key,table in POINTS.items():
        x = fractions(table)
        dets = {str(p):interpolate_jacobian(x,p) for p in (1009,1013)}
        assert all(dets.values())
        data[key] = {"jacobian_determinant_mod_p":dets,
                     "single_matching_scale_tree_observables":physical(x)}
    return data


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    check=parser.parse_args().check
    data=evaluate()
    path=HERE/"independent_results.json"
    if check:
        old=json.loads(path.read_text(encoding="utf-8"))
        for key, val in data.items():
            assert val["jacobian_determinant_mod_p"]==old[key]["jacobian_determinant_mod_p"]
            for name, number in val["single_matching_scale_tree_observables"].items():
                assert np.isclose(number,old[key]["single_matching_scale_tree_observables"][name],rtol=2e-8,atol=1e-12)
        print("1554: exact interpolated rank at two published hierarchical inputs")
    else:
        with path.open("x", encoding="utf-8") as out:
            json.dump(data,out,indent=2)
            out.write("\n")
