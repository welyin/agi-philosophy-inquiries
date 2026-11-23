"""Independent check with interval intersection, polynomial roots and complex congruences."""
import argparse
import itertools
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def get(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8"))


def spread_by_bisection(u, l, relative):
    # Intersect the three intervals for the real compensator.
    def feasible(radius):
        err = radius*l if relative else np.full(3, radius)
        return max((u-err)/l) <= min((u+err)/l)
    low, high = 0., max(u/l if relative else u)
    for _ in range(84):
        mid = (low+high)/2
        if feasible(mid):
            high = mid
        else:
            low = mid
    return high


def verify():
    book = get("results.json")
    pars = get("inputs.json")
    rng = np.random.default_rng(20261120)
    u = np.array([.07, .25, 1.3])
    l = np.array([.0001,.07,.27])
    hard = spread_by_bisection(u,l,False)
    eps = spread_by_bisection(u,l,True)
    # Diagonal aligned realization + exact formula from all six endpoint inequalities.
    pairs = itertools.combinations(range(len(l)),2)
    assert math.isclose(hard,max(abs(u[i]*l[j]-u[j]*l[i])/(l[i]+l[j]) for i,j in pairs))
    assert math.isclose(2*eps,max(u/l)-min(u/l))
    to_unitary = lambda z: np.linalg.qr(z)[0]
    for trial in range(80):
        q, v = [to_unitary(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))) for j in (0,1)]
        U, L = q @ np.diag(u) @ q.T, v @ np.diag(l) @ v.T
        c = (trial+1)*np.exp(1j*trial)
        invL = v.conj() @ np.diag(1/l) @ v.conj().T
        D = U-c*L
        assert np.linalg.norm(D,2) >= hard-1e-9
        assert np.linalg.norm(D@invL,2) >= eps-1e-9
        correction = D@invL@D.T
        lower = max(u[0]**2/l[0]-l[-1]*abs(c)**2,l[1]*abs(c)**2-u[-1]**2/l[1])-2*u[-1]*abs(c)
        assert np.linalg.norm(correction,2) >= lower-1e-8
        six = np.block([[L,D],[D.T,100*L]])
        assert np.linalg.norm(six,2) >= hard-1e-9

    up = get("../1386/inputs.json")["yukawa"]["up"]
    v, dv = pars["matched_vev_GeV"]
    u2 = (up[1][0]+2*up[1][1])*(v+2*dv)/math.sqrt(2)
    u3 = (up[2][0]-2*up[2][1])*(v-2*dv)/math.sqrt(2)
    rho = math.sqrt((7.388e-5-2*.078e-5)/(.002509+2*.000027))
    for row in book["spectral_bounds"]:
        t = row["spectral_transfer"]
        r = (rho-t)/(1+t)
        vals = [(u3*r-u2)/(1+r),(u3-u2/r)/(2*pars["mixing_expansion_budget"])]
        assert np.allclose(vals,[row["full_six_largest_mass_floor_GeV"],row["bare_heaviest_for_mixing_budget_GeV"]])
    win = book["unified_up_triplet_and_neutral_blocks_witness"]
    # Stable 2x2 light eigenvalue via determinant over heavy root, with arbitrary precision.
    with localcontext() as ctx:
        ctx.prec = 65
        dec = lambda value: Decimal(str(value))
        L = list(map(dec,win["triplet_diagonal_eV"]))
        uv = [dec(row[0])*dec(v)*Decimal(10)**9/Decimal(2).sqrt() for row in up]
        scale = dec(win["maximum_R_GeV"])*Decimal(10)**9/max(L)
        c = uv[-1]/L[-1]
        roots = []
        for left, top in zip(L,uv):
            R, d = scale*left, top-c*left
            heavy = (R+left+((R-left)**2+4*d**2).sqrt())/2
            roots.append(float((R*left-d**2)/heavy))
        assert np.allclose(roots,win["light_eV"],rtol=3e-13,atol=1e-18)
        assert np.allclose([roots[1]**2-roots[0]**2,roots[2]**2-roots[0]**2],
                           [7.388e-5,.002509],rtol=3e-13)
    return {"complex_six_block_trials":80,"absolute_minimax_test":round(hard,12),
            "ratio_minimax_test":round(eps,12),"stable_root_check_eV":[float(f"{r:.12g}") for r in roots],
            "oscillation_ratio_lower":round(rho,12)}


if __name__ == "__main__":
    args=argparse.ArgumentParser()
    cmd=args.add_mutually_exclusive_group(required=True)
    cmd.add_argument("--write-new",action="store_true")
    cmd.add_argument("--check",action="store_true")
    opt=args.parse_args()
    result=verify()
    if opt.write_new:
        with (ROOT/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result==get("independent_results.json")
    print("PASS independent minimax, six-block norms, near-singular correction, and high-precision light masses")
