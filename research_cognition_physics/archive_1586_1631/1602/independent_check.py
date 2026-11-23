"""Bound WH times bb from the same two canonical doublet directions and the 4f width."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def make_budget(k, widths=2):
    masses = read(HERE.parent/"1386/inputs.json")["yukawa"]
    def norm_bounds(rows):
        return math.sqrt(sum((m+widths*e)**2 for m,e in rows))
    bmin = masses["down"][-1][0] - widths*masses["down"][-1][1]
    anorm = (3*norm_bounds(masses["down"])+norm_bounds(masses["charged_leptons"])/k)/(4*bmin)
    bnorm = (1+(masses["charged_leptons"][-1][0]+widths*masses["charged_leptons"][-1][1])/(k*bmin))/4
    return anorm, bnorm


def envelope(k, epsilon, cells=4096):
    anorm,bnorm = make_budget(k)
    ell = .98*k*k/(3*1.30)
    g = .195/.61
    largest = 0.
    for i in range(cells):
        lo, hi = i/cells, (i+1)/cells
        V = min(1., epsilon*(1-lo)+math.sqrt(hi))
        # x=1/d; if t may vanish, admit an arbitrarily large down amplitude.
        xmin = 0. if i==0 else 1/(epsilon*anorm+bnorm/math.sqrt(lo))
        p = 4*epsilon*anorm
        gterm = g*V**2/1.30
        xstar = 3*ell*p/(ell*p*p+gterm)
        x = max(xmin,xstar)
        den = 1+ell*max(0.,3-p*x)**2+gterm*x*x
        largest = max(largest,V*V/(.55*den))
    return largest + .02


def results():
    rate = read(HERE.parent/"1470/inputs.json")["associated_VHbb"]["WH"]
    minimum = rate[0]-2*rate[1]
    out = {}
    for k in (.7,.85,1.,1.15):
        cases = {str(r):envelope(k,r) for r in (0.,.02,.05,.1,.2,.4)}
        low,high = 0.,1.
        if cases["0.0"] < minimum:
            for _ in range(28):
                mid=(low+high)/2
                if envelope(k,mid,1024) < minimum:
                    low=mid
                else:
                    high=mid
        out[str(k)] = {"WH_upper_by_absolute_5_radial_response_budget":cases,
                       "abs_r_necessary_from_coarser_outer_envelope":low}
    assert out["1.0"]["WH_upper_by_absolute_5_radial_response_budget"]["0.0"] < minimum
    return {"two_width_WH_floor":minimum, "cases":out,
            "scope":"SM branching inputs bracketed by [.55,.61] for bb, 4f >=.195; k identifies matched lepton/down amplitudes. Positive extra total width. Two doublet vevs, optional singlet admixture."}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args=parser.parse_args()
    out=results()
    dest=HERE/"independent_results.json"
    if args.check:
        def close(a,b):
            if isinstance(a,dict):
                return a.keys()==b.keys() and all(close(v,b[k]) for k,v in a.items())
            if isinstance(a,str):
                return a==b
            return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)
        assert close(out, read(dest))
    else:
        with dest.open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    print(out["cases"]["1.0"])
