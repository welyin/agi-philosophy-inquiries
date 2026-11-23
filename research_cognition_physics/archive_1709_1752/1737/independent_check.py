"""Profile the original lattice covariance for equal-charge mass ratios."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
# Fermilab/MILC/TUMQCD 1802.04248, Table IV, charm-mass coordinate marginalized.
# Entries are b/c, s/c, d/c, u/c in the authors' four-sea-flavor pure QCD.
CENTER = np.array([4.578, .08487, .004291, .001955])
ERROR = np.array([.008, .00018, .000039, .000037])
CORR = np.array([
    [1., .45502225, .04855992, .23627864],
    [.45502225, 1., .43609054, .47252309],
    [.04855992, .43609054, 1., -.32724921],
    [.23627864, .47252309, -.32724921, 1.],
])
COV = np.outer(ERROR, ERROR)*CORR


def read(n, filename="results.json"):
    return json.loads((ARCHIVE/str(n)/filename).read_text(encoding="utf-8"))


def linear_row(indices, value):
    numerator, denominator = indices
    row = np.eye(4)[numerator].copy()
    rhs = value
    if denominator is not None:
        row[denominator] -= value
        rhs = 0.
    return row, rhs


def fit(probes, relative=0.):
    """Find the closest covariance point to the exact ratio slab intersection."""
    constraints = []
    for indices, high in probes:
        bounds = high/(1+relative), high/(1-relative)
        bounds_ineq = []
        for sign,limit in zip((-1,1),bounds):
            row,rhs = linear_row(indices,limit)
            norm = np.linalg.norm(row*ERROR)
            # A z <= b, z the displacement in each coordinate's published width.
            bounds_ineq.append((sign*row*ERROR/norm, sign*(rhs-row@CENTER)/norm))
        constraints.append(bounds_ineq)
    best = None
    for faces in itertools.product((-1,0,1), repeat=len(probes)):
        indices = [constraints[i][side == 1] for i,side in enumerate(faces) if side]
        disp = np.zeros(4)
        distance = 0.
        if indices:
            rows, rhs = zip(*indices)
            rows, rhs = np.array(rows), np.array(rhs)
            mult = np.linalg.solve(rows @ CORR @ rows.T, rhs)
            disp = CORR @ rows.T @ mult
            distance = math.sqrt(max(0.,rhs @ mult))
        if any(row@disp > rhs+2e-9 for pair in constraints for row,rhs in pair):
            continue
        if best is None or distance < best[0]:
            best = distance, CENTER+ERROR*disp
    assert best is not None and min(best[1]) > 0
    actual = [(best[1][i]/(1. if j is None else best[1][j])) for (i,j),_ in probes]
    shifts = [high/lo-1. for (_,high),lo in zip(probes,actual)]
    assert max(abs(x) for x in shifts) < relative+1.e-9
    return best[0], best[1], shifts


def budget(probes, radius):
    if fit(probes)[0] <= radius:
        return 0.
    lo,hi = 0.,.1
    assert fit(probes,hi)[0] <= radius
    for _ in range(35):
        mid=(lo+hi)/2
        if fit(probes,mid)[0] > radius:
            lo=mid
        else:
            hi=mid
    return hi


def number(value):
    return 0. if abs(value)<5e-13 else float(f"{value:.9g}")


def evidence(up, down, charm):
    d,s,b=down
    probes = [((3,None),up/charm),((2,1),d/s)]
    R,x,_ = fit(probes)
    both = np.array([linear_row(ids,val)[0] for ids,val in probes])
    rhs = np.array([linear_row(ids,val)[1] for ids,val in probes])
    res = rhs-both@CENTER
    assert abs(R**2 - res @ np.linalg.solve(both@COV@both.T,res)) < 1e-10
    A=both[0]
    B=both[1]
    rho = A@COV@B/math.sqrt((A@COV@A)*(B@COV@B))
    record = {}
    for name,obs in [("ds_only",probes[1:]),("uc_ds",probes),("uc_ds_bs",probes+[((0,1),b/s)])]:
        record[name] = {
            "radius_by_relative_matching":{str(e):number(fit(obs,e)[0]) for e in (0.,.003,.005,.01,.02)},
            "min_relative_matching_for_radius": {str(rad):number(budget(obs,rad)) for rad in (1.,2.,math.sqrt(5.991464547),3.)}}
    record["equal_charge_point"] = {"uc":number(up/charm),"ds":number(d/s),"bs":number(b/s),
                                    "light_test_corr_of_residuals":number(rho),
                                    "light_test_best_published_vector":list(map(number,x)),
                                    "bs_at_light_only_best":number(x[0]/x[1])}
    # HPQCD 2102.09609 directly includes quenched QED at 3 GeV for b/c.
    # QCD+QED EFT: d ln(mb/mc)/d ln(mu)=alpha/(2 pi) at leading QED.
    # Unknown mixed, matching, sea-QED or full-SM effects go into the residual.
    log_range=math.log(91.1876/3)/(2*math.pi)
    run_lo,run_hi=math.exp(log_range/137.036), math.exp(log_range/127.)
    record["HPQCD_bc_3GeV_to_MZ"]={
        "experimental_lattice_ratio":[4.586,.012], "predicted_bc_MZ":number(b/charm),
        "one_loop_QED_transport_range":[number(run_lo),number(run_hi)],
        "relative_finite_factor_minus_one_at_one_width":[
            number((b/charm)/(4.598*run_hi)-1), number((b/charm)/(4.574*run_lo)-1)]}
    return record


def compute():
    previous = read(1736)["witnesses"]
    charm=read(1386,"inputs.json")["yukawa"]["up"][1][0]
    results = {key:evidence(val["up_MZ"],val["down_MZ"],charm) for key,val in previous.items()}
    # Published pure-QCD heavy/strange ratio sensitivity for the old length solutions.
    old=read(1732)["solutions"]
    alternate={}
    for tag,data in old.items():
        by_study={}
        for label,(central,sigma) in {"MILC_2018":(53.94,.12),"ETMC_2016":(51.4,1.4)}.items():
            by_study[label]={}
            for e in (0.,.01,.02):
                for n in (1,2):
                    admissible=[]
                    for pt in data:
                        d,s,b=pt["down_strange_bottom_MZ"]
                        if b/s/(1+e) <= central+n*sigma and b/s/(1-e) >= central-n*sigma:
                            admissible.append(pt["k"])
                    by_study[label][f"e={e},n={n}"]=admissible
        alternate[tag]=by_study
    # Same common-scale zeta_m as MILC Eq. 5.14 and light-quark Eq. 4 at mu=mb(5).
    # This illustrative alpha fixes the size of an individual conversion, not the ratio.
    a=.225/math.pi
    zeta=1+.2060*a**2+(1.8476+.0247*4)*a**3+(6.850-1.466*4+.05616*16)*a**4
    assert abs((4.578/zeta)/(.08487/zeta)/(4.578/.08487)-1.)<1e-15
    return {"ordering":["b/c","s/c","d/c","u/c"],"lattice_center":CENTER.tolist(),
            "lattice_marginal_errors":ERROR.tolist(),"published_corr":CORR.tolist(),
            "common_mu_four_loop_zeta_m_alpha_s_.225":number(zeta),
            "same_charge_K_reference":1., "fixed_yc_MZ":charm,
            "five_witnesses":results, "alternative_b_s_length_subsets":alternate}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    file=HERE/"independent_results.json"
    results=compute()
    if args.check:
        if results!=json.loads(file.read_text(encoding="utf-8")):
            raise SystemExit("1737 independent covariance check failed")
        print("1737 independent: five exact ratio profiles, bottom threshold convention and two lattice length tests verified")
    else:
        file.write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
