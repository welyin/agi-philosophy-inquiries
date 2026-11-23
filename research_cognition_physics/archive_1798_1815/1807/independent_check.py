"""Match the same vectorlike fermions in an equal positive gauge trace and in threshold running."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PI = math.pi
BASE_B = np.array([41/10, -19/6, -7.])
DOUBLET = np.array([3/5, 1., 0.])
TRIPLET = np.array([2/5, 0., 1.])
ORTHOGONAL = np.cross(np.ones(3), DOUBLET)


def read(path):
    return json.loads((HERE/path).read_bytes())


def coefficients(low):
    time = 2*PI*(ORTHOGONAL @ low)/(ORTHOGONAL @ BASE_B)
    at = low - BASE_B*time/(2*PI)
    x0, k = at[2], at[1]-at[2]
    np.testing.assert_allclose(at,x0+k*DOUBLET,atol=4e-14)
    return time, x0, k


def trial(low):
    t,x0,k=coefficients(low)
    max_log_weight=t/(3*PI)
    # Positive margins exclude a lone doublet, two doublets, and one completed (D,T) pair.
    direction_log = 3*PI*(x0-7*k)/8
    return {"log_scale":t,"X":x0,"K":k,"C_minus_K_for_one_D":x0/6-k,
            "single_complete_pair_missing_K":k-max_log_weight,
            "one_D_two_T_missing_K":k-(1.25*max_log_weight-x0/8),
            "two_D_one_T_log_min_geometric_mean_D_over_T":direction_log,
            "two_complete_pairs_log_mean_T_over_D":3*PI*k/2}


def example(low,reference,ld,lt):
    t,x0,k=coefficients(low)
    n_d,n_t=len(ld),len(lt)
    c=(x0-sum(lt)/(3*PI))/(6+n_t)
    assert c>0 and all(0<=v<=t for v in [*ld,*lt])
    trace=c*(6+n_d*DOUBLET+n_t*TRIPLET)
    threshold=(DOUBLET*sum(ld)+TRIPLET*sum(lt))/(3*PI)
    reconstructed=trace+threshold+BASE_B*t/(2*PI)
    np.testing.assert_allclose(reconstructed,low,atol=2e-13)
    scale=reference*math.exp(t)
    return {"M_D_GeV":[scale*math.exp(-l) for l in ld],
            "M_T_GeV":[scale*math.exp(-l) for l in lt],
            "common_trace_coefficient":c,"matching_scale_GeV":scale,
            "inverse_alphas_at_matching":trace.tolist(),
            "geometric_mean_T_over_D":math.exp(sum(ld)/n_d-sum(lt)/n_t)}


def envelopes(rows):
    names=list(rows[0])
    return {k:[min(row[k] for row in rows),max(row[k] for row in rows)] for k in names}


def run():
    old=read("../1369/inputs.json")
    data=read("../1769/inputs.json")
    ref=old["mu_ref_GeV"]
    g=np.array([old[k] for k in ("g1_GUT","g2","g3")])
    lattice=data["strong_lattice_5f"]["ALPHA_2026"]
    top=data["top_MS_matching_GeV"]

    def gauges(values,lattice_source):
        a=4*PI/(g[:,0]+g[:,1]*np.array(values[:3]))**2
        if lattice_source:
            a[2]=1/(lattice[0]+values[2]*lattice[1])+math.log((top[0]+values[3]*top[1])/ref)/(3*PI)
        return a

    result={}
    for key,use_lattice in [("common_MZ",False),("independent_lattice",True)]:
        points=[trial(gauges(ind,use_lattice)) for ind in itertools.product((-2,2),repeat=4)]
        outer=envelopes(points)
        for name in ("C_minus_K_for_one_D","single_complete_pair_missing_K",
                     "one_D_two_T_missing_K","two_D_one_T_log_min_geometric_mean_D_over_T"):
            assert outer[name][0]>0
        central=gauges([0]*4,use_lattice)
        t,x0,k=coefficients(central)
        ld=[3.,3.]
        lt=[3*PI*(x0-7*k)/8+7*sum(ld)/8]
        split=example(central,ref,ld,lt)
        # Two complete pairs also work; choose equal doublet masses and equal triplet masses.
        physical_ld=math.log(ref*math.exp(t)/1000.)
        equal_pairs=example(central,ref,[physical_ld]*2,[physical_ld-3*PI*k/2]*2)
        result[key]={"two_marginal_box":outer,"central":trial(central),
                     "two_doublets_one_triplet":split,"two_split_complete_pairs":equal_pairs,
                     "three_minimum_matter_witness_across_box":all(
                         0<=row["two_D_one_T_log_min_geometric_mean_D_over_T"]+7*6/8<=row["log_scale"]
                         and row["X"]-((row["two_D_one_T_log_min_geometric_mean_D_over_T"]+7*6/8)/(3*PI))>0
                         for row in points)}
    # Work out the minimum max additive inverse-coupling shift that allows the one-D hypothesis.
    # The affine separating row gives the exact unconstrained l-infinity distance to the half-space.
    covector=np.array([trial(row)["C_minus_K_for_one_D"] for row in np.eye(3)])
    central=gauges([0]*4,False)
    gap=covector@central
    correction=-gap/np.linalg.norm(covector,1)*np.sign(covector)
    t,x0,k=coefficients(central+correction)
    np.testing.assert_allclose(x0/6,k,atol=1e-13)
    assert t>0 and x0>0
    result["one_D_boundary_distance_in_inverse_alphas"]={
        "separating_row":covector.tolist(),"central_infimum":max(abs(correction)),
        "min_over_both_boxes":min(result[n]["two_marginal_box"]["C_minus_K_for_one_D"][0]
                                  for n in ("common_MZ","independent_lattice"))/sum(abs(covector)),
        "correction_for_central":correction.tolist()}
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    raw=json.dumps(run(),indent=2,allow_nan=False)+"\n"
    path=HERE/"independent_results.json"
    if parser.parse_args().check:
        assert json.loads(raw)==json.loads(path.read_bytes())
    else:
        path.write_text(raw,encoding="utf-8")
    print("1807 independent",hashlib.sha256(raw.encode()).hexdigest()[:12])
