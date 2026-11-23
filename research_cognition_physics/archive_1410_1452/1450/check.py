"""Deterministic flavor-volume bounds on two axion properties at fixed matching."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf8"))


def ckm_from_angles(angles):
    a,b,c,phase = angles
    x,y,z = np.sin([a,b,c])
    ca,cb,cc = np.cos([a,b,c])
    v = np.exp(1j*phase)
    return np.array([[ca*cb,x*cb,y/v],
                     [-x*cc-ca*z*y*v, ca*cc-x*z*y*v, z*cb],
                     [x*z-ca*cc*y*v,-ca*z-x*cc*y*v,cc*cb]])


def volume_ranges(sectors, widths):
    logs = []
    for sector in sectors:
        y = np.array(sector)
        assert np.all(y[:,0]-widths*y[:,1] > 0)
        logs.append((-float(np.log(y[:,0]+widths*y[:,1]).sum()),
                     -float(np.log(y[:,0]-widths*y[:,1]).sum())))
    return logs


def ratio(zu,zd,ze):
    return 8/3 + 2*(ze-zd)/(zu+zd)


def bounds(sectors, cfg, qcd):
    logs = volume_ranges(sectors, cfg["mass_widths"])
    h = 3*math.log(cfg["coefficient_singular_ceiling"])
    intervals = [(lo-h,hi+h) for lo,hi in logs]
    assert all(lo>0 for lo,_ in intervals)
    possible = [ratio(*z) for z in itertools.product(*intervals)]
    mix = qcd["photon_subtraction"]
    err = qcd["photon_subtraction_width"]
    c = 8/3 - (mix+err)
    min_h = (c*logs[0][0]+2*logs[2][0]-(2-c)*logs[1][1])/4
    exact_edge = ratio(logs[0][0]-min_h,logs[1][1]+min_h,logs[2][0]-min_h)
    assert abs(exact_edge-(mix+err)) < 1e-12
    min_one = math.exp(min(logs[0][0],logs[1][0])/3)
    factor = 1/(2*math.pi*qcd["inverse_alpha_zero_rounded"])
    mf_lo = qcd["m_times_f_GeV2"]-sum(qcd["m_times_f_errors_GeV2"])
    mf_hi = qcd["m_times_f_GeV2"]+sum(qcd["m_times_f_errors_GeV2"])
    return {"log_volume_ranges":logs,"E_over_N":[min(possible),max(possible)],
            "necessary_B_for_photon_cancellation":math.exp(min_h/3),
            "necessary_B_for_unit_domain_wall":min_one,
            "abs_g_gamma_over_mass_GeV_minus2":[factor*(min(possible)-mix-err)/mf_hi,
                                                factor*(max(possible)-mix+err)/mf_lo]}


def make_witness(sectors, V, epsilon, spectators):
    spectra = [np.array(rows)[:,0] for rows in sectors]
    charges = [np.rint(np.log(vals)/math.log(epsilon)).astype(int) for vals in spectra]
    # Q and L doublets have zero global charges. CKM is in the down left core.
    models = []
    for shift in (0,1):
        exps = [q.copy() for q in charges]
        exps[2][0] += shift
        cores = [frame @ np.diag(vals/epsilon**q)
                 for frame,vals,q in zip([np.eye(3),V,np.eye(3)],spectra,exps)]
        ys = [core@np.diag(epsilon**q) for core,q in zip(cores,exps)]
        for y, target in zip(ys,spectra):
            assert np.allclose(np.linalg.svd(y,compute_uv=False)[::-1],target,atol=1e-14)
        assert np.allclose(ys[1],V@np.diag(spectra[1]))
        dd = sum(int(x.sum()) for x in exps[:2])
        e_over_n = ratio(*(float(x.sum()) for x in exps))
        vals_c = np.concatenate([np.linalg.svd(c,compute_uv=False) for c in cores])
        assert vals_c.min() > 1/3 and vals_c.max() < 3
        models.append({"right_charges":[x.tolist() for x in exps],
                       "coefficient_singular_extrema":[float(vals_c.min()),float(vals_c.max())],
                       "D":dd,"E_over_N":e_over_n,
                       "with_extra_charged_pairs":e_over_n-2*sum(q*q for q in spectators)/dd})
    return models


def run():
    cfg = read("inputs.json")
    low = read(cfg["low_masses"])["yukawa"]
    high = read(cfg["high_quarks"])["high_scale"]
    high_rows = [high["up"],high["down"],read(cfg["high_charged_leptons"])["scales"]["1e16"][3:]]
    qcd = read(cfg["qcd"])["qcd"]
    values = {"MZ":bounds([low[k] for k in ("up","down","charged_leptons")],cfg,qcd),
              "1e16":bounds(high_rows,cfg,qcd)}
    v = ckm_from_angles([row[0] for row in high["angles_theta12_theta13_theta23_delta"]])
    assert np.allclose(v@v.conj().T,np.eye(3))
    values["high_scale_witnesses"] = make_witness(high_rows,v,cfg["witness_epsilon"],
                                                  cfg["witness_spectators_em_charges"])
    assert abs(values["high_scale_witnesses"][0]["with_extra_charged_pairs"]-1.92)<0.04
    # Normalize tuples and numpy-derived values to the on-disk representation.
    return json.loads(json.dumps(values))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=p.parse_args()
    results=run()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf8") as out:
            json.dump(results,out,indent=2)
            out.write("\n")
    else:
        assert results==read("results.json")
    print("1450 checked",results["1e16"]["E_over_N"])
