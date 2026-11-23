"""A common gauge threshold and the neutral current of two type-III triplets."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def predictions(cfg, gauges, box):
    plo, phi = gauges["whole_triplet_geometric_mean_GeV"]
    dlo,dhi = box["abs31_eV2"]
    rlo,rhi = box["solar_over_abs31"]
    light = {"NO": (dlo*math.sqrt(rlo), dhi*math.sqrt(rhi)),
             "IO": (dlo*math.sqrt(1+rlo), dhi*math.sqrt(1+rhi))}
    mmax_GeV = math.sqrt(dhi*(1+rhi))*1e-9
    data = cfg["atlas_muon"]
    pad = data["half_rounding_MeV"]
    lower = data["value_MeV"] - pad - data["width_count"]*math.hypot(data["stat_MeV"]+pad,data["syst_MeV"]+pad)
    ratio = lower/data["leptonic_reference_upper_MeV"]
    lepton_gamma = .5*(1-4*(cfg["muon_mass_upper_GeV"]/90)**2)**1.5
    nc_lower = data["minimum_ratio_transfer"]*ratio*lepton_gamma
    assert 2 < nc_lower < 3
    res = {}
    for label, delta in (("ATLAS_single_muon",3-nc_lower),
                         ("stronger_calibration_benchmark",cfg["normalized_width_deficit_example"])):
        hmax = 1-math.sqrt(1-delta)
        variants = []
        for f in cfg["matching_factors"]:
            floor = cfg["formal_floor_GeV"]/f
            assert (cfg["muon_mass_upper_GeV"]/floor)**2 < .5
            pole_lo = (plo/f)**2
            # Product bound from each bare singular value and the sterile-light block.
            pole_hi = (phi*f*(1+hmax*mmax_GeV/floor)/(1-hmax))**2
            bands = {order:[lo*1e-18/pole_hi,hi*1e-18/pole_lo] for order,(lo,hi) in light.items()}
            assert bands["NO"][1] < bands["IO"][0]
            variants.append(dict(factor=f, pole_product_GeV2=[pole_lo,pole_hi], Z_invariant=bands,
                                 ordering_separation_ratio=bands["IO"][0]/bands["NO"][1]))
        res[label]=dict(width_deficit_upper=delta,hmax=hmax,windows=variants)
    return dict(rounded_lower_width_MeV=lower,inferred_published_ratio_lower=ratio,
                muon_width_ratio_lower=lepton_gamma, corrected_light_Z_trace_square_lower=nc_lower,
                scenarios=res)


def physical_currents():
    rng = np.random.default_rng(1630)
    clebsch = np.diag([1.,1.,1.,math.sqrt(2),math.sqrt(2)])
    active = np.diag([1.,1.,1.,0.,0.])
    two = 2*np.eye(5)
    residual = 0.
    negative = 0
    for _ in range(36):
        b = rng.normal(size=(3,2)) + 1j*rng.normal(size=(3,2))
        sterile = np.diag([6.,13.])
        neutral = np.block([[np.zeros((3,3)), b],[b.T,sterile]])
        charged = np.block([[np.diag([.001,.02,.07]), math.sqrt(2)*b],[np.zeros((2,3)),sterile]])
        # For the neutral convention M=U diag(m) U^T use left
        # singular vectors. Column rephasings leave all checked spectra unchanged.
        un,m,_ = np.linalg.svd(neutral)
        un, m = un[:,::-1],m[::-1]
        uc,_,_ = np.linalg.svd(charged)
        uc=uc[:,::-1]
        v=uc.conj().T @ clebsch @ un
        zn=un.conj().T @ active @ un
        zc=uc.conj().T @ (two-active) @ uc
        c=zn[:3,:3]
        h=np.linalg.eigvalsh(np.eye(3)-c)[1:]
        prod=m[1]*m[2]/(m[3]*m[4])
        relative=abs(np.prod(h)/np.linalg.det(c).real/prod-1)
        error=max(relative,np.max(abs(zn-(two-v.conj().T@v))),
                  np.max(abs(zn-v.conj().T@(two-zc)@v)))
        residual=max(residual,error)
        assert np.all(m[3:] >= np.diag(sterile)-1.e-10)
        assert np.all(np.diag(sterile) >= (1-h[-1])*m[3:] - h[-1]*m[2] - 1.e-10)
        vlight=v[:3,:3]
        if min(np.linalg.eigvalsh(np.eye(3)-vlight@vlight.conj().T)) < -1e-5:
            negative+=1
    assert residual < 1.e-8 and negative > 0
    return dict(trials=36,projector_and_product_errors_below=1e-8,negative_naive_defects=negative)


def run():
    cfg = json.loads((HERE/"inputs.json").read_bytes())
    gauges = json.loads((HERE/cfg["gauge"]).read_bytes())
    box = json.loads((HERE/cfg["light"]).read_bytes())["box"]
    return {"matching":predictions(cfg,gauges,box),"matrix_verification":physical_currents()}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    data=run()
    p=HERE/"results.json"
    if args.check:
        assert data==json.loads(p.read_bytes())
    else:
        p.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
    print("1630 PASS: triplet weak currents, threshold-to-pole bounds, single-muon transport and ordering windows")
