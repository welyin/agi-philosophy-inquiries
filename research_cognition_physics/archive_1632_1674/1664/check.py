"""Scalar mass splittings, determinant bounds, and a bounded invariant potential."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
Q = Fraction


def read(rel):
    return json.loads((HERE/rel).read_text(encoding="utf-8"))


def exact_spectrum(central, scale):
    # Rationals are essential: the light squared mass is ~10^-28 of the heavy squared mass.
    def fraction_mass(value): return (Q(str(value))/Q(str(scale)))**2
    trip = fraction_mass(central["common_triplet_mass_GeV"])
    sext = fraction_mass(central["common_sextet_mass_GeV"])
    light = list(map(fraction_mass, central["doublet_masses_GeV"]))
    rows = []
    for mid in light:
        cross = (trip+sext-2*mid)/25
        slope = (sext-mid)/5-2*cross
        offset = mid+slope+6*cross
        def mass(x,y): return offset+(x+y)*slope+x*y*cross
        assert [mass(-3,-3),mass(-3,2),mass(2,2)] == [trip,mid,sext]
        rows.append((offset,slope,cross))
    return rows


def numeric_spotchecks():
    rng=np.random.default_rng(1664)
    def positive(n):
        x=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        return x.conj().T@x+np.eye(n)
    def logdet(a):
        sign,val=np.linalg.slogdet(a)
        assert abs(sign-1)<2e-13
        return val
    samples=0
    for n in (1,2,4,6):
        for _ in range(6):
            trip,mid,sext = [positive(n) for _ in range(3)]
            mix=(trip+sext)/2
            assert 2*logdet(mix)+1e-12 >= logdet(trip)+logdet(sext)
            cross=(trip+sext-2*mid)/25
            slope=(sext-mid)/5-2*cross
            offset=mid+slope+6*cross
            for x,y,target in ((2,2,sext),(-3,2,mid),(-3,-3,trip)):
                assert np.allclose(offset+(x+y)*slope+x*y*cross,target)
            eig,u=np.linalg.eigh(mid)
            inv=(u/np.sqrt(eig))@u.conj().T
            rel=np.linalg.eigvalsh(inv@(25*cross)@inv)
            d,s=[math.exp((logdet(a)-logdet(mid))/n) for a in (trip,sext)]
            assert rel[-1] >= d+s-2-1e-12
            samples+=1
    return samples


def report():
    cfg=read("inputs.json")
    saved=read(cfg["thresholds"])
    central=saved["central_witness"]
    mz=read(cfg["gauge_input"])["mu_ref_GeV"]
    scale=read(cfg["prior_contract"])["matching_GeV"]
    log_range=math.log(scale/mz)
    # Conservative outward lower bounds from the prior analytically bounded data rectangle.
    x,y=122.3516,105.2576
    assert x < saved["sum_log_triplet_over_doublet"][0]
    assert y < saved["sum_log_sextet_over_doublet"][0]
    relative_one=math.exp(x/2)+math.exp(y/2)-2
    relative_all=math.exp(2*(x-3*log_range))+math.exp(2*(y-3*log_range))-2
    rows=exact_spectrum(central,scale)
    # Each rational coefficient above uses MU units; transfer to units of the chosen vacuum v.
    g2=4*math.pi/central["common_inverse_alpha"]
    v=scale/(5*math.sqrt(g2))
    norm=25*g2
    scaled=np.array(rows,dtype=float)*norm
    slope_norm=max(abs(scaled[:,1]))
    cross_norm=max(abs(scaled[:,2]))
    lip=(4/5)*(slope_norm+3*cross_norm)+cross_norm
    quartic=cfg["benchmark_stabilizer_lambda_15"]
    assert 4*g2*quartic > lip**2
    masses=central["doublet_masses_GeV"]
    ratios=[25*float(row[2])*scale**2/m**2 for row,m in zip(rows,masses)]
    assert min(ratios)>relative_all and max(ratios)>relative_one
    # Adjoint block derivatives of h^2+h-6 at the high-scale vacuum.
    vac=[2,2,2,-3,-3]
    assert sum(vac)==0
    slope_sq=[(a+b+1)**2 for a in vac for b in vac]
    assert sorted(set(slope_sq))==[0,25]
    assert slope_sq.count(0)==12 and slope_sq.count(25)-1==12
    assert math.isclose(5*math.sqrt(g2)*v, scale)
    return {"mixed_matrix_samples":numeric_spotchecks(),"four_copy_relative_cross_lower":{
                "one_direction":relative_one,"every_direction":relative_all},
            "if_v_equals_MU_absolute_C_lower":{"one_direction":relative_one*mz**2/(25*scale**2),
                                                 "every_direction":relative_all*mz**2/(25*scale**2)},
            "central_potential":{"unification_scale_GeV":scale,"v_GeV":v,"lambda_Sigma":g2,
                   "lambda_Phi":quartic,"stability_Lipschitz_bound":lip,
                   "stabilizer_margin":4*g2*quartic-lip**2,
                   "four_dimensionless_cross_coefficients":[float(c)*norm for _,_,c in rows],
                   "ratios_25v2C_over_B":ratios,"rounded_data_reconstruction":"exact fractions",
                   "physical_adjoint_mass_GeV":5*math.sqrt(g2)*v,
                   "broken_vector_mass_GeV":5*math.sqrt(g2)*v}}


def tidy(obj):
    if isinstance(obj,dict): return {k:tidy(v) for k,v in obj.items()}
    if isinstance(obj,list): return list(map(tidy,obj))
    if isinstance(obj,float): return float(f"{obj:.11g}")
    return obj


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    opts=parser.add_mutually_exclusive_group(required=True)
    opts.add_argument("--check",action="store_true")
    opts.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=tidy(report())
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as dest:
            dest.write(json.dumps(data,indent=2)+"\n")
    else:
        assert data == json.loads(path.read_bytes())
        print("PASS 1664 separability, forced Hermitian cross term, exact spectrum and stable adjoint-breaking vacuum")
