"""Coloured completion of a cubic bottom mass: determinant, finite loops, and holdout."""
import argparse
import json
from math import sqrt
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8"))


def run():
    cfg=read("inputs.json")
    old=read(cfg["production"])
    radial=read(cfg["radial"])
    # Reuse the finite one-loop convention and the exact masses of the four-lepton comparison.
    form,masses=runpy.run_path(str(HERE.parent/"1377/check.py"))["inherited_spectrum"]()
    mh,v=cfg["h_GeV"],cfg["v_GeV"]
    sm=sum(form(mh,m)[1] for m in masses)
    bottom=form(mh,masses[-1])[1]
    assert np.allclose([sm.real,sm.imag],old["SM_quark_loop_real_imag"])
    mstar=cfg["min_heavy_singular_GeV"]
    z=(mh/(2*mstar))**2
    error=(2*cfg["max_insertion_over_v"]*v/mstar)*(14/45)*z/(1-z)
    f0=4/3
    k=radial["bb"]["kappa_upper"]
    amp_real=sm.real+(k-1)*bottom.real+(3-k)*f0-error
    bb_fraction=old["central_upper_bb_fraction"]
    transfer=cfg["min_gg_ratio_transfer"]
    lower=(1-bb_fraction)*transfer*amp_real**2/abs(sm)**2
    observed=old["marginal_rectangles"]["2"]["upper_ratio"]
    quark_gg_max=abs(sm)*sqrt(observed/((1-bb_fraction)*transfer))
    degree_max=k+(quark_gg_max-sm.real-(k-1)*bottom.real+error)/f0
    assert lower>observed and degree_max<3
    # Nonempty hold-in model: fit the running bottom singular value, never using the ggF data.
    def point(edge):
        yuk=np.array([[0.,edge,0],[0,0,700.],[edge,-700.,0.]])
        bare=np.diag([0.,300.,300.])
        left,m,right=np.linalg.svd(bare+yuk)
        kappas=np.diag(left.T@yuk@right.T)/m
        return m,kappas,yuk,bare
    lo,hi=0.,60.
    for _ in range(80):
        mid=(lo+hi)/2
        if point(mid)[0][-1] < masses[-1]:
            lo=mid
        else:
            hi=mid
    mass,kap,yuk,bare=point(hi)
    assert np.isclose(mass[-1],masses[-1]) and np.all(mass[:2]>=mstar)
    assert np.isclose(sum(kap),3) and 0<kap[-1]<k
    assert np.linalg.norm(yuk,2)<=cfg["max_insertion_over_v"]*v
    shift=sum(kr*form(mh,m)[1] for kr,m in zip(kap[:2],mass[:2]))
    example=abs(sm+(kap[-1]-1)*bottom+shift)**2/abs(sm)**2
    for t in [.8,1.15]:
        assert np.isclose(np.linalg.det(bare+t*yuk)/np.linalg.det(bare+yuk),t**3)
    # Tree-level gauge-invariant extra H†H G² amplitude can erase the real residual.
    restored=sm+(kap[-1]-1)*bottom+shift
    cancelling_contact=sm.real-restored.real
    return {"bottom_scalar_upper":k,"finite_heavy_remainder":error,"real_amplitude_floor":amp_real,
            "production_ratio_lower_after_transfer":lower,"two_width_data_upper":observed,
            "maximum_total_log_derivative":degree_max,
            "linear_determinant_fraction_lower":(3-degree_max)/2,
            "admitted_example":{"running_singular_GeV":list(mass),"kappas":list(kap),
                                "outer_mass_insertion_GeV":hi,"norm_insertion_over_v":np.linalg.norm(yuk,2)/v,
                                "pure_gg_rate_over_SM":example,
                                "extra_real_gluon_contact_for_SM_real_part":cancelling_contact,
                                "restored_rate_over_SM":abs(restored+cancelling_contact)**2/abs(sm)**2}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    result=run()
    if parser.parse_args().write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as file:
            json.dump(result,file,indent=2)
            file.write("\n")
    else:
        previous=read("results.json")
        def verify(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for key in a:
                    verify(a[key],b[key])
            else:
                np.testing.assert_allclose(a,b,atol=1e-11,rtol=1e-10)
        verify(result,previous)
        print("1439 determinant index, finite coloured loops and four-lepton holdout checked")
