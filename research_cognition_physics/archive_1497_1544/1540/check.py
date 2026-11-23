"""Oscillation, charged-current and two-loop bounds for the shared electron source."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
osc=json.loads((HERE/cfg["oscillation_inputs"]).read_text(encoding="utf-8"))


def limit(transfer, mass_transfer, weak_transfer):
    n=cfg["widths"]
    amp,_,d_amp=osc["Daya_Bay_A13"]
    p3=(1-math.sqrt(1-(amp+n*d_amp)))/2
    s,ds=cfg["IO_solar_angle"]
    atm,d_atm,_=osc["T2K_mu_disapp_dm13_IO_eV2"]
    mass_eV=(1-p3)*(1-2*(s+n*ds))*math.sqrt(atm-n*d_atm)
    tau=cfg["leptonic_tau_over_e_HFLAV"]
    mu=cfg["leptonic_mu_over_e_HFLAV"]
    lower=[(v-n*dv)/transfer for v,dv in [tau,mu]]
    denominator=sum(1/r**2 for r in lower)
    assert 2<denominator<4
    coeff=1-denominator/4
    s_upper=(denominator-2)/(1+math.sqrt(1+coeff*(denominator-2)))
    contact=math.sqrt(2)*weak_transfer*cfg["muon_Fermi_constant_GeVm2"]*s_upper
    unit=cfg["electron_pole_mass_GeV"]**2/(48*math.pi**2)
    min_product=mass_eV*1e-9/(mass_transfer*unit*contact)
    cost=mass_eV*1e-9/(mass_transfer*unit*math.sqrt(contact))
    budget=(cfg["illustrative_max_f_norm"]*cfg["illustrative_max_g_ee"]*
            cfg["illustrative_mu_over_Meff"])
    assert abs(2+2*s_upper+s_upper**2-denominator*(1+s_upper**2/4))<1e-14
    return {"B_min_eV":mass_eV,"p3_max":p3,"r_tau_e_and_mu_e_lower":lower,
            "C_upper":denominator,"s_upper":s_upper,"F_squared_over_Meff_squared_GeVm2_max":contact,
            "mu_times_g_ee_GeV_lower":min_product,
            "required_F_g_ee_mu_over_Meff":cost,
            "illustrative_budget_B_upper_eV":mass_transfer*unit*math.sqrt(contact)*budget*1e9,
            "fraction_of_required_block_max_under_budget":budget/cost}


def radial_integral(nodes, weights, masses_sq, mixture, ksq):
    # Split the radial quadrant at y=x; the four-dimensional angular mean is analytic.
    x=(nodes/(1-nodes))[:,None]
    y=x*nodes[None,:]
    def prop(z):
        return sum(frac/(z+m2) for frac,m2 in zip(mixture,masses_sq))
    d=x+y+ksq
    disc=(x-y)**2+2*ksq*(x+y)+ksq**2
    angular=2/(d+np.sqrt(disc))
    jac=2*x/(1-nodes[:,None])**2
    return float(np.sum((weights[:,None]*weights[None,:])*jac*prop(x)*prop(y)*angular))


def run():
    benchmark={}
    for order in [96,192]:
        nodes,weights=np.polynomial.legendre.leggauss(order)
        nodes=(nodes+1)/2
        weights=weights/2
        single=radial_integral(nodes,weights,[1],[1],0.)
        mixed=radial_integral(nodes,weights,[.64,4.],[.35,.65],2.25)
        bound=(.35/.64+.65/4)*math.pi**2/3
        assert abs(single-math.pi**2/3)<.001 and 0<mixed<bound
        benchmark[str(order)]={"massless_central_single":single,"mixed_integral":mixed,
                               "rigorous_mixed_upper":bound}
    assert abs(benchmark["96"]["mixed_integral"]-benchmark["192"]["mixed_integral"])<.001
    all_budgets=limit(cfg["lower_to_true_ratio_transfer_upper"],
                      cfg["physical_mass_block_over_two_loop_bound_upper"],
                      cfg["G0_over_Gmu_upper"])
    # Direction-specific exact tree rates underlying the two-ratio inequality.
    for s in [0,.001,.02,1,100]:
        for x in np.linspace(0,1,51):
            y=1-x
            dtm=1+s*s*x*y
            dte=1+2*s*x+s*s*x
            dme=1+2*s*y+s*s*y
            assert dtm<=min(dte,dme)+1e-12
            np.testing.assert_allclose(dte+dme,2+2*s+s*s,rtol=1e-13)
    return {"ideal_matching":limit(1,1,1),"finite_transfers":all_budgets,
            "quadrature_in_units_of_16pi2_squared":benchmark}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=1e-9,atol=1e-12)


if __name__=="__main__":
    cmd=argparse.ArgumentParser()
    group=cmd.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=cmd.parse_args()
    values=run()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(values,stream,indent=2)
            stream.write("\n")
    else:
        compare(values,json.loads(path.read_text(encoding="utf-8")))
    print("1540 OK",values["finite_transfers"])
