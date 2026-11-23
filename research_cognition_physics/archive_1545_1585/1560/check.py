"""One-loop gluino threshold, two soft masses, and a color-breaking ray."""
import argparse
import json
from decimal import Decimal, localcontext
from math import log, pi, sqrt
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def integral(a, b, c):
    """Positive three-propagator zero-momentum integral; masses are squared."""
    with localcontext() as ctx:
        ctx.prec = 52
        a, b, c = [Decimal(str(x)) for x in (a, b, c)]
        if min(a,b,c) <= 0:
            raise ValueError("nonpositive squared mass")
        if a == b == c:
            return float(1/(2*a))
        if a == b:
            return float((a-c+c*(c/a).ln())/(a-c)**2)
        def divided(x):
            return x*(x/c).ln()/(x-c) if x != c else Decimal(1)
        return float((divided(a)-divided(b))/(a-b))


def solve():
    inp = load(HERE/"inputs.json")
    prior = load(HERE/inp["previous"])
    g3 = load(HERE/inp["sources"]["MZ_g3"])["g3"][0]
    M, glu, mu, tan = (inp[k] for k in ("scale_TeV","gluino_TeV","mu_TeV","tan_beta"))
    v, mz, sw = (inp[k] for k in ("vev_174_convention_TeV","MZ_TeV","sin2thetaW_test"))
    cosbeta = 1/sqrt(1+tan**2)
    cos2beta = (1-tan**2)/(1+tan**2)
    vd = v*cosbeta
    # Between MZ and the gluino the light color spectrum is that of the SM.
    alpha_inv = 4*pi/g3**2 + (7*log(glu/mz)+5*log(M/glu))/(2*pi)
    alpha = 1/alpha_inv
    prefactor = 2*alpha*glu/(3*pi)
    eta = prefactor*integral(M*M,M*M,glu*glu)
    uniform = 1+eta*mu*tan
    (d, ed),(s, es) = inp["down_strange_SM_table_100TeV"]
    widths = inp["widths"]
    ratio_lo = (d-widths*ed)/(s+widths*es)
    c = inp["analytic_all_gluino_bound"]
    a_cap = 3+c["abs_EW_vacuum_energy_ceiling_TeV4"]*c["max_Yukawa"]**2/(2*M**4)
    etaM_cap = c["alpha_s_upper"]/6
    hcap = etaM_cap*a_cap
    bfloor = 1-etaM_cap*c["abs_mu_tanbeta_over_M_upper"]
    analytic_max = (bfloor+hcap)/(bfloor-hcap)
    out = {"alpha_s_M_one_loop":alpha,"eta_per_TeV":eta,"M_eta":M*eta,
           "universal_gain_b":uniform,"ratio_100TeV_two_width_lower":ratio_lo,
           "all_gluino_absolute_stability_ratio_upper":analytic_max,
           "specified_point_ray_condition_ratio_upper":(uniform+eta*M*a_cap)/(uniform-eta*M*a_cap),
           "cases":{}}
    for name, row in prior["cases"].items():
        pre = row["full"]["low_ratios"][0]
        central = d/(s*pre)
        lower = ratio_lo/pre
        robust = lower/inp["additional_differential_running_window"][1]
        min_a = uniform*(central-1)/(eta*(central+1))
        # Equal and opposite a minimizes the larger of the two |T/Y|.
        gains = [uniform+eta*min_a,uniform-eta*min_a]
        assert np.isclose(gains[0]/gains[1],central)
        # A second choice preserves the measured strange Yukawa separately.
        Y = s*pre/cosbeta
        tree_mass = vd*Y
        def finite_gain(a_d):
            ll = M**2+tree_mass**2+(-.5+sw/3)*cos2beta*mz**2
            rr = M**2+tree_mass**2-sw/3*cos2beta*mz**2
            lr = tree_mass*(a_d-mu*tan)
            eigen = np.linalg.eigvalsh([[ll,lr],[lr,rr]])
            k = prefactor*integral(*eigen,glu**2)
            return 1+k*(mu*tan-a_d),eigen
        linear_a = (uniform-central)/eta
        a,b = 1.01*linear_a,.99*linear_a
        assert finite_gain(a)[0]>central>finite_gain(b)[0]
        for _ in range(52):
            middle = (a+b)/2
            if finite_gain(middle)[0]>central:
                a=middle
            else:
                b=middle
        a_d = (a+b)/2
        gain,eigen=finite_gain(a_d)
        trilinear=Y*a_d
        # At this finite field the negative cubic is twice the quadratic.
        field=3*M*M/abs(trilinear)
        potential=3*M*M*field**2-2*abs(trilinear)*field**3+3*Y**2*field**4
        assert abs((gain*Y*cosbeta)/d-1)<1e-12
        assert min(eigen)>0 and potential < -1
        assert robust>analytic_max
        out["cases"][name]={"input_d_over_s":pre,"MZ_task_gain_lower":row["minimum_down_singular_ratio_gain"],
            "local_required_gain_two_width_lower":lower,"local_gain_to_two_central":central,
            "gain_lower_with_10pct_additional_transfer":robust,
            "minimum_favorable_transport_for_all_gluino_stable_class":lower/analytic_max,
            "exact_infimum_max_abs_a_TeV_for_central_unbroken_one_loop":min_a,
            "strange_fixed":{"Yd_before":Y,"Ys_before":s/cosbeta,"a_s_TeV":mu*tan,
                 "a_d_linear_TeV":linear_a,"a_d_broken_mass_check_TeV":a_d,
                 "T_d_TeV":trilinear,"T_s_TeV":(s/cosbeta)*mu*tan,
                 "squark_masses_TeV":np.sqrt(eigen).tolist(),
                 "insertion_fraction":abs(tree_mass*(a_d-mu*tan))/M**2,
                 "after_d_and_s_SM_like":[gain*Y*cosbeta,s],
                 "finite_field_over_M":field/M,"V_at_field_TeV4":potential}}
    assert np.isclose(integral(1,1,1),.5)
    return out


def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a: compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-9,atol=1e-11)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    ans=solve()
    dest=HERE/"results.json"
    if args.write_new:
        with dest.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(ans,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(load(dest),ans)
    print("1560 flavor thresholds and finite ray comparison verified")
