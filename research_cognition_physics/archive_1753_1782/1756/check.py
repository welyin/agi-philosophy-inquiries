"""A mass-independent obstruction between single-Higgs hiding and the frozen QCD trajectory."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def read(file):
    return json.loads((HERE/file).read_bytes())


CFG=read("inputs.json")


def signatures():
    obs=read(CFG["higgs"])
    family=obs["family"]
    mu,md=family["masses_chiral_UDEN"][:2]
    masses=family["vector_poles_including_inactive_weak_partner"]
    joint=family["extra_down_singlet_precision_holdout"]
    m=family["hard_vector_mass"]
    second=[m*(1+offset*q) for q in (joint["q_up"],joint["q_down"]) for offset in (-1,1)]
    return {"one_return_pair":[mu,md,*masses],"two_return_pairs":[mu,md,*second]}


def slope_check(bound, min_kappa):
    # Direct singular perturbation in nonsymmetric real 2x2 blocks. The sampler checks the
    # analytic inequalities; the class bound below does not depend on a scan.
    rng=np.random.default_rng(1756)
    count=0
    fraction=0.
    for _ in range(1400):
        mq,ms=bound*np.exp(rng.uniform(-2.,2.,2))
        a,c=bound*rng.uniform(-1,1,2)
        y=np.array([[0.,a],[c,0.]])
        m=np.array([[mq,a],[c,ms]])
        u,v,vh=np.linalg.svd(m)
        s=v[::-1]
        slopes=(np.diag(u.T@y@vh.T)/v)[::-1]
        np.testing.assert_allclose(sum(s*s*slopes),a*a+c*c,atol=2e-8)
        np.testing.assert_allclose(sum(slopes),-2*a*c/np.linalg.det(m),rtol=3e-9,atol=3e-11)
        if -slopes[0]>=min_kappa:
            count+=1
            cap=bound**2*(min_kappa**-2+2/min_kappa)
            assert np.prod(s)<=cap*(1+1e-12)
            fraction=max(fraction,np.prod(s)/cap)
    assert count>0
    return dict(examples_with_sufficient_light_slope=count,max_area_to_analytic_cap_in_sample=fraction)


def calculate():
    src=read(CFG["gauge"])["new_spectrum"]
    U=src["cross_GeV"]
    a0=src["a_at_cross"][2]
    low=read(CFG["low_boundary"])["mu_ref_GeV"]
    atlas=read(CFG["higgs"])["down_compensator_bound"]
    ymax=4*math.pi/math.sqrt(CFG["Nc"])
    B=ymax*read(CFG["higgs_setup"])["v_GeV"]/math.sqrt(2)
    mmin=read(CFG["higgs_data"])["higgs_GeV"]/2
    max_add=max(edge[1] for edge in atlas["extra_gluon_real_shift_intervals"])
    k=(8/3-max_add)/2
    assert k>0
    area=B**2*(1/k**2+2/k)
    max_vector_mass=area/mmin
    assert max_vector_mass<U
    max_geom=(B**2*area*max_vector_mass)**.2
    # At a fixed measured low boundary count only the portion of each threshold log above it.
    corrected=max_geom*max(1.,low/mmin)
    required=lambda cap:U*math.exp(-3*math.pi*(a0-cap)/5.)
    best_inv=a0-5/(3*math.pi)*math.log(U/corrected)
    rows={}
    for label,ms in signatures().items():
        mu=np.array(ms)
        assert min(mu)>low and max(mu)<U
        logs=sum(np.log(U/mu))
        loss=logs/(3*math.pi)
        rows[label]=dict(masses_GeV=ms,mass_geometric_mean_GeV=math.exp(np.mean(np.log(mu))),
                         inverse_alpha3_shift_at_old_U=-loss,formal_inverse_alpha3_at_old_U=a0-loss,
                         minimal_extra_positive_inverse_matching_for_alpha_le_one=max(0,1-a0+loss))
    assert best_inv<0 and corrected<required(0)<required(1)
    return dict(U_GeV=U,inverse_color_at_old_U=a0,low_reference_GeV=low,yukawa_cap=ymax,
                higgs_entry_cap_GeV=B,min_mass_GeV=mmin,additional_real_gluon_amplitude_max=max_add,
                lighter_state_negative_Higgs_log_slope_required=k,upper_s1s2_GeV2=area,
                upper_each_vector_mass_GeV=max_vector_mass,
                upper_five_pole_geometric_mean_GeV=max_geom,
                conservative_upper_after_low_threshold_clipping_GeV=corrected,
                lower_five_mean_for_positive_one_loop_color_at_U_GeV=required(0),
                lower_five_mean_for_alpha_le_one_at_U_GeV=required(1),
                class_upper_formal_inverse_alpha3_at_old_U=best_inv,
                toy_SVD_crosscheck=slope_check(B,k),specific_examples=rows)


def native(x):
    if isinstance(x,dict):
        return {k:native(v) for k,v in x.items()}
    if isinstance(x,list):
        return [native(y) for y in x]
    if isinstance(x,(np.floating,float)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=native(calculate())
    path=HERE/"results.json"
    if args.check:
        assert out==json.loads(path.read_bytes())
        print("1756 Higgs slope mass bound and cumulative colored threshold budget verified")
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2,allow_nan=False)
            f.write("\n")
