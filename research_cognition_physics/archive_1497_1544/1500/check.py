"""Fixed-theory mass response, KM ambiguity, and a separately calibrated sea source."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def km(masses,k):
    u,d,s=masses
    return np.array([u+k*d*s,d+k*s*u,s+k*u*d])


def mass_comparisons(masses,lepton_ratio,target):
    u,d,s=masses
    alpha9=(target*s-lepton_ratio*d)/(u*(lepton_ratio*s-target*d))
    alpha0=-u/(d*s)
    return float(lepton_ratio*d/s),float(alpha9),float(alpha0)


def unitary(rng):
    q,r=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))
    return q/np.linalg.det(q)**(1/3)


def adj(a):
    # Only invert nonsingular random test matrices. Analytic identity extends polynomially.
    return np.linalg.det(a)*np.linalg.inv(a)


def action_test():
    rng=np.random.default_rng(1500)
    residual=0.
    for _ in range(40):
        chi=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/3
        u=unitary(rng)
        a=u.conj().T@chi
        x=np.trace(a)
        o6=(x+x.conjugate())**2
        o7=(x-x.conjugate())**2
        o8=np.trace(a@a)+np.conj(np.trace(a@a))
        # Set F0=B0=kappa=1 in the coefficient identity.
        change=.25*(np.trace(.5*adj(chi.conj().T)@u.conj().T)
                    +np.trace(.5*adj(chi)@u))
        residual=max(residual,abs(change-(o6+o7-2*o8)/32))
    assert residual<2e-14
    # Store the tolerance, rather than machine-specific LAPACK rounding.
    return {"complex_field_tests":40,"absolute_tolerance":2e-14}


def ck_projection(archive):
    config=read(archive)
    ref=read(archive.parent/config["standard_model_reference"])
    x12,x13,x23,delta=[pair[0] for pair in ref["quark_sines_and_phase"]]
    c12,c13,c23=[math.sqrt(1-v*v) for v in (x12,x13,x23)]
    phase=np.exp(1j*delta)
    v=np.array([[c12*c13,x12*c13,x13/phase],
                [-x12*c23-c12*x23*x13*phase,c12*c23-x12*x23*x13*phase,x23*c13],
                [x12*x23-c12*c23*x13*phase,-c12*x23-x12*c23*x13*phase,c23*c13]])
    assert np.linalg.norm(v@v.conj().T-np.eye(3))<1e-14
    # Named down direction b, one physical alignment test. General invariance is independent of it.
    return {"down_axis":"b","reduced_D_flavour_factor":float(abs(v[0,2]*v[1,2]))}


def run():
    cfg=read(HERE/"inputs.json")
    base=read(HERE/cfg["same_scale"])
    raw_masses=[cfg["up_mass_MeV"],base["lattice"]["md"],base["lattice"]["ms"]]
    masses=np.array([v[0]/1000 for v in raw_masses])
    muon=base["muonium"]["pole_ratio"]
    target=cfg["comparison_double_ratio"]
    D,k9,k0=mass_comparisons(masses,muon[0],target)
    widths=cfg["work_widths"]
    def edges(pair,scale=1):
        return [(pair[0]+side*widths*pair[1])/scale for side in (-1,1)]
    outcomes=[mass_comparisons(x[:3],x[3],target) for x in
              itertools.product(*[edges(y,1000) for y in raw_masses],edges(muon))]
    boxes=[[min(xs),max(xs)] for xs in zip(*outcomes)]
    m9=km(masses,k9)
    m0=km(masses,k0)
    assert 0 < m9[0] < masses[0] and abs(m0[0])<1e-16
    assert abs(muon[0]*m9[1]/m9[2]-target)<1e-12
    val,stat,sys=cfg["sea_slope_GeV_inv"]
    zbox=[val-widths*(stat+sys),val+widths*(stat+sys)]
    width=zbox[1]-zbox[0]
    low_s,high_s=cfg["sea_pair_central_GeV"]
    def pair_ratio(z):
        return (1+high_s*z)/(1+low_s*z)
    jac=cfg["jacobian_operator_infty_example"]
    slope_min=2-3*jac/(1-jac)
    inverse_tangent_error=jac/(1-jac)
    example_shift=cfg["local_up_shift_example_MeV"]/1000
    scalar_example={}
    for name,k in (("nine",k9),("zero",k0)):
        # The difference between a linear scalar source on the second theory
        # and the push-forward of the first source at the same physical point.
        scalar_example[name]=(1000*(km(masses,k)-masses)).tolist()
    return {"fixed_QCD":{"absolute_error_in_d_and_s_tangents":inverse_tangent_error,
                         "min_abs_K0_squared_slope_in_B0":slope_min,
                         "conditional_delta_K0_squared_GeV2":slope_min*cfg["mass_scale_reference_B0_GeV"]*example_shift},
            "EFT_action":action_test(),"CKM_reduced_witness":ck_projection(HERE/cfg["tree_charged"]),
            "double_ratio_physical":[D,boxes[0]],"kappa9_GeV_inv":[k9,boxes[1]],"kappa0_GeV_inv":[k0,boxes[2]],
            "mass_at_9_MeV":list(1000*m9),"mass_at_zero_MeV":list(1000*m0),
            "linear_source_pushforward_minus_physical_tangent_MeV":scalar_example,
            "sea_slope_box_GeV_inv":zbox,"max_pair_displacement_GeV_inv":width,
            "minimal_unmodelled_difference_for_both":{"nine":-boxes[1][1]-width,"zero":-boxes[2][1]-width},
            "conditional_max_top_mass_fraction_if_other_slope_zero":zbox[1]/(-boxes[2][1]),
            "conditional_min_cancelling_other_slope_GeV_inv":-boxes[2][1]-zbox[1],
            "pion_mass_squared_ratio_104_over_94_in_linear_test":{
              "response_published_two_widths":[pair_ratio(v) for v in zbox],
              "full_topological_mass_hypothesis_at_ratio_5":pair_ratio(5)}}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    modes=p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-new",action="store_true")
    modes.add_argument("--check",action="store_true")
    args=p.parse_args()
    out=run()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    else:
        assert read(path)==out
        print("1500: conditional QCD identifiability, chiral action and independent strange source checked.")
