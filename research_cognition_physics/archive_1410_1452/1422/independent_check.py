"""Direct Z-pole ratio/asymmetry inference and a muon-DAR/IBD normalization check."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def unit_width(asym, alpha_inverse):
    ratio=asym/(1+math.sqrt(1-asym**2))
    return 2/((1+ratio**2)*(1+3/(4*math.pi*alpha_inverse)))


def recoil_example():
    # A full six-state collection of orthonormal weak ports; the first three are accessible.
    # An arbitrary complex rotation selects a non-maximal overlap, distinct from the author's 120 MeV example.
    eigen=np.array([.94,.987,.999])
    j=np.arange(3)
    rotation=np.exp(2j*math.pi*np.outer(j,j)/3)/math.sqrt(3)
    n=rotation@np.diag(np.sqrt(eigen))
    ports=np.column_stack([n,rotation@np.diag(np.sqrt(1-eigen))])
    assert np.allclose(ports@ports.conj().T,np.eye(3))
    a=n@n.conj().T
    e,mu=np.real(np.diag(a)[:2])
    gf0,vud=1.09,.974
    gfmu=gf0*math.sqrt(e*mu)
    beta_report=vud/math.sqrt(mu)
    from_lifetime=gf0**2*vud**2*e
    from_mu_beta=gfmu**2*beta_report**2
    source=n[1,:]/math.sqrt(mu)
    cross_section=gf0**2*vud**2*abs(np.vdot(n[0,:],source))**2
    normalized=abs(a[0,1])**2/(e*mu)
    assert np.isclose(from_lifetime,from_mu_beta) and np.isclose(cross_section/from_lifetime,normalized)
    return {"Aee":e,"Amumu":mu,"cross_section_normalized_by_beta":cross_section/from_lifetime,
            "normalized_overlap":normalized,"fixed_true_quark_input_instead_of_beta_factor":1/mu}


def calculate():
    inputs=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    # Complete LEP/SLD report eqs. (7.7), (3.15) and (1.39).
    r,err=5.943,.016
    asym,asym_err=.1513,.0021
    widths=inputs["matching_widths"]
    s0=unit_width(asym,128)
    max_single=unit_width(asym-widths*asym_err,128)
    robs=r-widths*err
    lsnd=(inputs["lsnd_percent"][0]-widths*math.hypot(*inputs["lsnd_percent"][1:]))/100
    arc=1.26693268*inputs["span_eVsquared"]*inputs["example_max_L_over_E_m_per_MeV"]
    initial_min=math.sin(math.asin(math.sqrt(lsnd))-arc)
    max_u=2+((1-initial_min)/(1+initial_min))**2
    phases={}
    for k in [1.,1.01,1.02]:
        lower=robs/(max_single*k)
        zero=(1-math.sqrt(lower-2))/(1+math.sqrt(lower-2))
        upper=math.sin(math.asin(zero)+arc)**2
        phases[str(k)]={"U_lower":lower,"finite_P_upper":upper}
    return {"pole_Rinv":r,"SLD_Aell":asym,"single_flavour_response_at_k_one":s0,
            "inferred_U_at_k_one_central":r/s0,"published_full_SM_width_ratio":1.99125,
            "effective_k_matching_original_SM_width_at_central":1.99125/s0,
            "LSND_two_width_floor":lsnd,"critical_U_for_assumed_finite_window":max_u,
            "minimum_k_to_allow_appearance_no_other_invisible":robs/(max_single*max_u),
            "minimum_extra_Rinv_in_lepton_width_units_at_k_one":robs-max_single*max_u,
            "rho_sensitivity":phases,"cross_calibration_example":recoil_example()}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=parser.parse_args()
    actual=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(actual)
    elif json.loads(path.read_text(encoding="utf-8"))!=json.loads(actual):
        raise SystemExit("1422 independent replay differs")
    print("1422 raw pole normalization, charged response and calibrated DAR event ratio checked")
