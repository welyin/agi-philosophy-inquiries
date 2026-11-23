"""Separate real shape curvature from the transverse-muon-polarization readout."""
import argparse
import json
from itertools import product
from math import exp, log, log1p, sqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_text(encoding="utf-8"))


def evaluate():
    cfg = read(HERE/"inputs.json")
    masses = cfg["kinematic_GeV"]
    muon = masses["m_mu"]
    strange = read(HERE.parent/"1480"/"inputs.json")["lattice"]["ms"]
    up = read(HERE.parent/"1500"/"inputs.json")["up_mass_MeV"]
    q = (strange[0]-up[0])/1000
    gap_hi = (strange[0]+2*strange[1]-up[0]+2*up[1])/1000
    channels = {}
    for tag,k,pion in [("charged","M_Kplus","M_pi0"),("neutral","M_K0","M_piplus")]:
        delta = masses[k]**2-masses[pion]**2
        tmin,tmax = muon**2,(masses[k]-masses[pion])**2
        # Effective real log C after dividing out the dispersive phase.
        eps = 0.001
        sigma = eps/(muon*q)
        shift0 = delta*sigma
        shift_min = delta*log1p(sigma*tmin)/tmin
        shift_max = delta*log1p(sigma*tmax)/tmax
        assert shift0 > shift_min > shift_max
        # Directly compare an exact scalar contact with a shifted standard shape.
        residue = (1+sigma*tmax)*exp(-sigma*tmax)-1
        channels[tag]={"logC_shift_linear_for_real_epsilon_0p001":shift0,
                       "two_endpoint_effective_logC_difference":shift_min-shift_max,
                       "max_fractional_amplitude_discrepancy_vs_shifted_SM":abs(residue),
                       "scalar_to_Im_xi_per_unit_f0_over_fplus":delta/(muon*q)}
    delta = masses["M_Kplus"]**2 - masses["M_pi0"]**2
    # The external experiment reports |Im xi| for its calibrated event response.
    limit_xi = 0.016
    factors = [1.0,0.85]
    limits = {str(b):limit_xi*muon*gap_hi/(delta*b) for b in factors}
    ktev,ktev_err = cfg["shape_KTeV_neutral"]
    charged,charged_err = cfg["shape_NA48_2_charged"]
    delta_coeff = (masses["M_K0"]**2-masses["M_piplus"]**2-delta)/(muon*q)
    def edges(key):
        x, s = cfg[key]
        return [x-2*s,x+2*s]
    limits_ratio = []
    for axial, nuclear, vector, fkpi, ff in product(
            edges("inclusive_axial_product"),edges("nuclear_vector_product"),
            edges("neutral_electron_vector_product"),edges("FK_over_Fpi_isosymmetric"),
            edges("fplus_neutral")):
        limits_ratio.append(0.5*log(axial*nuclear*ff/(vector*fkpi)))
    ratio_box = [min(limits_ratio),max(limits_ratio)]
    if (HERE/"results.json").exists():
        saved = read(HERE/"results.json")["half_log_Z_over_B"]
        assert max(abs(i-j) for i,j in zip(ratio_box,saved)) < 1e-10
    return {"ms_minus_mu_2GeV_GeV":q,"channels":channels,
            "independent_half_log_rate_lattice_window":ratio_box,
            "E246_abs_Imxi_90_CL":limit_xi,
            "Im_epsilon_nominal_unit_response":limit_xi*muon*q/delta,
            "Im_epsilon_upper_doubled_mass_widths_given_response_lower":limits,
            "response_needed_to_hide_Im_epsilon_0p001":limit_xi*muon*gap_hi/(delta*.001),
            "neutral_minus_charged_Central_quoted_logC":ktev-charged,
            "marginal_widths_neutral_charged":[ktev_err,charged_err],
            "two_channels_delta_logC_linear_contrast_for_0p001":delta_coeff*.001}


def check(x,y):
    if isinstance(x,dict):
        assert x.keys()==y.keys()
        for k in x:
            check(x[k],y[k])
    elif isinstance(x,list):
        for a,b in zip(x,y,strict=True):
            check(a,b)
    else:
        assert abs(x-y) <= 3e-10*max(abs(x),abs(y),1e-10)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    flag=parser.parse_args().check
    result=evaluate()
    target=HERE/"independent_results.json"
    if flag:
        check(result,read(target))
        print("1563 independent scalar shape curvature, charge contrast and E246 spin response checked")
    else:
        target.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
