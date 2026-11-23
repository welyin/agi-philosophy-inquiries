"""Two finite seesaw blocks; collider-independent normalization and kaon yield ceilings."""

import argparse
import json
from math import cosh, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def box(row, width):
    return row[0]-width*row[1], row[0]+width*row[-1]


def phase_factor(parent, lepton, neutral):
    x, y = (lepton/parent)**2, (neutral/parent)**2
    return ((x+y)-(x-y)**2)*sqrt((1-x-y)**2-4*x*y)/(x*(1-x)**2)


def calculate():
    c = json.loads((HERE/"inputs.json").read_text(encoding="utf8"))
    j = json.loads((HERE/c["oscillation_JUNO"]).read_text(encoding="utf8"))["juno_207_days"]
    h = json.loads((HERE/c["oscillation_reactor_and_accelerator"]).read_text(encoding="utf8"))
    widths, M_MeV = c["half_widths"], c["reference_bare_heavy_mass_MeV"]
    M = M_MeV*1e6
    t = sqrt(2)-1
    r = t*t
    def light_pair(gap31):
        m3 = sqrt(gap31)
        m2 = 2*r*m3*(1+m3/M)/(1+sqrt(1+4*r*(m3/M)*(1+m3/M)))
        assert abs(m2*(M+m2)/(m3*(M+m3))-r)<1e-14
        return np.array([m2,m3])
    solar = box(j["solar_21"],widths)
    def atmospheric_given_solar(gap21):
        m2=sqrt(gap21)
        m3=(2*m2*(1+m2/M)/r)/(1+sqrt(1+4*(m2/M)*(1+m2/M)/r))
        return m3*m3
    atm=box(j["normal_atmospheric_31"],widths)
    joint=(max(atm[0],atmospheric_given_solar(solar[0])),
           min(atm[1],atmospheric_given_solar(solar[1])))
    assert joint[0]<=joint[1]
    pair_lo,pair_hi=map(light_pair,joint)
    xlo = pair_lo/(M+2*pair_lo)
    xhi = pair_hi/(M+2*pair_hi)
    s12lo,s12hi = box(c["JUNO_sin2_theta12_207d"],widths)
    s13lo,s13hi = [(1-sqrt(1-a))/2 for a in box(h["DayaBay_amplitude"],widths)]
    s23lo,s23hi = box(h["NO"]["atmospheric_sine_squared"],widths)
    e_lo=np.array([s12lo*(1-s13hi),s13lo])
    e_hi=np.array([s12hi*(1-s13lo),s13hi])
    mu2_lo=max(0,sqrt((1-s12hi)*(1-s23hi))-sqrt(s12hi*s23hi*s13hi))**2
    mu2_hi=min(1,(sqrt((1-s12lo)*(1-s23lo))+sqrt(s12hi*s23hi*s13hi))**2)
    mu_lo=np.array([mu2_lo,s23lo*(1-s13hi)])
    mu_hi=np.array([mu2_hi,s23hi*(1-s13lo)])
    # The relative block weights are fixed much more tightly than either absolute mass.
    q_weight_lo=r*(1-4*pair_hi[1]/M)
    q_weight_hi=r*(1+4*pair_hi[1]/M)
    assert e_lo[0] > e_lo[1] and e_hi[0]>e_hi[1]
    electron_share=[(q_weight_lo*e_lo[0]+e_lo[1])/(1+q_weight_lo),
                    (q_weight_hi*e_hi[0]+e_hi[1])/(1+q_weight_hi)]
    transfer=c["additional_rate_transfer_relative"]
    bounds={}
    for channel, mixlo, mixhi in (("electron",xlo@e_lo,xhi@e_hi),("muon",xlo@mu_lo,xhi@mu_hi)):
        brlo,brhi=box(c[channel]["normalizing_Kenu" if channel=="electron" else "normalizing_Kmunu"],widths)
        kaons=box(c[channel]["effective_kaons"],widths)[1]
        rho=phase_factor(c["kaon_mass_MeV_reference"],c[channel+"_mass_MeV_reference"],M_MeV)
        lo=brlo*rho*mixlo/(1-mixlo)*(1-transfer)
        hi=brhi*rho*mixhi/(1-mixhi)*(1+transfer)
        bounds[channel]={"weak_mixing":[float(mixlo),float(mixhi)],
                         "rho":rho, "branch":[float(lo),float(hi)],
                         "ceiling_at_unit_acceptance":float(kaons*hi)}
    e=c["electron"]
    fluct=sqrt(e["observed_bin"]+e["background_standard_error"]**2)
    n_e=bounds["electron"]["ceiling_at_unit_acceptance"]
    el_br,mu_br=bounds["electron"]["branch"],bounds["muon"]["branch"]
    spectrum=light_pair(j["normal_atmospheric_31"][0])
    # A complex orthogonal rotation retains C5 while changing C6 at leading order.
    z=2.
    rotation=np.array([[np.cosh(z),1j*np.sinh(z)],[-1j*np.sinh(z),np.cosh(z)]])
    d = np.diag(np.sqrt(spectrum*M)) @ rotation
    assert np.allclose(d@d.T/M,np.diag(spectrum),atol=1e-13)
    extra=np.trace(d@d.conj().T).real/M/sum(spectrum)
    assert abs(extra-cosh(2*z))<1e-12
    return {"mixing_and_kaon_ranges":bounds,"electron_fraction":electron_share,
            "electron_to_muon_branch_ratio":[el_br[0]/mu_br[1],el_br[1]/mu_br[0]],
            "physical_heavy_splitting_eV_at_central":float(spectrum[1]-spectrum[0]),
            "max_e_count_over_published_fluctuation_scale":n_e/fluct,
            "enhancement_needed_to_equal_central_bin_excess":(e["observed_bin"]-e["background"])/n_e,
            "same_C5_trace_C6_enhancement_example":extra}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".11g"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=rounded(calculate())
    target=HERE/"results.json"
    if args.write:
        with target.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert result==json.loads(target.read_text(encoding="utf8"))
    print(json.dumps(result))
