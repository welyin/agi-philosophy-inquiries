"""Axial charged-meson holdout for the common B neutrino/muon amplitudes."""
import argparse
from math import pi, sqrt
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def solve():
    p=read("inputs.json")
    q=read(p["inherited"]["rare"])
    m=read(p["inherited"]["masses"])
    hp=p["HPQCD_tableI_eq15"]
    template=q["BelleII"]
    bform=template["reference_SD_in_1e-6"]*1e-6
    nref=hp["lambda_ref"]*hp["Xt_ref"]/(pi*hp["alpha_inverse_ref"]*hp["sine_squared_ref"])
    hnu=bform/nref**2
    c=m["const_GeV"]
    hbar=read(p["inherited"]["kinematics"])["hbar_GeV_ps"]
    tH,dt=m["tau_ps"]["s_heavy"]
    f,*err=p["fBs_GeV"]
    norm=c["GF"]**2*c["m_mu"]**2*sqrt(c["M_Bs"]**2-4*c["m_mu"]**2)/(16*pi*hbar)
    run=m["LHCb_run3_2024_25"]["s_in_1e9"]
    cms=q["CMS"]
    obs=[("CMS",cms["branching_in_1e-9"],sum(cms["upper_widths"]),.75),
         ("LHCb_Run3",run[0],run[2]+run[4],p["run3_working_light_response"])]
    kmin=q["sm_coefficient_ratio"]["range"][0]
    result={}
    for label,center,width,w in obs:
        for j in (1,2):
            h_min=norm*(tH-j*dt)*(f-j*sum(err))**2*(1-p["mu_rate_transfer_width"])
            n_min=sqrt((template["mu"]-j*template["width"])*bform/(hnu*(1+j*hp["hadronic_rate_width"])))
            p_max=sqrt((center+j*width)*1e-9/(w*h_min))
            ratio=((1-kmin+2*p_max/n_min)/(1+kmin))**2
            result[label+str(j)]=ratio
    N=nref*sqrt(template["mu"])
    P=sqrt(cms["branching_in_1e-9"]*1e-9/(norm*tH*f*f))
    k=q["sm_coefficient_ratio"]["center"]
    signed=[]
    for s in [1,-1]:
        F=(k*N-s*P)/(1+k)
        L=N-F
        signed.append(abs((L-F)/N)**2)
        assert abs((k*L-F)**2-P**2)<1e-18
    result["central_real_sign_branches"]=signed
    return json.loads(json.dumps(result),parse_float=lambda x:float(f"{float(x):.11g}"))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=solve()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as fp:
            json.dump(out,fp,indent=2)
            fp.write("\n")
    else:
        assert out==read("independent_results.json")
        print("1464: separate absolute templates and longitudinal strange-vector holdout checked")
