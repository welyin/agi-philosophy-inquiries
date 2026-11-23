"""Eliminate the complex left amplitude in physical Z-coupling units."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def interval(nlo, nhi, mmax, k):
    return max(0.,(k*nlo-mmax)/(1+k)), (k*nhi+mmax)/(1+k)


def solve():
    p = read("inputs.json")
    a = read(p["inherited"]["rare"])
    b = read(p["inherited"]["masses"])
    r = p["HPQCD_tableI_eq15"]
    kaon = a["BelleII"]
    nref = r["lambda_ref"]*r["Xt_ref"]/(math.pi*r["alpha_inverse_ref"]*r["sine_squared_ref"])
    hnu = kaon["reference_SD_in_1e-6"]*1e-6/nref**2
    c = b["const_GeV"]
    hbar = read(p["inherited"]["kinematics"])["hbar_GeV_ps"]
    beta = math.sqrt(1 - 4*(c["m_mu"]/c["M_Bs"])**2)
    factor = c["GF"]**2*c["m_mu"]**2*c["M_Bs"]*beta/(16*math.pi*hbar)
    th,dth=b["tau_ps"]["s_heavy"]
    f,*df=p["fBs_GeV"]
    hmu = factor*f*f*th
    k0=a["sm_coefficient_ratio"]["center"]
    kmin,kmax=a["sm_coefficient_ratio"]["range"]
    cms=a["CMS"]
    intercept,slope=cms["lifetime_BF_factor"]
    weight=(cms["light_ps_example"]/cms["heavy_ps"]
            *(intercept+slope*cms["heavy_ps"])/(intercept+slope*cms["light_ps_example"]))
    rows={}
    run3=b["LHCb_run3_2024_25"]["s_in_1e9"]
    run3_observation={"branching_in_1e-9":run3[0],"upper_widths":[run3[2],run3[4]]}
    run3_weight=cms["light_ps_example"]/cms["heavy_ps"]/(1+b["LHCb_run3_2024_25"]["A_s_rate_rel_increase_at_0_minus1"][1])
    assert p["run3_working_light_response"] < run3_weight
    for key,obs,w in [("CMS",cms,cms["accepted_light_to_heavy_weight_range"][0]),
                      ("LHCb",a["LHCb"],1.),
                      ("LHCb_R3",run3_observation,p["run3_working_light_response"])]:
        for n in [1,2]:
            half=r["hadronic_rate_width"]*n
            low_rate=(kaon["mu"]-n*kaon["width"])*4.97e-6
            high_rate=(kaon["mu"]+n*kaon["width"])*4.97e-6
            nlo=math.sqrt(low_rate/(hnu*(1+half)))
            nhi=math.sqrt(high_rate/(hnu*(1-half)))
            hmlo=factor*(f-n*sum(df))**2*(th-n*dth)*(1-p["mu_rate_transfer_width"])
            mu_max=(obs["branching_in_1e-9"]+n*sum(obs["upper_widths"]))*1e-9
            mmax=math.sqrt(mu_max/(w*hmlo))
            bounds=[interval(nlo,nhi,mmax,k) for k in (kmin,kmax)]
            lower,upper=min(v[0] for v in bounds),max(v[1] for v in bounds)
            rows[f"{key}_{n}"]={"abs_N":[nlo,nhi],"abs_P_max":mmax,"abs_F":[lower,upper],
                                  "minimum_extra_coherent_if_F_zero":max(0.,kmin*nlo-mmax),
                                  "weak_diagonal_sum_min":2*lower,
                                  "acceptance_weight_above_which_zero_excluded":mu_max/(hmlo*(kmin*nlo)**2),
                                  "minimum_Belle_mu_to_exclude_zero_at_fixed_mu_upper":hnu*(1+half)*(mmax/kmin)**2/4.97e-6}
    central_n=nref*math.sqrt(kaon["mu"])
    central_p=math.sqrt(cms["branching_in_1e-9"]*1e-9/hmu)
    centralF=[(k0*central_n+s*central_p)/(1+k0) for s in (-1,1)]
    # Make the more permissive two-width example in the experimentally sampled, pure-light lifetime template.
    exn=math.sqrt((kaon["mu"]-2*kaon["width"])*4.97e-6/(hnu*(1+2*r["hadronic_rate_width"])))
    expmu=1e9*hmu*weight*(kmin*exn)**2
    assert cms["branching_in_1e-9"]-2*sum(cms["lower_widths"]) <= expmu <= cms["branching_in_1e-9"]+2*sum(cms["upper_widths"])
    for value,sgn in zip(centralF,(1,-1)):
        left=central_n-value
        assert abs(k0*left-value-sgn*central_p)<2e-18
    assert rows["CMS_1"]["abs_F"][0]>0 and rows["CMS_2"]["abs_F"][0]==0
    return {"Hnu":hnu,"Hmu_heavy":hmu,"nu_reference_amplitude":nref,
            "calculated_light_response":weight,"run3_central_light_response":run3_weight,"working_boxes":rows,
            "central_heavy_real_two_signs_abs_F":centralF,
            "CMS_two_width_zero_witness":{"abs_L":exn,"k":kmin,"F":0,"mu_Belle":2.0,
                    "relative_Hnu":1+2*r["hadronic_rate_width"],"light_template_CMS_e9":expmu}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    results=solve()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as f:
            json.dump(results,f,indent=2)
            f.write("\n")
    else:
        assert results==read("results.json")
        print("1464 author OK")
