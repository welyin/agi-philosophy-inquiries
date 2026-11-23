"""Two neutral-B clocks and two phases versus an independently labeled matter texture."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())
BASE = runpy.run_path(str(HERE / CFG["charged_family"]))
CLOCKS = BASE["CLOCKS"]
CS = BASE["CFG"]["mixing_coefficient_ps_inv"][0]
S = CLOCKS["golden_S"][0]
PHASE = BASE["read"](BASE["CFG"]["secondary_tree_angle"])["time_integrated_gamma_deg"]
QCD = BASE["STRONG"]
TAU = runpy.run_path(str(HERE / CFG["leptonic_B"]))
BF_WINDOW = TAU["band"](TAU["CFG"]["BR_Belle_combined"],2)
BF_COEFF = TAU["factors"](2)["br_tau_unit"]


def unitary_roots(gamma_deg, xi, beta):
    gamma = math.radians(gamma_deg)
    ps = math.sqrt(CLOCKS["delta_ms"][0] / CS)
    ratio = xi * math.sqrt(CLOCKS["delta_md"][0]/CLOCKS["delta_ms"][0]
                           * CLOCKS["mass_s"][0]/CLOCKS["mass_d"][0])
    pd = ps*ratio
    a = pd*math.sin(beta)/math.sin(gamma)
    b = pd*math.sin(gamma+beta)/math.sin(gamma)
    v = pd**2+ps**2
    small = 2*v/(1+math.sqrt(1-4*v))
    out = []
    for t in [small, 1-small]:
        total = 1-t
        k = 1-pd**2/t
        trace = k*total+a*a-b*b
        disc = trace**2-4*k*a*a*total
        if disc < 0:
            continue
        x1 = 2*a*a*total/(trace+math.sqrt(disc))
        for x in [x1,(trace+math.sqrt(disc))/(2*k)]:
            if not 0 < x < total or k <= 0:
                continue
            bottom = np.sqrt([x,total-x,t])
            products = np.array([a*np.exp(1j*gamma), -b, pd*np.exp(-1j*beta)])
            down = products / bottom
            strange = np.cross(down.conj(),bottom)
            ckm = np.column_stack([down,strange,bottom])
            np.testing.assert_allclose(ckm.conj().T@ckm,np.eye(3),atol=2e-10)
            check_gamma = np.angle(-products[0]/products[1],deg=True)
            check_beta = np.angle(-products[1]/products[2],deg=True)
            assert abs(check_gamma-gamma_deg) < 1e-10 and abs(check_beta-math.degrees(beta))<1e-10
            assert abs(abs(ckm[2,1]*ckm[2,2]) - ps) < 1e-10
            br_band=BF_COEFF*x/pd**2
            passes = (bottom[2] >= CFG["top_decay"]["Vtb_lower_95"] and
                      br_band[0]<=BF_WINDOW[1] and BF_WINDOW[0]<=br_band[1])
            out.append({"Vtb":float(bottom[2]),"Vus":abs(strange[0]),"Vub":float(bottom[0]),
                        "Vcb":float(bottom[1]),"beta":check_beta,"gamma":check_gamma})
            out[-1].update(Btau_two_widths=br_band.tolist(),top_and_Btau_compatible=int(passes))
    return out


def fit(gamma, start, step=.25):
    def prediction(x):
        row = BASE["evaluate"](x,step)
        return row,np.log([row["gamma"]/gamma,row["sin2beta"]/S,
                           row["tbs_squared"]*CS/CLOCKS["delta_ms"][0]])
    x=np.log(start)
    row,res=prediction(x)
    h=1.e-4
    jac=np.column_stack([(prediction(x+d)[1]-res)/h for d in np.eye(3)*h])
    for _ in range(16):
        if max(abs(res)) < 2e-10:
            break
        dx=np.linalg.solve(jac,-res)
        x+=dx
        row,new=prediction(x)
        jac+=np.outer(new-res-jac@dx,dx)/(dx@dx)
        res=new
    assert max(abs(res)) < 2e-10
    return row,x


def compare(row):
    sd=1/row["ds"]
    old=BASE["compare"](row)
    xi=row["td_ts"] / math.sqrt(CLOCKS["delta_md"][0]/CLOCKS["delta_ms"][0]
                               *CLOCKS["mass_s"][0]/CLOCKS["mass_d"][0])
    return {**row,"sd":sd,"QCDSF_fixed_point_pull":(sd-QCD["QCDSF_s_d"][0])/QCD["QCDSF_s_d"][1],
            "xi_needed_to_match_Bd":xi,"xi_HP_fixed_point_pull":(xi-CLOCKS["HPQCD"][0])/CLOCKS["HPQCD"][1],
            "HPQCD_delta_md_pred":CLOCKS["delta_md"][0]*(xi/CLOCKS["HPQCD"][0])**2,
            "heldouts":old}


def competing_branch_bound():
    # Outward box for the newer gamma and both published xi choices.
    def extrema(d):
        return [d[0]-math.hypot(*d[1:]),d[0]+math.hypot(*d[1:])]
    md,ms=extrema(CLOCKS["delta_md"]),extrema(CLOCKS["delta_ms"])
    mbd,mbs=extrema(CLOCKS["mass_d"]),extrema(CLOCKS["mass_s"])
    q=max(extrema(CLOCKS[l])[1] for l in ("HPQCD","FNAL_MILC"))*math.sqrt(md[1]/ms[0]*mbs[1]/mbd[0])
    ps2=ms[1]/(CS-BASE["CFG"]["mixing_coefficient_ps_inv"][1])
    h=(1+q*q)*ps2
    small_top=2*h/(1+math.sqrt(1-4*h))
    kmin=1-q*q*ps2/(1-small_top)
    gam_min=math.radians(PHASE["68"][0])
    assert PHASE["68"][1]<90
    midpoint=4*q*q/((1+q*q)*math.sin(gam_min)**2)
    beta_min=(math.pi-math.asin(S+CLOCKS["golden_S"][1]))/2
    reflected=BF_COEFF[0]*math.sin(beta_min)**2
    large_ub=BF_COEFF[0]*(1+1/q**2)/2
    assert small_top < CFG["top_decay"]["Vtb_lower_95"]**2
    assert midpoint < kmin and min(reflected,large_ub)>BF_WINDOW[1]
    return dict(Pd_over_Ps_upper=q,small_top_sq_upper=small_top,k_lower=kmin,
                midpoint_left_upper=midpoint,reflected_beta_Btau_lower=reflected,
                large_ub_Btau_lower=large_ub)


def run():
    beta=math.asin(S)/2
    roots={}
    for g in [PHASE["central"],63.6]:
        for label in ["HPQCD","FNAL_MILC"]:
            surviving=0
            for phase,b in [("small",beta),("large",math.pi/2-beta)]:
                rows=unitary_roots(g,CLOCKS[label][0],b)
                roots[f"{g}_{label}_{phase}"]=rows
                surviving+=sum(row["top_and_Btau_compatible"] for row in rows)
            assert surviving==1
    central,x=fit(PHASE["central"],[.22426,.003839,.04251])
    alt,_=fit(63.6,np.exp(x))
    fine=BASE["evaluate"](x,.125)
    fine_diff=[fine[key]-central[key] for key in ("gamma","sin2beta","ds","bs","tbs_squared")]
    assert max(abs(v) for v in fine_diff)<2e-6
    return {"Belle_tau_comparison_two_widths":BF_WINDOW.tolist(),"branch_exclusion_box":competing_branch_bound(),
            "unitary_competitors":roots,"charged_family":{"updated":compare(central),
            "2024":compare(alt)},"fine_absolute_differences":fine_diff}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=BASE["rounded"](run())
    file=HERE/"results.json"
    if args.check:
        if result != json.loads(file.read_bytes()):
            raise SystemExit("saved 1766 result differs")
        print("1766 two-clock reconstruction and material comparison verified")
    else:
        file.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
