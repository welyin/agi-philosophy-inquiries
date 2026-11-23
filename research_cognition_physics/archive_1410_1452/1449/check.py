"""Check signed HS branches, mass inequalities, a true aligned vacuum and the eta rate."""
import argparse
import itertools
import json
from pathlib import Path
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def neutral(x):
    u, d, s = x
    return np.array([[u+d, (u-d)/math.sqrt(3)],
                     [(u-d)/math.sqrt(3), (u+d+4*s)/3]])


def evaluate():
    data = json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    branches = json.loads((HERE.parent/"1445"/"results.json").read_text(encoding="utf-8"))["branches"]
    signs = []
    for branch in branches:
        middle, up, down = branch["signs_middle_down_and_light_up_over_charm_down_over_strange"]
        j = branch["response_mc_mt_abs_ub_over_cb_alpha_beta_gamma_J"][-1]
        signs.append({"middle_upRatio_downRatio":[middle,up,down], "theta_in_pi":int(up!=down), "J":j})
    assert sum(v["theta_in_pi"] for v in signs)==4 and all(v["J"] for v in signs)
    m = data["mass_GeV"]
    combos = []
    for corner in itertools.product(*(m[k] for k in ("pi0","pic","kc","k0"))):
        p,a,c,n = np.square(corner)
        combos.append([2*p-a+c-n, a+n-c, c+n-a])
    intervals = np.stack([np.min(combos, axis=0),np.max(combos,axis=0)],axis=1)
    assert np.all(intervals[:,0]>data["max_residual_I_GeV2"])
    x = np.array(data["example_X_GeV2"])
    ell = math.sqrt(x[0]*x[1]+x[1]*x[2]+x[2]*x[0])
    h = np.hypot(x,ell)
    angles = np.arctan2(ell,x)
    nmat = neutral(x)
    eig, rot = np.linalg.eigh(nmat)
    assert np.isclose(sum(angles),np.pi) and np.allclose(h*np.cos(angles),x)
    assert np.isclose(4*ell**2,3*np.linalg.det(nmat))
    cstar = (h[1]**2+h[2]**2-(h[1]*h[2]/h[0])**2)/(2*h[1]*h[2])
    assert -1<cstar<1 and np.isclose(cstar,math.cos(angles[0]))
    coarse = np.linspace(-1,1,2001)
    energy = h[0]*coarse+np.sqrt(h[1]**2+h[2]**2-2*h[1]*h[2]*coarse)
    assert sum(x) >= np.max(energy)-1e-12
    # Coefficient of p+ p- eta8 in Tr(Phi^3), directly by polarization.
    plus=np.eye(3,k=1); plus[1,2]=0
    minus=plus.T
    octet=np.diag([1.,1.,-2.])/math.sqrt(6)
    t3=np.diag([1.,-1.,0.])/math.sqrt(2)
    def cubic(a,b,c):
        return sum(np.trace(v[0]@v[1]@v[2]) for v in itertools.permutations((a,b,c)))
    assert np.isclose(cubic(plus,minus,octet),math.sqrt(6)) and np.isclose(cubic(plus,minus,t3),0)
    transfer=data["LO_transfer"]
    min_amp=transfer["mpi_min_GeV"]*transfer["meta_min_GeV"]*transfer["octet_projection_min"]/transfer["F_max_GeV"]
    beta=math.sqrt(1-(2*m["pic"][1]/m["eta"][0])**2)
    width=data["eta_gammagamma_width_upper_GeV"]/data["eta_gammagamma_BR_lower"]
    max_amp=math.sqrt(16*math.pi*m["eta"][1]*width*data["eta_to_pipi_BR_90"]/beta)
    cc=data["charged_current_check"]
    pion=m["pic"][0]; mu=cc["mu_mass_GeV"]
    fmax=math.sqrt(4*math.pi*cc["hbar_GeVs"]/(cc["charged_pion_lifetime_lower_s"]*
        cc["GF_GeVminus2"]**2*cc["Vud_lower"]**2*cc["net_current_amplitude_factor_lower"]**2*
        pion*mu**2*(1-mu**2/pion**2)**2))
    assert fmax<transfer["F_max_GeV"] and max_amp<min_amp*transfer["net_amplitude_retention_min"]
    return {"branch_table":signs,"Iu_Id_Is_box_GeV2":intervals.tolist(),
        "two_vacua_same_spectrum":{"X_theta0_GeV2":x.tolist(),"Bm_theta_pi_GeV2":h.tolist(),
        "phi_at_pi_rad":angles.tolist(),"L_abs_GeV2":ell,"global_max_cos_phi_u":cstar,
        "neutral_masses_GeV":np.sqrt(eig).tolist(),"abs_heavy_octet":abs(rot[1,1])},
        "conditional_F_max_from_total_charged_pion_width_GeV":fmax,
        "eta_width_upper_GeV":width, "eta_amplitude_upper_GeV":max_amp,
        "LO_eta_amplitude_floor_GeV":min_amp,"max_allowed_retained_amplitude":max_amp/min_amp,
        "predicted_BR_floor_if_retention_0p001":beta*(min_amp*transfer["net_amplitude_retention_min"])**2/(16*math.pi*m["eta"][1]*width)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=json.dumps(evaluate(),ensure_ascii=False,indent=2,allow_nan=False)+"\n"
    target=HERE/"results.json"
    if args.write:
        with target.open("x",encoding="utf-8") as stream: stream.write(out)
    elif target.read_text(encoding="utf-8")!=out:
        raise SystemExit("1449 result mismatch")
    else: print("1449 OK")
