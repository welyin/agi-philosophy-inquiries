"""Transport scalar amplitudes between two different physical parity projections."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def box(data, n=2):
    return [data[0]-n*sum(data[1:]), data[0]+n*sum(data[1:])]


def branch_windows(masses, angle, already_sine):
    out = {}
    for sign in (-1,1):
        values = []
        for d,s,b,t in itertools.product(*(box(row) for row in masses),box(angle)):
            sn = t if already_sine else math.sin(t)
            values.append(sn/math.sqrt(1-sn**2)*abs((b+sign*s)/(b-sign*d)))
        out[str(sign)] = [min(values),max(values)]
    return out


def verify_vertices():
    # Rank-one physical port, positive masses, arbitrary nontrivial phases.
    u = [math.sqrt(.946),math.sqrt(.053)*cmath.exp(.31j),math.sqrt(.001)*cmath.exp(-.47j)]
    for sign in (-1,1):
        p = [abs(x)**2 for x in u]
        ys,yb = 0.028,0.90
        yd = (p[1]*ys + sign*p[2]*yb)/p[0]
        assert yd>0
        eta = [-sign,sign,1]
        lam = [a*b for a,b in zip(eta,[yd,ys,yb])]
        assert abs(sum(a*b for a,b in zip(p,lam))) < 1e-16
        K = [(lam[i]+lam[2])*u[i].conjugate()*u[2] for i in (0,1)]
        for c in (1,1j,.61+.72j):
            gR = [k*c for k in K]
            gL = [k*eta[i]*c.conjugate() for i,k in enumerate(K)]
            for parity in (-1,1):
                assert abs((gR[1]+parity*gL[1])*K[0]-(gR[0]-parity*gL[0])*K[1]) < 1e-16


def integral_upper(inp):
    kin = inp["kinematics_GeV"]
    pars = inp["HPQCD_table_V"]
    mb,mk,mmu = [kin[k] for k in ("Bplus","Kplus","mu")]
    tplus = (pars["B_average_GeV"]+pars["K_average_GeV"])**2
    z = lambda t: (math.sqrt(tplus-t)-math.sqrt(tplus))/(math.sqrt(tplus-t)+math.sqrt(tplus))
    mass_min = inp["common_low_scale_mass_work_GeV"]["mb"][0]-inp["common_low_scale_mass_work_GeV"]["ms"][1]
    lo,hi = inp["q2_bin_GeV2"]
    nb = inp["integration_cells"]
    h = (hi-lo)/nb
    total=0.0
    total_mid=0.0
    max_ff=0.0
    for j in range(nb):
        left,right = lo+j*h, lo+(j+1)*h
        # z is negative and decreases; coefficient signs are fixed only for a0,a1.
        poly_max = box(pars["a0"])[1]+box(pars["a1"])[0]*z(left)+box(pars["a2"])[1]*z(right)**2
        poly_min = box(pars["a0"])[0]+box(pars["a1"])[1]*z(right)+box(pars["a2"])[0]*z(right)**2
        assert poly_min>0 and box(pars["a1"])[0]>0
        fhi = box(pars["L"])[1]*poly_max/(1-right/box(pars["pole_mass_GeV"])[0]**2)
        max_ff=max(max_ff,fhi)
        lam = lambda t: ((mb+mk)**2-t)*((mb-mk)**2-t)
        beta = lambda t: math.sqrt(1-4*mmu**2/t)
        common=(mb**2-mk**2)**2/(128*math.pi**3*mb**3)
        total+=h*common*math.sqrt(lam(left))*right*beta(right)**3*fhi**2/mass_min**2
        t=(left+right)/2
        f=pars["L"][0]*sum(pars["a"+str(i)][0]*z(t)**i for i in range(3))/(1-t/pars["pole_mass_GeV"][0]**2)
        total_mid+=h*common*math.sqrt(lam(t))*t*beta(t)**3*f**2/(4.18-.08)**2
    return total*(1+1e-12),total_mid,max_ff


def compute():
    inp=load("inputs.json")
    low,high=load(inp["flavor_MZ"]),load(inp["flavor_1e16"])["high_scale"]
    ratios={"MZ":branch_windows(low["yukawa"]["down"],low["quark_sines_and_phase"][0],True),
            "1e16":branch_windows(high["down"],high["angles_theta12_theta13_theta23_delta"][0],False)}
    verify_vertices()
    den,mid,fmax=integral_upper(inp)
    k=inp["kinematics_GeV"]
    mbd=k["Bd"]
    beta=math.sqrt(1-4*(k["mu"]/mbd)**2)
    mass=inp["common_low_scale_mass_work_GeV"]
    bd_floor=inp["work_fBd_floor_GeV"]**2*mbd**5*beta**3/(8*math.pi*(mass["mb"][1]+mass["md"][1])**2)
    beta_d=inp["CMS_Bd_rate_upper_95"]
    eff=inp["work_CMS_efficiency_floor"]
    transfer=inp["work_cross_amplitude_transfer_upper"]
    Gd=math.sqrt(beta_d/eff*inp["hbar_GeV_ps"]/inp["work_tau_Bd_floor_ps"]/bd_floor)
    out={}
    for key,ran in ratios["MZ"].items():
        Gs=Gd*transfer*ran[1]
        bcap=Gs**2*den*inp["work_tau_Bplus_ceiling_ps"]/inp["hbar_GeV_ps"]
        q0,q1=inp["q2_bin_GeV2"]
        denom=box(inp["LHCb_K_differential_per_GeV2"])[0]*(q1-q0)
        out[key]={"abs_Gplus_s_cap_GeV_minus2":Gs,"scalar_K_bin_branching_cap":bcap,
                  "scalar_fraction_of_observed_bin_cap":bcap/denom}
    return {"ratio_windows_using_marginals_before_return":ratios,
            "Bd_width_per_abs_Gminus_d_squared_floor_GeV5":bd_floor,
            "abs_Gminus_d_cap_GeV_minus2":Gd,
            "K_bin_width_per_abs_Gplus_s_squared_upper_GeV5":den,
            "K_bin_central_check_GeV5":mid,"f0_grid_upper":fmax,
            "MZ_cross_caps":out,
            "coefficient_convention":"Gplus,Gminus=(C_R+/-C_L)/2; double each bound if quoting full sum/difference"}


def rounded(x):
    if isinstance(x,dict): return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list): return list(map(rounded,x))
    if isinstance(x,float): return float(format(x,".12g"))
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--write-new",action="store_true")
    g.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=rounded(compute())
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    elif result!=load("results.json"):
        raise SystemExit("Saved parity-transfer certificate differs")
    print(json.dumps(result))
