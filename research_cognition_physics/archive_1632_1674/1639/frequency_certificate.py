"""Phase-uniform analytic lower bound on the Bs clock from the two absolute rate disks."""
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def calculate():
    states=read("../1635/results.json")
    inputs=read("../1636/inputs.json")
    clock=read("../1636/results.json")
    rare=read("../1464/results.json")["working_boxes"]
    cal=read("../1464/inputs.json")["HPQCD_tableI_eq15"]
    k=read("../1406/inputs.json")["sm_coefficient_ratio"]["center"]
    weak=read("../1404/inputs.json")
    lattice=read("../1407/inputs.json")
    hat,dh=lattice["HPQCD_inputs"]["fBs_sqrtB_GeV"]
    scale=weak["const_GeV"]
    gf,mw=scale["GF"],scale["MW"]
    norm=gf**2*mw**2/(4*math.pi**2)
    pref=2*scale["M_Bs"]/(3*lattice["hbar_GeV_ps"])
    c=cal["Xt_ref"]/(math.pi*cal["alpha_inverse_ref"]*cal["sine_squared_ref"])
    bag_data=inputs["LR_F2B4_B5_GeV2"]
    bag_weights=inputs["LR_B4_B5_weights"]
    s0,ds0=weak["short_distance"]["S0_published"]
    eta,de=weak["short_distance"]["etaB"]
    xt=weak["short_distance"]["xt"][0]
    h1=-(xt-7)/(4*(xt-1))-(xt**2-2*xt+4)*math.log(xt)/(2*(xt-1)**2)
    emin,emax=read("../1468/inputs.json")["eta_Z_RGI_range"]
    result={}
    for key in ("common_witness","high_mass_witness"):
        ref=states[key]
        M=ref["M_Y_GeV"]
        ms,mb=ref["physical_masses_GeV"][1:3]
        r=math.sqrt(ref["Tbb"])
        a=ms*mb/M**2
        light=(1-(ms/M)**2)*(1-(mb/M)**2)
        dist=(ms/M)**2+ms*mb*r/M**2
        g0_b=1-abs(complex(*ref["left_S"][2]))**2-(mb*r/M)**2
        assert g0_b-dist>.98
        # Compare any allowed F to the frozen one using the square root of the light Gram matrix.
        lt_radius=2*dist/math.sqrt(.98)
        lt=abs(complex(*clock[key]["lambda_t"]))
        assert lt_radius<2.e-7
        lo,hi=lt-2.e-7,lt+2.e-7
        h,j=1+2*a,1-2*a
        hard=clock[key]["Qd_hard_to_tree_RR_after_RG_abs"]
        v=gf*c**2/(math.sqrt(2)*norm)*(1+4*a*a+hard)
        by_width={}
        for n in (1,2):
            box=rare[f"LHCb_R3_{n}"]
            Nlo=box["abs_N"][0]/(c*hi)
            Phi=box["abs_P_max"]/(c*lo)
            x=(Nlo**2-1-(h/j)**2*(Phi**2-k**2))/(2*h*(1+h*k/j))
            assert x>0
            y=x*x-(Phi/j)**2
            q=(hat+n*dh)**2*emax if y<0 else (hat-n*dh)**2*emin
            lr=sum(w*(b-n*db) for w,(b,db) in zip(bag_weights,bag_data))
            extent=((hat-n*dh)**2*(s0-n*ds0)*(eta-n*de)
                    +q*v*y+c*xt*(math.log(M/mw)+h1)*light*x*lr)
            assert extent>0
            freq=pref*lo**2*norm*extent
            by_width[str(n)]={"normalized_real_F_lower":float(f"{x:.10g}"),
                              "normalized_real_F_squared_lower":float(f"{y:.10g}"),
                              "clock_lower_ps_inverse":float(f"{freq:.10g}")}
        assert by_width["1"]["clock_lower_ps_inverse"]>17.78
        result[key]=by_width
    return result


if __name__=="__main__":
    print(json.dumps(calculate(),indent=2))
