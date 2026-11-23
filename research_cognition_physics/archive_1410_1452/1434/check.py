"""Single universal left-handed leptoquark: closure and experimental envelopes."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def mixing(s12,s13,s23,delta):
    c12,c13,c23 = [math.sqrt(1-x*x) for x in (s12,s13,s23)]
    phase = np.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,s13/phase],
                     [-s12*c23-c12*s23*s13*phase,c12*c23-s12*s23*s13*phase,s23*c13],
                     [s12*s23-c12*c23*s13*phase,-c12*s23-s12*c23*s13*phase,c23*c13]])


def normalized_rates(v,w,u,coupling):
    # Rows j = u,c,t and columns ell = e,mu,tau. The outgoing neutrino is unresolved.
    first = v[:,2,None,None]*u[None,:,:]
    second = (v@w@u)[:,None,:]*w[2,:].conj()[None,:,None]
    result = np.sum(abs(first+coupling*second)**2,axis=2)
    np.testing.assert_allclose(result.sum(axis=0)-1,(2*coupling+3*coupling**2)*abs(w[2])**2,atol=3e-14)
    np.testing.assert_allclose(result.sum(axis=1),(3+2*coupling)*abs(v[:,2])**2+coupling**2,atol=3e-14)
    return result


def high_q_bounds(inputs,number):
    threshold = inputs["q2_cut_GeV2"]
    masses = inputs["masses_GeV"]
    helicity = [(1-m*m/threshold)**2*(1+m*m/(2*threshold)) for m in masses]
    lifetime,stat,systematic = inputs["B0_lifetime_ps"]
    lifetime -= number*math.hypot(stat,systematic)
    norm,error = inputs["HPQCD_partial_massless_per_ps"]
    norm -= number*error
    def bin_sum(key):
        return 1e-4*sum(rate+number*math.hypot(stat,syst) for rate,stat,syst in inputs[key])
    upper_rates = [bin_sum("branch_e_high_bins_1e4"),bin_sum("branch_mu_high_bins_1e4"),
                   inputs["branch_tau_full_upper_90"]]
    lower_response=inputs["min_relative_rate_response"]
    strength = [upper/(lower_response*life_factor*lifetime*norm)
                for upper,life_factor in zip(upper_rates,helicity)]
    strength_total=sum(strength)
    upper=math.sqrt(strength_total)
    return {"widths":number, "unit_helicity_floors":helicity, "measured_branch_caps":upper_rates,
            "normalized_strength_caps":strength, "summed_upper":strength_total,
            "lambda_upper":upper,"M_over_g_lower_GeV":inputs["v_GeV"]/math.sqrt(2*upper)}


def charged_tau_control(inputs, light_caps):
    """Substitute the measured purely leptonic tau channel for the neutral-B tau limit."""
    cfg = inputs["charged_tau_control"]
    gf = read(HERE/cfg["fermi_source"])["GF"]
    hbar = read(HERE/cfg["hbar_source"])["hbar_GeV_ps"]
    def edge(row, sign):
        return row[0]+sign*2*math.hypot(*row[1:])
    branch_upper = edge(cfg["branch_tau_1e4"],1)*1e-4
    decay_constant = edge(cfg["f_Bplus_GeV"],-1)
    lifetime = edge(cfg["tau_Bplus_ps"],-1)
    mass = edge(cfg["mass_Bplus_GeV"],-1)
    mt = min(cfg["mass_tau_GeV_window"])
    assert max(cfg["mass_tau_GeV_window"]) < mass/math.sqrt(3)
    unit = (gf**2/(8*math.pi)*decay_constant**2*lifetime/hbar*mass*mt**2
            *(1-(mt/mass)**2)**2)
    tau_upper = branch_upper/(inputs["min_relative_rate_response"]*unit)
    total = sum(light_caps[:2])+tau_upper
    coupling = math.sqrt(total)
    return {"charged_tau_branch_upper":branch_upper,"branch_per_unit_strength_floor":unit,
            "charged_tau_strength_upper":tau_upper,"three_strengths_upper":total,
            "lambda_upper":coupling,"M_over_g_floor_GeV":inputs["v_GeV"]/math.sqrt(2*coupling)}


def run():
    inputs = read(HERE/"inputs.json")
    tag = read(HERE/inputs["down_b_tag"])
    calibrated = read(HERE/inputs["ckm_example"])
    v = mixing(*[x[0] for x in calibrated["quark_sines_and_phase"]])
    u = mixing(.55,.149,.72,1.24)  # Any unitary works, this one contains nontrivial phases.
    p = tag["central_b_mu_sharp_floor"]
    c,s=math.sqrt(1-p),math.sqrt(p)
    w = np.array([[1,0,0],[0,c,-s],[0,s,c]],dtype=complex)
    lam=inputs["phase_example_lambda"]
    examples=[]
    z=float(sum(abs(v[:2,2])**2))
    for sign in (1,-1):
        # Different relative s-quark left orientations at the same fixed CKM and mass data.
        w_direction=np.diag([1,sign,1])@w
        rates=normalized_rates(v,w_direction,u,lam)
        in_light=float(sum(rates[:2,1]))
        examples.append({"relative_s_sign":sign,"open_mu_sum":in_light,"open_mu_shift":in_light-z,
                         "full_mu_shift":float(sum(rates[:,1])-1),
                         "sum_b_to_u_three_charged_leptons":float(sum(rates[0]))})
    assert examples[0]["open_mu_shift"]*examples[1]["open_mu_shift"] < 0
    full_shift=(2*lam+3*lam**2)*tag["sharp_floor_range_two_widths"][0]
    rng=np.random.default_rng(1434)
    for _ in range(30):
        def random_unitary():
            return np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
        vv,ww,uu = random_unitary(), random_unitary(), random_unitary()
        g=float(rng.uniform(.001,.3))
        rates=normalized_rates(vv,ww,uu,g)
        assert np.all(rates[:2].sum(axis=0)+1e-14 >= g*g*abs(ww[2])**2)
        np.testing.assert_allclose(rates,normalized_rates(vv,ww,np.eye(3),g),atol=3e-14)
    bounds=[high_q_bounds(inputs,n) for n in (0,1,2)]
    assert bounds[-1]["lambda_upper"] < .019
    charged=charged_tau_control(inputs,bounds[-1]["normalized_strength_caps"])
    assert charged["lambda_upper"] < bounds[-1]["lambda_upper"]
    return {"experimental_work_boxes":bounds,"charged_B_tau_replacement":charged,
            "SM_open_mu_baseline_for_example":z,
            "minimum_full_mu_shift_for_lambda_0p005":full_shift,
            "two_equal_spectra_and_moduli_realizations":examples}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=run()
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8") as stream:
            json.dump(data,stream,indent=2)
            stream.write("\n")
    elif args.check:
        def close(x,y):
            if isinstance(x,dict):
                return x.keys()==y.keys() and all(close(x[k],y[k]) for k in x)
            if isinstance(x,list):
                return len(x)==len(y) and all(close(a,b) for a,b in zip(x,y))
            return math.isclose(x,y,rel_tol=3e-11,abs_tol=3e-15)
        assert close(data,read(target))
        print("1434 closure, truncated sign counterexample and data bound verified")
    else:
        print(json.dumps(data,indent=2))
