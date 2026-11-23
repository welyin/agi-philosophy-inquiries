"""Reconstruct gauge poles and dimension-six lifetimes on the triplet/sextet curve."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent
TWOPI = 2*math.pi
SM = np.array([4.1, -19/6, -7])
MSSM = np.array([6.6, 1, -3])
T = np.array([3.6, 4, 0])
S = np.array([3.2, 0, 5])
Z = 7-T-S
H = np.array([.4, 0, 1])
ADJ = np.array([0, 2, 3])
GAMMA = np.array([[11/20, 9/4, 2], [11/30, 3/2, 4/3]])


def read(name):
    return json.loads((HERE/name).read_bytes())


def low_boundary(data, calibration):
    ew = read(data["electroweak_boundary"])
    low = 4*math.pi / np.array([ew[k][0] for k in ("g1_GUT", "g2", "g3")])**2
    if calibration == "ALPHA":
        cf = read(data["lattice_coupling"])
        low[2] = 1/cf["strong_lattice_5f"]["ALPHA_2026"][0]
        low[2] += math.log(cf["top_MS_matching_GeV"][0]/ew["mu_ref_GeV"])/(3*math.pi)
    return ew["mu_ref_GeV"], low


def baseline_events(base):
    ms = base["soft_GeV"]
    events = [(ms, MSSM-SM, -ADJ/(12*math.pi))]
    for label, charges in (("colored", H), ("weak", 1-H)):
        for spin, weight in (("fermions", 2/3), ("scalars", 1/3)):
            for m in base[f"{label}_{spin}_GeV"]:
                events.append((m, weight*charges, np.zeros(3)))
    return events


def integrate(start, a, soft, events, stop):
    """Solve inverse couplings and the scalar duql Wilson RGE at each physical pole."""
    a = a.copy()
    b = SM.copy()
    exposure = np.zeros(3)
    extra_hypercharge = 0.
    min_inv = min(a)
    for m, step, offset in sorted([*events, (stop,np.zeros(3),np.zeros(3))],
                                  key=lambda event:event[0]):
        if m > stop:
            break
        dt = math.log(m/start)/TWOPI
        slope = b != 0
        area = dt/a
        area[slope] = -np.log1p(-b[slope]*dt/a[slope])/b[slope]
        exposure += GAMMA[int(start>=soft)]*area
        extra_hypercharge += (.4 if start>=soft else .6)*area[0]
        a = a-b*dt+offset
        min_inv = min(min_inv, min(a))
        b += step
        start = m
    return a, math.exp(sum(exposure)), math.exp(sum(exposure)+extra_hypercharge), min_inv


def split(base, low_shift, mt):
    u0 = base["U_GeV"]
    length = math.log(u0/mt)
    h,s = np.linalg.solve([[5.2,3.2],[8,-5]],
                          [.4*length, -4*length-TWOPI*low_shift])
    u = u0*math.exp(h)
    return dict(U=u, T=mt, S=u*math.exp(-s), Z=u)


def physical(data, base, calibration, poles, radius):
    start, a = low_boundary(data, calibration)
    u = poles["U"]
    soft = base["soft_GeV"]
    events = baseline_events(base)
    for m, charge in ((poles["T"],T),(poles["S"],S),(poles["Z"],Z),
                      (u,H),(u/radius**2,ADJ)):
        events.append((m, charge, np.zeros(3)))
    end, gain, gain_10, min_inv = integrate(start, a, soft, events, u*radius)
    assert np.ptp(end) < 1e-6 and min_inv > 5
    a5 = np.mean(end)
    A = gain*data["long_distance_MZ_to_2GeV"]
    g = math.sqrt(4*math.pi/a5)
    return {"inverse_alpha_X":a5, "duql_A":A,
            "qque_A":gain_10*data["long_distance_MZ_to_2GeV"], "gX":g,
            "contact_scale_GeV":u*radius/g/math.sqrt(A),
            "min_inverse_coupling":min_inv}


def scale_limits(data, row, error):
    p = read(data["pion_inputs"])
    k = read(data["kaon_inputs"])
    quark_d,quark_s = row[:2]
    mp,mn,pk,nk = (k["mass_GeV"][key] for key in ["p","n","K+","K0"])
    def lower(arr):
        return max(0, abs(arr[0])-error*sum(arr[1:]))
    w = lower(p["pion_W0_GeV2"])
    wa,wb,wc = [lower(k["W0_GeV2"][key]) for key in ["us_d","ud_s","ds_u"]]
    channels = [
        ("p_pi_plus",mp,p["pion_mass_GeV"],(w*quark_d)**2,p["pion_lifetime_years"]),
        ("n_pi_zero",mn,data["neutron_neutral_pion_mass_GeV"],(w*quark_d)**2/2,
         data["neutron_neutral_pion_limit_years"]),
        ("p_K_plus",mp,pk,(wa*quark_d)**2+(wb*quark_s)**2,k["proton"]["partial_lifetime_lower_yr"]),
        ("n_K_zero",mn,nk,(wc*quark_d)**2+(wb*quark_s)**2,k["neutron"]["partial_lifetime_lower_yr"])]
    return {key:(mass/(32*math.pi)*(1-(meson/mass)**2)**2*amplitude_sq*lifetime*
                 p["year_seconds"]/p["hbar_GeV_s"])**.25
            for key,mass,meson,amplitude_sq,lifetime in channels}


def required(data,base,label,poles,limits):
    f=lambda r:physical(data,base,label,poles,r)
    contact=max(limits.values())
    if f(1)["contact_scale_GeV"]>=contact:
        return 1.
    a,b=data["MX_over_U_window"]
    assert f(b)["contact_scale_GeV"]>contact
    for _ in range(47):
        mid=(a+b)/2
        if f(mid)["contact_scale_GeV"]<contact:
            a=mid
        else:
            b=mid
    return b


def summary(data,base,label,poles,limits):
    p=physical(data,base,label,poles,1.)
    return dict(poles_GeV=poles,r1=p,
                r1_rate_over_search_by_channel={k:(v/p["contact_scale_GeV"])**4 for k,v in limits[0].items()},
                required_r_in_declared_window=[required(data,base,label,poles,l) for l in limits],
                possible_amplitude_cancellation_needed_r1=max(0,1-(p["contact_scale_GeV"]/max(limits[0].values()))**2))


def marked_electron(data,base,fermi,shape,labels,comparisons,endpoints):
    api=runpy.run_path(str(HERE/data["electron_mode_basis_and_kinematics"]))
    p=read(data["pion_inputs"])
    ew=read("../1725/inputs.json")
    factor=ew["electron_mode_partial_lifetime_lower_years"]*p["year_seconds"]/p["hbar_GeV_s"]
    mp,mpi,m=api["M_P"],api["M_PI"],api["LEPTONS"][0]
    e=(mp*mp-mpi*mpi+m*m)/(2*mp)
    momentum=math.sqrt(e*e-m*m)
    w0,w1=map(abs,api["FORM"][0])
    w0err,w1err=sum(api["ERR0"][0]),sum(api["ERR1"][0])
    rows={}
    for key,x in labels.items():
        *_,le,_,d=fermi["matrices"](np.array(x),shape)
        a=np.vdot(le[:,0],d[:,0])
        b=d[0,0]*le[0,0].conj()
        cases={}
        points=comparisons | {f"{cal}_{branch}_upper":(cal,row[key][branch]["poles_GeV"])
                              for cal,row in endpoints.items() for branch in row[key]}
        for name,(label,poles) in points.items():
            def ratios(r,error):
                phys=physical(data,base,label,poles,r)
                norm=(phys["gX"]/(r*poles["U"]))**4*factor
                norm*=abs(phys["duql_A"]*a)**2+abs(phys["qque_A"]*(a+b))**2
                # Conservative two-sided response for every relative chiral phase.
                common=e*momentum/(16*math.pi*mp)
                return [norm*common*(1-m/e)*(w0-error*w0err-m/mp*(w1+error*w1err))**2,
                        norm*common*(1+m/e)*(w0+error*w0err+m/mp*(w1+error*w1err))**2]
            def r_needed(error):
                if ratios(1,error)[0]<=1:
                    return 1.
                lo,hi=data["MX_over_U_window"]
                assert ratios(hi,error)[0]<1
                for _ in range(40):
                    mid=(lo+hi)/2
                    if ratios(mid,error)[0]>1:
                        lo=mid
                    else:
                        hi=mid
                return hi
            cases[name]=dict(r1_central_rate_sandwich=ratios(1,0),
                             necessary_r_central_and_conservative_one_width=[r_needed(n) for n in (0,1)],
                             comparison_r=(6.5 if "split" in name or "upper" in name else 2.0),
                             comparison_rate_sandwich=ratios(6.5 if "split" in name or "upper" in name else 2.0,0))
        rows[key]=dict(electron_duql_and_qque_modulus=[abs(a),abs(a+b)],cases=cases)
    return rows

def main():
    data=read("inputs.json")
    base=read(data["ordinary_split_spectrum"])["shared_witness"]
    fits=read(data["triplet_limits"])["five_cases"]
    fermi=runpy.run_path(str(HERE/data["symmetric_weak_frame"]))
    shape=fermi["backend"]()
    row=abs(shape[3][0])
    limits=[scale_limits(data,row,n) for n in [0,1,2]]
    mtmin=data["triplet_window_lower_GeV"]
    thresholds={key:{label:case[label]["triplet_mass_at_product_GeV"]["cap_product"]
                      for label in ["old_absolute_matching","restored_common_running"]} for key,case in fits.items()}
    low0=low_boundary(data,"original")[1]
    results={}
    global_inverse_lower=[]
    for label in ["original","ALPHA"]:
        shift=low_boundary(data,label)[1][2]-low0[2]
        at=lambda t:split(base,shift,t)
        mtmax=max(max(v.values()) for v in thresholds.values())
        # Concavity in log(mu),log(MT),log(r) bounds all points by the
        # endpoint vertices; extending the finite scheme step to MZ lowers a.
        lower=min(low_boundary(data,label)[1]-ADJ/(12*math.pi))
        for m in (mtmin,mtmax):
            for r in data["MX_over_U_window"]:
                lower=min(lower,physical(data,base,label,at(m),r)["inverse_alpha_X"]-1e-6)
        global_inverse_lower.append(lower)
        endpoints={}
        for key,table in thresholds.items():
            endpoints[key]={q:summary(data,base,label,at(t),limits) for q,t in table.items()}
        cross={}
        for n,l in enumerate(limits):
            a,b=mtmin,max(max(v.values()) for v in thresholds.values())
            measure=lambda t:physical(data,base,label,at(t),1)["contact_scale_GeV"]
            if measure(a)>=max(l.values()):
                cross[str(n)]=dict(all_pass=True)
            elif measure(b)<max(l.values()):
                cross[str(n)]=dict(all_fail=True)
            else:
                for _ in range(47):
                    m=math.sqrt(a*b)
                    if measure(m)<max(l.values()):
                        a=m
                    else:
                        b=m
                cross[str(n)]=dict(min_triplet_mass_for_r1_GeV=b)
        results[label]=dict(at_lower_end=summary(data,base,label,at(mtmin),limits),
                            at_case_upper_ends=endpoints,r1_threshold_vs_MT=cross)
    u0=base["U_GeV"]
    complete=summary(data,base,"original",dict(U=u0,T=mtmin,S=mtmin,Z=mtmin),limits)
    a=read(data["adjoint_comparison"])["threshold_cases"]["1698_two_pairs"]["experimental_working_boxes"]["0"]
    t,s,z=a["lowest_scale_positive_mass_witness_GeV_T_S_Z"]
    adj=summary(data,base,"ALPHA",dict(U=a["witness_U_GeV_and_inverse_alpha_U"][0],T=t,S=s,Z=z),limits)
    comparisons={"split_at_1e13":("ALPHA",results["ALPHA"]["at_lower_end"]["poles_GeV"]),
                           "degenerate_15":("original",complete["poles_GeV"]),
                           "adjoint_close":("ALPHA",adj["poles_GeV"])}
    labels=read("../1799/inputs.json")
    marked=marked_electron(data,base,fermi,shape,
                           labels["mixed_type_I_II"]|labels["fixed_delta_minus70_mass_readout"],comparisons,
                           {k:v["at_case_upper_ends"] for k,v in results.items()})
    finite=read(data["finite_gram"])["cases"]
    amin=min(global_inverse_lower)
    rmax=data["MX_over_U_window"][1]
    gamma_max=GAMMA[1]+[.4,0,0]
    dlogA=sum(gamma_max)/(TWOPI*amin)+2*float(gamma_max@ADJ)*3*math.log(rmax)/(TWOPI*amin)**2
    # Derivatives bound ln(g^2 A/M_X^2); a negative sign guarantees one crossing.
    dr=16/(TWOPI*amin)+dlogA-2
    dmt=(9/43)*(sum(gamma_max)/(TWOPI*amin)-2)
    assert dr<0 and dmt<0
    return dict(weak_first_row_moduli=row.tolist(),monotonicity_certificate={
                "whole_curve_inverse_alpha_lower":amin, "dlog_amplitude_dlog_r_upper":dr,
                "dlog_amplitude_dlog_MT_upper":dmt},
                max_finite_neutrino_missing_gram=
                3*max(c["sterile_norm_upper_per_light_vector"]**2 for c in finite.values()),
                minimum_contact_scale_GeV_by_error=limits,split=results,
                degenerate_15_original=complete,adjoint_split_ALPHACentral=adj,
                marked_single_X_electron=marked)


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    answer=tidy(main())
    blob=json.dumps(answer,indent=2,allow_nan=False)+"\n"
    target=HERE/"results.json"
    if args.check:
        assert json.loads(target.read_bytes())==answer
    else:
        target.write_bytes(blob.encode())
    print("1804",hashlib.sha256(blob.encode()).hexdigest()[:12])
