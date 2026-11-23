"""Bound the inverted Dirac phase by eliminating the unknown Majorana phase."""
import argparse
import cmath
import json
import math
from pathlib import Path


ROOT=Path(__file__).resolve().parent
ARCH=ROOT.parent


def load(rel):
    return json.loads((ARCH/rel).read_text(encoding="utf-8-sig"))


def dipole_cap(m,j,t,br,eta):
    norm=.0072973525205055605*j**3*t*(1-(m/j)**2)**3/(4*6.582119569e-25)
    return math.sqrt(br/(eta*norm))


def evaluate(scenario,ranges):
    s0,s1,a0,a1,t0,t1,r0,r1,m20,m21=ranges  # s,a,t are squared sines
    me,mm,mt=.00051099895,.1056583755,1.775
    eta=scenario["rate_floor"]*(1-scenario["cancel"])**2
    em=dipole_cap(me,mm,2.1969803e-6,1.5e-13,eta)
    et=dipole_cap(me,mt,288.e-15,3.3e-8,eta)
    mu_tau=dipole_cap(mm,mt,288.e-15,4.4e-8,eta)
    # N=U_2 diag(m1,m2) U_2^T. The other entries of each flavor row
    # have squared length 1-|U_i3|^2, even for complex phases.
    ee=1-a0
    mu=1-(1-a1)*t0
    mu_tau_entry=(1+a1)/2
    e_tau_entry=math.sqrt((1-a0)*(1-(1-a1)*(1-t1)))
    gap=((1.176e-13-scenario["other_ae"])/(2*me**2)
         -(1.64e-9+scenario["other_amu"])/(2*mm**2))
    h=(em*math.hypot(ee/me,mu/mm)+(et*mu_tau_entry+mu_tau*e_tau_entry)/mt
       +scenario["defect_per_m2_GeV_minus2"])/gap
    # Work with A=m1 |Ue1 Um1|/m2 and B=|Ue2 Um2|.
    # |A²-B²|/(1-a)=|d0+d1 cosδ| <= h(A+B)/(1-a).
    AplusB_div_c13sq=math.sqrt(mu/(1-a1))
    d0_lo=-s1*(1-s1)*(1-t0)*(1-r0**2)+t0*a0*(r0**2*(1-s1)**2-s1**2)
    d0_hi=-s0*(1-s0)*(1-t1)*(1-r1**2)+t1*a1*(r1**2*(1-s0)**2-s0**2)
    assert r0**2*(1-s1)**2>s1**2
    d1_lo=(2*math.sqrt(s0*(1-s0)*min(t0*(1-t0),t1*(1-t1))*a0)
           *(r0**2+(1-r0**2)*s0))
    d1_hi=2*math.sqrt(s1*(1-s1)*.25*a1)*(r1**2+(1-r1**2)*s1)
    lim=h*AplusB_div_c13sq
    quotients=[n/d for n in [-lim-d0_hi,lim-d0_lo] for d in [d1_lo,d1_hi]]
    cmin=max(-1.,min(quotients)); cmax=min(1.,max(quotients))
    deg_lo=math.degrees(math.acos(cmax));deg_hi=math.degrees(math.acos(cmin))
    # Complementary phasor identity for |Mee|, fixing no Majorana phase.
    bmin=r0**2+(1-r0**2)*s0
    phi0=math.atan(math.sqrt(t1*a1/(1-t1)))
    denom=math.sqrt((1-a1)*bmin*(1-t1+t1*a0))
    phi=phi0+math.asin(min(1.,h/denom))
    floor=m20*(1-a1)*math.sqrt(bmin)*max(0.,math.cos(phi))
    return {"z_gap_GeV-2":gap,"h_max":h,"cos_delta_interval":[cmin,cmax],
            "abs_delta_interval_degrees":[deg_lo,deg_hi],"Mee_lower_eV":floor,
            "generic_endpoint_mbeta_min_eV":m20*math.sqrt((1-a1)*bmin),
            "tube_needed_for_exactly_CP_even": min(abs(d0_lo+d1_lo),abs(d0_hi-d1_lo))
            / AplusB_div_c13sq if max(abs(d0_lo),abs(d0_hi))<d1_lo else 0.,
            "Tyrs_upper_at_weak_122meV":3.8e26*(.122/floor)**2}


def physical_point(io):
    s,a,t=[math.sqrt(x) for x in (io["s12"][0],io["s13"][0],
                                  load("1580/inputs.json")["theta23_sin2"]["IO"][0])]
    c,ca,u=[math.sqrt(1-x*x) for x in (s,a,t)]
    m2=math.sqrt(io["atm_abs"][0]); m1=math.sqrt(m2*m2-io["solar"][0])
    numerator=m2**2*s*s*(c*c*u*u+s*s*t*t*a*a)-m1**2*c*c*(s*s*u*u+c*c*t*t*a*a)
    denom=2*s*c*t*u*a*(m1**2*c*c+m2**2*s*s)
    delta=-math.acos(numerator/denom)
    d=cmath.exp(1j*delta)
    e1,e2=c*ca,s*ca
    m1v,m2v=-s*u-c*t*a*d,c*u-s*t*a*d
    phase=-m1*e1*m1v/(m2*e2*m2v)
    assert abs(abs(phase)-1)<1e-14 and abs(m1*e1*m1v+m2*phase*e2*m2v)<1e-17
    mee=abs(m1*e1**2+m2*phase*e2**2)
    mmu=abs(m1*m1v**2+m2*phase*m2v**2)
    me,mu,tau=.00051099895,.1056583755,1.77709
    old=load("1594/independent_results.json")["IO_witness"]
    rates={}
    for name,m,weight,gain in zip(["tau_e","tau_mu"],[me,mu],[a,t*ca],old["positive_k_GeV-2"]):
        norm=.0072973525205055605*tau**3*290.17e-15*(1-(m/tau)**2)**3/(4*6.582119569e-25)
        rates[name]=norm*(gain*m*weight/(u*ca))**2
    assert rates["tau_e"] < 3.3e-8 and rates["tau_mu"] < 4.4e-8
    return {"delta_deg":math.degrees(delta),"relative_mass_coefficient_phase_deg":math.degrees(cmath.phase(phase)),
            "Mee_eV":mee,"M_mumu_magnitude_eV":mmu,"masses_eV":[m1,m2,0.],"LFV_rates":rates}


def report():
    io=load("1370/inputs.json")["IO"]
    s0,s1=io["s12"][1:]; a0,a1=io["s13"][1:]
    t0,t1=load("1580/inputs.json")["theta23_sin2"]["IO"][1:]
    m20,m21=map(math.sqrt,io["atm_abs"][1:])
    d0,d1=io["solar"][1:]
    ranges=s0,s1,a0,a1,t0,t1,math.sqrt(1-d1/m20**2),math.sqrt(1-d0/m21**2),m20,m21
    result={case["label"]:evaluate(case,ranges) for case in load("1596/inputs.json")["scenarios"]}
    point=physical_point(io)
    mee=point["Mee_eV"]
    for label,d in result.items():
        dlow,dhi=d["abs_delta_interval_degrees"]
        assert dlow<=abs(point["delta_deg"])<=dhi and d["Mee_lower_eV"]<=mee
    # Literature endpoints are different matrix-element dictionaries for one measurement.
    xenon=load("1370/inputs.json")["kamland_zen"]
    nme=xenon["m_bb_upper_eV_endpoints_by_nme"]
    T=xenon["half_life_lower_years"]
    return {"round":1596,"cases":result,"matched_IO_point":point,
            "IO_example_absolute_delta_deg":abs(point["delta_deg"]),"IO_example_Mee_eV":mee,
            "IO_example_KZ_half_life_yrs":[T*(edge/mee)**2 for edge in nme],
            "IO_example_min_cancelling_amplitude_fraction_for_28meV":max(0.,1-.028/mee),
            "IO_published_delta_interval_degrees":
               [180*v for v in load("1416/inputs.json")["phase_credible_arcs_pi"]["IO_99.73"]]}


def compare(actual,saved):
    if isinstance(actual,dict):
        assert saved.keys()==actual.keys()
        for k in actual: compare(actual[k],saved[k])
    elif isinstance(actual,list):
        assert len(actual)==len(saved)
        for a,b in zip(actual,saved): compare(a,b)
    else:
        assert math.isclose(actual,saved,rel_tol=2.e-11,abs_tol=2.e-25)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    cli.add_argument("--write",action="store_true")
    args=cli.parse_args()
    out=report()
    if args.write:
        with (ROOT/"independent_results.json").open("x",encoding="utf-8",newline="\n") as f:
            json.dump(out,f,indent=2,ensure_ascii=False)
            f.write("\n")
    elif args.check:
        compare(out,load("1596/independent_results.json"))
    print(json.dumps(out,ensure_ascii=False))
