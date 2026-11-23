"""A charged B width, neutral B_d oscillation and dimuon width on the actual flavor endpoints."""
import argparse
import cmath
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
def read(n):return json.loads((HERE/n).read_bytes())
CFG=read("inputs.json")
WEAK=read(CFG["input_scenarios"])
SD=read(CFG["dimuon_bag_shortdistance"])
TAU=read(CFG["tau"])
MIX=read(CFG["oscillation_measurements"])["measurements_in_order_as_source_urls"]
CP=read(CFG["golden_and_gamma"])


def band(values,k):
    width=math.hypot(*values[1:])
    return np.array([values[0]-k*width,values[0]+k*width])


def matrix(us,ub,cb,gamma_deg):
    c13=math.sqrt(1-ub**2)
    s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12**2),math.sqrt(1-s23**2)
    g=math.radians(gamma_deg)
    phase=g+math.asin(c12*s23*ub*math.sin(g)/(s12*c23))
    e=cmath.exp(1j*phase)
    V=np.array([[c12*c13,s12*c13,ub*e.conjugate()],
                [-s12*c23-c12*s23*ub*e,c12*c23-s12*s23*ub*e,s23*c13],
                [s12*s23-c12*c23*ub*e,-c12*s23-s12*c23*ub*e,c23*c13]])
    np.testing.assert_allclose(V@V.conj().T,np.eye(3),atol=4e-16)
    return V


def invariants(V):
    a,b,c=V[:,0]*V[:,2].conjugate()
    beta=cmath.phase(-b/c)
    gamma=cmath.phase(-a/b)
    alpha=cmath.phase(-c/a)
    ratio=abs(V[0,2]/c)
    np.testing.assert_allclose(ratio*abs(V[0,0]),math.sin(beta)/math.sin(gamma),rtol=2e-13)
    return dict(ud=abs(V[0,0]),ub=abs(V[0,2]),cb=abs(V[1,2]),bd=abs(c),
                ratio_squared=ratio**2,alpha=math.degrees(alpha),beta=math.degrees(beta),gamma=math.degrees(gamma))


def factors(k):
    # Coefficient boxes, simultaneously favorable signs. Isospin uses the published difference.
    me=TAU["meson_masses_GeV"]["B"]
    mt=TAU["mu_tau_masses_GeV"][1]
    massbd=band(MIX["mass_d_MeV"],k)/1000
    short=SD["short_distance"]
    S=band(short["S0_published"],k)
    eta=band(short["etaB"],k)
    bag=band(SD["bag"]["d"],k)
    fplus=band(TAU["decay_constants_MeV"]["B"],k)
    diff=band(CFG["fneutral_minus_fcharged_MeV"],k)
    f=1/(1+diff[::-1]/fplus)
    dmass=band(MIX["dm_d_inverse_ps"],k)
    taup=band(TAU["B_lifetime_ps"],k)
    taud=band(SD["tau_ps"]["d"],k)
    c=SD["const_GeV"]
    mw=c["MW"]; mu=c["m_mu"]; gf=c["GF"]
    axial=band(short["CA_at5GeV"],k)
    beta=np.sqrt(1-(2*mu/massbd)**2)
    tau_pref=3*math.pi/4*(mt/mw)**2*(1-(mt/me)**2)**2*me
    mu_pref=3*(gf*mw*mu)**2/math.pi**3
    ktau=tau_pref*f**2/(massbd*bag*eta*S)[::-1]
    kmu=mu_pref*axial**2*beta/(bag*eta*S)[::-1]
    # Direct tau/mumu removes the mixing inputs before any uncertainty propagation.
    both=tau_pref/mu_pref*f**2*taup/(taud*massbd*axial**2*beta)[::-1]
    return dict(f_ratio=f,qtau=ktau,br_tau_unit=ktau*dmass*taup,
                br_mumu=kmu*dmass*taud,br_tau_to_mumu_unit=both)


def free_family(ud,beta,gamma):
    z=math.sin(beta)/math.sin(beta+gamma)*cmath.exp(1j*gamma)
    rows=[]
    for t in (.025,.045,.065):
        a=abs(z)**2*t*t*(1-t*t)/abs(1-z*t*t)**2
        c12sq=(ud**2+a)/(1+a)
        sx=z*t*math.sqrt((1-c12sq)*(1-t*t)/c12sq)/(1-z*t*t)
        ub=abs(sx); c13=math.sqrt(1-ub**2)
        V=matrix(math.sqrt(1-c12sq)*c13,ub,t*c13,math.degrees(gamma))
        inv=invariants(V)
        np.testing.assert_allclose([inv["ud"],inv["beta"],inv["gamma"]],
                                   [ud,math.degrees(beta),math.degrees(gamma)],rtol=4e-14)
        rows.append(dict(cb=inv["cb"],ub=inv["ub"],us=abs(V[0,1]),leptonic_to_mixing_ratio_squared=inv["ratio_squared"]))
    return rows


def beta_branches(k):
    S=band(CP["S_d"],k)
    C=band(CP["C_d"],k)
    abs_c=np.array([0 if C[0]<=0<=C[1] else min(abs(C)),max(abs(C))])
    small=.5*np.arcsin(S/np.sqrt(1-abs_c**2))
    return [small,np.pi/2-small[::-1]]


def add_absolute_decay_constant(br,ud,beta,gamma):
    """Use the absolute tree normalization as extra information beyond the width ratio."""
    mass=TAU["meson_masses_GeV"]["B"]
    tau_mass=TAU["mu_tau_masses_GeV"][1]
    f=TAU["decay_constants_MeV"]["B"][0]/1000
    time=TAU["B_lifetime_ps"][0]/CFG["hbar_GeV_ps"]
    norm=time*SD["const_GeV"]["GF"]**2*mass*tau_mass**2*f**2*(1-(tau_mass/mass)**2)**2/(8*math.pi)
    ub=math.sqrt(br/norm)
    z=math.sin(beta)/math.sin(beta+gamma)*cmath.exp(1j*gamma)
    left=(1-ub**2-ud**2)*abs(z)**2
    right=ub**2*ud**2
    roots=np.roots([left+right*abs(z)**2,-left-2*right*z.real,right])
    rows=[]
    for v in sorted(roots):
        if not np.isreal(v) or not 0<v<1:
            continue
        cb=math.sqrt(v*(1-ub**2))
        us=math.sqrt(1-ub**2-ud**2)
        # Reconstruct delta from the complex triangle, including the near-maximal s23 branch.
        base=1-z*v
        delta=cmath.phase(z/base)
        c12=ud/math.sqrt(1-ub**2)
        s12=us/math.sqrt(1-ub**2)
        t=math.sqrt(v)
        z_check=c12*ub*cmath.exp(1j*delta)/(t*(s12*math.sqrt(1-v)+c12*t*ub*cmath.exp(1j*delta)))
        np.testing.assert_allclose(z_check,z,rtol=1e-8)
        rows.append(dict(cb=cb,ub=ub,delta_deg=math.degrees(delta)))
    return rows


def inference():
    experiments={"Belle_2015_combination":CFG["BR_Belle_combined"],"Belle_II_2025":TAU["b_branching"]}
    rows={}
    for exp,sample in experiments.items():
        results={}
        for k in (0,1,2):
            F=factors(k)
            raw=band(sample,k)
            ratio=np.sqrt(raw/F["br_tau_unit"][::-1])
            gp=CP["time_integrated_gamma_deg"]
            gamma=np.radians([gp["central"]]*2 if k==0 else gp["68" if k==1 else "95"])
            beta=beta_branches(k)
            alternatives=[]
            for bi in beta:
                sn=np.sin(bi)/np.sin(gamma[::-1])
                ud=sn/ratio[::-1]
                alternatives.append(dict(beta_deg=np.degrees(bi),ud_interval_before_unitarity=ud,
                                    consistent_unitary=bool(ud[1]>0 and ud[0]<1)))
            rec=dict(weak_ratio=ratio,angles_and_ud=alternatives)
            if k==0:
                ud=alternatives[0]["ud_interval_before_unitarity"][0]
                if ud<1:
                    rec["same_rates_angles_different_cb"]=free_family(ud,beta[0][0],gamma[0])
                    rec["using_extra_absolute_f_examples"]=add_absolute_decay_constant(sample[0],ud,beta[0][0],gamma[0])
            results[str(k)]=rec
        rows[exp]=results
    return rows


def comparison(V):
    inv=invariants(V); r2=inv["ratio_squared"]
    F=factors(0)
    bpred=F["br_tau_unit"][0]*r2
    mu=F["br_mumu"][0]
    rare=SD["LHCb_run3_2024_25"]["d_in_1e10"]
    expmu=rare[0]*1e-10
    levels={}
    for exp,data in [("Belle",CFG["BR_Belle_combined"]),("Belle_II",TAU["b_branching"])]:
        expt=data[0]
        closure=(expt/expmu)/(r2*F["br_tau_to_mumu_unit"][0])
        np.testing.assert_allclose(closure,(expt/bpred)/(expmu/mu),rtol=5e-15)
        tauBox=band(data,1)
        dimuBox=np.array([rare[0]-math.hypot(rare[1],rare[3]),rare[0]+math.hypot(rare[2],rare[4])])*1.e-10
        hbox=(tauBox/dimuBox[::-1])/(r2*factors(1)["br_tau_to_mumu_unit"][::-1])
        levels[exp]=dict(branch_difference_in_exp_widths=(bpred-expt)/math.hypot(*data[1:]),
                       tau_over_mixing_relative_to_SM=expt/bpred,
                       tau_vs_dimuon_closure=closure,closure_marginal_box=hbox,
                       relative_decay_amplitudes_needed_even_if_mixing_free=math.sqrt(closure))
    return dict(CKM=inv,tau_BR=bpred,dimuon_BR=mu,
                BR_tau_divided_by_BR_dimuon=bpred/mu,rate_tests=levels,
                tau_BR_nuisance_box=factors(1)["br_tau_unit"]*r2,
                direct_BR_ratio_nuisance_box=factors(1)["br_tau_to_mumu_unit"]*r2)


def dataset():
    engine=runpy.run_path(str(HERE/CFG["joint"]))
    def resolve(u):
        row=engine["shoot"](HERE,1/WEAK["QCDSF_s_d"][0],u,WEAK["inclusive_cb"][0],.25)
        m=row["low_moduli"]
        return matrix(m["Vus"],m["Vub"],m["Vcb"],row["gamma_deg"])
    orig=read("../1755/results.json")["full_two_threshold_shooting"]["QCDSF"]
    reactor=read(CFG["reactor_calibrations"])
    sets={}
    for label,p in orig.items():
        m=p["low_moduli"]
        sets["QCDSF_B_"+label]=matrix(m["Vus"],m["Vub"],m["Vcb"],p["gamma_deg"])
    for label,p in {"central":reactor["central"],"reactor_high":reactor["only_reactor_at_upper_one_width"],
                    **reactor["eight_one_width_vertices"]}.items():
        sets["QCDSF_reactor_"+label]=matrix(p["us"],p["ub"],p["cb"],p["gamma_deg"])
    inc=resolve(WEAK["inclusive_ub"][0])
    np.testing.assert_allclose(inc,sets["QCDSF_B_central"],atol=6.e-10)
    profile=read(CFG["clock_recalibration"])["combined_central_solution"]
    sets["QCDSF_reprofiled_two_clocks"]=resolve(profile["Vub"])
    return {label:comparison(V) for label,V in sets.items()}


def precision(a,b,key,width=.0):
    x,y=sorted([a[key],b[key]],key=lambda s:s[0])
    bottom=y[0]*(1-width)
    top=x[1]*(1+width)
    return (bottom-top)/(bottom+top)


def calc():
    points=dataset()
    high=points["QCDSF_B_central"]
    low=points["QCDSF_reprofiled_two_clocks"]
    n=CFG["stress_relative_rate_budget"]
    tests={}
    for name,box in (("tau_over_mixing","tau_BR_nuisance_box"),("tau_over_dimuon","direct_BR_ratio_nuisance_box")):
        tests[name]=dict(separating_fractional_experimental_budget=precision(low,high,box),
                         threshold_with_optional_extra_relative_rate_budget=precision(low,high,box,n))
    return dict(points=points,coefficients={str(k):factors(k) for k in (0,1,2)},
                from_tau_mixing_and_phase_only=inference(),resolution=tests,
                input_point_spans={prefix:[min(vals),max(vals)] for prefix in
                    ("QCDSF_B_","QCDSF_reactor_") for vals in
                    [[v["tau_BR"] for key,v in points.items() if key.startswith(prefix)]]},
                fractional_splitting_high_to_low=high["tau_BR"]/low["tau_BR"]-1)


def tidy(z):
    if isinstance(z,dict):return {k:tidy(v) for k,v in z.items()}
    if isinstance(z,(tuple,list)):return [tidy(v) for v in z]
    if isinstance(z,np.ndarray):return tidy(z.tolist())
    if isinstance(z,(float,np.floating)):return float(f"{z:.10g}")
    return z


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=tidy(calc())
    target=HERE/"results.json"
    if args.check:
        assert result==json.loads(target.read_bytes())
        print("1763 OK: physical CKM, charged-neutral isospin and dimuon cancellation, phase inversion and CKM continuum.")
    else:
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

