"""Dirac charged-singlet loop: sharp right-handed flavor fraction and full-box tests."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def extent(entry,n):
    return entry[0]-n*entry[1],entry[0]+n*entry[-1]


def mixing(s12sq,s13sq,s23sq,phase):
    s12,s13,s23=np.sqrt([s12sq,s13sq,s23sq])
    c12,c13,c23=np.sqrt(1-np.array([s12sq,s13sq,s23sq]))
    e=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/e],
                     [-s12*c23-c12*s23*s13*e,c12*c23-s12*s23*s13*e,s23*c13],
                     [s12*s23-c12*c23*s13*e,-c12*s23-s12*c23*s13*e,c23*c13]])


def reactor(amplitude):
    return (1-math.sqrt(1-amplitude))/2


def full_box(cfg,data,ordering):
    n=cfg["widths"]
    osc=data[ordering]
    lo13,hi13=map(reactor,extent(data["DayaBay_amplitude"],n))
    _,hi12=extent(osc["solar_sine_squared"],n)
    lo23,hi23=extent(osc["atmospheric_sine_squared"],n)
    _,solar_max=extent(osc["solar_gap"],n)
    atm_min,_=extent(osc["atmospheric_gap"],n)
    e,mu,tau=cfg["charged_masses_GeV"]
    r=(tau/cfg["scalar_window_GeV_min"])**2
    loop_loss=-r*math.log(r)/(1-r)**2
    lo,hi=cfg["charged_mass_transport_each"]
    relative=hi/(lo*(1-loop_loss))
    if ordering=="NO":
        e0_min=(1-hi12)*(1-hi13)
        hratio=(hi23*(1-lo13)+solar_max/(solar_max+atm_min))/((1-hi23)*(1-hi13))
        left_min=e0_min/2
    else:
        e0_min=lo13
        hratio=(1-lo23*(1-hi13)+solar_max/atm_min)/(1-(1-lo23)*(1-lo13))
        left_min=(1-(1-lo23)*(1-lo13))/2
    k=mu/tau*relative*math.sqrt(hratio)
    ne_ratio=math.sqrt((1-e0_min)/e0_min)
    rho_mu=e/mu*relative*ne_ratio
    rho_tau=e/tau*relative*ne_ratio
    slope2=k*k+(rho_tau+k*rho_mu)**2
    rmax=slope2/(1+slope2)
    lower=left_min*(1-r)**2
    # Phase-space suppression cannot increase the tau share. For L the maximum
    # after weighting is also <= its massless value, hence at most 1/2.
    frac=cfg["illustration_inclusive_tau_fraction"]
    return {"max_loop_transfer_loss":loop_loss,"min_zero_e_fraction":e0_min,
            "H_mu_mu_over_H_tau_tau_upper":hratio,"R_tau_upper":rmax,
            "L_tau_lower_with_phase_space":lower,
            "min_L_fraction_for_registered_inclusive_readout":(frac-rmax)/(.5-rmax)}


def flavor_strength(g):
    return np.sum(abs(g)**2,axis=1)


def benchmark(cfg,data,ordering):
    source=cfg["benchmark"]
    osc=data[ordering]
    e,mu,tau=ms=np.array(cfg["charged_masses_GeV"])
    m1,m2=source["charged_poles_GeV"]
    sine=source["sin_theta"]
    cosine=math.sqrt(1-sine**2)
    def log_kernel(m):
        r=(ms/m)**2
        return -np.log(r)/(1-r)
    c0=2*sine*cosine*math.log((m2/m1)**2)/(16*math.pi**2)
    d=ms*(log_kernel(m2)-log_kernel(m1))/math.log((m2/m1)**2)
    u=mixing(osc["solar_sine_squared"][0],reactor(data["DayaBay_amplitude"][0]),
             osc["atmospheric_sine_squared"][0],source["delta_rad"])
    z=0 if ordering=="NO" else 2
    n=u[:,z].conj()
    f=source["f_scale"]*np.array([[0,n[2],-n[1]],[-n[2],0,n[0]],[n[1],-n[0],0]])
    A=c0*f@np.diag(d)
    solar,atm=osc["solar_gap"][0],osc["atmospheric_gap"][0]
    sq=([0,solar,solar+atm] if ordering=="NO" else [atm,atm+solar,0])
    M=u@np.diag(np.sqrt(sq)*1e-9)
    G=np.linalg.pinv(A,rcond=1e-12)@M
    p=n/d
    p/=np.linalg.norm(p)
    assert np.linalg.norm(A@G-M)/np.linalg.norm(M)<1e-12
    assert np.linalg.norm(p.conj()@G)<1e-12
    Gsq=np.linalg.norm(G)**2
    sharp=abs(p)**2+flavor_strength(G)/Gsq
    # For any requested row solve for the null component that saturates its fraction.
    witnesses={}
    for i,label in enumerate(("e","mu","tau")):
        row=G[i]
        free=p[i].conjugate()*Gsq*row/(np.linalg.norm(row)**2)
        basis=G+np.outer(p,free)
        shares=flavor_strength(basis)/np.linalg.norm(basis)**2
        assert abs(shares[i]-sharp[i])<1e-12
        assert np.linalg.norm(A@basis-M)/np.linalg.norm(M)<1e-10
        if i==2:
            witnesses["max_tau"]=basis
    witnesses["without_e"]=G-np.outer(p,G[0]/p[0])
    witnesses["without_tau"]=G-np.outer(p,G[2]/p[2])
    left=flavor_strength(f)
    phase_space=(1-(ms/m1)**2)**2
    assert np.max(abs(left/left.sum()-(1-abs(n)**2)/2))<1e-14
    out={}
    for name,g in witnesses.items():
        right=flavor_strength(g)
        left_rate=4*cosine**2*left*phase_space
        right_rate=sine**2*right*phase_space
        br_mu_e=(1/source["inverse_alpha"]/(768*math.pi*source["fermi_GeV_minus2"]**2)
                 *(16*abs((f@f.conj().T)[1,0])**2*(cosine**2/m1**2+sine**2/m2**2)**2
                   +abs((g@g.conj().T)[1,0])**2*(sine**2/m1**2+cosine**2/m2**2)**2))
        out[name]={"R_tau":right[2]/sum(right),"max_abs_g":float(np.max(abs(g))),
                   "pole1_inclusive_tau":(left_rate[2]+right_rate[2])/sum(left_rate+right_rate),
                   "leading_mu_to_e_gamma":br_mu_e}
    J=float(np.imag(u[0,0]*u[1,1]*u[0,1].conjugate()*u[1,0].conjugate()))
    return {"charged_loop_effective_mass_GeV":d.tolist(),"C_abs":c0,
            "left_fractions":(left/sum(left)).tolist(),"sharp_R_share_suprema":sharp.tolist(),
            "sum_of_separately_maximized_shares":sum(sharp),"J_example":J,"witnesses":out}


def run():
    cfg=read(HERE/"inputs.json")
    data=read(HERE/cfg["oscillation_input"])
    limit=read(HERE/cfg["left_and_radiative_data"])["MEG_II"]["BR_90CL_upper"]
    result={}
    for name in ("NO","IO"):
        bounds=full_box(cfg,data,name)
        central=benchmark(cfg,data,name)
        assert abs(central["sum_of_separately_maximized_shares"]-2)<1e-12
        assert central["sharp_R_share_suprema"][2]<bounds["R_tau_upper"]<.01
        assert bounds["L_tau_lower_with_phase_space"]>.19
        assert central["witnesses"]["max_tau"]["leading_mu_to_e_gamma"]<limit
        result[name]={"uniform_bound":bounds,"benchmark":central}
    return result


def verify(x,y):
    if isinstance(x,dict):
        return x.keys()==y.keys() and all(verify(v,y[k]) for k,v in x.items())
    return np.allclose(x,y,atol=1e-30,rtol=2e-8)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cli.parse_args()
    path=HERE/"results.json"
    values=run()
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(values,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        assert verify(values,read(path))
        print("1538: full-box bounds, attainable sharp flavor extrema, common scalar widths and dipole checks pass")
