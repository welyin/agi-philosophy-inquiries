"""Independent simultaneous dipole and Majorana test at the matched scale."""
import argparse
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
E, MU, TAU = .00051099895, .1056583755, 1.775
ALPHA, HBAR, VEV = .0072973525205055605, 6.582119569e-25, 246.22
UNITS = {"me": (E,MU,2.1969803e-6,1.5e-13),
         "te": (E,TAU,288.e-15,3.3e-8), "tm":(MU,TAU,288.e-15,4.4e-8)}


def reading(path):
    return json.loads((ARCH/path).read_text(encoding="utf-8-sig"))


def rate(i,j,life):
    return ALPHA*j**3*(1-(i/j)**2)**3*life/(4*HBAR)


def caps():
    return {name:math.sqrt(br/rate(i,j,t)) for name,(i,j,t,br) in UNITS.items()}


def radiative_angle(osc):
    return [math.sqrt((1-math.sqrt(1-amp))/2) for amp in
            (osc["DayaBay_amplitude"][0]+i*osc["DayaBay_amplitude"][1]
             for i in (-2,2))]


def angles(osc, order):
    def endpoints(values):
        x=values[0]
        return math.sqrt(x-2*values[1]),math.sqrt(x+2*values[-1])
    return (*endpoints(osc[order]["solar_sine_squared"]),
            *endpoints(osc[order]["atmospheric_sine_squared"]),*radiative_angle(osc))


def electron_bounds(osc,order,b):
    sl,sh,tl,th,al,ah=angles(osc,order)
    if order=="NO":
        e_min=math.sqrt((1-sh**2)*(1-ah**2))
        tau_max=sh*th+math.sqrt((1-sl**2)*(1-tl**2))*ah
        mu_max=sh*math.sqrt(1-tl**2)+math.sqrt(1-sl**2)*th*ah
    else:
        e_min=al
        tau_max=math.sqrt((1-tl**2)*(1-al**2))
        mu_max=th*math.sqrt(1-al**2)
    # A row defect |n†C[:,e]| is in GeV^-1. Signed diagonals may have
    # imaginary parts; taking only their measured real part is conservative.
    cases={}
    for name,ae in {"Cs_two_width":3.024e-13,"Rb_one_width":1.176e-13}.items():
        for efficiency,rho,diagonal in [(1.,0.,0.),(.8,.2,1.e-13)]:
            spill=e_min*max(0.,ae-diagonal)/(2*E)-mu_max*b["me"]/(math.sqrt(efficiency)*(1-rho))
            amplitude=max(0.,spill)/tau_max
            cases[f"{name}_q_{efficiency}"]={
                "projected_e_column_defect_to_erase_GeV-1":max(0.,spill),
                "visible_Br_tau_e_min":efficiency*(1-rho)**2*rate(E,TAU,288e-15)*amplitude**2}
    return cases


def normal_bound(osc,b):
    sl,sh,tl,th,al,ah=angles(osc,"NO")
    cs,ct=math.sqrt(1-sl**2),math.sqrt(1-tl**2)
    g21=osc["NO"]["solar_gap"];g32=osc["NO"]["atmospheric_gap"]
    sol=g21[0]+2*g21[1]; atm=g32[0]-2*g32[1]
    q=math.sqrt(sol/(sol+atm))
    # Enclose |N_ij|/m3 in a rectangular oscillation box.
    em=math.sqrt(1-ah**2)*(al*tl-q*sh*(cs*ct+sh*th*ah))
    ee=ah**2+q*sh**2
    mm=th**2+q*(cs*ct+sh*th*ah)**2
    et=ah*ct+q*sh*(cs*th+sh*ct*ah)
    mt=th*ct+q*(cs*ct+sh*th*ah)*(cs*th+sh*ct*ah)
    # Both mu->egamma chiralities share a single experimental upper bound.
    rhs=b["me"]*math.hypot(ee/E,mm/MU)+(b["te"]*mt+b["tm"]*et)/TAU
    zmin=1.176e-13/(2*E**2)-1.64e-9/(2*MU**2)
    qeff=.8*(1-.2)**2
    return {"q_hi":q,"N_em_over_m3_min":em,"N_ee_mm_et_mt_over_m3_max":[ee,mm,et,mt],
            "z_difference_needed_GeV-2":zmin,"z_difference_max_GeV-2":rhs/em,
            "stress_z_difference_max_GeV-2":rhs/(em*math.sqrt(qeff)),
            "minimum_stress_extra_electron_anomaly":
                2*E**2*max(0.,zmin-rhs/(em*math.sqrt(qeff))),
            "minimum_stress_skew_residual_over_m3_GeV-2":
                max(0.,em*zmin-rhs/math.sqrt(qeff)),
            "rate_efficiency_needed_at_rho_02":(rhs/(em*zmin*.8))**2}


def pmns(s,c,t,ct,a,ca,delta):
    d=np.exp(1j*delta)
    return np.array([[c*ca,s*ca,a*d.conjugate()],
          [-s*ct-c*t*a*d,c*ct-s*t*a*d,t*ca],
          [s*t-c*ct*a*d,-c*t-s*ct*a*d,ct*ca]])


def mass_function(x):
    return x*x.ln()/(x-1) if x!=1 else Decimal(1)


def physical_loops(gain, mass, neutral_coeff):
    # Same mass ratios, one dark parity per species, independent overall scales.
    with localcontext() as precision:
        precision.prec=60
        x1,x2=Decimal("1.21"),Decimal("2.56")
        def mag(x):
            d=1-x
            return (3*x-1)/d**2+2*x*x*x.ln()/d**3+2*(1/d+x*x.ln()/d**2)
        df = mass_function(x2)-mass_function(x1)
        response = float((mag(x1)-mag(x2))/(2*df))
        fermion=math.sqrt(response/gain)
        sin2=min(.35, .6*VEV/fermion)
        sin_sq=(1-math.sqrt(1-sin2*sin2))/2
        x0=x1+(x2-x1)*Decimal(str(sin_sq))
        ma=fermion*sin2*float(df)/(32*math.pi**2)
        yuk_left=.07 if mass==E else .4
        yuk_right=mass/(ma*yuk_left)
        target=Decimal(str(abs(neutral_coeff)*1.e-9*32*math.pi**2/(fermion*yuk_left**2)))
        lo,hi=Decimal(0),Decimal(".1")
        for _ in range(155):
            h=(lo+hi)/2
            if mass_function(x0+h/2)-mass_function(x0-h/2)<target:
                lo=h
            else:
                hi=h
        split=(lo+hi)/2
        reconstructed=fermion*yuk_left**2*float(mass_function(x0+split/2)-mass_function(x0-split/2))/(32*math.pi**2)*1e9
        assert math.isclose(reconstructed,abs(neutral_coeff),rel_tol=1e-12)
        return {"M_Majorana_GeV":fermion,"charged_poles_GeV":[1.1*fermion,1.6*fermion],
                "sin_2theta":sin2,"H_eta_chi_Trilinear_over_M":
                    sin2*float(x2-x1)*fermion/(math.sqrt(2)*VEV),
                "left_Yukawa_multiplier":yuk_left,"right_Yukawa_magnitude":yuk_right,
                "lambda5":float(split)*fermion**2/VEV**2,
                "Ma_diagonal_magnitude_eV":reconstructed}


def inverted_example(osc):
    o=osc["IO"]
    s,t=math.sqrt(o["solar_sine_squared"][0]), math.sqrt(o["atmospheric_sine_squared"][0])
    a=math.sqrt((1-math.sqrt(1-osc["DayaBay_amplitude"][0]))/2)
    c,ct,ca=[math.sqrt(1-x*x) for x in (s,t,a)]
    m1=math.sqrt(o["atmospheric_gap"][0])
    m2=math.sqrt(o["atmospheric_gap"][0]+o["solar_gap"][0])
    num=m2**2*s*s*(c*c*ct*ct+s*s*t*t*a*a)-m1**2*c*c*(s*s*ct*ct+c*c*t*t*a*a)
    den=2*s*c*t*ct*a*(m1*m1*c*c+m2*m2*s*s)
    delta=-math.acos(num/den)
    u=pmns(s,c,t,ct,a,ca,delta)
    v,w,n=u.T
    ph=-(m1*v[0]*v[1])/(m2*w[0]*w[1])
    assert abs(abs(ph)-1)<1e-13
    nu=m1*np.outer(v,v)+m2*ph*np.outer(w,w)
    S=np.eye(3,dtype=complex)[:,:2]
    S[2]=-n[:2].conj()/n[2].conjugate()
    assert np.linalg.norm(nu-S@np.diag(np.diag(nu)[:2])@S.T)<1.e-15
    tau_pole=1.77709
    m=np.diag([E,MU,tau_pole])
    a_e,a_mu=3.38e-13,.25e-9
    k=np.array([a_e/(2*E**2),a_mu/(2*MU**2)])
    mass_loop=S@m[:2,:]
    tree=m-mass_loop
    assert np.linalg.matrix_rank(tree, tol=1.e-12)==1
    dip=S@np.diag(k)@m[:2,:]
    x=dip@np.diag(1/np.diag(m))
    assert np.linalg.norm(x@nu-nu@x.T)<1.e-22 and np.linalg.norm(n.conj()@dip)<1.e-23
    rates={}
    for label,(i,j) in {"me":(0,1),"te":(0,2),"tm":(1,2)}.items():
        mi,mj,tau,limit=UNITS[label]
        if j==2: mj,tau=tau_pole,290.17e-15
        rates[label]=rate(mi,mj,tau)*(abs(dip[i,j])**2+abs(dip[j,i])**2)
        assert rates[label]<limit
    return {"delta_deg":math.degrees(delta),"Majorana_relative_coefficient_deg":np.angle(ph,deg=True),
            "mass_spectrum_eV":[m1,m2,0.],"anomalies":[a_e,a_mu],"rates":rates,
            "positive_k_GeV-2":k.tolist(),"tree_frobenius_GeV":float(np.linalg.norm(tree)),
            "channels":[physical_loops(ki,mi,nu[i,i]) for i,(ki,mi) in enumerate(zip(k,[E,MU]))]}


def contrast_cosets():
    # A single possible complex population and two lower needs on the two
    # chiralities of e-mu. Both pass their individual cap, yet jointly fail.
    ne,nm,nt=.8,.36,.48
    assert math.isclose(ne**2+nm**2+nt**2,1.)
    e_diag,mu_diag=.5,2.
    e_tau,mu_tau=(ne*e_diag-.8*nm)/nt,(nm*mu_diag-.8*ne)/nt
    needed=math.hypot((ne*e_diag-nt*e_tau)/nm,(nm*mu_diag-nt*mu_tau)/ne)
    assert e_tau>0 and mu_tau>0 and needed>1.
    return {"null_moduli":[ne,nm,nt],"fraction_of_MEG_amplitude_required":needed}


def calculate():
    osc=reading("1386/inputs.json")
    b=caps()
    return {"r":1594,"LFV_caps_GeV-1":b,"electron_column":{
        order:electron_bounds(osc,order,b) for order in ("NO","IO")},
        "NO_two_matched_mediators":normal_bound(osc,b),"joint_phase_space_example":contrast_cosets(),"IO_witness":inverted_example(osc)}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): compare(x,y)
    else:
        assert math.isclose(a,b,rel_tol=2.e-11,abs_tol=1.e-32),(a,b)


if __name__=="__main__":
    cmd=argparse.ArgumentParser()
    cmd.add_argument("--check",action="store_true")
    cmd.add_argument("--write",action="store_true")
    args=cmd.parse_args()
    obj=calculate()
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf8",newline="\n") as buffer:
            json.dump(obj,buffer,ensure_ascii=False,indent=2)
            buffer.write("\n")
    elif args.check:
        compare(obj,json.loads(path.read_text(encoding="utf8")))
    print(json.dumps(obj,ensure_ascii=False))

