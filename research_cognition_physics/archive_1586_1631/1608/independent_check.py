"""Independent reverse running, high-field interpolation, and real-singlet control."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial as Poly

OUT = Path(__file__).with_name("independent_results.json")
REFERENCE = 173.34
GAUGE = (.35830,.64779,1.16660)
RUN = (41/6,-19/6,-7)
K = 16*math.pi**2
LAMBDA = .12604


def gauge(t,initial):
    return [g*g/(1-2*b*g*g*t/K) for g,b in zip(initial,RUN)]


def reverse(t_end,n,initial):
    gy,gw,_ = gauge(t_end,initial)
    u = .25*math.sqrt(2*gw*gw+(gw+gy)**2)
    def flow(t,z):
        gy,gw,gs=gauge(t,initial)
        q,l=z
        return (q*(9*q-17*gy/6-9*gw/2-16*gs)/K,
                (24*l*l+(12*q-9*gw-3*gy)*l-6*q*q+3/8*(2*gw*gw+(gw+gy)**2))/K)
    z=(u,0.)
    dt=-t_end/n
    for i in range(n):
        time=t_end*(1-i/n)
        def step(q,h):
            return tuple(c+h*v for c,v in zip(z,q))
        a=flow(time,z)
        b=flow(time+dt/2,step(a,dt/2))
        c=flow(time+dt/2,step(b,dt/2))
        d=flow(time+dt,step(c,dt))
        z=tuple(c0+dt*(a0+2*b0+2*c1+d0)/6 for c0,a0,b0,c1,d0 in zip(z,a,b,c,d))
    return math.sqrt(z[0]),z[1]


def critical(n,ref=REFERENCE,initial=GAUGE,quartic=LAMBDA):
    lo,hi=(math.log(m/ref) for m in [1e10,1.22e19])
    endpoints=[reverse(t,n,initial)[1] for t in (lo,hi)]
    assert (endpoints[0]-quartic)*(endpoints[1]-quartic)<0
    for _ in range(40):
        m=(lo+hi)/2
        if (reverse(m,n,initial)[1]-quartic)*(endpoints[0]-quartic)>0:
            lo=m
        else:
            hi=m
    t=(lo+hi)/2
    y,l=reverse(t,n,initial)
    return {"mu0_GeV":ref,"low_y":y,"low_quartic":l,"critical_log":t,
            "formal_scale_GeV":ref*math.exp(t),
            "required_quartic_at_log_endpoints":endpoints}


def high_field_jets():
    a,b,c=3e-10,-2e-10,1e-10
    w=.125
    y=Poly([0,1])
    d=y-1
    f=y**5*(a+(b-5*a)*d+(c/2-5*b+15*a)*d**2)
    got=[f.deriv(i)(1.) for i in range(3)]
    low=[f.deriv(i)(0.) for i in range(5)]
    np.testing.assert_allclose(got,[a,b,c],atol=4e-22)
    assert low==[0]*5
    # Potential in units Lambda^4; z=(H†H-v²/2)/Lambda² and y=z/w.
    wilsons={str(2*i):coef/w**i for i,coef in enumerate(f.coef) if coef}
    assert max(map(abs,wilsons.values()))<.1
    return {"field_fraction_z":w,"low_jet":low,"high_jet_y":got,
            "dimensionless_coeff_in_shifted_basis_by_operator_dimension":wilsons}


def singlet():
    f=1e12
    v=246.22
    h_mass=125.25
    lam=h_mass**2/(2*v*v)
    k=.25
    b=.2
    h_high=math.sqrt((v/f)**2+2*k)
    ss=Poly([0,1])
    records={}
    for eps in [-.002,0,.002]:
        valley=b*ss**2*(ss-1)**2+eps*ss**3*(4-3*ss)
        critical_points=sorted(valley.deriv().roots())
        saddle=b/(2*(b-3*eps))
        np.testing.assert_allclose(critical_points,[0,saddle,1],atol=1e-13)
        np.testing.assert_allclose(valley(1.),eps,atol=1e-15)
        high_hess=2*lam*np.outer([h_high,-2*k],[h_high,-2*k])
        high_hess[1,1]+=valley.deriv(2)(1.)
        high_eigen=np.linalg.eigvalsh(high_hess)
        assert min(high_eigen)>0 and valley.deriv(2)(saddle)<0
        assert valley(saddle)>max(0,eps)
        np.testing.assert_allclose(valley.deriv(2)(0.),2*b)
        records[str(eps)]={"low_s_mass_GeV":math.sqrt(2*b)*f,"lambda":lam,
            "high_h_over_f":h_high,"high_S_over_f":1.,"gap_over_f4":valley(1.),
            "saddle_s_over_f":saddle,"barrier_above_low_over_f4":valley(saddle),
            "high_radial_mass2_over_f2":high_eigen.tolist(),
            "ssss_contact_at_low_vacuum":24*(b+lam*k*k-3*eps)}
    return {"F_GeV":f,"kappa":k,"b":b,"v_GeV":v,"tree_higgs_mass_GeV":h_mass,
            "epsilon_family":records}


def evaluate():
    slow=critical(300)
    fast=critical(150)
    assert abs(slow["low_y"]-fast["low_y"])<3e-8
    assert abs(slow["critical_log"]-fast["critical_log"])<1e-5
    newer=critical(400,200.,(.358852,.646832,1.15251),.123533)
    return {"reverse_SM":slow,"reverse_SM_200GeV":newer,
            "coarse_fine_log_difference":abs(slow["critical_log"]-fast["critical_log"]),
            "single_Higgs_considered_EFT":high_field_jets(),"renormalizable_control":singlet()}


def verify(x,y):
    if isinstance(y,dict):
        assert x.keys()==y.keys()
        for k in y:
            verify(x[k],y[k])
    else:
        np.testing.assert_allclose(x,y,rtol=1e-9,atol=5e-12)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group()
    g.add_argument("--check",action="store_true")
    g.add_argument("--write",action="store_true")
    args=p.parse_args()
    data=evaluate()
    if args.write:
        with OUT.open("x",encoding="utf8") as f:
            json.dump(data,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        verify(data,json.loads(OUT.read_text(encoding="utf8")))
        print("1608 independent: reverse critical flow, Hermite data, renormalizable vacuum family passed")
