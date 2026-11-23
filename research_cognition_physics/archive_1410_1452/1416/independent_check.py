"""Two real, named flavor seeds, circular intervals and paired T-conjugate vacuum rates."""
import argparse
import json
from math import cos, isclose, pi, sin, sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
CONVERSION=1.26693268  # Delta_ij = conversion * (dm^2/eV^2) * (L/km)/(E/GeV)


def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))


def rot(i,j,sine):
    x=np.eye(3)
    x[i,i]=x[j,j]=sqrt(1-sine**2)
    x[i,j],x[j,i]=sine,-sine
    return x


def angles(U):
    s13=abs(U[0,2])
    c13=sqrt(1-s13*s13)
    s12,s23=abs(U[0,1])/c13,abs(U[1,2])/c13
    c12,c23=sqrt(1-s12*s12),sqrt(1-s23*s23)
    frac=(abs(U[1,0])**2-s12**2*c23**2-c12**2*s23**2*s13**2)/(2*s12*c12*s23*c23*s13)
    quartet=(U[0,0]*U[1,1]*U[0,1].conj()*U[1,0].conj()).imag
    J0=s12*c12*s13*c13*c13*s23*c23
    return np.array([s12**2,s13**2,s23**2]),quartet/J0,frac


def radians_on_circle(a,ends):
    return any(ends[0]*pi<=a+2*k*pi<=ends[1]*pi for k in (-2,-1,0,1,2))


def calculate():
    hist=read("../1386/inputs.json")
    splits=read("../1370/inputs.json")
    alpha,err=hist["quark_sines_and_phase"][3]
    amp=hist["DayaBay_amplitude"][0]
    x13=(1-sqrt(1-amp))/2
    phases={"minus_alpha":-alpha,"plus_alpha":alpha,"pi_minus_alpha":pi-alpha,"pi_plus_alpha":pi+alpha}
    periods={"marginal68":[-.81,-.26],"NO_9973":[-1.38,.30],"IO_9973":[-.92,-.04]}
    memberships={p:{name:radians_on_circle(a,ends) for name,ends in periods.items()}
                 for p,a in phases.items()}
    def row_tests(order):
        par=hist[order]
        x12,x23=par["solar_sine_squared"][0],par["atmospheric_sine_squared"][0]
        s12,s13,s23=map(sqrt,[x12,x13,x23])
        c12,c13,c23=map(sqrt,[1-x12,1-x13,1-x23])
        da=splits[order]["solar"][0]
        db=splits[order]["atm_abs"][0] if order=="NO" else da-splits[order]["atm_abs"][0]
        phi=750*CONVERSION*np.array([0.,da,db])
        E=np.exp(-2j*phi)
        X=c13*s12*c12*c23*(E[1]-1)
        Y=c13*s13*s23*(E[2]-c12*c12-x12*E[1])
        Q=X.conjugate()*Y
        base=abs(X)**2+abs(Y)**2
        outputs={}
        for label,delta in phases.items():
            second=label in ("plus_alpha","pi_plus_alpha")
            flip=label in ("pi_minus_alpha","pi_plus_alpha")
            seed_e_T=rot(1,2,s23)@rot(0,2,-s13 if flip else s13)
            seed_n=rot(0,1,s12)
            if second:
                swap=np.eye(3)[[1,0,2]]
                seed_e_T=seed_e_T@swap.T
                seed_n=swap@seed_n
            U=seed_e_T@np.diag([np.exp(-1j*alpha),1,1])@seed_n
            xs,sd,cd=angles(U)
            np.testing.assert_allclose(xs,[x12,x13,x23])
            np.testing.assert_allclose([sd,cd],[sin(delta),cos(delta)])
            S=(U*E)@U.conj().T
            e_mu,mu_e=abs(S[1,0])**2,abs(S[0,1])**2
            odd,even=e_mu-mu_e,(mu_e+e_mu)/2-base
            np.testing.assert_allclose(odd,-4*Q.imag*sin(delta),atol=2e-14)
            np.testing.assert_allclose(even,2*Q.real*cos(delta),atol=2e-14)
            outputs[label]=[float(e_mu),float(mu_e),float(odd),float(even)]
        return {"CP_even_baseline":base,"phase_odd_multiplier":-4*Q.imag,
                "phase_even_multiplier":2*Q.real,"ordered_rates_difference_and_subtracted_mean":outputs}

    # Outward scope deliberately wider than the point example. Delta m^2_31 is used for both orderings.
    lx=(742.5,757.5)
    den=(.0021,.0029)
    small=(lx[0]*CONVERSION*.000060,lx[1]*CONVERSION*.000090)
    large=(lx[0]*CONVERSION*den[0],lx[1]*CONVERSION*den[1])
    phase_diff=(large[0]-small[1],large[1]+small[1])
    assert pi/2<phase_diff[0]<phase_diff[1]<pi
    sb=min(sin(v) for v in large)
    sd=min(sin(v) for v in phase_diff)
    cdelta=abs(cos(phase_diff[0]))
    p12=(.27,.35)
    p23=(.38,.65)
    p13=(.018,.025)
    j0=sqrt(min(v*(1-v) for v in p12)*min(v*(1-v) for v in p23)*p13[0])*(1-p13[1])
    lower_D=16*j0*sin(small[0])*sb*sd*sin(alpha-2*err)
    re_term=sb*cdelta-p12[1]*sin(small[1])
    assert re_term>0
    lower_E=8*j0*sin(small[0])*re_term*cos(alpha+2*err)
    return {"memberships_at_central_alpha":memberships,
            "normal_order_vacuum":row_tests("NO"),"inverted_order_vacuum":row_tests("IO"),
            "window_certificates":{"min_abs_difference":lower_D,"min_abs_subtracted_mean":lower_E,
                                   "absolute_phase21":list(small),"absolute_phase31":list(large),
                                   "absolute_phase32_envelope":list(phase_diff)}}


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    return bool(np.allclose(a,b,rtol=1e-10,atol=1e-14))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert same(result,read("independent_results.json"))
        print("1416 same vacuum, four physical orientations and two-channel propagation verified")
