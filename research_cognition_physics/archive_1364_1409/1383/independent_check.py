"""Independent 1383 check of tagged lepton normalizations and a finite-threshold control."""
import argparse
import json
from math import log, sqrt
from pathlib import Path

import numpy as np


def two_lepton_fit():
    values=np.array([.2155,.2167,.2153,.2152,.2132])
    errors=np.array([.0007,.0009,.0014,.0013,.0015])
    corr=np.array([[1,.28,.16,.07,.01],[.28,1,.08,.01,.18],[.16,.08,1,.04,.01],
                   [.07,.01,.04,1,.67],[.01,.18,.01,.67,1]])
    cov=corr*np.outer(errors,errors)
    whiten=np.linalg.cholesky(cov)
    design=np.array([[1.,0.],[0.,1.],[1.,0.],[1.,0.],[0.,1.]])
    x=np.linalg.solve(whiten,design)
    y=np.linalg.solve(whiten,values)
    pars=np.linalg.lstsq(x,y,rcond=None)[0]
    pcov=np.linalg.inv(x.T@x)
    cweights=np.linalg.solve(pcov,np.ones(2))
    cweights/=sum(cweights)
    assert min(cweights)>0
    total=pars@cweights
    assert abs(total-.21567820258)<1e-10
    ratio=pars[0]**2/pars[1]**2
    return {"separate_e_mu_product":pars.tolist(),"covariance":pcov.tolist(),
            "common_fit_flavor_weights":cweights.tolist(),"combined_product":total,
            "two_role_ratio_Ae_over_Amu":ratio}


def data_width(n):
    amax=.97373+n*.00031
    bmax=(.0351+n*.0045)/.80*(.083+n*.008)
    fmin=.9704-n*.0032
    row=[]
    for p,dp in [(.2157,.0006),(.2155,.0007),(.2167,.0009)]:
        kmax=(p+n*dp)/fmin
        row.append(amax**2+bmax**2+kmax**2)
    pmax=(1.2344+n*.0042)/(1.2352-n*.0002)
    bmax_e=bmax*max(1,sqrt(pmax))
    ke_high=(.2155+n*.0007)/fmin
    tau_low=1.0016-n*.0014
    qmax=sqrt(amax**2+bmax_e**2+ke_high**2)/tau_low
    assert (max(row)<1)==(n==1) and (qmax<1)==(n==1)
    return {"row_upper_average_e_mu":row,"q_upper_no_tau_contact":qmax,
            "true_b_required_by_beta_Kmu_without_b_semileptonic":
                sqrt(max(0.,1-amax**2-((.2167+n*.0009)/fmin)**2))}


def threshold_control():
    mass,mu=61.,105.6584
    x=(mass/mu)**2
    # t=q^2/m_mu^2, q is the sum of the two massless momenta. Set t=tmax*(1-z^2)
    # to make the square-root endpoint smooth.
    nodes,weights=np.polynomial.legendre.leggauss(48)
    z=(nodes+1)/2
    tmax=(1-sqrt(x))**2
    t=tmax*(1-z*z)
    rad=(tmax*z*z)*((1+sqrt(x))**2-t)
    phase=12*sum(weights/2*(2*tmax*z)*t*(1+x-t)*np.sqrt(rad))
    exact=1-8*x+8*x**3-x**4-12*x**2*log(x)
    assert abs(phase-exact)<5e-15
    a=.974
    s=.001
    mu_den=1-s+s*phase
    visible=(1-s*a*a)/mu_den
    assert visible<1 and phase>1-a*a
    # Conditional absolute experimental ceiling from PIENU, conservatively using the
    # least restrictive published envelope 1e-7 over its 60--135 MeV search range.
    cap=1e-7
    largest=cap*(phase+a*a-1)/(1-cap*(1-phase))
    assert 0<largest<1e-8
    return {"mass_MeV":mass,"three_body_fraction":phase,"true_ud":a,
            "trial_e_mixing":s,"apparent_row_at_trial":visible,
            "missing_probability_envelope_from_peak_search_at_this_ud":largest}


def run():
    first,second=data_width(1),data_width(2)
    return {"table7":two_lepton_fit(),"one_width":first,"two_widths":second,
            "PIENU_over_SM":1.2344/1.2352,"finite_mass_control":threshold_control()}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".10g")) if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--check",action="store_true")
    out=Path(__file__).with_name("independent_results.json")
    result=rounded(run())
    if ap.parse_args().write:
        with out.open("x",encoding="utf-8") as fp:
            json.dump(result,fp,indent=2,ensure_ascii=False)
            fp.write("\n")
    else:
        assert result==json.loads(out.read_text(encoding="utf-8"))
    print("1383 independent flavor-resolved fit, tau comparison and 61 MeV threshold verified.")
