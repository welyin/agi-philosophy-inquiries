"""Bound charged pion-neutral proton modes using the actual charged flavor flags."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
M_P=0.938272088
M_PI=0.1349768
LEPTONS=[0.00051099895, 0.105658376]
# Table VIII continuum, mixed chirality, pi+, in GeV^2; the pi0 amplitude is /sqrt(2).
FORM=[(-.159,.169),(-.161,.170)]
ERR0=[(.015,.020,.025),(.015,.020,.026)]
ERR1=[(.014,.018,.029),(.014,.018,.028)]


def kinematic(j, w0, w1):
    m=LEPTONS[j]
    energy=(M_P**2-M_PI**2+m**2)/(2*M_P)
    momentum=math.sqrt(energy**2-m**2)
    return sorted(energy*momentum*(1-s*m/energy)*(w0-s*m/M_P*w1)**2
                  /(16*math.pi*M_P) for s in [-1,1])


def response(j, envelope=False):
    w0,w1=FORM[j]
    if not envelope:
        return kinematic(j,w0,w1)
    d0,d1=math.sqrt(sum(t*t for t in ERR0[j])),math.sqrt(sum(t*t for t in ERR1[j]))
    vals=[value for x in [-d0,d0] for y in [-d1,d1] for value in kinematic(j,w0+x,w1+y)]
    assert abs(w0)-d0 > LEPTONS[j]/M_P*(abs(w1)+d1)
    return min(vals),max(vals)


def extrema_ratio(a,b):
    # Quotient of quadratics in r >= 0, including r=infinity by continuity.
    def coeff(x,y):
        return np.array([abs(x)**2,2*(x.conjugate()*y).real,abs(y)**2])
    n=coeff(a[1],b[1]);d=coeff(a[0],b[0])
    assert min(d[0],d[2])>0 and (d[1]>=0 or 4*d[0]*d[2]>d[1]**2)
    polynomial=[n[2]*d[1]-n[1]*d[2],2*(n[2]*d[0]-n[0]*d[2]),n[1]*d[0]-n[0]*d[1]]
    rs=[0.,1.]+[z.real for z in np.roots(polynomial) if abs(z.imag)<1e-10 and z.real>0]
    values=[(abs((a[1]+r*b[1])/(a[0]+r*b[0]))**2,r) for r in rs]
    values.append((abs(b[1]/b[0])**2,None))
    return min(values,key=lambda p:p[0]),max(values,key=lambda p:p[0])


def windows(ratios, envelope=False):
    e=response(0,envelope);m=response(1,envelope)
    return [float(m[0]/e[1]*min(ratios)),float(m[1]/e[0]*max(ratios))]


def witness(x,data,old):
    _,_,_,_,_,E,_,D=old['matrices'](x,data)
    # Physical Takagi columns can instead be conjugated and phase-normalized. Only
    # the absolute value of each row's a+r*b enters; phases/conjugation drop out.
    u=np.array([1,0,0])
    d=D[:,0]
    a=E.conj().T@d
    b=(u@d)*(E.conj().T@u)
    a2=a+b
    assert np.isclose(np.linalg.norm(a),1) and np.isclose(np.vdot(a,b),abs(u@d)**2)
    pure=abs(a2[1]/a2[0])**2
    other=abs(a[1]/a[0])**2
    rlo,rhi=extrema_ratio(a,b)
    return dict(cc_squared_mu_over_e=float(pure),pentabar_squared_mu_over_e_single_X=float(other),
                second_chiral_ratio_extrema_and_arg=[rlo,rhi],
                single_X=windows([pure,other]),equal_X_Xprime=windows([pure]),
                single_X_equal_chiral_transport=windows([(abs(a2[1])**2+abs(a[1])**2)/(abs(a2[0])**2+abs(a[0])**2)]),
                Xprime_dominant_equal_transport_limit=windows([abs(b[1]/b[0])**2]),
                any_X_Xprime=windows([pure,rlo[0],rhi[0]]),
                wide_lattice_any_X_Xprime=windows([pure,rlo[0],rhi[0]],True))


def compute():
    folder=HERE.parent/'1799'
    old=runpy.run_path(str(folder/'check.py'))
    data=old['backend']()
    cfg=json.loads((folder/'inputs.json').read_bytes())
    points=cfg['mixed_type_I_II'] | cfg['fixed_delta_minus70_mass_readout']
    cases={name:witness(np.array(x),data,old) for name,x in points.items()}
    bounds=[2.4e34,1.6e34]
    all_upper=max(s['wide_lattice_any_X_Xprime'][1] for s in cases.values())
    tolerance=math.sqrt((bounds[0]/bounds[1])/all_upper)
    return dict(pion_phasespace_masses_GeV=[M_P,M_PI]+LEPTONS,
                pi_plus_form_factors=FORM,form_factor_quadrature_inputs=[ERR0,ERR1],
                response_eig_central_GeV5=[response(0),response(1)],
                response_eig_1_marginal_rect_GeV5=[response(0,True),response(1,True)],
                sk_partial_lifetime_bounds_years=bounds,
                all_five_wide_upper_mu_over_e=all_upper,
                implied_mu_partial_lifetime_over_years=bounds[0]/all_upper,
                max_extra_relative_amplitude_to_improve_direct_mu_limit=(tolerance-1)/(tolerance+1),
                cases=cases)


def verify_values(observed, saved):
    if isinstance(saved, dict):
        assert isinstance(observed, dict) and observed.keys() == saved.keys()
        for key, value in saved.items():
            verify_values(observed[key], value)
    elif isinstance(saved, list):
        assert isinstance(observed, list) and len(observed) == len(saved)
        for result, value in zip(observed, saved, strict=True):
            verify_values(result, value)
    elif isinstance(saved, float):
        assert math.isclose(observed, saved, rel_tol=2e-11, abs_tol=2e-15), (observed, saved)
    else:
        assert observed == saved


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    raw=json.dumps(compute(),indent=2,allow_nan=False)+'\n'
    file=HERE/'independent_results.json'
    if args.check:
        reference = file.read_bytes()
        verify_values(json.loads(raw), json.loads(reference))
    else:
        reference = raw.encode()
        file.write_bytes(reference)
    print('1804 independent',hashlib.sha256(reference).hexdigest()[:12])
