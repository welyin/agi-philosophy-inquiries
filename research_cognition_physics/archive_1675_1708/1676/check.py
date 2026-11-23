"""A three-channel one-loop identity and a finite-mass bound on two width ratios."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def reference(cfg):
    loops = runpy.run_path(str(HERE/cfg['loop_convention']))
    ref = loops['prior']
    mw,mz,fermions = ref.common_inputs(read('../1373/inputs.json'))
    mh = cfg['m_h_GeV']
    cw = mw/mz
    sw = math.sqrt(1-cw*cw)
    tw2 = (sw/cw)**2
    i1,i2 = loops['integrals'](mh,mz,mw)
    aw = (2*mw/mh)**2
    az = cw/sw*(4*(3-tw2)*i2+((1+2/aw)*tw2-5-2/aw)*i1)
    for m,nc,q in fermions:
        a,b = loops['integrals'](mh,mz,m)
        az += 2*nc*q*(math.copysign(.5,q)-2*q*sw**2)/(sw*cw)*(a-b)
    def forms(m):
        return ref.loop_functions(mh,m)[0],2*loops['integrals'](mh,mz,m)[0]
    return az,sw,cw,forms


def interval(row,n):
    return [row[0]-n*row[1],row[0]+n*row[2]]


def compute():
    cfg = read('inputs.json')
    data = read(cfg['production_input'])
    xcap = max(row[1] for row in data['gluon_amplitude_shift_sets']['production'])+1e-10
    az,sw,cw,forms = reference(cfg)
    photon = complex(*data['SM_gluon_and_photon'][1])
    charges = np.array([2/3,-1/3])
    isospin = np.array([.5,-.5])
    coeff = 3*charges**2
    kz = 3*charges*(isospin-sw**2*charges)/(sw*cw)
    np.testing.assert_allclose(kz,((.5-sw**2)*coeff+1/3)/(sw*cw))
    a,b = forms(cfg['scalar_mass_floor_GeV'])
    rmax = b/a+1e-10
    xmax = kz.max()*rmax*xcap
    ymax = coeff.max()*xcap
    assert az.real-xmax>abs(az.imag) and -photon.real-ymax>abs(photon.imag)
    def normalized(a,shift):
        return abs(1+shift/a)**2
    def slope(x,b):
        return x/(x*x+b*b)
    # Monotonicity on the entire positive simplex, including all mass ratios r.
    inc = coeff[0]*slope(-photon.real,photon.imag)-kz[0]*rmax*slope(az.real-xmax,az.imag)
    dec = kz[1]*slope(az.real,az.imag)-coeff[1]*slope(-photon.real-ymax,photon.imag)
    assert min(inc,dec)>0
    def double(i,r):
        return normalized(az,-kz[i]*r*xcap)/normalized(photon,coeff[i]*xcap)
    lo,hi = double(1,rmax),double(0,1.)
    sm,unc = cfg['cms_sm_ratio']
    absolute = [(sm-unc)*lo,(sm+unc)*hi]
    zlo = normalized(az,-xmax)
    cms = cfg['cms_Zgamma_over_gammagamma_absolute']
    atlas = cfg['atlas_Zgamma_over_ZZ_normalized']
    assert absolute[1]<interval(cms,1)[0] and interval(cms,2)[0]<absolute[0]
    assert 1<interval(atlas,1)[0] and interval(atlas,2)[0]<zlo
    masses = cfg['sign_example_masses_GeV']
    pairs = [forms(m) for m in masses]
    f = cfg['positive_first_fraction']
    g = f*pairs[0][0]/pairs[1][0]
    def example(fractions):
        x = sum(fi*a for fi,(a,_) in zip(fractions,pairs))
        z = -sum(kz)*sum(fi*b for fi,(_,b) in zip(fractions,pairs))
        return [2*x,sum(coeff)*x,z]
    a,b = example([f,0,0,0]),example([0,g,0,0])
    np.testing.assert_allclose(a[:2],b[:2],rtol=1e-13)
    assert abs(a[2]-b[2])>.004
    return {'finite_form_ratio_max':rmax,'positive_gg_shift_max':xcap,
            'weak_angle_sine_squared':sw**2,'positive_Zgamma_ZZ_normalized':[zlo,1.],
            'derivative_sign_margins': [inc,dec],
            'positive_Zgamma_gammagamma_normalized_outer':[lo,hi],
            'absolute_CMS_ratio_with_reference_band':absolute,
            'CMS_one_and_two_widths':[interval(cms,1),interval(cms,2)],
            'ATLAS_one_and_two_widths':[interval(atlas,1),interval(atlas,2)],
            'example_gg_gamma_Z_amplitudes':[a,b],'second_portal':g}


def rounded(o):
    if isinstance(o,dict):
        return {k:rounded(v) for k,v in o.items()}
    if isinstance(o,list):
        return [rounded(v) for v in o]
    return round(float(o),10)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    out=rounded(compute())
    if args.check:
        assert out==read('results.json')
        print('PASS: finite scalar loops, two width-ratio bounds and positive competitors')
    else:
        with (HERE/'results.json').open('x',encoding='utf8') as f:
            f.write(json.dumps(out,indent=2)+'\n')
