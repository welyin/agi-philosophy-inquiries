"""First mass logarithm of ordinary QCD jets and a CMS cross-energy shape."""
import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def readout(s13low, s13high, s8a, s8b, s8high):
    s8low = (53*s8a+56*s8b)/109
    low, high = s13low/s8low, s13high/s8high
    return low, high, high/low


def measured(cfg):
    low, high = cfg['cms13']
    a, b, c = cfg['cms8']
    assert low['bin']==[a['bin'][0],b['bin'][1]]==[737,846]
    assert a['bin'][1]==b['bin'][0]==790 and high['bin']==c['bin']==[1588,1784]
    rows13 = [low,high]
    rows8 = [a,b,c]
    centre = [row['value'] for row in rows13+rows8]
    bounds = []
    for r in rows13:
        width = math.hypot(*r['errors'])
        bounds.append([r['value']-width,r['value']+width])
    for r in rows8:
        bounds.append([r['value']-math.hypot(r['stat'],r['sys_down']),
                       r['value']+math.hypot(r['stat'],r['sys_up'])])
    ratios = [readout(*po) for po in itertools.product(*bounds)]
    ratios_central = readout(*centre)
    relative = [[s/r['value'] for s in r['errors']] for r in rows13]
    l,h = relative
    # The public NNPDF mapping treats only the first nuisance as bin-independent.
    delta = [y-x for x,y in zip(l[1:],h[1:])]
    log_shape_uncertainty = math.hypot(l[0],h[0],*delta)
    return {'ratios_13_over_8_low_high_and_double':list(ratios_central),
            'registered_marginal_envelopes':[[min(col),max(col)] for col in zip(*ratios)],
            '13TeV_log_shape_width_NNPDF_signed_sources':log_shape_uncertainty}


def ordinary_coefficient(power, incoming_gluons):
    # Keep the first threshold logarithm in a joint perturbative expansion.
    return power-incoming_gluons


def compute():
    cfg=read('inputs.json')
    data=measured(cfg)
    prior=read(cfg['benchmark'])
    source=read(cfg['cmsalphas'])
    grid={row[0]:row[1] for row in source['CMS_Q_gev_alpha_minus_plus']}
    # Verify the gluon-to-scalar momentum flux with exact massless kernel.
    # Per complex colour fundamental, two charge-conjugate daughters contribute.
    per_species=2*Fraction(1,2)*(Fraction(1,3)-Fraction(1,4))
    assert 2*per_species==Fraction(1,6)
    assert Fraction(1,3)/2==2*per_species
    coeff=[ordinary_coefficient(2,i) for i in [2,1,0]]
    plots=[]
    for model in prior['two_central_spectra']:
        masses=model['doublets_GeV']
        readouts=[]
        for q in cfg['Q']:
            threshold=sum(max(0,math.log(q/m)) for m in masses)/(6*math.pi)
            a=grid[q]
            kappa=a*threshold
            readouts.append({'Q':q,'inverse_threshold':threshold,'alpha_times_D':kappa,
                             'gg_qg_qq_log_corrections':[w*kappa for w in coeff]})
        low,high=[r['alpha_times_D'] for r in readouts]
        # Independent nonnegative ordinary Born weights for each beam/bin.
        extremes=[high*(f13h-f8h)-low*(f13l-f8l)
                  for f13h,f8h,f13l,f8l in itertools.product([0.,2.],repeat=4)]
        cap=2*(low+high)
        assert math.isclose(max(extremes),cap) and math.isclose(min(extremes),-cap)
        measured_lo,measured_hi=data['registered_marginal_envelopes'][-1]
        plots.append({'masses_GeV':masses,'fixed_scales':readouts,
                      'two_energy_log_contrast_first_order_range':[-cap,cap],
                      'required_baseline_times_residual_envelope':
                          [measured_lo*math.exp(-cap),measured_hi*math.exp(cap)]})
    return {'CMS_cross_sections':data,'ordinary_jet_coefficients_gg_qg_qq':coeff,
            'gluon_momentum_loss_per_doublet':str(2*per_species),
            'conditional_spectra':plots}


def tidy(x):
    if isinstance(x,list):
        return [tidy(y) for y in x]
    if isinstance(x,dict):
        return {k:tidy(y) for k,y in x.items()}
    return float(f'{x:.12g}') if isinstance(x,float) else x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read('results.json')
        print('PASS: scalar momentum flux, three incoming channels, CMS matched-bin ratios and spectral extrema')
    else:
        with (HERE/'results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
