"""Higgs production and decay read out the same scalar mass derivative."""
import argparse
import importlib.util
import itertools
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def loop_inputs(cfg):
    spec=importlib.util.spec_from_file_location('old_loops', HERE/cfg['reference_amplitudes'])
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    old=module.read('inputs.json')
    w,_,fermions=module.common_inputs(old)
    mh=cfg['higgs_GeV']
    def form(m, spin):
        return module.loop_functions(mh,m)[spin]
    gluon=sum(form(m,1) for m,nc,_ in fermions if nc==3)
    photon=form(w,2)+sum(nc*q*q*form(m,1) for m,nc,q in fermions)
    return gluon,photon,lambda m:float(form(m,0)),old['v_MZ_GeV']


def real_shift_intervals(amplitude, coefficient, rates):
    # Real x reaches points on a line through the complex amplitude.
    radii=[r*abs(amplitude)**2-amplitude.imag**2 for r in rates]
    if radii[1]<0:
        return []
    middle=-amplitude.real/coefficient
    inner,outer=[math.sqrt(max(0,d))/coefficient for d in radii]
    if inner==0:
        return [[middle-outer,middle+outer]]
    return [[middle-outer,middle-inner],[middle+inner,middle+outer]]


def overlap(left,right):
    return [[max(a,c), min(b,d)] for (a,b),(c,d) in itertools.product(left,right)
            if max(a,c)<=min(b,d)]


def window(row,width):
    x,lower,upper=row
    return [x-width*lower,x+width*upper]


def compute():
    cfg=read('inputs.json')
    g,photon,scalar,v=loop_inputs(cfg)
    production=window(cfg['VBF_over_gluon'],cfg['widths'])
    decay=window(cfg['diphoton_over_ZZ'],cfg['widths'])
    pure_rates=[(1/r-b)/(1-b) for r,b in itertools.product(production,cfg['bb_fraction'])]
    prod_x=real_shift_intervals(g,1,[min(pure_rates),max(pure_rates)])
    photon_x=real_shift_intervals(photon,5/6,decay)
    both=overlap(prod_x,photon_x)
    assert len(prod_x)==len(photon_x)==2 and len(both)==1

    mass_cap=math.ceil(read(cfg['mass_box'])['four_copies_geometric_mean_upper_GeV']*1e4)/1e4
    full=[8*scalar(mass_cap),8*scalar(cfg['mass_floor_GeV'])]
    assert full[0]>both[0][1] and full[1]<-photon.real/(5/6)
    def gamma_rate(x):
        return abs(1+5*x/(6*photon))**2
    def gluon_rate(x):
        return abs(1+x/g)**2
    bb=cfg['bb_fraction'][1]
    holdout=overlap(prod_x,[[0.,full[1]]])[0]
    predicted_photons=[gamma_rate(holdout[1]),gamma_rate(holdout[0])]
    exact_fraction=cfg['example_signed_fraction']
    mass=cfg['example_doublet_mass_GeV']
    f=[exact_fraction,-exact_fraction,0,0]
    lam=[2*m_i*mass**2/v**2 for m_i in f]
    x=2*sum(a*scalar(mass) for a in f)
    linear=sum(f)
    curvature=sum(a-2*a*a for a in f)
    assert abs(x)<1e-15 and curvature<0 and gamma_rate(x)==gluon_rate(x)==1
    return {'SM_gluon_and_photon':[[g.real,g.imag],[photon.real,photon.imag]],
            'two_width_windows':[production,decay],
            'gluon_amplitude_shift_sets':{'production':prod_x,'decay':photon_x,'joint':both},
            'sum_positive_f_upper':1.5*both[0][1],
            'common_fraction_upper_in_mass_box':both[0][1]/full[0],
            'predicted_gamma_over_ZZ_given_positive_production':predicted_photons,
            'all_four_pure_Higgs':{'shift_envelope':full,
                'normalized_VBF_over_grouped_ggF_at_most':1/(bb+(1-bb)*gluon_rate(full[0])),
                'normalized_gamma_over_ZZ_at_most':gamma_rate(full[0]),
                'max_gluon_rate_transfer_to_enter_production_box':
                    max(pure_rates)/gluon_rate(full[0])},
            'signed_example':{'mass_GeV':mass,'f':f,'lambda':lam,'one_higgs_shift':x,
                'log_mass_sum_linear':linear,'log_mass_sum_second_derivative':curvature}}


def tidy(x):
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [tidy(v) for v in x]
    return float(f'{x:.11g}') if isinstance(x,float) else x


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read('results.json')
        print('PASS: finite loops, two amplitude branches, Higgs mass bound, and signed competitor')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
