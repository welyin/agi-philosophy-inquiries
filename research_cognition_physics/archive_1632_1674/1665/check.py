"""Compare the unavoidable colour thresholds with the scale reach of the CMS jet fit."""
import argparse
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
N = 4
NORMALIZATION = 6*math.pi


def read(path):
    return json.loads((HERE/path).read_bytes())


def interval_witness(q1, q2, lower, mean_upper):
    entry = math.log(q1/lower)
    width = math.log(q2/q1)
    budget = N*math.log(mean_upper/lower)
    choices = [(min(r*width, budget-r*entry), r)
               for r in range(N+1) if r*entry <= budget]
    saved, count = max(choices)
    # Realize the best saving using equal paid masses; any remaining budget
    # is optional because the scalar constraint is an upper bound on the product.
    logs = [entry+saved/count]*count if count else []
    logs.extend([0.]*(N-count))
    masses = [lower*math.exp(x) for x in logs]
    total = sum(max(math.log(q2/m),0)-max(math.log(q1/m),0) for m in masses)
    assert math.isclose(total,N*width-saved,abs_tol=1e-12)
    assert sum(logs)<=budget+1e-12
    return {'minimal_increment':total/NORMALIZATION,'realizing_doublets_GeV':masses,
            'simple_upper':N*width/NORMALIZATION}


def compute():
    data=read('inputs.json')
    earlier=read(data['bounds'])
    # Round up for a conservative bound, including the precision of the source JSON.
    cap = math.ceil(earlier['four_copies_geometric_mean_upper_GeV']*1e5)/1e5
    fit=data['CMS_Q_gev_alpha_minus_plus']
    lower=data['MZ_GeV']
    rows=[]
    for q, alpha, down, up in fit:
        threshold=N*max(math.log(q/cap),0)/NORMALIZATION
        inverse_decrease=1/alpha-1/(alpha+up)
        inverse_increase=1/(alpha-down)-1/alpha
        rows.append({'Q_GeV':q,'inverse_coupling_decrease_minimum':threshold,
                     'reported_inverse_width_toward_larger_alpha':inverse_decrease,
                     'reported_inverse_width_toward_smaller_alpha':inverse_increase,
                     'ratio_to_width_toward_larger_alpha':threshold/inverse_decrease})
    q1,q2=data['selected_pair']
    pair=interval_witness(q1,q2,lower,cap)
    # At one actual charged-flavour point, exhibit spectra with the same X,Y
    # and different running inside the interval; the unification scale is untouched.
    _,x,y=earlier['central_gain_and_log_ratios']
    central_cap=1e16*math.exp(-x/N)
    central_interval=interval_witness(q1,q2,lower,central_cap)
    def complete(masses):
        geometric_mean=math.exp(sum(math.log(m) for m in masses)/N)
        m_delta=geometric_mean*math.exp(x/N)
        m_sextet=geometric_mean*math.exp(y/N)
        assert q2<m_delta<=1e16*(1+1e-12) and q2<m_sextet<=1e16
        return {'doublets_GeV':masses,'common_weak_triplet_GeV':m_delta,
                'common_colour_sextet_GeV':m_sextet}
    central_examples=[complete(central_interval['realizing_doublets_GeV']),
                      complete([240]*N)]
    assert pair['minimal_increment']>0
    return {'four_copy_geometric_mean_upper_used_GeV':cap,
            'cumulative_and_published_widths':rows,'interbin_envelope':pair,
            'central_interbin_envelope':central_interval,'two_central_spectra':central_examples}


def tidy(x):
    if isinstance(x,list):
        return [tidy(y) for y in x]
    if isinstance(x,dict):
        return {k:tidy(v) for k,v in x.items()}
    return float(f'{x:.11g}') if isinstance(x,float) else x


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read('results.json')
        print('PASS: mass-product bound, full two-bin optimization, two spectra, and CMS uncertainty conversion')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
