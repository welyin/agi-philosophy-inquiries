"""Finite York reference diagnostic. Default replay is read-only."""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, math
import numpy as np

HERE = Path(__file__).resolve().parent

def calculate():
    tlo, thi = F(99,500), F(101,500)
    lower_r2, upper_r2 = F(256,125)/thi**3, F(256,125)/tlo**3
    assert F(64) < lower_r2 <= upper_r2 < F(289)
    c = F(1)-F(9,200)*F(22,7)
    assert c == F(601,700)
    plo, phi = c/(1+c), 1/(1+c)
    assert 1-2*plo-plo**2 < 0
    lower = phi*(1-phi)/(1+phi)
    assert lower == F(420700,2603301) and lower > F(4,25)
    assert F(4,25)-2*F(1,100) == F(7,50)
    # Constant-free exponent difference: avoid subtracting two ~10^4 terms.
    beta = 320*math.pi/math.sqrt(5)
    root0 = math.sqrt(.2)
    def relf(m):
        d = m-40
        # E(m)-E(40)=2d/(sqrt(1-m/50)+sqrt(.2)).
        return d*(4*math.pi*(m+40)-2*beta/(np.sqrt(1-m/50)+root0))
    rows=[]
    for n in [1,8,32,128,512]:
        low=39.9+(np.arange(n)+.5)*(.1/n)
        high=40+(np.arange(n)+.5)*(.1/n)
        gl,gh=np.exp(relf(low)),np.exp(relf(high))
        assert gl.min()>float(c) and gh.min()>float(c)
        assert gl.max()<=1+1e-12 and gh.max()<=1+1e-12
        l,h=gl.sum(),gh.sum()
        p0=h/(l+h); p1=2*h/(l+2*h)
        assert float(plo)<p0<float(phi) and p1-p0>float(lower)
        assert abs(p1-2*p0/(1+p0))<1e-14
        # Finite quantum diagnostic: weights are explicit and not E-only Gibbs.
        diag0=np.r_[gl,gh];diag0/=diag0.sum()
        diag1=np.r_[gl,gh,gh];diag1/=diag1.sum()
        assert abs(diag0[n:].sum()-p0)<1e-14
        assert abs(diag1[n:].sum()-p1)<1e-14
        rows.append({'N_per_half':n,'p0':float(p0),'p1':float(p1),
                     'gap':float(p1-p0),'dimension0':2*n,'dimension1':3*n})
    # Algebraic midpoint and normalization bounds, not numeric quadrature proof.
    eta_n = F(99,14000)
    denom_lower = 2*F(1,10)*c
    tight = eta_n/denom_lower
    stated = F(297,4808)
    assert tight == F(99,2404) and tight < stated
    return {'status':'PASS_BOUNDED_REFERENCE_MEASURE_SCREEN',
      'scope':'Formal candidate probabilities, not an actual instrument, TR implementation, or GR derivation.',
      'exact':{'r_squared_lower':str(lower_r2),'r_squared_upper':str(upper_r2),
        'g_lower':str(c),'p0_lower':str(plo),'p0_upper':str(phi),
        'gap_lower':str(lower),'simple_gap_lower':'4/25',
        'conditional_gap_two_1pct_errors':'7/50',
        'midpoint_probability_error_tight_times_N':str(tight),
        'midpoint_probability_error_stated_times_N':str(stated),
        'reweighted_probability_error_stated_times_N':'297/2404'},
      'diagnostics_not_proof':rows,'new_formal_rounds':0,'new_adopted_axioms':0}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); result=calculate(); path=HERE/'results.json'
    if args.write:
        assert not path.exists(),'Do not overwrite frozen results'
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('WROTE_FIRST_RESULTS')
    else:
        prior=json.loads(path.read_text(encoding='utf-8'))
        assert prior['exact']==result['exact'] and prior['status']==result['status']
        assert len(prior['diagnostics_not_proof'])==len(result['diagnostics_not_proof'])
        for a,b in zip(prior['diagnostics_not_proof'],result['diagnostics_not_proof']):
            assert a['N_per_half']==b['N_per_half']
            assert all(abs(a[k]-b[k])<1e-12 for k in ['p0','p1','gap'])
        print('READ_ONLY_MATCH_BOUNDED_REFERENCE_SCREEN')
    print(json.dumps({'status':result['status'],'gap_lower':result['exact']['gap_lower']},ensure_ascii=False))

if __name__=='__main__':main()
