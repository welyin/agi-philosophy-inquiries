"""Conditional arithmetic and Wick algebra diagnostics, never field constants."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import numpy as np

def calculate():
    # Entirely specified conditional inputs. They are NOT estimates for the QFT.
    g,alpha,v0,v1,delta=F(1,8),F(1,64),F(9,2),F(5,4),F(1,16)
    y=F(6) # s^2: both sufficient bounds evaluate below or at six.
    assert y*y*delta*delta >= 2*g*g*v0
    assert y*delta*delta >= 2*g*alpha*v1
    variance=g*g*v0/(y*y)+g*alpha*v1/y
    assert variance<=delta*delta
    # A=(a+a^*)^2/2 in a finite Fock matrix. The vacuum polynomial reaches
    # at most |4>, so dimensions >=5 reproduce these moments exactly.
    rows=[]
    for n in (6,10):
        a=np.diag(np.sqrt(np.arange(1,n,dtype=float)),1)
        x=(a+a.T)/np.sqrt(2)
        q=x@x
        z=np.zeros(n);z[0]=1
        for eta in (0.,0.5,2.):
            shifted=q+2*eta*x+eta*eta*np.eye(n)
            mean=float(z@shifted@z)
            var=float(z@(shifted@shifted)@z-mean*mean)
            # Gaussian vacuum: Var(X^2)=1/2 and Var(2X)=2.
            predicted=0.5+2*eta*eta
            assert abs(var-predicted)<1e-12
            rows.append({'dimension':n,'eta':eta,'variance':var,'predicted':predicted})
    return {'status':'PASS_CONDITIONAL_ARITHMETIC_AND_ALGEBRA_ONLY',
        'actual_field_noise_constants_computed':False,
        'actual_geometry_or_instrument_error_certified':False,
        'conditional_inputs':{k:str(v) for k,v in locals().copy().items()
            if k in ('g','alpha','v0','v1','delta','y')},
        'conditional_variance':str(variance),'budget':str(delta*delta),
        'finite_polynomial_diagnostics':rows}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args(); path=Path(__file__).with_name('results.json')
    result=calculate()
    if args.write:
        assert not path.exists(),'Refuse overwrite'
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(path.read_text(encoding='utf-8'))
        assert saved==result,'Saved result mismatch'
    print(json.dumps({'status':result['status'],'mode':'first_write' if args.write else 'read_only_match'}))

if __name__=='__main__':main()
