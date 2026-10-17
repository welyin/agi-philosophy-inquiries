"""Exact finite-window certificates for a conditional York thermal boundary.

No quantum-gravity path integral, device, or microscopic reservoir is simulated.
First --write creates the result exclusively; default rechecks without mutation.
"""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import argparse
import json

DEN=10**50

def interval(x):
    return x if isinstance(x,tuple) else (F(x),F(x))

def add(a,b):
    a,b=interval(a),interval(b)
    return a[0]+b[0],a[1]+b[1]

def neg(a):
    a=interval(a)
    return -a[1],-a[0]

def sub(a,b):return add(a,neg(b))

def mul(a,b):
    a,b=interval(a),interval(b)
    values=[x*y for x in a for y in b]
    return min(values),max(values)

def sqrt_interval(v):
    v=F(v); k=isqrt((v.numerator*DEN*DEN)//v.denominator)
    lo,hi=F(k,DEN),F(k+1,DEN)
    assert lo*lo<=v<=hi*hi
    return lo,hi

def out(a):
    a=interval(a)
    return {'exact':[str(a[0]),str(a[1])],
            'display_only':[float(a[0]),float(a[1])]}

def calculate():
    # Ebar=G E/R, bbar=beta/(4 pi R), Sbar=G S/(pi R^2).
    x0,x1=F(9,10),F(91,100)
    assert F(8,9)<x0<x1<1  # negative classical action and positive heat capacity
    q0,q1=sqrt_interval(1-x0),sqrt_interval(1-x1)
    e0,e1=sub(1,q0),sub(1,q1)
    b0,b1=mul(x0,q0),mul(x1,q1)
    de,db=sub(e1,e0),sub(b1,b0)
    ds=x1*x1-x0*x0
    residual=sub(ds,mul(4,mul(b0,de)))
    dem=(x1-x0)/2
    wrong=sub(ds,mul(4,mul(b0,dem)))
    magnitude=neg(residual)
    de2=mul(de,de)
    lower=mul(2*(3*x0-2),de2)
    upper=mul(2*(3*x1-2),de2)
    assert de[0]>0 and db[1]<0
    assert magnitude[0]>lower[1] and magnitude[1]<upper[0]
    assert magnitude[1]<F(4,10000)
    assert wrong[0]>F(12,1000)
    gap=sub(wrong,residual)
    assert gap[0]>F(12,1000)
    # Direct boundary integral vs beta E - area/4G at both endpoints.
    action_checks=[]
    for x,q,e,b in [(x0,q0,e0,b0),(x1,q1,e1,b1)]:
        action_ghy=mul(4*x,add(sub(q,1),3*x/4))
        action_thermo=sub(mul(4,mul(b,e)),x*x)
        err=sub(action_ghy,action_thermo)
        assert err[0]<=0<=err[1] and err[1]-err[0]<F(1,10**45)
        assert action_ghy[1]<0
        heat_capacity=4*x*x*(1-x)/(3*x-2)  # G C/(pi R^2)
        assert heat_capacity>0
        action_checks.append({'x':str(x),'Ibar':out(action_ghy),
          'same_action_identity_residual':out(err),'Cbar':str(heat_capacity)})
    return {'status':'PASS_CONDITIONAL_YORK_FINITE_IDENTITIES',
      'scope':'leading stationary EH+GHY canonical saddle; not complete quantum/process qualification',
      'normalizations':{'Ebar':'G*E/R','bbar':'beta/(4*pi*R)','Sbar':'G*S/(pi*R^2)'},
      'x_interval':[str(x0),str(x1)],'E_change':out(de),'beta_change':out(db),'S_change':str(ds),
      'fixed_initial_beta_residual':out(residual),'quadratic_magnitude_lower':out(lower),
      'quadratic_magnitude_upper':out(upper),'wrong_ADM_residual':out(wrong),
      'wrong_vs_correct_gap':out(gap),'action_checks':action_checks,
      'physical_loop_error_certified':False,'actual_reservoir_source_certified':False,
      'new_adopted_axioms':0,'new_numbered_round':False}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    dest=Path(__file__).with_suffix('.json');result=calculate()
    if args.write:
        with dest.open('x',encoding='utf-8') as stream:json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
        action='CREATED_ONCE'
    else:
        assert json.loads(dest.read_text('utf-8'))==result
        action='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'action':action,
      'residual':result['fixed_initial_beta_residual']['display_only'],
      'wrong_ADM':result['wrong_ADM_residual']['display_only'],'physical_completion':False}))

if __name__=='__main__':main()
