"""Finite scalar entropy-meter screen; default read-only replay, no device claim."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def exact_checks():
    mi, mf, d, tau = F(10), F(101, 10), F(1, 10), F(1000)
    nu = 15*d/(8*tau)
    lo, hi, zlo, zhi = F(49, 500), F(51, 500), F(1, 5), F(9, 5)
    z0 = [F(1), F(10**6, 101**3)]
    zdiff = 3*d/(4*lo)
    assert F(1)/mf-nu > lo and F(1)/mi+nu < hi
    assert min(z0)-zdiff > zlo and max(z0)+zdiff < zhi
    stable_bound = 4*zhi*(F(1, 10)/lo)**3
    assert stable_bound < 8 < mi**2
    width = 2*nu*d/lo**2
    assert width == F(75, 19208) and width < F(1, 250)
    work = d/4 + zhi*(d/(mi*mf)+2*nu)
    assert work < F(7, 250) < 1
    assert (mf**2-mi**2)/2 == F(201, 200)
    # At fixed total energy E=1, signs of d(Tgas-TH)/dM.
    signs = []
    for m in [F(3, 5), F(9, 10)]:
        th = 1/m
        derivative = th*(1/m-1/(4*(1-m)))
        signs.append({'M_over_E':str(m), 'derivative':str(derivative),
                      'stable':derivative<0})
    assert [r['stable'] for r in signs] == [False, True]
    # Same-energy non-grey radiation screen, not a quantum thermal state.
    moments = [F(6), F(1, 16)*F(6)/(F(1, 2)**4)]
    absorptions = [F(6)/2**4, F(1, 16)*F(6)/(F(3, 2)**4)]
    assert moments[0] == moments[1]
    assert absorptions[0]-absorptions[1] == F(65, 216)
    return {'nu':str(nu), 'T_bounds':[str(lo),str(hi)],
            'Z_bounds':[str(zlo),str(zhi)],
            'stability_margin_bound':str(stable_bound),
            'entropy_width_bound':str(width), 'work_variation_bound':str(work),
            'delta_SW':str(F(201,200)), 'stable_branch_examples':signs,
            'radiation_equal_energy':list(map(str,moments)),
            'radiation_Q_difference':str(absorptions[0]-absorptions[1])}


def path(t, start, end):
    s = t/1000.
    f = 10*s**3-15*s**4+6*s**5
    fp = 30*s*s*(1-s)*(1-s)
    fpp = 60*s-180*s*s+120*s**3
    mass = start+(end-start)*f
    q = (end-start)*fp/1000.
    dq = (end-start)*fpp/1000.**2
    th = 1./mass
    temp = th+q
    dt = -q/mass**2+dq
    return mass, q, th, temp, dt


def trajectory(start, end, steps):
    z0 = 1000./start**3
    state = np.array([z0, 1., 0.])
    initial_energy = start+z0/start+1.
    h = 1000./steps
    max_balance = 0.
    ranges = np.array([[np.inf,-np.inf],[np.inf,-np.inf],[np.inf,-np.inf]])

    def deriv(t, y):
        mass,q,th,temp,dt = path(t,start,end)
        return np.array([-.75*q/temp, -q/4-y[0]*dt, q*q/(temp*th)])

    for i in range(steps+1):
        t = i*h
        mass,q,th,temp,dt = path(t,start,end)
        z,ec,sigma = state
        x = z/temp**3
        max_balance = max(max_balance, abs(mass+z*temp+ec-initial_energy))
        vals = np.array([temp,z,ec])
        ranges[:,0] = np.minimum(ranges[:,0],vals)
        ranges[:,1] = np.maximum(ranges[:,1],vals)
        assert .098 < temp < .102 and .2 < z < 1.8 and x > 0 and ec > 0
        assert 4*x < mass**5  # neighbouring equilibrium's stability
        if i == steps:
            break
        k1 = deriv(t,state)
        k2 = deriv(t+h/2,state+h*k1/2)
        k3 = deriv(t+h/2,state+h*k2/2)
        k4 = deriv(t+h,state+h*k3)
        state += h*(k1+2*k2+2*k3+k4)/6
    j = -4*(state[0]-z0)/3
    dsw = (end*end-start*start)/2
    balance = abs(dsw-j-state[2])
    assert balance < 2e-10 and max_balance < 2e-10
    for t in [0.,1000.]:
        m,q,th,temp,dt = path(t,start,end)
        assert abs(q) < 1e-15 and abs(dt) < 1e-15 and abs(temp-th) < 1e-15
    return {'steps':steps,'start':start,'end':end,'signed_meter_integral':j,
            'entropy_production':float(state[2]),'entropy_identity_error':balance,
            'max_energy_balance_error':max_balance,'final_Z':float(state[0]),
            'final_work_balance':float(state[1]),'final_X':float(state[0]/temp**3),
            'ranges_T_Z_EC':ranges.tolist()}


def calculate():
    exact = exact_checks()
    scans=[]
    for n in [512,1024]:
        plus=trajectory(10.,10.1,n)
        minus=trajectory(10.1,10.,n)
        jplus=plus['signed_meter_integral']
        jminus=-minus['signed_meter_integral']
        assert jplus < 1.005 < jminus
        assert jminus-jplus < float(F(75,19208))
        scans.append({'forward':plus,'reverse':minus,'J_plus':jplus,
                      'J_minus':jminus,'width':jminus-jplus})
    diff=max(abs(scans[0][k]-scans[1][k]) for k in ['J_plus','J_minus','width'])
    assert diff < 2e-10
    return {'status':'PASS_FINITE_SCALAR_ENTROPY_METER_ONLY','exact':exact,
            'numerical_diagnostics':scans,'refinement_max_difference':diff,
            'scope':{'actual_device_certified':False,'covariant_source_certified':False,
                     'all_state_entropy_certified':False,'GR_generated':False,
                     'new_numbered_round':False}}


def close(a,b):
    if isinstance(a,dict): return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    if isinstance(a,list): return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    if isinstance(a,float): return abs(a-b) <= 1e-11*max(1.,abs(a),abs(b))
    return a==b


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    target=HERE/'results.json'
    if args.write:
        assert not target.exists(),'Refuse overwrite'
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        mode='CREATED'
    else:
        assert close(result,json.loads(target.read_text('utf-8-sig'))),'Saved result mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'mode':mode,
          'J_plus':result['numerical_diagnostics'][-1]['J_plus'],
          'J_minus':result['numerical_diagnostics'][-1]['J_minus'],
          'width':result['numerical_diagnostics'][-1]['width']},ensure_ascii=False))
