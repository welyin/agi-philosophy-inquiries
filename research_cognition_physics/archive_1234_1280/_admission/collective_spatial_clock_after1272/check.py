"""Common constrained spatial symbol and clock; no actual quantum instrument.

The Fraction certificates check exact algebraic coefficients and conservative
bounds. NumPy grids diagnose the same full equations and Taylor remainder.
Default execution is read-only against the stored results.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import numpy as np

HERE = Path(__file__).resolve().parent


def symbol(k):
    omega = np.sum(1-np.cos(k), axis=-1)
    return 7-2*np.exp(-5*omega)


def initial(u):
    c = np.full_like(u, 0.25)
    d = 10/u-1.25
    total = c+d
    t1sq = d*u/(total)-0.25
    t1 = np.sqrt(t1sq)
    t2 = -c*t1/d
    return c, d, t1, t2


def state(u, tau):
    c, d, t10, t20 = initial(u)
    bb = 2*t10
    aa1 = bb*bb+1
    aa2 = (bb*bb+1)/(4*d)
    p1, p2 = c+bb*tau+aa1*tau*tau, d-bb*tau+aa2*tau*tau
    v1, v2 = bb+2*aa1*tau, -bb+2*aa2*tau
    acc1, acc2 = 2*aa1, 2*aa2
    den = 1+0.4*u*u*tau*tau
    denp = 0.8*u*u*tau
    uu = u/den
    ss = -0.2*u*u*tau/den
    up = -u*denp/(den*den)
    sp = (-0.2*u*u*den+0.2*u*u*tau*denp)/(den*den)
    t1, t2 = ss/2+v1/(8*p1), ss/2+v2/(8*p2)
    tv1 = sp/2+(acc1*p1-v1*v1)/(8*p1*p1)
    tv2 = sp/2+(acc2*p2-v2*v2)/(8*p2*p2)
    return p1, p2, v1, v2, t1, t2, tv1, tv2, uu, ss, up, sp


def maxabs(x):
    return float(np.max(np.abs(x)))


def diagnostic(u, tau):
    p1,p2,v1,v2,t1,t2,tv1,tv2,uu,ss,up,sp = state(u,tau)
    qc = 1+p1+p2
    uc = ss*ss+4*p1*t1*t1+4*p2*t2*t2+1/(16*p1)+1/(16*p2)
    res = {
        'U_definition': maxabs(uc-uu),
        'Hamiltonian_constraint': maxabs(-uu+0.1*uu*qc*uu),
        'momentum_constraint': maxabs(ss+2*(p1*t1+p2*t2)),
        'U_ode': maxabs(up-4*ss*uu),
        's_ode': maxabs(sp+0.2*uu*uu-2*ss*ss),
    }
    for name,p,v,t,tv in [('1',p1,v1,t1,tv1),('2',p2,v2,t2,tv2)]:
        res['p'+name+'_ode']=maxabs(v-8*p*t+4*p*ss)
        res['T'+name+'_ode']=maxabs(tv+4*t*t+0.1*uu*uu-1/(16*p*p)-4*ss*t)
    assert max(res.values()) < 3e-12, res
    assert float(np.min(p1)) >= .25-1e-15
    assert float(np.min(p2)) > .125
    assert 3 < float(np.min(v1)) <= float(np.max(v1)) < 5
    return {'max_equation_residual':max(res.values()),
            'p1_min':float(np.min(p1)), 'p2_min':float(np.min(p2)),
            'site_clock':float(np.mean(p1)), 'site_clock_rate':float(np.mean(v1)),
            'clock_spectral_rate_range':[float(np.min(v1)),float(np.max(v1))]}


def encode(x):
    return {'fraction':str(x),'decimal':float(x)}


def run():
    # Exact algebra: no square roots required for these identities.
    rational_checks=[]
    for u in [F(5),F(11,2),F(6),F(13,2),F(7)]:
        c,d=F(1,4),10/u-F(5,4)
        total=c+d
        t1sq=d*u/total-F(1,4)
        bsq=4*t1sq
        a1,a2=(bsq+1)/(4*c),(bsq+1)/(4*d)
        assert F(8,3)<=t1sq<4
        assert 4*c*t1sq+4*c*c*t1sq/d+1/(16*c)+1/(16*d)==u
        assert a1+a2==4*u
        assert 4*c*a1-bsq==4*d*a2-bsq==1
        assert 1+total==10/u
        rational_checks.append({'u':str(u),'p2_initial':str(d),'T1_squared':str(t1sq),
                                'p_quadratic_coefficients':[str(a1),str(a2)]})
    # Rational enclosure constants used in the analytic proof, not samples.
    p2_lower=F(5,28)-F(4,100)
    clock_upper=4+F(34,100)
    metric_lower=60*(1-F(30,10000))/(1+F(10,10000))**3
    v_upper=F(2,5)*49/F(10000)
    assert p2_lower>F(1,8) and clock_upper<5
    assert metric_lower>59 and v_upper<F(1,500)
    assert 1-14*F(1,500)>0
    taylor_bound=F(90,8)+F(14,24)
    relative_bound=taylor_bound/(F(59,12))*F(1,100)
    assert taylor_bound==F(71,6) and relative_bound==F(142,5900)<F(1,40)

    grids=[]
    for ell in (8,16,32,64):
        line=2*np.pi*np.fft.fftfreq(ell)
        k=np.stack(np.meshgrid(line,line,line,indexing='ij'),axis=-1)
        u=symbol(k)
        assert float(u.min())>=5 and float(u.max())<=7
        times=[{'tau':tau,**diagnostic(u,tau)} for tau in (0.,.005,.01)]
        grids.append({'ell':ell,'sites':ell**3,'u_range':[float(u.min()),float(u.max())],
                      'checks':times})
    # Smooth-symbol Taylor diagnostic: includes axes and non-axial directions.
    rng=np.random.default_rng(1272)
    directions=rng.normal(size=(128,3))
    directions/=np.linalg.norm(directions,axis=1)[:,None]
    directions=np.concatenate([np.eye(3),directions])
    k=np.concatenate([radius*directions for radius in (.02,.05,.1)])
    q=np.sum(k*k,axis=1)
    u=symbol(k)
    max_scaled=0.;max_relative=0.
    for tau in np.linspace(0,.01,21):
        vv=.4*u*u*tau*tau
        f=.1*u*u*(1-vv)/(1+vv)**2
        f0=2.5*(1-10*tau*tau)/(1+10*tau*tau)**2
        aa=2.5*(1-30*tau*tau)/(1+10*tau*tau)**3
        remainder=np.abs(f-f0-2*aa*q)
        assert np.all(remainder<=float(taylor_bound)*q*q+2e-14)
        max_scaled=max(max_scaled,float(np.max(remainder/(q*q))))
        max_relative=max(max_relative,float(np.max(remainder/(2*aa*q))))
    return {
        'status':'PASS_COMMON_SPATIAL_COLLECTIVE_CLOCK_SCOPE',
        'exact_rational_identity_diagnostics':rational_checks,
        'uniform_certificates':{
            'time_window':['0','1/100'],'p1_lower':'1/4',
            'p2_lower':encode(p2_lower),'site_clock_derivative_lower_strict':'3',
            'site_clock_derivative_upper':encode(clock_upper),
            'inverse_spatial_metric_coefficient_lower':encode(metric_lower),
            'inverse_spatial_metric_coefficient_upper':'60',
            'v_upper':encode(v_upper),'spatial_taylor_absolute_factor':encode(taylor_bound),
            'momentum_norm_upper':'1/10','spatial_taylor_relative_upper':encode(relative_bound)},
        'finite_periodic_grid_diagnostics':grids,
        'smooth_symbol_taylor_diagnostic':{'wave_vectors':len(k),'times':21,
             'largest_observed_quartic_factor':max_scaled,
             'largest_observed_relative_remainder':max_relative},
        'scope':{
            'three_dimensional_lattice_and_initial_symbol_are_inputs':True,
            'same_positive_collective_initial_data_and_flow':True,
            'geometry_uses_original_analytic_symbol_long_wave_dictionary':True,
            'finite_lattice_metric_alias_error_certified':False,
            'full_TT_spacetime_remainder_certified':False,
            'normal_quantum_state_or_actual_probe_certified':False,
            'cognitive_axioms_select_parent_constraints':False,
            'new_numbered_round':False,'new_axiom':False}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=run()
    dest=HERE/'results.json'
    if args.write:
        with dest.open('x',encoding='utf-8') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
        print('SAVED_COMMON_SPATIAL_COLLECTIVE_CLOCK_SCOPE')
    else:
        old=json.loads(dest.read_text(encoding='utf-8'))
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for key in a: compare(a[key],b[key])
            elif isinstance(a,list):
                assert len(a)==len(b)
                for x,y in zip(a,b): compare(x,y)
            elif isinstance(a,float):
                assert np.isclose(a,b,rtol=2e-10,atol=2e-12),(a,b)
            else: assert a==b,(a,b)
        compare(result,old)
        print('READ_ONLY_MATCH: common full constraints, positive clocks and spatial Taylor bound')
