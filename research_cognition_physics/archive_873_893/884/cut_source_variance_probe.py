"""Working884: positive current variance in883 cut strip.
Classical uniform flat-holonomy average, NOT original full dynamical Gauss Gibbs.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'883'))
import spectral_source_bridge as old
TARGET=HERE/'cut_source_variance_probe_results.json'

def step_derivative(z):
    z=np.asarray(z);c=old.smooth_step(z);out=np.zeros_like(c)
    m=(z>0)&(z<1)
    out[m]=c[m]*(1-c[m])*(1/z[m]**2+1/(1-z[m])**2)
    return out

def source(N,p):
    a=2*np.pi/N
    th=np.angle(np.exp(1j*a*np.asarray(p)))
    pp=th/a
    z=(np.pi-abs(th))/(a*old.ETA)
    return old.smooth_step(z)-abs(pp)*step_derivative(z)/old.ETA

def solve_u(N,w):
    lo=0.;hi=old.ETA
    for i in range(65):
        m=(lo+hi)/2
        if (N/2-m)*old.smooth_step(np.array([m/old.ETA]))[0]<w:lo=m
        else:hi=m
    u=(lo+hi)/2;z=np.array([u/old.ETA])
    derivative=-old.smooth_step(z)[0]+(N/2-u)*step_derivative(z)[0]/old.ETA
    assert derivative>0
    return u,float(derivative)

def variance(H,G):
    ev,v=np.linalg.eigh(H);n=.5-.5*np.tanh(old.BETA*ev/2)
    g=v.conj().T@G@v
    return float(np.sum(n[:,None]*(1-n[None,:])*abs(g)**2).real)

def run():
    m=old.load();q,_=old.charges(m)
    x,weights=np.polynomial.legendre.leggauss(32)
    ws=.75*x+1.25;weights=.75*weights
    rows=[]
    for N in (17,33,65,129,257,1025,4097):
        value=0.;continuum=0.
        for w,weight in zip(ws,weights):
            u,wp=solve_u(N,w);alpha=-(.5-u)/6
            variance_sum=0.;base_sum=0.
            for k in (N//2,-(N//2)):
                hn=old.H(m,q,N,[k,0,0],alpha)
                gn=m.GAMMA[0]@np.diag(q*source(N,k+alpha*q))
                hc=old.H(m,q,N,[k,0,0],alpha,kind='continuum')
                gc=m.GAMMA[0]@np.diag(q)
                variance_sum+=variance(hn,gn)/2
                base_sum+=variance(hc,gc)/2
            value+=weight*variance_sum/(6*wp)
            continuum+=weight*base_sum/(6*wp)
        low,_=solve_u(N,.5);high,_=solve_u(N,2.)
        rows.append(dict(N=N,holonomy_strip_width=(high-low)/6,
                         integrated_positive_current_variance=value,
                         variance_over_logN_squared=value/np.log(N)**2,
                         continuum_same_strip_variance=continuum))
    assert all(y['integrated_positive_current_variance']>x['integrated_positive_current_variance']
               for x,y in zip(rows,rows[1:]))
    return dict(working_round=884,formal_round_still=883,cumulative_groups_still=3668,
       previous_goal_turn_classification='progress',original_Nambu_components=64,
       average='Lebesgue d alpha over a shrinking subinterval near -1/12 of the uniform flat-hypercharge parameter circle; no Gauss Gibbs measure asserted.',
       variance='1/2 sum over paired k and -k of Tr[n G (1-n) G], all original mass and charge blocks retained',
       strip_defined_by='Right-electron spectral momentum between0.5 and2; q=-6. k=(N//2,0,0).',
       rows=rows,full_dynamical_gauge_model_refuted=False,
       status='working quantitative evidence; asymptotic proof and scope audit pending')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
