"""Exploratory joint flow; normalization still under independent audit."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import json
import numpy as np
from protected_pair_rg import inputs,gauge_squared,LOOP

old,b=inputs(); T=b['T']; r=b['r']; xs=old['matched_inverse_couplings']

def initial(q):
    q=np.asarray(q); S=(T-3*q)/r
    return np.array([q,S,np.full_like(q,T/(2*r)),(3*q*q+r*S*S)/T,
                     S,np.full_like(q,T/r),np.ones_like(q),np.ones_like(q)])

def rhs(u,a):
    q,S,z,lh,p,ls,x,y=a
    gy,gw,gc=gauge_squared(-u,xs)
    B=3*gy+9*gw; C=3/8*(gy*gy+2*gy*gw+3*gw*gw)
    return -np.array([
        2*q*(4.5*q+S-17*gy/12-9*gw/4-8*gc),
        S*(6*q+5*S+z-B/2),2*z*(5*z+S),
        24*lh*lh-B*lh+C+4*(3*q+S)*lh-2*(3*q*q+S*S)+2*p*p,
        p*(12*lh+6*ls+8*p+6*q+2*S+4*z-B/2)-4*z*S,
        18*ls*ls+8*ls*z+8*p*p-8*z*z,
        (12*lh+6*q+2*S-B/2)*x+2*p*y,(6*ls+4*z)*y+8*p*x])/LOOP

def flow(q,u,steps=800):
    a=initial(q);du=u/steps
    for j in range(steps):
        t=j*du;k1=rhs(t,a);k2=rhs(t+du/2,a+du*k1/2)
        k3=rhs(t+du/2,a+du*k2/2);k4=rhs(t+du,a+du*k3)
        a+=du*(k1+2*k2+2*k3+k4)/6
    return a

if __name__=='__main__':
    qs=np.linspace(0,T/3,101)
    for u in (0,1,5,10,15,18):
        a=flow(qs,u); q,S,z,lh,p,ls,x,y=a;D=lh*ls-p*p
        hn=ls*x-p*y;sn=lh*y-p*x
        mask=(hn>0)&(sn>0)&(D>0)
        print(json.dumps(dict(u=u,Hnumerator_range=[float(hn.min()),float(hn.max())],
          Snumerator_range=[float(sn.min()),float(sn.max())],Dmin=float(D.min()),
          valid_q0=qs[mask].tolist(),max_ratio=float(np.max(sn[mask]/hn[mask])) if mask.any() else None,
          at_zero=a[:,0].tolist(),at_upper=a[:,-1].tolist())))
