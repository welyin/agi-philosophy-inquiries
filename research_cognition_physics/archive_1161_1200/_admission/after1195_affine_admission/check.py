"""Finite algebraic affine admission; no quantum-process or GR generation claim."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
import numpy as np

def data(x,kind):
    if kind=='sqrt':
        f=np.sqrt(1-x);p=-.5/f;pp=-.25/f**3
    elif kind=='linear':
        f=1-x/2;p=-.5;pp=0.
    else: raise ValueError(kind)
    b=f-2*x*p
    N=(f-x*p)*(f-3*x*p)-x*x*f*pp
    Y=x*np.sqrt(b)/f**1.5
    Yp=N/(f**2.5*np.sqrt(b))
    U=(f-x*p)/(f**1.5*np.sqrt(b))
    Uy=-p/b
    Uyy=-(f*pp+p*p)/(b*b*Yp)
    return dict(f=f,p=p,pp=pp,b=b,N=N,Y=Y,Yp=Yp,U=U,Uy=Uy,Uyy=Uyy)
def lag(P,v,kind):
    x=float(v@np.linalg.solve(P,v))
    return np.sqrt(abs(np.linalg.det(P)))*data(x,kind)['f']
def matrix_check(kind,x):
    P=np.diag([-1.,1.,1.,1.]);v=np.array([np.sqrt(-x),0.,0.,0.])
    d=data(x,kind);s=np.sqrt(abs(np.linalg.det(P)))
    inv=np.linalg.inv(P);w=inv@v
    Q=s*(d['f']*inv-2*d['p']*np.outer(w,w))
    vol=np.sqrt(-np.linalg.det(Q));gi=Q/vol;g=np.linalg.inv(gi)
    Y=float(v@gi@v)
    U=(.5*np.sum(Q*P)-lag(P,v,kind))/vol
    recovered=g/np.sqrt(d['f']*d['b'])-2*d['p']/d['b']*np.outer(v,v)
    A=np.array([[1.,.1,0.,0.],[0.,1.2,.1,0.],[0.,0.,.9,.15],[.05,0.,0.,1.1]])
    density_error=abs(lag(A.T@P@A,A.T@v,kind)-abs(np.linalg.det(A))*lag(P,v,kind))
    H=np.array([[.2,.07,-.02,0.],[.07,-.3,.04,.01],[-.02,.04,.1,-.02],[0.,.01,-.02,.05]])
    eps=2.e-5
    diff=(lag(P+eps*H,v,kind)-lag(P-eps*H,v,kind))/(2*eps)
    xeps=1.e-5
    left,right=data(x-xeps,kind),data(x+xeps,kind)
    derivative=(right['U']-left['U'])/(right['Y']-left['Y'])
    second=(right['Uy']-left['Uy'])/(right['Y']-left['Y'])
    errors={'metric_Y':abs(Y-d['Y']),'Legendre_U':abs(U-d['U']),
      'inverse_P':float(np.max(abs(recovered-P))),'density_covariance':density_error,
      'L_P_pairing':abs(diff-.5*np.sum(Q*H)),'U_Y':abs(derivative-d['Uy']),
      'U_YY':abs(second-d['Uyy'])}
    assert max(errors.values())<2.e-8,errors
    if kind=='sqrt':
        assert abs(d['U']-(1+d['Y']/2))<1.e-12
        assert abs(d['Uyy'])<1.e-12
        K=P-np.outer(v,v)
        assert abs(lag(P,v,kind)-np.sqrt(abs(np.linalg.det(K))))<1.e-12
    return {'kind':kind,'X':x,'errors':errors}
def run():
    samples=[matrix_check(k,x) for k in ('sqrt','linear') for x in (-.2,-.5)]
    x=F(-1,2);f=1-x/2;b=1+x/2
    Uy=1/(2*b)
    cs2=(1+x)/(1+x*x/2)
    time_coeff=Uy/cs2
    assert (f,b,Uy,time_coeff,cs2)==(F(5,4),F(3,4),F(2,3),F(3,2),F(4,9))
    d=data(float(x),'linear')
    assert abs(d['Uy']/(d['Uy']+2*d['Y']*d['Uyy'])-float(cs2))<1.e-14
    rho=d['U']-2*d['Y']*d['Uy']
    assert abs(rho-8*np.sqrt(15)/25)<1.e-14 and rho>0
    return {'status':'PASS_WITH_SCOPE','new_formal_rounds':0,'samples':samples,
     'finite_counterexample':{'signature':'-+++','c':1,'M_squared':1,
      'f':'1-X/2','X':str(x),'f_value':str(f),'b':str(b),'N':'1/2',
      'gradient_coefficient_U_Y':str(Uy),'time_coefficient':str(time_coeff),
      'sound_speed_squared':str(cs2),'energy_density_exact':'8*sqrt(15)/25','energy_density':rho,
      'nondegenerate_Y_prime':d['Yp']},
     'interpretation':'All regular members are classically Einstein plus k-essence; common-cone condition selects the scalar kinetic function only.',
     'not_certified':['cognitive origin of Ricci-only action class','actual clocks or source instruments',
        'quantum process and arbitrary reference','full GR goal','UV front velocity','dark-sector identity']}
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true');args=ap.parse_args()
    o=run();p=Path(__file__).with_name('results.json')
    if args.write_results:
        with p.open('x',encoding='utf-8') as h:json.dump(o,h,ensure_ascii=False,indent=2);h.write('\n')
    elif p.exists():assert o==json.loads(p.read_text(encoding='utf-8'))
    print(json.dumps(o,ensure_ascii=False,indent=2))
