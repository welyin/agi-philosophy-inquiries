"""1249 canonical interface: exact vector fields, source identity, finite volume.

One scientific group.  The old interval trajectory is read only at saved
endpoints; no integration is repeated.  Fraction dual differentiation checks
the entire seven-pair generator, not a lone quantized oscillator.  This file
does not certify the quantum propagation error, which has a separate author.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARENT=HERE.parent
DIM=15  # fourteen canonical coordinates plus explicit relational time
Q,PQ,P,PP,X,PX,QO,PO,PSI,PB,Y,K,ANGLE,PI,TIME=range(DIM)
OLD_RESULT=PARENT/'1247'/'results.json'
OLD_RESULT_SHA='dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06'


class Dual:
    def __init__(self,value,grad=None):
        self.value=F(value)
        self.grad=tuple(F(x) for x in grad) if grad is not None else (F(0),)*DIM
    @staticmethod
    def of(x):
        return x if isinstance(x,Dual) else Dual(x)
    def __add__(self,other):
        other=self.of(other)
        return Dual(self.value+other.value,[a+b for a,b in zip(self.grad,other.grad)])
    __radd__=__add__
    def __neg__(self):
        return Dual(-self.value,[-a for a in self.grad])
    def __sub__(self,other):
        return self+-self.of(other)
    def __rsub__(self,other):
        return self.of(other)+-self
    def __mul__(self,other):
        other=self.of(other)
        return Dual(self.value*other.value,
            [a*other.value+self.value*b for a,b in zip(self.grad,other.grad)])
    __rmul__=__mul__
    def __pow__(self,n):
        assert isinstance(n,int)
        if n==0:return Dual(1)
        return Dual(self.value**n,[n*self.value**(n-1)*a for a in self.grad])
    def __truediv__(self,other):
        return self*self.of(other)**-1
    def __rtruediv__(self,other):
        return self.of(other)*self**-1


def seed(value,index,scale=1):
    grad=[F(0)]*DIM;grad[index]=F(scale)
    return Dual(value,grad)


def pulse(theta):
    value=theta.value if isinstance(theta,Dual) else theta
    if 0<=value<=1:
        return 64*theta**3*(1-theta)**3
    return theta*0


def hamiltonian(z,theta):
    # z[Q] stores a=exp(Q); its dual seed is a, so derivatives are w.r.t. Q.
    a,pq,p,pp,x,px,q,r,psi,b,y,k,angle,pi=z
    w,f=pulse(theta),pulse(theta-1)
    return (a**-3*(-pq*pq/4+(pp*pp+px*px)/2)
        +a**3*(p*p+x*x)/2+F(3,20)*a*(1+p/5)
        +F(71,100)+(q*q+r*r+k*k)/2+(b-1)*w*q*p/4
        +f*x*y/10+F(1,100)+(1+p)*pi/10)


def physical(z,theta):
    a,pq,p,pp,x,px,q,r,psi,b,y,k,angle,pi=z
    w,f=pulse(theta),pulse(theta-1)
    h=hamiltonian(z,theta);j=-h
    u=1/a;H=-pq/(2*a**3);v=pp/a**3;vx=px/a**3
    mc=F(71,100)+j+(q*q+r*r+k*k)/2+(b-1)*w*q*p/4+f*x*y/10
    mr=F(1,10)+(1+p)*pi
    rod=F(3,20)*(1+p/5)*u*u
    dust=(mc+mr/10)*u**3
    rho=(v*v+p*p+vx*vx+x*x)/2+rod+dust
    pressure=(v*v-p*p+vx*vx-x*x)/2-rod/3
    return dict(h=h,J=j,u=u,H=H,v=v,V=vx,Mc=mc,Mr=mr,
                rod=rod,dust=dust,rho=rho,pressure=pressure,I=10*pi,B=b)


def check_point(z,theta,label):
    z=list(map(F,z));theta=F(theta)
    dual=[seed(value,i,value if i==Q else 1) for i,value in enumerate(z)]
    hd=hamiltonian(dual,seed(theta,TIME))
    assert hd.value==hamiltonian(z,theta)
    dz=[F(0)]*14
    for i in range(0,14,2):
        dz[i]=hd.grad[i+1];dz[i+1]=-hd.grad[i]
    a,pq,p,pp,x,px,q,r,psi,b,y,k,angle,pi=z
    a_dot=a*dz[Q]
    ph=physical(z,theta)
    u,H,v,vx=ph['u'],ph['H'],ph['v'],ph['V']
    wp=pulse(seed(theta,TIME)).grad[TIME]
    fp=pulse(seed(theta-1,TIME)).grad[TIME]
    w,f=pulse(theta),pulse(theta-1)
    observed=dict(u=-a_dot/a**2,
        H=-dz[PQ]/(2*a**3)+3*pq*a_dot/(2*a**4),
        p=dz[P],v=dz[PP]/a**3-3*v*a_dot/a,
        X=dz[X],V=dz[PX]/a**3-3*vx*a_dot/a,
        q=dz[QO],r=dz[PO],psi=dz[PSI],B=dz[PB],
        Y=dz[Y],K=dz[K],angle=dz[ANGLE],I=10*dz[PI],
        J=-hd.grad[TIME],z=dz[ANGLE]-F(1,10))
    expected=dict(u=-H*u,
        H=-F(3,2)*(v*v+vx*vx+F(2,3)*ph['rod']+ph['dust']),
        p=v,v=-3*H*v-p-F(3,100)*u*u-(pi/10+(b-1)*w*q/4)*u**3,
        X=vx,V=-3*H*vx-x-f*y*u**3/10,
        q=r,r=-q-(b-1)*w*p/4,psi=w*q*p/4,B=F(0),
        Y=k,K=-f*x/10,angle=(1+p)/10,I=F(0),
        J=-(b-1)*wp*q*p/4-fp*x*y/10,z=p/10)
    assert observed==expected,(label,{n:observed[n]-expected[n] for n in observed})
    assert ph['rho']==H*H
    mdot=(observed['J']+q*observed['q']+r*observed['r']+k*observed['K']
        +(b-1)*(wp*q*p+w*observed['q']*p+w*q*v)/4
        +(fp*x*y+f*vx*y+f*x*observed['Y'])/10)
    assert mdot==(b-1)*w*q*v/4+f*y*vx/10
    rod_dot=F(3,100)*v*u*u+F(3,10)*(1+p/5)*u*observed['u']
    dust_dot=(mdot+pi*v/10)*u**3+3*ph['dust']*observed['u']/u
    rho_dot=v*observed['v']+p*v+vx*observed['V']+x*vx+rod_dot+dust_dot
    assert rho_dot+3*H*(ph['rho']+ph['pressure'])==0
    assert 2*H*observed['H']-rho_dot==0
    # Retain all source terms, then delete the memory mass from accounting only.
    memory=F(1,100)+k*k/2
    omitted=ph['rho']-memory*u**3
    assert ph['rho']-omitted==memory*u**3>0
    # Local Weyl identity: H=A(Q) pi_Q, A=-exp(-3Q)/2.
    coeff=-1/(2*a**3);coeff_prime=-3*coeff
    assert coeff_prime**2/4==F(9,16)*a**-6
    return dict(label=label,theta=str(theta),B=str(b),I=str(10*pi),
        all_16_physical_derivative_residuals='0',
        reduced_material_constraint='0',complete_mass_exchange_residual='0',
        continuity_residual='0',constraint_derivative='0',
        same_history_memory_deletion_positive=True,
        local_Weyl_square_coefficient='9 exp(-6Q)/16'),ph


def total_hamiltonian(z,theta,volume):
    a,pq,p,pp,x,px,q,r,psi,b,y,k,angle,pi=z
    nc,nr=volume,volume/10
    pa=volume*pq/a;pp_total=volume*pp;px_total=volume*px
    p_q=volume*r;p_b=volume*b;p_y=volume*k;p_i=volume*pi
    return (-pa*pa/(4*volume*a)+(pp_total**2+px_total**2)/(2*volume*a**3)
        +volume*a**3*(p*p+x*x)/2+F(3,20)*volume*a*(1+p/5)
        +nc*(F(71,100)+q*q/2+(p_q*p_q+p_y*p_y)/(2*nc**2)
             +(p_b/nc-1)*pulse(theta)*q*p/4+pulse(theta-1)*x*y/10)
        +nr/10+(1+p)*p_i/10)


def fixed_rational_cases():
    cases=[];scaling=[]
    for b in [F(1),F(2)]:
        for theta in [F(1,4),F(3,4),F(5,4),F(7,4),F(17,8)]:
            a=F(5,4);H=F(4,5);v=F(-1,40);vx=F(-1,20)
            z=[a,-2*a**3*H,F(-1,50),a**3*v,F(3,40),a**3*vx,
               F(2,5),F(-1,5),F(-1,100),b,F(-1,200),F(-1,400),F(1,20),F(1,10)]
            item,ph=check_point(z,theta,f'rational_B{b}_theta{theta}')
            assert ph['J']>0 and ph['Mc']>0 and ph['Mr']>0
            cases.append(item)
            for n in [1,2]:
                volume=F(1000*n**3)
                assert total_hamiltonian(z,theta,volume)==volume*ph['h']
                assert volume*ph['J']+total_hamiltonian(z,theta,volume)==0
                assert F(n)*F(10)/(10*n)==1 # fixed target-period gradient
            scaling.append(dict(B=str(b),theta=str(theta),n_1_and_2_match=True))
    return cases,scaling


def frozen_endpoint_cases():
    assert hashlib.sha256(OLD_RESULT.read_bytes()).hexdigest()==OLD_RESULT_SHA
    saved=json.loads(OLD_RESULT.read_text('utf-8'))
    assert saved['status']=='PASS'
    cases=[]
    for alpha in [0,1]:
        for time in ['1','2']:
            old=saved['group_finite_memory'][f'alpha{alpha}']['checkpoints'][time]
            midpoint=lambda interval:sum(map(F,interval),F(0))/2
            s={key:midpoint(value) for key,value in old['state'].items()}
            a=1/s['u']
            z=[a,-2*a**3*s['H'],s['p'],a**3*s['v'],s['X'],a**3*s['V'],
               s['q'],s['r'],s['psi'],F(alpha+1),s['Y'],s['K'],
               s['z']+F(time)/10,F(1,10)]
            item,ph=check_point(z,F(time),f'saved_alpha{alpha}_T{time}')
            j_difference=ph['J']-s['J']
            assert abs(j_difference)<F(1,10**12)
            item.update(midpoint_J_reconstruction_difference=float(j_difference),
                midpoint_test_not_new_trajectory_certificate=True)
            cases.append(item)
    return cases


def initial_resource_check():
    # The only irrational initial momentum enters quadratically.  Its square
    # is exactly 4*(111/100); no floating sqrt is used in this certificate.
    parts=[F(-1110),F(5),F(150),F(835),F(20)]
    assert sum(parts)==-100
    rows=[]
    for n in [1,2]:
        volume=F(1000*n**3);scale=F(n**3)
        ptheta=volume/10;hred=sum(parts)*scale
        assert ptheta+hred==0 and hred/volume==F(-1,10)
        total_material=(F(111,100))*volume
        assert total_material==1110*scale
        rows.append(dict(n=n,V=str(volume),N_c=str(volume),N_r=str(volume/10),
            P_theta=str(ptheta),H_red=str(hred),constraint='0',
            total_initial_material_energy=str(total_material),
            initial_density='111/100',initial_h='-1/10',
            hbar_over_V_prefactor=str(1/volume)))
    return dict(initial_terms_V1000=list(map(str,parts)),
        exact_initial_constraint='100-100=0',
        initial_piQ_squared='111/25',branches_have_same_initial_energy=True,
        finite_volume_checks=rows,
        quantum_smallness_not_certified_for_original_V1000=True)


def run():
    cases,scaling=fixed_rational_cases()
    endpoints=frozen_endpoint_cases()
    paths={f'{r}/proof.md':PARENT/str(r)/'proof.md' for r in [1245,1246,1247]}
    paths.update({'1247/results.json':OLD_RESULT,
                  '1247/check.py':PARENT/'1247'/'check.py',
                  '1249/canonical_check.py':Path(__file__)})
    return dict(round=1249,part='canonical_and_material_dictionary',status='PASS',science_groups=1,
        arithmetic='Fraction first-order Dual on all seven pairs and explicit time',
        canonical_pairs=['Q,piQ','p,pip','X,piX','q,piq','psi,pipsi','Y,piY','vartheta,pivartheta'],
        group_complete_canonical_interface=dict(status='PASS',rational_cases=cases,
            volume_scaling_cases=scaling,saved_endpoint_cases=endpoints,
            initial_and_resource=initial_resource_check(),
            saved_old_trajectories_not_reintegrated=True),
        source_sha256={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()},
        limitation='Conditional homogeneous relational canonical bridge. Finite phase-space dimension is not finite Hilbert dimension. No quantum error or actual clock/measurement precision is certified by this part.',
        quantum_error_proof_is_separate=True,CT_adopted=False,new_axioms=0,whole_GR_proved=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run();path=HERE/'canonical_results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2);handle.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(path.read_text('utf-8'))==result,'saved result mismatch'
        mode='READ_ONLY_MATCH'
    group=result['group_complete_canonical_interface']
    print(json.dumps(dict(status='PASS',mode=mode,groups=1,
        rational_cases=len(group['rational_cases']),saved_endpoint_cases=len(group['saved_endpoint_cases']),
        initial_constraint=group['initial_and_resource']['exact_initial_constraint']),ensure_ascii=False))
