"""565 provisional Galerkin diagnosis of the full-matter large-b limit.

The matrices are finite form compressions. They diagnose, not prove,
convergence of the untruncated model or its physical measurement.
"""
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent;BASE=HERE.parent
if str(BASE) not in sys.path:sys.path.insert(0,str(BASE))
import joint_finite_time_gauge_probe as inherited
TARGET=HERE/'large_b_matter_probe_results.json'


def laguerre_polys(n):
    return [np.array([(-1.)**j*math.comb(m+1,m-j)/math.factorial(j) for j in range(m+1)])/math.sqrt(m+1)
            for m in range(n)]


def radial_moment(n,power):
    polys=laguerre_polys(n)
    return np.array([[sum(c*d*math.gamma(i+j+2+power) for i,c in enumerate(P) for j,d in enumerate(Q))
                      for Q in polys] for P in polys])


def build(nr=2,ns=3,nang=4):
    a=.5;k=.8;L,u,_=inherited.constants()
    Ir=np.eye(nr);Is=np.eye(ns);Ic=np.eye(nr*ns)
    R=radial_moment(nr,.5);R2=radial_moment(nr,1);R4=radial_moment(nr,2);Rinv2=radial_moment(nr,-1)
    assert np.max(abs(radial_moment(nr,0)-Ir))<1e-12
    ext=ns+4
    q=np.diag(np.sqrt(np.arange(1,ext)/2),1)+np.diag(np.sqrt(np.arange(1,ext)/2),-1)
    S2=(q@q)[:ns,:ns];S4=(q@q@q@q)[:ns,:ns]
    Tr=a*np.diag(4*np.arange(nr)+4)-a*R2
    Ts=a*np.diag(2*np.arange(ns)+1)-a*S2
    V=(L[0,0]*np.kron(R4-2*u[0]*R2+u[0]**2*Ir,Is)
       +2*L[0,1]*np.kron(R2-u[0]*Ir,S2-u[1]*Is)
       +L[1,1]*np.kron(Ir,S4-2*u[1]*S2+u[1]**2*Is))/4
    Hcell=np.kron(Tr,Is)+np.kron(Ir,Ts)+V+k*np.kron(R2,Is)/2
    HZ=np.kron(Hcell,Ic)+np.kron(Ic,Hcell)
    Rcell=np.kron(R,Is);R2cell=np.kron(R2,Is);Invcell=np.kron(Rinv2,Is)
    rr=np.kron(Rcell,Rcell)
    sumsq=np.kron(R2cell,Ic)+np.kron(Ic,R2cell)
    readR=sumsq/2-rr
    A0=a*(np.kron(Invcell,Ic)+np.kron(Ic,Invcell))
    Cas=np.diag(np.arange(nang)*(np.arange(nang)+2))
    Z=np.diag(np.full(nang-1,.5),1)+np.diag(np.full(nang-1,.5),-1)
    I0=np.eye(HZ.shape[0]);Ia=np.eye(nang)
    Hrest=np.kron(HZ,Ia)+np.kron(A0,Cas)-k*np.kron(rr,Z)
    penalty=np.kron(I0,Cas)/4
    Hs=lambda b:Hrest+b*penalty
    source=np.zeros(Hrest.shape[0]);source[0]=1.
    sourceZ=np.zeros(HZ.shape[0]);sourceZ[0]=1.
    return dict(H=Hs,HZ=HZ,R=np.kron(readR,Ia),RZ=readR,source=source,sourceZ=sourceZ,
        nang=nang,penalty=penalty,norm_moment_endpoint_correction=float(np.linalg.norm(R4-R2@R2)),
        coordinate_square_compression_difference=float(np.linalg.norm(R2-R@R)),
        E0=float(source@Hrest@source),Hcell=Hcell,dimension=Hrest.shape[0])


def exp_apply(H,psi,t):
    ev,U=np.linalg.eigh(H)
    return U@(np.exp(-1j*t*ev)*(U.T.conj()@psi))


def embed(phi,nang):
    out=np.zeros(len(phi)*nang,dtype=complex);out[::nang]=phi
    return out


def dynamics(m,b):
    H=m['H'](b);E=m['E0'];rows=[]
    for t in (.2,.6):
        psi=exp_apply(H,m['source'],t);target=embed(exp_apply(m['HZ'],m['sourceZ'],t),m['nang'])
        leakage=float(np.sum(abs(psi.reshape(-1,m['nang'])[:,1:])**2))
        assert leakage<=4*E/(3*b)+1e-12
        rows.append(dict(time=t,state_distance=float(np.linalg.norm(psi-target)),
            angular_leakage=leakage,full_H_energy=float(np.vdot(psi,H@psi).real)))
        assert abs(rows[-1]['full_H_energy']-E)<1e-10
    compressed=m['H'](b)[::m['nang'],::m['nang']]
    assert np.max(abs(compressed-m['HZ']))<1e-13
    return dict(b=b,rows=rows,all_time_leakage_bound=4*E/(3*b),
                smallest_compressed_eigenvalue=float(np.linalg.eigvalsh(H)[0]))


def apparatus(m,b,n=20):
    # After an unmeasured source time, keep source H during finite tau.
    # Free pointer spreading cancels only with the declared compensated reader.
    source_time=.25;tau=.2;g=1.;sigmaQ=.5
    H=m['H'](b)
    psi=exp_apply(H,m['source'],source_time)
    phi=exp_apply(m['HZ'],m['sourceZ'],source_time)
    x,w=np.polynomial.hermite.hermgauss(n);p=math.sqrt(2)*x;w/=math.sqrt(math.pi)
    squared=0.
    for momentum,weight in zip(p,w):
        actual=exp_apply(tau*H+g*momentum*m['R'],psi,1.)
        effective=embed(exp_apply(tau*m['HZ']+g*momentum*m['RZ'],phi,1.),m['nang'])
        squared+=weight*np.linalg.norm(actual-effective)**2
    return dict(b=b,source_time=source_time,tau=tau,g=g,sigmaQ=sigmaQ,
        pointer_momentum_quadrature_nodes=n,joint_state_distance=math.sqrt(float(squared)),
        any_compensated_pointer_event_difference_upper=min(1.,math.sqrt(float(squared))))


def run():
    m=build();d=[dynamics(m,b) for b in (8.,32.,128.,512.)]
    apparatus_rows=[apparatus(m,b,20) for b in (16.,64.,256.)]
    refined=apparatus(m,256.,32)
    qdiff=abs(refined['joint_state_distance']-apparatus_rows[-1]['joint_state_distance'])
    assert qdiff<1e-7
    assert d[-1]['rows'][-1]['state_distance']<d[0]['rows'][-1]['state_distance']
    assert apparatus_rows[-1]['joint_state_distance']<apparatus_rows[0]['joint_state_distance']
    # One angular cutoff change, still only a diagnostic of the form compression.
    m5=build(nang=5);d5=dynamics(m5,512.)
    cutoff=max(abs(r['state_distance']-s['state_distance']) for r,s in zip(d[-1]['rows'],d5['rows']))
    assert cutoff<1e-7
    return dict(status='preliminary_565_not_completed',Galerkin_dimension=m['dimension'],
        radial_modes=2,singlet_modes=3,class_modes=4,initial_source_energy=m['E0'],
        radial_fourth_moment_endpoint_correction=m['norm_moment_endpoint_correction'],
        radial_coordinate_square_compression_difference=m['coordinate_square_compression_difference'],
        source_dynamics=d,finite_duration_apparatus=apparatus_rows,
        pointer_integral_refinement_difference=qdiff,angular_cutoff_refinement_difference=cutoff,
        scope=dict(finite_form_compression_only=True,not_untruncated_dynamical_error_certificate=True,
            no_b_dependent_rescaling_of_k_or_W=True,no_spacetime_or_unification_claim=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
