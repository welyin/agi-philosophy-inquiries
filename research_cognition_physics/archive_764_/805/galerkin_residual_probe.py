"""805 working: continuous PDE residual bound, with analytic omitted modes.

Diagnostic Dirac equation on a circle, not the original 3+1 joint background.
Galerkin trajectories are exact in each trial space; their residuals are tested
in the full Fourier space. Midpoint integration errors use Lipschitz bounds.
Floating point checks calibrate estimates; they are not interval arithmetic.
"""
from pathlib import Path
import argparse, json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'galerkin_residual_probe_results.json'
S1=np.array([[0,1],[1,0]],complex)
S2=np.array([[0,-1j],[1j,0]],complex)
S3=np.diag([1.,-1.])
G=.23;M=.4;T=.6

def matrices(n):
    modes=np.arange(-n,n+1);size=2*len(modes)
    h=np.zeros((size,size),complex);d=h.copy()
    for a,k in enumerate(modes):
        for b,l in enumerate(modes):
            hi=np.zeros((2,2),complex);di=hi.copy()
            if k==l:hi=k*S3+M*S1
            if abs(k-l)==1:
                hi=G*S1;di=.05*(k+l)/2*S2
            if k-l==2:di+=.07/(2j)*S1
            if k-l==-2:di-=.07/(2j)*S1
            h[2*a:2*a+2,2*b:2*b+2]=hi
            d[2*a:2*a+2,2*b:2*b+2]=di
    assert np.max(abs(h-h.conj().T))<1e-14
    assert np.max(abs(d-d.conj().T))<1e-14
    u=np.zeros(size,complex);z=u.copy()
    u[2*(n+0)]=1;z[2*(n+1)+1]=1
    # Missing output Fourier modes -(n+1) and n+1. All other residuals vanish.
    outside=np.zeros((4,size),complex)
    outside[:2,:2]=G*S1;outside[2:,-2:]=G*S1
    return h,d,u,z,outside

def trajectory(h,v,times):
    e,q=np.linalg.eigh(h)
    return ((q@((q.conj().T@v)[:,None]*np.exp(-1j*e[:,None]*times))).T)

def calculate(n,steps=2048):
    h,d,u0,z0,b=matrices(n)
    times=(np.arange(steps)+.5)*T/steps
    u=trajectory(h,u0,times);z=trajectory(h,z0,times)
    integrand=np.einsum('ti,ij,tj->t',u.conj(),d,z)
    value=complex(T*np.mean(integrand))
    ru=np.linalg.norm(u@b.T,axis=1);rz=np.linalg.norm(z@b.T,axis=1)
    # Lipschitz norms of ||B exp(-it H)v||, valid through zeros as well.
    lip=float(np.linalg.norm(b@h,2))
    qres=lip*T*T/(4*steps)
    eu=float(T*np.mean(ru)+qres)
    rz1=float(np.sqrt(1+(n+1)**2)*(T*np.mean(rz)+qres))
    ez=float(np.exp(2*G*T)*rz1)
    z_h1=float(np.exp(2*G*T)*np.sqrt(2.))
    # D: H1 -> L2 norm <= .1 + (.05 + .07) = .22.
    wave_error=T*.22*(eu*(z_h1+ez)+ez)
    quad_error=float(np.linalg.norm(h@d-d@h,2)*T*T/(4*steps))
    return dict(modes_each_side=n,time_midpoints=steps,
        integral=[value.real,value.imag],u_L2_error_bound=eu,z_H1_error_bound=ez,
        full_PDE_wave_pairing_error_bound=wave_error,
        quadrature_error_bound=quad_error,total_error_bound=wave_error+quad_error,
        omitted_residual_integral_u=float(T*np.mean(ru)),
        omitted_residual_integral_z=float(T*np.mean(rz)))

def run():
    rows=[calculate(n) for n in (2,4,6)]
    reference=calculate(18,4096)
    ref=complex(*reference['integral'])
    for row in rows:
        discrepancy=abs(complex(*row['integral'])-ref)
        assert discrepancy<=row['total_error_bound']+reference['total_error_bound']+1e-13
        row['larger_discretization_difference']=float(discrepancy)
    assert rows[0]['omitted_residual_integral_z']>1e-5
    assert rows[-1]['omitted_residual_integral_z']<rows[0]['omitted_residual_integral_z']
    return dict(working_round=805,all_checks_passed=True,
        model='i dt psi = (-i sigma3 dx + .4 sigma1 + .46 cos(x) sigma1) psi on S1',
        diagnostic_not_original_background=True,full_omitted_Fourier_residual_included=True,
        trials=rows,larger_discretization_check=reference,
        discrete_residual_would_be_zero=True,
        floating_point_is_not_interval_certificate=True,
        original_PDE_response_evaluated=False,formal_round_completed=False,new_numbered_test_groups=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
