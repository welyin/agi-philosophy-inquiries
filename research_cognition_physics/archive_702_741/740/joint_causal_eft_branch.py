"""740: the original stationary response at a consistent finite loop order.

Keeps the complete memory and the matched potential at the SAME loop order.
Inverse Laplace sampling is a diagnostic of the causal analytic construction,
not an approximation claim about the unstable exactly-resummed equation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
import joint_covariant_response_closure as original
import joint_mixed_neutral_response as mixed
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_causal_eft_branch_results.json'
sys.path.insert(0,str(HERE/'round740_drafts'))
import total_pole_entry as pole_entry


def data():
    old=mixed.old;x,G,V=original.local_data()
    phi=np.array([0.,np.sqrt(old.U0[0]),0.,0.,np.sqrt(old.U0[1])])
    F=old.matter.original.F(phi)
    J=np.linalg.inv(original.massmap.jacobian(phi))[np.ix_([1,4],[1,4])]
    D=np.diag(2*np.sqrt(old.U0))
    V0=J.T@(D@(old.L/2)@D/F**2)@J
    return x,G,V0,V-V0


def kernel_batch(s,channel,n=160):
    x,q=mixed.quad(n);terms=original.channel_terms(channel)
    dim=terms[0][-1].shape[0];out=np.zeros((len(s),dim,dim),complex)
    for _,a,r,kind,C in terms:
        weights=q*original.profile(x,r,kind)
        integral=(s[:,None]**2*x[None,:]/(a*a+s[:,None]**2*x[None,:]**2))@weights
        out+=integral[:,None,None]*C[None,:,:]
    return out


def pulse_run(N=4096,nquad=160):
    period=128.;dt=period/N;gamma=.25;k=.3
    freq=2*np.pi*np.fft.fftfreq(N,dt);s=gamma+1j*freq;w=s*s+k*k
    Q=np.sqrt(w);x,G,V0,V1=data()
    Im=kernel_batch(Q,'neutral',nquad);I4=kernel_batch(Q,'shear',nquad)[:,0,0]
    H0=w[:,None,None]*G+V0
    def leading(b):
        out=np.empty_like(b)
        out[:,:2]=np.einsum('nij,nj->ni',H0,b[:,:2])
        out[:,2]=-6*w*b[:,2];out[:,3]=w*b[:,3]
        return out
    def inverse0(f):
        out=np.empty_like(f)
        out[:,:2]=np.linalg.solve(H0,f[:,:2,None])[:,:,0]
        out[:,2]=-f[:,2]/(6*w);out[:,3]=f[:,3]/w
        return out
    def correction(b):
        z=b[:,:2]+b[:,2,None]*x
        response=w[:,None]*np.einsum('nij,nj->ni',Im,z)
        out=np.empty_like(b)
        out[:,:2]=b[:,:2]@V1-response
        out[:,2]=-response@x
        out[:,3]=w*w*I4*b[:,3]
        return out
    degree=12;rate=2.;v=np.array([.4,-.3,0.,.15])
    B0=(rate/(s+rate))[:,None]**(degree+1)*v
    forcing=leading(B0)
    Q0=correction(B0)
    B1=-inverse0(Q0)
    Q1=correction(B1)
    # Check the full original total symbol, not only the residual formula.
    rows=[]
    for epsilon in (.2,.1,.05):
        finite=B0+epsilon*B1
        residual=leading(finite)+epsilon*correction(finite)-forcing
        error=float(np.max(abs(residual-epsilon**2*Q1)))
        assert error<3e-14
        rows.append(dict(loop_parameter=epsilon,full_symbol_residual_identity_error=error))
    def invert(B,indices):
        times=np.where(indices<=N//2,indices,indices-N)*dt
        values=np.fft.ifft(B,axis=0)[indices]/dt*np.exp(gamma*times[:,None])
        return times,values
    ids=np.arange(int(16/dt)+1)
    t,b0=invert(B0,ids);_,b1=invert(B1,ids)
    _,q0=invert(Q0,ids);_,q1=invert(Q1,ids)
    analytic=np.zeros_like(t)
    mask=t>0
    analytic[mask]=np.exp((degree+1)*np.log(rate)+degree*np.log(t[mask])-rate*t[mask]-math.lgamma(degree+1))
    baseline_error=float(np.max(abs(b0-analytic[:,None]*v)))
    # Independent finite-difference ODE check of the reconstructed correction.
    interior=(t[1:-1]>=.5)&(t[1:-1]<=15.5)
    b=b1.real
    dd=(b[2:]-2*b[1:-1]+b[:-2])/dt**2
    wave=dd+k*k*b[1:-1]
    residual=np.empty_like(wave)
    residual[:,:2]=wave[:,:2]@G+b[1:-1,:2]@V0+q0.real[1:-1,:2]
    residual[:,2]=-6*wave[:,2]+q0.real[1:-1,2]
    residual[:,3]=wave[:,3]+q0.real[1:-1,3]
    fd_error=float(np.max(abs(residual[interior])))
    neg=np.arange(N-int(4/dt),N)
    _,negative=invert(B1,neg)
    source_norm=float(np.max(np.linalg.norm(b0.real,axis=1)))
    correction_norm=float(np.max(np.linalg.norm(b1.real,axis=1)))
    residual_norm=float(np.max(np.linalg.norm(q1.real,axis=1)))
    assert baseline_error<1e-10 and np.max(abs(negative))<1e-10
    return dict(points=N,quadrature=nquad,dt=dt,spatial_k=k,Laplace_line=gamma,
                period=period,observed_window=16.,pulse_degree=degree,pulse_rate=rate,
                pulse_direction=v.tolist(),baseline_inverse_Laplace_error=baseline_error,
                finite_difference_correction_equation_error=fd_error,
                maximum_pre_pulse_alias=float(np.max(abs(negative))),
                inverse_transform_imaginary_error=float(np.max(abs(b1.imag))),
                baseline_window_norm=source_norm,one_loop_correction_window_norm=correction_norm,
                correction_to_baseline_window_norm=correction_norm/source_norm,
                second_order_residual_coefficient_norm=residual_norm,
                residual_rows=rows),b1.real


def causal_pulse_check():
    coarse,bc=pulse_run(4096,384)
    fine,bf=pulse_run(8192,384)
    refined,br=pulse_run(8192,512)
    time_error=float(np.max(abs(bf[::2]-bc)))
    spectral_error=float(np.max(abs(br-bf)))
    ratio=coarse['finite_difference_correction_equation_error']/fine['finite_difference_correction_equation_error']
    assert 3.5<ratio<4.5 and time_error<1e-9 and spectral_error<2e-8
    return dict(coarse=coarse,fine=fine,spectral_refinement=refined,
                time_sampling_refinement_error=time_error,spectral_integral_refinement_error=spectral_error,
                independent_ODE_error_ratio=ratio,
                initial_failed_resolution_check=dict(quadratures=[160,256],
                    observed_difference=7.534834884410178e-8,required_difference=2e-8,
                    resolution_increased_without_relaxing_tolerance=True),
                analytic_retarded_support_not_proved_by_FFT=True,
                source_is_a_mathematical_calibration_not_a_built_record_apparatus=True,
                residual_bound_not_exact_resummed_solution_error=True)


def run():
    poles=pole_entry.run()
    assert poles==json.loads(pole_entry.TARGET.read_text('utf8'))
    pulse=causal_pulse_check();x,G,V0,V1=data()
    assert np.min(np.linalg.eigvalsh(G))>0 and np.min(np.linalg.eigvalsh(V0))>0
    deps=('research_note_601.md','research_note_630.md','research_note_632.md',
          'research_note_735.md','research_note_738.md','research_note_739.md',
          'joint_covariant_response_closure.py','joint_mixed_neutral_response.py',
          'round740_drafts/total_pole_entry.py','round740_drafts/total_pole_entry_results.json')
    return dict(round=740,tests_run=2,failures=0,errors=0,
                checks=['original_total_poles_check','common_causal_pulse_check'],
                results=dict(original_total_poles=poles,common_causal_pulse=pulse,
                    loop_counting_data=dict(x=x.tolist(),G=G.tolist(),V_tree=V0.tolist(),
                        V_one_loop=V1.tolist(),tree_mass_squared_eigenvalues=np.linalg.eigvals(np.linalg.solve(G,V0)).tolist())),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original all-frequency one-loop resummation has gauge-invariant growing poles. A consistently counted finite-order retarded response retains full memory and linear constraints with a derivative-budgeted O(lambda^2) equation residual. No convergence to the unstable exact inverse, microscopic theory, or nonlinear actual-record solution is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
