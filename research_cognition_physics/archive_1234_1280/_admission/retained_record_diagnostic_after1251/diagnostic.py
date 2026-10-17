"""Read-only-model diagnostic, NOT a new numbered scientific certificate.

Fixed 1247 alpha=0. Only H' is forced by eta*w(theta), eta=1e-3.
This is off shell: the complete original material stress stays conserved,
but C'=2 H eta w.  Two fixed RK4 resolutions are diagnostics, not bounds.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in [HERE, *HERE.parents]
            if (p/'research_cognition_physics').is_dir())
OLD = ROOT/'research_cognition_physics/archive_1234_/1247'
ETA = 1e-3


def pulses(t, phase):
    x = t if phase == 0 else t-1
    pulse = 64*x**3*(1-x)**3
    derivative = 192*x*x*(1-x)**2*(1-2*x)
    zero = np.zeros_like(t)
    if phase == 0:
        return pulse, derivative, zero, zero
    return zero, zero, pulse, derivative


def source(y, phase):
    u,H,p,v,z,light,theta,J,q,r,psi,X,V,Y,K = y
    w,wp,f,fp = pulses(theta, phase)
    Mc = .71+J+(q*q+r*r+K*K)/2+.1*f*X*Y
    rod = .15*(1+.2*p)*u*u
    dust = (Mc+.02+.01*p)*u**3
    rho = (v*v+p*p+V*V+X*X)/2+rod+dust
    pressure = (v*v-p*p+V*V-X*X)/2-rod/3
    return Mc,rod,dust,rho,pressure


def rhs(y, phase):
    u,H,p,v,z,light,theta,J,q,r,psi,X,V,Y,K = y
    w,wp,f,fp = pulses(theta, phase)
    Mc,rod,dust,rho,pressure = source(y,phase)
    return np.array([-H*u,-1.5*(v*v+V*V+2*rod/3+dust),v,
        -3*H*v-p-.03*u*u-.01*u**3,.1*p,u,np.ones_like(u),
        -.1*fp*X*Y,r,-q,.25*w*q*p,V,-3*H*V-X-.1*f*Y*u**3,K,-.1*f*X])


def jac(y, phase):
    # Polynomial RHS: simultaneous complex-step Jacobian, no subtraction.
    z = y[:,None].astype(complex)+1j*1e-30*np.eye(15)
    return np.imag(rhs(z,phase))/1e-30


def integrate(steps):
    y0 = np.array([1,np.sqrt(1.11),0,0,0,0,0,.1,.5,0,0,.1,0,0,0])
    # background, full fundamental matrix, tangent to eta, forced history
    state = np.r_[y0,np.eye(15).reshape(-1),np.zeros(15),y0]
    h = 1/steps
    history = [state.copy()]
    def fun(a,phase):
        y=a[:15]; U=a[15:240].reshape(15,15); Z=a[240:255]; f=a[255:270]
        A=jac(y,phase); forcing=np.zeros(15)
        forcing[1]=float(pulses(y[6],phase)[0])
        finite=rhs(f,phase); finite[1]+=ETA*float(pulses(f[6],phase)[0])
        return np.r_[rhs(y,phase),(A@U).reshape(-1),A@Z+forcing,finite]
    for phase in [0,1]:
        for j in range(steps):
            k1=fun(state,phase); k2=fun(state+h*k1/2,phase)
            k3=fun(state+h*k2/2,phase); k4=fun(state+h*k3,phase)
            state+=h*(k1+2*k2+2*k3+k4)/6
            state[6]=state[261]=phase+(j+1)*h
            history.append(state.copy())
    states=np.asarray(history)
    times=np.arange(2*steps+1)*h
    final_row=state[15:240].reshape(15,15)[14]
    adj=[]
    for item in states:
        adj.append(np.linalg.solve(item[15:240].reshape(15,15).T,final_row))
    adj=np.asarray(adj)
    weights=np.ones(2*steps+1);weights[1:-1:2]=4;weights[2:-1:2]=2
    w=np.where(times<=1,64*times**3*(1-times)**3,0)
    adj_projection=h*np.dot(weights,adj[:,1]*w)/3
    direct=state[269]-state[14]
    tangent=state[254]
    # Initial X perturbation remains on the original Hamiltonian constraint
    # by adjusting H0 along dH0/dX0 = X0/(2H0).
    initial_X_constrained=final_row[11]+final_row[1]*y0[11]/(2*y0[1])
    constraint=state[256]**2-source(state[255:270],1)[3]
    forced_H=states[:,256]
    integrated_constraint=ETA*h*np.dot(weights,2*forced_H*w)/3
    source_chain=[]
    for j in [steps//4,steps//2,steps+steps//2]:
        phase=0 if j<steps else 1
        yf=states[j,255:270]
        dy=rhs(yf,phase);dy[1]+=ETA*float(pulses(yf[6],phase)[0])
        rdot=np.imag(source(yf.astype(complex)+1j*1e-30*dy,phase)[3])/1e-30
        _,_,_,rho,P=source(yf,phase)
        source_chain.append(float(rdot+3*yf[1]*(rho+P)))
    return dict(steps_per_phase=steps,method='piecewise binary64 RK4 and complex-step full 15D Jacobian; diagnostic only',
        background_K=float(state[14]),forced_K=float(state[269]),delta_K=float(direct),
        tangent_K_per_unit_residual_amplitude=float(tangent),
        adjoint_H_w_integral=float(adj_projection),
        adjoint_vs_tangent_error=float(adj_projection-tangent),
        finite_minus_linear=float(direct-ETA*tangent),
        ratio_nonlinear_to_signal=float(abs(direct-ETA*tangent)/abs(direct)),
        final_constraint_C=float(constraint),integrated_2H_residual=float(integrated_constraint),
        maximum_state_difference=float(np.max(np.abs(states[:,255:270]-states[:,:15]))),
        terminal_read_budget_two_histories=2e-10,
        diagnostic_gap_after_terminal_read=float(abs(direct)-2e-10),
        initial_X_constrained_sensitivity=float(initial_X_constrained),
        initial_X_linear_offset_matching_signal=float(direct/initial_X_constrained),
        initial_K_sensitivity=float(final_row[14]),
        original_full_source_continuity_diagnostics=source_chain,
        adjoint_max_by_component=[float(x) for x in np.max(np.abs(adj),axis=0)],
        adjoint_at_initial=[float(x) for x in adj[0]],
        kernel_samples=[dict(t=float(times[j]),zeta_H=float(adj[j,1]),
            zeta_X=float(adj[j,11]),zeta_K=float(adj[j,14])) for j in range(0,2*steps+1,steps//8)],
        exact_new_permission='none; model off-shell forcing only',
        rigorous_error_bound=False)


def run():
    spec=importlib.util.spec_from_file_location('frozen1247',OLD/'check.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    maximum=0.
    for phase in [0,1]:
        y=np.array([.7,.5,-.01,-.02,0,.2,phase+.25,.1,.3,-.2,0,.08,-.03,-.001,-.002])
        maximum=max(maximum,float(np.max(np.abs(rhs(y,phase)-np.asarray(old.rhs(y,0,phase),dtype=float)))))
    assert maximum<1e-14
    runs=[integrate(512),integrate(1024)]
    result=dict(status='DIAGNOSTIC_ONLY',numbered_round=None,
        direction='fixed original alpha=0; r_H=1e-3*w(theta) on [0,1], zero afterwards',
        finite_model='complete original matter equations and initial values, modified H equation only',
        not_a_covariant_candidate=True,old_rhs_implementation_difference=maximum,
        runs=runs,
        resolution_K_difference=float(runs[1]['delta_K']-runs[0]['delta_K']),
        resolution_adjoint_difference=float(runs[1]['adjoint_H_w_integral']-runs[0]['adjoint_H_w_integral']),
        source_sha256={n:hashlib.sha256((OLD/n).read_bytes()).hexdigest() for n in ['proof.md','check.py','results.json']},
        stop='No interval, preparation, full physical source-matching, or nonlinear error certificate is supplied; numerical nonzero is not admission.')
    return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=run()
    output=HERE/'results.json'
    if args.write:
        with output.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('DIAGNOSTIC_SAVED_EXCLUSIVELY',output)
    else:
        saved=json.loads(output.read_text('utf-8'))
        assert result==saved,'saved diagnostic mismatch'
        print('READ_ONLY_MATCH: DIAGNOSTIC_ONLY; no rigorous numerical error bound')
    runs=result['runs']
    for r in runs:
        print({k:r[k] for k in ['steps_per_phase','delta_K','adjoint_H_w_integral','finite_minus_linear','initial_X_linear_offset_matching_signal','final_constraint_C']})


if __name__=='__main__':main()
