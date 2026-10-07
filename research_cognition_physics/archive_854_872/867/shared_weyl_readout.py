"""867: original finite-loop menu -> shared physical Weyl readout.
Exact Fourier transfer is inherited from866. Numerical noncommuting matrices
only calibrate the normalization/Taylor identity; they are NOT the original
continuum CCR or a truncation of its physical Hilbert space.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'866'))
import joint_loop_instrument as prior
probe=prior.probe
TARGET=HERE/'shared_weyl_readout_results.json'

def finite_menu():
    n=12;a,b,c=np.meshgrid(np.arange(n),np.arange(n),np.arange(n),indexing='ij')
    w=(-a-b)%n;mask=(a!=b)&(a!=w)&(b!=w)
    ry,ay,weights=probe.menu_grid(n)
    return ry[mask],ay[mask],weights[mask]

def herm_fun(M,fun):
    d,v=np.linalg.eigh(M)
    return (v*fun(d)[...,None,:])@v.conj().swapaxes(-1,-2)

def readout(p0,Q,t,delta=.5):
    I=np.eye(2,dtype=complex)
    scaled=t*Q/(delta*p0[:,None,None])
    sine=herm_fun(scaled,np.sin)
    A=p0[:,None,None]*(I+delta*sine)
    B=A.sum(axis=0);F=herm_fun(B,lambda x:1/np.sqrt(x))
    L=herm_fun(A,np.sqrt)@F
    E=L.conj().swapaxes(-1,-2)@L
    return dict(L=L,E=E,normalization=float(np.max(abs(E.sum(axis=0)-I))),
        B_without_normalization_error=float(np.max(abs(B-I))),
        relative_effect_lower=float(np.min(np.linalg.eigvalsh(E)/p0[:,None])))

def channel(L,X):
    return np.sum(L.conj().swapaxes(-1,-2)@X@L,axis=0)

def run():
    assert prior.exact_certificate()==json.loads(prior.TARGET.read_text('utf-8'))['exact_certificate']
    ry,ay,w=finite_menu()
    with probe.ResearchRuntime(probe.Layout()).installed():
        import joint_reference_constraint_strata as old
        U,dU,_,_=probe.original.matrices(old,probe.original.inherited.constants(),.4)
    # Central-anchor calibration. The theorem uses the full smooth anchor
    # average and full source kernels; no delta anchor is asserted physical.
    phase=np.exp(-.74j);r0=phase*np.trace(U);dr=phase*np.trace(dU)
    a0=abs(r0)**2-1;da=2*np.real(np.conj(r0)*dr)
    base=np.array([r0.real,r0.imag,a0])
    tangent=np.array([dr.real,dr.imag,da])
    C=w[:,None]*np.stack([2*probe.ETA_R*ry.real,2*probe.ETA_R*ry.imag,probe.ETA_A*ay],axis=-1)
    p0=w+C@base
    values=np.stack([ry.real/probe.ETA_R,ry.imag/probe.ETA_R,ay/probe.ETA_A,1+ay/probe.ETA_A],axis=0)
    target=np.array([[1,0,0],[0,1,0],[0,0,1],[0,0,1]],float)
    transfer_error=float(np.max(abs(values@C-target)))
    assert transfer_error<2e-13 and np.max(abs(C.sum(axis=0)))<1e-15
    assert np.min(p0/w)>.2 and abs(sum(p0)-1)<1e-13
    G=C.T@(C/p0[:,None])
    eig=np.linalg.eigvalsh(G);assert min(eig)>0
    fisher=float(tangent@G@tangent)
    assert fisher>0
    paulis=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    Q=np.einsum('ja,amn->jmn',C,paulis)
    target_observables=np.einsum('ka,amn->kmn',target,paulis)
    X=.3*paulis[0]-.2*paulis[1]+.4*paulis[2]
    comm=Q@X-X@Q
    second=-np.sum((Q@comm-comm@Q)/p0[:,None,None],axis=0)/8
    rows=[]
    for t in (.4,.2,.1,.05):
        plus=readout(p0,Q,t);minus=readout(p0,Q,-t)
        derivative=np.einsum('kj,jmn->kmn',values,(plus['E']-minus['E'])/(2*t))
        first_error=float(np.max(abs(derivative-target_observables)))
        got_second=(channel(plus['L'],X)+channel(minus['L'],X)-2*X)/(2*t*t)
        second_error=float(np.max(abs(got_second-second)))
        assert plus['normalization']<1e-13 and minus['normalization']<1e-13
        assert plus['relative_effect_lower']>(1-.5)/(1+.5)-1e-13
        rows.append(dict(t=t,complete_normalization_error=plus['normalization'],
            omitted_B_normalization_error=plus['B_without_normalization_error'],
            relative_effect_lower=plus['relative_effect_lower'],
            first_all_menu_source_error=first_error,
            universal_second_backaction_error=second_error))
    for key in ('first_all_menu_source_error','universal_second_backaction_error'):
        assert all(3.7<rows[i][key]/rows[i+1][key]<4.3 for i in range(3))
    tail={str(n):float((2+np.sqrt(2))*2.**(-n)) for n in (12,24,40)}
    return dict(round=867,date='2026-10-06',formal_reports=867,
        cumulative_numbered_groups=3652,fresh_numbered_groups=1,all_checks_passed=True,
        original_finite_menu_labels=len(w),all_menu_linear_transfer_error=transfer_error,
        source_coefficient_Gram_eigenvalues=eig.tolist(),
        central_anchor_physical_tangent=tangent.tolist(),
        central_anchor_probability_Fisher=fisher,
        inherited_companion_variance_t2_coefficient=fisher/4,
        physical_nonzero_for_smooth_anchor='analytic continuity and full-rank source transfer',
        matrix_calibration_scope='arbitrary noncommuting self-adjoint source matrices, not continuum CCR',
        matrix_rows=rows,
        weighted_Weyl_series_operator_tail_upper_per_sqrt_p=tail,
        full_source_Ward_and_common_domain='analytic proof using original809/863',
        free_CP_instrument_for_finite_real_readout_parameter=True,
        original_interacting_transport='fixed formal bulk coupling only',
        original_graph_matching='zeroth probability and complete first sources, not nonlinear all-orders',
        leading_same_state_backreaction_and_Ward=True,
        local_probe_implementation_and_strong_causality_proved=False,
        finite_bulk_coupling_convergence=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
