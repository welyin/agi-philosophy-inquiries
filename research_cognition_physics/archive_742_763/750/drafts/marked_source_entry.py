"""750 entry: actual conditional source family of the original occupation record.

Full original128 Nambu local-symbol calibration. Continuum smoothness, Ward
identities and linear gravitational response are analytic, not simulated here.
"""
import sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'round742_drafts'))
import joint_dynamic_continuum_reference as native
import joint_relative_source_development as sources
import gaussian_record_entry as old_pair
TARGET=HERE/'marked_source_entry_results.json'
def mx(a):return float(np.max(abs(a)))
def conditional(P,e,dP=None):
    n=len(P)//2;C=np.block([[np.zeros((n,n)),np.eye(n)],[np.eye(n),np.zeros((n,n))]])
    ce=C@e.conj();E=np.outer(e,e.conj());Ec=np.outer(ce,ce.conj());F=np.eye(2*n)-E-Ec
    p=float(np.vdot(e,P@e).real);v=F@P@e;w=F@P@ce
    D=np.outer(v,v.conj())-np.outer(w,w.conj());base=F@P@F
    P1=base-D/p+E;P0=base+D/(1-p)+Ec
    if dP is None:return p,(P0,P1)
    dp=float(np.vdot(e,dP@e).real);dv=F@dP@e;dw=F@dP@ce
    dD=np.outer(dv,v.conj())+np.outer(v,dv.conj())-np.outer(dw,w.conj())-np.outer(w,dw.conj())
    dbase=F@dP@F
    dP1=dbase-dD/p+D*dp/p**2
    dP0=dbase+dD/(1-p)+D*dp/(1-p)**2
    return p,(P0,P1),dp,(dP0,dP1)
def full_sources(B):
    bg=sources.background(0.);g,phi,a,a0=bg
    zero=[np.zeros_like(x) for x in bg];labels=['energy'];matrices=[B]
    for i in range(3):
        for j in range(i,3):
            args=[x.copy() for x in zero];args[0][i,j]=1;args[0][j,i]=1
            labels.append('metric_'+str(i)+str(j));matrices.append(sources.matrix_and_source(*bg,*args)[1])
    for i in range(5):
        args=[x.copy() for x in zero];args[1][i]=1
        labels.append('scalar_'+str(i));matrices.append(sources.matrix_and_source(*bg,*args)[1])
    for i in range(3):
        for j in range(3):
            args=[x.copy() for x in zero];args[2][i,j]=1
            labels.append('weak_'+str(i)+str(j));matrices.append(sources.matrix_and_source(*bg,*args)[1])
    for i in range(3):
        args=[x.copy() for x in zero];args[3][i]=1
        labels.append('circle_'+str(i));matrices.append(sources.matrix_and_source(*bg,*args)[1])
    return labels,matrices
def fock_check(P,branches,p,matrices,e):
    n=len(P)//2;G=old_pair.covariance(P);f=e[:n].real
    u=np.r_[f,np.zeros(n)];v=np.r_[np.zeros(n),f]
    outside=np.eye(2*n)-np.outer(u,u)-np.outer(v,v)
    a=float(u@G@v);kappa=np.sqrt(1-a*a)
    W=np.array([u,v,u@G@outside/kappa,v@G@outside/kappa]);small=W@G@W.T
    c,d=old_pair.destroy(0),old_pair.destroy(1)
    ws=[c+c.conj().T,-1j*(c-c.conj().T),d+d.conj().T,-1j*(d-d.conj().T)]
    H=-.25j*sum(small[i,j]*ws[i]@ws[j] for i in range(4) for j in range(4))
    _,V=np.linalg.eigh(H);rho=np.outer(V[:,0],V[:,0].conj());N=c.conj().T@c
    lr=[np.eye(4)-N,N];rhos=[l@rho@l for l in lr]
    probs=[float(np.trace(x).real) for x in rhos];rhos=[x/pb for x,pb in zip(rhos,probs)]
    error=max(abs(probs[1]-p),mx(W@W.T-np.eye(4)))
    for r,Pr in enumerate(branches):
        gg=np.array([[float((.5j*np.trace(rhos[r]@(ws[i]@ws[j]-ws[j]@ws[i]))).real)
                      for j in range(4)] for i in range(4)])
        error=max(error,mx(gg-W@old_pair.covariance(Pr)@W.T))
    T=np.block([[np.eye(n),np.eye(n)],[-1j*np.eye(n),1j*np.eye(n)]])
    source_error=0.
    for J in matrices:
        aa=(T@J@T.conj().T).imag/2;aa=W@aa@W.T
        O=.25j*sum(aa[i,j]*ws[i]@ws[j] for i in range(4) for j in range(4))
        for rr,Pr in zip(rhos,branches):
            small_diff=float(np.trace((rr-rho)@O).real)
            full_diff=float(np.trace(J@(Pr-P)).real/2)
            source_error=max(source_error,abs(small_diff-full_diff))
    assert error<1e-10 and source_error<1e-10
    return dict(independent_Fock_conditional_covariance_error=error,
                all_24_source_compressions_Fock_error=source_error)
def run():
    P,dP,B,dB=native.flow(.08,96,derivative=True)
    e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2)
    p,branches,dp,derivatives=conditional(P,e,dP)
    weights=np.array([1-p,p]);dw=np.array([-dp,dp])
    Q=np.outer(e,e.conj());ce=np.r_[e[64:].conj(),e[:64].conj()];Q+=np.outer(ce,ce.conj())
    R=np.eye(128)-2*Q;bar=(P+R@P@R)/2
    mixture=weights[0]*branches[0]+weights[1]*branches[1]
    errors=[mx(mixture-bar)]
    C=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    for Pr in branches:
        errors.extend([mx(Pr@Pr-Pr),mx(C@Pr.conj()@C+Pr-np.eye(128))])
    labels,Js=full_sources(B)
    delta=np.array([[np.trace(J@(Pr-P)).real/2 for J in Js] for Pr in branches])
    average=weights@delta;gap=delta[1]-delta[0]
    covariance=p*(1-p)*np.outer(gap,gap)
    reconstructed=sum(w*np.outer(q-average,q-average) for w,q in zip(weights,delta))
    assert mx(covariance-reconstructed)<1e-11 and max(errors)<2e-11
    pair=fock_check(P,branches,p,Js,e)
    # The same geometry perturbation changes BOTH reference and actual probabilities.
    conditional_Eprime=np.array([np.trace(dB@(Pr-P)+B@(dPr-dP)).real/2
                                 for Pr,dPr in zip(branches,derivatives)])
    weighted=dw@delta[:,0];within=weights@conditional_Eprime
    direct=np.trace(dB@(bar-P)+B@((R@dP@R-dP)/2)).real/2
    assert abs(weighted+within-direct)<1e-10
    h=2e-5;side=[]
    for gamma in (h,-h):
        PP,_,BB,_=native.flow(.08,96,gamma=gamma)
        pp,rr=conditional(PP,e)
        ee=np.array([np.trace(BB@(Pr-PP)).real/2 for Pr in rr])
        side.append((pp,rr,ee))
    der_err=max(mx((side[0][1][i]-side[1][1][i])/(2*h)-derivatives[i]) for i in range(2))
    prob_err=abs((side[0][0]-side[1][0])/(2*h)-dp)
    Eerr=mx((side[0][2]-side[1][2])/(2*h)-conditional_Eprime)
    assert max(der_err,prob_err,Eerr)<3e-7
    return dict(entry_round=750,latest_completed_round=749,formal_test_count_unchanged=3436,
        Nambu_dimension=128,probabilities=weights.tolist(),probability_one_derivative=dp,
        positivity_purity_reality_and_mixture_error=max(errors),source_labels=labels,
        conditional_source_differences=delta.tolist(),source_difference=gap.tolist(),
        nonselective_source_differences=average.tolist(),
        marked_source_covariance=covariance.tolist(),
        marked_source_covariance_rank=int(np.linalg.matrix_rank(covariance,tol=1e-10)),
        Fock_check=pair,
        energy_response=dict(probability_weight_term=float(weighted),conditional_source_term=float(within),
                             sum=float(weighted+within),old_nonselective_response=float(direct),
                             covariance_derivative_error=der_err,probability_derivative_error=prob_err,
                             conditional_energy_derivative_error=Eerr),
        continuum_sources_not_numerically_simulated=True,
        all_24_symbols_are_from_one_original_model=True,
        no_claim_record_label_is_physical_collapse=True,
        no_actual_local_detector_or_nonlinear_geometry_realization=True)
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:r[k] for k in ('entry_round','probabilities','marked_source_covariance_rank','Fock_check','energy_response')},ensure_ascii=False))
