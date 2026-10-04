"""703 entry: original record effect, bounded source derivatives and time cut.

Reuses702's declared radial/neutral diagnostic; not a new completed round or
full graph thermal evaluation. Analytic full-model bounds are in the note.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_gibbs_preserving_transfer as previous
TARGET=HERE/'bounded_source_entry_results.json'


def multiply(x,y):
    return (x[0]@y[0],x[1]@y[0]+x[0]@y[1],
            x[2]@y[0]+2*x[1]@y[1]+x[0]@y[2])


def power_derivatives(s,n):
    out=(np.eye(len(s[0])),np.zeros_like(s[0]),np.zeros_like(s[0]))
    while n:
        if n&1:out=multiply(out,s)
        n//=2
        if n:s=multiply(s,s)
    return out


def run():
    assert previous.run()==json.loads(previous.TARGET.read_text('utf8'))
    K,W,B,effects=previous.interacting_radial_diagnostic()
    O=effects[0]@effects[0];o=np.diag(O).real
    beta=1.2;theta=.5;lambda0=.23
    A=K+theta*W;D=(1-theta)*W+B;H=A+D+lambda0*O
    assert np.linalg.norm(D@O-O@D)<1e-12
    e,v,rho,Z,_,_=previous.thermal(H,beta)
    weight=np.exp(-beta*e);ot=v.conj().T@O@v
    z=beta*(e[None,:]-e[:,None])/2
    sinhc=np.divide(np.sinh(z),z,out=np.ones_like(z),where=abs(z)>1e-10)
    kernel=np.exp(-beta*(e[:,None]+e[None,:])/2)*sinhc
    deriv=v@(-beta*kernel*ot)@v.conj().T
    first=float(np.trace(deriv).real)
    second=beta**2*float(np.sum(kernel*abs(ot)**2))
    rhoprime=deriv/Z-rho*first/Z
    logfirst=first/Z;logsecond=second/Z-logfirst**2
    target_prob=float(np.trace(O@rho).real)
    target_response=float(np.trace(O@rhoprime).real)
    assert abs(target_prob+logfirst/beta)<1e-12
    assert abs(target_response+logsecond/beta)<1e-12
    c=max(0.,-float(np.linalg.eigvalsh(D)[0]));ZA=float(np.exp(-beta*np.linalg.eigvalsh(A)).sum())
    rows=[]
    for n in (4,8,16,32,64):
        a=beta/n;half=previous.exp_h(A,a/2);middle=previous.exp_h(D,a)*np.exp(-a*lambda0*o)[None,:]
        pieces=tuple(half@(middle*((-a*o)**j)[None,:])@half for j in range(3))
        p,p1,p2=power_derivatives(pieces,n)
        zn=float(np.trace(p).real);z1=float(np.trace(p1).real);z2=float(np.trace(p2).real)
        rn=p/zn;rnprime=p1/zn-p*z1/zn**2
        lnfirst=z1/zn;lnsecond=z2/zn-lnfirst**2
        prob=float(np.trace(O@rn).real);response=float(np.trace(O@rnprime).real)
        rows.append(dict(time_slices=n,step=a,
            partition_error=abs(zn-Z),partition_first_derivative_error=abs(z1-first),
            partition_second_derivative_error=abs(z2-second),
            Gibbs_derivative_trace_norm_error=previous.norm1(rnprime-rhoprime),
            original_record_probability=prob,original_record_response=response,
            original_record_response_error=abs(response-target_response),
            partition_score_probability=-lnfirst/beta,
            partition_score_response=-lnsecond/beta,
            finite_cut_probability_gap=abs(prob+lnfirst/beta),
            finite_cut_response_gap=abs(response+lnsecond/beta)))
    for key in ('partition_first_derivative_error','partition_second_derivative_error',
                'Gibbs_derivative_trace_norm_error','original_record_response_error'):
        assert rows[-1][key]<rows[0][key]/100,key
    assert rows[0]['finite_cut_response_gap']>1e-9
    complex_checks=[]
    for lam in (.2+.4j,-.3+.1j):
        n=8;a=beta/n;half=previous.exp_h(A,a/2)
        s=half@(previous.exp_h(D,a)*np.exp(-a*lam*o)[None,:])@half
        p=np.linalg.matrix_power(s,n)
        actual=float(np.linalg.svd(p,compute_uv=False).sum())
        bound=float(np.exp(beta*(c+abs(lam)*.75))*ZA)
        assert actual<=bound*(1+1e-12)
        complex_checks.append(dict(lambda_real=lam.real,lambda_imag=lam.imag,
            actual_trace_norm=actual,Holder_bound=bound))
    deps=('research_note_603.md','research_note_624.md','research_note_643.md',
          'research_note_655.md','research_note_662.md','research_note_665.md',
          'research_note_702.md','joint_gibbs_preserving_transfer.py',
          'joint_gibbs_preserving_transfer_results.json','research_round_702_checks.json')
    return dict(date='2026-10-02',entry_round=703,latest_formal_round=702,new_formal_round=False,
        source='original L_plus(s)^2, bounded scalar multiplication',beta=beta,theta=theta,lambda0=lambda0,
        target=dict(partition=Z,partition_first=first,partition_second=second,
            record_probability=target_prob,record_response=target_response),
        rows=rows,complex_trace_checks=complex_checks,
        scope=dict(original_all_coupling_analytic_bound=True,
            diagnostic_reuses_702_radial_neutral_sector=True,
            source_deformation_is_declared_probe_not_autonomous_instrument=True,
            exact_partition_derivatives_not_finite_difference=True,
            bounded_source_response_only=True,
            general_geometry_or_real_time_source_derivatives_not_proved=True,
            finite_step_score_not_assumed_equal_bare_record_effect=True,
            no_new_dimension_or_GR_claim=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry=703,formal=702,last=result['rows'][-1],all_checks_passed=True)))
