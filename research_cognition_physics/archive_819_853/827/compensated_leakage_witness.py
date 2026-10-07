"""827: finite CAR calibration of an uncancellable external response component.

The boson is a two-dimensional ordered-covariance diagnostic. This is not a
numerical evaluation of the original curved-background kernel or its bound.
"""
from pathlib import Path
import argparse,itertools,json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'compensated_leakage_witness_results.json'

def run():
    I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.])
    def kron(parts):
        out=np.ones((1,1),complex)
        for part in parts:out=np.kron(out,part)
        return out
    gamma=[kron([Z]*i+[p]+[I]*(3-i)) for i in range(4) for p in (X,Y)]
    eye=np.eye(16,dtype=complex)
    tr=lambda a:np.trace(a)/16
    def comm(a,b):return 1j*(a@b-b@a)
    monomials=[];even=[]
    for degree in range(5):
        for indices in itertools.combinations(range(4),degree):
            m=eye.copy()*(1j)**(degree*(degree-1)//2)
            for i in indices:m=m@gamma[i]
            monomials.append(m)
            if degree in (2,4):even.append(m)
    def cond(a):return sum((tr(m.conj().T@a)*m for m in monomials),np.zeros_like(a))
    outer=lambda a:a-cond(a)
    logical=[1j*gamma[1]@gamma[2],1j*gamma[0]@gamma[2],1j*gamma[0]@gamma[1]]
    assert np.max(abs(logical[0]@logical[1]-1j*logical[2]))<1e-13
    Hcross=sum((.5j*t*gamma[i]@gamma[i+4] for i,t in enumerate((1.,1.5,.75,1.25))),np.zeros_like(eye))
    HIcross=.125j*gamma[0]@gamma[5]-.0625j*gamma[2]@gamma[7]
    Hin=.5*logical[2]+.25*logical[0]
    HIin=.375*logical[1]-.125*logical[2]
    H=Hcross+Hin;HI=HIcross+HIin;K=H+1j*HI
    residuals=[outer(comm(a,K)) for a in even]
    G=np.array([[tr(a.conj().T@b).real for b in residuals] for a in residuals])
    eig=np.linalg.eigvalsh(G);assert eig.min()>.1
    projection_error=max(float(np.max(abs(cond(r)))) for r in residuals)
    assert projection_error<1e-12
    # The leading trace state is precisely the finite support contract of 817.
    state=np.kron(np.diag([1.,0.]),eye/16)
    sig=[np.kron(I,a) for a in logical];ident=np.eye(32)
    def e_log(a):return sum((s@a@s for s in [ident,*sig]),np.zeros_like(a))/4
    def noise(v):
        b=[e_log(s@v) for s in sig]
        assert max(abs(np.trace(state@x)) for x in b)<1e-13
        c=np.array([[np.trace(state@x@y) for y in b] for x in b])
        d=[comm(a,v) for a in sig]
        actual=np.array([[np.trace(state@x@y).real for y in d] for x in d])
        assert np.max(abs(actual-4*(np.trace(c).real*np.eye(3)-c.real)))<1e-12
        return c,b,actual
    V=np.kron(X,H)+np.kron(Y,HI)
    Cold,Bold,oldG=noise(V)
    rows=[]
    for strength in (0.,1.,2.,-3.,10.):
        R=-strength*Hin;RI=-strength*HIin
        correction=R+1j*RI
        Vh=np.kron(X,R)+np.kron(Y,RI)
        C,B,totalG=noise(V+Vh);Ch,Bh,_=noise(Vh)
        cross=np.array([[np.trace(state@(Bold[r]@Bh[s]+Bh[r]@Bold[s])) for s in range(3)] for r in range(3)])
        cross_error=float(np.max(abs(C-Cold-Ch-cross)))
        shifted=[comm(a,K+correction) for a in even]
        gram=np.array([[tr(a.conj().T@b).real for b in shifted] for a in shifted])
        projection_change=max(float(np.max(abs(outer(x)-y))) for x,y in zip(shifted,residuals))
        remainder_min=float(np.linalg.eigvalsh(gram-G).min())
        assert cross_error<1e-12 and projection_change<1e-12 and remainder_min>-1e-12
        assert np.trace(C).real>=3*eig.min()/8-1e-12
        if strength==1.:assert np.max(abs(gram-G))<1e-12
        rows.append(dict(compensation_strength=strength,noise_trace=float(np.trace(C).real),
            old_plus_preparation_without_cross=float(np.trace(Cold+Ch).real),
            cross_trace=float(np.trace(cross).real),full_cross_identity_error=cross_error,
            external_projection_change=projection_change,seven_record_Gram_min=float(np.linalg.eigvalsh(gram).min()),
            Gram_minus_external_Gram_min=remainder_min))
    assert rows[1]['noise_trace']<rows[0]['noise_trace']
    # Allowing forbidden external fermion control would cancel the remaining
    # witness; this tests the scope rather than asserting a universal no-go.
    illegal=-K
    outside_control=float(np.sqrt(tr(outer(illegal).conj().T@outer(illegal)).real))
    erased=max(float(np.linalg.norm(comm(a,K+illegal))) for a in even)
    assert outside_control>.1 and erased==0.
    return dict(round=827,all_checks_passed=True,finite_CAR_modes=4,
        same_finite_trace_state_used_for_all_compensations=True,
        external_response_Gram_eigenvalues=eig.tolist(),external_projection_error=projection_error,
        fixed_noise_trace_lower_bound=float(3*eig.min()/8),rows=rows,
        outside_fermion_control_norm=outside_control,noise_response_after_allowed_scope_violated=erased,
        original_continuum_noise_constant_computed=False,finite_boson_is_covariance_diagnostic=True,
        scope='The outside-CAR component survives every internal-record pure-boson compensation; internal cross terms can cancel. Outside fermion controls invalidate this restriction.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
