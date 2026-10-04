"""Round 532: a common-weight constraint across gauge selection and dynamics.

Finite representation audit, not a complete Lorentzian spectral model.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'shared_trace_unimodularity_results.json'


def blockdiag(*blocks):
    n=sum(b.shape[0] for b in blocks)
    out=np.zeros((n,n),dtype=complex)
    offset=0
    for b in blocks:
        k=b.shape[0]
        out[offset:offset+k,offset:offset+k]=b
        offset+=k
    return out


def rep(N, lam, quat, color):
    a=blockdiag(np.diag([lam,np.conjugate(lam)]),quat)
    return blockdiag(a,lam*np.eye(4),np.kron(a,np.eye(N)),np.kron(np.eye(4),color))


def structures(N,wq,wl,t=1.0):
    dim=8*(N+1)
    Z=np.diag([wl]*8+[wq]*(8*N)).astype(complex)
    J=np.zeros((dim,dim),complex)
    for a,b,size in ((0,4,4),(8,8+4*N,4*N)):
        J[a:a+size,b:b+size]=np.eye(size)
        J[b:b+size,a:a+size]=np.eye(size)
    gamma=blockdiag(np.diag([-1,-1,1,1]),np.diag([1,1,-1,-1]),
                    np.kron(np.diag([-1,-1,1,1]),np.eye(N)),
                    np.kron(np.diag([1,1,-1,-1]),np.eye(N)))
    def yukawa(a,b):
        s=np.zeros((4,4),complex)
        s[0,2]=s[2,0]=t*a
        s[1,3]=s[3,1]=t*b
        return s
    sl,sq=yukawa(2,3),yukawa(5,7)
    D=blockdiag(sl,sl.conj(),np.kron(sq,np.eye(N)),np.kron(sq.conj(),np.eye(N)))
    D[0,4]=D[4,0]=11  # gauge-neutral Majorana term
    return Z,J,gamma,D


def raw_abelian(N,h,v):
    a=np.diag([h,-h,0,0])
    return blockdiag(a,h*np.eye(4),np.kron(a,np.eye(N)),v*np.eye(4*N))


def lhs_charges(N,h,q):
    return dict(Q=(2*N,q),uc=(N,-q-h),dc=(N,-q+h),
                L=(2,-h),ec=(1,2*h),nuc=(1,F(0)))


def anomalies(N,h,q):
    fs=lhs_charges(N,h,q)
    return {'weak':N*q-h, 'cubic':sum(m*y**3 for m,y in fs.values()),
            'mixed_gravity':sum(m*y for m,y in fs.values()),
            'color':2*q+fs['uc'][1]+fs['dc'][1]}


def close(actual,expected):
    assert np.max(np.abs(np.asarray(actual)-np.asarray(expected)))<1e-9


def run():
    checks=[]
    for N in (3,5):
        Z,J,gamma,D=structures(N,1,2)
        pauli=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]))
        basis=[rep(N,z,np.zeros((2,2)),np.zeros((N,N))) for z in (1,1j)]
        basis.extend(rep(N,0,q,np.zeros((N,N))) for q in [np.eye(2)]+[1j*s for s in pauli])
        for i in range(N):
            for j in range(N):
                e=np.zeros((N,N),complex); e[i,j]=1
                basis.extend(rep(N,0,np.zeros((2,2)),z*e) for z in (1,1j))
        for A in basis+[D,gamma]:
            close(Z@A,A@Z)
        close(J@Z.conj()@J,Z)
        close(D,D.conj().T)
        close(J@D.conj()@J,D)
        close(gamma@D+D@gamma,0)
    checks.append('unequal_positive_weights_commute_with_full_representation_D_J_grading')
    lifted_rows=[]
    for N in (3,5):
        for r in (F(1),F(2),F(7,3)):
            h=F(1,2); v=-r*h/N
            Z,J,_,_=structures(N,1,float(r))
            A=raw_abelian(N,float(h),float(v))
            close(np.trace(Z@A),4*(float(r*h+N*v)))
            B=A-J@A.conj()@J
            # Particle charges R,R,L,L. Conjugate right fields when computing anomalies.
            close(np.diag(B)[:4],[0,-2*float(h),-float(h),-float(h)])
            expected=np.repeat([float(h-v),float(-h-v),float(-v),float(-v)],N)
            close(np.diag(B)[8:8+4*N],expected)
            # Trace of lifted B vanishes identically; it is NOT unimodularity on A.
            close(np.trace(Z@B),0)
            lifted_rows.append({'N':N,'r':str(r),'q':str(-v),'h':str(h)})
    checks.append('raw_weighted_unimodularity_and_actual_J_lift')
    for N in range(3,20):
        for h in (F(1,2),F(3,7),F(-2)):
            for r in (F(1),F(2),F(7,3),F(1,5)):
                a=anomalies(N,h,r*h/N)
                assert a['weak']==(r-1)*h
                assert a['cubic']==6*h**3*(1-r)
                assert a['mixed_gravity']==a['color']==0
                assert (a['weak']==a['cubic']==0)==(r==1)
    checks.append('unweighted_physical_anomalies_force_equal_weights_for_nonzero_h')
    # Dropping the nonzero-U(1) hypothesis is a genuine degenerate branch.
    assert all(x==0 for x in anomalies(3,F(0),F(0)).values())
    # Three identical generations multiply, not cancel, this anomaly.
    assert 3*anomalies(3,F(1,2),F(1,3))['weak']==F(3,2)
    checks.append('zero_h_branch_and_generation_multiplicity_controls')
    for N in (3,5):
        h=F(1,2);q=h/N
        assert all(a==0 for a in anomalies(N,h,q).values())
        Z,J,_,_=structures(N,1,2)
        A=raw_abelian(N,float(h),float(-q))
        close(np.trace(A),0)
        close(np.trace(Z@A),2)  # ordinary condition works; same-Z one does not.
    checks.append('separate_unimodularity_and_action_traces_avoid_the_obstruction')
    moments=[]
    for N in (3,5):
        for wq,wl in ((1.0,1.0),(1.0,2.0)):
            for t in (0.0,0.3,1.0):
                Z,J,_,D=structures(N,wq,wl,t)
                a=wl*(2**2+3**2)+N*wq*(5**2+7**2)
                b=wl*(2**4+3**4)+N*wq*(5**4+7**4)
                moment0=8*(N*wq+wl)
                moment2=4*a*t*t+2*wl*11**2
                moment4=4*b*t**4+8*wl*11**2*2**2*t*t+2*wl*11**4
                close(np.trace(Z),moment0)
                close(np.trace(Z@D@D),moment2)
                # Larger fourth moments need relative tolerance for float t.
                assert np.isclose(np.trace(Z@np.linalg.matrix_power(D,4)),moment4,rtol=1e-13,atol=1e-8)
            moments.append({'N':N,'wq':wq,'wl':wl,'TrZ':moment0,'aZ':a,'bZ':b})
    checks.append('same_weights_change_identity_Higgs_and_Majorana_spectral_moments')
    previous=json.loads((HERE/'joint_gauge_matter_constraints_results.json').read_text('utf8'))
    r=previous['historical_running']['positive_weight_ratio_lepton_over_quark']
    gap=(r-1)/2
    assert r>2 and gap>0.6
    checks.append('round531_historical_weight_cannot_satisfy_same_Z_anomaly_condition')
    names=['research_note_531.md','joint_gauge_matter_constraints_results.json',
           'unified_physics_condition_ledger_531.md']
    return dict(round=532,tests_run=len(checks),failures=0,errors=0,checks=checks,
        exact_same_trace_counterexample={k:str(v) for k,v in anomalies(3,F(1,2),F(1,3)).items()},
        lifted_examples=lifted_rows,spectral_moment_examples=moments,
        round531_historical_r=r,round531_historical_weak_anomaly_if_same_trace=round(gap,10),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope=dict(same_weight_for_unimodularity_is_extra_hypothesis=True,
            unequal_weights_compatible_with_listed_finite_operators=True,
            equal_weights_forced_by_this_added_hypothesis_and_anomaly=True,
            all_weighted_actions_ruled_out=False,full_spectral_or_quantum_gravity_model=False,
            three_spatial_dimensions_or_standard_model_derived=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
