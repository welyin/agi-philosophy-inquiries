"""697 entry: actual all-channel centre histories and the remaining full holonomy.
Checks covariance/normalization only, never samples a claimed Haar integral.
"""
import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent;sys.path.insert(0,str(ARCHIVE))
import joint_gauss_boundary_functional as base
TARGET=HERE/'center_holonomy_entry_results.json'

def kernel(first,second,a,c):
    r=base.prior.rep(np.eye(3),-np.eye(2),1);n=128
    ts=np.zeros((n,n),complex)
    ts[:64,:64]=np.kron(np.eye(4),r if first else np.eye(16))
    ts[64:,64:]=np.kron(np.eye(4),r if second else np.eye(16))
    tt=np.zeros_like(ts);tt[:64,64:]=np.kron(np.eye(4),a);tt[64:,:64]=-np.kron(np.eye(4),c)
    gs=np.kron(np.kron(np.eye(2),base.internal.spin.G5),np.eye(16))
    g0=np.kron(np.kron(np.eye(2),base.internal.spin.GAMMA[3]),np.eye(16))
    x=np.eye(n)-(ts+ts.conj().T)/2-(tt+tt.conj().T)/2+g0@(tt-tt.conj().T)/2
    h=gs@x;ev,v=np.linalg.eigh(h);gap=float(min(abs(ev)));assert gap>1e-8
    d=(np.eye(n)+gs@(v*np.sign(ev))@v.conj().T)/2
    return v[:,ev<0],v[:,ev>0],d,gap

def weight(first,second,a,c,e):
    u,v,d,gap=kernel(first,second,a,c)
    mat=base.fixed_matrices(e,np.zeros((2,5)))
    return base.regular(u,v,d,mat,lam=0)['weight'],gap

def weyl_normalization():
    # Constant terms of actual SU3 and SU2 Vandermonde absolute squares.
    def sign(p):return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
    exponent=[(1,0),(0,1),(-1,-1)];delta={}
    for p in itertools.permutations(range(3)):
        key=tuple(sum(p[i]*exponent[i][j] for i in range(3)) for j in range(2))
        delta[key]=delta.get(key,0)+sign(p)
    coeff={}
    for a,x in delta.items():
        for c,y in delta.items():
            key=tuple(i-j for i,j in zip(a,c));coeff[key]=coeff.get(key,0)+x*y
    assert coeff[(0,0)]==6
    weak={2:-1,0:2,-2:-1};assert weak[0]==2
    assert coeff[(0,0)]*weak[0]==12
    # Original Z6 generator is represented trivially, so normalized product
    # Haar pushes forward to the quotient with no extra factor6.
    zeta=np.exp(2j*np.pi/3);r=base.prior.rep(zeta*np.eye(3),-np.eye(2),np.exp(1j*np.pi/3))
    err=float(np.max(abs(r-np.eye(16))));assert err<2e-14
    return dict(SU3_Weyl_square_constant=6,SU2_Weyl_square_constant=2,
        full_Weyl_normalizing_factor=12,torus_rank=4,quotient_generator_identity_error=err,
        no_extra_factor6_for_normalized_quotient_Haar=True)

def run():
    rng=np.random.default_rng(69701);e=rng.normal(size=(2,10));e/=np.linalg.norm(e,axis=1)[:,None]
    a=base.prior.rep(*base.prior.group(69711,.34));c=base.prior.rep(*base.prior.group(69712,.31))
    h=a@c;g=base.prior.rep(*base.prior.group(69713,.27))
    ea=e.copy();ea[1]=base.prior.vector_rotation(a)@e[1]
    eg=np.array([base.prior.vector_rotation(g)@v for v in ea]);rows=[]
    for first,second in ((False,False),(True,False),(False,True),(True,True)):
        f,gap=weight(first,second,a,c,e)
        sewn,gap2=weight(first,second,np.eye(16),h,ea)
        conjugated,gap3=weight(first,second,np.eye(16),g@h@g.conj().T,eg)
        err=float(abs(f-sewn)/max(abs(f),abs(sewn),1e-30))
        err2=float(abs(sewn-conjugated)/max(abs(sewn),abs(conjugated),1e-30))
        assert err<3e-10 and err2<3e-10
        wrong,_=weight(first,second,np.eye(16),h,e)
        rows.append(dict(half_centres=[first,second],Wilson_gap_min=min(gap,gap2,gap3),
            original_full_weight=[float(f.real),float(f.imag)],
            complete_sewing_relative_error=err,simultaneous_conjugation_relative_error=err2,
            omit_E_transport_relative_error=float(abs(f-wrong)/max(abs(f),1e-30))))
    assert max(x['omit_E_transport_relative_error'] for x in rows)>1e-3
    deps=('research_note_673.md','research_note_675.md','research_note_696.md',
          'joint_dynamic_auxiliary_integral_results.json','joint_gauss_boundary_functional.py')
    return dict(date='2026-10-02',entry_round=697,new_formal_round=False,
        original_matrix_sewing=rows,Weyl_normalization=weyl_normalization(),
        scope=dict(original_full_group_and16channels_retained=True,
            class_function_only_after_full_E_average=True,remaining_Haar_integral_not_computed=True,
            fixed_holonomy_samples_not_used_as_average=True,physical_RP_undecided=True,
            old_space_interfaces_inherited=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps},
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(entry_round=697,all_checks_passed=True,
        max_sewing_relative=max(x['complete_sewing_relative_error'] for x in r['original_matrix_sewing']))))
