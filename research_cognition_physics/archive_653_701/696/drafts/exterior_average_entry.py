"""696 entry only: exact signed sector of the original local S9 pairing tensor.

Integer arithmetic certifies the exterior-square spectrum. This is a specific
realignment of the original local tensor, not the physical averaged RP form.
"""
import argparse
import hashlib
import itertools
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_gauss_boundary_functional as base
TARGET=HERE/'exterior_average_entry_results.json'

def integer(a):
    r=np.rint(a.real).astype(np.int64);i=np.rint(a.imag).astype(np.int64)
    assert np.array_equal(a,r+1j*i)
    return r,i

def mul(a,b):
    # Every multiplication here has length<=120 and components<=32.
    # Bound proves int64 cannot overflow; no floating product decides a claim.
    bound=2*a[0].shape[1]*max(int(abs(x).max()) for x in a)*max(int(abs(x).max()) for x in b)
    assert bound<2**62
    return a[0]@b[0]-a[1]@b[1],a[0]@b[1]+a[1]@b[0]

def adj(a):return a[0].T,-a[1].T
def eq(a,b):return all(np.array_equal(x,y) for x,y in zip(a,b))
def ident(n):return np.eye(n,dtype=np.int64),np.zeros((n,n),dtype=np.int64)

def wedge2(a):
    pairs=list(itertools.combinations(range(a[0].shape[0]),2))
    out=[np.zeros((len(pairs),len(pairs)),dtype=np.int64) for _ in range(2)]
    def scalar(i,k,j,l):
        x,y=int(a[0][i,k]),int(a[1][i,k]);z,w=int(a[0][j,l]),int(a[1][j,l])
        return x*z-y*w,x*w+y*z
    for p,(i,j) in enumerate(pairs):
        for q,(k,l) in enumerate(pairs):
            x=scalar(i,k,j,l);y=scalar(i,l,j,k)
            out[0][p,q]=x[0]-y[0];out[1][p,q]=x[1]-y[1]
    return tuple(out)

def run():
    ts=[integer(t) for t in base.internal.T]
    us=[mul(adj(ts[0]),t) for t in ts]
    assert eq(us[0],ident(16))
    for a in us:assert eq(mul(adj(a),a),ident(16))
    for a in us[1:]:assert eq(adj(a),(-a[0],-a[1]))
    for i in range(1,10):
        for j in range(1,10):
            ab=mul(us[i],us[j]);ba=mul(us[j],us[i])
            assert eq((ab[0]+ba[0],ab[1]+ba[1]),(-2*int(i==j)*ident(16)[0],ident(16)[1]))
    wedges=[wedge2(a) for a in us]
    a=tuple(sum((w[k] for w in wedges),np.zeros((120,120),dtype=np.int64)) for k in range(2))
    assert eq(a,adj(a))
    square=mul(a,a);assert eq(square,(16*ident(120)[0],ident(120)[1]))
    trace=int(np.trace(a[0]));assert trace==192 and int(np.trace(a[1]))==0
    positive=(120+trace//4)//2;negative=120-positive
    assert (negative,positive)==(36,84)
    # Exact negative eigenvector w=(4I-A)e_j; this retains the actual basis.
    j=next(j for j in range(120) if a[0][j,j]!=4)
    w=(4*ident(120)[0][:,j:j+1]-a[0][:,j:j+1],-a[1][:,j:j+1])
    assert eq(mul(a,w),(-4*w[0],-4*w[1]))
    norm=int(mul(adj(w),w)[0][0,0]);assert norm>0
    # Normalized measure on S9: E_a E_b=delta_ab/10. Thus A/10 is
    # exactly the integrated exterior-square matrix, not numerical quadrature.
    # The original spin pairing splits into two disjoint16+16 blocks.
    # B is transpose-antisymmetric, not anti-Hermitian in this convention.
    b=integer(base.internal.B);assert eq((b[0].T,b[1].T),(-b[0],-b[1]))
    assert eq(mul(adj(b),b),ident(4))
    support=[]
    for i in range(4):
        for j2 in range(i+1,4):
            if b[0][i,j2] or b[1][i,j2]:support.append([i,j2,int(b[0][i,j2]),int(b[1][i,j2])])
    assert len(support)==2 and len({x for row in support for x in row[:2]})==4
    deps=('research_note_615.md','research_note_657.md','research_note_675.md','research_note_695.md',
          'joint_gauss_boundary_functional.py')
    return dict(date='2026-10-02',entry_round=696,new_formal_round=False,
        exact_tensor=dict(original_internal_dimension=16,sphere_dimension=9,
            sphere_second_moment='delta_ab/10',exterior_square_dimension=120,
            scaled_Hermitian_matrix_square=16,scaled_trace=trace,
            eigenvalues=['-2/5','2/5'],multiplicities=[negative,positive],
            negative_witness_real=w[0][:,0].tolist(),negative_witness_imag=w[1][:,0].tolist(),
            negative_witness_squared_norm=norm,negative_witness_Rayleigh_quotient='-2/5',
            spin_pair_support=support,
            actual64mode_tensor_has_two16_by16_pairings=True),
        scope=dict(specific_local_exterior_realignment_not_positive=True,
            all_original_local_S9_moments_retained=True,
            not_a_physical_state_or_original_RP_counterexample=True,
            full_auxiliary_gauge_Hb_average_undecided=True,
            no_HF_continuum_or_GR_completion=True,old_space_interfaces_inherited=True),
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in deps},
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(entry_round=696,all_checks_passed=True,spectrum=r['exact_tensor']['eigenvalues'],
        multiplicities=r['exact_tensor']['multiplicities'])))
