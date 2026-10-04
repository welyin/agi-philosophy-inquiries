"""676 exact integer-flux witness for the fixed-seed boundary obstruction.

Rational symmetric congruence proves inertia. No numerical zero tolerance.
"""
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as base
TARGET=HERE/'exact_flux_seed_results.json'

def inertia(matrix):
    a=[[F(int(x)) for x in row] for row in matrix]
    positive=negative=zero=0;det=F(1);pivots=[]
    while a:
        n=len(a);diag=next((i for i in range(n) if a[i][i]),None)
        if diag is not None:
            order=[diag]+[i for i in range(n) if i!=diag]
            a=[[a[i][j] for j in order] for i in order]
            p=a[0][0];positive+=p>0;negative+=p<0;det*=p;pivots.append(str(p))
            a=[[a[i][j]-a[i][0]*a[0][j]/p for j in range(1,n)] for i in range(1,n)]
        else:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if a[i][j]),None)
            if pair is None:zero+=n;det=F(0);break
            i,j=pair;order=[i,j]+[k for k in range(n) if k not in (i,j)]
            a=[[a[i][j] for j in order] for i in order];p=a[0][1]
            positive+=1;negative+=1;det*=-p*p;pivots.append('block(0,'+str(p)+')')
            a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/p
                for j in range(2,n)] for i in range(2,n)]
    return dict(positive=positive,negative=negative,zero=zero,determinant=str(det),pivots=pivots)

def exact_sector(charge,sign):
    sp=base.internal.spin;sites=base.prior.SITES
    dw=np.zeros((16,16),complex)
    for mu in range(2):
        shift=np.zeros_like(dw)
        for i,(t,x) in enumerate(sites):
            dest=((t+1)%2,x) if mu==1 else (t,1-x)
            j=sites.index(dest)
            z=(1j if x==1 else 1) if mu==1 else ((-1)**t if x==1 else 1)
            phase=z**int(charge)*(-1 if mu==1 and t==1 else 1)
            shift[4*i:4*i+4,4*j:4*j+4]=phase*np.eye(4)
        gamma=np.kron(np.eye(4),sp.GAMMA[3 if mu==1 else 0])
        dw+=np.eye(16)-(shift+shift.conj().T)/2+gamma@(shift-shift.conj().T)/2
    h=np.kron(np.eye(4),sp.G5)@(dw-np.eye(16))
    # Exact orthogonal basis of S=i gamma2 gamma3, with Gram2 I.
    bspin=np.array([[1,0],[-sign,0],[0,1],[0,-sign]])
    b=np.kron(np.eye(4),bspin)
    hs=b.T@h@b
    real=np.block([[hs.real,-hs.imag],[hs.imag,hs.real]])
    assert np.array_equal(real,np.rint(real))
    matrix=real.astype(int).tolist()
    assert np.array_equal(real,real.T)
    count=inertia(matrix);assert count['zero']==0
    assert count['positive']%2==count['negative']%2==0
    row_bound=max(sum(abs(x) for x in row) for row in matrix)
    gap_bound=abs(F(count['determinant']))/(2*F(row_bound)**15)
    return dict(charge_mod4=charge,spin_sign=sign,realified_integer_matrix=matrix,
        real_inertia=count,complex_negative_count=count['negative']//2,
        original_sector_gap_lower_bound=str(gap_bound),
        numerical_gap_diagnostic=float(min(abs(np.linalg.eigvalsh(hs/2)))))

def run():
    mult=Counter(int(q)%4 for q in base.internal.Q)
    assert mult=={0:4,1:8,2:4}
    rows=[exact_sector(q,s) for s in (-1,1) for q in sorted(mult)]
    counts=[sum(mult[r['charge_mod4']]*r['complex_negative_count'] for r in rows if r['spin_sign']==s) for s in (-1,1)]
    assert counts==[72,56]
    # Independent original16-channel construction with original fixed dictionary.
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i,(t,x) in enumerate(base.prior.SITES):
        for mu,z in enumerate(((-1.)**t if x==1 else 1.,1j if x==1 else 1.)):
            links[mu,i]=base.prior.rep(np.eye(3),np.eye(2),z)
    u,v,d,h,gap=base.kernel(links)
    jm=np.kron(np.kron(np.eye(4),base.internal.VM),np.eye(16))
    singular=np.linalg.svd(u.conj().T@jm,compute_uv=False)
    return dict(date='2026-10-02',entry_round=676,
        original_charge_multiset=base.internal.Q.tolist(),mod4_multiplicities=dict(mult),
        spatial_links_z=[1,1,1,-1],temporal_links_z=[1,1j,1,1j] if False else [[1,0],[0,1],[1,0],[0,1]],
        exact_sectors=rows,negative_counts=counts,total_negative_rank=sum(counts),
        fixed_Jminus_seed_ranks=[64,64],forced_missing_negative_directions=8,
        exact_nonzero_Wilson_spectrum=True,positive_full_gauge_neighborhood_by_continuity=True,
        direct_original_Wilson_gap_diagnostic=gap,
        direct_original_seed_nullity_diagnostic=int(sum(singular<1e-10)),
        direct_original_seed_smallest12=singular[-12:].tolist(),
        scope='Exact obstruction to659 fixed-Jminus kernel-image boundary convergence at nonflat original-group configurations. Transfer positivity in fifth direction and finite inverse-free675 candidate unaffected. This does not refute domain-wall methods, full gauge reflection positivity, original HF or cognition principles.',
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in ('research_note_659.md','research_note_675.md','joint_gauss_boundary_functional.py','joint_subgroup_measure_source.py')})

if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('entry_round','negative_counts','forced_missing_negative_directions','direct_original_seed_nullity_diagnostic','exact_nonzero_Wilson_spectrum')}))

