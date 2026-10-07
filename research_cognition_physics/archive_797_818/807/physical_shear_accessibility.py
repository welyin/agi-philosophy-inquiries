"""807: symbol checks supporting the analytic physical-shear accessibility proof.

This does NOT solve the original coupled PDE or evaluate its switched response.
The exact constrained solution and fixed-state accessibility are proved in the note.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'802'))
import original_bff_vertex as original
TARGET=HERE/'physical_shear_accessibility_results.json'

def tf(a):return a-np.trace(a)*np.eye(3)/3
def basis():
    out=[np.diag([1.,-1.,0.])/np.sqrt(2),np.diag([1.,1.,-2.])/np.sqrt(6)]
    for i,j in ((0,1),(0,2),(1,2)):
        a=np.zeros((3,3));a[i,j]=a[j,i]=1/np.sqrt(2);out.append(a)
    return np.array(out)

def run():
    rng=np.random.default_rng(807)
    es=basis();I=np.eye(3);errors=[];ranks=[]
    for _ in range(12):
        n=rng.normal(size=3);n/=np.linalg.norm(n)
        B=np.column_stack([np.einsum('aij,ij->a',es,tf(np.outer(n,x)+np.outer(x,n))) for x in I])
        gram=B.T@B
        q=np.eye(5)-B@np.linalg.solve(gram,B.T)
        expected=2*I+(2/3)*np.outer(n,n)
        errors.extend([np.max(abs(gram-expected)),np.max(abs(q@B)),np.max(abs(q@q-q))])
        ranks.append(int(np.linalg.matrix_rank(q,tol=1e-12)))
    assert max(errors)<1e-12 and ranks==[2]*12
    idx=[30,31]
    G=[a[np.ix_(idx,idx)] for a in original.original.GAMMA]
    def gamma(n):return sum((x*g for x,g in zip(n,G)),np.zeros((2,2),complex))
    def cross(k,n):
        p=(np.eye(2)+gamma(n))/2
        return p@gamma(k@n)@(np.eye(2)-p)/2
    directions=[*I,*(-I),*( (I[i]+I[j])/np.sqrt(2) for i,j in ((0,1),(0,2),(1,2)) )]
    fullbasis=[*es,I/np.sqrt(3)]
    mat=np.column_stack([np.concatenate([cross(k,n).ravel() for n in directions]) for k in fullbasis])
    sv=np.linalg.svd(mat,compute_uv=False)
    assert np.linalg.matrix_rank(mat,tol=1e-12)==5 and sv[-1]<1e-13
    assert np.linalg.norm(mat[:,-1])<1e-13
    tfsv=np.linalg.svd(mat[:,:5],compute_uv=False)
    assert tfsv[-1]>.1
    hermitian_residual=0.;sample_min=1e100
    for _ in range(12):
        coeff=rng.normal(size=5)+1j*rng.normal(size=5)
        k=np.einsum('a,aij->ij',coeff,es)
        responses=[]
        for n in directions:
            a=cross(k,n);ell=(a-a.conj().T)/1j
            hermitian_residual=max(hermitian_residual,float(np.max(abs(ell-ell.conj().T))))
            responses.append(np.linalg.norm(ell))
        sample_min=min(sample_min,float(max(responses)/np.linalg.norm(k)))
        assert max(responses)>1e-10
    invisible=np.diag([1.,-1.,0.])
    assert np.linalg.norm(cross(invisible,I[2]))==0
    detected=max(float(np.linalg.norm(cross(invisible,n))) for n in directions)
    assert detected>.1
    return dict(round=807,all_checks_passed=True,
        calibration_scope='Conformal-Killing symbol quotient and original sterile frequency cross; not a coupled PDE solution.',
        original_Nambu_dimension=64,sterile_indices=idx,
        spatial_symbol_trials=12,TT_quotient_ranks=ranks,
        max_projection_and_ellipticity_residual=float(max(errors)),
        complex_symmetric_tensor_dimension=6,finite_direction_count=9,
        cross_map_complex_rank=5,cross_map_singular_values=sv.tolist(),
        tracefree_min_singular_value=float(tfsv[-1]),
        conformal_null_residual=float(np.linalg.norm(mat[:,-1])),
        complex_tensor_trials=12,hermiticity_residual=hermitian_residual,
        minimum_sample_max_response_per_tensor_norm=sample_min,
        single_direction_invisible_tensor_other_direction_norm=detected,
        actual_boson_accessibility_proof_in_note=True,
        original_PDE_numerically_solved=False,switched_response_nonzero_proven=False,
        new_numbered_test_groups=1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
