"""All real radial responses at fixed diagonal Higgs pseudoscalars."""
import argparse
import cmath
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
EDGES=[(0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2)]
C=np.array([[0,0,0,0,0,-1,1],[0,0,-1,1,-1,0,1]],dtype=float)
PAIRS=[(0,1),(1,0),(0,2),(2,0),(1,2),(2,1)]


def read(path):
    return json.loads((HERE/path).read_bytes())


def projected_norms(g,yb):
    return [float((abs(g[i,j])**2+abs(g[j,i])**2)/yb**2)
            for i,j in [(0,1),(0,2),(1,2)]]


def latched(slopes,derivatives,m):
    dg=np.einsum('a,aij->ij',slopes,derivatives)
    # Add a universal real slope to set the scalar bb coupling to its reference.
    common=-dg[2,2].real/m[2]
    locked=slopes+common
    g=np.diag(m)+np.einsum('a,aij->ij',locked,derivatives)
    assert abs(g[2,2].real/m[2]-1)<1e-12
    p=np.diag(g).imag/m
    rates=projected_norms(g,m[2])
    return {'slopes_at_scalar_b_equal_one':locked.tolist(),'ds_bd_bs_over_SM_bb':rates,
            'total_FC_over_SM_bb':sum(rates),'p':p.tolist()}


def calculate():
    cfg=read('inputs.json')
    ref=read(cfg['reconstruction'])['reconstruction']
    a,b,c,e,f=ref['normalized_real_edges_A_B_C_E_F']
    r,d=ref['two_real_coefficients']
    vals=[a,b,c,e,f,r,d*cmath.exp(1j*np.radians(ref['eta_example_deg']))]
    terms=np.zeros((7,3,3),dtype=complex)
    for j,((s,t),amp) in enumerate(zip(EDGES,vals,strict=True)):
        terms[j,s,t]=amp
    Y=terms.sum(axis=0)
    ul,m,vh=np.linalg.svd(Y)
    idx=np.argsort(m)
    m,ul,vr=m[idx],ul[:,idx],vh.conj().T[:,idx]
    deriv=np.array([ul.conj().T@t@vr for t in terms])
    def real_operator(pairs):
        # Scale each row by yb; a squared Euclidean norm equals a tree width ratio.
        mat=np.array([deriv[:,i,j]/m[2] for i,j in pairs])
        return np.concatenate([mat.real,mat.imag])
    A=real_operator(PAIRS)
    K=real_operator(PAIRS[:2])
    U,sv,V=np.linalg.svd(C,full_matrices=True)
    inv=V[:2].T@np.diag(1/sv)@U.T
    N=V[2:].T
    assert np.linalg.matrix_rank(A)==6 and np.linalg.matrix_rank(K@N)==4
    assert np.linalg.norm(A@np.ones(7))<1e-12
    # The factors below are shared by every choice of the two contrasts.
    least=inv-N@np.linalg.pinv(A@N,rcond=1e-12)@A@inv
    cancel=inv-N@np.linalg.pinv(K@N,rcond=1e-12)@K@inv
    B=A@least
    Q=B.T@B
    np.testing.assert_allclose((A@N).T@B,0,atol=1e-15)
    np.testing.assert_allclose(K@cancel,0,atol=1e-14)
    assert min(np.linalg.eigvalsh(Q))>0
    old=read(cfg['radial'])['sample']
    scenarios={}
    for model in ['phase_path_dim6','real_path_dim6','both_nonlinear','universal_linear']:
        z=old[model]['two_contrasts']
        minimum=latched(least@z,deriv,m)
        kaon_free=latched(cancel@z,deriv,m)
        for solution in [minimum,kaon_free]:
            np.testing.assert_allclose(solution['p'],old[model]['diagonal_p'],atol=2e-10)
        assert abs(minimum['total_FC_over_SM_bb']-z@Q@z)<1e-14
        assert kaon_free['ds_bd_bs_over_SM_bb'][0]<1e-23
        bottom=read('../1386/inputs.json')['yukawa']['down'][2]
        yukawa_low=(bottom[0]-2*bottom[1])/np.sqrt(2)
        kaon_free['physical_bs_coupling_lower_from_pair_norm']=yukawa_low*np.sqrt(kaon_free['ds_bd_bs_over_SM_bb'][2]/2)
        scenarios[model]={'global_minimum':minimum,'cancel_both_kaon_vertices':kaon_free}
    return {'quadratic_form':Q.tolist(),'eigenvalues':np.linalg.eigvalsh(Q).tolist(),
            'single_kaon_rank_on_invisible_space':4,'full_offdiag_rank_on_all_slopes':6,'scenarios':scenarios}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):
        return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return 0. if abs(x)<1e-16 else float(f'{x:.10g}')
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    data=clean(calculate())
    out=HERE/'results.json'
    if args.check:
        assert data==json.loads(out.read_bytes())
        print('PASS: all real slopes, two independently cancelled kaon vertices, sharp inclusive FC minimum')
    else:
        with out.open('x',encoding='utf-8') as f:
            f.write(json.dumps(data,indent=2)+'\n')
