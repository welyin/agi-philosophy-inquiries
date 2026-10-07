"""839: consistent finite CAR inclusions, even restriction and support cases.

The matrix tower is a diagnostic of the analytic original-region construction,
not a discretization or a numerical solution of the curved background.
"""
from pathlib import Path
import argparse,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'nested_region_entropy_split_results.json'
from gaussian_reference_entropy_split import kron_all,state,entropy,relative

def marginal(rho,keep,total):
    a=2**keep;b=2**(total-keep)
    return np.einsum('aibi->ab',rho.reshape(a,b,a,b))
def remainder(rho,drop,total):
    a=2**drop;b=2**(total-drop)
    return np.einsum('aiaj->ij',rho.reshape(a,b,a,b))
def extended_relative(a,b):
    ev,u=np.linalg.eigh(b);support=ev>1e-13
    outside=u[:,~support]
    if outside.size and float(np.trace(outside.conj().T@a@outside).real)>1e-12:return math.inf
    aa=u[:,support].conj().T@a@u[:,support];bb=np.diag(ev[support])
    ew=np.linalg.eigvalsh(aa);positive=ew>1e-13
    return float(np.dot(ew[positive],np.log(ew[positive]))-np.trace(aa@np.diag(np.log(ev[support]))).real)
def run():
    total=6;record=3;bias=.6
    unit=np.eye(2);x=np.array([[0,1],[1,0]],complex);y=np.array([[0,-1j],[1j,0]]);z=np.diag([1.,-1.])
    gamma=[]
    for j in range(total):
        for a in (x,y):gamma.append(kron_all([z]*j+[a]+[unit]*(total-1-j)))
    rng=np.random.default_rng(83901);a=rng.normal(size=(12,12));a=(a-a.T)*.08
    k=sum(1j*a[i,j]*gamma[i]@gamma[j] for i in range(12) for j in range(i+1,12))
    ref=state(k);nu=remainder(ref,record,total)
    sigma=(np.eye(8)+bias*kron_all([z]*3))/8
    rho=np.kron(sigma,nu);g=np.kron(np.eye(8)/8,nu)
    c=float(np.log(8)-entropy(sigma));rows=[];previous=(-1.,-1.)
    for modes in range(3,7):
        rn=marginal(rho,modes,total);gn=marginal(g,modes,total);fn=marginal(ref,modes,total)
        d=relative(rn,fn);b=relative(gn,fn)
        assert abs(d-b-c)<2e-12
        assert d>=previous[0]-2e-12 and b>=previous[1]-2e-12
        parity=np.diag(kron_all([z]*modes)).real
        # The even algebra is a direct sum of parity sectors. Its relative
        # entropy includes central probability weights; no block normalization.
        de=0.
        for sign in (1,-1):
            ids=np.flatnonzero(parity==sign)
            de+=relative(rn[np.ix_(ids,ids)],fn[np.ix_(ids,ids)])
        assert abs(de-d)<2e-12
        if modes<6:
            consistent=marginal(marginal(ref,modes+1,total),modes,modes+1)
            assert np.linalg.norm(consistent-fn)<1e-13
        rows.append(dict(modes=modes,full_relative_entropy=d,
                         common_gaussian_cost=b,constant_record_cost=c,
                         split_residual=abs(d-b-c),even_algebra_residual=abs(de-d),
                         reference_minimum_eigenvalue=float(np.linalg.eigvalsh(fn).min())))
        previous=d,b
    # Nonfaithful common exterior with correct support: no artificial noise.
    vacuum=np.diag([1.,0.]);cold=np.diag([.7,.3])
    ref8=kron_all([cold,np.eye(2)/2,np.eye(2)/2])
    rs=np.kron(sigma,vacuum);gs=np.kron(np.eye(8)/8,vacuum);fs=np.kron(ref8,vacuum)
    ds=extended_relative(rs,fs);bs=extended_relative(gs,fs)
    assert math.isfinite(ds) and abs(ds-bs-c)<1e-12
    bad=np.zeros((16,16));bad[0,0]=1
    assert math.isinf(extended_relative(rs,bad)) and math.isinf(extended_relative(gs,bad))
    # A fixed finite number of modes does not bound the Gaussian background
    # cost uniformly as the reference approaches a forbidden support.
    approaching=[]
    for probability in (1e-2,1e-4,1e-8):
        f=kron_all([np.diag([1-probability,probability]),unit/2,unit/2])
        b=relative(np.eye(8)/8,f)
        expected=-.5*np.log(4*probability*(1-probability))
        assert abs(b-expected)<1e-12
        approaching.append(dict(reference_occupation=probability,background_cost=b))
    return dict(round=839,all_checks_passed=True,fresh_test_groups=1,
                nested_CAR_tower_rows=rows,same_original_pair_restricted_at_every_level=True,
                parity_center_probability_weights_retained=True,
                nonfaithful_common_support=dict(full_relative_entropy=ds,gaussian_cost=bs,
                                                record_cost=c,support_dimension=8),
                support_failure=dict(full_relative_entropy='infinite',gaussian_cost='infinite',
                                     both_share_support_failure=True),
                reference_boundary_examples=approaching,
                boundary_examples_are_not_a_single_regional_tower=True,
                original_continuum_embedding_and_limit_are_analytic=True,
                original_regional_relative_entropy_finite_proven=False,
                area_or_Einstein_response_computed=False,
                scope='Literal nested restrictions preserve the entropy split and its finite record term, including the even algebra with central weights. Nonfaithful references are treated by common support or infinity. The original-region result uses the analytic dense Cauchy-mode construction and Araki monotone convergence, not an inference from this finite numerical tower.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
