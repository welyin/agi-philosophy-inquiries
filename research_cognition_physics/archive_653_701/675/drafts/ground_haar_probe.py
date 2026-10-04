"""675 entry: full-group/full-S9 pilot for the fixed-lambda ground-kernel limit.

This is exploratory integration, not a proof of its sign or of RP.
Original 16 channels, original fields, no fixed E and no weight division.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_gauss_boundary_functional as model
old=model.old
TARGET=HERE/'ground_haar_probe_results.json'

def pack(z):
    z=np.asarray(z)
    return dict(real=z.real.tolist(),imag=z.imag.tolist())

def haar_su(rng,n):
    z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    q,r=np.linalg.qr(z)
    # Positive R diagonal fixes QR phase, producing Haar U(n).
    d=np.diag(r)
    q=q*(d/np.abs(d))
    q[:,0]*=np.linalg.det(q).conjugate()
    assert old.err(q.conj().T@q-np.eye(n))<3e-14
    assert abs(np.linalg.det(q)-1)<3e-14
    return q

def haar_group(rng):
    return haar_su(rng,3),haar_su(rng,2),np.exp(2j*np.pi*rng.random())

def run():
    seed=675101;samples=8;rng=np.random.default_rng(seed)
    spatial=np.array([[model.prior.rep(*model.prior.group(67515+2*i+x,s))
        for x in range(2)] for i,s in enumerate((.18,.31))])
    phis=model.mass.car.PHI.reshape(2,2,5)
    raw=[];gaps=[];ranks=[];norms=[];nearzeros=[]
    for sample in range(samples):
        temporal=np.array([model.prior.rep(*haar_group(rng)) for _ in range(4)])
        e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
        block=np.zeros((2,2),complex)
        gapblock=np.zeros((2,2));rankblock=np.zeros((2,2),int)
        nullblock=np.zeros((2,2),int)
        for i in range(2):
            for j in range(2):
                links=np.array([np.r_[spatial[i],spatial[j]],temporal])
                u,v,d,h,gap=model.kernel(links)
                mat=model.fixed_matrices(e,np.r_[phis[i],phis[j]])
                data=model.regular(u,v,d,mat,lam=.37)
                block[i,j]=data['weight']
                gapblock[i,j]=gap;rankblock[i,j]=u.shape[1]
                sv=np.linalg.svd(data['N'],compute_uv=False)
                nullblock[i,j]=int(sum(sv<1e-10))
        raw.append(block);gaps.append(gapblock.tolist());ranks.append(rankblock.tolist())
        nearzeros.append(nullblock.tolist())
        norms.append(float(np.linalg.norm(block)))
    raw=np.array(raw)
    # By674 and measure-preserving exchange of the two E fields and inverse
    # seams, this is an unbiased antithetic estimator of the Hermitian integral.
    # It is not pointwise positivity or projection onto the PSD cone.
    paired=(raw+raw.conj().transpose(0,2,1))/2
    mean=paired.mean(axis=0)
    re_stderr=paired.real.std(axis=0,ddof=1)/np.sqrt(samples)
    im_stderr=paired.imag.std(axis=0,ddof=1)/np.sqrt(samples)
    absolute=np.abs(paired)
    sums=absolute.sum(axis=0)
    largest=np.divide(absolute.max(axis=0),sums,out=np.zeros_like(sums),where=sums!=0)
    spread=np.divide(sums**2,(absolute**2).sum(axis=0),out=np.zeros_like(sums),
        where=(absolute**2).sum(axis=0)!=0)
    scale=float(abs(mean).max())
    assert np.isfinite(raw).all() and scale>0
    eigenvalues=np.linalg.eigvalsh(mean/scale)*scale
    deps=('research_note_591.md','research_note_603.md','research_note_643.md',
          'research_note_673.md','research_note_674.md','joint_gauss_boundary_functional.py')
    return dict(date='2026-10-02',entry_round=675,seed=seed,iid_base_samples=samples,
        lambda_fixed=.37,time_parameter_not_identified_with_original_HF=True,
        spatial_seeds=[67515,67516,67517,67518],spatial_scales=[.18,.31],
        original_fields=phis.tolist(),Haar_method='complex Ginibre QR positive R phase; determinant correction to SU3 and SU2; uniform circle; pushforward to original Z6 quotient',
        E_method='independent normalized real Gaussian at each of four sites: full S9',
        auxiliary_E_integrated_not_declared_observable=True,
        raw_matrices=[pack(v) for v in raw],antithetic_mean=pack(mean),
        nominal_standard_error_real=re_stderr.tolist(),nominal_standard_error_imag=im_stderr.tolist(),
        standard_errors_not_confidence_or_sign_certificates=True,
        largest_absolute_sample_share=largest.tolist(),
        absolute_contribution_spread_not_MCMC_ESS=spread.tolist(),
        sample_mean_eigenvalues=eigenvalues.tolist(),sample_matrix_norms=norms,
        Wilson_gaps=gaps,negative_chiral_ranks=ranks,numerical_N_nullities=nearzeros,
        weights_not_divided_or_absolutized=True,no_PSD_projection=True,
        no_positive_or_negative_integral_claim=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps},
        scope='Pilot for full double-Haar and full auxiliary sphere average at fixed two physical half configurations; mathematical ground-projection branch only. No evaluated original Hb heat kernel, no sign certificate, no fermionic RP, no original HF or continuum identification.')

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_round','iid_base_samples','sample_mean_eigenvalues','largest_absolute_sample_share','no_positive_or_negative_integral_claim')}))

