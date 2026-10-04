"""Numerical full-measure diagnostic with explicit positive proposal densities.
Independent pilot fixes a control coefficient. Standard errors are NOT a
deterministic integration certificate or a theorem about physical RP.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np
import full_average_probe as base
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
TARGET=HERE/'importance_full_average_results.json'
CENTERS=np.array([[2*np.pi*k/3,2*np.pi*k/3,np.pi*k,np.pi*k/3] for k in range(6)])
KAPPA=np.array([1.,1.,1.,15.])

def reference():
    data=json.loads((ARCHIVE/'joint_gauss_support_marginal_results.json').read_text('utf8'))
    rows=data['exact_projected_support']['rows']
    exact=sum(r['Gauss_multiplicity']*F(r['mu'])**2 for r in rows)/2**64
    old=next(r for r in data['CAR_marginal_certificate']['rows'] if r['N']==2) if 'CAR_marginal_certificate' in data else None
    if old:assert exact==F(old['raw_projected_partition'])
    return exact

def samples(rng,n):
    uniform=rng.uniform(0,2*np.pi,(n,4));component=rng.integers(6,size=n)
    vm=rng.vonmises(0,KAPPA,size=(n,4))+CENTERS[component]
    choice=rng.random(n)<.9;angles=np.where(choice[:,None],vm,uniform)
    offsets=angles[:,None,:]-CENTERS[None,:,:]
    logp=np.sum(KAPPA*np.cos(offsets)-np.log(np.i0(KAPPA)),axis=2)
    pg=.1+.9*np.mean(np.exp(logp),axis=1)
    h,haar=base.torus(angles)
    e=rng.normal(size=(n,10));e/=np.linalg.norm(e,axis=1)[:,None]
    z=rng.beta(12.5,4.5,size=n);t=2*z-1
    v=rng.normal(size=(n,10));v-=np.sum(v*e,axis=1)[:,None]*e;v/=np.linalg.norm(v,axis=1)[:,None]
    correlated=t[:,None]*e+np.sqrt(1-t*t)[:,None]*v
    ordinary=rng.normal(size=(n,10));ordinary/=np.linalg.norm(ordinary,axis=1)[:,None]
    choice=rng.random(n)<.9;f=np.where(choice[:,None],correlated,ordinary)
    z=(1+np.sum(e*f,axis=1))/2
    norm=F(1)
    for k in range(8):norm*=F(9,2)+k;norm/=9+k
    pe=.1+.9*z**8/float(norm)
    fields=np.stack((e,f),axis=1)
    vals=np.stack([base.values(h,fields,*flags) for flags in ((False,False),(True,False),(True,True))],axis=1)
    vals*=haar[:,None]/(pg*pe)[:,None]
    assert np.max(abs(vals.imag))<1e-15
    return vals.real

def accumulate(n,seed,batch=512):
    rng=np.random.default_rng(seed);total=np.zeros(3);second=np.zeros((3,3));minimum=np.full(3,np.inf);maximum=np.full(3,-np.inf)
    for start in range(0,n,batch):
        x=samples(rng,min(batch,n-start));total+=x.sum(axis=0);second+=x.T@x
        minimum=np.minimum(minimum,x.min(axis=0));maximum=np.maximum(maximum,x.max(axis=0))
    mean=total/n;cov=(second-n*np.outer(mean,mean))/(n-1)
    return mean,cov,minimum,maximum

def run():
    original_error=base.crosscheck();ref=reference()
    pilot_n=16384;pilot=accumulate(pilot_n,69761)
    beta=pilot[1][0,:]/pilot[1][0,0];beta[0]=1
    count=524288;mean,cov,mi,ma=accumulate(count,69762)
    transform=np.eye(3);transform[:,0]-=beta
    adjusted=mean-beta*(mean[0]-float(ref));adjusted_cov=transform@cov@transform.T/count
    coef=np.array([1.,-2/1000,1/10**6]);test=float(coef@adjusted);se=float(np.sqrt(coef@adjusted_cov@coef))
    raw_se=np.sqrt(np.diag(cov)/count)
    return dict(date='2026-10-02',entry_round=697,new_formal_round=False,
        full_original_group_and_two_S9_measures=True,original_full_matrix_relative_error=original_error,
        inherited_exact_B00=str(ref),inherited_exact_B00_decimal=float(ref),
        independent_pilot=dict(samples=pilot_n,seed=69761,control_coefficients=beta.tolist()),
        main=dict(samples=count,seed=69762,raw_mean=mean.tolist(),raw_standard_error=raw_se.tolist(),
            controlled_mean=adjusted.tolist(),controlled_covariance_of_mean=adjusted_cov.tolist(),
            sampled_range=[mi.tolist(),ma.tolist()]),
        fixed_rational_test=dict(vector=[1,'-1/1000'],estimated_value=test,estimated_standard_error=se,
            standard_errors_below_zero=-test/se,
            negative_deterministic_certificate=False,confidence_interval_not_rigorously_coverage_certified=True),
        proposals=dict(torus='0.1 uniform +0.9 six-Z6-centre product-vonMises mixture',kappa=KAPPA.tolist(),
            spheres='uniform E;0.1 uniform F+0.9 conditional((1+E.F)/2)^8',
            both_densities_bounded_below_by=.1,self_normalized_importance_not_used=True),
        limitation='Sample standard errors do not prove exact sign; original Hb and finite-tau physical RP remain separate.',
        dependency_hashes={p:hashlib.sha256((ARCHIVE/p).read_bytes()).hexdigest() for p in
            ('joint_gauss_support_marginal_results.json','round697_drafts/full_average_probe.py','round697_drafts/center_holonomy_entry.py')})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['fixed_rational_test']))
