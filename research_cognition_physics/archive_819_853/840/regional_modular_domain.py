"""840: finite covariance calibration of the original-region domain criterion.

The reference compression tower and fixed record projection are consistent.
No comparison between physical energy and modular energy is assumed.
"""
from pathlib import Path
import argparse,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'regional_modular_domain_results.json'

def spectral(a,f):
    w,u=np.linalg.eigh(a);return (u*f(w))@u.conj().T
def logm(a):return spectral(a,np.log)
def entropy_cov(a):
    w=np.linalg.eigvalsh(a);assert w.min()>0 and w.max()<1
    return float(-.5*np.sum(w*np.log(w)+(1-w)*np.log1p(-w)))
def run():
    rng=np.random.default_rng(840);dim=48;q=20;logd=(q/2)*np.log(2)
    o,_=np.linalg.qr(rng.normal(size=(dim,dim)))
    blocks=np.zeros((dim,dim))
    for k,t in enumerate(np.linspace(.1,.97,dim//2)):
        blocks[2*k,2*k+1]=t;blocks[2*k+1,2*k]=-t
    full=(np.eye(dim)+1j*o@blocks@o.T)/2
    full_h=logm(np.eye(dim)-full)-logm(full)
    full_L=-logm(full)-logm(np.eye(dim)-full)
    zfull=float(np.trace(spectral(full_h,np.abs)[:q,:q]).real)
    ellfull=float(np.trace(full_L[:q,:q]).real)
    assert ellfull<=zfull+2*q*np.log(2)+1e-12
    rows=[];prev=(-1.,-1.,-1.)
    for n in (20,24,28,36,48):
        a=full[:n,:n];qmat=np.diag([1.]*q+[0.]*(n-q));r=np.eye(n)-qmat
        g=qmat/2+r@a@r
        h=logm(np.eye(n)-a)-logm(a)
        z=float(np.trace(spectral(h,np.abs)[:q,:q]).real)
        ell=-logm(a)-logm(np.eye(n)-a)
        jensen_error=float(np.linalg.eigvalsh(full_L[:n,:n]-ell).min())
        assert jensen_error>-1e-12
        delta_s=entropy_cov(g)-entropy_cov(a)
        delta_k=float(.5*np.trace(h@(g-a)).real)
        b=delta_k-delta_s
        independent=float(.5*np.trace(g@(logm(g)-logm(a))+
                           (np.eye(n)-g)@(logm(np.eye(n)-g)-logm(np.eye(n)-a))).real)
        hq=h[:q,:q];aq=a[:q,:q]
        block_k=float(.5*(np.trace(aq@hq)-2*np.trace((a@h)[:q,:q])).real)
        assert abs(block_k-delta_k)<1e-12
        assert abs(b-independent)<1e-12
        assert -1e-12<=delta_s<=2*logd+1e-12
        assert b>=-1e-12 and delta_k<=.75*z+1e-12
        lower=.25*z-q/math.e-2*logd
        assert b>=lower-1e-12
        assert z<=ellfull+1e-12
        assert all(v>=p-1e-12 for v,p in zip((b,delta_s,delta_k),prev))
        prev=b,delta_s,delta_k
        rows.append(dict(self_dual_modes=n,background_relative_entropy=b,
                         entropy_increase=delta_s,modular_energy_increase=delta_k,
                         record_absolute_modular_moment=z,compression_log_domain_bound=ellfull,
                         relative_entropy_lower_bound=lower,relative_entropy_upper_bound=.75*z,
                         independent_relative_entropy_residual=abs(b-independent),
                         self_dual_block_formula_residual=abs(block_k-delta_k),
                         operator_Jensen_minimum_difference=jensen_error))
    # Dimension-two direct factor: proves the lower estimate really permits
    # arbitrarily large cost. Only log-odds are used, avoiding underflow.
    large=[]
    for t in (1.,8.,40.,120.):
        a=1/(1+math.exp(t))
        b=math.log1p(math.exp(-t))+t/2-math.log(2)
        z=2*t;deltak=(.5-a)*t
        deltas=math.log(2)+a*math.log(a)+(1-a)*math.log1p(-a)
        assert abs(b-deltak+deltas)<1e-12
        assert z/4-2/math.e-2*math.log(2)<=b<=.75*z
        large.append(dict(log_odds=t,background_relative_entropy=b,absolute_modular_moment=z))
    # An abstract fixed countable reference: a_j = const*exp(-j),
    # physical test energy E_j=j, modular log-odds h_j=exp(3j).
    # The analytic series gives all polynomial energy moments finite but
    # sum |a_j|^2 h_j divergent. This is not the original spacetime reference.
    c=math.exp(2)-1
    partials=[]
    for cutoff in (4,8,12):
        moment=2*c*sum(math.exp(j) for j in range(1,cutoff+1))
        partials.append(dict(cutoff=cutoff,absolute_modular_moment_partial_sum=moment))
    return dict(round=840,all_checks_passed=True,fresh_test_groups=1,
                original_record_complex_modes=10,self_dual_record_dimension=q,
                literal_compression_tower=rows,
                original_full_reference_record_absolute_moment=zfull,
                original_full_reference_log_sum_moment=ellfull,
                large_modular_cost_calibration=large,
                abstract_energy_domain_separation_partial_sums=partials,
                abstract_example_claimed_to_be_original_Hadamard_background=False,
                exact_region_criterion_is_analytic=True,
                original_selected_modes_log_domain_verified=False,
                physical_H_equals_regional_modular_generator_assumed=False,
                unification_or_area_response_proven=False,
                scope='For the declared trace-reset Gaussian preparation, regional reference-relative entropy is finite exactly when the fixed record modes have finite absolute one-particle modular first moment. The proof uses the finite conditional entropy bound, a two-sided modular-energy estimate, compression Jensen inequality and the 839 increasing-algebra limit. The criterion is not automatically implied by physical Sobolev or energy regularity.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
