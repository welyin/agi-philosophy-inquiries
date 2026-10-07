"""839 working: correlated Gaussian reference, organized record entropy split.

Full four-mode Fock diagnostics include pairing and cross-region couplings.
The original ten-mode entropy is evaluated analytically. Continuum limits
remain a separate proof obligation.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'gaussian_reference_entropy_split_results.json'
sys.path.insert(0,str(HERE.parent/'837'))
from organized_reference_modular_audit import entropy as code_entropy

def kron_all(items):
    out=np.array([[1.]],complex)
    for item in items:out=np.kron(out,item)
    return out
def logm(a):
    w,u=np.linalg.eigh(a);assert w.min()>0
    return (u*np.log(w))@u.conj().T
def state(k):
    w,u=np.linalg.eigh(k);ex=np.exp(-(w-w.min()))
    return (u*(ex/ex.sum()))@u.conj().T
def entropy(rho):return float(-np.trace(rho@logm(rho)).real)
def relative(rho,sigma):return float(np.trace(rho@(logm(rho)-logm(sigma))).real)

def run():
    x=np.array([[0,1],[1,0]],complex);y=np.array([[0,-1j],[1j,0]]);z=np.diag([1.,-1.]);unit=np.eye(2)
    gamma=[]
    for j in range(4):
        for a in (x,y):gamma.append(kron_all([z]*j+[a]+[unit]*(3-j)))
    rng=np.random.default_rng(839);m=rng.normal(size=(8,8));m=(m-m.T)*.13
    k=sum(1j*m[a,b]*gamma[a]@gamma[b] for a in range(8) for b in range(a+1,8))
    reference=state(k)
    nu=np.einsum('aiaj->ij',reference.reshape(8,2,8,2))
    gaussianized=np.kron(np.eye(8)/8,nu)
    parity=kron_all([z,z,z]);pc=np.kron(parity,unit)
    background_cost=relative(gaussianized,reference)
    assert np.linalg.norm(reference-gaussianized)>1e-3
    assert np.linalg.norm(reference@gaussianized-gaussianized@reference)>1e-4
    rows=[]
    for bias in (0.,.2,.6,-.7):
        sigma=(np.eye(8)+bias*parity)/8;rho=np.kron(sigma,nu)
        covariance_error=max(abs(np.trace((rho-gaussianized)@gamma[a]@gamma[b]))
                             for a in range(8) for b in range(8))
        modular_error=abs(np.trace((rho-gaussianized)@logm(reference)))
        nong=np.log(8)-entropy(sigma)
        d=relative(rho,reference)
        residual=abs(d-background_cost-nong)
        assert covariance_error<1e-13 and modular_error<1e-12 and residual<1e-12
        rows.append(dict(record_bias=bias,shared_covariance_difference=float(covariance_error),
                         common_gaussian_background_relative_entropy=background_cost,
                         nonGaussian_record_entropy_deficit=float(nong),
                         full_relative_entropy=d,additive_identity_residual=residual,
                         gaussian_modular_energy_difference=float(modular_error)))
    eta=.31;bias=.6
    nonGaussian_reference=state(k+eta*pc)
    rho=np.kron((np.eye(8)+bias*parity)/8,nu)
    difference=relative(rho,nonGaussian_reference)-relative(gaussianized,nonGaussian_reference)
    missing=difference-(np.log(8)-entropy((np.eye(8)+bias*parity)/8))
    assert abs(missing-eta*bias)<1e-12
    original=[]
    for epsilon in (.1,.4,511/512):
        for r in (0.,.6):
            original.append(dict(epsilon=epsilon,r=r,
                                 ten_mode_nonGaussian_cost_nats=float(np.log(1024)-code_entropy(epsilon,r)),
                                 added_cost_from_bias=float(code_entropy(epsilon,0)-code_entropy(epsilon,r))))
    return dict(round=839,status='working_not_formal',formal_test_groups_added=0,
                all_working_checks_passed=True,
                full_Fock_modes=4,includes_pairing_and_cross_block_quadratic_terms=True,
                gaussian_reference_product_error=float(np.linalg.norm(reference-gaussianized)),
                gaussian_reference_commutator_with_gaussianized_norm=float(np.linalg.norm(reference@gaussianized-gaussianized@reference)),
                correlated_gaussian_reference_rows=rows,
                nonGaussian_reference_counterexample=dict(
                    added_six_Majorana_modular_coefficient=eta,
                    record_bias=bias,unaccounted_modular_response=float(missing),
                    expected=eta*bias),
                original_ten_mode_analytic_costs=original,
                original_continuum_relative_entropy_limit_proven=False,
                original_regional_relative_entropy_finite_proven=False,
                original_spatial_modular_flow_computed=False,
                native_preparation_or_gravity_proven=False,
                scope='For finite CAR restrictions of a fixed faithful Gaussian reference, an organized state and its Gaussianization share quadratic modular energy. Relative entropy separates exactly into the nonGaussian record deficit and the Gaussian covariance/background cost even for correlated, paired and noncommuting references. NonGaussian references need additional modular terms. Original region inclusion, limits and finiteness remain to be proved.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
