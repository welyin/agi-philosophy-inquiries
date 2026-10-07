"""879: full self-dual covariance, graded regional product, and parity audit.
A four-CAR finite matrix calibration with pairing and hopping.
This is not a discretization of the original space-dependent field theory.
"""
from pathlib import Path
import argparse,json
import numpy as np
TARGET=Path(__file__).with_name('graded_regional_covariance_results.json')
def adj(a):return a.conj().T
def herm(a):return (a+adj(a))/2
def fun(a,f):
    v,u=np.linalg.eigh(herm(a))
    return (u*f(v))@adj(u)
def norm(a):return float(np.linalg.norm(a,2))
def sqrt(a):
    assert np.linalg.eigvalsh(herm(a)).min()>-1e-12
    return fun(a,lambda x:np.sqrt(np.maximum(0,x)))
def annihilator(n,mode):
    z=np.diag([1,-1]);a=np.array([[0,1],[0,0]])
    out=np.array([[1.]])
    for i in range(n):out=np.kron(out,z if i<mode else a if i==mode else np.eye(2))
    return out.astype(complex)
def selfdual_cov(rho,ops):
    return np.array([[np.trace(rho@a@adj(b)) for b in ops] for a in ops])
def run():
    n=4;beta=0.8
    c=[annihilator(n,k) for k in range(n)]
    number=[adj(x)@x for x in c];identity=np.eye(2**n)
    h=np.array([[.7,.11+.08j,.23,0],[.11-.08j,-.2,0,.17j],
                [.23,0,.45,-.09],[0,-.17j,-.09,.9]],dtype=complex)
    pairing=np.array([[0,.12,.04,.15],[-.12,0,-.08,.02],
                      [-.04,.08,0,.1],[-.15,-.02,-.1,0]],dtype=complex)
    b=np.block([[h,pairing],[-pairing.conj(),-h.T]])
    ph=np.block([[np.zeros((n,n)),np.eye(n)],[np.eye(n),np.zeros((n,n))]])
    p=fun(b,lambda x:1/(1+np.exp(-beta*x)))
    ops=c+[adj(x) for x in c]
    ham=sum(adj(c[i])@c[j]*h[i,j] for i in range(n) for j in range(n))
    create=sum(adj(c[i])@adj(c[j])*pairing[i,j]/2 for i in range(n) for j in range(n))
    ham=herm(ham+create+adj(create))
    rho=fun(ham,lambda x:np.exp(-beta*x));rho/=np.trace(rho)
    actual=selfdual_cov(rho,ops)
    slots=[np.array([0,1,4,5]),np.array([2,3,6,7])]
    product_p=np.zeros_like(p)
    for indices in slots:product_p[np.ix_(indices,indices)]=p[np.ix_(indices,indices)]
    four=rho.reshape(4,4,4,4)
    rho_a=np.einsum('abcb->ac',four)
    rho_b=np.einsum('abad->bd',four)
    rho_product=np.kron(rho_a,rho_b)
    product_actual=selfdual_cov(rho_product,ops)
    pa=(identity-2*number[0])@(identity-2*number[1])
    pb=(identity-2*number[2])@(identity-2*number[3])
    pt=pa@pb
    cross=adj(c[0])@c[2]+adj(c[2])@c[0]
    even_each=(cross+pa@cross@pa+pb@cross@pb+pt@cross@pt)/4
    pair=c[0]@c[3]
    delta=p-product_p
    eigdelta=np.linalg.eigvalsh(herm(delta))
    trace_norm=float(np.abs(eigdelta).sum())
    hs2=float(np.linalg.norm(sqrt(p)-sqrt(product_p),'fro')**2)
    mean=lambda op:float(np.trace(rho@op).real)
    cov_n=mean(number[0]@number[2])-mean(number[0])*mean(number[2])
    result=dict(round=879,date='2026-10-06',fresh_numbered_groups=1,
        cumulative_numbered_groups=3664,fixture=dict(complex_CAR_modes=n,selfdual_dimension=2*n,Fock_dimension=2**n,beta=beta),
        all_checks_passed=True,
        selfdual_checks=dict(B_charge_conjugation=norm(ph@b.conj()@ph+b),
            covariance_reality=norm(ph@p.conj()@ph-(np.eye(2*n)-p)),
            product_reality=norm(ph@product_p.conj()@ph-(np.eye(2*n)-product_p)),
            actual_Fock_covariance_defect=norm(actual-p),
            product_Fock_covariance_defect=norm(product_actual-product_p),
            product_covariance_min=float(np.linalg.eigvalsh(product_p).min()),
            product_covariance_max=float(np.linalg.eigvalsh(product_p).max())),
        regional_checks=dict(
            normal_state_trace_norm_difference=float(np.abs(np.linalg.eigvalsh(herm(rho-rho_product))).sum()),
            local_marginal_a_preserved=norm(np.einsum('abcb->ac',rho_product.reshape(4,4,4,4))-rho_a),
            local_marginal_b_preserved=norm(np.einsum('abad->bd',rho_product.reshape(4,4,4,4))-rho_b),
            original_cross_number_connected=cov_n,
            original_cross_pairing_abs=float(abs(np.trace(rho@pair))),
            original_cross_hopping_mean=mean(cross),
            product_cross_hopping_mean=float(np.trace(rho_product@cross).real),
            total_even_operator_commutator=norm(pt@cross-cross@pt),
            local_even_projection_norm=norm(even_each),
            lost_cross_even_operator_norm=norm(cross-even_each)),
        quasi_equivalence_calibration=dict(sqrt_difference_HS_squared=hs2,covariance_difference_trace_norm=trace_norm),
        analytic_scope='The full original finite-matrix Dirac/Majorana free CAR on a fixed background has a graded normal regional product for finitely many separated regions. Together with878, the original free bosonic algebra and regional even CAR admit a common normal tensor presentation. Actual correlated and non-Gaussian normal preparations are transported, not reset. Global even cross-odd products are not identified with the join of regional even algebras.',
        continuum_quasi_equivalence_proved_by_finite_fixture=False,
        global_type_I_collar_or_energy_control_proved=False,full_goal_completed=False)
    checks=result['selfdual_checks']
    assert max(v for k,v in checks.items() if k not in ('product_covariance_min','product_covariance_max'))<2e-12
    assert 0<checks['product_covariance_min']<=checks['product_covariance_max']<1
    rr=result['regional_checks']
    assert max(rr[k] for k in ('local_marginal_a_preserved','local_marginal_b_preserved','total_even_operator_commutator','local_even_projection_norm'))<2e-12
    assert abs(rr['product_cross_hopping_mean'])<1e-12
    assert abs(rr['original_cross_hopping_mean'])>0.01 and rr['original_cross_pairing_abs']>0.01
    assert abs(cov_n)>1e-4 and abs(rr['lost_cross_even_operator_norm']-1)<1e-12
    assert hs2<=trace_norm+1e-12
    return result
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
