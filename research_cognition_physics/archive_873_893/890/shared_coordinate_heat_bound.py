"""890: one coordinate kinetic budget controls the full cut Fock thermal sum.
All original charges/masses in the analytic bound. The direct finite joint
Gibbs experiment is explicitly a two-mode calibration, not the full Gauss model.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'889'))
import dynamic_cut_spectrum as prior
old=prior.old
TARGET=HERE/'shared_coordinate_heat_bound_results.json'
KAPPA=prior.KAPPA;ETA=prior.ETA;D=prior.HALF_INTERVAL
BETA=old.BETA;SHIFT=6*(abs(prior.CENTER)+D);CUT_FACTOR=.25

def entropy_bound(N,m2):
    k,r,c=prior.prior.last.weight.radial_counts(N)
    out=0.
    for x in k[1:-1]:
        e=np.sqrt(max(abs(x)-SHIFT,0)**2+m2[:,None]+r)
        out+=.5*float(np.sum(c*np.log1p(np.exp(-BETA*e/2))))
    e=np.sqrt(m2[:,None]+r)
    out+=float(np.sum(c*np.log1p(np.exp(-BETA*e/2))))
    return out

def boundary_energy(N,q,k1,alpha,m2,perp2=0.):
    b=old.buffered(N,2*np.pi/N*(k1+alpha*q))
    return np.sqrt(b*b+perp2+m2)

def original_weight_checks(q,m2):
    worst=0.;worstcut=0.
    mmin=float(m2.min())
    for N in (9,17,33):
        L=N//2
        for x in np.linspace(-D,D,41):
            alpha=prior.CENTER+x
            for k in range(-L,L+1):
                e=boundary_energy(N,q,k,alpha,m2)
                w=np.sqrt((0 if abs(k)==L else max(abs(k)-SHIFT,0))**2+m2)
                worst=min(worst,float(np.min(e-w)))
                if abs(k)==L:
                    v=CUT_FACTOR*float(prior.energy(N,x,mmin))
                    worstcut=min(worstcut,float(np.min(e-v)))
    assert worst>-2e-12 and worstcut>-2e-12
    return dict(minimum_weight_slack=worst,minimum_common_cut_slack=worstcut,
                q_values=sorted(set(q.tolist())),original_full_components=len(q))

def scalar_cut_lower(N,m2):
    # Dirichlet Sobolev/sublevel-set bound, not a finite-difference estimate.
    # For normalized psi with energy lambda: 1<=4a sqrt(lambda/kappa)+lambda/V0.
    ss=np.log(float(N))
    candidates=np.r_[np.geomspace(D/10000,D*.999,256),
                     ETA/(6*ss)*np.linspace(1.05,8,160)]
    candidates=candidates[(candidates>0)&(candidates<D)]
    u=6*candidates
    b=(float(N)/2-6*D)*old.smooth_step(u/ETA)
    V0=CUT_FACTOR*np.hypot(b,np.sqrt(m2))
    A=4*candidates/np.sqrt(KAPPA)
    lb=(2/(A+np.sqrt(A*A+4/V0)))**2
    j=int(np.argmax(lb))
    baseline=KAPPA*(np.pi/(2*D))**2+CUT_FACTOR*np.sqrt(m2)
    return dict(N=str(N),certified_formula_lower_evaluated=float(max(lb[j],baseline)),
                optimized_halfwidth=float(candidates[j]),exterior_potential_lower=float(V0[j]))

def heat_coordinate(beta):
    n=np.arange(1,2001,dtype=float)
    e=KAPPA*(np.pi*n/(2*D))**2
    z=np.exp(-beta*e)
    return float(np.log(z.sum()))

def exact_two_mode_gibbs(m,q,m2):
    # The same coordinate T appears once in the 4-sector Fock Hamiltonian.
    N=17;n=31;h=2*D/(n+1)
    xs=np.linspace(-D+h,D-h,n)
    T=np.diag(np.full(n,2*KAPPA/h**2))+np.diag(np.full(n-1,-KAPPA/h**2),1)+np.diag(np.full(n-1,-KAPPA/h**2),-1)
    e1=prior.energy(N,xs,m2)
    e2=np.sqrt((-6*(prior.CENTER+xs))**2+1+m2)
    potentials=np.array([np.zeros(n),e1,e2,e1+e2])
    # Explicit coordinate x Fock matrix with ordinary 2-mode fermion occupations.
    H=np.kron(T,np.eye(4))+np.diag(potentials.T.reshape(-1))
    e,V=np.linalg.eigh(H);w=np.exp(-BETA*(e-e.min()));z=w.sum()
    rho=(V*w)@V.T/z
    individual=[]
    for potential in potentials:
        ev=np.linalg.eigvalsh(T+np.diag(potential))
        individual.append(float(np.exp(-BETA*(ev-e.min())).sum()))
    assert abs(sum(individual)-z)<1e-9
    joint_cut=(individual[1]+individual[3])/z
    center=n//2
    diagonal=np.diag(rho).reshape(n,4)[center]
    conditional_cut=float((diagonal[1]+diagonal[3])/diagonal.sum())
    frozen_cut=float(1/(1+np.exp(BETA*e1[center])))
    assert abs(conditional_cut-frozen_cut)>.1
    assert joint_cut<frozen_cut
    return dict(N=N,coordinate_interior_points=n,total_matrix_dimension=4*n,
        full_vs_sector_partition_error=float(abs(sum(individual)-z)),
        joint_Gibbs_cut_probability=float(joint_cut),
        actual_joint_conditional_cut_probability_at_center=conditional_cut,
        frozen_pointwise_free_cut_probability_at_center=frozen_cut,
        exact_finite_matrix_min_eigenvalue=float(e.min()),
        experiment_scope='Two original dispersion modes calibrate shared-coordinate Gibbs; full species controlled by analytic bounds, not diagonalized here.')

def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.last.em_charges(m,matter)
    qq,m2=prior.prior.last.joint_spectrum(m,q)
    audit=original_weight_checks(qq,m2)
    entropy=[dict(N=N,full_original_Fock_log_product_bound=entropy_bound(N,m2)) for N in (9,17,33,65)]
    assert 0<entropy[-1]['full_original_Fock_log_product_bound']-entropy[-2]['full_original_Fock_log_product_bound']<.01
    # Closed analytic uniform ceiling; uses <=64 labels, conservatively.
    t=BETA/(2*np.sqrt(3))
    A=1+2*np.exp(t*SHIFT)/np.expm1(t)
    C3=A/np.tanh(t/2)**2
    C2=1/np.tanh(BETA/(4*np.sqrt(2)))**2
    uniform=64*C3+128*C2
    assert entropy[-1]['full_original_Fock_log_product_bound']<uniform
    logs=[]
    coordinate_ratio=heat_coordinate(BETA/2)-heat_coordinate(BETA)
    for power in (4,12,30,120,300,600):
        N=2**power+1
        lb=scalar_cut_lower(N,float(m2.min()))
        lb['log_probability_upper_ceiling']=float(min(0.,-BETA*lb['certified_formula_lower_evaluated']/2+uniform+coordinate_ratio))
        logs.append(lb)
    assert logs[-1]['log_probability_upper_ceiling']<-100
    electron_mass=float(abs(np.exp(old.XI)*m.MASS[26,28])**2)
    direct=exact_two_mode_gibbs(m,q,electron_mass)
    return dict(round=890,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3675,
        argument_scope='In889 explicitly transported and pointwise vacuum-subtracted free branch, one shared coordinate kinetic term suppresses every occupation sector containing a boundary momentum mode. A uniform full-species Fock entropy product and heat trace bound give vanishing total boundary thermal probability, without assigning one coordinate kinetic term per mode. No local Gauss or common physical-source bridge is proved.',
        beta=BETA,kappa=KAPPA,coordinate_half_interval=D,maximum_charge_shift=SHIFT,
        full_original_weight_checks=audit,entropy_product_rows=entropy,
        explicit_analytic_uniform_log_product_ceiling=float(uniform),
        scalar_cut_analytic_lower_evaluations=logs,two_mode_actual_joint_Gibbs=direct,
        full_Fock_thermal_boundary_probability_tends_to_zero=True,
        single_shared_coordinate_kinetic_preserved=True,
        own_joint_Gibbs_not_pointwise_free=True,
        full_original_mass_and_species_in_entropy_bound=True,
        local_Gauss_constraint_solved=False,original_bare_kinetic_recovered=False,
        original_vacuum_source_matched=False,full_continuous_thermal_state_convergence=False,
        source_derivatives_and_real_time_records_matched=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
