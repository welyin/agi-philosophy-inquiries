"""Working887: c-number thermal repair versus the same EM dynamic source.
All original matter is retained in the matrix check. No dynamic gauge claim.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'886'))
import unbroken_holonomy_locality as last
weight=last.weight;cut=weight.cut;old=last.old
sys.path.insert(0,str(HERE.parent/'884'))
import cut_time_smeared_lower_bound as previous
TARGET=HERE/'renormalized_weight_variance_probe_results.json'

def integral(N,m2,order=48):
    x,a=np.polynomial.legendre.leggauss(order);ws=.75*x+1.25;a=.75*a
    total=0.
    for w,aw in zip(ws,a):
        _,wp=cut.solve_u(N,w)
        e=np.sqrt(w*w+m2);r=np.exp(-old.BETA*e);v=r/(1+r)**2
        total+=aw*24*v*w*w/(e*e)*wp
    return float(total)

def run():
    m=old.load();_,matter=old.charges(m);q=last.em_charges(m,matter)
    mass=np.exp(old.XI)*m.MASS[26,28];m2=float(abs(mass)**2)
    audit=[]
    for N in (17,65,257):
        for w in (.5,1.,2.):
            u,wp=cut.solve_u(N,w);alpha=-(.5-u)/6
            e=np.sqrt(w*w+m2);r=np.exp(-old.BETA*e);v=r/(1+r)**2
            selected=144*v*w*w/(e*e)*wp*wp
            actual=0.
            for k in (N//2,-N//2):
                h=old.H(m,q,N,[k,0,0],alpha)
                g=m.GAMMA[0]@np.diag(q*cut.source(N,k+alpha*q))
                actual+=previous.diagonal_lower(h,g)/2
            relative=abs(actual-selected)/max(1.,selected)
            # Other original modes can add; only the selected lower bound is asserted.
            assert actual>=selected*(1-2e-10)
            audit.append(dict(N=N,w=w,full64_diagonal_lower=actual,
                exact_electron_lower=selected,relative_difference=relative))
    rows=[]
    for N in (17,65,257,4097,1048577):
        val=integral(N,m2)
        rows.append(dict(N=N,unweighted_lower_integral=val,scaled_by_logN_squared=val/np.log(N)**2))
    x,a=np.polynomial.legendre.leggauss(64);w=.75*x+1.25;a=.75*a
    e=np.sqrt(w*w+m2);r=np.exp(-old.BETA*e);v=r/(1+r)**2
    coefficient=float(24/old.ETA*np.dot(a,w**3/(w*w+m2)*v))
    qerror=abs(integral(257,m2,32)-integral(257,m2,64))
    assert qerror<1e-8 and coefficient>0
    # Actual finite Hamiltonian c-number shift: prepared Gibbs unchanged exactly.
    h=old.H(m,q,17,[0,0,0],-.04)
    en,vec=np.linalg.eigh(h);p=np.exp(-old.BETA*(en-en.min()));p/=p.sum()
    shift=31.2
    pshift=np.exp(-old.BETA*((en+shift)-(en+shift).min()));pshift/=pshift.sum()
    assert np.max(abs(p-pshift))<2e-14
    # This last Gibbs identity is matrix algebra calibration; the Fock identity is analytic.
    return dict(working_round=887,formal_round_still=886,cumulative_groups_still=3671,
        prior_goal_turn_classification='progress',
        original_mass_squared=m2,matrix_lower_bound_checks=audit,rows=rows,
        coefficient=coefficient,quadrature32_vs64_error=qerror,
        c_number_shift_preparation_calibration=float(np.max(abs(p-pshift))),
        exact_phase_cancellation_is_analytic=True,
        matching_weights_requires_uniform_positive_density_on_cut_band=True,
        weak_measure_convergence_alone_is_insufficient=True,
        full_dynamic_gauge_counterterm_refuted=False,
        status='Working same-background tradeoff; formal scope and source-contract review pending.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
