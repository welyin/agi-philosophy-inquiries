"""885: full original64 static fermion weight, no dynamic Gauss claim.
The original ordering is retained: logZ contains the vacuum energy.
Quadratures calibrate explicit constants; asymptotic proofs are in the note.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'884'))
import cut_source_variance_probe as cut
old=cut.old
TARGET=HERE/'joint_static_thermal_weight_results.json'

def radial_counts(N):
    k=np.arange(-(N//2),N//2+1,dtype=int)
    r,c=np.unique((k[:,None]**2+k[None,:]**2).ravel(),return_counts=True)
    return k,r.astype(float),c.astype(float)

def x_squares(m,q,N,alpha):
    vals=[]
    for k in range(-(N//2),N//2+1):
        eig=np.linalg.eigvalsh(old.H(m,q,N,[k,0,0],alpha))
        vals.append(np.sort(eig*eig))
    return np.array(vals)

def log_trace_difference(a,b,rad,count):
    # a,b are all original64 x-Hamiltonian squared eigenvalues.
    vacuum=0.;thermal=0.
    for aa,bb in zip(a,b):
        ea=np.sqrt(aa[:,None]+rad[None,:])
        eb=np.sqrt(bb[:,None]+rad[None,:])
        de=np.divide((aa-bb)[:,None],ea+eb,
            out=np.zeros_like(ea),where=(ea+eb)>0)
        vacuum+=old.BETA/4*float(np.sum(de*count))
        thermal+=.5*float(np.sum((np.log1p(np.exp(-old.BETA*ea))-
                np.log1p(np.exp(-old.BETA*eb)))*count))
    return dict(logZ_difference=vacuum+thermal,vacuum_contribution=vacuum,
                excitation_contribution=thermal)

def original_matrix_checks(m,q):
    M=np.exp(old.XI)*m.MASS
    errors={}
    errors['mass_clifford_anticommutation']=max(float(np.max(abs(g@M+M@g))) for g in m.GAMMA)
    errors['charge_clifford_commutator']=max(float(np.max(abs(g@np.diag(q)-np.diag(q)@g))) for g in m.GAMMA)
    errors['source_mass_cross_trace']=float(abs(np.trace(m.GAMMA[0]@np.diag(q)@M)))
    square=0.;direct=0.
    for N,alpha in ((5,.02),(17,-.08)):
        hx=old.H(m,q,N,[N//2,0,0],alpha)
        h=old.H(m,q,N,[N//2,1,-2],alpha)
        square=max(square,float(np.max(abs(h@h-hx@hx-5*np.eye(64)))))
        ex=np.linalg.eigvalsh(hx)
        e=np.linalg.eigvalsh(h)
        predicted=.5*np.sum(np.logaddexp(old.BETA*np.sqrt(ex*ex+5)/2,
                                         -old.BETA*np.sqrt(ex*ex+5)/2))
        measured=.5*np.sum(np.logaddexp(old.BETA*e/2,-old.BETA*e/2))
        direct=max(direct,float(abs(predicted-measured)))
    errors['perpendicular_square_error']=square
    errors['direct64_logcosh_error']=direct
    assert max(errors.values())<2e-11
    return errors

def constants(q,order=64):
    x,w=np.polynomial.legendre.leggauss(order)
    y=x/2;weights=w[:,None]*w[None,:]/4
    r2=y[:,None]**2+y[None,:]**2
    # Patch constant has a smooth face-integral representation.
    patch=old.BETA*np.sum(q*q)/8*float(np.sum(weights/np.sqrt(.25+r2)))
    # Split quadrant and integrate in polar coordinates to avoid a cusp at0.
    theta=(x+1)*np.pi/8;wt=w*np.pi/8
    radius=.5/np.cos(theta)
    # Integral r [sqrt(1/4+r^2)-r] dr, from0 to radial boundary.
    primitive=((.25+radius*radius)**1.5-.125-radius**3)/3
    cutconst=old.BETA*8*float(np.dot(wt,primitive)) # multiplicity2 times beta/2
    return dict(patch_N_squared_coefficient=patch,cut_N_cubed_loss_coefficient=cutconst)

def run():
    m=old.load();q,_=old.charges(m)
    checks=original_matrix_checks(m,q)
    c32=constants(q,32);c64=constants(q,64)
    assert max(abs(c32[k]-c64[k]) for k in c32)<1e-10
    rows=[]
    for N in (17,33,65,129):
        k,rad,count=radial_counts(N)
        base=x_squares(m,q,N,0.)
        zero_mass=np.repeat((k*k)[:,None],64,axis=1)
        baseline_mass=log_trace_difference(base,zero_mass,rad,count)['logZ_difference']
        patch=[]
        for alpha in (.02,.04):
            a=x_squares(m,q,N,alpha)
            d=log_trace_difference(a,base,rad,count)
            patch.append(dict(alpha=alpha,**d,
                scaled_coefficient=d['logZ_difference']/(alpha*alpha*N*N)))
        u,wp=cut.solve_u(N,1.)
        alpha=-(.5-u)/6
        a=x_squares(m,q,N,alpha)
        d=log_trace_difference(a,base,rad,count)
        # For uniform prior on[-.1,.1], Z_total >= |J| Z_massless(0), J=[-1/30,1/30].
        density_log_bound=np.log(15.)+d['logZ_difference']+baseline_mass
        rows.append(dict(N=N,original_components=64,total_momentum_points=N**3,
            baseline_mass_logZ_shift=baseline_mass,patch=patch,
            cut_alpha=alpha,cut_w=1.,cut_slope=wp,cut_difference=d,
            cut_logZ_difference_over_N_cubed=d['logZ_difference']/N**3,
            cut_log_density_upper_bound=density_log_bound))
    assert rows[0]['cut_log_density_upper_bound']>0  # coarse lower bound is inconclusive at17
    assert all(row['cut_log_density_upper_bound']<0 for row in rows[1:])
    assert all(row['patch'][0]['logZ_difference']>0 for row in rows)
    return dict(round=885,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3670,
        argument_scope='Static classical flat-parameter integration with the unmodified original64 fermion Fock weight, fixed bounded prior positive near0, and883 smooth-cut hopping. Its884 cut band is exponentially suppressed, while an O(N^2) vacuum bias remains already inside the exact continuum patch. This is not the dynamic Gauss Gibbs marginal or a completed Q/E bridge.',
        beta=old.BETA,xi=old.XI,eta=old.ETA,
        original_charge_square_trace=float(np.sum(q*q)),matrix_checks=checks,
        constants=c64,quadrature32_vs64={k:abs(c32[k]-c64[k]) for k in c32},
        rows=rows,
        original_Majorana_and_all_species_retained=True,
        original_Fock_ordering_and_vacuum_energy_retained=True,
        static_cut_suppression_proved=True,
        fixed_patch_quadratic_vacuum_bias_proved=True,
        full_dynamic_Gauss_measure_derived=False,
        full_averaged_source_convergence_proved=False,
        nonlocal_counterterm_added=False,
        full_interacting_Q_E_bridge_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
