"""721: original velocity-reference instruments preserve full fixed-graph energy.

The full energy-form proof is analytic. Numerical groups check target flows,
all original mass coefficients, and the inherited64D radial diagnostic.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quantum_reference_forms as old
import joint_fermion_gauss_completion as matter

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_velocity_record_energy_results.json'
original=old.original;M=old.M
L,VAC,_=original.lattice.scalar.parameters()


def data(p,label,w=1.):
    F=original.F(p);K=original.inverse(p)
    df=np.zeros(5);ddf=np.zeros((5,5))
    if label=='T':
        df[:4]=p[:4];ddf[:4,:4]=np.eye(4)
        div=F*(4-p[:4]@p[:4]/(2*M))
        ddiv=-p/3*(4-p[:4]@p[:4]/(2*M))-F/M*np.r_[p[:4],0.]
        ell=2*M;beta=7*np.sqrt(6)*M/3;c=M/2
    else:
        df[4]=1.;div=-F*p[4]/(3*M)
        ddiv=-(F*df-p[4]*p/3)/(3*M)
        ell=np.sqrt(6*M)/3;beta=np.sqrt(M);c=2*np.sqrt(2*M)/9
    X=K@df
    DY=np.zeros((5,5))
    for j in range(5):
        e=np.eye(5)[j]
        dk=-p[j]/3*(np.eye(5)-np.outer(p,p)/(6*M))
        dk-=F*(np.outer(e,p)+np.outer(p,e))/(6*M)
        DY[:,j]=(dk@df+K@ddf[:,j])/w
    hessian=ddf-(np.outer(p,df)+np.outer(df,p))/(6*F)
    return X/w,DY,div/w,hessian,ddiv,ell/w,beta/w,c/w


def flow(p,label,t,w=1.,steps=128):
    y=np.r_[p,np.eye(5).ravel(),0.]
    def rhs(v):
        X,D,div,*_=data(v[:5],label,w)
        return np.r_[X,(D@v[5:30].reshape(5,5)).ravel(),div]
    h=t/steps
    for _ in range(steps):
        a=rhs(y);b=rhs(y+h*a/2);c=rhs(y+h*b/2);d=rhs(y+h*c)
        y+=h*(a+2*b+2*c+d)/6
    return y[:5],y[5:30].reshape(5,5),y[-1]


def norm(a):
    return float(np.linalg.norm(a,2))


def root(a):
    e,u=np.linalg.eigh(a)
    return (u*np.sqrt(e))@u.T


def target_check():
    rng=np.random.default_rng(721)
    rows=[];calculus=[]
    for label in ('s','T'):
        for frac in (.2,.7,.94):
            p=rng.normal(size=5);p*=np.sqrt(6*M)*frac/np.linalg.norm(p)
            Y,D,div,hess,ddiv,ell,beta,c=data(p,label,1.2)
            K=original.inverse(p)
            assert norm(root(K)@hess@root(K))<=ell*1.2+1e-12
            assert np.sqrt(ddiv@K@ddiv)<=beta*1.2+1e-12
            eps=2e-6
            fd=np.column_stack([(data(p+eps*e,label,1.2)[0]-data(p-eps*e,label,1.2)[0])/(2*eps) for e in np.eye(5)])
            calculus.append(norm(fd-D))
            t=.37
            q,J,logj=flow(p,label,t,1.2);q2,_,_=flow(p,label,t,1.2,256)
            gram=root(K)@J.T@original.metric(q)@J@root(K)
            stretch=float(np.linalg.eigvalsh(gram).max())
            actualj=np.linalg.slogdet(J)[1]-3*np.log(original.F(q)/original.F(p))
            inv=flow(q,label,-t,1.2)[0]
            assert stretch<=np.exp(2*ell*abs(t))*(1+1e-11)
            assert abs(np.log(original.F(q)/original.F(p)))<=c*abs(t)+1e-12
            assert abs(logj-actualj)<2e-10 and norm(inv-p)<2e-10
            # Half-density Jacobian gradient is not discarded.
            grad=np.array([(flow(p+eps*e,label,t,1.2)[2]-flow(p-eps*e,label,t,1.2)[2])/(2*eps) for e in np.eye(5)])
            gradnorm=float(np.sqrt(grad@K@grad))
            assert gradnorm<=beta*abs(t)*np.exp(ell*abs(t))+1e-9
            rows.append(dict(reference=label,radius_fraction=frac,
                             flow_refinement_error=norm(q2-q),inverse_error=norm(inv-p),
                             stretch_squared=stretch,stretch_bound=float(np.exp(2*ell*abs(t))),
                             measure_jacobian_error=abs(logj-actualj),
                             log_jacobian_gradient=gradnorm,
                             log_jacobian_gradient_bound=float(beta*abs(t)*np.exp(ell*abs(t)))))
    assert max(calculus)<2e-9
    return dict(rows=rows,max_cartesian_variational_error=max(calculus),
                analytic_global_bounds_not_inferred_from_samples=True)


def coefficient_check():
    lam=np.linalg.det(L)/np.trace(L);D=6*M-sum(VAC)
    A=32/(lam*D*D);B=144/(D*D)
    numerator=norm(L)/4*sum((6*M+VAC)**2)
    rng=np.random.default_rng(7211);rows=[];coordinate_errors=[]
    for label in ('s','T'):
        ratios=[];diffs=[]
        for frac in (.2,.7,.98):
            p=rng.normal(size=5);p*=np.sqrt(6*M)*frac/np.linalg.norm(p)
            t=.21;c=data(p,label)[-1]
            F=original.F(p);x=p/np.sqrt(F);Y=data(p,label)[0]
            dx=Y/np.sqrt(F)+p*(p@Y)/(6*F**1.5)
            desired=np.r_[F*x[:4],0.] if label=='T' else np.r_[np.zeros(4),np.sqrt(F)]
            coordinate_errors.append(norm(dx-desired))
            pp=flow(p,label,t)[0];pm=flow(p,label,-t)[0]
            pot=float(original.node_potential(p))
            bound=numerator*np.exp(2*c*t)*(A*pot+B)
            ratios.extend([float(original.node_potential(q))/bound for q in (pp,pm)])
            h,d=matter.mass_matrices(p)
            hp,dp=matter.mass_matrices(pp);hm,dm=matter.mass_matrices(pm)
            diffs.append(dict(radius_fraction=frac,Dirac_feedback=norm((hp+hm)/2-h),
                              Majorana_feedback=norm((dp+dm)/2-d)))
            assert max(ratios)<=1+1e-10
            # Nontrivial original weak/U1 transformation; covariance retained.
            W=matter.gauge.group_exp(np.array([.3,-.2,.4]),2);z=np.exp(.23j)
            X=z**3*W@(p[:2]+1j*p[2:4])
            gp=np.r_[X.real,X.imag,p[4]]
            flowed=flow(gp,label,t)[0]
            Xq=z**3*W@(pp[:2]+1j*pp[2:4])
            assert norm(flowed-np.r_[Xq.real,Xq.imag,pp[4]])<1e-11
        rows.append(dict(reference=label,max_potential_bound_ratio=max(ratios),mass_rows=diffs))
    assert max(r['Dirac_feedback'] for r in rows[0]['mass_rows'])<1e-11
    assert max(r['Majorana_feedback'] for r in rows[1]['mass_rows'])<1e-11
    assert max(r['Majorana_feedback'] for r in rows[0]['mass_rows'])>1e-5
    assert max(r['Dirac_feedback'] for r in rows[1]['mass_rows'])>1e-3
    assert max(coordinate_errors)<1e-13
    return dict(rows=rows,confinement_A=A,confinement_B=B,numerator_bound=numerator,
                complementary_mass_flow_identity_error=max(coordinate_errors),
                all_original32_mode_masses_retained=True)


def record_fixture(g,label,t=.13):
    d=old.diagnostic()
    _,basis=np.linalg.eigh(d['kinetic']+d['potential'])
    T=basis.T@d['kinetic']@basis;W=basis.T@d['potential']@basis
    H=np.exp(-6*g)*T+np.exp(6*g)*W
    G=-6*np.exp(-6*g)*T+6*np.exp(6*g)*W
    V=np.exp(-6*g)*d['V'][0 if label=='T' else 1]
    e,u=np.linalg.eigh(H);probs=np.exp(-old.old625.BETA*(e-e[0]));probs/=sum(probs)
    rho=(u*probs)@u.conj().T
    ev,uv=np.linalg.eigh(V);U=(uv*np.exp(-.5j*t*ev/old.HBAR))@uv.conj().T
    Ud=3j*t/old.HBAR*V@U
    Ks=[(U+1j*U.conj().T)/2,(U-1j*U.conj().T)/2]
    Kds=[(Ud+1j*Ud.conj().T)/2,(Ud-1j*Ud.conj().T)/2]
    B=sum(k.conj().T@H@k for k in Ks)-H
    energy=float(np.trace(rho@B).real)
    source=sum(k.conj().T@G@k for k in Ks)-G
    instrument=sum(kd.conj().T@H@k+k.conj().T@H@kd for k,kd in zip(Ks,Kds))
    return locals()


def actual_history_check():
    rows=[]
    for label in ('s','T'):
        d=record_fixture(.017,label);h=1e-5
        plus=record_fixture(.017+h,label);minus=record_fixture(.017-h,label)
        drho=(plus['rho']-minus['rho'])/(2*h)
        fd=(plus['energy']-minus['energy'])/(2*h)
        pieces=[float(np.trace(d['rho']@d[k]).real) for k in ('source','instrument')]
        pieces.append(float(np.trace(drho@d['B']).real))
        err=abs(fd-sum(pieces))
        assert err<2e-6 and abs(pieces[1])>.01
        # Full records and poststates with same original waiting Hamiltonian.
        eh,uh=np.linalg.eigh(d['H']);wait=(uh*np.exp(-.31j*eh/old.HBAR))@uh.conj().T
        leaves=[]
        for k1 in d['Ks']:
            for k2 in d['Ks']:
                a=k2@wait@k1;leaves.append(a@d['rho']@a.conj().T)
        total=sum(leaves);probs=[float(np.trace(x).real) for x in leaves]
        assert abs(sum(probs)-1)<1e-12
        A=d['H']+(1-min(eh))*np.eye(len(eh))
        inv=np.linalg.inv(root(A))
        average=sum(k.conj().T@A@k for k in d['Ks'])
        C=float(np.linalg.eigvalsh(inv@average@inv).max())
        budget=C*C*np.trace(A@d['rho']).real
        actual=np.trace(A@total).real
        assert actual<=budget+1e-10
        rows.append(dict(reference=label,energy=d['energy'],geometry_total_fd=fd,
                         source_instrument_preparation_terms=pieces,derivative_error=err,
                         two_record_probabilities=probs,actual_final_shifted_energy=float(actual),
                         finite_matrix_two_step_budget=float(budget)))
    return dict(rows=rows,conditional64D_radial_diagnostic=True,
                no_full_Gauss_spectrum_or_continuum_simulation=True)


def run():
    results=dict(target_flow=target_check(),original_coefficients=coefficient_check(),
                 records_and_sources=actual_history_check())
    names=('joint_quantum_reference_forms.py','joint_quantum_reference_forms_results.json',
           'joint_fermion_gauss_completion.py','joint_operator_domain_completion.py',
           'round721_drafts/material_velocity_entry_results.json','research_note_652.md','research_note_720.md')
    return dict(round=721,tests_run=3,failures=0,errors=0,results=results,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        scope='Full original fixed finite graph, fixed positive external geometry: velocity-flow instruments preserve the energy form; actual records and fixed-instrument source forms retained. Optional dynamical theta, uniform spatial refinement, autonomous apparatus and arbitrary finite-energy total source differentiability not proved.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
