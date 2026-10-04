"""722: original-energy failure of a phase readout and a resolvent alternative.

Domain theorems are analytic, not inferred from matrix cutoffs. Numerical
groups check the original 5D measure/state lift, actual rational instruments,
and geometry derivatives on the inherited radial diagnostic. The nonzero
spatial term there is a conditioned-neighbour reference, not a full graph.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_velocity_record_energy as previous
from round722_drafts import squared_reference_entry as entry

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_squared_reference_readout_results.json'
old=previous.old
M=old.M


def radius(xi,z,w=1.):
    lo=np.full_like(np.asarray(xi+z,dtype=float),-15.)
    hi=lo+20.
    A=1+z*z/6
    for _ in range(90):
        q=(lo+hi)/2
        value=w/M*(A*q+np.exp(2*q)/12)
        lo=np.where(value<xi,q,lo);hi=np.where(value>=xi,q,hi)
    return np.exp((lo+hi)/2)


def gl(n,a,b):
    t,w=np.polynomial.legendre.leggauss(n)
    return (a+b)/2+(b-a)*t/2,w*(b-a)/2


def lift_quadrature(n,radial=False):
    # Normalized constant S3 harmonic and the original empty CAR vacuum.
    # b(z)=(1-4z²)^2 has compact H1 support, sufficient for the energy form.
    z,wz=gl(n,-.5,.5);bnorm=np.sqrt(np.dot(wz,(1-4*z*z)**4))
    Z=z[:,None];B=(1-4*Z*Z)**2/bnorm
    Bz=-16*Z*(1-4*Z*Z)/bnorm
    norm=kin=pot=hardy=0.
    for a,b in ((-1.,0.),(0.,1.)):
        t,wt=gl(n,0,1)
        if radial:
            r0=radius(a,Z);r1=radius(b,Z)
            R=r0+(r1-r0)*t
            xi=((1+Z*Z/6)*np.log(R)+R*R/12)/M
            dxi=((1+Z*Z/6)/R+R/6)/M
            weights=wz[:,None]*wt*(r1-r0)*dxi
        else:
            xi=a+(b-a)*t+0*Z;R=radius(xi,Z)
            weights=wz[:,None]*wt*(b-a)
        D=1+(R*R+Z*Z)/6
        U=np.sqrt(1.5)*(1-abs(xi))*B
        Uxi=-np.sqrt(1.5)*np.sign(xi)*B
        Uz=np.sqrt(1.5)*(1-abs(xi))*Bz
        xiR=((1+Z*Z/6)/R+R/6)/M
        xiz=Z*np.log(R)/(3*M)
        # sqrt(m) * derivatives of the physical five-dimensional wavefunction.
        AR=Uxi*xiR-.5*U*(4/R-R/(2*D))
        Az=Uz+Uxi*xiz+.25*U*Z/D
        kval=AR*AR*(1+R*R/6)+Az*Az*(1+Z*Z/6)+AR*Az*R*Z/3
        phi=np.zeros(R.shape+(5,));phi[...,0]=R*np.sqrt(M/D)
        phi[...,4]=Z*np.sqrt(M/D)
        norm+=np.sum(weights*U*U)
        kin+=old.HBAR**2/2*np.sum(weights*kval)
        pot+=np.sum(weights*U*U*old.original.node_potential(phi))
        hardy+=np.sum(weights*U*U/R**2)
    return dict(norm=float(norm),original_kinetic=float(kin),
                original_potential=float(pot),weighted_Hardy_norm=float(hardy))


def measure_and_failure_check():
    errors=[]
    for x in (np.array([.2,.3,-.4,.5,.1]),np.array([1.,2.,-.1,.7,-.4])):
        D=1+x@x/6;F=M/D;p=x*np.sqrt(F)
        J=np.sqrt(F)*(np.eye(5)-np.outer(x,x)/(6*D))
        pull=J.T@old.original.metric(p)@J
        gx=np.eye(5)-np.outer(x,x)/(6+x@x)
        errors.append(float(np.max(abs(pull-gx))))
        volume=np.sqrt(M)/F**3*np.linalg.det(J)
        errors.append(abs(volume-D**-.5))
    assert max(errors)<3e-15
    # Classical four-dimensional radial Hardy identity, not the H5 spectral one.
    r,w=gl(48,0,1);v=(1-r)**3;dv=-3*(1-r)**2
    lhs=np.dot(w,r**3*dv**2-r*v**2)
    rhs=np.dot(w,r**3*(dv+v/r)**2)
    assert abs(lhs-rhs)<1e-14 and lhs>0
    lift=lift_quadrature(56);other=lift_quadrature(56,True)
    discrepancy=max(abs(lift[k]-other[k]) for k in lift)
    assert discrepancy<2e-10 and abs(lift['norm']-1)<1e-13
    # Reuse (not re-count) the already saved scalar-tail certificate.
    scalar=entry.run()
    assert scalar==json.loads(entry.TARGET.read_text('utf8'))
    # Analytic interval lower bound: integral over [-L,-L+1].
    cmin=M/(1+.5**2/6)
    amp2constant=6/np.pi*np.sin(.25)**2
    tails=[]
    for L in (80.,160.,320.):
        assert (L-1)/2>=scalar['lower_threshold']
        tails.append(dict(L=L,log_integral_lower_bound=float(
            np.log(amp2constant)+2*cmin*(L-1)-4*np.log(L))))
    assert all(x['log_integral_lower_bound']>200 for x in tails)
    return dict(original_measure_metric_error=max(errors),Hardy_identity_error=abs(lhs-rhs),
                finite_energy_Gauss_lift=lift,radial_vs_flow_quadrature_error=discrepancy,
                analytic_tail_interval_bounds=tails,
                infinity_proved_analytically_not_by_cutoff=True,
                all32_CAR_modes_retained_empty_vacuum_used=True)


def fixture(g,label,lam=16.):
    d=old.diagnostic();_,basis=np.linalg.eigh(d['kinetic']+d['potential'])
    T=np.exp(-6*g)*basis.T@d['kinetic']@basis
    P=np.exp(6*g)*basis.T@d['potential']@basis
    H=T+P;G=-6*T+6*P
    index=0 if label=='T' else 1;V=np.exp(-6*g)*d['V'][index]
    I=np.eye(len(H));f=d['refs'][index]
    Wbase=(f-.2*I)@(f-.2*I)
    # One conditioned-neighbour difference, unit grid spacing, same f.
    W=np.exp(-4*g)*Wbase
    b=2*np.linalg.norm(Wbase,2)+1
    R=V@V+b*I-W;Rprime=-12*V@V+4*W
    assert np.linalg.eigvalsh(R).min()>0
    Z=np.linalg.inv(lam*I+1j*R)
    K0=lam*Z;Ks=[K0,I-K0]
    K0prime=-1j*lam*Z@Rprime@Z;Kprimes=[K0prime,-K0prime]
    eh,uh=np.linalg.eigh(H)
    p=np.exp(-old.old625.BETA*(eh-eh[0]));p/=sum(p)
    rho=(uh*p)@uh.conj().T
    B=sum(k.conj().T@H@k for k in Ks)-H
    energy=float(np.trace(rho@B).real)
    source=sum(k.conj().T@G@k for k in Ks)-G
    instrument=sum(kp.conj().T@H@k+k.conj().T@H@kp for k,kp in zip(Ks,Kprimes))
    return locals()


def readout_check():
    rows=[]
    for label in ('T','s'):
        d=fixture(.017,label)
        I=d['I'];R=d['R'];Ks=d['Ks'];lam=d['lam'];rho=d['rho']
        E0=lam**2*np.linalg.inv(lam**2*I+R@R)
        complete=sum(k.conj().T@k for k in Ks)
        errors=[np.linalg.norm(complete-I,2),np.linalg.norm(Ks[0].conj().T@Ks[0]-E0,2)]
        probe=R@np.linalg.inv(I+R)
        errors.append(np.linalg.norm(sum(k.conj().T@probe@k for k in Ks)-probe,2))
        assert max(errors)<2e-12
        comm=np.linalg.norm(d['W']@d['V']-d['V']@d['W'],2)
        assert comm>.01
        Rwrong=d['V']@d['V']+d['b']*I
        wrong=lam**2*np.linalg.inv(lam**2*I+Rwrong@Rwrong)
        diff=abs(np.trace(rho@(wrong-E0)).real);assert diff>1e-4
        A=d['H']+(1-min(d['eh']))*I;root=previous.root(A);inv=np.linalg.inv(root)
        norms=[float(np.linalg.norm(root@k@inv,2)) for k in Ks]
        actual=sum(float(np.trace(A@k@rho@k.conj().T).real) for k in Ks)
        bound=sum(v*v for v in norms)*float(np.trace(A@rho).real)
        assert actual<=bound+1e-10
        rows.append(dict(reference=label,instrument_spectral_errors=[float(x) for x in errors],
                         spatial_velocity_commutator=float(comm),
                         dropping_spatial_term_probability_error=float(diff),
                         record0_probability=float(np.trace(rho@E0).real),
                         actual_shifted_energy=actual,finite_matrix_energy_bound=bound))
    # Numerical identity supporting the analytic Neumann step, with arbitrary
    # noncommuting bounded D. This is not a numerical infinite-domain proof.
    d=fixture(.017,'T');lam=256.;C=2.;omega=.7;hbar=.7
    eta=C/(np.sqrt(lam)*(np.sqrt(lam/2)-hbar*omega))
    D=d['b']*d['I']-d['W'];V=d['V'];I=d['I']
    dnorm=np.linalg.norm(D,2);assert dnorm*eta<.5
    Z0=np.linalg.inv(V@V-1j*lam*I)
    Z=Z0@np.linalg.inv(I+D@Z0)
    error=float(np.linalg.norm(Z-np.linalg.inv(V@V+D-1j*lam*I),2))
    assert error<1e-14
    return dict(rows=rows,resolvent_identity_error=error,
                actual_phases_not_Luders=True,conditioned_neighbour64D_diagnostic=True,
                full_graph_energy_bound_is_analytic=True)


def source_and_history_check():
    rows=[]
    for label in ('T','s'):
        g=.017;d=fixture(g,label);step=1e-5
        plus=fixture(g+step,label);minus=fixture(g-step,label)
        kerror=max(float(np.linalg.norm((p-m)/(2*step)-kp,2))
                   for p,m,kp in zip(plus['Ks'],minus['Ks'],d['Kprimes']))
        rho_prime=(plus['rho']-minus['rho'])/(2*step)
        parts=[float(np.trace(d['rho']@d[k]).real) for k in ('source','instrument')]
        parts.append(float(np.trace(rho_prime@d['B']).real))
        total=(plus['energy']-minus['energy'])/(2*step)
        error=abs(total-sum(parts))
        assert kerror<1e-7 and error<1e-6 and abs(parts[1])>.001
        wait=(d['uh']*np.exp(-.31j*d['eh']/old.HBAR))@d['uh'].conj().T
        leaves=[]
        for k in d['Ks']:
            for j in d['Ks']:
                word=j@wait@k;leaves.append(word@d['rho']@word.conj().T)
        probs=[float(np.trace(r).real) for r in leaves]
        assert abs(sum(probs)-1)<2e-12 and min(probs)>0
        A=d['H']+(1-min(d['eh']))*d['I']
        inv=np.linalg.inv(previous.root(A))
        average=sum(k.conj().T@A@k for k in d['Ks'])
        C=float(np.linalg.eigvalsh(inv@average@inv).max())
        actual=float(np.trace(A@sum(leaves)).real)
        budget=C*C*float(np.trace(A@d['rho']).real)
        assert actual<=budget+1e-10
        rows.append(dict(reference=label,energy=d['energy'],total_derivative=total,
                         source_instrument_preparation=parts,derivative_error=error,
                         Kraus_derivative_error=kerror,full_two_record_probabilities=probs,
                         actual_final_shifted_energy=actual,finite_matrix_budget=budget))
    return dict(rows=rows,original_waiting_not_reset=True,
                no_autonomous_apparatus_or_continuum_claim=True)


def run():
    results=dict(measure_and_domain_failure=measure_and_failure_check(),
                 complete_readout=readout_check(),geometry_and_history=source_and_history_check())
    names=('joint_velocity_record_energy.py','research_note_721.md',
           'joint_quantum_reference_forms.py','research_note_652.md',
           'joint_operator_domain_completion.py','research_note_704.md',
           'round722_drafts/squared_reference_entry.py','round722_drafts/squared_reference_entry_results.json')
    return dict(round=722,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original fixed finite graph, external positive geometry, no dynamical theta. One-node phase-of-square readout can leave the original full energy form. A sufficiently broad rational instrument for full Q including spatial W preserves that form and has bounded energy-space first geometry derivative. Not an autonomous apparatus, coordinate reconstruction, uniform continuum or gravity-generation proof.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
