"""601: one consistent EFT order for matter, lapse and constraint evolution.

General consistency is analytic; numerics use the original pure-singlet
classical sector and the entire known 581 scalar-loop reduced coefficient.
Not a full quantum constraint calculation or continuum-limit proof.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_eft_constraint_order_results.json'
A4=-5/18  # Declared Lorentz continuation of the 581 scalar coefficient.
MAT,VAC,_=original.lattice.scalar.parameters()


def tensor_potential(phi):
    """Holomorphic formulas; no complex conjugation in analytic derivatives."""
    f=original.F(phi)
    d=np.array([np.sum(phi[:4]**2)-VAC[0],phi[4]**2-VAC[1]])
    ld=MAT@d;num=d@MAT@d
    group=np.array([0,0,0,0,1])
    n1=4*phi*ld[group]
    n2=4*np.diag(ld[group])+8*np.outer(phi,phi)*MAT[group[:,None],group[None,:]]
    f1=-phi/3;f2=-np.eye(5)/3
    U=num/(4*f*f)
    grad=n1/(4*f*f)-num*f1/(2*f**3)
    hess=(n2/(4*f*f)-(np.outer(n1,f1)+np.outer(f1,n1))/(2*f**3)
          +3*num*np.outer(f1,f1)/(2*f**4)-num*f2/(2*f**3))
    hcov=hess-(np.outer(grad,phi)+np.outer(phi,grad))/(6*f)
    inv=original.inverse(phi);mixed=inv@hcov
    V=U*U-2*U*np.trace(mixed)/3+np.trace(mixed@mixed)/2+grad@inv@grad/6
    return U,grad,hess,V


def coeff_raw(chi):
    phi=np.zeros(5,dtype=np.result_type(chi,float))
    phi[4]=np.sqrt(6*original.M)*np.tanh(chi/np.sqrt(6))
    U,grad,hess,V=tensor_potential(phi)
    return np.array([U,U/9,V])


def coefficients(chi):
    val=coeff_raw(chi).real
    deriv=coeff_raw(chi+1e-24j).imag/1e-24
    return tuple(float(z) for z in (*val,*deriv))


def original_coefficient_check():
    rows=[]
    for chi in (.4,.8,1.1):
        phi=np.zeros(5);phi[4]=np.sqrt(6*original.M)*np.tanh(chi/np.sqrt(6))
        U,g,H,V=tensor_potential(phi)
        estimates=[]
        for eps in (1e-4,5e-5):
            numerical=np.zeros((5,5))
            for i in range(5):
                for j in range(5):
                    ei=np.eye(5)[i]*eps;ej=np.eye(5)[j]*eps
                    numerical[i,j]=(original.node_potential(phi+ei+ej)-original.node_potential(phi+ei-ej)
                                    -original.node_potential(phi-ei+ej)+original.node_potential(phi-ei-ej))/(4*eps**2)
            estimates.append(numerical)
        numerical=(4*estimates[1]-estimates[0])/3
        assert abs(U-original.node_potential(phi))<1e-14
        he=float(np.max(abs(H-numerical)));assert he<2e-8
        *_,up,bp,vp=coefficients(chi)
        derivative_errors=[]
        for first_step in (1e-4,5e-5,2.5e-5):
            fd=(coeff_raw(chi+first_step)-coeff_raw(chi-first_step))/(2*first_step)
            derivative_errors.append(float(np.max(abs(fd-np.array([up,bp,vp])))))
        derivative_error=derivative_errors[-1]
        assert derivative_error<2e-8
        assert derivative_errors[0]>12*derivative_error
        rows.append(dict(chi=chi,U=float(U),V1=float(V),hessian_error=he,
                         hessian_raw_errors=[float(np.max(abs(H-e))) for e in estimates],
                         coefficient_derivative_error=derivative_error,
                         derivative_step_halving_errors=derivative_errors))
    return dict(rows=rows,full_five_field_hessian_retained=True,
                pure_singlet_background_not_one_field_loop_replacement=True)


def derivative(x):
    k=np.fft.fftfreq(len(x),1/len(x))
    return np.fft.ifft(1j*k*np.fft.fft(x)).real


def hamiltonian_check():
    x=np.arange(512)*2*np.pi/512
    chi=.7+.08*np.sin(x);sp=.08*np.cos(x);p=.18+.025*np.cos(2*x);y=sp*sp
    N=1+.2*np.cos(x);M=1+.15*np.sin(2*x)
    co=np.array([coefficients(c) for c in chi]);U,b,V,up,bp,vp=co.T
    vv=N*derivative(M)-M*derivative(N)
    mean=lambda f:float(2*np.pi*np.mean(f))
    D=mean(vv*p*sp);aa=A4*(p*p-y)+b
    rows=[]
    for lam in (.1,.05,.025):
        hp=p*(1-lam*aa);hy=(1+lam*aa)/2
        hc=up-lam*(bp*(p*p-y)/2-vp)
        dn=N*hc-derivative(N*2*hy*sp)
        dm=M*hc-derivative(M*2*hy*sp)
        bracket=mean(dn*M*hp-N*hp*dm)
        expected=-lam*lam*mean(vv*p*aa*aa*sp)
        residual=bracket-D
        assert abs(residual-expected)<3e-16
        rows.append(dict(eft_parameter=lam,bracket=bracket,momentum_constraint=D,
                         algebra_defect=residual,expected_order_two_defect=expected))
    assert abs(rows[0]['algebra_defect'])>1e-10
    # Cubic inverse has a regular small branch and additional cutoff-scale roots.
    U0,b0,V0,*_=coefficients(.7);pp=.2;yy=.01;lam=.05
    raw_roots=np.roots([lam*A4,0.,1+lam*(b0-A4*yy),-pp])
    assert np.max(abs(raw_roots.imag))<1e-12
    roots=np.sort(raw_roots.real)
    branches=[]
    for u in roots:
        z=1+lam*(b0+A4*(3*u*u-yy))
        branches.append(dict(velocity=float(u),legendre_derivative=float(z),
                             quartic_expansion_size=float(abs(lam*A4*u*u))))
    low=min(branches,key=lambda row:abs(row['velocity']-pp))
    assert low['legendre_derivative']>.99
    assert all(r['quartic_expansion_size']>.9 for r in branches if r is not low)
    return dict(rows=rows,cubic_branches=branches,continuum_functional_bracket_only=True,
                exact_graph_or_full_fermionic_HDA_not_claimed=True)


def local(chi,u,H,lam,scheme):
    U,b,V,up,bp,vp=coefficients(chi);X=u*u/2
    px=1+lam*(2*A4*X+b);z=1+lam*(6*A4*X+b)
    rho=X+U+lam*(3*A4*X*X+b*X+V)
    f0=-3*H*u-up
    e1=(6*A4*X+b)*f0+3*H*(2*A4*X+b)*u+bp*X+vp
    if scheme=='mixed':acc=f0
    elif scheme=='reduced':acc=f0-lam*e1
    elif scheme=='resummed':acc=-(3*H*px*u+up+lam*(bp*X+vp))/z
    else:raise ValueError(scheme)
    hdot=-X*px
    cdot=6*H*hdot-(up+lam*(bp*X+vp))*u-u*z*acc
    expected={'mixed':-lam*u*e1,'reduced':lam*lam*u*(6*A4*X+b)*e1,'resummed':0}[scheme]
    assert abs(cdot-expected)<1e-14
    return np.array([u,acc,hdot]),rho,z,cdot


def evolve(lam,scheme,steps=160):
    chi=.8;u=.2;U,b,V,*_=coefficients(chi);X=u*u/2
    rho=X+U+lam*(3*A4*X*X+b*X+V);assert rho>0
    state=np.array([chi,u,np.sqrt(rho/3)]);dt=.5/steps
    zmin=1e9
    for _ in range(steps):
        def f(q):
            nonlocal zmin
            rhs,_,z,_=local(*q,lam,scheme);zmin=min(zmin,z)
            return rhs
        k1=f(state);k2=f(state+dt*k1/2);k3=f(state+dt*k2/2);k4=f(state+dt*k3)
        state+=dt*(k1+2*k2+2*k3+k4)/6
    _,rho,_,_=local(*state,lam,scheme)
    assert zmin>.99
    return dict(final_state=state.tolist(),final_constraint=float(3*state[2]**2-rho),
                minimum_legendre_derivative=float(zmin),steps=steps)


def constraint_evolution_check():
    rows=[]
    for lam in (.04,.02,.01,.005,.0025):
        row=dict(eft_parameter=lam)
        for scheme in ('mixed','reduced','resummed'):row[scheme]=evolve(lam,scheme)
        assert abs(row['resummed']['final_constraint'])<2e-13
        rows.append(row)
    ratios={}
    for scheme,target in (('mixed',2),('reduced',4)):
        ratios[scheme]=[abs(rows[j][scheme]['final_constraint']/rows[j+1][scheme]['final_constraint'])
                        for j in range(len(rows)-1)]
        # Finite-lambda trajectories also change; test convergence to the
        # analytically derived order, retaining the non-asymptotic samples.
        deviations=[abs(r-target) for r in ratios[scheme]]
        assert all(deviations[j+1]<deviations[j] for j in range(len(deviations)-1))
        assert all(d<.1 for d in deviations[-2:])
    fine=evolve(.04,'reduced',320)
    refinement=abs(fine['final_constraint']-rows[0]['reduced']['final_constraint'])
    assert refinement<1e-12
    return dict(rows=rows,parameter_halving_ratios=ratios,step_refinement_error=refinement,
                time_interval=.5,initial_constraint_satisfied=True,
                resummed_first_gradient_action_is_comparison_not_all_orders_physics=True)


def run():
    coefficients_check=original_coefficient_check()
    hamiltonian=hamiltonian_check();evolution=constraint_evolution_check()
    names=('research_note_344.md','research_note_351.md','research_note_352.md','research_note_581.md',
           'research_note_585.md','research_note_588.md','research_note_600.md',
           'joint_curved_quantum_source.py')
    return dict(round=601,tests_run=3,failures=0,errors=0,coefficients=coefficients_check,
                hamiltonian=hamiltonian,evolution=evolution,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope=dict(original_known_scalar_EFT_sector=True,Lorentz_continuation_declared=True,
                           one_common_perturbative_order=True,full_five_field_scalar_loop_potential_used=True,
                           finite_window_classical_constraint_checks=True,
                           general_Noether_compatibility_not_global_existence=True,
                           no_full_quantum_HDA_or_unified_GR_completion_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=601,tests=3,all_passed=True,
                         halving_ratios=result['evolution']['parameter_halving_ratios'])))
