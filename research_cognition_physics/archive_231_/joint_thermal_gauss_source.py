"""603: thermal trace and form-source certificates for the full finite graph.

Analytic extension uses the original confining potential, including fermions
and Gauss restriction. Numerics validate its tail bound and the response
inequality on the exact inherited neutral Fock factor; not full-graph spectra.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_quantum_response_matching as previous
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_thermal_gauss_source_results.json'

def constants():
    L,u,_=original.lattice.scalar.parameters()
    ell=float(np.linalg.det(L)/np.trace(L));D=6*original.M-u.sum()
    a=ell*D*D/(512*original.M**2)
    b=4/np.sqrt(6)
    R0=np.sqrt(6)*np.arctanh(np.sqrt((6*original.M+u.sum())/(12*original.M)))
    return float(a),float(b),float(R0)

def integrate(f,lo,hi,n):
    x,w=np.polynomial.legendre.leggauss(n)
    r=(lo+hi)/2+(hi-lo)*x/2
    return float((hi-lo)/2*np.sum(w*f(r)))

def thermal_tail_check():
    a,b,R0=constants();omega=8*np.pi**2/3
    density=lambda r:omega*(np.sqrt(6)*np.sinh(r/np.sqrt(6)))**4
    core=integrate(density,0,R0,192);rows=[]
    directions=[np.array([0.,0.,0.,0.,1.]),np.array([1.,0.,0.,0.,0.]),
                np.array([1.,0.,0.,0.,1.])/np.sqrt(2)]
    for alpha in (.01,.2,2.):
        c=alpha*a;Rmax=max(R0+.1,np.log(60/c)/b)
        integrand=lambda r:density(r)*np.exp(-c*np.exp(b*r))
        coarse=integrate(integrand,R0,Rmax,96);fine=integrate(integrand,R0,Rmax,192)
        tail=omega*9/(4*b*c)*np.exp(-c*np.exp(b*R0))
        omitted=omega*9/(4*b*c)*np.exp(-c*np.exp(b*Rmax))
        assert fine+omitted<tail and abs(fine-coarse)<1e-7*max(1,fine)
        rays=[]
        for direction in directions:
            def actual(r):
                p=np.sqrt(6*original.M)*np.tanh(r/np.sqrt(6))[:,None]*direction
                return density(r)*np.exp(-alpha*original.node_potential(p))
            value=integrate(actual,R0,Rmax,192)
            assert value<=fine*(1+1e-10)
            rays.append(value)
        rows.append(dict(alpha=alpha,core_volume=core,rigorous_tail_upper_bound=float(tail),
                         envelope_quadrature=fine,envelope_omitted_upper_bound=float(omitted),
                         quadrature_change=abs(fine-coarse),original_potential_ray_integrals=rays,
                         full_angular_integral_upper_bound=float(core+tail)))
    # Known H5 heat kernel, curvature -1/6, recurrence from H3.
    heat=[]
    for t in (.2,1.,3.):
        diag=(4*np.pi*t)**-2.5*np.exp(-2*t/3)*(1+t/9)
        errors=[]
        for r in (.04,.02,.01):
            # Scale to unit curvature: x=r/sqrt(6), tau=t/6.
            x=r/np.sqrt(6);tau=t/6
            p3=(4*np.pi*tau)**-1.5*np.exp(-tau-x*x/(4*tau))*x/np.sinh(x)
            log_derivative=1/x-1/np.tanh(x)-x/(2*tau)
            p5=-np.exp(-3*tau)/(2*np.pi*np.sinh(x))*p3*log_derivative/6**2.5
            errors.append(float(abs(p5-diag)/diag))
        assert errors[0]>12*errors[-1] and errors[-1]<2e-4
        heat.append(dict(t=t,diagonal_kernel=float(diag),radial_relative_errors=errors))
    return dict(a=a,b=b,R0=R0,rows=rows,heat_kernel_scaling=heat,
                no_full_graph_partition_function_evaluated=True,
                ray_integrals_are_diagnostics_not_angular_proof=True)

def response_bound_check():
    _,pairs=previous.matrices(previous.S0,[24,25,30,31])
    H,G,_=[previous.fock(pair) for pair in pairs]
    e,v,p,rho=previous.thermal(H,2.)
    aa=e-e.min()+1;g=v.conj().T@G@v
    normalized=g/np.sqrt(aa[:,None]*aa[None,:])
    L=float(np.linalg.norm(normalized,2))
    de=e[None,:]-e[:,None];same=abs(de)<previous.TOL
    weight=np.zeros_like(de);np.divide(p[:,None]-p[None,:],de,out=weight,where=~same)
    mean=float(np.sum(p*np.diag(g)).real);centered=g-mean*np.eye(len(e))
    D=float(np.sum(p[:,None]*same*abs(centered)**2))
    chi0=float(np.sum(weight*abs(g)**2));chie=chi0+2*D
    moment1=float(p@aa);moment2=float(p@(aa*aa))
    Dbound=L*L*moment2
    bound=L*L*(4*moment1+5*2*moment2)
    assert 0<D<Dbound and 0<chi0<chie<bound
    rows=[]
    for eta in (.2,.1,.05,.025):
        kernel=np.zeros_like(de);np.divide(de*de,de*de+eta*eta,out=kernel)
        value=float(np.sum(weight*kernel*abs(g)**2))
        assert 0<value<chi0
        rows.append(dict(Abel_parameter=eta,retarded_susceptibility_imaginary_axis=value,
                         difference_from_zero=chi0-value))
    assert all(rows[j+1]['difference_from_zero']<rows[j]['difference_from_zero'] for j in range(3))
    return dict(exact_inherited_Fock_dimension=16,beta=2.,form_constant=L,
                first_shifted_energy_moment=moment1,second_shifted_energy_moment=moment2,
                conserved_variance=D,variance_bound=Dbound,
                static_susceptibility=chie,isolated_zero_susceptibility=chi0,
                universal_form_susceptibility_bound=bound,Abel_rows=rows,
                full_graph_extension_is_analytic_not_a_matrix_substitution=True)

def run():
    a=thermal_tail_check();b=response_bound_check()
    deps=('research_note_579.md','research_note_589.md','research_note_590.md','research_note_591.md',
          'research_note_598.md','research_note_602.md','joint_curved_quantum_source.py',
          'joint_quantum_response_matching.py')
    return dict(round=603,tests_run=2,failures=0,errors=0,thermal_tail=a,response_bound=b,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(full_nonlinear_fixed_graph_Gibbs_existence_proved=True,
                           finite_CAR_and_Gauss_restriction_retained=True,
                           source_zero_frequency_form_bounds_proved=True,
                           instantaneous_noise_not_claimed=True,
                           equilibrium_preparation_and_temperature_are_inputs=True,
                           no_cutoff_uniform_or_quantum_GR_or_continuum_completion=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=603,tests=2,all_passed=True,form_bound=result['response_bound']['universal_form_susceptibility_bound'])))
