"""724: original material-selected regions, gauge gluing and moving records.

Three diagnostics accompany the full fixed-graph analytic theorem. Finite
Peter-Weyl coefficient tests are not a new physical factorization. The radial
matrix is the old 652 diagnostic, not the continuum Hamiltonian's spectrum.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
import joint_material_boundary_transport as wall
import joint_quantum_reference_forms as refs
import joint_region_energy_gluing as glue
from joint_reference_process_transport import gibbs_with_derivative

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_material_region_records_results.json'


def cdf(x):
    return .5*(1+np.vectorize(math.erf,otypes=[float])(np.asarray(x)/np.sqrt(2)))


def bin_kernel(values,lo,hi,sigma):
    u=np.asarray(values)
    mean=(u[:,None]+u[None,:])/2
    overlap=np.exp(-(u[:,None]-u[None,:])**2/(8*sigma**2))
    upper=np.ones_like(mean) if hi==np.inf else cdf((hi-mean)/sigma)
    lower=np.zeros_like(mean) if lo==-np.inf else cdf((lo-mean)/sigma)
    return overlap*(upper-lower)


def boundary_kernel(values,b,sigma):
    k=(2*np.pi*sigma**2)**(-.25)*np.exp(-(b-values)**2/(4*sigma**2))
    return np.outer(k,k)


def wall_value(phi,lam):
    return float(phi[4]-lam*np.linalg.norm(phi[:4]))


def material_energy_check():
    _,lam,_=wall.base_data();M=refs.M
    rng=np.random.default_rng(724);points=[]
    for frac in (.1,.4,.8,.96):
        p=rng.normal(size=5);points.append(p/np.linalg.norm(p)*np.sqrt(6*M)*frac)
    points += [np.array([r,0.,0.,0.,.5]) for r in (1e-2,1e-4,1e-6)]
    rows=[];epsbound=0.
    for p in points:
        h=np.linalg.norm(p[:4]);value=wall_value(p,lam)
        grad=np.r_[-lam*p[:4]/h,1.]
        exact=float(grad@refs.original.inverse(p)@grad)
        formula=float(refs.original.F(p)*(1+lam**2-value**2/(6*M)))
        step=min(2e-6,h*1e-3)
        fd=np.array([(wall_value(p+step*e,lam)-wall_value(p-step*e,lam))/(2*step)
                     for e in np.eye(5)])
        error=float(np.linalg.norm(fd-grad)/(1+np.linalg.norm(grad)))
        assert abs(exact-formula)<2e-12 and error<2e-7
        for eps in (.01,.1,.5):
            regular=np.r_[-lam*p[:4]/np.sqrt(h*h+eps*eps),1.]
            squared=float(regular@refs.original.inverse(p)@regular)
            epsbound=max(epsbound,squared)
            assert 0<=squared<=M*(1+lam**2)+1e-10
        rows.append(dict(h=float(h),field=value,gradient_squared=exact,
                         algebra_error=abs(exact-formula),relative_gradient_fd_error=error))
    sigma=.3;x,w=np.polynomial.legendre.leggauss(220)
    z=10*sigma*x;weights=10*sigma*w
    k=(2*np.pi*sigma*sigma)**(-.25)*np.exp(-z*z/(4*sigma*sigma))
    kp=z*k/(2*sigma*sigma)
    moments=[float(weights@(k*k)),float(weights@(k*kp)),float(weights@(kp*kp))]
    assert abs(moments[0]-1)<5e-14 and abs(moments[1])<5e-14
    assert abs(moments[2]-1/(4*sigma*sigma))<2e-13
    return dict(lambda_from651=float(lam),original_inverse_metric=True,rows=rows,
                weak_gradient_bound=float(M*(1+lam**2)),regularized_gradient_max=epsbound,
                gaussian_integrals=moments,sigma=sigma,
                per_node_energy_bound_for_original_w08=float(refs.HBAR**2*M*(1+lam**2)/(8*refs.VOLUME*sigma**2)),
                continuum_energy_identity_proved_analytically=True)


def nested_region_check():
    c,lam,Y=wall.base_data();sigma=.3;b0=float(Y[1])
    field=np.array([c['S']-.01-lam*(c['H']-.015),
                    c['S']+.01-lam*(c['H']+.015)])
    configurations=np.array(list(itertools.product(range(2),repeat=2)))
    values=[field[configurations[:,v]] for v in range(2)]
    cuts=[-np.inf,b0-.35,b0+.4,np.inf]
    bins=[[bin_kernel(v,cuts[j],cuts[j+1],sigma) for j in range(3)] for v in values]
    rng=np.random.default_rng(72401)
    z=rng.normal(size=(32,32))+1j*rng.normal(size=(32,32))
    rho=z@z.conj().T;rho/=np.trace(rho)
    Jlink=glue.cut_matrix(2)
    unitary=glue.gauge.group_exp(np.array([.2,-.4,.3]),2)
    action=np.kron(np.kron(np.kron(np.eye(2),unitary.conj().T),unitary.T),np.eye(2))
    gauss_error=float(np.max(abs(action@Jlink-Jlink)))
    assert gauss_error<1e-13
    def J(label):
        return np.kron(np.kron(np.eye(4),Jlink if label[0]!=label[1] else np.eye(4)),np.eye(2))
    # A declared finite energy/source probe, not a simulation of all original CAR.
    a=rng.normal(size=(32,32))+1j*rng.normal(size=(32,32))
    H=a.conj().T@a/32+np.eye(32)
    a=rng.normal(size=(32,32))+1j*rng.normal(size=(32,32));G=(a+a.conj().T)/2
    fine={};outputs={};iso_error=0.;energy_error=0.;source_error=0.;luders_error=0.
    min_kernel=0.;min_output=0.
    for label in itertools.product(range(3),repeat=2):
        m=bins[0][label[0]]*bins[1][label[1]]
        min_kernel=min(min_kernel,float(np.linalg.eigvalsh(m).min()))
        y=np.kron(m,np.ones((8,8)))*rho;fine[label]=y
        j=J(label);out=j@y@j.conj().T;outputs[label]=out
        iso_error=max(iso_error,float(np.max(abs(j.conj().T@j-np.eye(32)))))
        min_output=min(min_output,float(np.linalg.eigvalsh(out).min()))
        energy_error=max(energy_error,float(abs(np.trace((j@H@j.conj().T)@out)-np.trace(H@y))))
        source_error=max(source_error,float(abs(np.trace((j@G@j.conj().T)@out)-np.trace(G@y))))
        # Square-root effect implements another instrument, despite same probabilities.
        wrong=np.sqrt(np.outer(np.diag(m),np.diag(m)))
        luders_error+=float(np.linalg.norm(np.kron(wrong-m,np.ones((8,8)))*rho))
    pi=(0,0,1);coarse_error=0.;glue_error=0.;probabilities=[]
    bounds=(-np.inf,cuts[2],np.inf)
    for coarse in itertools.product(range(2),repeat=2):
        j=J(coarse)
        m=bin_kernel(values[0],bounds[coarse[0]],bounds[coarse[0]+1],sigma)
        m*=bin_kernel(values[1],bounds[coarse[1]],bounds[coarse[1]+1],sigma)
        direct=np.kron(m,np.ones((8,8)))*rho
        accumulated=np.zeros_like(j@direct@j.conj().T);pulled=np.zeros_like(rho)
        for label,out in outputs.items():
            if tuple(pi[x] for x in label)==coarse:
                r=j@J(label).conj().T
                accumulated+=r@out@r.conj().T;pulled+=fine[label]
        coarse_error=max(coarse_error,float(np.max(abs(accumulated-j@direct@j.conj().T))))
        glue_error=max(glue_error,float(np.max(abs(pulled-direct))))
        probabilities.append(float(np.trace(direct).real))
    total=sum(fine.values())
    expected=rho.copy()
    for v in values:expected*=np.kron(np.exp(-(v[:,None]-v[None,:])**2/(8*sigma*sigma)),np.ones((8,8)))
    completeness_error=float(abs(sum(probabilities)-1))
    assert max(iso_error,energy_error,source_error,coarse_error,glue_error,completeness_error)<1e-12
    assert min(min_kernel,min_output)>-1e-12 and luders_error>1e-3
    assert np.max(abs(total-expected))<1e-13
    return dict(nodes=2,scalar_configurations=4,weak_link_coefficient_dimension=4,
                passive_reference_dimension=2,fine_records=9,coarse_records=4,
                coarse_probabilities=probabilities,cut_Gauss_error=gauss_error,
                isometry_error=iso_error,energy_transport_error=energy_error,
                source_transport_error=source_error,nested_record_gluing_error=coarse_error,
                pulled_back_coarse_instrument_error=glue_error,trace_error=completeness_error,
                minimum_bin_kernel_eigenvalue=min_kernel,minimum_output_eigenvalue=min_output,
                actual_vs_Luders_poststate_gap=luders_error,
                same_single_measurement_not_remeasurement=True,
                full_endpoint_Gauss_and_CAR_claim_is_analytic_not_this_fixture=True)


def moving_record_check():
    d=refs.diagnostic();_,lam,Y=wall.base_data();n=8
    hs=np.linspace(.15,1.55,n+2)[1:-1];ss=np.linspace(-1.1,1.1,n+2)[1:-1]
    h,s=np.meshgrid(hs,ss,indexing='ij');values=(s-lam*h).ravel()
    sigma=1.1;nu=7.;g0=.017;b0=float(Y[1])
    def at(g):
        H=np.exp(-6*g)*d['kinetic']+np.exp(6*g)*d['potential']
        G=-6*np.exp(-6*g)*d['kinetic']+6*np.exp(6*g)*d['potential']
        rho,rhop=gibbs_with_derivative(H,G);b=b0+nu*g
        kernels=[bin_kernel(values,-np.inf,b,sigma),bin_kernel(values,b,np.inf,sigma)]
        flux=nu*boundary_kernel(values,b,sigma)*rho
        outputs=[m*rho for m in kernels]
        primes=[m*rhop+sign*flux for m,sign in zip(kernels,(1,-1))]
        p=np.array([np.trace(y).real for y in outputs])
        ep=np.array([np.trace(G@y+H@yp).real for y,yp in zip(outputs,primes)])
        E=np.array([np.trace(H@y).real for y in outputs])
        pp=np.array([np.trace(yp).real for yp in primes])
        return dict(H=H,G=G,rho=rho,rhop=rhop,kernels=kernels,outputs=outputs,
                    primes=primes,p=p,E=E,ep=ep,pp=pp,conditional=E/p,
                    conditional_prime=ep/p-E*pp/p**2,flux=flux)
    center=at(g0);step=2e-6;plus=at(g0+step);minus=at(g0-step)
    derivative_error=max(float(np.linalg.norm((yp-ym)/(2*step)-analytic))
                         for yp,ym,analytic in zip(plus['outputs'],minus['outputs'],center['primes']))
    energy_error=float(np.max(abs((plus['E']-minus['E'])/(2*step)-center['ep'])))
    conditional_error=float(np.max(abs((plus['conditional']-minus['conditional'])/(2*step)-center['conditional_prime'])))
    omitted_flux=float(np.trace(center['H']@center['flux']).real)
    probability_flux=float(np.trace(center['flux']).real)
    total=bin_kernel(values,-np.inf,np.inf,sigma)*center['rho']
    assert np.max(abs(total-sum(center['outputs'])))<1e-14
    assert derivative_error<1e-7 and energy_error<1e-6 and conditional_error<2e-5
    assert abs(omitted_flux)>.01 and probability_flux>.01
    totalprime=sum(m*center['rhop'] for m in center['kernels'])
    cancellation_error=float(np.max(abs(totalprime-sum(center['primes']))))
    assert cancellation_error<1e-13
    return dict(original652_dimension=64,boundary_speed=nu,sigma=sigma,
                probabilities=center['p'].tolist(),energy_numerators=center['E'].tolist(),
                energy_derivatives=center['ep'].tolist(),
                conditional_energy_derivatives=center['conditional_prime'].tolist(),
                state_derivative_fd_error=derivative_error,energy_derivative_fd_error=energy_error,
                conditional_energy_derivative_fd_error=conditional_error,
                omitted_branch_energy_flux=omitted_flux,omitted_branch_probability_flux=probability_flux,
                pulled_back_boundary_flux_cancellation_error=cancellation_error,
                record_boundary_flux_not_GHY_or_Einstein_equation=True)


def run():
    results=dict(original_wall_energy=material_energy_check(),nested_gauge_regions=nested_region_check(),
                 moving_record_sources=moving_record_check())
    names=('joint_material_boundary_transport.py','joint_quantum_reference_forms.py',
           'joint_region_energy_gluing.py','joint_reference_process_transport.py',
           'research_note_563.md','research_note_617.md','research_note_647.md',
           'research_note_649.md','research_note_651.md','research_note_652.md','research_note_723.md',
           'round724_drafts/relational_wall_entry.py','round724_drafts/relational_wall_entry_results.json')
    return dict(round=724,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in names},
                scope='Same original material wall and full fixed-graph Gauss Hamiltonian: finite-resolution field records select vertex regions, original matching isometries transport all sectors and energy; coarsening the same record commutes with gluing. Moving record thresholds contribute branch source flux. No autonomous apparatus, physical graph refinement, sharp zero-cost boundaries, spacetime causal-type reconstruction or quantum gravity gluing.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
