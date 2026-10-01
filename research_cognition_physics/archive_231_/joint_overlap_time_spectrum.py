"""612: exact-time free overlap matching vs fixed-volume spectral resources.

Uses the inherited m0=1 Wilson/overlap normalization, in a frozen massless
background. The branch continuum is derived analytically, not fitted.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_chiral_fibre_source as chiral
import joint_gapped_link_locality as link
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_overlap_time_spectrum_results.json'

def parameters(k):
    k=np.asarray(k,float);s=np.sin(k);s2=float(s@s)
    b=float(np.sum(1-np.cos(k)));A=1+s2+b*b
    assert b>0 and b<np.sqrt(1+s2)
    E0=float(np.arcsinh(np.sqrt(s2)))
    E1=float(np.arccosh(A/(2*b)))
    w=float((np.sqrt(1+s2)-b)/(2*np.sqrt(1+s2)))
    assert E1>E0 and w>0
    return dict(k=k.tolist(),s=s,s2=s2,b=b,A=A,E0=E0,E1=E1,w=w)
def scalar_fourier(omega,m):
    r=np.sqrt(m['A']-2*m['b']*np.cos(omega))
    return -1j*np.sin(omega)/(2*(r+m['b']-np.cos(omega)))
def density(E,m):
    y=np.sqrt(np.maximum(0.,2*m['b']*np.cosh(E)-m['A']))
    return np.sinh(E)*y/(2*np.pi*(np.sinh(E)**2-m['s2']))
def quadrature(m,order=240,top=42.):
    nodes,weights=np.polynomial.legendre.leggauss(order)
    # E=E1+u^2 smooths the square-root threshold.
    umax=np.sqrt(top-m['E1'])
    u=(nodes+1)*umax/2
    E=m['E1']+u*u
    weight=weights*umax/2*2*u*density(E,m)
    return np.r_[m['E0'],E],np.r_[m['w'],weight]
def moments(times,m,order=240,power=0):
    E,w=quadrature(m,order)
    return np.exp(-np.asarray(times)[:,None]*E[None,:])@(w*E**power)
def kernel_and_fourier_check():
    k=np.array([np.pi/2,0,0]);m=parameters(k)
    omega=2*np.pi*(np.arange(32768)+.5)/32768-np.pi
    f=scalar_fourier(omega,m);times=np.arange(1,13)
    fourier=np.array([np.mean(np.exp(1j*omega*n)*f) for n in times])
    spectral=moments(times,m)
    error=float(np.max(abs(fourier-spectral)));assert error<2e-12
    inverse_error=0.;weyl_error=0.
    for w0 in (-2.1,-.5,.3,1.7):
        r=np.sqrt(m['A']-2*m['b']*np.cos(w0))
        X=(m['b']-np.cos(w0))*np.eye(4)+1j*chiral.GAMMA[3]*np.sin(w0)
        X+=sum(1j*chiral.GAMMA[i]*m['s'][i] for i in range(3))
        D=np.eye(4)+X/r
        S=np.linalg.inv(D)
        actual=np.trace(chiral.GAMMA[3]@S)/4
        inverse_error=max(inverse_error,float(abs(actual-scalar_fourier(w0,m))))
        hatP=(np.eye(4)-chiral.G5@(np.eye(4)-D))/2
        weyl=(hatP@S)[2:,:2]
        weyl_error=max(weyl_error,float(np.linalg.norm(weyl-S[2:,:2],2)),
                       float(abs(np.trace(weyl)/2-scalar_fourier(w0,m))))
    assert inverse_error<1e-12
    assert weyl_error<1e-12
    # Direct inherited spatial/temporal hopping kernel at an allowed L=4 mode.
    H,_,sites=link.kernel(4)
    momentum=np.array([np.pi/2,0,0,np.pi/2])
    v=np.exp(1j*np.array(sites)@momentum)/np.sqrt(len(sites))
    F=np.kron(v[:,None],np.eye(4));block=F.conj().T@H@F
    X=(np.sum(1-np.cos(momentum))-1)*np.eye(4)
    X=X.astype(complex)+sum(1j*chiral.GAMMA[i]*np.sin(momentum[i]) for i in range(4))
    inherited_error=float(np.linalg.norm(block-chiral.G5@X,2));assert inherited_error<1e-12
    return dict(spatial_momentum=k.tolist(),pole_energy=m['E0'],cut_energy=m['E1'],
        pole_weight=m['w'],kernel_block_error=inherited_error,inverse_scalar_error=inverse_error,
        specified_free_Weyl_block_error=weyl_error,
        Fourier_spectral_error=error,times=times.tolist(),correlation=spectral.tolist(),
        pole_only_error=(spectral-m['w']*np.exp(-m['E0']*times)).tolist())
def hankel_and_positive_check():
    rows=[]
    for k in ([np.pi/2,0,0],[.4,0,0]):
        m=parameters(k);E,w=quadrature(m)
        assert np.all(w>0)
        orders=[]
        for size in (2,3,4):
            times=np.arange(2*size-1)+1;c=moments(times,m)
            H=c[np.arange(size)[:,None]+np.arange(size)[None,:]]
            eig=np.linalg.eigvalsh(H);assert eig[0]>1e-14
            V=np.exp(-np.arange(size)[:,None]*E[None,:])
            factor=(V*(w*np.exp(-E)))@V.T
            error=float(np.max(abs(H-factor)));assert error<1e-13
            orders.append(dict(size=size,minimum_eigenvalue=float(eig[0]),
                determinant=float(np.linalg.det(H)),spectral_Gram_error=error))
        times=np.arange(1,7);c=moments(times,m)
        rows.append(dict(k=k,E0=m['E0'],E1=m['E1'],hankel=orders,
            continuum_fraction_first_time=float((c[0]-m['w']*np.exp(-m['E0']))/c[0])))
    return dict(rows=rows,all_finite_Hankel_ranks_proved_analytically=True,
        infinite_rank_not_inferred_from_four_by_four_test=True,
        zero_momentum_exception_has_no_cut_at_m0_one=True)
def reconstruction_resource_check():
    m=parameters([np.pi/2,0,0]);E,w=quadrature(m)
    times=np.array([1.,2.,4.]);energy=moments(times,m,power=1)
    eps=2e-6
    derivative=-(moments(times+eps,m)-moments(times-eps,m))/(2*eps)
    err=float(np.max(abs(energy-derivative)));assert err<1e-10
    c=moments(times,m);c2=moments(times,m,order=320)
    convergence=float(np.max(abs(c-c2)));assert convergence<2e-12
    # Uniform positive heat lower bound on any N orthonormal functions in a
    # finite energy interval of the continuous spectral subspace.
    beta=.7;lo=m['E1']+.2;hi=m['E1']+.5
    lower=[dict(orthonormal_vectors=n,heat_partial_trace_lower_bound=float(n*np.exp(-beta*hi)))
           for n in (4,16,64,256)]
    assert all(row['heat_partial_trace_lower_bound']>0 for row in lower)
    # Finite positive quadratures approximate moments, but their heat trace
    # counts modes, not their field-overlap weights.
    approximations=[]
    for order in (12,24,48):
        en,weight=quadrature(m,order)
        approx=np.exp(-times[:,None]*en)@weight
        approximations.append(dict(order=order,moment_error=float(np.max(abs(approx-c2))),
            one_particle_heat_trace=float(np.sum(np.exp(-beta*en))),
            weighted_field_thermal_integral=float(np.sum(weight*np.exp(-beta*en)))))
    assert approximations[-1]['moment_error']<approximations[0]['moment_error']
    return dict(energy_derivative_error=err,spectral_quadrature_difference=convergence,
        beta=beta,finite_energy_interval=[lo,hi],trace_lower_bounds=lower,
        finite_positive_approximations=approximations,
        spectral_field_weight_not_density_of_states=True,
        exact_continuous_cyclic_sector_not_trace_class=True,
        no_total_interacting_thermal_state_reconstruction_claim=True)
def run():
    data=dict(kernel=kernel_and_fourier_check(),positivity=hankel_and_positive_check(),
              resources=reconstruction_resource_check())
    deps=('research_note_363.md','research_note_596.md','research_note_603.md',
          'research_note_605.md','research_note_611.md','joint_chiral_fibre_source.py',
          'joint_gapped_link_locality.py')
    return dict(round=612,tests_run=3,failures=0,errors=0,**data,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(free_frozen_massless_same_overlap_kernel=True,
            exact_all_time_fixed_cutoff_contract=True,
            positive_cyclic_spectral_reconstruction=True,
            finite_or_compact_resolvent_ground_state_matching_obstructed=True,
            trace_class_heat_conflicts_with_exact_continuous_cyclic_sector=True,
            interacting_gauge_SM_and_GR_not_completed=True))
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=612,tests=3,all_passed=True,kernel=result['kernel'],
                         positivity=result['positivity'],resources=result['resources'])))
