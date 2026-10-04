"""725: the same nonsmooth material wall has an energy-controlled square form.

The continuous Hardy and Gaussian tests use the original curved target
measure. The 64-dimensional matrix test retains652's original H and uses a
declared differential-form discretization for C* C, not a claimed full
Gauss spectrum or a self-adjoint signed wall velocity.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_material_boundary_transport as wall
import joint_quantum_reference_forms as old
import joint_material_region_records as records
from joint_reference_process_transport import gibbs_with_derivative

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_relational_normal_form_results.json'
M=old.M;HBAR=old.HBAR;W=old.VOLUME


def radial_data(order,alpha=.3,beta=.4,phase=.15):
    # s=R*t, h=sqrt(R^2-s^2)*r; integrate the full radial h^3 F^-3 measure.
    x,q=np.polynomial.legendre.leggauss(order)
    r=(x+1)/2;wr=q/2;t=x;wt=q
    rr,tt=np.meshgrid(r,t,indexing='ij')
    s=np.sqrt(6*M)*tt;cap=6*M-s*s;h=np.sqrt(cap)*rr
    F=M-(h*h+s*s)/6
    weight=np.sqrt(M)*cap**2*rr**3*np.sqrt(6*M)/F**3*wr[:,None]*wt[None,:]
    psi=F**3*np.exp(-alpha*h*h-beta*s*s+1j*phase*h*h*s)
    norm=float(np.sum(weight*abs(psi)**2));psi/=np.sqrt(norm)
    dh=(-h/F-2*alpha*h+2j*phase*h*s)*psi
    ds=(-s/F-2*beta*s+1j*phase*h*h)*psi
    khh=F*(1-h*h/(6*M));kss=F*(1-s*s/(6*M));khs=-F*h*s/(6*M)
    _,lam,_=wall.base_data();a=-lam;f=s+a*h
    xh=a*khh+khs;xs=a*khs+kss
    lap=F*(-f/(3*M)+3*a/h)
    Cpsi=-1j*HBAR/W*(xh*dh+xs*ds+lap*psi/2)
    gamma=F*(1+a*a-f*f/(6*M))
    return locals()


def integrate(d,value):
    return float(np.sum(d['weight']*value).real)


def quadratic(d,u,v):
    return d['khh']*abs(u)**2+d['kss']*abs(v)**2+2*d['khs']*(u.conj()*v).real


def hardy_energy_check():
    rows=[]
    for alpha,beta,phase in ((.3,.4,.15),(.05,.1,.3),(.8,.15,-.2)):
        values=[]
        for order in (60,90):
            d=radial_data(order,alpha,beta,phase)
            grad=integrate(d,quadratic(d,d['dh'],d['ds']))
            hardy=integrate(d,d['F']/d['h']**2*abs(d['psi'])**2)
            remainder=integrate(d,quadratic(d,d['dh']+d['psi']/d['h'],d['ds']))
            c2=integrate(d,abs(d['Cpsi'])**2)
            qkin=HBAR**2/(2*W)*grad
            bound=2*M*(2+11*d['a']**2)/W*qkin+2*M*HBAR**2*(1+d['a']**2)/(3*W**2)
            error=abs(grad-hardy-remainder)
            assert error<2e-10 and c2<=bound and hardy>0
            values.append([grad,hardy,remainder,c2,bound])
        convergence=float(np.max(abs(np.array(values[1])-values[0])/(1+abs(np.array(values[1])))))
        assert convergence<1e-11
        rows.append(dict(parameters=[alpha,beta,phase],gradient_form=values[-1][0],
                         original_Hardy_weight=values[-1][1],positive_remainder=values[-1][2],
                         squared_wall_derivative=values[-1][3],full_energy_sufficient_bound=values[-1][4],
                         normalized_quadrature_convergence=convergence))
    return dict(rows=rows,angular_and_CAR_extension_analytic=True,
                source_family='Normalized F^3 exp(-alpha h^2-beta s^2+i phase h^2 s), original measure.',
                signed_velocity_selfadjointness_not_used=True)


def normal_record_check():
    d=radial_data(90);sigma=.7
    x,q=np.polynomial.legendre.leggauss(140);z=10*x;qw=10*q
    k=(2*np.pi*sigma*sigma)**(-.25)*np.exp(-z*z/4)
    normal0=integrate(d,abs(d['Cpsi'])**2)
    measured=0.;prob=0.;cross=0.
    # At fixed physical y, C K_y psi = K_y C psi - i*hbar/w Gamma K_y' psi.
    # Then substitute y=f+sigma*z for this independent Gaussian quadrature.
    for t,weight,amplitude in zip(z,qw,k):
        out=amplitude*d['Cpsi']-1j*HBAR/W*d['gamma']*t/(2*sigma)*amplitude*d['psi']
        measured+=weight*sigma*integrate(d,abs(out)**2)
        prob+=weight*sigma*amplitude**2
        cross+=weight*sigma*t*amplitude**2
    mean_gamma2=integrate(d,d['gamma']**2*abs(d['psi'])**2)
    predicted=HBAR**2/(4*W*W*sigma*sigma)*mean_gamma2
    mf=integrate(d,d['f']*abs(d['psi'])**2)
    mf2=integrate(d,d['f']**2*abs(d['psi'])**2)
    spatial=2*(mf2-mf*mf)  # two identical independent nodes, epsilon=1, gamma^11=1
    error=abs((measured-normal0)-predicted)
    assert error/(1+predicted)<5e-13 and abs(prob-1)<2e-14 and abs(cross)<2e-14
    return dict(sigma=sigma,initial_normal_square=normal0,postrecord_normal_square=measured,
                exact_positive_backaction=predicted,independent_integral_error=error,
                same_two_node_spatial_difference_mean=spatial,
                initial_Q_mean=spatial-normal0,postrecord_Q_mean=spatial-measured,
                probability_normalization=prob,gaussian_cross_integral=cross,
                same_original_two_node_product_Gauss_state=True,
                sign_of_mean_not_a_classical_causal_signature=True)


def matrix_fixture():
    d=old.diagnostic();n=8
    hs=np.linspace(.15,1.55,n+2)[1:-1];ss=np.linspace(-1.1,1.1,n+2)[1:-1]
    h,s=np.meshgrid(hs,ss,indexing='ij');h=h.ravel();s=s.ravel()
    F=M-(h*h+s*s)/6;mu=np.sqrt(M)*h**3/F**3;root=np.sqrt(mu)
    def der(step):return (np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/(2*step)
    Dh=np.kron(der(hs[1]-hs[0]),np.eye(n))
    Ds=np.kron(np.eye(n),der(ss[1]-ss[0]))
    khh=F*(1-h*h/(6*M));kss=F*(1-s*s/(6*M));khs=-F*h*s/(6*M)
    def first(xh,xs,lap):
        return -1j*HBAR/W*root[:,None]*(xh[:,None]*Dh+xs[:,None]*Ds+np.diag(lap/2))/root[None,:]
    Ch=first(khh,khs,F*(3/h-h/(3*M)))
    Cs=first(khs,kss,-F*s/(3*M))
    _,lam,Y=wall.base_data();C=Cs-lam*Ch;N=C.conj().T@C
    f=s-lam*h;b0=float(Y[1]);spatial=np.diag((f-b0)**2)
    return dict(H0=d,Ch=Ch,Cs=Cs,C=C,N=N,f=f,b0=b0,spatial=spatial)


def joint_source_check():
    d=matrix_fixture();sigma=1.1;nu=7.;g0=.017;step=2e-6
    def at(g):
        H=np.exp(-6*g)*d['H0']['kinetic']+np.exp(6*g)*d['H0']['potential']
        G=-6*np.exp(-6*g)*d['H0']['kinetic']+6*np.exp(6*g)*d['H0']['potential']
        rho,rhop=gibbs_with_derivative(H,G)
        N=np.exp(-12*g)*d['N'];spatial=np.exp(-4*g)*d['spatial']
        Q=spatial-N;Qp=-4*spatial+12*N
        b=d['b0']+nu*g
        kernels=[records.bin_kernel(d['f'],-np.inf,b,sigma),
                 records.bin_kernel(d['f'],b,np.inf,sigma)]
        flux=nu*records.boundary_kernel(d['f'],b,sigma)*rho
        outputs=[m*rho for m in kernels]
        primes=[m*rhop+sign*flux for m,sign in zip(kernels,(1,-1))]
        moments=np.array([np.trace(Q@y).real for y in outputs])
        derivatives=np.array([np.trace(Qp@y+Q@yp).real for y,yp in zip(outputs,primes)])
        probs=np.array([np.trace(y).real for y in outputs])
        return locals()
    center=at(g0);plus=at(g0+step);minus=at(g0-step)
    fd=(plus['moments']-minus['moments'])/(2*step)
    relative=float(np.max(abs(fd-center['derivatives'])/(1+abs(center['derivatives']))))
    assert relative<2e-8
    omitted_boundary=float(np.trace(center['Q']@center['flux']).real)
    omitted_geometry=[float(np.trace(center['Qp']@y).real) for y in center['outputs']]
    assert abs(omitted_boundary)>.1 and max(abs(x) for x in omitted_geometry)>1
    # The old mixed term must be retained even in the declared finite form.
    _,lam,_=wall.base_data()
    no_mixed=d['Cs'].conj().T@d['Cs']+lam*lam*d['Ch'].conj().T@d['Ch']
    mixed=-lam*(d['Cs'].conj().T@d['Ch']+d['Ch'].conj().T@d['Cs'])
    identity=float(np.max(abs(d['N']-no_mixed-mixed)))
    omitted_mixed_norm=float(np.linalg.norm(mixed,2))
    assert identity<1e-10 and omitted_mixed_norm>1
    eig=np.linalg.eigvalsh(d['N'])
    assert eig[0]>-1e-10
    return dict(original652_H_dimension=64,finite_normal_form_min_eigenvalue=float(eig[0]),
                form_discretization_not_signed_selfadjoint_velocity=True,
                branch_probabilities=center['probs'].tolist(),branch_normal_moments=center['moments'].tolist(),
                joint_source_derivatives=center['derivatives'].tolist(),
                independent_source_differences=fd.tolist(),source_relative_error=relative,
                omitted_branch_boundary_source=omitted_boundary,omitted_explicit_geometry_sources=omitted_geometry,
                original_mixed_identity_error=identity,omitted_mixed_operator_norm=omitted_mixed_norm,
                W_from_declared_fixed_neighbor_diagnostic=True,
                no_continuum_backaction_identity_asserted_for_finite_difference_H=True)


def run():
    result=dict(original_global_Hardy=hardy_energy_check(),
                actual_region_normal_backaction=normal_record_check(),joint_normal_sources=joint_source_check())
    deps=('research_note_556.md','research_note_563.md','research_note_652.md','research_note_722.md',
          'research_note_723.md','research_note_724.md','joint_quantum_reference_forms.py',
          'joint_material_region_records.py','joint_reference_process_transport.py',
          'round725_drafts/wall_domain_entry.py','round725_drafts/wall_domain_entry_results.json')
    return dict(round=725,tests_run=3,failures=0,errors=0,results=result,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original f=s-lambda h on the complete curved target: original-core closed derivative yields a positive square form controlled by the full fixed-graph energy. Same actual region records retain finite normal moments and an exact negative Q backaction; original matching, finite process approximation and first external-geometry sources retain these terminal forms. Does not prove a self-adjoint signed wall velocity, an energy-preserving Q-read instrument, sharp classical signature, continuum quantum normal or Einstein gluing.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
