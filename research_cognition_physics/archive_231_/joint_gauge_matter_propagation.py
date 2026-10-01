"""568: the same Gauss-Higgs matter, transverse propagation and common scale.

Classical quadratic audit of the actual graph H; magnetic faces and matching
remain explicit inputs. This is not a nonperturbative quantum continuum limit.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_scalar_propagation_matching as scalar

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauge_matter_propagation_results.json'


def mul(q,r):
    return np.r_[q[0]*r[0]-q[1:]@r[1:], q[0]*r[1:]+r[0]*q[1:]+np.cross(q[1:],r[1:])]


def conj(q):
    return q*np.array([1.,-1.,-1.,-1.])


def exp_half(x):
    length=float(np.linalg.norm(x))
    return np.r_[np.cos(length/2),(.5 if length==0 else np.sin(length/2)/length)*x]


def square():
    edges=[(0,1),(1,2),(3,2),(0,3)]
    D=np.zeros((4,4))
    for e,(i,j) in enumerate(edges):D[e,i]=1;D[e,j]=-1
    faces=[[(0,1),(1,1),(2,-1),(3,-1)]]
    return D,np.array([[1.,1.,-1.,-1.]]),edges,faces


def torus(length,dimension):
    points=list(itertools.product(range(length),repeat=dimension));lookup={p:i for i,p in enumerate(points)}
    def shift(p,a):
        q=list(p);q[a]=(q[a]+1)%length;return tuple(q)
    N=len(points);edges=[]; labels={}
    for p in points:
        for a in range(dimension):
            labels[(p,a)]=len(edges);edges.append((lookup[p],lookup[shift(p,a)]))
    D=np.zeros((len(edges),N));faces=[]
    for e,(i,j) in enumerate(edges):D[e,i]=1;D[e,j]=-1
    for p in points:
        for a in range(dimension):
            for b in range(a+1,dimension):
                faces.append([(labels[p,a],1),(labels[shift(p,a),b],1),
                    (labels[shift(p,b),a],-1),(labels[p,b],-1)])
    C=np.zeros((len(faces),len(edges)))
    for f,face in enumerate(faces):
        for e,sign in face:C[f,e]+=sign
    return D,C,edges,faces,points


def full_potential(fields,links,edges,faces,L,u,volume,kh,ks,beta):
    X,s=fields[:,:4],fields[:,4]
    d=np.column_stack((np.sum(X*X,axis=1)-u[0],s*s-u[1]))
    value=volume*np.einsum('ni,ij,nj->',d,L,d)/4
    for e,(i,j) in enumerate(edges):
        value+=kh*np.sum((X[i]-mul(links[e],X[j]))**2)/2+ks*(s[i]-s[j])**2/2
    for face in faces:
        q=np.array([1.,0.,0.,0.])
        for e,sign in face:q=mul(q,links[e] if sign==1 else conj(links[e]))
        value+=beta*(1-q[0])
    return float(value)


def nonlinear_hessian_check():
    L,u,mass=scalar.parameters();D,C,edges,faces=square();N=4;E=4
    volume=1.;kh=ks=.8;beta=1.6;h0=np.sqrt(u[0]);s0=np.sqrt(u[1])
    def fun(z):
        radial=z[:2*N].reshape(N,2)+[h0,s0]
        phi=z[2*N:5*N].reshape(N,3);theta=z[5*N:].reshape(E,3)
        fields=np.column_stack((np.array([radial[i,0]*exp_half(phi[i]) for i in range(N)]),radial[:,1]))
        links=np.array([exp_half(t) for t in theta])
        return full_potential(fields,links,edges,faces,L,u,volume,kh,ks,beta)
    B=np.column_stack((-D,np.eye(E)))
    angular=kh*h0*h0/4*np.kron(B.T@B,np.eye(3))
    angular[3*N:,3*N:]+=beta/4*np.kron(C.T@C,np.eye(3))
    radial=volume*np.kron(np.eye(N),mass)+np.kron(D.T@D,np.diag([kh,ks]))
    expected=np.zeros((5*N+3*E,5*N+3*E))
    expected[:2*N,:2*N]=radial;expected[2*N:,2*N:]=angular
    errors=[]
    for step in (1e-3,5e-4):
        found=scalar.hessian(fun,np.zeros(len(expected)),step)
        errors.append(float(np.max(abs(found-expected))))
    assert max(errors)<3e-7 and errors[-1]<errors[0]/2
    return dict(Cartesian_and_exact_SU2_full_potential_variables=len(expected),
        Hessian_refinement_errors=errors,
        radial_angle_cross_block_max=float(np.max(abs(found[:2*N,2*N:]))),
        convention='half-angle Lie basis, beta*(1-ReTr(U_face)/2)')


def symplectic_check():
    _,u,_=scalar.parameters();h2=u[0];v=1.3;kh=.8;b=.7;beta=1.1
    D,C,_,_=square();N=D.shape[1];E=D.shape[0]
    rng=np.random.default_rng(568);phi=rng.normal(size=N);theta=rng.normal(size=E);pi=rng.normal(size=E)
    p=-D.T@pi;a=theta-D@phi;mu=kh*h2/4
    A=2*b*np.eye(E)+4/(v*h2)*(D@D.T);K=mu*np.eye(E)+beta/4*(C.T@C)
    full=2/(v*h2)*(p@p)+b*(pi@pi)+mu/2*(a@a)+beta/8*np.sum((C@theta)**2)
    reduced=.5*pi@A@pi+.5*a@K@a
    dphi=rng.normal(size=N);da=rng.normal(size=E);dtheta=da+D@dphi
    symplectic=float(abs(p@dphi+pi@dtheta-pi@da))
    pdot=mu*D.T@a;pidot=-mu*a-beta/4*C.T@C@theta
    gauss=float(np.max(abs(pdot+D.T@pidot)))
    adot=2*b*pi-D@(4/(v*h2)*p)
    reduced_eq=max(float(np.max(abs(adot-A@pi))),float(np.max(abs(pidot+K@a))))
    assert max(abs(full-reduced),symplectic,gauss,reduced_eq)<2e-13
    return dict(energy_residual=abs(full-reduced),canonical_one_form_residual=symplectic,
        Gauss_derivative_residual=gauss,reduced_Hamilton_equation_residual=reduced_eq)


def hodge_check():
    rows=[]
    for d in (2,3):
        D,C,edges,faces,points=torus(3,d);N=len(points);E=len(edges)
        L0=D.T@D;L1=D@D.T+C.T@C
        identity=float(np.max(abs(L1-np.kron(L0,np.eye(d)))))
        assert np.max(abs(C@D))==0 and identity==0
        harmonic=int(np.sum(np.linalg.eigvalsh(L1)<1e-9))
        rankD=int(np.linalg.matrix_rank(D));rankC=int(np.linalg.matrix_rank(C))
        assert rankD==N-1 and harmonic==d and rankD+rankC+harmonic==E
        rows.append(dict(dimension=d,sites=N,links=E,faces=len(faces),
            rank_gradient=rankD,rank_curl=rankC,harmonic_modes=harmonic,
            one_form_Hodge_identity_residual=identity))
    return dict(tori=rows,dimension_not_selected=True)


def dispersion_check():
    _,u,_=scalar.parameters();h2=u[0];v=1.;kh=.8;b=1.;mass2=b*kh*h2/2
    D,C,edges,faces,points=torus(3,3);E=len(edges);N=len(points)
    A=2*b*np.eye(E)+4/(v*h2)*D@D.T;mu=kh*h2/4
    lapvals=np.linalg.eigvalsh(D.T@D);positive=lapvals[lapvals>1e-10]
    rows=[]
    for beta in (0.,.8,1.6):
        K=mu*np.eye(E)+beta/4*C.T@C;omega=A@K
        formula=mass2*np.eye(E)+(kh/v)*D@D.T+(b*beta/2)*C.T@C
        residual=float(np.max(abs(omega-formula)))
        symmetry=float(np.max(abs(omega-omega.T)))
        expected=np.sort(np.r_[np.repeat(mass2,3),mass2+(kh/v)*positive,
            np.repeat(mass2+(b*beta/2)*positive,2)])
        spectral=float(np.max(abs(np.linalg.eigvalsh((omega+omega.T)/2)-expected)))
        assert residual<1e-13 and symmetry<1e-13 and spectral<3e-13
        rows.append(dict(beta=beta,mass_squared=mass2,longitudinal_lattice_slope=kh/v,
            transverse_lattice_slope=b*beta/2,dynamic_matrix_residual=residual,spectral_residual=spectral,
            number_at_mass_gap=int(np.sum(abs(expected-mass2)<1e-9))))
    return dict(one_colour_modes=E,three_colours_have_same_quadratic_operator=True,branches=rows)


def nonlinear_gauge_check():
    L,u,_=scalar.parameters();D,C,edges,faces=square();rng=np.random.default_rng(5681)
    errors=[];positivity=[]
    for _ in range(12):
        fields=rng.normal(size=(4,5));q=rng.normal(size=(8,4));q/=np.linalg.norm(q,axis=1)[:,None]
        g,U=q[:4],q[4:];new=fields.copy()
        for i in range(4):new[i,:4]=mul(g[i],fields[i,:4])
        Ug=np.array([mul(mul(g[i],U[e]),conj(g[j])) for e,(i,j) in enumerate(edges)])
        original=full_potential(fields,U,edges,faces,L,u,1.,.8,.8,1.6)
        changed=full_potential(new,Ug,edges,faces,L,u,1.,.8,.8,1.6)
        errors.append(abs(original-changed));positivity.append(original)
    assert max(errors)<4e-13 and min(positivity)>0
    return dict(full_noncommuting_SU2_gauge_residual=max(errors),minimum_sampled_positive_potential=min(positivity),
        positivity_and_invariance_proved_algebraically=True)


def scales_check():
    _,u,_=scalar.parameters();c2=.8;g2=2.;rows=[]
    for d in (2,3,4):
        for eps in (1.,.5,.25):
            v=eps**d;kh=c2*eps**(d-2);b=g2/(2*eps**(d-2));beta=4*c2*eps**(d-4)/g2
            ch2=kh*eps**2/v;cg2=b*beta*eps**2/2;mass2=b*kh*u[0]/2
            assert abs(ch2-cg2)<1e-14 and abs(ch2-c2)<1e-14
            assert abs(mass2-c2*g2*u[0]/4)<1e-14
            rows.append(dict(d=d,epsilon=eps,volume=v,k_h=kh,b=b,beta=beta,
                common_speed_squared=ch2,vector_frequency_gap_squared=mass2))
    return dict(rows=rows,coupling_is_normalization_input_not_measured_weak_coupling=True,
        fixed_b_fast_link_limit_is_not_this_scale_family=True)


def joint_boson_check():
    _,u,mass=scalar.parameters();c2=.8;v=1.;kh=ks=.8;b=1.;beta=1.6;mA2=b*kh*u[0]/2
    radial_m2=np.linalg.eigvalsh(mass);rows=[]
    for momentum in ([0.,0.,0.],[.3,.7,1.1],[1.,.5,.2]):
        p=np.array(momentum);dvec=1-np.exp(1j*p);lam=float((dvec.conj()@dvec).real)
        # Direct curl rows: C_ij theta = d_j theta_i - d_i theta_j.
        C=np.zeros((3,3),dtype=complex)
        for f,(i,j) in enumerate(((0,1),(0,2),(1,2))):C[f,i]=dvec[j];C[f,j]=-dvec[i]
        D=dvec[:,None];A=2*b*np.eye(3)+4/(v*u[0])*(D@D.conj().T)
        K=kh*u[0]/4*np.eye(3)+beta/4*C.conj().T@C
        vec=A@K;rad=mass+lam*np.diag([kh/v,ks/v])
        assert np.max(abs(vec-(mA2+c2*lam)*np.eye(3)))<2e-14
        assert np.max(abs(np.linalg.eigvalsh(rad)-(radial_m2+c2*lam)))<2e-14
        rows.append(dict(dimensionless_momentum=momentum,lattice_Laplacian=lam,
            radial_omega_squared=np.linalg.eigvalsh(rad).tolist(),
            each_colour_vector_omega_squared=np.linalg.eigvalsh((vec+vec.conj().T)/2).tolist()))
    return dict(rows=rows,common_slope=c2,total_modes_per_momentum_in_d3=11,
        scalar_and_three_massive_SU2_vectors_only=True,
        no_photon_hypercharge_fermions_or_gravity_added=True)


def smooth_energy_check():
    # A fixed smooth periodic family. Each A_i is constant along its own link,
    # so the link exponential is the exact parallel transport for this family.
    L,u,_=scalar.parameters();h,s0=np.sqrt(u);g=np.sqrt(2.);c2=.8
    def product(q,r):
        return np.concatenate(((q[...,0]*r[...,0]-np.sum(q[...,1:]*r[...,1:],axis=-1))[...,None],
            q[...,0,None]*r[...,1:]+r[...,0,None]*q[...,1:]+np.cross(q[...,1:],r[...,1:])),axis=-1)
    rows=[]
    for count in (8,16,32):
        eps=2*np.pi/count;v=eps**3;k=c2*eps;beta=4*c2/(g*g*eps)
        axis=np.arange(count)*eps;x,y,z=np.meshgrid(axis,axis,axis,indexing='ij')
        r=h*(1+.12*np.cos(x)+.07*np.sin(y));s=s0*(1+.1*np.sin(z))
        amplitudes=[.11*(1+.2*np.sin(y)),.09*(1+.15*np.cos(z)),.08*(1+.1*np.cos(x))]
        X=np.zeros(r.shape+(4,));X[...,0]=r
        links=[]
        for i,amplitude in enumerate(amplitudes):
            U=np.zeros(r.shape+(4,));U[...,0]=np.cos(g*eps*amplitude/2)
            U[...,i+1]=np.sin(g*eps*amplitude/2);links.append(U)
        d=np.stack((r*r-u[0],s*s-u[1]),axis=-1)
        potential=np.einsum('...i,ij,...j->...',d,L,d)/4
        lattice=float(v*np.sum(potential));magnetic=0.
        for i in range(3):
            diff=X-product(links[i],np.roll(X,-1,axis=i))
            lattice+=k/2*float(np.sum(diff*diff))+k/2*float(np.sum((s-np.roll(s,-1,axis=i))**2))
            for j in range(i+1,3):
                loop=product(product(product(links[i],np.roll(links[j],-1,axis=i)),
                    np.roll(links[i],-1,axis=j)*[1,-1,-1,-1]),links[j]*[1,-1,-1,-1])
                magnetic+=beta*float(np.sum(1-loop[...,0]))
        lattice+=magnetic
        scalar_gradient=(.12*h*np.sin(x))**2+(.07*h*np.cos(y))**2+(.1*s0*np.cos(z))**2
        covariant_part=r*r*g*g/4*sum(a*a for a in amplitudes)
        curvature=(.022*np.cos(y))**2+(.008*np.sin(x))**2+(.0135*np.sin(z))**2
        curvature+=g*g*sum((amplitudes[i]*amplitudes[j])**2 for i in range(3) for j in range(i+1,3))
        # All terms are trigonometric polynomials of degree <=4 in each variable.
        # These uniform periodic quadratures integrate the continuum density exactly.
        continuum=float(v*np.sum(potential+c2/2*(scalar_gradient+covariant_part+curvature)))
        rows.append(dict(points_per_axis=count,epsilon=eps,full_nonlinear_lattice_energy=lattice,
            continuum_density_integral=continuum,absolute_difference=abs(lattice-continuum),
            magnetic_energy=magnetic))
    assert max(abs(r['continuum_density_integral']-rows[0]['continuum_density_integral']) for r in rows)<1e-12
    assert rows[-1]['absolute_difference']<rows[0]['absolute_difference']/12
    return dict(rows=rows,noncommuting_gauge_directions=True,
        scope='smooth classical energy consistency, not evolution or quantum convergence')


def run():
    checks=[nonlinear_hessian_check,symplectic_check,hodge_check,dispersion_check,
            nonlinear_gauge_check,scales_check,joint_boson_check,smooth_energy_check]
    evidence={f.__name__:f() for f in checks}
    deps=['joint_scalar_propagation_matching.py','joint_finite_time_gauge_probe.py',
          'joint_singlet_common_mass_rg_results.json','research_note_566.md','research_note_567.md']
    return dict(round=568,tests_run=len(checks),failures=0,errors=0,checks=[f.__name__ for f in checks],
        evidence=evidence,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(classical_quadratic_Gauss_reduction=True,
            original_longitudinal_matter_mediated_propagation_preserved=True,
            magnetic_plaquettes_and_common_speed_are_added_inputs=True,
            same_radial_potential_vacuum_and_mass_matrix=True,
            no_quantum_phase_or_continuum_limit_or_full_Lorentz_proof=True,
            no_dimension_GR_or_standard_model_completion=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False))
