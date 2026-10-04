"""575: an explicit autonomous neutral pointer in the common matter/geometry model.
New circle field and target connection are inputs, not deductions of the SM or GR.
Classical constraints and finite-graph quantum response have separate scopes.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as old
import joint_gravity_material_coordinates as chart

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_autonomous_pointer_geometry_results.json'
geo=old.geometry
lat=old.lattice


def connection(phi,g=1.):
    out=np.zeros_like(phi);out[...,:4]=g*phi[...,4,None]*phi[...,:4]
    return out


def target_blocks(phi,g=1.):
    k=old.metric(phi);ki=old.inverse(phi);a=connection(phi,g);v=np.einsum('...ij,...j->...i',ki,a)
    shape=phi.shape[:-1]+(6,6)
    G=np.zeros(shape);inv=np.zeros(shape)
    G[...,:5,:5]=k+a[..., :,None]*a[...,None,:]
    G[...,:5,5]=G[...,5,:5]=a;G[...,5,5]=1
    inv[...,:5,:5]=ki;inv[...,:5,5]=inv[...,5,:5]=-v
    inv[...,5,5]=1+np.sum(a*v,axis=-1)
    return G,inv


def edge_integral(phi,chi):
    dx=chi[...,:4]-phi[...,:4]
    ds=chi[...,4]-phi[...,4];sm=(chi[...,4]+phi[...,4])/2
    dt=(np.sum(chi[...,:4]**2,axis=-1)-np.sum(phi[...,:4]**2,axis=-1))/2
    return sm*dt+ds*np.sum(dx*dx,axis=-1)/12


def source(N,g=1.):
    q=geo.make_source(N)
    x,y,_=np.moveaxis(q['grid'],-1,0)
    hs=np.sqrt(lat.PAR['h2']);_,u,_=lat.scalar.parameters();ss=np.sqrt(u[1])
    dh=np.zeros_like(q['grid']);ds=np.zeros_like(dh)
    dh[...,0]=dh[...,1]=.05*hs*np.cos(x+y);ds[...,1]=-.06*ss*np.sin(y)
    C=g*q['f']['s'][...,None]*q['f']['h'][...,None]*dh
    extra=np.sum(C*C,axis=-1)
    new=dict(q);new['B']=q['B']+extra
    psi,tensor,stats=geo.solve_hamiltonian(new)
    velocity=np.einsum('...ij,...j->...i',old.inverse(q['phi']),q['p'])*psi[...,None]**-6
    vh,vs=velocity[...,1],velocity[...,4]
    vchi=-g*q['f']['s']*q['f']['h']*vh
    rh=psi**-4*np.sum(dh*dh,axis=-1)-vh*vh
    rs=psi**-4*np.sum(ds*ds,axis=-1)-vs*vs
    return dict(q=new,original=q,psi=psi,tensor=tensor,stats=stats,dh=dh,ds=ds,C=C,extra=extra,
                v=velocity,vchi=vchi,rh=rh,rs=rs,hstar=hs,sstar=ss)


def metric_and_canonical_check():
    rng=np.random.default_rng(575);phi=rng.normal(size=(30,5))*.3
    G,inv=target_blocks(phi)
    err=float(np.max(abs(G@inv-np.eye(6))))
    deter=float(np.max(abs(np.linalg.det(G)/np.linalg.det(old.metric(phi))-1)))
    velocity=rng.normal(size=(30,6));p=np.einsum('...ij,...j->...i',G,velocity)
    split=p[...,:5]-connection(phi)*p[...,5,None]
    Hsplit=(np.einsum('...i,...ij,...j->...',split,old.inverse(phi),split)+p[...,5]**2)/2
    L=np.einsum('...i,...ij,...j->...',velocity,G,velocity)/2
    assert max(err,deter,np.max(abs(Hsplit-L)))<1e-12
    assert np.min(np.linalg.eigvalsh(G))>0
    # Canonical momentum flux: horizontal matter plus pointer equals p_A Dphi + Pchi dchi.
    spatial=rng.normal(size=(30,6));C=spatial[...,5]+np.sum(connection(phi)*spatial[...,:5],axis=-1)
    full=np.sum(split*spatial[...,:5],axis=-1)+p[...,5]*C
    canonical=np.sum(p*spatial,axis=-1)
    assert np.max(abs(full-canonical))<1e-13
    return dict(inverse_error=err,volume_relative_error=deter,Legendre_error=float(np.max(abs(Hsplit-L))),
                canonical_flux_error=float(np.max(abs(full-canonical))),target_dimension=6)


def edge_connection_check():
    rng=np.random.default_rng(5751);a=rng.normal(size=(20,5))*.3;b=rng.normal(size=(20,5))*.3
    exact=edge_integral(a,b);x,w=np.polynomial.legendre.leggauss(5);t=(x+1)/2;w=w/2
    z=a[:,None,:]+t[None,:,None]*(b-a)[:,None,:]
    integral=np.sum(w[None,:]*z[...,4]*np.sum(z[...,:4]*(b-a)[:,None,:4],axis=-1),axis=-1)
    err=float(np.max(abs(exact-integral)));rev=float(np.max(abs(edge_integral(a,b)+edge_integral(b,a))))
    U=lat.old.group_exp(np.array([.3,-.4,.2]),2);phase=np.exp(.23j)
    def gauge(phi):
        X=phi[...,:2]+1j*phi[...,2:4];X=phase**3*np.einsum('ij,...j->...i',U,X)
        return old.real_phi(X,phi[...,4])
    invariance=float(np.max(abs(edge_integral(gauge(a),gauge(b))-exact)))
    assert max(err,rev,invariance)<1e-14
    # Rectangle in target (T,s): integral s dT = minus oriented area.
    t0,t1,s0,s1=.1,.2,.4,.7
    points=[np.array([0.,np.sqrt(2*t),0.,0.,s]) for t,s in ((t0,s0),(t1,s0),(t1,s1),(t0,s1))]
    loop=sum(edge_integral(points[i],points[(i+1)%4]) for i in range(4))
    assert abs(loop+(t1-t0)*(s1-s0))<1e-14
    return dict(line_integral_error=err,reversal_error=rev,gauge_error=invariance,
                target_loop_integral=float(loop),connection_not_an_exact_relabeling=True)


def resource_obstruction_check():
    f=source(16);a=f['C'];N=16;dx=f['q']['dx'];k=geo.waves(N);k2=np.sum(k*k,axis=-1)
    safe=k2.copy();safe[0,0,0]=1
    ah=np.fft.fftn(a,axes=(0,1,2));longitudinal=k*np.sum(k*ah,axis=-1)[...,None]/safe[...,None]
    transverse=np.fft.ifftn(ah-longitudinal,axes=(0,1,2)).real
    norm=float(dx**3*np.sum(transverse**2))
    alpha,beta=.05,.06;H,S=f['hstar'],f['sstar']
    formula=geo.VOL*(S*H*H*alpha)**2*beta**2*(3/20+9*alpha**2/1040)
    assert abs(norm-formula)<1e-16 and formula>0
    harmonic=float(np.max(abs(np.mean(a,axis=(0,1,2)))))
    assert harmonic<1e-15
    lower=f['stats']['lower_barrier']**2*formula/2
    actual=float(dx**3*np.sum(f['psi']**2*np.sum(a*a,axis=-1))/2)
    assert actual>lower>0
    return dict(flat_transverse_norm_squared=norm,analytic_Fourier_norm_squared=formula,
        harmonic_max=harmonic,all_smooth_circle_configs_gradient_energy_lower=lower,
       _chi_zero_initial_gradient_energy=actual,
        winding_harmonics_cannot_cancel_transverse_part=True)


def constraints_and_reference_check():
    oldbounds=chart.rational_bounds()
    extra_bound=(Q('1.06')*Q('.55')*Q('1.05')*Q('.67'))**2*2*(Q('.05')*Q('.67'))**2
    Bmax=Q(oldbounds['B_upper']['exact'])+extra_bound
    assert Bmax<1
    rows=[]
    for N in (16,24,32):
        f=source(N);q=f['q'];psi=f['psi'];base,_,_=geo.solve_hamiltonian(f['original'])
        c=chart.point_coefficients(f);i=(0,N//4,N//8)
        z=geo.derivative(psi,2)
        vh=c['Ch']*psi[i]**-6
        det=8*c['Ch']**2*(-c['b'])*psi[i]**-31*z[i]*c['Chx']*c['b']**2*(c['R']-psi[i]**8)
        rx,rz=geo.derivative(f['rh'],0)[i],geo.derivative(f['rh'],2)[i]
        sx,sz=geo.derivative(f['rs'],0)[i],geo.derivative(f['rs'],2)[i]
        direct=vh*c['b']*(rx*sz-rz*sx)
        relative=float(abs(direct-det)/abs(det))
        assert det>0 and relative<.012 and np.max(z[:,:,1:N//4])<0
        assert np.min(psi-base)>0 and np.max(psi)<1.5
        assert f['stats']['original_Hamiltonian_residual']<3e-8
        assert np.max(abs(q['mom']-f['original']['mom']))==0
        rows.append(dict(N=N,max_psi_change=float(np.max(psi-base)),min_psi_change=float(np.min(psi-base)),
            constraint_residual=f['stats']['original_Hamiltonian_residual'],psi_at_reference=float(psi[i]),
            h_velocity=float(vh),pointer_velocity=float(f['vchi'][i]),clock_norm=float(-vh*vh),
            reference_determinant=float(det),independent_minor_relative_error=relative))
    assert rows[-1]['independent_minor_relative_error']<1e-5
    return dict(extra_B_rational_upper=str(extra_bound),total_B_rational_upper=str(Bmax),rows=rows,
        same_canonical_source_momentum=True,new_geometry_includes_pointer_gradient=True)


def graph_limit_check():
    rows=[]
    for N in (8,12,20):
        f=source(N);psi=f['psi'];q=old.shared_source(N);eps=q['eps'];X=q['f']['X']
        W=lat.su2(q['a'],eps);z=np.exp(1j*eps*q['a0']);graph=0.
        for mu in range(3):
            NX=z[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],np.roll(X,-1,axis=mu))
            phi1=old.real_phi(NX,np.roll(q['f']['s'],-1,axis=mu))
            I=edge_integral(q['phi'],phi1);pm=(psi+np.roll(psi,-1,axis=mu))/2
            graph+=float(eps*np.sum(pm**2*2*np.sin(I/2)**2))
        continuous=float(eps**3*np.sum(psi**2*f['extra'])/2)
        rows.append(dict(N=N,probe_edge_energy=graph,continuum_quadrature=continuous,
                         error=abs(graph-continuous)))
    assert rows[-1]['error']<rows[0]['error']/4
    return dict(rows=rows,new_circle_edge_potential_nonnegative=True,
                no_quantum_continuum_limit_claim=True)


def radial_density(hbar,nodes=70):
    H=np.sqrt(lat.PAR['h2']);mat,u,_=lat.scalar.parameters();S=np.sqrt(u[1])
    center=np.array([1.05*H,S]);p=np.array([0.,.025-.0056/(.06*S)])
    x,w=np.polynomial.hermite.hermgauss(nodes);xx=np.stack(np.meshgrid(x,x,indexing='ij'),axis=-1)
    d=np.sqrt(hbar)*xx;mask=np.all(abs(d)<.18,axis=-1)
    q=(center+d)[mask];d=d[mask];weight=(w[:,None]*w[None,:])[mask]
    h,s=q.T;F=2-(h*h+s*s)/6
    chi=np.exp(-np.sum(d*d/(.18**2-d*d),axis=-1))
    measure=h**3*np.sqrt(2)*F**-3
    prob=weight*chi**2*measure;prob/=np.sum(prob)
    log1=-d/hbar-2*d*.18**2/(.18**2-d*d)**2
    log2=-1/hbar-2*.18**2*(.18**2+3*d*d)/(.18**2-d*d)**3
    ki=F[:,None,None]*(np.eye(2)-q[:,:,None]*q[:,None,:]/12)
    drift=np.stack((3*F/h-F*h/6,-F*s/6),axis=-1)
    a=np.stack((s*h,np.zeros_like(h)),axis=-1)
    mixed=-np.einsum('...ij,...j->...i',ki,a)
    vertical=1-np.sum(a*mixed,axis=-1)
    bchi=-(s*h*drift[:,0]+s*ki[:,0,0]+h*ki[:,1,0])
    delta=q*q-u
    U=np.einsum('...i,ij,...j->...',delta,mat,delta)/(4*F*F)
    return q,p,prob,log1,log2,ki,drift,mixed,vertical,bchi,U


def pointer_packet(hbar,nodes=70,angles=4096):
    q,p,prob,L1,L2,ki,drift,mixed,vertical,bchi,U=radial_density(hbar,nodes)
    theta=np.arange(angles)*2*np.pi/angles;kap=1/hbar
    weight=np.exp(kap*(np.cos(theta)-1));weight/=np.sum(weight)
    lp=-kap*np.sin(theta)/2;lp2=-kap*np.cos(theta)/2
    cosmean=float(weight@np.cos(theta));p2=hbar*hbar*float(weight@(lp*lp))
    dl=L1+1j*p/hbar
    lapbase=np.einsum('...i,...ij,...j->...',dl,ki,dl)+np.sum(np.diagonal(ki,axis1=-2,axis2=-1)*L2,axis=-1)+np.sum(drift*dl,axis=-1)
    Hbase=-hbar*hbar*lapbase/2+U
    cross=2*np.sum(mixed*dl,axis=-1)+bchi
    # Full three-coordinate Laplacian applied to the product wavefunction.
    # Tensor-product quadrature is contracted without a huge dense array.
    meanH=np.sum(prob*Hbase)-hbar*hbar/2*(
        np.sum(prob*vertical)*(weight@(lp*lp+lp2))+np.sum(prob*cross)*(weight@lp))
    meanHsin=np.sum(prob*Hbase)*(weight@np.sin(theta))-hbar*hbar/2*(
        np.sum(prob*vertical)*(weight@((lp*lp+lp2)*np.sin(theta)))+
        np.sum(prob*cross)*(weight@(lp*np.sin(theta))))
    record_dot=float(2*meanHsin.imag/hbar)
    current=float(np.sum(prob*np.sum(mixed*p,axis=-1)))
    expected=cosmean*current
    added=float(p2*np.sum(prob*vertical)/2)
    assert abs(record_dot-expected)<1e-12
    assert abs(p2-hbar*hbar*kap*cosmean/4)<1e-12
    assert abs(float(meanH.real-np.sum(prob*Hbase).real)-added)<1e-12
    center=np.array([1.05*np.sqrt(lat.PAR['h2']),np.sqrt(lat.scalar.parameters()[1][1])])
    f0=2-np.dot(center,center)/6;k0=f0*(np.eye(2)-np.outer(center,center)/12)
    classical=-center[1]*center[0]*float((k0@p)[0])
    return dict(hbar=hbar,cos_mean=cosmean,probe_momentum_variance=p2,
        quantum_sin_response=record_dot,classical_unit_cell_response=classical,
        response_error=abs(record_dot-classical),extra_node_kinetic_energy=added,
        full_Hamiltonian_commutator_error=abs(record_dot-expected),
        binary_effect_response=record_dot/4)


def quantum_record_check():
    rows=[pointer_packet(h) for h in (.002,.0005,.000125)]
    for r in rows:
        assert r['quantum_sin_response']<0 and r['probe_momentum_variance']>0 and r['extra_node_kinetic_energy']>0
    assert rows[-1]['response_error']<rows[0]['response_error']/8
    repeat=pointer_packet(.0005,90,5120)
    err=abs(repeat['quantum_sin_response']-rows[1]['quantum_sin_response'])
    assert err<1e-11
    return dict(rows=rows,independent_quadrature_error=err,
        unit_cell_diagnostic_not_full_graph_evolution=True,
        full_graph_short_time_record_from_bounded_effect_derivative=True,
        uniform_circle_state_has_no_absolute_sin_response=True)


def stress_exchange_check():
    # Independent finite-difference divergence and exact off-shell identity.
    eta=np.array([-1.,1.,1.,1.])
    def fields(x):
        t,a,b,_=x
        kh=np.array([1.,1.,0.,0.]);ks=np.array([1.,0.,-1.,0.]);kc=np.array([1.,1.,1.,0.])
        h=.7+.03*np.sin(t+a);s=.55+.04*np.cos(t-b)
        dh=.03*np.cos(t+a)*kh;ds=-.04*np.sin(t-b)*ks
        dc=.1*np.cos(t+a+b)*kc
        ddh=-.03*np.sin(t+a)*np.outer(kh,kh);ddc=-.1*np.sin(t+a+b)*np.outer(kc,kc)
        C=dc+s*h*dh
        dC=ddc+np.outer(ds*h+s*dh,dh)+s*h*ddh
        tensor=np.outer(C,C)-np.diag(eta)*np.dot(eta*C,C)/2
        return C,dC,tensor
    point=np.array([.2,.4,-.3,.1]);C,dC,_=fields(point);step=1e-5
    divergence=np.zeros(4)
    for mu in range(4):
        e=np.eye(4)[mu]*step
        divergence+=eta[mu]*(fields(point+e)[2][mu]-fields(point-e)[2][mu])/(2*step)
    rhs=float(np.sum(eta*np.diag(dC)))*C+np.einsum('i,ij->j',eta*C,dC-dC.T)
    err=float(np.max(abs(divergence-rhs)))
    assert err<1e-10
    exchange=np.einsum('i,ij->j',eta*C,dC-dC.T)
    assert np.linalg.norm(exchange)>0
    return dict(off_shell_stress_identity_error=err,nonzero_exchange=exchange.tolist(),
                source_and_pointer_stress_must_be_combined=True)


def run():
    checks=('metric_and_canonical_check','edge_connection_check','resource_obstruction_check',
            'constraints_and_reference_check','graph_limit_check','quantum_record_check','stress_exchange_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_curved_quantum_source.py','joint_gauss_einstein_initial_data.py',
          'joint_gravity_material_coordinates.py','research_round_574_checks.json')
    return dict(round=575,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='one explicitly added neutral pointer and autonomous target connection; same matter canonical initial source, joint classical Einstein constraints and finite-graph bounded quantum record; no finite-hbar self-consistent quantum geometry, full apparatus preparation or dimension/GR/SM derivation')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
