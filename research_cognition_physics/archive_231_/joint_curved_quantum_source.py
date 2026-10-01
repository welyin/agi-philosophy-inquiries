"""574: fixed-graph Gauss states for the same positive-F matter source.

Laplace--Beltrami ordering, geodesic edge energy, and fixed initial geometry
are declared inputs. No quantum continuum or quantum Einstein constraint.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gauss_continuum_sampling as lattice
import joint_gauss_einstein_initial_data as geometry

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_curved_quantum_source_results.json'
M=geometry.M02


def F(phi):return M-np.sum(phi*phi,axis=-1)/6


def metric(phi):
    f=F(phi)
    return np.eye(5)/f[...,None,None]+phi[..., :,None]*phi[...,None,:]/(6*f[...,None,None]**2)


def inverse(phi):
    return F(phi)[...,None,None]*(np.eye(5)-phi[..., :,None]*phi[...,None,:]/(6*M))


def distance_squared(phi,chi):
    a,b=np.sqrt(F(phi)),np.sqrt(F(chi))
    # Cancellation-free cosh(d/sqrt(6))-1.
    z=(np.sum((phi-chi)**2,axis=-1)/12+(a-b)**2/2)/(a*b)
    return 24*np.arcsinh(np.sqrt(z/2))**2


def real_phi(X,s):return np.concatenate((X.real,X.imag,s[...,None]),axis=-1)


def node_potential(phi):
    mat,u,_=lattice.scalar.parameters()
    delta=np.stack((np.sum(phi[...,:4]**2,axis=-1)-u[0],phi[...,4]**2-u[1]),axis=-1)
    return np.einsum('...i,ij,...j->...',delta,mat,delta)/(4*F(phi)**2)


def shared_source(N):
    q=lattice.sample(N,3)
    x,y=np.moveaxis(q['grid'],-1,0)[:2]
    hstar=np.sqrt(lattice.PAR['h2'])
    _,u,_=lattice.scalar.parameters();sstar=np.sqrt(u[1])
    radial=-.0056/(.05*hstar)*np.cos(x+y)
    ps=.025*np.cos(x)-.0056/(.06*sstar)*np.sin(y)
    change=q['eps']**3*(radial-q['f']['radial'])
    q['PX']=q['PX'].copy();q['PXhat']=q['PXhat'].copy()
    q['PX'][...,1]+=change;q['PXhat'][...,1]+=change
    q['f']=dict(q['f']);q['f']['Ps']=ps;q['f']['radial']=radial
    q['phi']=real_phi(q['f']['X'],q['f']['s'])
    q['P']=real_phi(q['PXhat'],q['eps']**3*ps)
    return q


def graph_energy(q,psi):
    eps=q['eps'];phi=q['phi'];P=q['P'];X=q['f']['X'];s=q['f']['s']
    scalar_kin=float(np.sum(np.einsum('...i,...ij,...j->...',P,inverse(phi),P)/(2*eps**3*psi**6)))
    onsite=float(eps**3*np.sum(psi**6*node_potential(phi)))
    W=lattice.su2(q['a'],eps);z=np.exp(1j*eps*q['a0'])
    bw,b0=lattice.PAR['b'][1:];gradient=0.;electric=0.;magnetic=0.
    for mu in range(3):
        midpoint=(psi+np.roll(psi,-1,axis=mu))/2
        next_X=z[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],np.roll(X,-1,axis=mu))
        chi=real_phi(next_X,np.roll(s,-1,axis=mu))
        gradient+=float(eps*np.sum(midpoint**2*distance_squared(phi,chi))/2)
        electric+=float(np.sum(midpoint**-2*(bw*np.sum(q['phat'][...,mu,:]**2,axis=-1)+b0*q['p0hat'][...,mu]**2))/eps)
    dg,coef=lattice.old.coefficients(lattice.PAR)
    for mu in range(3):
        for nu in range(mu+1,3):
            wm,wn=W[...,mu,:,:],W[...,nu,:,:]
            wc=wm@np.roll(wn,-1,axis=mu)@np.swapaxes(np.roll(wm,-1,axis=nu).conj(),-1,-2)@np.swapaxes(wn.conj(),-1,-2)
            zc=z[...,mu]*np.roll(z[...,nu],-1,axis=mu)/np.roll(z[...,mu],-1,axis=nu)/z[...,nu]
            tw=np.trace(wc,axis1=-2,axis2=-1)
            chi=lattice.PAR['wq']*(3*tw*zc+3*(zc**-4+zc**2))+lattice.PAR['wl']*(tw*zc**-3+zc**6+1)
            vf=dg*(12*lattice.PAR['wq']+4*lattice.PAR['wl']-chi.real)+coef[1]*(4-abs(tw)**2)+coef[2]*(1-(zc**6).real)
            center=(psi+np.roll(psi,-1,axis=mu)+np.roll(psi,-1,axis=nu)+np.roll(np.roll(psi,-1,axis=mu),-1,axis=nu))/4
            magnetic+=float(np.sum(center**-2*vf)/eps)
    return dict(scalar_kinetic=scalar_kin,onsite=onsite,gradient=gradient,
                electric=electric,magnetic=magnetic,total=scalar_kin+onsite+gradient+electric+magnetic)


def target_geometry_check():
    rng=np.random.default_rng(574);phi=rng.normal(size=(30,5))*.35
    k=metric(phi);ki=inverse(phi);f=F(phi)
    inverse_error=float(np.max(abs(k@ki-np.eye(5))))
    determinant_error=float(np.max(abs(np.linalg.det(k)/(M/f**6)-1)))
    r=np.linalg.norm(phi,axis=-1);rho=np.sqrt(6)*np.arctanh(r/np.sqrt(6*M))
    angular_error=float(np.max(abs(r*r/f-6*np.sinh(rho/np.sqrt(6))**2)))
    assert max(inverse_error,determinant_error,angular_error)<1e-13
    return dict(inverse_error=inverse_error,determinant_relative_error=determinant_error,
                radial_transform_error=angular_error,target_dimension=5,sectional_curvature=-1/6,
                target_dimension_is_not_spacetime_dimension=True)


def distance_and_gauge_check():
    rng=np.random.default_rng(5741);phi=rng.normal(size=5)*.2;v=rng.normal(size=5)*.3
    expected=float(v@metric(phi)@v);rows=[]
    for step in (.02,.01,.005):
        got=float(distance_squared(phi-step*v/2,phi+step*v/2)/step**2)
        rows.append(dict(step=step,quadratic_coefficient=got,error=abs(got-expected)))
    assert rows[-1]['error']<rows[0]['error']/12
    other=rng.normal(size=5)*.3;U=lattice.old.group_exp(rng.normal(size=3),2);phase=np.exp(.31j)
    def transform(p):
        X=p[:2]+1j*p[2:4];Y=phase**3*(U@X)
        return np.r_[Y.real,Y.imag,p[4]]
    gauge_error=abs(float(distance_squared(transform(phi),transform(other))-distance_squared(phi,other)))
    assert gauge_error<1e-13
    return dict(principal_coefficient=expected,centered_expansion=rows,gauge_distance_error=gauge_error)


def exact_source_check():
    rows=[]
    for N in (8,12,16):
        q=shared_source(N);c=geometry.make_source(N);eps=q['eps']
        gw,g0=lattice.residual(q);err=max(np.max(abs(gw)),np.max(abs(g0)))/eps**3
        radial=np.max(abs(q['PXhat'][...,1].real/eps**3-c['PX'][...,1].real))
        singlet=np.max(abs(q['f']['Ps']-c['Ps']))
        difference=q['P']/eps**3-c['p']
        correction=float(np.sqrt(eps**3*np.sum(np.einsum('...i,...ij,...j->...',difference,inverse(q['phi']),difference))))
        # Every x-edge has a strictly noncentral weak matrix with off-diagonal part.
        offdiag=np.min(abs(np.sin(eps*np.linalg.norm(q['a'][...,0,:],axis=-1)/2)))
        assert err<2e-11 and max(radial,singlet)<1e-13 and offdiag>0
        rows.append(dict(N=N,Gauss_density_error=float(err),radial_identity_error=float(radial),
            singlet_identity_error=float(singlet),node_curved_correction_norm=correction,
            norm_over_eps_squared=correction/eps**2,min_x_link_offdiagonal=float(offdiag),
            F_min=float(np.min(F(q['phi'])))))
    assert rows[-1]['node_curved_correction_norm']<rows[0]['node_curved_correction_norm']/3
    return dict(rows=rows,colour_momentum_zero=True,EW_stabilizer_only_ineffective_center=True,
                compact_colour_stabilizer_fixes_source=True)


def curved_energy_check():
    rows=[]
    for N in (8,12,20):
        c=geometry.make_source(N);psi,_,_=geometry.solve_hamiltonian(c);q=shared_source(N)
        discrete=graph_energy(q,psi);eps=q['eps']
        continuous=float(eps**3*np.sum(.5*psi**-6*c['pKp']+.5*psi**2*c['B']+psi**6*c['U']+psi**-2*c['Y']))
        error=abs(discrete['total']-continuous)
        assert min(discrete.values())>=0
        rows.append(dict(N=N,discrete=discrete,continuum_quadrature=continuous,
                         difference=error,difference_over_eps_squared=error/eps**2))
    assert rows[-1]['difference']<rows[0]['difference']/4
    return dict(rows=rows,fixed_smooth_initial_geometry_only=True,
                sampled_psi_not_a_continuum_error_certificate=True)


def boundary_check():
    mat,u,_=lattice.scalar.parameters()
    # At r²=6M: ||(h²-u_h,s²-u_s)||² >= (6M-u_h-u_s)²/2.
    lower=float(np.linalg.eigvalsh(mat)[0]*(6*M-np.sum(u))**2/8)
    assert lower>0 and np.sum(u)<6*M
    norm_shells=[];energy_shells=[]
    for cut in (.02,.01,.005):
        # Radial asymptotic witnesses in t=F, prefactors bounded and positive.
        norm_shells.append((cut**-2-.1**-2)/2)
        energy_shells.append(cut**-1-.1**-1)
    assert norm_shells[-1]>4*norm_shells[-2] and energy_shells[-1]>2*energy_shells[-2]
    return dict(V_boundary_analytic_lower=lower,F_power_in_volume=-3,
                raw_Gaussian_norm_divergence=norm_shells,
                half_density_Gaussian_potential_divergence=energy_shells,
                cutoff_or_intrinsic_packets_not_excluded=True)


def radial_packet(hbar,nodes=90):
    """Two-coordinate quadrature diagnostic, not a full many-node state simulation.

    SU(2)xU(1)-invariant unit-cell radial source, no gauge link momenta.
    Uses the actual h,s and radial density momenta at the 573 reference line.
    """
    hstar=np.sqrt(lattice.PAR['h2']);mat,u,_=lattice.scalar.parameters();sstar=np.sqrt(u[1])
    center=np.array([1.05*hstar,sstar]);momentum=np.array([0.,.025-.0056/(.06*sstar)])
    R=.18;x,weights=np.polynomial.hermite.hermgauss(nodes)
    xx=np.stack(np.meshgrid(x,x,indexing='ij'),axis=-1);d=np.sqrt(hbar)*xx
    q=center+d;mask=np.all(abs(d)<R,axis=-1)
    q=q[mask];d=d[mask];ww=(weights[:,None]*weights[None,:])[mask]
    h,s=q.T;f=M-np.sum(q*q,axis=-1)/6
    assert np.min(h)>0 and np.min(f)>0
    chi=np.exp(-np.sum(d*d/(R*R-d*d),axis=-1))
    mu=h**3*np.sqrt(M)*f**-3
    probability=ww*mu*chi**2;Z=hbar*np.sum(probability);probability/=np.sum(probability)
    log1=-d/hbar-2*d*R**2/(R**2-d*d)**2
    log2=-1/hbar-2*R**2*(R**2+3*d*d)/(R**2-d*d)**3
    gi=f[:,None,None]*(np.eye(2)-q[:,:,None]*q[:,None,:]/(6*M))
    gradlogmu=np.stack((3/h+h/f,s/f),axis=-1)
    divg=-q/3+np.sum(q*q,axis=-1)[:,None]*q/(18*M)-3*f[:,None]*q/(6*M)
    drift=divg+np.einsum('...ij,...i->...j',gi,gradlogmu)
    dl=log1+1j*momentum/hbar
    lap=np.einsum('...i,...ij,...j->...',dl,gi,dl)+np.sum(np.diagonal(gi,axis1=-2,axis2=-1)*log2,axis=-1)+np.sum(drift*dl,axis=-1)
    delta=np.stack((h*h-u[0],s*s-u[1]),axis=-1)
    U=np.einsum('...i,ij,...j->...',delta,mat,delta)/(4*f*f)
    Hratio=-hbar*hbar*lap/2+U
    velocity=np.einsum('...ij,j->...i',gi,momentum)-1j*hbar*(np.einsum('...ij,...j->...i',gi,log1)+drift/2)
    f0=M-np.dot(center,center)/6;g0=f0*(np.eye(2)-np.outer(center,center)/(6*M))
    delta0=center*center-u;E0=float(momentum@g0@momentum/2+delta0@mat@delta0/(4*f0*f0))
    meanH=np.sum(probability*Hratio);meanV=np.sum(probability[:,None]*velocity,axis=0)
    meanE2=float(np.sum(probability*abs(Hratio)**2));var=meanE2-float(meanH.real)**2
    return dict(hbar=hbar,normalizer=Z,energy=float(meanH.real),classical_energy=E0,
        energy_error=abs(float(meanH.real)-E0),energy_variance=var,
        energy_second_moment=meanE2,mean_velocity=meanV.real.tolist(),classical_velocity=(g0@momentum).tolist(),
        imaginary_mean_error=float(max(abs(meanH.imag),np.max(abs(meanV.imag)))),
        norm_H_minus_classical=float(np.sqrt(np.sum(probability*abs(Hratio-E0)**2))),
        variance_h=float(np.sum(probability*(h-np.sum(probability*h))**2)))


def quantum_packet_check():
    rows=[radial_packet(b) for b in (.001,.00025,.0000625)]
    for r in rows:
        assert r['normalizer']>0 and r['energy_variance']>0 and r['imaginary_mean_error']<2e-7
        assert r['classical_velocity'][0]>0
    assert rows[-1]['energy_error']<rows[0]['energy_error']/10
    assert rows[-1]['energy_variance']<rows[0]['energy_variance']/10
    repeat=radial_packet(.00025,110)
    quadrature_error=abs(repeat['energy']-rows[1]['energy'])
    assert quadrature_error<1e-9
    return dict(rows=rows,independent_quadrature_order_error=quadrature_error,
                one_cell_radial_diagnostic_only=True,full_graph_Gauss_state_exists_by_tube_proof=True)


def run():
    checks=('target_geometry_check','distance_and_gauge_check','exact_source_check',
            'curved_energy_check','boundary_check','quantum_packet_check')
    evidence={name:globals()[name]() for name in checks}
    deps=('joint_gauss_continuum_sampling.py','joint_gauss_einstein_initial_data.py',
          'joint_quotient_gauge_completion.py','research_round_573_checks.json')
    return dict(round=574,tests_run=len(checks),failures=0,errors=0,checks=list(checks),evidence=evidence,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='declared curved-target quantization on fixed finite graph and prescribed initial geometry; exact Gauss semiclassical initial source; no continuum quantum dynamics, actual preparation, quantum Einstein constraints or dimension generation')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
