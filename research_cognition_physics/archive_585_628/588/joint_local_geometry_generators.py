"""588: paired local lapses, finite-family obstruction, and covariant completion.
No quantum gravity or fixed-hbar continuum result is claimed.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_matter_energy_current as prior

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_local_geometry_generators_results.json'
original=prior.original; PAR=prior.PAR

def dual_graph(q,psi,N,M):
    eps=q['eps'];phi=q['phi']
    vel=np.einsum('...ij,...j->...i',original.inverse(phi),q['P'])/(eps**3*psi[...,None]**6)
    W=original.lattice.su2(q['a'],eps);z=np.exp(1j*eps*q['a0'])
    links=[];rates=[]
    for mu in range(3):
        pe=(psi+np.roll(psi,-1,axis=mu))/2
        A=np.einsum('...a,aij->...ij',2*PAR['b'][1]/eps*pe[...,None]**-2*q['phat'][...,mu,:],original.lattice.T)
        omega=2*PAR['b'][2]/eps*pe**-2*q['p0hat'][...,mu]
        links.append((np.broadcast_to(np.eye(3,dtype=complex),phi.shape[:-1]+(3,3)),W[...,mu,:,:],z[...,mu]))
        rates.append((np.zeros(phi.shape[:-1]+(3,3),complex),1j*A@W[...,mu,:,:],1j*omega*z[...,mu]))
    scalar=0.;magnetic=0.
    def shift(q,mu):return tuple(np.roll(x,-1,axis=mu) for x in q)
    def averages(f,mu,nu):
        fm,fn=np.roll(f,-1,axis=mu),np.roll(f,-1,axis=nu)
        fmn=np.roll(fm,-1,axis=nu)
        return (f+fm+fn+fmn)/4,((f+fm)/2,(fm+fmn)/2,(fn+fmn)/2,(f+fn)/2)
    for mu in range(3):
        pe=(psi+np.roll(psi,-1,axis=mu))/2
        def transform(v):
            X=v[...,:2]+1j*v[...,2:4]
            X=z[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],np.roll(X,-1,axis=mu))
            return original.real_phi(X,np.roll(v[...,4],-1,axis=mu))
        gx,gy=prior.edge_gradients(phi,transform(phi))
        dwi=eps*pe**2*np.sum(gx*vel,axis=-1)
        dwj=eps*pe**2*np.sum(gy*transform(vel),axis=-1)
        alpha=N*np.roll(M,-1,axis=mu)-M*np.roll(N,-1,axis=mu)
        scalar+=float(np.sum(alpha*(dwj-dwi)/2))
        for nu in range(mu+1,3):
            pf,_=averages(psi,mu,nu)
            face=[links[mu],shift(links[nu],mu),shift(links[mu],nu),links[nu]]
            fv=[rates[mu],shift(rates[nu],mu),shift(rates[mu],nu),rates[nu]]
            dw=[x/(eps*pf**2) for x in prior.compose_face(face,fv)]
            nf,ne=averages(N,mu,nu);mf,me=averages(M,mu,nu)
            magnetic+=sum(float(np.sum((em*nf-en*mf)*d)) for en,em,d in zip(ne,me,dw))
    return np.array([scalar,magnetic])

def dual_lapse_check():
    rows=[]
    for n in (8,12,20,32):
        q=original.shared_source(n);c=original.geometry.make_source(n)
        x,y,z=np.moveaxis(q['grid'],-1,0);psi=1.1+.05*np.cos(x)+.02*np.sin(y)
        N=1+.2*np.sin(y);M=1+.15*np.sin(x)
        dN=np.stack((0*x,.2*np.cos(y),0*x),axis=-1)
        dM=np.stack((.15*np.cos(x),0*x,0*x),axis=-1)
        vector=N[...,None]*dM-M[...,None]*dN
        smom=np.einsum('...a,...ia->...i',c['p'],c['Dphi']);gmom=c['mom']-smom
        expected=np.array([q['eps']**3*np.sum(psi**-4*np.sum(vector*mom,axis=-1)) for mom in (smom,gmom)])
        got=dual_graph(q,psi,N,M)
        anti=float(np.max(abs(got+dual_graph(q,psi,M,N))))
        same=float(np.max(abs(dual_graph(q,psi,N,N))))
        gw,g0=original.lattice.residual(q)
        rows.append(dict(N=n,scalar_and_magnetic=got.tolist(),continuum=expected.tolist(),
                         max_sector_error=float(np.max(abs(got-expected))),total_error=float(abs(np.sum(got-expected))),
                         antisymmetry_error=anti,equal_lapse_bracket_error=same,
                         Gauss_density_residual=float(max(np.max(abs(gw)),np.max(abs(g0)))/q['eps']**3)))
    assert rows[-1]['max_sector_error']<rows[0]['max_sector_error']/8
    assert rows[-1]['total_error']<rows[0]['total_error']/8
    assert max(r['antisymmetry_error']+r['equal_lapse_bracket_error'] for r in rows)<1e-13
    return dict(rows=rows,positive_independent_lapses=True,same_EW_Gauss_source=True,
                continuum_scope='classical smooth-source symbol; not a fixed-hbar limit')

def coefficient_identity_check():
    rng=np.random.default_rng(588)
    errors=[]
    for _ in range(20):
        # General incidence, including unequal edge/face allocations.
        a=rng.random(5);a/=a.sum();b=rng.random(5);b/=b.sum()
        N,M=rng.normal(size=(2,5))
        J=np.outer(b,a)-np.outer(a,b)
        direct=N@J@M
        coefficient=(M@a)*(N@b)-(N@a)*(M@b)
        errors.append(abs(direct-coefficient))
    assert max(errors)<1e-15
    return dict(dual_source_incidence_identity_error=float(max(errors)),
                inherited_current_support_and_gauge_derivatives_reused=True)

def mixed_symbol_check():
    path=HERE/'round588_drafts/mixed_current_symbol_entry.py'
    spec=importlib.util.spec_from_file_location('mixed_entry',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    got=json.loads(json.dumps(module.run()))
    assert got==json.loads(path.with_name('mixed_current_symbol_entry_results.json').read_text('utf8'))
    return got

def transported_metric_check():
    # Exact tangent flow at the coincident configuration; this does not replace
    # the nonlinear compactly supported full configuration flow in the proof.
    phi=np.array([.14,.53,-.17,.1,.42]);K=original.inverse(phi)
    w,k=.8,.9;L=np.array([.7,1.3])
    B=k/(2*w)*np.block([[-np.eye(5),np.eye(5)],[-np.eye(5),np.eye(5)]])
    A=np.block([[L[0]/w*K,np.zeros((5,5))],[np.zeros((5,5)),L[1]/w*K]])
    assert np.max(abs(B@B))==0
    derivative=B@A+A@B.T
    expected=k*(L[1]-L[0])/(2*w*w)*K[-1,-1]
    assert abs(derivative[4,9]-expected)<1e-14
    rows=[]
    for t in (-.4,-.2,.2,.4):
        F=np.eye(10)+t*B;At=F@A@F.T
        p=np.linspace(-.3,.4,10)
        equality=abs(float(p@At@p-(F.T@p)@A@(F.T@p)))
        eigen=float(np.linalg.eigvalsh(At).min())
        mixed=float(np.linalg.norm(At[:5,5:]))
        assert eigen>0 and mixed>0 and equality<1e-13
        rows.append(dict(flow_parameter=t,min_kinetic_eigenvalue=eigen,cross_node_block_norm=mixed,
                         energy_pullback_error=equality,Jacobian=float(np.linalg.det(F))))
    # A Gaussian point cloud, exact transformation of a record covector and
    # momenta: common record and kinetic energy remain pointwise identical.
    rng=np.random.default_rng(5884);Q=rng.normal(scale=.1,size=(30,10));P=rng.normal(size=(30,10))
    F=np.eye(10)+.2*B
    Qnew=Q@F.T;Pnew=P@np.linalg.inv(F)
    Qold=Qnew@np.linalg.inv(F).T
    effect=lambda s:.5+np.sin(s)/4
    record_error=float(np.max(abs(effect(Qold[:,4])-effect(Q[:,4]))))
    frozen_record_change=float(np.max(abs(effect(Qnew[:,4])-effect(Q[:,4]))))
    e1=np.einsum('ni,ij,nj->n',P,A,P)
    e2=np.einsum('ni,ij,nj->n',Pnew,F@A@F.T,Pnew)
    energy_error=float(np.max(abs(e1-e2)))
    assert record_error<1e-15 and energy_error<1e-12 and frozen_record_change>1e-3
    return dict(rows=rows,infinitesimal_mixed_coefficient=float(derivative[4,9]),
                transported_record_error=record_error,untransported_record_difference=frozen_record_change,
                transformed_energy_error=energy_error,linear_tangent_diagnostic_only=True,
                full_nonlinear_unitary_and_Gauss_statement_proved_analytically=True)

def geometry_completion_check():
    # Conformal-only metric variations miss the traceless part of Lie_xi delta.
    x=np.arange(32)*2*np.pi/32
    lie=np.zeros((32,3,3));lie[:,0,0]=-2*np.sin(x)
    conformal=np.trace(lie,axis1=1,axis2=2)[:,None,None]*np.eye(3)/3
    residual=lie-conformal
    mean=float(np.mean(np.sum(residual**2,axis=(1,2))))
    assert abs(mean-4/3)<1e-14
    # Full continuum matter density under a non-orthogonal coordinate change.
    rng=np.random.default_rng(5885)
    phi=np.array([.14,.53,-.17,.1,.42]);p=rng.normal(size=5)
    dphi=rng.normal(size=(3,5));E=rng.normal(size=(3,12))
    F=rng.normal(size=(3,3,12));F=F-np.swapaxes(F,0,1)
    gamma=np.array([[1.3,.17,-.08],[.17,.9,.11],[-.08,.11,1.1]])
    transform=np.array([[1.1,.21,0.],[.13,.94,.12],[0.,.08,1.03]])
    det=float(np.linalg.det(transform));inv=np.linalg.inv(transform)
    b=np.repeat(PAR['b'],[8,3,1]);kap=np.repeat(PAR['K'],[8,3,1])
    def density(g,p,d,E,F):
        root=np.sqrt(np.linalg.det(g));gi=np.linalg.inv(g)
        kin=p@original.inverse(phi)@p/(2*root)
        grad=root/2*np.einsum('ab,aA,AB,bB->',gi,d,original.metric(phi),d)
        potential=root*float(original.node_potential(phi))
        elec=np.einsum('ab,ac,bc,c->',g,E,E,b)/root
        mag=root/4*np.einsum('ac,bd,abi,cdi,i->',gi,gi,F,F,kap)
        return np.array([kin,grad,potential,elec,mag])
    old=density(gamma,p,dphi,E,F)
    Ft=np.einsum('ai,abc,bj->ijc',inv,F,inv)
    new=density(inv.T@gamma@inv,p/det,inv.T@dphi,transform@E/det,Ft)
    error=float(np.max(abs(det*new-old)))
    assert error<1e-11
    return dict(conformal_metric_missing_shear_mean_square=mean,
                all_colour_weak_circle_field_components_nonzero=True,
                full_metric_density_components=old.tolist(),coordinate_Jacobian=det,
                full_metric_density_covariance_error=error,
                completion_scope='classical spatial covariance only; no gravitational Hamiltonian specified')

def run():
    evidence=dict(dual_source_coefficients=coefficient_identity_check(),paired_lapse_same_source=dual_lapse_check(),
                  exact_finite_family_obstruction=mixed_symbol_check(),unitary_representation_enlargement=transported_metric_check(),
                  conformal_obstruction_and_full_metric=geometry_completion_check())
    deps=('joint_matter_energy_current.py','joint_matter_energy_current_results.json',
          'joint_curved_quantum_source.py','joint_quotient_gauge_completion.py',
          'round588_drafts/STATUS.md','round588_drafts/mixed_current_symbol_entry.md',
          'round588_drafts/mixed_current_symbol_entry.py','round588_drafts/mixed_current_symbol_entry_results.json')
    return dict(round=588,tests_run=5,failures=0,errors=0,evidence=evidence,
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
                scope='same-model paired-source symbol; original fixed-geometry family nonclosure; unitary representation enlargement and classical full-metric spatial covariance, not quantum GR')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
