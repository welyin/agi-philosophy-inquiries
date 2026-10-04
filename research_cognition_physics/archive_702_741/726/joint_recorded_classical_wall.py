"""726: actual position records and the original classical-wall task.

The radial quadrature reuses574's compact packet recipe with actual651
canonical normalization. Neighbor fields are fixed for this numerical
conditional fixture; the full-graph Gauss argument is analytic.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_material_boundary_transport as wall
import joint_curved_quantum_source as original
import joint_fermion_gauss_completion as matter

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_recorded_classical_wall_results.json'
M=original.M


def source():
    N=32;d=wall.old.fields(N);c,lam,_=wall.base_data();idx=(0,N//4,N//8)
    eps=2*np.pi/N;psi=float(d['psi'][idx]);w=eps**3*psi**6
    center=np.array([c['H'],c['S']]);p=eps**3*np.array([0.,c['p']])
    field=d['q']['phi'][...,4]-lam*np.linalg.norm(d['q']['phi'][...,:4],axis=-1)
    neighbor=np.array([np.roll(field,-1,axis=j)[idx] for j in range(3)])
    return dict(N=N,c=c,lam=lam,a=-lam,eps=eps,psi=psi,w=w,center=center,p=p,neighbor=neighbor)


def packet(hbar,order=60):
    src=source();center=src['center'];p=src['p'];w=src['w'];a=src['a']
    R=.18;x,weights=np.polynomial.hermite.hermgauss(order)
    xx=np.stack(np.meshgrid(x,x,indexing='ij'),axis=-1);delta=np.sqrt(hbar)*xx
    keep=np.all(abs(delta)<R,axis=-1);q=(center+delta)[keep];delta=delta[keep]
    weights=(weights[:,None]*weights[None,:])[keep]
    h,s=q.T;F=M-np.sum(q*q,axis=1)/6
    chi=np.exp(-np.sum(delta*delta/(R*R-delta*delta),axis=1))
    probability=weights*h**3*np.sqrt(M)/F**3*chi**2;probability/=sum(probability)
    l1=-delta/hbar-2*delta*R**2/(R**2-delta*delta)**2+1j*p/hbar
    l2=-1/hbar-2*R**2*(R**2+3*delta*delta)/(R**2-delta*delta)**3
    K=F[:,None,None]*(np.eye(2)-q[:,:,None]*q[:,None,:]/(6*M))
    v=np.array([a,1.]);f=q@v;X=np.einsum('nij,j->ni',K,v)
    Gamma=X@v;lap=F*(-f/(3*M)+3*a/h)
    dX=np.empty((len(q),2,2))
    for j in range(2):
        for k in range(2):
            dX[:,j,k]=-q[:,k]/3*(v[j]-q[:,j]*f/(6*M))-F*((j==k)*f+q[:,j]*v[k])/(6*M)
    accel=np.einsum('njk,nk->nj',dX,X)
    dlap=-q/3*(-f/(3*M)+3*a/h)[:,None]+F[:,None]*np.column_stack((
        np.full(len(q),-a/(3*M))-3*a/h**2,np.full(len(q),-1/(3*M))))
    r=np.sum(X*l1,axis=1)+lap/2
    curv=np.sum(X*X*l2,axis=1)+np.sum(accel*l1,axis=1)+np.sum(X*dlap,axis=1)/2
    L0=-hbar*hbar/w**2*(r*r+curv)
    spatial=src['psi']**-4*np.sum((src['neighbor'][None,:]-f[:,None])**2,axis=1)/src['eps']**2
    Q0=spatial-L0
    drift=np.column_stack((F*(3/h-h/(3*M)),-F*s/(3*M)))
    lapratio=np.einsum('ni,nij,nj->n',l1,K,l1)+np.sum(np.diagonal(K,axis1=1,axis2=2)*l2,axis=1)+np.sum(drift*l1,axis=1)
    phi=np.column_stack((np.zeros_like(h),h,np.zeros_like(h),np.zeros_like(h),s))
    potential=w*original.node_potential(phi);T0=-hbar*hbar/(2*w)*lapratio
    H0=T0+potential;G0=-6*T0+6*potential
    fstar=float(center@v);Kstar=original.inverse(np.array([0.,center[0],0.,0.,center[1]]))[np.ix_([1,4],[1,4])]
    vstar=float(v@Kstar@p/w)
    Wstar=src['psi']**-4*sum((src['neighbor']-fstar)**2)/src['eps']**2
    Qstar=Wstar-vstar*vstar
    phistar=np.array([0.,center[0],0.,0.,center[1]])
    Tstar=float(p@Kstar@p/(2*w));Pstar=float(w*original.node_potential(phistar))
    Hstar=Tstar+Pstar;Gstar=-6*Tstar+6*Pstar
    sigma=hbar**.4
    # Polynomial in t=(y-f)/sigma. Differentiation is at fixed y first.
    r1=Gamma/(2*sigma)
    Lcoef=np.stack((
        L0+hbar*hbar/w**2*Gamma**2/(2*sigma*sigma),
        -hbar*hbar/w**2*(2*r*r1+(accel@v)/(2*sigma)),
        -hbar*hbar/w**2*r1*r1),axis=1)
    Qcoef=-Lcoef;Qcoef[:,0]+=spatial
    Tcoef=np.stack((T0+hbar*hbar*Gamma/(4*w*sigma*sigma),
                   -hbar*hbar*r/(2*w*sigma),-hbar*hbar*Gamma/(8*w*sigma*sigma)),axis=1)
    Hcoef=Tcoef.copy();Hcoef[:,0]+=potential
    Gcoef=-6*Tcoef;Gcoef[:,0]+=6*potential
    gram=np.array([[1.,0.,1.],[0.,1.,0.],[1.,0.,3.]])
    def stats(coef,star):
        residual=coef.copy();residual[:,0]-=star
        mse=float(probability@np.einsum('ni,ij,nj->n',residual.conj(),gram,residual).real)
        mean=probability@(coef[:,0]+coef[:,2])
        return dict(mean=float(mean.real),classical_value=float(star),mean_error=float(abs(mean.real-star)),
                    mean_imaginary=float(abs(mean.imag)),record_weighted_square_residual=mse)
    Qstats=stats(Qcoef,Qstar);Hstats=stats(Hcoef,Hstar);Gstats=stats(Gcoef,Gstar)
    position=float(probability@(f-fstar)**2+sigma*sigma)
    deltaQ=float((probability@(Qcoef[:,0]+Qcoef[:,2]-Q0)).real)
    formulaQ=float(-hbar*hbar/(4*w*w*sigma*sigma)*(probability@(Gamma*Gamma)))
    assert abs(deltaQ-formulaQ)<2e-12
    # Independent Gaussian integration of the full polynomials, not the moment table.
    xx,ww=np.polynomial.hermite.hermgauss(12);t=np.sqrt(2)*xx;ww/=np.sqrt(np.pi)
    direct=0.
    for z,weight in zip(t,ww):
        direct+=weight*float(probability@abs(Qcoef[:,0]+z*Qcoef[:,1]+z*z*Qcoef[:,2]-Qstar)**2)
    assert abs(direct-Qstats['record_weighted_square_residual'])<1e-12
    # Raw zero-radius convention is not used; packets stay in h>0,F>0.
    mf2=float(probability@(abs(matter.Y['s'])**2*s*s/F))
    return dict(hbar=hbar,sigma=sigma,position_mean_square_error=position,
                Q=Qstats,H_bosonic_local=Hstats,G_bosonic_local=Gstats,
                exact_Q_backaction=deltaQ,predicted_Q_backaction=formulaQ,
                independent_gaussian_residual_error=float(abs(direct-Qstats['record_weighted_square_residual'])),
                Q_wrong_sign_spectral_bound=min(1.,Qstats['record_weighted_square_residual']/Qstar**2),
                retained_Majorana_square=mf2,
                full_H_local_square_residual=Hstats['record_weighted_square_residual']+mf2)


def recorded_packet_check():
    rows=[packet(hbar) for hbar in (1e-10,2.5e-11,6.25e-12)]
    repeat=packet(2.5e-11,80);mid=rows[1]
    error=max(abs(repeat[k]['record_weighted_square_residual']-mid[k]['record_weighted_square_residual'])
              for k in ('Q','H_bosonic_local','G_bosonic_local'))
    assert error<1e-10
    assert rows[-1]['position_mean_square_error']<rows[0]['position_mean_square_error']/5
    for key in ('Q','H_bosonic_local','G_bosonic_local'):
        assert rows[-1][key]['record_weighted_square_residual']<rows[0][key]['record_weighted_square_residual']/4
    return dict(rows=rows,quadrature_order_check=error,
                original574_compact_packet_recipe=True,actual_K_y_poststate_not_hidden_phase_labels=True,
                conditional_radial_fixture_not_full_graph_wavefunction=True)


def prescribed_wall_check():
    s=source();b=s['c']['b'];u=.05*s['c']['H']/1.05;a=s['a'];psi=s['psi']
    target=psi**-4*b*b;rows=[]
    for eps in (2*np.pi/32,2*np.pi/64,2*np.pi/128):
        t=a*u*(1-np.cos(eps))/eps
        value=psi**-4*(t*t+(b*np.sin(eps)/eps+t)**2)
        bound=psi**-4*(b*a*u*eps+(a*a*u*u/2+b*b/3)*eps*eps)
        assert value>0 and abs(value-target)<=bound+1e-13
        rows.append(dict(epsilon=float(eps),normal_symbol=float(value),continuous_normal=target,
                         actual_error=float(abs(value-target)),analytic_error_bound=float(bound)))
    f0=float(s['center']@np.array([a,1.]))
    direct=s['psi']**-4*sum((s['neighbor']-f0)**2)/s['eps']**2
    assert abs(direct-rows[0]['normal_symbol'])<1e-13
    assert rows[-1]['actual_error']<rows[0]['actual_error']/4
    return dict(rows=rows,original32_grid_formula_error=float(abs(direct-rows[0]['normal_symbol'])),
                same_positive_geometry_at_original_point=True,
                not_a_quantum_grid_refinement_or_uniform_packet_estimate=True)


def retained_fermion_check():
    s=source();h,z=s['center'];phi=np.array([0.,h,0.,0.,z])
    hd,delta=matter.mass_matrices(phi)
    # Full32-mode vacuum action: number-conserving part annihilates, pair norm.
    norm2=float(np.sum(abs(delta)**2)/2);expected=abs(matter.Y['s']*z/np.sqrt(original.F(phi)))**2
    # Independent exact two-mode CAR restriction for the only nonzero vacuum output.
    c0=np.array([[0,1],[0,0]],complex);parity=np.diag([1,-1])
    c1=np.kron(c0,np.eye(2));c2=np.kron(parity,c0)
    m=delta[30,31];pair=m*c1.conj().T@c2.conj().T
    B=pair+pair.conj().T;vac=np.eye(4)[:,0]
    direct=float(np.linalg.norm(B@vac)**2)
    assert abs(norm2-expected)<1e-14 and abs(direct-expected)<1e-14 and expected>0
    return dict(original_majorana_square=float(expected),full32_pair_norm=norm2,
                independent_CAR_norm=direct,original_Dirac_block_norm=float(np.linalg.norm(hd)),
                empty_reference_not_ground_state=True,
                same_reading_keeps_vacuum_Fock_factor_but_not_full_H_stationarity=True,
                no_global_no_semiclassical_limit_claim=True)


def run():
    results=dict(actual_recorded_classical_packet=recorded_packet_check(),
                 original_wall_spatial_mapping=prescribed_wall_check(),
                 retained_fermion_fluctuation=retained_fermion_check())
    deps=('research_note_574.md','research_note_579.md','research_note_651.md','research_note_719.md',
          'research_note_725.md','joint_curved_quantum_source.py','joint_fermion_gauss_completion.py',
          'joint_material_boundary_transport.py','round726_drafts/common_wall_limit_entry.py',
          'round726_drafts/common_wall_limit_entry_results.json')
    return dict(round=726,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Same original574 Gauss packets and724 Gaussian region instrument: fixed-graph recorded positions and bosonic wall-normal/energy/source menu jointly concentrate when sigma->0 and hbar/sigma->0. Original651 forward-difference wall retained. Original598 empty CAR reference keeps zero first mass mean but a strictly positive Majorana full-energy variance floor; no full-source fluctuation-free limit, finite-hbar quantum continuum, autonomous detector or Einstein gluing claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');a=parser.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
