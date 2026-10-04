"""758: an initial Gauss matrix bridge, not a physical-time limit.

The full fixed-graph proof is in research_note_758.md. The calculations below
check its original color stabilizer, finite CAR moments and a curved radial
matrix packet. The radial calculation is not a full graph simulation.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_correlated_hadamard_noise as old
import joint_matter_ground_source as bands

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_matrix_background_results.json'

def unitary(a):
    e,v=np.linalg.eigh(a)
    return (v*np.exp(1j*e))@v.conj().T

def color_stabilizer_check():
    col=old.background.color(1.)
    A,E,T=col['A'],col['E'],old.background.T
    rows=[]
    for N in (4,7,12):
        eps=2*np.pi/N
        links=[unitary(eps*a) for a in A]
        loops=[np.linalg.matrix_power(u,N) for u in links]
        loop_error=max(float(np.max(abs(w-unitary(2*np.pi*a)))) for w,a in zip(loops,A))
        gauss=float(np.max(abs(sum(e-u.conj().T@e@u for e,u in zip(E,links)))))
        # Joint commutant on su(3), tested using complete noncontractible loops.
        columns=np.array([np.concatenate([(t@w-w@t).ravel() for w in loops]) for t in T]).T
        gram=(columns.conj().T@columns).real
        eig=np.linalg.eigvalsh(gram)
        assert loop_error<3e-14 and gauss<1e-14 and eig[0]>.05
        rows.append(dict(N=N,loop_identity_error=loop_error,color_Gauss_error=gauss,
                         joint_commutant_smallest_gram_eigenvalue=float(eig[0])))
    # Analytically: each loop has three distinct eigenvalues; commuting with
    # T1,T2 gives diag(z,z,w), then T4 forces z=w. SU(3) leaves its center.
    gaps=[]
    for a in A:
        phases=np.exp(2j*np.pi*np.linalg.eigvalsh(a))
        gaps.append(float(min(abs(phases[i]-phases[j]) for i in range(3) for j in range(i))))
    assert min(gaps)>.1
    return dict(rows=rows,loop_eigenvalue_minimum_separations=gaps,
                discrete_center_not_removed=True,nonzero_color_momenta_fixed_by_stabilizer=True,
                old_electroweak_Gauss_and_stabilizer_argument_reused=True)

def fiber_check():
    q,psi,*_=old.background.completed(12,1.)
    point=q['phi'][0,3,2];h,d=old.matter_source.matter.mass_matrices(point)
    ix=np.ix_(np.arange(24,32),np.arange(24,32))
    B=old.fock_quadratic(h[ix],d[ix],old.car(8))
    # Full eight-lepton block, original sterile bits 6,7; other six in vacuum.
    vacuum=np.zeros((64,64));vacuum[0,0]=1
    rho=np.kron(old.tau(.5,.5),vacuum)
    ref=np.kron(old.tau(.5,.25),vacuum);difference=rho-ref
    F=float(old.matter_source.matter.original.F(point));Ys=old.matter_source.matter.Y['s']
    second=float(np.trace(difference@B@B).real)
    expected=.5*abs(Ys)**2*point[4]**2/F
    assert abs(second-expected)<2e-14 and abs(np.trace(difference@B))<1e-14
    eig,v=np.linalg.eigh(B);weights=np.diag(v.conj().T@difference@v)
    rows=[]
    for u in (.1,.05,.025):
        gap=np.sum(weights*np.exp(1j*u*eig))
        rows.append(dict(spectral_parameter=u,gap_real=float(gap.real),gap_imaginary=float(gap.imag),
                         gap_over_u_squared=float(gap.real/u**2)))
    assert abs(rows[-1]['gap_over_u_squared']+second/2)<2e-6
    # Both q states are positive with the same complete two-point data.
    assert min(np.linalg.eigvalsh(rho))>=0 and min(np.linalg.eigvalsh(ref))>=0
    return dict(Fock_dimension=256,second_moment_gap=second,analytic_second_gap=expected,
                mean_gap=float(np.trace(difference@B).real),characteristic_rows=rows,
                original_point=point.tolist(),single_point_not_full_graph_spectrum=True,
                spectral_parameter_not_physical_time=True)

def packet(hbar,order=32):
    src,psi,*_=old.background.completed(12,1.)
    idx=(0,3,2);phi0=src['phi'][idx];center=np.array([np.linalg.norm(phi0[:4]),phi0[4]])
    eps=2*np.pi/12;w=float(eps**3*psi[idx]**6)
    momentum=eps**3*src['p'][idx][[1,4]]
    radius=.045
    x,weights=np.polynomial.hermite.hermgauss(order)
    xx=np.stack(np.meshgrid(x,x,indexing='ij'),axis=-1)
    delta=np.sqrt(hbar)*xx;keep=np.all(abs(delta)<radius,axis=-1)
    delta=delta[keep];r=center+delta;h,s=r.T
    weights=(weights[:,None]*weights[None,:])[keep]
    F=2-np.sum(r*r,axis=1)/6
    cutoff=np.exp(-np.sum(delta*delta/(radius**2-delta*delta),axis=1))
    probability=weights*h**3*np.sqrt(2)/F**3*cutoff**2;probability/=probability.sum()
    first=-delta/hbar-2*delta*radius**2/(radius**2-delta*delta)**2+1j*momentum/hbar
    second=-1/hbar-2*radius**2*(radius**2+3*delta*delta)/(radius**2-delta*delta)**3
    metric=F[:,None,None]*(np.eye(2)-r[:,:,None]*r[:,None,:]/12)
    drift=np.column_stack((F*(3/h-h/6),-F*s/6))
    lap=np.einsum('ni,nij,nj->n',first,metric,first)+np.sum(np.diagonal(metric,axis1=1,axis2=2)*second,axis=1)+np.sum(drift*first,axis=1)
    phi=np.column_stack((np.zeros_like(h),h,np.zeros_like(h),np.zeros_like(h),s))
    V=w*bands.old.original.node_potential(phi)
    Hb=-hbar*hbar*lap/(2*w)+V
    K0=old.matter_source.matter.original.inverse(phi0)[np.ix_([1,4],[1,4])]
    hb0=float(momentum@K0@momentum/(2*w)+w*bands.old.original.node_potential(phi0))
    pair,der=bands.radial(*center,bands.IDS)
    B0=bands.response.fock(pair);S0=bands.response.fock(der);I=np.eye(len(B0))
    vac=np.zeros((4,4));vac[0,0]=1
    diff=np.kron(old.tau(.5,.5)-old.tau(.5,.25),vac)
    Hnorm=Snorm=source_gap=0.
    for j,(hj,sj) in enumerate(r):
        pair,der=bands.radial(hj,sj,bands.IDS)
        B=bands.response.fock(pair);S=bands.response.fock(der)
        D=(Hb[j]-hb0)*I+B-B0
        # Frobenius norm is an upper bound on operator norm, without diagonalizing.
        Hnorm+=probability[j]*float(np.vdot(D,D).real)
        Snorm+=probability[j]*float(np.vdot(S-S0,S-S0).real)
        source_gap+=probability[j]*float(np.trace(diff@S@S).real)
    target_source_gap=float(np.trace(diff@S0@S0).real)
    c=center[1];F0=float(old.matter_source.matter.original.F(phi0))
    analytic_source_gap=.5*abs(old.matter_source.matter.Y['s'])**2*(F0**-.5+c*c/(6*F0**1.5))**2
    assert abs(target_source_gap-analytic_source_gap)<2e-14
    # Original bounded two-output field instrument, at the initial instant.
    L=np.sqrt(.5+np.sin(s)/4);L0=np.sqrt(.5+np.sin(c)/4)
    record_error=float(probability@abs(L-L0)**2)
    return dict(hbar=hbar,normalization=float(probability.sum()),
                H_matrix_squared_Frobenius_residual=Hnorm,source_matrix_squared_Frobenius_residual=Snorm,
                retained_source_variance_gap=source_gap,limiting_source_variance_gap=target_source_gap,
                original_initial_record_squared_residual=record_error,
                background_point=center.tolist(),cell_volume=w,canonical_radial_momentum=momentum.tolist())

def curved_matrix_packet_check():
    rows=[packet(h) for h in (1e-4,2.5e-5,6.25e-6)]
    fine=packet(6.25e-6,48)
    keys=('H_matrix_squared_Frobenius_residual','source_matrix_squared_Frobenius_residual',
          'retained_source_variance_gap','original_initial_record_squared_residual')
    quad=max(abs(fine[k]-rows[-1][k]) for k in keys)
    assert quad<2e-10
    for k in keys[:2]:assert rows[-1][k]<rows[0][k]/8
    assert abs(rows[-1]['retained_source_variance_gap']-rows[-1]['limiting_source_variance_gap'])<1e-6
    return dict(rows=rows,quadrature_refinement_error=quad,
                retained_quantum_matrix_not_one_band_or_mean_field=True,
                original_curved_radial_fixture_not_full_graph_or_full_Gauss_simulation=True)

def run():
    results=(color_stabilizer_check(),fiber_check(),curved_matrix_packet_check())
    deps=('research_note_574.md','research_note_643.md','research_note_704.md','research_note_706.md',
          'research_note_726.md','research_note_727.md','research_note_753.md','research_note_756.md',
          'research_note_757.md','joint_reference_constraint_strata.py','joint_matter_ground_source.py',
          'joint_correlated_hadamard_noise.py')
    return dict(round=758,tests_run=3,failures=0,errors=0,color_stabilizer=results[0],
                original_correlated_fiber=results[1],curved_matrix_packet=results[2],
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='An initial fixed-graph Gauss isometry retaining the stabilizer-invariant finite Fock matrix, full principal-symbol source matrices, initial smooth instruments and bounded spectral tests. No gap or one-band replacement. Numerical checks are color loop algebra and local matter/radial fixtures. No finite physical-time, fixed-hbar continuum, dynamic Einstein constraint, autonomous preparation or quantum-gravity theorem.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
