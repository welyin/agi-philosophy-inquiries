"""554: finite-cutoff h,s reference menu, joint records and discarded-mode bounds.

No continuum error bound, autonomous detector, full gauge QFT, or quantum
gravity is claimed. The Gaussian POVM is an explicitly added instrument.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_matter_reference_quantum_readout_results.json'


def geometry():
    xyz=np.array(list(itertools.product((-.5,.5),repeat=3)))
    N=len(xyz); D=[]
    for axis in range(3):
        mat=np.zeros((N,N))
        for i,x in enumerate(xyz):
            y=x.copy();y[axis]*=-1
            j=np.flatnonzero(np.all(xyz==y,axis=1))[0]
            mat[i,j]=1/(y[axis]-x[axis]);mat[i,i]=-mat[i,j]
        assert np.allclose(mat@xyz[:,axis],1)
        assert np.allclose(mat@np.ones(N),0)
        D.append(mat)
    return xyz,D


def menu(N,v,D):
    dim=4*N; K=sum(d.T@d for d in D)
    O=np.block([[np.zeros((2*N,2*N)),np.eye(2*N)],
                [-np.eye(2*N),np.zeros((2*N,2*N))]])
    ls=[np.zeros(dim) for _ in range(4)]
    As=[np.zeros((dim,dim)) for _ in range(4)]
    ls[0][:N]=1/N;ls[1][N:2*N]=1/N
    for field,j in ((0,2),(1,3)):
        q=slice(field*N,(field+1)*N);p=slice((field+2)*N,(field+3)*N)
        As[j][q,q]=K/N;As[j][p,p]=-np.eye(N)/(N*v*v)
    return O,ls,As


def gaussian_moments(m,S,ls,As):
    means=np.array([l@m+m@A@m+np.trace(A@S) for l,A in zip(ls,As)])
    grads=[l+2*A@m for l,A in zip(ls,As)]
    cov=np.array([[grads[i]@S@grads[j]+2*np.trace(As[i]@S@As[j]@S)
                   for j in range(4)] for i in range(4)])
    return means,cov


def independent_variance(m,S,l,A):
    # Rotate Gaussian quadratic to independent normal coordinates; exact
    # 5-node Hermite integration of each degree-four one-dimensional moment.
    root=np.diag(np.sqrt(np.diag(S)))
    eig,U=np.linalg.eigh(root@A@root)
    b=U.T@root@(l+2*A@m)
    nodes,weights=np.polynomial.hermite.hermgauss(5)
    nodes*=math.sqrt(2);weights/=math.sqrt(math.pi)
    value=0.
    for lam,lin in zip(eig,b):
        samples=lam*nodes**2+lin*nodes
        value+=weights@samples**2-(weights@samples)**2
    return float(value)


def energy(m,G,N,v,D,C,L,V0):
    qh,qs,ph,ps=np.split(m,4)
    ah=np.diag(G)[:N];ass=np.diag(G)[N:2*N]
    kinetic=(ph@ph+ps@ps+np.trace(G[2*N:,2*N:]))/(2*v)
    gradient=0.
    for d in D:
        for q,block in ((qh,G[:N,:N]),(qs,G[N:2*N,N:2*N])):
            gradient+=v*(np.linalg.norm(d@q)**2+np.trace(d@block@d.T))/2
    X=qh*qh+ah;Y=qs*qs+ass
    H4=qh**4+6*qh*qh*ah+3*ah*ah
    S4=qs**4+6*qs*qs*ass+3*ass*ass
    potential=v*np.sum(V0-(C[0]*X+C[1]*Y)/2+
                       (L[0,0]*H4+2*L[0,1]*X*Y+L[1,1]*S4)/4)
    return dict(kinetic=float(kinetic),gradient=float(gradient),
                potential=float(potential),total=float(kinetic+gradient+potential))


def run():
    checks=[];hbar=1.;v=1.;xyz,D=geometry();N=len(xyz);V=N*v
    O,ls,As=menu(N,v,D)
    coeff=2*ls[0]@O@As[2]
    wanted=np.zeros(4*N);wanted[2*N:3*N]=-2/(V*V)
    assert np.max(abs(coeff-wanted))<1e-15
    assert np.max(abs(2*ls[1]@O@As[3]-np.roll(wanted,N)))<1e-15
    # The spatial-gradient matrices do not enter either commutator.
    assert np.max(abs(coeff[:2*N]))==0
    checks.append('exact_finite_cell_composite_menu_commutator_including_all_spatial_modes')

    source=HERE/'joint_singlet_common_mass_rg_results.json'
    saved=json.loads(source.read_text('utf8'));row=saved['examples'][2]['state']
    L=np.array([[row['lambda_H'],row['p']],[row['p'],row['lambda_s']]])
    C=.25*np.array([row['x'],row['y']]);u=np.linalg.solve(L,C)
    hs=np.sqrt(u);V0=C@u/4
    omega=32.;a=hbar/(2*v*omega);b=hbar*v*omega/2
    G=np.diag([a]*(2*N)+[b]*(2*N));seed=G.copy();S=G+seed
    assert np.min(np.linalg.eigvalsh(G+1j*hbar*O/2))>-1e-12
    assert np.min(np.linalg.eigvalsh(seed+1j*hbar*O/2))>-1e-12
    cases=[];maxquad=0.;maxjac=0.
    for eps in (.25,1.,4.,16.):
        m=np.r_[np.full(N,hs[0]),hs[1]+eps*(xyz[:,0]+xyz[:,0]*xyz[:,2]),
                v*eps*(1+xyz[:,1]),np.zeros(N)]
        J=np.diag([eps,eps,-2*eps**2,2*eps**2]);inv=np.linalg.inv(J)
        # Four specified Cauchy-data perturbations. These are a finite menu
        # calibration, not a proof of continuum spacetime-coordinate accuracy.
        dm=np.zeros((4*N,4));dm[:N,0]=eps;dm[N:2*N,1]=eps
        dm[2*N:3*N,2]=v*eps;dm[N:2*N,3]=eps*xyz[:,0]
        measuredJ=np.vstack([l+2*A@m for l,A in zip(ls,As)])@dm
        maxjac=max(maxjac,float(np.max(abs(measuredJ-J))))
        assert np.allclose(measuredJ,J,atol=1e-12,rtol=0)
        classical=np.array([l@m+m@A@m for l,A in zip(ls,As)])
        outmean,outcov=gaussian_moments(m,S,ls,As)
        seedbias=np.array([np.trace(A@seed) for A in As])
        meanerror=inv@(outmean-seedbias-classical)
        sourcebias=inv@np.array([np.trace(A@G) for A in As])
        assert np.allclose(meanerror,sourcebias,atol=1e-12,rtol=0)
        recordcov=inv@outcov@inv.T
        MSE=np.diag(recordcov)+meanerror**2
        for i in range(4):
            for j in range(i,4):
                vij=independent_variance(m,S,ls[i]+ls[j],As[i]+As[j])
                vi=independent_variance(m,S,ls[i],As[i])
                vj=independent_variance(m,S,ls[j],As[j])
                maxquad=max(maxquad,abs((vij-vi-vj)/2-outcov[i,j]))
        _,wigner_cov=gaussian_moments(m,G,ls,As)
        quantum_cov=wigner_cov+np.array([[hbar*hbar*np.trace(A@O@B@O)/2
                                         for B in As] for A in As])
        intrinsic=inv@quantum_cov@inv.T
        bound=hbar/(2*V*eps**2)
        assert intrinsic[0,0]*intrinsic[2,2]>=bound**2-1e-12
        assert recordcov[0,0]*recordcov[2,2]>=(2*bound)**2-1e-12
        # A finite four-number record: clip beyond radius, retaining overflow.
        alpha=.05;radius=math.sqrt(float(np.sum(MSE))/alpha);step=radius/32
        limit=radius+step;bins=math.ceil(2*limit/step)+2
        bits=4*math.ceil(math.log2(bins))
        sourceE=energy(m,G,N,v,D,C,L,V0)
        assert sourceE['total']>0 and all(np.isfinite(list(sourceE.values())))
        # Union bound for Q_h crossing the positive radial patch at this instant.
        patch_bound=N*.5*math.erfc(hs[0]/math.sqrt(2*a))
        cases.append(dict(epsilon=eps,V=V,source_covariance_q=a,source_covariance_p=b,
            classical_menu=classical.tolist(),source_bias_after_meter_debiasing=meanerror.tolist(),
            decoded_record_variances=np.diag(recordcov).tolist(),decoded_MSE=MSE.tolist(),
            intrinsic_std_product=float(math.sqrt(intrinsic[0,0]*intrinsic[2,2])),
            intrinsic_Robertson_lower_bound=bound,
            record_std_product=float(math.sqrt(recordcov[0,0]*recordcov[2,2])),
            Gaussian_record_lower_bound=2*bound,
            source_energy=sourceE,ancilla_seed_oscillator_energy=2*N*hbar*omega/2,
            h_negative_event_union_bound=patch_bound,
            finite_record=dict(confidence_at_anchor_at_least=1-alpha,radius=radius,
                quantization_step=step,decoded_error_at_most=radius+step/2,
                bins_per_channel_including_overflow=bins,total_bits=bits)))
    assert maxquad<1e-7
    checks.append('same_frozen_matter_parameters_and_four_reference_calibration_directions')
    checks.append('positive_Gaussian_joint_POVM_full_quadratic_moments_and_independent_quadrature')
    checks.append('source_bias_intrinsic_vs_record_uncertainty_and_finite_record_resource_accounting')

    # Exact missing-modes identity: averages of squares are not square averages.
    rng=np.random.default_rng(554)
    P=rng.normal(size=(40,N));Ptot=P.sum(axis=1)
    perpendicular=P-P.mean(axis=1,keepdims=True)
    Korth=np.sum(perpendicular**2,axis=1)/(2*v)
    defect=np.mean((P/v)**2,axis=1)-(Ptot/V)**2
    assert np.max(abs(defect-2*Korth/V))<1e-14
    checks.append('omitted_momentum_modes_exact_energy_defect_not_a_smearing_identity')

    # Mixtures of coherent momentum displacements in ONE orthogonal mode.
    # Same uniform-mode reduced state, same q distribution and full H mean energy,
    # but no common fourth-moment / squared-readout risk bound.
    basevar=Qvar=1/2;extra=1.;tails=[]
    for amplitude in (2.,4.,16.,64.):
        weight=extra/amplitude**2
        second=basevar+extra
        fourth=3*basevar**2+6*basevar*extra+extra*amplitude**2
        variance=fourth-second**2
        nodes,weights=np.polynomial.hermite.hermgauss(5)
        normal=math.sqrt(2*basevar)*nodes;weights=weights/math.sqrt(math.pi)
        q2=q4=0.
        for prob,shift in ((1-weight,0.),(weight/2,amplitude),(weight/2,-amplitude)):
            q2+=prob*(weights@(normal+shift)**2)
            q4+=prob*(weights@(normal+shift)**4)
        assert abs(q2-second)<1e-12 and abs(q4-fourth)<1e-9
        tails.append(dict(momentum_displacement=amplitude,mixture_weight=weight,
            orthogonal_mode_p_second=second,p_fourth=fourth,
            mean_kinetic_energy=second/(2*v),
            variance_of_missing_A=variance/(N*N*v**4)))
    assert len({row['mean_kinetic_energy'] for row in tails})==1
    assert all(tails[j+1]['variance_of_missing_A']>tails[j]['variance_of_missing_A'] for j in range(3))
    checks.append('fixed_energy_same_retained_state_unbounded_omitted_quadratic_readout_variance')
    deps=('research_note_548.md','research_note_549.md','research_note_553.md',
          'joint_singlet_common_mass_rg_results.json')
    return dict(round=554,tests_run=len(checks),failures=0,errors=0,checks=checks,
        N=N,v=v,finite_cutoff_coordinates=xyz.tolist(),
        commutator_coefficient_of_total_P=-2/(V*V),
        calibration_residual=maxjac,independent_quadrature_covariance_residual=maxquad,
        same_matter_parameters=dict(L=L.tolist(),C=C.tolist(),V0=float(V0),vevs=hs.tolist()),
        examples=cases,fixed_energy_omitted_mode_counterfamily=tails,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope=dict(finite_canonical_cutoff_of_reduced_two_scalar_sector=True,
            full_Higgs_gauge_quantization_not_proved=True,Gaussian_POVM_is_extra_instrument_input=True,
            source_and_meter_initially_independent=True,no_autonomous_control_claim=True,
            calibration_parameters_not_proved_spacetime_coordinates=True,
            no_controlled_continuum_or_truncation_error=True,no_quantum_gravity_or_dimension_generation=True,
            entire_unification_completed=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors',
        'calibration_residual','independent_quadrature_covariance_residual')}))
