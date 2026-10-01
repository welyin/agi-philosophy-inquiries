"""617: original compact-group cutting with the same energy and sources.

Finite Peter-Weyl fixtures test the analytic full-Haar matching map. They
are not independent physical regions or a cutoff continuum simulation.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_full_spatial_metric as geometry
from joint_chiral_fibre_source import annihilators

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_region_energy_gluing_results.json'
gauge=matter.gauge
original=matter.original

def sample(rng):
    return (gauge.group_exp(rng.normal(size=8),3),
            gauge.group_exp(rng.normal(size=3),2),np.exp(1j*rng.normal()))

def product(A,B):
    return A[0]@B[0],A[1]@B[1],A[2]*B[2]

def inverse(A):
    return A[0].conj().T,A[1].conj().T,1/A[2]

def rep(A,name):
    R=matter.representation(*A)
    block=R[matter.SLICES[name],matter.SLICES[name]]
    return block[::2,::2]  # one fixed spin; gauge multiplicity unchanged

def cut_matrix(d):
    J=np.zeros((d**4,d*d),complex)
    for i,j,k in itertools.product(range(d),repeat=3):
        J[np.ravel_multi_index((i,k,k,j),(d,)*4),i*d+j]=1/np.sqrt(d)
    return J

def axis_apply(A,X,axis,dims):
    tensor=X.reshape(tuple(dims)+(X.shape[-1],))
    moved=np.moveaxis(tensor,axis,0)
    out=np.tensordot(A,moved,axes=(1,0))
    return np.moveaxis(out,0,axis).reshape(X.shape)

def tensor_product(items):
    out=np.ones((1,1))
    for item in items:out=np.kron(out,item)
    return out

def matching_check():
    rng=np.random.default_rng(617)
    center=(np.exp(2j*np.pi/3)*np.eye(3),-np.eye(2),np.exp(1j*np.pi/3))
    rows=[];quotient_error=0.;product_error=0.;transport_error=0.;cut_error=0.
    T=gauge.generators(2)
    for name in matter.SLICES:
        d=rep(center,name).shape[0];J=cut_matrix(d)
        iso=float(np.max(abs(J.conj().T@J-np.eye(d*d))))
        quotient_error=max(quotient_error,float(np.max(abs(rep(center,name)-np.eye(d)))))
        rows.append(dict(original_module=name,dimension=d,isometry_error=iso))
    for _ in range(9):
        A,B,h=sample(rng),sample(rng),sample(rng)
        Ah=product(A,inverse(h));hB=product(h,B)
        for name in matter.SLICES:
            RA,RB=rep(A,name),rep(B,name)
            product_error=max(product_error,float(np.max(abs(rep(product(A,B),name)-RA@RB))))
            cut_error=max(cut_error,float(np.max(abs(rep(Ah,name)@rep(hB,name)-RA@RB))))
        # Original weak doublet, with its actual charge -3.
        RA,RB=rep(A,'L'),rep(B,'L')
        O=2*np.einsum('aij,jk,bkl,li->ab',T,A[1],T,A[1].conj().T).real
        for a in range(3):
            transported=sum(O[a,b]*RA@T[b]@RB for b in range(3))
            transport_error=max(transport_error,float(np.max(abs(transported-T[a]@RA@RB))))
    d=2;J=cut_matrix(d)
    first=np.kron(J,np.eye(d*d))@J
    second=np.kron(np.eye(d*d),J)@J
    associativity=float(np.max(abs(first-second)))
    # Reference preservation at the coefficient level, including entanglement.
    x=rng.normal(size=(4,3))+1j*rng.normal(size=(4,3));x/=np.linalg.norm(x)
    y=J@x
    ref_error=float(np.max(abs(y.conj().T@y-x.conj().T@x)))
    assert max(quotient_error,product_error,transport_error,cut_error,
               associativity,ref_error,max(r['isometry_error'] for r in rows))<2e-13
    return dict(modules=rows,quotient_error=quotient_error,product_error=product_error,
        cut_gauge_error=cut_error,transported_momentum_error=transport_error,
        three_segment_associativity_error=associativity,reference_error=ref_error,
        full_Haar_onto_and_domain_claims_are_analytic=True)

def electric_check():
    # Three original weak links out of one source node. Each PW right index j=0.
    # This fixture tests coefficients before endpoint Gauss; full Gauss is analytic.
    T=gauge.generators(2);I=np.eye(2)
    J=tensor_product([cut_matrix(2)[:,[0,2]]]*3)  # 4096 x 8
    P=[];D=[]
    for mu in range(3):
        pm=[];dm=[]
        for a in range(3):
            pm.append(tensor_product([T[a].T if k==mu else I for k in range(3)]))
            dm.append(axis_apply(T[a].T,J,4*mu,[2]*12))
        P.append(pm);D.append(dm)
    intertwining=max(float(np.max(abs(D[m][a]-J@P[m][a])))
                     for m in range(3) for a in range(3))
    G=np.empty((3,3,8,8),complex);GC=np.empty_like(G)
    for m,n in itertools.product(range(3),repeat=2):
        G[m,n]=sum(P[m][a].conj().T@P[n][a] for a in range(3))
        GC[m,n]=sum(D[m][a].conj().T@D[n][a] for a in range(3))
    S=np.array([[.24,.31,-.12],[.31,-.07,.17],[-.12,.17,-.17]])
    L=np.array([[.3,.06,0],[.06,.6,.04],[0,.04,.8]])
    eps=.73;b=float(geometry.PAR['b'][1])
    def matrices(t,sigma):
        shape=geometry.shape_exp(S,t)
        K=b/eps*np.exp(-2*sigma)*shape
        # A non-proportional positive split. K=KA+KB at every geometry.
        vals,U=np.linalg.eigh(K);root=(U*np.sqrt(vals))@U.T
        KA=root@L@root;KB=K-KA
        h=np.einsum('mn,mnij->ij',K,G)
        hc=np.einsum('mn,mnij->ij',KA,GC)+np.einsum('mn,mnij->ij',KB,GC)
        return K,KA,KB,h,hc
    t=.8;sigma=.12;K,KA,KB,H,HC=matrices(t,sigma)
    step=2e-6
    ds=-2*H
    dt=np.einsum('mn,mnij->ij',K@S,G)
    def fd(index,cut):
        args=[t,sigma];args[index]+=step
        plus=matrices(*args)[4 if cut else 3]
        args[index]-=2*step
        minus=matrices(*args)[4 if cut else 3]
        return (plus-minus)/(2*step)
    source_errors=[float(np.max(abs(fd(0,True)-dt))),float(np.max(abs(fd(1,True)-ds)))]
    diagonal=np.einsum('mn,mnij->ij',np.diag(np.diag(K)),G)
    copied=np.einsum('mn,mnij->ij',2*K,GC)
    energy_error=float(np.max(abs(H-HC)))
    omission=float(np.linalg.norm(H-diagonal,2))
    duplicate=float(np.linalg.norm(copied-H,2))
    source_duplicate=float(np.linalg.norm(-2*copied-ds,2))
    beta=.9;e,U=np.linalg.eigh(H);heat=(U*np.exp(-beta*e))@U.conj().T
    shear_source=float((np.trace(heat@dt)/np.trace(heat)).real)
    def free(tt):
        return -np.log(np.exp(-beta*np.linalg.eigvalsh(matrices(tt,sigma)[4])).sum())/beta
    thermal_fd=(free(t+step)-free(t-step))/(2*step)
    assert min(np.linalg.eigvalsh(KA))>0 and min(np.linalg.eigvalsh(KB))>0
    assert max(intertwining,energy_error,max(source_errors),abs(shear_source-thermal_fd))<2e-8
    assert omission>1e-3 and duplicate>1e-2 and source_duplicate>1e-2
    return dict(coefficient_matrix=K.tolist(),split_A=KA.tolist(),split_B=KB.tolist(),
        frozen_weak_electric_coefficient=b,epsilon=eps,shape_t=t,log_conformal_factor=sigma,
        cut_isometry_shape=list(J.shape),momentum_intertwining_error=intertwining,
        energy_error=energy_error,shear_and_conformal_source_errors=source_errors,
        discarded_shear_operator_defect=omission,duplicated_energy_operator_defect=duplicate,
        duplicated_conformal_source_defect=source_duplicate,
        finite_sector_ordinary_thermal_shear_source=shear_source,thermal_source_fd_error=float(abs(shear_source-thermal_fd)),
        actual_energy_eigenvalues=e.tolist(),fixture_before_endpoint_Gauss=True,
        full_matched_Gauss_thermal_equivalence_is_analytic=True)

def interaction_check():
    rng=np.random.default_rng(6173);c=annihilators(4)
    parity=np.diag([(-1)**int(i).bit_count() for i in range(16)])
    errors=[];hop_norms=[];gradient_values=[];magnetic_values=[]
    def hop(R):
        h=sum(R[i,j]*c[i].conj().T@c[2+j] for i,j in itertools.product(range(2),repeat=2))
        return h+h.conj().T
    for _ in range(6):
        A,B,h=sample(rng),sample(rng),sample(rng)
        U=product(A,B);Ah=product(A,inverse(h));hB=product(h,B)
        single=hop(rep(U,'L'));split=hop(rep(A,'L')@rep(B,'L'))
        recut=hop(rep(Ah,'L')@rep(hB,'L'))
        errors.extend([np.max(abs(single-split)),np.max(abs(single-recut)),
                       np.max(abs(single@parity-parity@single))])
        hop_norms.append(float(np.linalg.norm(single,2)))
        phi=rng.normal(size=5)*.18;chi=rng.normal(size=5)*.18
        def transformed(R):
            X=R@(chi[:2]+1j*chi[2:4])
            return np.r_[X.real,X.imag,chi[4]]
        RU=U[2]**3*U[1];RA=A[2]**3*A[1];RB=B[2]**3*B[1]
        d1=float(original.distance_squared(phi,transformed(RU)))
        d2=float(original.distance_squared(phi,transformed(RA@RB)))
        errors.append(abs(d1-d2));gradient_values.append(d1)
        # Actual non-Abelian plaquette with one cut edge and three unchanged edges.
        V,W,X=sample(rng),sample(rng),sample(rng)
        face=product(product(product(U,V),inverse(W)),inverse(X))
        facecut=product(product(product(product(A,B),V),inverse(W)),inverse(X))
        v1=gauge.potential(*face,geometry.PAR);v2=gauge.potential(*facecut,geometry.PAR)
        errors.append(abs(v1-v2));magnetic_values.append(v1)
    # Noncommuting factors witness that dropping a half-link changes matter.
    A,B=sample(rng),sample(rng)
    defect=float(np.linalg.norm(hop(rep(A,'L'))-hop(rep(product(A,B),'L')),2))
    assert max(errors)<1e-12 and defect>1e-2
    return dict(max_multiplication_and_even_CAR_error=float(max(errors)),
        original_hopping_norms=hop_norms,Higgs_distance_squared_values=gradient_values,
        original_magnetic_potential_values=magnetic_values,
        omitted_half_link_hopping_defect=defect,
        Fock_modes_unchanged=True,new_boundary_matter_modes=0,
        whole_unbounded_scalar_and_fermion_domains_transported_analytically=True)

def run():
    deps=('research_note_360.md','research_note_362.md','research_note_365.md',
          'research_note_574.md','research_note_589.md','research_note_598.md',
          'research_note_603.md','research_note_608.md','research_note_616.md',
          'joint_quotient_gauge_completion.py','joint_curved_quantum_source.py',
          'joint_full_spatial_metric.py','round589_drafts/full_metric_magnetic_entry.py',
          'joint_fermion_gauss_completion.py','joint_chiral_fibre_source.py')
    return dict(round=617,tests_run=3,failures=0,errors=0,
        matching=matching_check(),electric=electric_check(),interaction=interaction_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_quotient_group_and_geometry_coefficients_retained=True,
            matching_unitarily_preserves_full_existing_process_and_sources=True,
            representation_cut_not_new_physical_refinement=True,
            independent_region_preparability_not_proved=True,
            chiral_continuum_and_GR_still_open=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))

