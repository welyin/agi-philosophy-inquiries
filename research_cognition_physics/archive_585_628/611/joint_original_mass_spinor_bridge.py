"""611: original 598 masses in the explicit 606 zero-momentum spinor fibres.

Rest-frame spinor dictionary is input. Off-background results are mass-block
changes, not a computation of physical pole masses or a chiral Hamiltonian.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as original
import joint_chiral_fibre_source as chiral
import joint_gapped_link_locality as link

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_original_mass_spinor_bridge_results.json'
CHARGES=dict(Q=1,u=4,d=-2,L=-3,e=-6,nu=0)
LEFT={'Q','L'}
VL=np.vstack([np.zeros((2,2)),np.eye(2)]).astype(complex)
VR=chiral.GAMMA[3]@VL
BETA=chiral.GAMMA[3]
EPS=np.array([[0,1],[-1,0]],complex)
EPS4=np.kron(np.eye(2),EPS)
PHI=np.array([.4,.7,-.2,.1,.5])
def op(A):return float(np.linalg.norm(A,2))
def dagger(A):return A.conj().T
def module_for_mode(i):
    return next(name for name,s in original.SLICES.items() if s.start<=i<s.stop)
def embed(theta=0.,mu=0):
    J=np.zeros((64,32),complex);dJ=np.zeros_like(J)
    for i in range(16):
        name=module_for_mode(2*i);q=CHARGES[name]
        V=VL if name in LEFT else VR
        U=np.cos(q*theta/2)*np.eye(4)+1j*np.sin(q*theta/2)*chiral.GAMMA[mu]
        J[4*i:4*i+4,2*i:2*i+2]=U@V
        dJ[4*i:4*i+4,2*i:2*i+2]=1j*q*chiral.GAMMA[mu]@U@V/2
    return J,dJ
def ambient(phi,spin=BETA):
    h,d=original.mass_matrices(phi)
    normal=np.zeros((64,64),complex);pair=np.zeros_like(normal)
    for i in range(16):
        for j in range(16):
            block=h[2*i:2*i+2,2*j:2*j+2]
            scalar=np.trace(block)/2
            assert op(block-scalar*np.eye(2))<1e-12
            normal[4*i:4*i+4,4*j:4*j+4]=scalar*spin
    nu=original.SLICES['nu'].start//2
    pair[4*nu:4*nu+4,4*nu:4*nu+4]=d[-2,-1]*EPS4
    return normal,pair
def bdg(h,d):
    return np.block([[h,d],[dagger(d),-h.T]])
def spin_rep(R):
    out=np.zeros((64,64),complex)
    for i in range(16):
        for j in range(16):
            block=R[2*i:2*i+2,2*j:2*j+2]
            scalar=np.trace(block)/2
            assert op(block-scalar*np.eye(2))<1e-12
            out[4*i:4*i+4,4*j:4*j+4]=scalar*np.eye(4)
    return out
def rest_dictionary_check():
    J,_=embed();p=J@dagger(J);q=np.eye(64)-p
    h,d=original.mass_matrices(PHI);b,delta=ambient(PHI)
    naive,_=ambient(PHI,np.eye(4))
    compressed=dagger(J)@b@J;paired=dagger(J)@delta@J.conj()
    normal_error=op(compressed-h);pair_error=op(paired-d)
    assert max(normal_error,pair_error,op(dagger(J)@J-np.eye(32)))<1e-12
    assert op(dagger(J)@naive@J)<1e-12 and op(h)>.5
    assert op(p@b@q)<1e-12
    qq=op(q@delta@q.T);assert qq>.1
    chosen=p@delta@p.T
    assert op(dagger(J)@chosen@J.conj()-d)<1e-12
    source=bdg(compressed,paired);target=bdg(h,d)
    assert op(source-target)<1e-12
    rng=np.random.default_rng(611);covariance=0.
    for _ in range(4):
        C=original.gauge.group_exp(rng.normal(size=8),3)
        W=original.gauge.group_exp(rng.normal(size=3),2);z=np.exp(.3j*rng.normal())
        R=original.representation(C,W,z);RR=spin_rep(R)
        X=z**3*W@(PHI[:2]+1j*PHI[2:4])
        phi=np.r_[X.real,X.imag,PHI[4]]
        bt,dt=ambient(phi)
        covariance=max(covariance,op(RR@J-J@R),op(bt-RR@b@dagger(RR)),
                       op(dt-RR@delta@RR.T))
    assert covariance<1e-12
    return dict(original_normal_mass_norm=op(h),naive_identity_compressed_norm=op(dagger(J)@naive@J),
        normal_mass_error=normal_error,Majorana_error=pair_error,
        BdG_original_spectrum_error=float(np.max(abs(np.linalg.eigvalsh(source)-np.linalg.eigvalsh(target)))),
        complement_pair_block_removed_norm=qq,gauge_dictionary_error=covariance,
        no_Fock_dimension_counted_as_species=True)

def expected(phi,theta,mu):
    h,d=original.mass_matrices(phi);out=h.copy();derivative=np.zeros_like(h)
    pairs=[('Q','u'),('Q','d'),('L','nu'),('L','e')]
    for left,right in pairs:
        charge=CHARGES[left]+CHARGES[right] if mu<3 else CHARGES[right]-CHARGES[left]
        factor=np.cos(charge*theta/2);df=-charge*np.sin(charge*theta/2)/2
        a=original.SLICES[left];b=original.SLICES[right]
        out[a,b]=factor*h[a,b];out[b,a]=factor*h[b,a]
        derivative[a,b]=df*h[a,b];derivative[b,a]=df*h[b,a]
    return out,d,derivative
def varying_fibre_check():
    rows=[];actual=0.
    b,delta=ambient(PHI)
    for mu in (0,3):
        for theta in (0.,.08,.2):
            J,_=embed(theta,mu);h=dagger(J)@b@J;d=dagger(J)@delta@J.conj()
            expected_h,expected_d,_=expected(PHI,theta,mu)
            assert max(op(h-expected_h),op(d-expected_d))<1e-12
            old_h,old_d=original.mass_matrices(PHI)
            shift=float(np.max(abs(np.linalg.eigvalsh(bdg(h,d))-np.linalg.eigvalsh(bdg(old_h,old_d)))))
            factors={}
            for left,right in [('Q','u'),('Q','d'),('L','nu'),('L','e')]:
                q=CHARGES[left]+CHARGES[right] if mu<3 else CHARGES[right]-CHARGES[left]
                factors[right]=float(np.cos(q*theta/2))
            rows.append(dict(direction=mu,theta=theta,Dirac_factors=factors,
                compressed_BdG_mass_block_spectrum_change=shift,
                Majorana_change=op(d-old_d)))
        # Extract full-kernel zero momentum projectors, all original charges.
        theta=.08;F=np.kron(np.ones((16,1))/4,np.eye(4))
        for name,charge in CHARGES.items():
            phases={(x,mu):theta for x in chiral.SITES}
            H,_,_=link.kernel(2,phases,charge=charge)
            E,V=np.linalg.eigh(H);select=(E>0) if name in LEFT else (E<0)
            p=(V*select)@dagger(V)
            U=np.cos(charge*theta/2)*np.eye(4)+1j*np.sin(charge*theta/2)*chiral.GAMMA[mu]
            V0=VL if name in LEFT else VR
            actual=max(actual,op(dagger(F)@p@F-U@V0@dagger(V0)@dagger(U)))
    assert actual<1e-12
    return dict(rows=rows,actual_full_Wilson_projection_error=actual,
        original_Y_and_scalar_background_not_refitted=True,
        compressed_mass_block_not_physical_pole_mass=True,
        time_direction_choice_is_input=True)

def source_and_adjustment_check():
    theta=.17;mu=0;eps=1e-6
    b,delta=ambient(PHI);J,dJ=embed(theta,mu)
    h=dagger(J)@b@J;dh=dagger(dJ)@b@J+dagger(J)@b@dJ
    analytic=expected(PHI,theta,mu)[2]
    Jp,_=embed(theta+eps,mu);Jm,_=embed(theta-eps,mu)
    fd=(dagger(Jp)@b@Jp-dagger(Jm)@b@Jm)/(2*eps)
    assert op(dh-analytic)<1e-12 and op(dh-fd)<1e-8 and op(dh)>.1
    scalar_errors=[]
    for a in range(5):
        delta_phi=np.eye(5)[a]*eps
        hp,dp=original.mass_matrices(PHI+delta_phi);hm,dm=original.mass_matrices(PHI-delta_phi)
        bp,pairp=ambient(PHI+delta_phi);bm,pairm=ambient(PHI-delta_phi)
        left=dagger(J)@((bp-bm)/(2*eps))@J
        right=(expected(PHI+delta_phi,theta,mu)[0]-expected(PHI-delta_phi,theta,mu)[0])/(2*eps)
        scalar_errors.append(op(left-right))
        assert op(left-right)<1e-8
    # Restoring every original mass by explicitly co-rotating its ambient
    # extension is possible on this slice, but changes the physical candidate.
    U=np.zeros((64,64),complex)
    for i in range(16):
        q=CHARGES[module_for_mode(2*i)]
        U[4*i:4*i+4,4*i:4*i+4]=np.cos(q*theta/2)*np.eye(4)+1j*np.sin(q*theta/2)*chiral.GAMMA[mu]
    rotated=U@b@dagger(U);rotpair=U@delta@U.T
    original_h,original_d=original.mass_matrices(PHI)
    restored=max(op(dagger(J)@rotated@J-original_h),op(dagger(J)@rotpair@J.conj()-original_d))
    assert restored<1e-12 and op(rotated-b)>.1
    return dict(link_source_norm=op(dh),link_source_analytic_error=op(dh-analytic),
        link_source_difference_error=op(dh-fd),max_scalar_source_error=max(scalar_errors),
        naive_zero_ambient_derivative_misses_link_source=True,
        corotated_mass_restore_error=restored,ambient_mass_change=op(rotated-b),
        corotation_is_new_field_dependent_mass_not_a_free_basis_change=True)
def run():
    result=dict(rest=rest_dictionary_check(),varying=varying_fibre_check(),
                sources=source_and_adjustment_check())
    deps=('research_note_598.md','research_note_604.md','research_note_606.md',
          'research_note_609.md','research_note_610.md','joint_fermion_gauss_completion.py',
          'joint_chiral_fibre_source.py','joint_gapped_link_locality.py')
    return dict(round=611,tests_run=3,failures=0,errors=0,**result,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_32_mode_mass_used=True,explicit_rest_spinor_dictionary=True,
            field_dependent_projection_changes_mass_blocks_and_sources=True,
            no_global_spinor_frame_or_full_gauge_configuration_extension=True,
            no_Euclidean_transfer_reconstruction_or_chiral_measure=True,
            no_pole_mass_continuum_or_GR_completion=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=611,tests=3,all_passed=True,rest=result['rest'],sources=result['sources'])))
