"""666: original reference grading, Majorana stabilizer and transported weights.

Conditional reconstruction inside the already given614 Spin(10) carrier,
not a derivation of the carrier, colour split, neutrino or physical gauge group.
Full32 conditional auxiliary traces are independently computed as four8-mode
Fock blocks for a diagonal-colour but nontrivial weak/hypercharge closure.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_spinor_subgroup_mass as dictionary
import joint_connection_matter_matching as contact
import joint_gauss_fermion_influence as car
import joint_fermion_gauss_completion as matter

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_reference_mass_stabilizer_results.json'
BLOCKS=[list(range(4*c,4*c+4))+[12+2*c,13+2*c,18+2*c,19+2*c] for c in range(3)]+[list(range(24,32))]


def norm(x):return float(np.linalg.norm(x,'fro'))
def expm(h,z):
    ev,u=np.linalg.eigh(h)
    return (u*np.exp(z*ev))@u.conj().T


def real_columns(arrays):return np.column_stack([np.r_[x.real.ravel(),x.imag.ravel()] for x in arrays])


def joint_stabilizer_check():
    generators,annih=dictionary.clifford();pairs=list(itertools.combinations(range(10),2))
    lam=-contact.current_matrices()[0][::2,::2];j=dictionary.dictionary();lc=j@lam@j.conj().T
    parity=np.diag([-(-1)**((m&7).bit_count()) for m in dictionary.MASKS])
    assert norm(lc-parity)==0
    b0=np.zeros((16,16),complex);b0[0,0]=1
    comms=[t@lc-lc@t for t in generators]
    masses=[t@b0+b0@t.T for t in generators]
    mixed=[i for i,(a,b) in enumerate(pairs) if (a<6)!=(b<6)]
    same=[i for i in range(45) if i not in mixed]
    assert len(mixed)==24 and len(same)==21 and max(norm(comms[i]) for i in same)==0
    mc=real_columns([comms[i] for i in mixed]);gram=mc.T@mc
    gram_error=float(np.max(abs(gram-16*np.eye(24))))
    assert gram_error==0
    rank_mass=int(np.linalg.matrix_rank(real_columns(masses),tol=1e-12))
    joint=np.vstack([real_columns(comms),real_columns(masses)])
    rank_joint=int(np.linalg.matrix_rank(joint,tol=1e-12))
    # Explicit s(u(3)+u(2)) basis, lifted from exactly the same5 oscillators.
    basis5=[]
    for ids in (range(3),range(3,5)):
        ids=list(ids)
        for a,b in itertools.combinations(ids,2):
            x=np.zeros((5,5),complex);x[a,b]=x[b,a]=1;basis5.append(x)
            y=np.zeros((5,5),complex);y[a,b]=-1j;y[b,a]=1j;basis5.append(y)
        for a in ids[:-1]:
            x=np.zeros((5,5));x[a,a]=1;x[ids[-1],ids[-1]]=-1;basis5.append(x)
    basis5.append(np.diag([-2.,-2.,-2.,3.,3.]))
    lifted=[]
    for h in basis5:
        t=sum(h[a,b]*annih[a].conj().T@annih[b] for a in range(5) for b in range(5))
        lifted.append(t[np.ix_(dictionary.MASKS,dictionary.MASKS)])
    residual=max(max(norm(t@lc-lc@t),norm(t@b0+b0@t.T)) for t in lifted)
    basis_rank=int(np.linalg.matrix_rank(real_columns(lifted),tol=1e-12))
    assert (rank_mass,rank_joint,basis_rank)==(21,33,12) and residual==0
    return dict(internal_dimension=16,reference_parity_identity_error=0.,
                reference_preserving_generators=21,mixed_generators=24,mixed_Gram_error=gram_error,
                Majorana_stabilizer_Lie_dimension=45-rank_mass,
                joint_stabilizer_Lie_dimension=45-rank_joint,explicit_SM_basis_dimension=basis_rank,
                explicit_basis_residual=residual,extra_central_minus_one_preserves_pair=True,
                connected_group_only_is_SM_quotient=True,
                original_carrier_split_and_mass_line_are_inputs=True)


def transported_auxiliary_check():
    currents=contact.current_matrices();chi=currents[0];spin=currents[1:]/2
    h,d=matter.mass_matrices(car.PHI[0]);w=dictionary.ph_matrix();new_b=w@car.bdg(h,d)@w.conj().T
    delta=.19;a0=.035;fields=np.array([.31,-.27,.43,.18,-.36,.22])
    colour=np.diag(np.exp(1j*np.array([.13,-.07,-.06])))
    weak=matter.gauge.group_exp(np.array([.17,-.09,.21]),2)
    gauge=matter.representation(colour,weak,np.exp(.12j))
    zero=np.zeros((32,32),complex)
    old_heat=[];new_heat=[]
    for ids in BLOCKS:
        old_heat.append(expm(car.fock(h[np.ix_(ids,ids)],d[np.ix_(ids,ids)]),-delta))
        new_heat.append(expm(car.fock(new_b[np.ix_(ids,ids)],new_b[np.ix_(ids,[i+32 for i in ids])]),-delta))
    occupation=np.array([[bool(m&(1<<i)) for i in range(8)] for m in range(256)],int)
    def evaluate(theta,charge_average=False,points=80):
        a=a0*np.exp(-6*theta);t=np.sqrt(2*a);s=t*fields[0]
        rot=np.eye(32,dtype=complex)
        for k,(axis,c) in enumerate(((0,4),(1,4),(2,8),(1,4),(0,4)),start=1):
            rot=rot@expm(spin[axis],1j*np.sqrt(c*a)*fields[k])
        base=gauge@rot
        rb=(w@np.block([[base,zero],[zero,np.linalg.inv(base).T]])@w.conj().T)[:32,:32]
        physical=[];left=[];phase=1.+0j
        for ids,eo,en in zip(BLOCKS,old_heat,new_heat):
            olddiag=np.diag(car.exterior(base[np.ix_(ids,ids)])@eo)
            newdiag=np.diag(car.exterior(rb[np.ix_(ids,ids)])@en)
            n=occupation.sum(axis=1);ch=occupation@np.diag(chi)[ids].real
            relative=-ch
            physical.append((olddiag,n,ch));left.append((newdiag,n,relative))
            right=[i for i in ids if chi[i,i]==1]
            phase*=np.linalg.det(base[np.ix_(right,right)])
        def weights(s,shift=False):
            oldv=1.+0j;newv=1.+0j
            for (od,n,ch),(nd,nb,rel) in zip(physical,left):
                oldv*=np.sum(od*np.exp(2*a*n+s*ch))
                newv*=np.sum(nd*np.exp(2*a*rel-s*nb-(32*a*nb if shift else 0)))
            return oldv/2.**32,newv/2.**32
        oldv,bare=weights(s);factor=phase*np.exp(32*a+16*s);newv=factor*bare
        out=(oldv,newv,bare,factor)
        if charge_average:
            z,qw=np.polynomial.hermite.hermgauss(points);z*=np.sqrt(2);qw/=np.sqrt(np.pi)
            av_old=sum(q*weights(t*x)[0] for x,q in zip(z,qw))
            av_new=sum(q*phase*np.exp(32*a+16*t*x)*weights(t*x)[1] for x,q in zip(z,qw))
            av_shift=sum(q*phase*np.exp(288*a)*weights(t*x,shift=True)[1] for x,q in zip(z,qw))
            return out,(av_old,av_new,av_shift)
        return out
    (oldv,newv,bare,factor),means=evaluate(0.,True)
    _,means64=evaluate(0.,True,64)
    relative=float(abs(oldv-newv)/abs(oldv));assert relative<1e-12
    step=2e-5;plus=evaluate(step);minus=evaluate(-step)
    source_old=(plus[0]-minus[0])/(2*step)/oldv
    source_new=(plus[1]-minus[1])/(2*step)/newv
    source_bare=(plus[2]-minus[2])/(2*step)/bare
    scalar_derivative=-192*a0-48*np.sqrt(2*a0)*fields[0]
    source_error=float(abs(source_old-source_new))
    missing_error=float(abs(source_old-source_bare-scalar_derivative))
    mean_errors=[float(abs(x-means[0])/abs(means[0])) for x in means[1:]]
    mean_errors64=[float(abs(x-means64[0])/abs(means64[0])) for x in means64[1:]]
    convergence=[float(abs(x-y)/abs(x)) for x,y in zip(means,means64)]
    assert source_error<2e-9 and missing_error<1e-6
    assert max(mean_errors+mean_errors64+convergence)<3e-11
    return dict(original_modes=32,independent_Fock_blocks=[8,8,8,8],time_step=delta,a=a0,
                normalized_original_trace=[float(oldv.real),float(oldv.imag)],
                exact_reference_factor=[float(factor.real),float(factor.imag)],
                particle_hole_trace_relative_error=relative,transported_source_error=source_error,
                omitted_reference_source=scalar_derivative,source_factor_identity_error=missing_error,
                normalized_charge_average_errors=mean_errors,
                Gaussian_quadrature_points=[64,80],charge_average_errors64=mean_errors64,
                quadrature_convergence_relative_errors=convergence,
                initial_32_point_failure_preserved='round666_drafts/transport_first_failure.json',
                shifted_Gaussian_mean=16*np.sqrt(2*a0),shifted_constant_exponent=288*a0,
                original_SM_closure_and_all_complex_masses_retained=True,
                auxiliary_spin_fields_held_fixed_in_charge_average=True,
                no_full_Gauss_bosonic_or_chiral_measure_equivalence=True)


def run():
    deps=('joint_spinor_subgroup_mass.py','joint_connection_matter_matching.py',
          'joint_fermion_gauss_completion.py','joint_gauss_fermion_influence.py',
          'research_note_614.md','research_note_653.md','research_note_661.md','research_note_665.md',
          'round666_drafts/particle_hole_contact_probe.py','round666_drafts/particle_hole_contact_probe_results.json')
    return dict(date='2026-10-02',round=666,tests_run=2,failures=0,errors=0,
                joint_reference_mass_stabilizer=joint_stabilizer_check(),
                actual_full_CAR_transport=transported_auxiliary_check(),
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
                scope='Within the given614 Spin(10) carrier, the original reference grading and fixed Majorana line have connected common stabilizer equal to the original SM quotient; central fermion parity remains. Actual32 auxiliary weights, reference constants and geometry sources transported by the same particle-hole map. No selection of carrier/split, full enlarged symmetry, continuum anomaly, original chiral auxiliary measure or quantum GR.',
                all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
