"""720: original massive Weyl propagation, mean spin, records and mass profiles.

Exact identities in the 633 constant-background neutral BdG branch. No claim
of an interacting continuum limit, local FW instrument, or emergent geometry.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_continuum_record_sources as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_momentum_spin_reference_results.json'
B0=old.B0
I=np.eye(8,dtype=complex)
ALPHA=np.array([old.bdg(e)-B0 for e in np.eye(3)])
S0=np.array([np.block([[np.kron(np.eye(2),s)/2,np.zeros((4,4))],
                       [np.zeros((4,4)),-np.kron(np.eye(2),s.T)/2]])
             for s in old.PAULI])
ev,vec=np.linalg.eigh(B0)
M=(vec*abs(ev))@vec.conj().T
MI=(vec/abs(ev))@vec.conj().T
BETA=B0@MI
MASSES=old.mass_data()[0]
MASS_P=[(M-m*I)/(MASSES[j]-m) for j,m in enumerate(MASSES[::-1])]


def norm(a):
    return float(np.linalg.norm(a,2))


def comm(a,b):
    return a@b-b@a


def block(*items):
    out=np.zeros((sum(len(x) for x in items),)*2,complex)
    pos=0
    for x in items:
        out[pos:pos+len(x),pos:pos+len(x)]=x
        pos+=len(x)
    return out


def fw(k):
    k=np.asarray(k,float)
    K=np.einsum('a,aij->ij',k,ALPHA)
    energies=np.sqrt(k@k+ev**2)
    E=(vec*energies)@vec.conj().T
    invden=(vec/np.sqrt(2*energies*(energies+abs(ev))))@vec.conj().T
    U=(E+M+BETA@K)@invden
    return U,BETA@E


def spin(k):
    U,_=fw(k)
    return np.array([U.conj().T@s@U for s in S0])


def free_check():
    rng=np.random.default_rng(720)
    errs=[]
    for k in np.vstack((np.zeros(3),rng.normal(size=(15,3)))):
        U,D=fw(k);B=old.bdg(k);s=spin(k)
        errs += [norm(U@U.conj().T-I),norm(U@B@U.conj().T-D),
                 norm(sum(a@a for a in s)-.75*I),
                 norm(old.JCH@fw(-k)[0].conj()@old.JCH-U)]
        for a in range(3):
            errs += [norm(comm(s[a],B)),norm(s[a]-s[a].conj().T),
                     norm(comm(s[a],s[(a+1)%3])-1j*s[(a+2)%3])]
        axis=rng.normal(size=3);axis/=np.linalg.norm(axis);angle=.71
        uu=np.cos(angle/2)*np.eye(2)-1j*np.sin(angle/2)*np.einsum('a,aij->ij',axis,old.PAULI)
        V=block(np.kron(np.eye(2),uu),np.kron(np.eye(2),uu.conj()))
        rotation=np.array([[np.trace(old.PAULI[a]@uu@old.PAULI[b]@uu.conj().T).real/2
                            for b in range(3)] for a in range(3)])
        errs.append(norm(fw(rotation@k)[0]-V@U@V.conj().T))
    deltas=[norm(spin(np.array([0.,0.,t]))[0]-S0[0]) for t in (.1,.4,1.)]
    assert max(errs)<2e-13 and min(deltas)>.1
    return dict(neutral_masses=MASSES.tolist(),max_identity_error=max(errs),
                spin_x_difference_at_kz_01_04_1=deltas,
                real_structure_preserved=True,
                exactly_original_Dirac_and_Majorana=True)


def record_data(lam):
    k=np.array([.27,.41,.19])*np.exp([lam,-lam,0.])
    U1,D1=fw(k);U2,D2=fw(-k)
    U=block(U1,U2);D=block(D1,D2)
    B=block(old.bdg(k),old.bdg(-k))
    C=np.block([[np.zeros((8,8)),old.JCH],[old.JCH,np.zeros((8,8))]])
    f=np.zeros(16,complex);f[2]=np.sqrt(.75);f[10]=np.exp(.31j)*np.sqrt(.25)
    cf=C@f.conj()
    Q=np.outer(f,f.conj())+np.outer(cf,cf.conj())
    R=np.eye(16)-2*Q
    P=block(old.spectral(old.bdg(k))[1],old.spectral(old.bdg(-k))[1])
    PF=U@P@U.conj().T
    RF=U@R@U.conj().T
    post=(P+R@P@R)/2
    postF=(PF+RF@PF@RF)/2
    naive=(PF+R@PF@R)/2
    energy=float(np.trace(B@(post-P)).real/2)
    errs=[norm(R@R-np.eye(16)),norm(C@U.conj()@C-U),
          norm(PF-block((I-BETA)/2,(I-BETA)/2)),
          norm(postF-U@post@U.conj().T),
          abs(energy-np.trace(D@(postF-PF)).real/2)]
    return locals()


def transport_check():
    lam=.17;d=record_data(lam);h=2e-5
    up=record_data(lam+h);dn=record_data(lam-h)
    Udot=(up['U']-dn['U'])/(2*h)
    Ddot=(up['D']-dn['D'])/(2*h)
    G=(up['B']-dn['B'])/(2*h)
    RFdot=(up['RF']-dn['RF'])/(2*h)
    Gamma=Udot@d['U'].conj().T
    source_err=norm(Ddot-d['U']@G@d['U'].conj().T-comm(Gamma,d['D']))
    record_err=norm(RFdot-comm(Gamma,d['RF']))
    source_only=float(np.trace(Ddot@(d['postF']-d['PF'])).real/2)
    instrument_term=float(np.trace(d['D']@(RFdot@d['PF']@d['RF']+
                                               d['RF']@d['PF']@RFdot)).real/4)
    fd=(up['energy']-dn['energy'])/(2*h)
    sdiff=spin(.5*np.array([0.,0.,1.]))[0]-spin(-.5*np.array([0.,0.,1.]))[0]
    # Exact Fourier block of commutator with cos(q.x), not a locality fit.
    localization_block=norm(sdiff)/2
    # Independent 8-mode Fock implementation checks the Nambu energy factor.
    order=list(range(4))+list(range(8,12))+list(range(12,16))+list(range(4,8))
    bs=d['B'][np.ix_(order,order)]
    HF=old.finite.fock((bs[:8,:8],bs[:8,8:]))
    ef,vf=np.linalg.eigh(HF);vac=vf[:,0]
    aa=old.finite.annihilators(8)
    mode=np.sqrt(.75)*aa[2]+np.exp(-.31j)*np.sqrt(.25)*aa[6]
    nn=mode.conj().T@mode
    branches=[(np.eye(256)-nn)@vac,nn@vac]
    exact_energy=sum(np.vdot(v,HF@v).real for v in branches)-ef[0]
    errs=d['errs']+[source_err,record_err,abs(fd-source_only-instrument_term),
                   abs(exact_energy-d['energy'])]
    assert max(errs)<2e-8
    assert norm(d['postF']-d['naive'])>.05
    assert abs(instrument_term)>1e-3 and localization_block>.1
    return dict(max_transport_error=max(errs),original_record_energy=d['energy'],
                unchanged_record_wrong_poststate_opnorm=norm(d['postF']-d['naive']),
                fixed_direction_localization_commutator_block=localization_block,
                geometry_derivative_finite_difference=fd,
                transported_energy_derivative_term=source_only,
                actual_instrument_derivative_term=instrument_term,
                omitted_connection_source_opnorm=norm(comm(Gamma,d['D'])),
                source_identity_error=source_err,
                genuine_two_momentum_CAR_record=True,
                full_fock_unitary_in_infinite_volume_not_assumed=True)


def profile_check():
    # Exact commutant calculation on original matrices; analytic proof in note.
    mats=[B0,*ALPHA]
    columns=[]
    for j in range(64):
        e=np.zeros((8,8),complex);e.flat[j]=1
        columns.append(np.concatenate([comm(e,a).ravel() for a in mats]))
    sv=np.linalg.svd(np.array(columns).T,compute_uv=False)
    nullity=int(np.count_nonzero(sv<1e-11))
    assert nullity==2 and sv[-3]>.1
    basis_errors=[]
    for P in MASS_P:
        vp,wp=np.linalg.eigh(P);V=wp[:,vp>.5]
        gamma=[V.conj().T@a@V for a in (*ALPHA,BETA)]
        cliff=[]
        for mask in range(16):
            a=np.eye(4,dtype=complex)
            for j in range(4):
                if mask>>j&1:a=a@gamma[j]
            cliff.append(a)
        gram=np.array([[np.trace(a.conj().T@b)/4 for b in cliff] for a in cliff])
        basis_errors.append(norm(gram-np.eye(16)))
    q=np.array([.21,.16,.13]);mom=[j*q for j in range(-2,3)]
    H=block(*[old.bdg(k) for k in mom]);S=block(*[spin(k)[0] for k in mom])
    eps=.17
    for j in range(4):
        H[j*8:(j+1)*8,(j+1)*8:(j+2)*8]+=eps*B0/2
        H[(j+1)*8:(j+2)*8,j*8:(j+1)*8]+=eps*B0/2
    witness=norm(comm(H,S))
    assert witness>.001 and max(basis_errors)<2e-13
    # c(x) B0 is an exact original radial field path, not an added coupling.
    targetM=old.bg.matter.original.M
    r0=float(np.linalg.norm(old.PHI));direction=old.PHI/r0
    tau0=r0/np.sqrt(old.bg.matter.original.F(old.PHI))
    radial_errors=[]
    for c in (.83,1.,1.17):
        tau=c*tau0
        phi=np.sqrt(targetM)*tau/np.sqrt(1+tau*tau/6)*direction
        hh,dd=old.bg.matter.mass_matrices(phi)
        h0=hh[np.ix_(old.IDS,old.IDS)];de=dd[np.ix_(old.IDS,old.IDS)]
        bb=np.block([[h0,de],[de.conj().T,-h0.T]])
        radial_errors.append(norm(bb-c*B0))
    assert max(radial_errors)<2e-14
    return dict(common_complex_commutant_dimension=nullity,
                smallest_nonzero_commutant_singular_value=float(sv[-3]),
                mass_sector_Clifford_gram_errors=basis_errors,
                profile_amplitude=eps,profile_wavevector=q.tolist(),
                exact_Fourier_compression_commutator_opnorm=witness,
                original_radial_mass_path_error=max(radial_errors),
                witness_is_matrix_element_of_continuum_operator=True,
                no_claim_all_adapted_or_orbital_references_excluded=True)


def run():
    results=dict(free_massive_reference=free_check(),
                 records_locality_and_geometry=transport_check(),
                 original_nonuniform_profiles=profile_check())
    names=('joint_continuum_record_sources.py','joint_continuum_record_sources_results.json',
           'round720_drafts/spin_propagation_entry.py','round720_drafts/spin_propagation_entry_results.json',
           'joint_local_relational_record.py','research_note_633.md','research_note_719.md')
    return dict(round=720,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope='Original constant-background neutral BdG: exact massive FW spin, transported CAR record and geometry; no fixed translation-invariant SU2 reference for every original radial mass profile. Not a full interacting or gravitational no-go.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
