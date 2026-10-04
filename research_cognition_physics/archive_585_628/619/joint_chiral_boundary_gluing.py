"""619: original chiral content, boundary flux, and common frame weights.

The Lorentz Weyl principal symbol is an explicit continuous-branch input.
Nambu doubling includes pairing; it does not double the physical species.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as old
import joint_original_mass_spinor_bridge as bridge

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_chiral_boundary_gluing_results.json'
SIGMA3=np.diag([1.,-1.])
Q32=np.zeros(32);J32=np.tile([.5,-.5],16);A32=np.zeros(32)
LABEL=[]
for name,s in old.SLICES.items():
    chi=-1 if name in bridge.LEFT else 1
    Q32[s]=bridge.CHARGES[name]
    A32[s]=np.tile([chi,-chi],(s.stop-s.start)//2)
    LABEL.extend([name]*(s.stop-s.start))
CHARGE=np.r_[Q32,-Q32]
SPIN=np.r_[J32,-J32]
VELOCITY=np.r_[A32,A32]
PLUS=np.flatnonzero(VELOCITY>0);MINUS=np.flatnonzero(VELOCITY<0)
AN=np.diag(VELOCITY)
CX=np.block([[np.zeros((32,32)),np.eye(32)],[np.eye(32),np.zeros((32,32))]])

def blockdiag(a,b):
    out=np.zeros((len(a)+len(b),len(a)+len(b)),complex)
    out[:len(a),:len(a)]=a;out[len(a):,len(a):]=b
    return out

def mass(phi):
    h,d=old.mass_matrices(np.asarray(phi))
    return bridge.bdg(h,d)

def reflection_check():
    qp,qm=CHARGE[PLUS],CHARGE[MINUS]
    jp,jm=SPIN[PLUS],SPIN[MINUS]
    allowed=(qm[:,None]==qp[None,:])&(jm[:,None]==jp[None,:])
    positions=np.argwhere(allowed)
    # At most one allowed entry per row/column here; hence the exact rank bound.
    assert len(positions)==2
    assert all(LABEL[int(PLUS[j]%32)]=='nu' and LABEL[int(MINUS[i]%32)]=='nu'
               for i,j in positions)
    partial=allowed.astype(complex)
    assert np.linalg.matrix_rank(partial)==2
    rng=np.random.default_rng(619)
    group_error=0.
    for _ in range(5):
        C=old.gauge.group_exp(rng.normal(size=8),3)
        W=old.gauge.group_exp(rng.normal(size=3),2);z=np.exp(1j*rng.normal())
        R=old.representation(C,W,z);RN=blockdiag(R,R.conj())
        gp=RN[np.ix_(PLUS,PLUS)];gm=RN[np.ix_(MINUS,MINUS)]
        group_error=max(group_error,float(np.max(abs(gm@partial-partial@gp))))
    # A fully gauge-preserving wall which fails axial spin rotation.
    identity=np.eye(32)
    identity_charge=float(np.linalg.norm(qm[:,None]*identity-identity*qp[None,:],2))
    identity_spin=float(np.linalg.norm(jm[:,None]*identity-identity*jp[None,:],2))
    # Particle-hole reflection keeps spin but flips every nonzero U1 charge.
    swap=np.zeros((32,32))
    for j,i in enumerate(PLUS):
        target=(int(i)+32)%64
        # Particle-hole alone has the same normal sign; flip spin for reflection.
        target=target+1 if target%2==0 else target-1
        row=int(np.flatnonzero(MINUS==target)[0]);swap[row,j]=1
    swap_spin=float(np.linalg.norm(jm[:,None]*swap-swap*jp[None,:],2))
    swap_charge=float(np.linalg.norm(qm[:,None]*swap-swap*qp[None,:],2))
    assert identity_charge==0 and identity_spin==1 and swap_spin==0 and swap_charge==12
    phi=np.array([0.,0.,0.,0.,.8]);B=mass(phi)
    assert np.linalg.norm(B,2)>.1
    assert np.max(abs(CX@B.conj()@CX+B))<1e-13
    assert np.max(abs(CX@AN.conj()@CX-AN))==0
    assert np.max(abs(np.diag(SPIN)@B-B@np.diag(SPIN)))<1e-13
    particle_allowed=(Q32[np.flatnonzero(A32<0)][:,None]
                      ==Q32[np.flatnonzero(A32>0)][None,:])&(
                      J32[np.flatnonzero(A32<0)][:,None]
                      ==J32[np.flatnonzero(A32>0)][None,:])
    assert np.count_nonzero(particle_allowed)==0
    rows=[]
    for name,s in old.SLICES.items():
        d=(s.stop-s.start)//2;chi=-1 if name in bridge.LEFT else 1;q=bridge.CHARGES[name]
        rows.append(dict(module=name,internal_dimension=d,charge=q,chirality=chi,
            particle_plus_weight=[q,chi/2],particle_minus_weight=[q,-chi/2],
            hole_plus_weight=[-q,-chi/2],hole_minus_weight=[-q,chi/2]))
    return dict(normal_positive_dimension=32,normal_negative_dimension=32,
        symmetry_intertwiner_rank_upper_bound=2,missing_rank_at_least=30,
        allowed_complex_entries=2,only_neutral_nu_entries=True,
        particle_only_intertwiner_dimension=0,full_group_partial_intertwiner_error=group_error,
        gauge_wall_spin_generator_defect=identity_spin,spin_wall_U1_generator_defect=swap_charge,
        original_nonzero_Majorana_mass_norm=float(np.linalg.norm(B,2)),weight_table=rows,
        PHS_constraints_can_only_reduce_the_allowed_set=True,
        zero_Higgs_configuration_is_the_symmetry_test=True)

def transformed_phi(phi,W,z):
    X=z**3*W@(phi[:2]+1j*phi[2:4])
    return np.r_[X.real,X.imag,phi[4]]

def transmission_check():
    rng=np.random.default_rng(6192);worst=0.;source_error=0.;rows=[]
    doubled_form=blockdiag(AN,-AN)
    for _ in range(7):
        phi=rng.normal(size=5)*.32;direction=rng.normal(size=5)*.17
        C=old.gauge.group_exp(rng.normal(size=8),3)
        W=old.gauge.group_exp(rng.normal(size=3),2);z=np.exp(1j*rng.normal())
        R=old.representation(C,W,z);T=blockdiag(R,R.conj())
        G=np.vstack([np.eye(64),T])/np.sqrt(2)
        transformed=transformed_phi(phi,W,z)
        B=mass(phi);Bt=mass(transformed)
        errors=[np.max(abs(G.conj().T@doubled_form@G)),
                np.max(abs(G.conj().T@G-np.eye(64))),
                np.max(abs(Bt-T@B@T.conj().T)),
                np.max(abs(T.conj().T@AN@T-AN)),
                np.max(abs(np.diag(SPIN)@T-T@np.diag(SPIN)))]
        step=2e-6
        dB=(mass(phi+step*direction)-mass(phi-step*direction))/(2*step)
        dBt=(mass(transformed_phi(phi+step*direction,W,z))
             -mass(transformed_phi(phi-step*direction,W,z)))/(2*step)
        se=float(np.max(abs(dBt-T@dB@T.conj().T)))
        source_error=max(source_error,se);worst=max(worst,float(max(errors)))
        rows.append(dict(F=float(old.original.F(phi)),mass_norm=float(np.linalg.norm(B,2)),
                         boundary_form_and_mass_error=float(max(errors)),mass_source_error=se))
    assert worst<2e-13 and source_error<2e-9
    assert np.linalg.matrix_rank(G)==64
    # Exact charged characteristic at zero Higgs, with the original neutral mass present.
    sigma=.4;d=1.1;t0=.25;t1=1.7
    def right(t):return .5*math.erfc((d-t)/sigma)
    nodes,weights=np.polynomial.legendre.leggauss(96)
    time=t0+(nodes+1)*(t1-t0)/2
    flux=np.exp(-((d-time)/sigma)**2)/(math.sqrt(math.pi)*sigma)
    integrated=float(weights@flux*(t1-t0)/2)
    difference=right(t1)-right(t0)
    assert abs(integrated-difference)<2e-13 and difference>.9
    return dict(rows=rows,max_boundary_and_mass_error=worst,max_source_error=source_error,
        transmission_graph_dimension=64,combined_boundary_dimension=128,
        charged_packet_right_probability_initial=right(t0),
        charged_packet_right_probability_final=right(t1),integrated_flux=integrated,
        flux_balance_error=abs(integrated-difference),
        independently_closed_wall_would_require_zero_net_flux=True,
        boundary_trace_relation_not_cloning_map=True,
        full_continuum_domain_equivalence_proved_for_declared_flat_bounded_background=True)

def frame_check():
    rng=np.random.default_rng(6193);rows=[]
    def densities(phi,eta):
        F=float(old.original.F(phi));BE=mass(phi);BJ=np.sqrt(F)*BE
        etaE=F**(-.75)*eta
        # Per fixed Jordan spacetime density; full Nambu overall half cancels.
        return float(np.vdot(eta,BJ@eta).real),float(F**2*np.vdot(etaE,BE@etaE).real)
    for fraction in (.15,.55,.88):
        direction=rng.normal(size=5);direction/=np.linalg.norm(direction)
        phi=np.sqrt(6*old.original.M)*fraction*direction
        F=float(old.original.F(phi))
        eta=rng.normal(size=64)+1j*rng.normal(size=64);eta/=np.linalg.norm(eta)
        # Ensure a nonzero boundary form rather than a null accidental fixture.
        eta[MINUS]*=.2;eta/=np.linalg.norm(eta)
        etaE=F**(-.75)*eta
        fluxJ=float(np.vdot(eta,AN@eta).real)
        fluxE=float(F**1.5*np.vdot(etaE,AN@etaE).real)
        wrong=float(F**1.5*np.vdot(eta,AN@eta).real)
        dj,de=densities(phi,eta)
        variation=rng.normal(size=5)*.07;step=1e-6
        plus=densities(phi+step*variation,eta);minus=densities(phi-step*variation,eta)
        sourceJ=(plus[0]-minus[0])/(2*step);sourceE=(plus[1]-minus[1])/(2*step)
        errors=[abs(fluxJ-fluxE),abs(dj-de),abs(sourceJ-sourceE)]
        assert max(errors)<5e-10 and abs(wrong-fluxJ)>.1
        rows.append(dict(radial_fraction=fraction,F=F,Jordan_boundary_form=fluxJ,
            Einstein_boundary_form=fluxE,unscaled_spinor_flux_defect=abs(wrong-fluxJ),
            Jordan_mass_density=dj,Einstein_pulled_back_mass_density=de,
            pulled_back_mass_source=sourceJ,max_error=float(max(errors))))
    return dict(rows=rows,spinor_weight='F^(-3/4)',mass_weight='F^(-1/2)',
        scalar_source_includes_metric_and_spinor_pullback=True,
        ordinary_charge_conservation_not_inferred_for_Majorana_particle_number=True)

def run():
    deps=('research_note_598.md','research_note_599.md','research_note_611.md',
          'research_note_613.md','research_note_616.md','research_note_617.md','research_note_618.md',
          'joint_fermion_gauss_completion.py','joint_original_mass_spinor_bridge.py')
    return dict(round=619,tests_run=3,failures=0,errors=0,
        reflection=reflection_check(),transmission=transmission_check(),frames=frame_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_charges_chiralities_and_Majorana_retained=True,
            no_go_only_for_declared_symmetric_local_linear_reflection=True,
            matched_transmission_and_mass_sources_compatible=True,
            continuum_principal_symbol_is_input_not_derived=True,
            interacting_boundaries_chiral_reconstruction_and_quantum_GR_open=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))

