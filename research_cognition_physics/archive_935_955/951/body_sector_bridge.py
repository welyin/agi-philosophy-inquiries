"""951: number-conserving material sectors preserve the actual shared Hamiltonian.

The numerical example tests the sector algebra with two actual (h,B) spectral
levels and the inherited two-mode portal. It is not a new continuum, detector,
SM or all-order gravity certificate. Full protocol transport is an analytic
sector equivalence applied to the frozen round-947 Hamiltonian.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'body_sector_bridge_results.json'
def read(p): return json.loads(p.read_text('utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a): return float(np.linalg.norm(a))
def unitary(h,t):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T
def destroy(n):
    return np.diag(np.sqrt(np.arange(1,n)),1)
def kron3(a,b,c): return np.kron(np.kron(a,b),c)

def run():
    old=read(STAGE/'947/protocol_field_transport_results.json')
    portal=read(STAGE/'946/portal_common_process_results.json')
    bridge=read(STAGE/'948/neutral_matter_bridge_results.json')
    assert all(x['all_scientific_checks_passed'] for x in (old,portal,bridge))
    K=np.array(bridge['parameters']['inherited_K'])
    mass2,O=np.linalg.eigh(K)
    omega=np.sqrt(mass2)
    cutoff=3
    a=destroy(cutoff)
    eye=np.eye(cutoff)
    b=[np.kron(a,eye),np.kron(eye,a)]
    number=sum(x.T@x for x in b)
    occupations=np.diag(number).astype(int)
    one=np.flatnonzero(occupations==1)
    vac=np.flatnonzero(occupations==0)
    assert len(one)==2 and len(vac)==1
    # Ordering of the one-body basis is inherited from occupation enumeration.
    # b_i^dag |0> fixes the basis map; do not assume sorted indices give it.
    W=np.column_stack([x.T[:,vac[0]] for x in b])
    def dgamma(mat):
        return sum(mat[i,j]*(b[i].T@b[j]) for i in range(2) for j in range(2))
    I2=np.eye(2)
    I9=np.eye(9)
    Z=np.diag([1.,-1.])
    Y=np.array([[0,-1j],[1j,0]])
    X=np.array([[0,1.],[1.,0]])
    kappa=old['same_physical_protocol']['internal_clock_scale']
    hint=np.diag([0.,kappa])
    ms,mr=20.,30. # diagnostic masses, NOT replacement of the frozen 947 masses
    M=ms*I2+hint
    mode_a=[np.kron(a,eye),np.kron(eye,a)]
    Hfield=sum(w*x.T@x for w,x in zip(omega,mode_a))
    coordinates=[(x+x.T)/math.sqrt(2*w) for x,w in zip(mode_a,omega)]
    pfield=sum(O[0,i]*coordinates[i] for i in range(2))
    etafield=sum(O[1,i]*coordinates[i] for i in range(2))
    gS,gR,G,u=.5,.5,2e-4,.5
    def full_h(g=gS,grav=G,source=0.):
        hs=dgamma(M+source*X)
        hr=mr*number
        return (kron3(hs,I9,I9)+kron3(I9,hr,I9)+kron3(I9,I9,Hfield)
                +g*kron3(dgamma(Z),I9,pfield)+gR*kron3(I9,dgamma(Y),pfield)
                -grav*u*kron3(hs,hr,I9))
    def first_h(g=gS,grav=G,source=0.):
        msource=M+source*X
        return (kron3(msource,I2,I9)+kron3(I2,mr*I2,I9)+kron3(I2,I2,Hfield)
                +g*kron3(Z,I2,pfield)+gR*kron3(I2,Y,pfield)
                -grav*u*kron3(msource,mr*I2,I9))
    H=full_h()
    h=first_h()
    W11=kron3(W,W,I9)
    W00=kron3(np.eye(9)[:,[0]],np.eye(9)[:,[0]],I9)
    Ns=kron3(number,I9,I9);Nr=kron3(I9,number,I9)
    sector_error=norm(H@W11-W11@h)
    vacuum_error=norm(H@W00-W00@Hfield)
    number_error=max(norm(H@Ns-Ns@H),norm(H@Nr-Nr@H))
    assert max(sector_error,vacuum_error,number_error)<1e-12
    U=unitary(H,.7);smallU=unitary(h,.7)
    full_U_error=norm(U@W11-W11@smallU)
    vacuum_U_error=norm(U@W00-W00@unitary(Hfield,.7))
    assert max(full_U_error,vacuum_U_error)<2e-12
    step=1e-4
    source_errors={}
    vacuum_source_errors={}
    for name,plus,minus,splus,sminus in (
        ('record_coupling',full_h(g=gS+step),full_h(g=gS-step),first_h(g=gS+step),first_h(g=gS-step)),
        ('Newton_coupling',full_h(grav=G+step),full_h(grav=G-step),first_h(grav=G+step),first_h(grav=G-step)),
        ('noncommuting_internal_source',full_h(source=step),full_h(source=-step),first_h(source=step),first_h(source=-step))):
        derivative=(plus-minus)/(2*step)
        small_derivative=(splus-sminus)/(2*step)
        source_errors[name]=norm(derivative@W11-W11@small_derivative)
        vacuum_source_errors[name]=norm(derivative@W00)
    assert max(source_errors.values())<2e-10 and max(vacuum_source_errors.values())<1e-12
    assert norm(M@X-X@M)>.1
    # Energy exchange is present at finite particle number. It is not
    # erased by the exact empty-body sector identity.
    A=kron3(Z,I2,pfield)
    J=1j*(h@A-A@h)
    current_eigs=np.linalg.eigvalsh(J)
    current_max=float(max(abs(current_eigs)))
    assert current_max>.1
    scalar_read=kron3(I2,I2,etafield)
    scalar_effect_velocity=norm(1j*(h@scalar_read-scalar_read@h))
    # A number-changing term is an explicit negative control. On the full
    # S Fock space, the body vacuum no longer stays vacuum.
    pair=b[0].T@b[1].T
    pair_term=.2*(pair+pair.T)
    Rvac_indices=np.array([s*81+f for s in range(9) for f in range(9)])
    H_Rvac=H[np.ix_(Rvac_indices,Rvac_indices)]
    bad=H_Rvac+np.kron(pair_term,I9)
    initial=np.eye(81)[:,0]
    bad_out=unitary(bad,.7)@initial
    bad_empty_probability=float(np.sum(abs(bad_out[:9])**2))
    pair_leak=1-bad_empty_probability
    assert pair_leak>1e-7
    # Explicit inherited data are copied for provenance, not recomputed as
    # new task tests. The analytic theorem transports the whole original H.
    inherited=old['finite_joint_control']
    files=[Path(__file__),STAGE/'947/protocol_field_transport_results.json',
           STAGE/'946/portal_common_process_results.json',STAGE/'948/neutral_matter_bridge_results.json',
           STAGE/'944/finite_field_communication.py',STAGE/'950/material_vacuum_matching_results.json']
    return dict(round=951,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(material_representation='normal-ordered, fixed-body-number effective theory',
            numerical_internal_levels=2,numerical_body_mode_max_occupation=2,scalar_mode_max_occupation=2,
            numerical_Fock_dimension=729,one_each_dimension=36,diagnostic_body_masses=[ms,mr],
            actual_947_internal_dimension=old['same_physical_protocol']['original_work_and_fault_and_clock_dimension'],
            source_K=K.tolist(),internal_level_spacing=kappa),
        checks=dict(H_intertwining_error=sector_error,empty_sector_generator_error=vacuum_error,
            material_number_commutator_error=number_error,full_unitary_intertwining_error=full_U_error,
            empty_sector_unitary_error=vacuum_U_error,physical_source_intertwining_errors=source_errors,
            empty_sector_material_source_errors=vacuum_source_errors,
            noncommuting_internal_source_commutator=norm(M@X-X@M),
            finite_material_field_exchange_current_norm=current_max,
            scalar_readout_velocity_norm=scalar_effect_velocity,
            number_changing_negative_control_vacuum_leak=pair_leak),
        inherited_by_exact_sector_identity=dict(
            full_cq_instrument_distance_bound=inherited['full_cq_instrument_distance_bound'],
            record_probability_contrast_lower=inherited['record_probability_contrast_lower'],
            higgs_probability_contrast_lower=inherited['higgs_probability_contrast_lower'],
            gravity_probability_contrast_lower=inherited['gravity_probability_contrast_lower'],
            bounds_are_frozen_947_not_new_numerical_estimates=True),
        scope=dict(exact_fixed_sector_representation_theorem=True,
            body_empty_vacuum_gives_no_additional_material_loops_in_this_EFT=True,
            high_energy_matching_coefficients_must_be_kept=True,
            all_internal_levels_reinterpreted_as_fundamental_Dirac_species=False,
            microscopic_binding_or_SM_material_construction_proved=False,
            number_changing_or_pair_creation_domain_covered=False,
            full_SM_Einstein_parent_matching_certified=False,
            new_space_or_gravity_derivation_claimed=False,
            diagnostic_oscillator_cutoff_is_physical_minimum_scale=False,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))

