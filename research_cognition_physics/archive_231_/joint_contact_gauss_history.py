"""665: original full-CAR auxiliary lift, normalized sources and spatial flow.

Whole-graph bounded-perturbation and path statements are proved in the note.
Numerics: original32 conditional lift, eight-lepton auxiliary integrations,
and the actual invariant two-node neutral sector. No full bosonic/Gauss integral.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_connection_matter_matching as old
import joint_gauss_fermion_influence as car
import joint_fermion_gauss_completion as matter

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_contact_gauss_history_results.json'


def expm(h,z):
    ev,u=np.linalg.eigh(h)
    return (u*np.exp(z*ev))@u.conj().T


def norm(a):return float(np.linalg.norm(a,'fro'))
def jet_product(x,y):return x[0]@y[0],x[1]@y[0]+x[0]@y[1]


def averaged_exponential(op,amplitude,points=24):
    """Normalized N(0,1). Return derivative with amplitude'=-3 amplitude.

    Independent physical-variable score includes the +3 normalization term.
    """
    z,w=np.polynomial.hermite.hermgauss(points);z*=np.sqrt(2);w/=np.sqrt(np.pi)
    ev,u=np.linalg.eigh(op)
    values=np.exp(amplitude*z[:,None]*ev[None,:])
    mean=w@values
    derivative=w@(-3*amplitude*z[:,None]*ev[None,:]*values)
    score=(w*(3-3*z*z))@values
    omitted=(w*(-3*z*z))@values
    matrix=lambda d:(u*d)@u.conj().T
    return (matrix(mean),matrix(derivative)),matrix(score),matrix(omitted)


def auxiliary_lift_check():
    currents=old.current_matrices();chi=currents[0];spin=currents[1:]/2
    phi=car.PHI[0];h,d=matter.mass_matrices(phi);delta=.19;a0=.035
    field=np.array([.31,-.27,.43,.18,-.36,.22])
    gauge=matter.representation(matter.gauge.group_exp(np.array([.13,-.06,.08,.09,-.12,.07,.04,.11]),3),
          matter.gauge.group_exp(np.array([.17,-.09,.21]),2),np.exp(.12j))
    zero=np.zeros_like(h);bdg=car.bdg(h,d);heat_n=expm(bdg,-delta)
    heat_l=expm(car.fock(h[24:,24:],d[24:,24:]),-delta);heat_q=expm(h[:24,:24],-delta)
    def evaluate(theta):
        a=a0*np.exp(-6*theta)
        r=np.exp(2*a)*expm(chi,np.sqrt(2*a)*field[0])
        for k,(axis,c) in enumerate(((0,4),(1,4),(2,8),(1,4),(0,4)),start=1):
            r=r@expm(spin[axis],1j*np.sqrt(c*a)*field[k])
        r=gauge@r
        fock_trace=np.linalg.det(np.eye(24)+r[:24,:24]@heat_q)*np.trace(car.exterior(r[24:,24:])@heat_l)
        # Nonunitary CAR lift has inverse transpose, NOT complex conjugate.
        rn=np.block([[r,zero],[zero,np.linalg.inv(r).T]])
        mn=rn@heat_n
        determinant=np.linalg.det(np.eye(64)+mn)/2.**64
        normalized=fock_trace/2.**32
        return normalized,determinant,mn,a
    value,determinant,mn,a=evaluate(0.)
    error=float(abs(value*value-np.exp(64*a)*determinant))
    relative=error/max(abs(value*value),1e-20)
    assert relative<3e-12
    step=2e-5;plus=evaluate(step);minus=evaluate(-step)
    derivative=(plus[0]-minus[0])/(2*step)/value
    dm=(plus[2]-minus[2])/(2*step)
    bare=.5*np.trace(np.linalg.solve(np.eye(64)+mn,dm))
    complete=bare-192*a
    source_error=float(abs(complete-derivative))
    assert source_error<3e-6 and abs(bare-derivative)>1
    return dict(original_modes=32,auxiliary_real_fields=field.tolist(),time_step=delta,a=a,
                normalized_trace=[float(value.real),float(value.imag)],
                lifted_trace_square_relative_error=relative,
                missing_lift_trace_factor=float(np.exp(32*a)),
                full_source_error=source_error,omitted_lift_source_error=float(abs(bare-derivative)),
                required_lift_geometric_derivative=-192*a,
                conditional_field_not_integrated_full_Gauss_process=True)


def auxiliary_measure_check():
    mats=old.current_matrices()[:,24:,24:];q,raw,js=old.contact_operator(mats)
    chi=js[0];spin=[x/2 for x in js[1:]]
    number=car.fock(np.eye(8),np.zeros((8,8)));a=.021
    scalar=expm(number,2*a);pair=(scalar,-12*a*number@scalar)
    wrong=pair;score_pair=pair
    factors=[(chi,np.sqrt(2*a))]+[(spin[i],1j*np.sqrt(c*a)) for i,c in ((0,4),(1,4),(2,8),(1,4),(0,4))]
    score_errors=[]
    for op,amp in factors:
        jet,score,omitted=averaged_exponential(op,amp)
        score_errors.append(norm(score-jet[1]))
        pair=jet_product(pair,jet);score_pair=jet_product(score_pair,(jet[0],score))
        wrong=jet_product(wrong,(jet[0],omitted))
    # Independent closed-factor construction (spin symmetric product).
    closed=expm(number,2*a)@expm(chi@chi,a)
    for i,c in ((0,2),(1,2),(2,4),(1,2),(0,2)):
        closed=closed@expm(spin[i]@spin[i],-c*a)
    factor_error=norm(pair[0]-closed)
    minimum=float(np.linalg.eigvalsh(pair[0])[0]);assert minimum>0
    h,d=matter.mass_matrices(car.PHI[1]);e=expm(car.fock(h[24:,24:],d[24:,24:]),-.17)
    trace=np.trace(e@pair[0]);source=np.trace(e@pair[1])/trace
    bad=np.trace(e@wrong[1])/trace
    assert max(score_errors)<3e-12 and factor_error<3e-12
    assert abs(bad-source+18)<1e-11
    assert norm(score_pair[1]-pair[1])<1e-10
    return dict(conditional_lepton_modes=8,Gaussian_points_per_factor=24,independent_auxiliary_fields=6,
                mean_factor_error=factor_error,max_physical_score_error=max(score_errors),
                positive_mean_step_minimum=minimum,
                source=[float(source.real),float(source.imag)],
                dropped_normalization_source_shift=[float((bad-source).real),float((bad-source).imag)],
                mean_spin_Stranger_error=norm(pair[0]-expm(q,-a)),
                finite_spin_splitting_not_exact_rotation_invariant=True,
                positive_average_does_not_imply_positive_path_integrands=True)


def neutral_spatial_check():
    phis=np.array([[0.,0.,0.,0.,.62],[0.,0.,0.,0.,-.41]])
    h=np.zeros((4,4),complex);d=np.zeros_like(h);h[:2,2:]=.23*np.eye(2);h[2:,:2]=.23*np.eye(2)
    onsite=[];contacts=[];spins=[];theta=np.array([.13,-.08]);zero=np.zeros((4,4),complex)
    for node,phi in enumerate(phis):
        hh,dd=matter.mass_matrices(phi);sl=slice(2*node,2*node+2)
        part=np.zeros((4,4),complex);part[sl,sl]=dd[30:,30:];d+=part
        onsite.append(car.fock(zero,part))
        mats=np.zeros((4,4,4),complex);mats[0,sl,sl]=np.eye(2)
        for i in range(3):mats[i+1,sl,sl]=old.PAULI[i]
        q,_,js=old.contact_operator(mats)
        contacts.append(3*q/(16*np.exp(6*theta[node])));spins.append(js[3]/2)
    hopping=car.fock(h,zero);h0=car.fock(h,d);v=sum(contacts);total=h0+v
    old_local=[onsite[i]+hopping/2 for i in range(2)]
    new_local=[old_local[i]+contacts[i] for i in range(2)]
    comm=lambda x,y:(x@y-y@x)/1j
    old_current=comm(*old_local);current=comm(*new_local)
    predicted=comm(contacts[0],hopping/2)+comm(hopping/2,contacts[1])
    error=norm(current-old_current-predicted)
    assert error<1e-12 and norm(predicted)>.1
    local_spin_comm=norm(spins[0]@hopping-hopping@spins[0])
    global_spin_comm=norm(sum(spins)@total-total@sum(spins))
    assert local_spin_comm>.1 and global_spin_comm<1e-12
    beta=.83;e=expm(total,-beta);rho=e/np.trace(e)
    source=6*beta*np.trace(contacts[0]@rho).real;step=2e-5
    def logz(t):
        ev=np.linalg.eigvalsh(h0+np.exp(-6*t)*contacts[0]+contacts[1])
        return float(np.log(np.exp(-beta*ev).sum()))
    fd=(logz(step)-logz(-step))/(2*step)
    assert abs(fd-source)<2e-7
    return dict(actual_original_neutral_modes=4,original_node_singlets=phis[:,4].tolist(),
                hopping=.23,geometry_independent_hopping_is_allowed_diagnostic_input=True,
                contact_energy_current_norm=norm(predicted),complete_current_identity_error=error,
                local_spin_hopping_commutator=local_spin_comm,global_spin_commutator=global_spin_comm,
                joint_thermal_logZ=logz(0),first_node_contact_geometry_source=float(source),
                independent_source_difference_error=abs(fd-source),
                stationary_current_mean=float(np.trace(rho@current).real),
                stationary_current_second_moment=float(np.trace(rho@current@current).real),
                no_full_bosonic_path_integral_claim=True)


def run():
    deps=('joint_connection_matter_matching.py','joint_gauss_fermion_influence.py',
          'joint_fermion_gauss_completion.py','research_note_603.md','research_note_623.md',
          'research_note_624.md','research_note_643.md','research_note_655.md','research_note_664.md',
          'round665_drafts/spin_charge_auxiliary_probe.py','round665_drafts/spin_charge_auxiliary_probe_results.json')
    return dict(date='2026-10-02',round=665,tests_run=3,failures=0,errors=0,
                full_CAR_auxiliary_lift=auxiliary_lift_check(),normalized_auxiliary_measure=auxiliary_measure_check(),
                original_spatial_current=neutral_spatial_check(),
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
                scope='Declared contact extension on the original positive-geometry finite graph. Whole-graph common domain, Gauss heat state, original scalar records and ordered auxiliary representation proved in note; numerics use actual32 conditional lift, eight-lepton auxiliary mean and invariant two-node neutral spatial sector. No sign-free integrand, unique torsion measure, continuum limit or quantum GR.',
                all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
