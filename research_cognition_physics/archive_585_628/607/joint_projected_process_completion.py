"""607: a declared form completion for the 606 varying chiral projector.
Diagnostics retain an actual flat-holonomy block of the previous Wilson kernel.
This changes the ambient Hamiltonian; it is not original dynamics equivalence.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import joint_chiral_fibre_source as prior
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_projected_process_completion_results.json'

def build(N,gamma=.12,half_width=1.):
    step=2*half_width/(N+1);theta=-half_width+step*np.arange(1,N+1)
    kappa=np.exp(-2*gamma)
    T=(2*np.eye(N)-np.eye(N,k=1)-np.eye(N,k=-1))/step**2
    ev,vec=np.linalg.eigh(-prior.G5);V0=vec[:,ev>0]
    J=np.zeros((4*N,2*N),complex)
    for i,t in enumerate(theta):
        U=np.cos(t/2)*np.eye(4)+1j*np.sin(t/2)*prior.GAMMA[0]
        J[4*i:4*i+4,2*i:2*i+2]=U@V0
    P=J@J.conj().T;Q=np.eye(4*N)-P
    H=kappa*np.kron(T,np.eye(4))
    completed=P@H@P+Q@H@Q
    h=J.conj().T@H@J
    return dict(theta=theta,step=step,kappa=kappa,T=T,J=J,P=P,H=H,Hb=completed,h=h)

def unitary(H,t):
    E,V=np.linalg.eigh(H)
    return (V*np.exp(-1j*t*E))@V.conj().T

def form_and_spectrum_check():
    rows=[]
    # This is the zero-momentum block of the actual 2^4 kernel, not an arbitrary rotating qubit.
    F=np.kron(np.ones((16,1))/4,np.eye(4))
    projector_error=0.
    for t in (-.3,.0,.2):
        H,_,_=prior.wilson(t,1)
        block=F.conj().T@H@F
        exact=-prior.G5*np.cos(t)+1j*prior.G5@prior.GAMMA[0]*np.sin(t)
        projector_error=max(projector_error,float(np.max(abs(block-exact))))
    assert projector_error<1e-13
    for N in (12,24,48):
        m=build(N);J=m['J'];P=m['P'];H=m['H'];Hb=m['Hb'];h=m['h']
        norm=lambda x:float(np.linalg.norm(x,2))
        intertwining=norm(Hb@J-J@h);leak=norm(H@J-J@h)
        formula=m['kappa']*np.cos(m['step']/2)*np.kron(m['T'],np.eye(2))
        phi=2*m['kappa']*(1-np.cos(m['step']/2))/m['step']**2
        formula+=phi*np.eye(2*N)
        assert intertwining<2e-11 and norm(h-formula)<2e-11
        assert norm(Hb@P-P@Hb)<2e-11 and leak>.1
        minimum=float(np.linalg.eigvalsh(Hb)[0]);assert minimum>0
        lowT=float(np.linalg.eigvalsh(m['T'])[0])
        low_shift=float(np.linalg.eigvalsh(h)[0]-m['kappa']*lowT)
        rows.append(dict(N=N,intertwining_error=intertwining,original_generator_leakage=leak,
            minimum_completed_energy=minimum,discrete_horizontal_potential=phi,
            continuum_horizontal_potential=m['kappa']/4,lowest_energy_shift=low_shift,
            compressed_formula_error=norm(h-formula)))
    assert abs(rows[-1]['lowest_energy_shift']-rows[-1]['continuum_horizontal_potential'])<.0007
    return dict(actual_Wilson_block_error=projector_error,rows=rows,
        completion_is_new_dynamics=True,all_prior_intrasection_forms_retained=True)

def instruments(theta,which,ambient=False):
    e=.5+.2*(np.cos(2*theta) if which==0 else np.sin(2*theta))
    spin=4 if ambient else 2
    return [np.kron(np.diag(np.sqrt(e)),np.eye(spin)),
            np.kron(np.diag(np.sqrt(1-e)),np.eye(spin))]

def full_history_check():
    m=build(18);J=m['J'];N=18
    x=np.arange(1,N+1)
    s1=np.sin(np.pi*x/(N+1));s1/=np.linalg.norm(s1)
    s2=np.sin(2*np.pi*x/(N+1));s2/=np.linalg.norm(s2)
    logical=np.zeros((2*N,2),complex)
    logical[::2,0]=s1/np.sqrt(2);logical[1::2,1]=s2/np.sqrt(2)
    full=J@logical
    U0=unitary(m['Hb'],.14);u0=unitary(m['h'],.14)
    U3=unitary(m['Hb'],.06);u3=unitary(m['h'],.06)
    first=instruments(m['theta'],0,True);first_l=instruments(m['theta'],0)
    error=0.;trace_norm=0.;total=0.;energy_error=0.
    for a in range(2):
        U2=unitary(m['Hb'],.09+.04*a);u2=unitary(m['h'],.09+.04*a)
        second=instruments(m['theta'],a,True);second_l=instruments(m['theta'],a)
        for b in range(2):
            A=U3@second[b]@U2@first[a]@U0
            B=u3@second_l[b]@u2@first_l[a]@u0
            error=max(error,float(np.linalg.norm(A@J-J@B,2)))
            psi=A@full;phi=J@B@logical
            aa=psi.ravel();bb=phi.ravel()
            diff=np.outer(aa,aa.conj())-np.outer(bb,bb.conj())
            trace_norm+=float(np.sum(abs(np.linalg.eigvalsh(diff))))
            total+=float(np.linalg.norm(psi)**2)
            ef=float(np.trace(psi.conj().T@m['Hb']@psi).real)
            el=float(np.trace((B@logical).conj().T@m['h']@(B@logical)).real)
            energy_error=max(energy_error,abs(ef-el))
    assert error<1e-11 and trace_norm<1e-11 and abs(total-1)<1e-11 and energy_error<1e-10
    old_wait=unitary(m['H'],.14)@full
    outside=(np.eye(4*N)-m['P'])@old_wait
    old_leak=float(np.linalg.norm(outside)**2);assert old_leak>1e-4
    return dict(reference_dimension=2,record_branches=4,operator_history_error=error,
        complete_record_state_trace_norm=trace_norm,total_probability=total,
        branch_energy_error=energy_error,old_wait_outside_probability=old_leak,
        arbitrary_reference_guarantee_is_analytic=True)

def common_source_check():
    m=build(16);eps=2e-6;J=m['J']
    dH=(build(16,.12+eps)['Hb']-build(16,.12-eps)['Hb'])/(2*eps)
    source=-2*m['Hb']
    derivative=float(np.max(abs(dH-source)))
    assert derivative<2e-7
    source_map=float(np.max(abs(J.conj().T@source@J+2*m['h'])))
    L=instruments(m['theta'],1,True);l=instruments(m['theta'],1)
    injection=sum(x@m['Hb']@x for x in L)-m['Hb']
    small=sum(x@m['h']@x for x in l)-m['h']
    injection_error=float(np.max(abs(J.conj().T@injection@J-small)))
    source_injection=sum(x@source@x for x in L)-source
    assert max(source_map,injection_error)<1e-11
    assert np.max(abs(source_injection+2*injection))<1e-11
    return dict(source_derivative_error=derivative,source_map_error=source_map,
        instrument_energy_map_error=injection_error,
        finite_grid_injection_not_claimed_positive_for_all_states=True,
        original_continuum_scalar_instrument_identity_restricted_in_note=True)

def run():
    data=dict(completion=form_and_spectrum_check(),history=full_history_check(),sources=common_source_check())
    deps=('research_note_363.md','research_note_459.md','research_note_506.md','research_note_592.md',
          'research_note_598.md','research_note_603.md','research_note_605.md','research_note_606.md',
          'joint_chiral_fibre_source.py','cognitive_foundation_bridge_605.md')
    return dict(round=607,tests_run=3,failures=0,errors=0,**data,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(declared_form_completion_not_original_H=True,
            conditional_Gauss_heat_and_source_closure=True,
            full_history_and_reference_mapping=True,
            actual_flat_GW_block_is_diagnostic=True,
            no_locality_continuum_chiral_measure_or_GR_completion=True,
            actual_preparation_and_full_FUCP_rights_open=True))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(dict(round=607,tests=3,all_passed=True,history=r['history'],sources=r['sources'])))
