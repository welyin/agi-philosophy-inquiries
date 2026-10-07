"""865: original quotient-group fusion constrains loop normalization.
The existing constrained 863 color family is used. Proper lengths are evaluated
on the original 859 constraint background, not on a replacement flat geometry.
The perimeter normalization is a declared finite candidate, not a computed
counterterm for the complete interacting theory.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'relational_loop_fusion_matching_results.json'
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(HERE.parent/'863'))
from research_layout import Layout,ResearchRuntime
import physical_relational_loop_bridge as inherited
import loop_self_contact_probe as contact

def matrices(old,c,L):
    data=inherited.color_at(old,0.,c);A=-1j*data['A']
    dA=-1j*np.array([c['a_prime']*old.T[0],np.zeros((3,3)),c['c']*old.T[3]])
    V=A[0]-A[1];Z=A[2];dV=dA[0]-dA[1];dZ=dA[2]
    U=np.eye(3,dtype=complex);dU=np.zeros_like(U)
    for B,dB in ((-L*V,-L*dV),(-L*Z,-L*dZ),(L*V,L*dV),(L*Z,L*dZ)):
        C,dC=inherited.exp_with_tangent(B,dB);dU=dC@U+C@dU;U=C@U
    adj=np.array([[2*np.trace(Ta@U@Tb@U.conj().T) for Tb in old.T] for Ta in old.T])
    dadj=np.array([[2*np.trace(Ta@(dU@Tb@U.conj().T+U@Tb@dU.conj().T)) for Tb in old.T] for Ta in old.T])
    return U,dU,adj,dadj

def real_trig_basis(points,N):
    freq=np.fft.fftfreq(N,d=1/N)
    answer=[np.exp(1j*points[:,i,None]*freq[None,:]) for i in range(3)]
    # Even grids have a single stored Nyquist coefficient. Split it equally
    # between +N/2 and -N/2; the real cosine interpolant agrees at every node.
    if N%2==0:
        for i in range(3):answer[i][:,N//2]=np.cos((N//2)*points[:,i])
    return answer

def interpolation_check(psi):
    N=psi.shape[0];coeff=np.fft.fftn(psi*psi)/N**3
    indices=np.array([[0,0,0],[1,2,3],[N-1,N//2,0],[N//3,N//2,N-2]])
    points=indices*(2*np.pi/N)
    interp=np.einsum('ijk,bi,bj,bk->b',coeff,*real_trig_basis(points,N),optimize=True)
    error=float(np.max(abs(interp-(psi*psi)[tuple(indices.T)])))
    assert error<1e-12
    return error

def proper_length(psi,L,order):
    N=psi.shape[0];coeff=np.fft.fftn(psi*psi)/N**3
    nodes,weights=np.polynomial.legendre.leggauss(order)
    us=(nodes+1)/2;weights=weights/2
    center=np.array([0.,np.pi/2,np.pi/4])
    v=L*np.array([1.,-1.,0.]);w=L*np.array([0.,0.,1.])
    vertices=np.array([center,center+v,center+v+w,center+w,center])
    answer=0.;imaginary=0.
    for left,right in zip(vertices,vertices[1:]):
        delta=right-left;points=left[None,:]+us[:,None]*delta[None,:]
        e=real_trig_basis(points,N)
        interp=np.einsum('ijk,bi,bj,bk->b',coeff,*e,optimize=True)
        assert np.min(interp.real)>0
        imaginary=max(imaginary,float(np.max(abs(interp.imag))))
        answer+=float(np.dot(weights,interp.real)*np.linalg.norm(delta))
    assert imaginary<1e-12
    return answer,imaginary

def run():
    old_contact=contact.run()
    assert old_contact==json.loads(contact.TARGET.read_text('utf-8'))
    c=inherited.constants();CF=Q(4,3);CA=Q(3)
    assert 2*CF-CA==Q(-1,3)
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        fundamental=sum(T@T for T in old.T)
        ad_generators=np.array([[[2*np.trace(Tb@(Ta@Tc-Tc@Ta)) for Tc in old.T] for Tb in old.T] for Ta in old.T])
        adjoint_casimir=sum(T@T for T in ad_generators)
        casimir_error=max(float(np.max(abs(fundamental-float(CF)*np.eye(3)))),
                          float(np.max(abs(adjoint_casimir-float(CA)*np.eye(8)))))
        assert casimir_error<1e-14
        geo=old.geo;N=24;q,psi0,tensor0,_,_=old.completed(N,1.)
        x,y,z=np.moveaxis(q['grid'],-1,0);epsilon=.02;probe=epsilon*np.sin(x)
        new=dict(q);new['B']=q['B']+epsilon**2*np.cos(x)**2
        new['U']=q['U']+.5*probe**2;new['C']=q['C']-probe**2
        psi,tensor,stats=geo.solve_hamiltonian(new,initial=psi0)
        assert stats['original_Hamiltonian_residual']<3e-8 and np.max(abs(tensor-tensor0))<1e-14
        interpolation_grid_error=interpolation_check(psi)
        rows=[]
        for L in (.4,.2,.1):
            U,dU,ad,dad=matrices(old,c,L)
            chi=np.trace(U);dchi=np.trace(dU);chi8=float(np.trace(ad).real);dchi8=float(np.trace(dad).real)
            error=abs(abs(chi)**2-(1+chi8))
            derr=abs(2*np.real(np.conj(chi)*dchi)-dchi8)
            global_errors=[]
            for theta in (0.,.37,1.14):
                chiR=np.exp(-2j*theta)*chi
                global_errors.append(abs(chiR*np.conj(chiR)-(1+chi8)))
            first=2*float(CF)*abs(chi)**2-float(CA)*chi8
            expected=(8-chi8)/3
            tangent=2*float(CF)*2*float(np.real(np.conj(chi)*dchi))-float(CA)*dchi8
            assert expected>0 and tangent>0
            assert error<2e-14 and derr<2e-14 and max(global_errors)<2e-14
            assert abs(first-expected)<3e-14 and abs(tangent+dchi8/3)<3e-14
            length32,imag32=proper_length(psi,L,32);length64,imag64=proper_length(psi,L,64)
            assert abs(length32-length64)<1e-12
            rows.append(dict(coordinate_side=L,original_proper_perimeter=length64,
                length_quadrature_difference=abs(length32-length64),
                real_interpolation_imaginary_residual=max(imag32,imag64),
                fundamental_trace_real=float(chi.real),fundamental_trace_imag=float(chi.imag),
                adjoint_trace=chi8,bare_fusion_residual=error,quotient_U1_phase_fusion_maximum=float(max(global_errors)),
                constrained_family_adjoint_derivative=dchi8,adjoint_derivative_crosscheck=derr,
                first_finite_normalization_fusion_defect_per_alpha=expected,
                defect_derivative_along_original_solution_per_alpha=tangent,
                existing_probe_mass=1.,actual_perimeter_weighted_defect=length64*expected,
                actual_perimeter_weighted_physical_derivative=length64*tangent))
    alpha=.01
    bad=3-3*np.exp(float(CF)*alpha);assert bad<0
    return dict(round=865,date='2026-10-06',formal_reports=865,fresh_numbered_groups=1,
        cumulative_numbered_groups=3650,all_checks_passed=True,
        representation='R=(3,1)_{-2}, conjugate Rbar; R tensor Rbar = 1 + (8,1)_0',
        CF=str(CF),CA=str(CA),original_Casimir_error=casimir_error,
        even_grid_Nyquist='real cosine splitting; exact at original grid nodes',
        interpolation_grid_error=interpolation_grid_error,original_background_N=N,original_background_probe_amplitude=epsilon,
        original_Hamiltonian_residual=stats['original_Hamiltonian_residual'],rows=rows,
        positive_observable_counterexample=dict(alpha=alpha,
            input='3 - Re chi_R >= 0 on original compact quotient group',
            rescaled_value_at_identity=float(bad),unital_positive_map_possible=False),
        normalized_anchor_products='same-anchor composite, not product of separately averaged characters',
        finite_normalization_candidate_not_computed_counterterm=True,
        shared_source_partner_required_if_preserving_fusion=True,
        nontriviality_on_original_constrained_solutions='analytic small-loop proof and original matrix/metric calibration',
        all_order_interacting_relational_loop_defined=False,
        full_graph_state_dynamics_matching_proved=False,full_goal_completed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
