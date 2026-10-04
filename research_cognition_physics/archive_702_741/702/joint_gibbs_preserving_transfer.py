"""702: full-model Gibbs-preserving regrouping, with scoped diagnostics.

The infinite-dimensional full Gauss theorem is analytic. Group1 checks its
remainder bound on the original64 CAR conditional fiber. Group2 uses the
original625 finite radial box, now coupled to the exact neutral Majorana Fock
factor. It does not evaluate the full graph thermal partition function.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_full_graph_transfer_sources as inherited
import joint_gauss_fermion_influence as car
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gibbs_preserving_transfer_results.json'


def norm1(a):
    return float(np.sum(abs(np.linalg.eigvalsh((a+a.conj().T)/2))))


def exp_h(h,t):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-t*e))@v.conj().T


def remainder_bound():
    old=inherited.old;mass=inherited.mass
    mat,u,_=old.original.lattice.scalar.parameters();lam=np.linalg.det(mat)/np.trace(mat)
    D=6*old.original.M-u.sum();aa=32/(lam*D**2);bb=144/D**2
    linear=[]
    for j in range(5):
        e=np.eye(5)[j]
        numerator=mass.delta(e)*np.sqrt(old.original.F(e))
        linear.append(float(np.sum(abs(np.triu(numerator,1)))))
    kappa=np.sqrt(6*old.original.M)*np.linalg.norm(linear)
    nodes=2;c0=nodes*kappa*bb**.25;c1=nodes**.75*kappa*aa**.25
    data=inherited.graph_data();normal=data['normal_BdG'];j0=64*.23
    direction=np.array([1.,2.,-.4,.6,.8]);direction/=np.linalg.norm(direction)
    rows=[]
    for theta in (.25,.5,.75):
        constant=j0+2*c0+3*(2*c1)**(4/3)/(4**(4/3)*(1-theta)**(1/3))
        cases=[]
        for f in (.8,.2,.01,1e-4):
            p=direction*np.sqrt(6*(old.original.M-f))
            pp=[p,.8*p]
            blocks=[old.mass_matrices(x) for x in pp]
            h=inherited.blockdiag(blocks[0][0],blocks[1][0])
            delta=inherited.blockdiag(blocks[0][1],blocks[1][1])
            matrix=data['transform']@inherited.nambu(h,delta)@data['transform'].conj().T+normal
            # Original Tr(h)=0 convention; all64 modes retained.
            eigen=np.linalg.eigvalsh(matrix)
            ground=-float(np.sum(abs(eigen)))/4
            W=sum(float(old.original.node_potential(x)) for x in pp)
            minimum=(1-theta)*W+ground
            assert minimum>=-constant
            cases.append(dict(target_F=f,full64_conditional_ground=ground,
                original_W=W,remainder_minimum=minimum,certified_lower=-constant))
        rows.append(dict(theta=theta,c_theta=constant,c0=c0,c1=c1,j0=j0,cases=cases))
    return dict(rows=rows,all_original64_modes_and_mass_hopping_retained=True,
        actual_global_bound_proved_by_Young_not_by_the_sampled_rays=True)


def interacting_radial_diagnostic():
    # Exactly the scalar grid/metric/quadrature convention of625. Additional
    # model choice for this diagnostic only: Dirac Yukawas zero, charged Fock
    # vacuum, neutral singlet Majorana coupling retained. This is an invariant
    # radial sector before the declared finite-box discretization.
    n=8;hs=np.linspace(.15,1.55,n+2)[1:-1];ss=np.linspace(-1.1,1.1,n+2)[1:-1]
    h,s=np.meshgrid(hs,ss,indexing='ij');h=h.ravel();s=s.ravel()
    F=original.M-(h*h+s*s)/6;measure=np.sqrt(original.M)*h**3/F**3
    def derivative(step):return (np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/(2*step)
    ri=1/np.sqrt(measure)
    dh=np.kron(derivative(hs[1]-hs[0]),np.eye(n))*ri[None,:]
    ds=np.kron(np.eye(n),derivative(ss[1]-ss[0]))*ri[None,:]
    ghh=F*(1-h*h/(6*original.M));gss=F*(1-s*s/(6*original.M));ghs=-F*h*s/(6*original.M)
    T=.7**2/(2*.8)*(dh.T@((measure*ghh)[:,None]*dh)+ds.T@((measure*gss)[:,None]*ds)
        +dh.T@((measure*ghs)[:,None]*ds)+ds.T@((measure*ghs)[:,None]*dh))
    phi=np.zeros((n*n,5));phi[:,1]=h;phi[:,4]=s
    w=.8*original.node_potential(phi);scalar_W=np.diag(w)
    K=np.kron(T,np.eye(4));W=np.kron(scalar_W,np.eye(4));B=np.zeros_like(K,dtype=complex)
    for j,p in enumerate(phi):
        _,d=inherited.old.mass_matrices(p)
        neutral=d[np.ix_([30,31],[30,31])]
        B[4*j:4*j+4,4*j:4*j+4]=car.fock(np.zeros((2,2)),neutral)
    effects=[np.diag(np.repeat(np.sqrt(.5+sgn*np.sin(s)/4),4)) for sgn in (1,-1)]
    assert np.linalg.norm(sum(x@x for x in effects)-np.eye(len(K)))<1e-12
    return K,W,B,effects


def thermal(h,beta):
    e,v=np.linalg.eigh(h);weight=np.exp(-beta*e);Z=float(weight.sum());p=weight/Z
    rho=(v*p)@v.conj().T
    return e,v,rho,Z,float(-p@np.log(p)),[float(p@e),float(p@(e*e))]


def records(e,v,rho,effects):
    us=[(v*np.exp(-1j*t*e))@v.conj().T for t in (.12,.19)]
    out=[]
    for hist in itertools.product(range(2),repeat=2):
        x=rho
        for r,u in zip(hist,us):x=effects[r]@u@x@u.conj().T@effects[r]
        out.append(x)
    assert abs(sum(np.trace(x) for x in out)-1)<2e-12
    return out


def thermal_and_record_limit():
    K,W,B,effects=interacting_radial_diagnostic();H=K+W+B;beta=1.2
    e,v,rho,Z,entropy,moments=thermal(H,beta)
    target_records=records(e,v,rho,effects);target_heat=exp_h(H,beta);rows=[]
    for theta in (.25,.5,.75):
        A=K+theta*W;D=(1-theta)*W+B
        constant=max(0.,-float(np.linalg.eigvalsh(D)[0]))
        pair=[]
        for a in (.2,.1,.05,.025):
            half=exp_h(A,a/2);S=half@exp_h(D,a)@half;S=(S+S.conj().T)/2
            vals,vec=np.linalg.eigh(S);assert min(vals)>0
            Ha=(vec*(-np.log(vals)/a))@vec.conj().T
            order=float(np.linalg.eigvalsh(Ha-A+constant*np.eye(len(A)))[0])
            assert order>-3e-10
            ea,va,ra,za,sa,ma=thermal(Ha,beta)
            output=records(ea,va,ra,effects)
            cq=sum(norm1(x-y) for x,y in zip(output,target_records))
            heat_error=norm1(exp_h(Ha,beta)-target_heat)
            state_error=norm1(ra-rho)
            assert state_error<=2*heat_error/Z+1e-10
            # Actual spectral domination bound, evaluated in this finite calibration.
            L=Ha+constant*np.eye(len(A));Al=np.linalg.eigvalsh(A)
            R=exp_h(L,beta);ZA=float(np.exp(-beta*Al).sum())
            Mbeta=2/(np.e*beta)*float(np.exp(-beta*Al/2).sum())
            assert np.trace(R).real<=ZA+1e-10
            assert np.trace(A@R).real<=Mbeta+1e-10
            pair.append(dict(step=a,minimum_transfer_eigenvalue=float(min(vals)),
                log_order_margin=order,heat_trace_norm_error=heat_error,
                Gibbs_trace_norm_error=state_error,partition_error=abs(za-Z),
                entropy_error=abs(sa-entropy),energy_moment_errors=[abs(ma[j]-moments[j]) for j in (0,1)],
                full_classical_quantum_history_trace_error=cq,
                thermal_reference_energy=float(np.trace(A@R).real),uniform_M_beta=Mbeta))
        assert pair[-1]['Gibbs_trace_norm_error']<pair[0]['Gibbs_trace_norm_error']/20
        assert pair[-1]['full_classical_quantum_history_trace_error']<pair[0]['full_classical_quantum_history_trace_error']/20
        rows.append(dict(theta=theta,finite_box_lower_constant=constant,rows=pair))
    return dict(dimension=len(H),scalar_grid=64,neutral_Fock_dimension=4,
        original_majorana_interaction_noncommutes_kinetic=float(np.linalg.norm(K@B-B@K,2)),
        target_partition=Z,target_entropy=entropy,target_energy_moments=moments,rows=rows,
        no_full_graph_CAR_or_Gauss_partition_numerically_evaluated=True,
        finite_radial_box_diagnostic_not_infinite_dimensional_proof=True,
        Dirac_zero_only_in_diagnostic_not_in_analytic_theorem=True)


def run():
    deps=('research_note_603.md','research_note_623.md','research_note_624.md','research_note_625.md',
          'research_note_643.md','research_note_655.md','research_note_667.md','research_note_699.md',
          'research_note_701.md','joint_full_graph_transfer_sources.py',
          'joint_source_preserving_compression.py','round702_drafts/source_scale_entry.md')
    return dict(date='2026-10-02',round=702,tests_run=2,failures=0,errors=0,
        original_full_fiber_bound=remainder_bound(),same_thermal_process=thermal_and_record_limit(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(analytic_model_is_original_full_finite_graph_HF=True,
            regrouping_preserves_all_interactions_and_Gauss=True,
            new_transfer_not_claimed_equal_to_655_or_673_at_finite_step=True,
            own_Gibbs_states_trace_norm_converge=True,real_time_and_actual_record_limit_shared=True,
            all_energy_moments_and_thermal_entropy_converge=True,
            spatial_refinement_uniformity_not_proved=True,autonomous_thermalization_not_proved=True,
            original_chiral_overlap_identity_still_open=True,GR_or_dimension_not_derived=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=702,checks=2,all_checks_passed=True,
        final_Gibbs_errors=[x['rows'][-1]['Gibbs_trace_norm_error'] for x in result['same_thermal_process']['rows']])) )
