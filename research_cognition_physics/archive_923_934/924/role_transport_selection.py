"""Round 924: active encoding contracts, transport selection and composition scope.
Finite matrices check analytic claims; no microscopic continuum or gauge dynamics input.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
HERE=Path(__file__).resolve().parent
I2=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)
PAULI=[I2,X,Y,Z]
def op(a):return float(np.linalg.norm(a,2))
def comm(a,b):return a@b-b@a
def hermbasis(d):
    out=[]
    for i in range(d):
        a=np.zeros((d,d),complex);a[i,i]=1;out.append(a)
    for i in range(d):
        for j in range(i+1,d):
            a=np.zeros((d,d),complex);a[i,j]=a[j,i]=1/np.sqrt(2);out.append(a)
            a=np.zeros((d,d),complex);a[i,j]=1j/np.sqrt(2);a[j,i]=-1j/np.sqrt(2);out.append(a)
    return out

def commutant_dimension(gens):
    d=gens[0].shape[0]
    cols=[]
    for a in hermbasis(d):
        q=np.concatenate([comm(a,g).reshape(-1) for g in gens])
        cols.append(np.r_[q.real,q.imag])
    vals=np.linalg.svd(np.column_stack(cols),compute_uv=False)
    rank=int(np.count_nonzero(vals>1e-10))
    return d*d-rank

def weyl(r):
    shift=np.roll(np.eye(r,dtype=complex),1,axis=0)
    phase=np.diag(np.exp(2j*np.pi*np.arange(r)/r))
    return shift,phase

def swaps(r):
    out=[]
    for j in range(r-1):
        u=np.eye(r,dtype=complex);u[[j,j+1]]=u[[j+1,j]];out.append(u)
    return out

def evolve(h,t):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T

def random_hermitian(d,rng):
    a=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
    return (a+a.conj().T)/2

def run():
    rng=np.random.default_rng(924)
    classified=[];twirl=[]
    for r in (2,3,4,5):
        shift,phase=weyl(r)
        dim_perm=commutant_dimension(swaps(r))
        dim_shift=commutant_dimension([shift])
        dim_full=commutant_dimension([shift,phase])
        dim_transport=commutant_dimension([np.kron(I2,g) for g in (shift,phase)])
        assert (dim_perm,dim_shift,dim_full,dim_transport)==(2,r,1,4)
        classified.append(dict(r=r,label_permutation_commutant=dim_perm,cyclic_only_commutant=dim_shift,
            two_encoding_commutant=dim_full,transport_real_parameters_before=4*r*r,
            transport_real_parameters_after=dim_transport))
        v=random_hermitian(2*r,rng)
        v=v/op(v)
        eff=np.kron(np.trace(v.reshape(2,r,2,r),axis1=1,axis2=3)/r,np.eye(r))
        avg=np.zeros_like(v)
        for a in range(r):
            for b in range(r):
                u=np.kron(I2,np.linalg.matrix_power(shift,a)@np.linalg.matrix_power(phase,b))
                avg+=u@v@u.conj().T/(r*r)
        residuals=[op(np.kron(I2,g)@v@np.kron(I2,g).conj().T-v) for g in (shift,phase)]
        error=op(v-eff);bound=(r-1)*sum(residuals)/2
        assert op(avg-eff)<2e-14 and error<=bound+2e-14
        twirl.append(dict(r=r,twirl_identity_error=op(avg-eff),actual_distance=error,
            encoding_residuals=residuals,analytic_distance_bound=bound))
    r=3;s=np.ones(r)/np.sqrt(r);q=np.array([1.,-1.,0])/np.sqrt(2)
    ps=np.outer(s,s);c=np.eye(r)+ps
    vz=np.kron(Z,c);eplus=np.array([1.,0.])
    velocities=[float((np.kron(eplus,a).conj()@vz@np.kron(eplus,a)).real) for a in (s,q)]
    classical=list(np.diag(c).real)
    perm_error=max(op(comm(c,u)) for u in swaps(r))
    assert perm_error==0 and np.allclose(velocities,[2,1]) and np.allclose(classical,[4/3]*3)
    u=np.kron(I2,weyl(r)[1]);rho=np.outer(np.kron(eplus,s),np.kron(eplus,s))
    passive=abs(np.trace((u@vz@u.conj().T)@(u@rho@u.conj().T))-np.trace(vz@rho))
    active=float(abs(np.trace(vz@(u@rho@u.conj().T))-np.trace(vz@rho)))
    assert passive<1e-14 and abs(active-1)<1e-14
    # Same geometry for two inequivalent allowed operation sets; no inference of a physical gauge group.
    su2=[X,Y,Z];finite=[X,Z]
    assert commutant_dimension(su2)==commutant_dimension(finite)==1
    # Requiring common U tensor U on a pair is weaker than independent replacements.
    s4=np.zeros((4,4),complex)
    for a in range(2):
        for b in range(2):s4[2*b+a,2*a+b]=1
    plus=(np.eye(4)+s4)/2;minus=(np.eye(4)-s4)/2
    collective=[np.kron(g,I2)+np.kron(I2,g) for g in (X,Y,Z)]
    independent=[np.kron(g,I2) for g in (X,Z)]+[np.kron(I2,g) for g in (X,Z)]
    c_pair=plus+2*minus
    singlet=np.array([0,1,-1,0],complex)/np.sqrt(2)
    triplet=np.array([1,0,0,0],complex)
    pairvel=[float(np.vdot(a,c_pair@a).real) for a in (triplet,singlet)]
    collective_error=max(op(comm(c_pair,g)) for g in collective)
    separate_error=max(op(comm(c_pair,g)) for g in independent)
    assert commutant_dimension(collective)==2 and commutant_dimension(independent)==1
    assert collective_error==0 and separate_error>.9 and np.allclose(pairvel,[1,2])
    assert op(comm(c_pair,s4))==0
    # Full Hamiltonian invariance would additionally remove internal mixing; principal-only does not.
    velocities_m=[np.kron(g,np.eye(2)) for g in (X,Y,Z)]
    mixing=np.kron(I2,Z)
    lower_order_error=op(comm(mixing,np.kron(I2,X)))
    assert max(op(comm(v,np.kron(I2,X))) for v in velocities_m)==0 and lower_order_error==2
    symbol_commutant=commutant_dimension(velocities_m)
    full_commutant=commutant_dimension(velocities_m+[mixing])
    assert symbol_commutant==4 and full_commutant==2
    # Finite time, finite momentum comparison for the SAME lower-order term and arbitrary old reference.
    r=2;eps=0.03
    vs=[velocities_m[0]+eps*np.kron(Z,X),velocities_m[1]+eps*np.kron(X,Z),velocities_m[2]]
    effs=[];etas=[]
    for v in vs:
        ve=np.kron(np.trace(v.reshape(2,r,2,r),axis1=1,axis2=3)/r,I2);effs.append(ve)
        eta=sum(op(np.kron(I2,g)@v@np.kron(I2,g).conj().T-v) for g in (X,Z))/2
        assert op(v-ve)<=eta+1e-14
        etas.append(eta)
    p=np.array([.2,-.4,.3]);t=.7
    h=mixing+sum(p[i]*vs[i] for i in range(3));he=mixing+sum(p[i]*effs[i] for i in range(3))
    unitary_error=op(evolve(h,t)-evolve(he,t));prob_bound=t*float(np.abs(p)@etas)
    assert unitary_error<=prob_bound+1e-14
    psi=rng.normal(size=8)+1j*rng.normal(size=8);psi/=np.linalg.norm(psi)
    ua=np.kron(evolve(h,t),I2)@psi;ub=np.kron(evolve(he,t),I2)@psi
    rho_a=np.outer(ua,ua.conj());rho_b=np.outer(ub,ub.conj())
    td=float(np.linalg.norm(rho_a-rho_b,ord='nuc')/2)
    assert td<=prob_bound+1e-14
    return dict(round=924,date='2026-10-06',all_scientific_checks_passed=True,
        classification=classified,weyl_twirl=twirl,
        permutation_counterexample=dict(all_classical_label_velocities=classical,coherent_velocities=velocities,
            permutation_commutator_error=perm_error,passive_covariance_error=float(passive),active_phase_change=active),
        composite_counterexample=dict(collective_commutant_dimension=2,independent_commutant_dimension=1,
            symmetric_antisymmetric_velocities=pairvel,collective_commutator_error=collective_error,
            independent_replacement_violation=separate_error,whole_swap_commutator_error=op(comm(c_pair,s4))),
        finite_operations_and_continuous_group_have_same_commutant=True,
        principal_only_allows_lower_order_mixing=True,full_invariance_mixing_violation=lower_order_error,
        identical_symbol_different_full_symmetry_commutant_dimensions=[symbol_commutant,full_commutant],
        finite_scope=dict(velocity_error_bounds=etas,momentum=p.tolist(),time=t,
            unitary_distance=unitary_error,reference_entangled_trace_distance=td,all_effect_probability_bound=prob_bound,
            same_lower_order_term_in_both_models=True,source_or_backreaction_budget_certified=False),
        active_contract_is_added_hypothesis=True,internal_factorization_is_added_input=True,
        irreversible_record_compression_used=False,physical_gauge_group_selected=False,
        gauge_dynamics_generated=False,mass_parameters_selected=False,common_geometry_for_all_matter_proved=False,
        closed_finite_claim='Encoding replacement selects a scalar internal transport coefficient iff its action is irreducible; synchronized composition and classical labels need not suffice.',
        whole_stage_completed=False,full_goal_completed=False,
        source_hashes={str(Path('924')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    data=run();target=HERE/'role_transport_selection_results.json'
    if args.write:
        assert not target.exists();target.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in data.items() if k not in ('weyl_twirl','source_hashes')},ensure_ascii=False,indent=2))
