"""975: readiness/correlation accounting on the unchanged 974 joint process.
The product Gibbs reference is a bookkeeping choice, not a prepared bath or
the equilibrium state of the interacting Hamiltonian. No new dynamics is added.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'record_resource_audit_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,p):
    spec=importlib.util.spec_from_file_location(name,p)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def entropy(rho):
    ev=np.linalg.eigvalsh((rho+rho.conj().T)/2)
    assert min(ev)>-1e-12 and abs(sum(ev)-1)<2e-11
    ev=ev[ev>1e-15];return float(-np.sum(ev*np.log(ev)))
def h2(x):
    if x==0 or x==1:return 0.
    return -x*math.log(x)-(1-x)*math.log1p(-x)
def thermal_entropy(n):return math.log1p(n)+n*math.log1p(1/n)
def reduce_pure(tensor,axis):
    v=np.moveaxis(tensor,axis,0).reshape(tensor.shape[axis],-1)
    return v@v.conj().T
def mean(op,block):return float(np.sum(block.conj()*(op@block)).real)
def norm(a):return float(np.linalg.norm(a,2))

def run():
    core=load('material965',STAGE/'965/material_field_window.py')
    prior=read(STAGE/'974/joint_energy_source_results.json')
    m,Q,Hm,S,W0,local=core.material();par=prior['parameters']
    n=par['computational_cutoff']+1;g=par['g'];omega=par['omega'];T=par['T']
    taus=par['proper_times'];pot=par['potentials'];m0=par['m0']
    eps=prior['certificate']['transported_bound'];beta=2*math.log(2)
    ann=np.diag(np.sqrt(np.arange(1,n)),1)
    N=np.diag(np.arange(n,dtype=float))
    Hbare=np.kron(m['H'],np.eye(6))+np.kron(np.eye(6),m['H'])
    Vcap=Hm-Hbare
    F=np.kron(np.eye(36),omega*N)
    C=g*np.kron(S,ann+ann.T);D=g*g/omega*np.kron(S@S,np.eye(n))
    Vlocal=np.kron(Vcap,np.eye(n))+C+D
    H=np.kron(Hm,np.eye(n))+F+C+D
    val,V=np.linalg.eigh(H)
    us=[(V*core.phases(val,t))@V.T for t in taus]
    phi=np.zeros(n,complex);phi[:2]=1/np.sqrt(2)
    p_singlet=np.eye(6)-m['single']@(np.eye(4)-m['ps'])@m['single'].T
    conserved_label=np.kron(p_singlet,np.eye(6))
    assert norm(conserved_label@Hm-Hm@conserved_label)<1e-12
    assert norm(conserved_label@S-S@conserved_label)<1e-12
    source_effect=np.kron(p_singlet,np.eye(4))
    targetB=(np.exp(1j*m['J']*T)*local[:,0]+local[:,1])/np.sqrt(2)
    receiver_effect=np.outer(targetB,targetB.conj())
    pure_in=[];pure_out=[];energy=[]
    for label in (0,1):
        psi=np.kron(np.kron(local[:,label],local@np.array([1.,1.])/np.sqrt(2)),phi)
        block=psi.reshape(6,4,6,4,n).transpose(0,2,4,1,3).reshape(36*n,16)
        inp=psi.reshape(24,24,n)
        pure_in.append(np.stack([inp,inp],axis=-1)/np.sqrt(2))
        branches=[];rows=[]
        for l,u in enumerate(us):
            x=u@block
            tensor=x.reshape(6,6,n,4,4).transpose(0,3,1,4,2).reshape(24,24,n)
            phase=core.phases(np.array([m0]),taus[l]-taus[0])[0]
            branches.append(phase*tensor/np.sqrt(2))
            rows.append(dict(branch=l,full_local_before=mean(H,block),full_local_after=mean(H,x),
                local_interaction_before=mean(Vlocal,block),local_interaction_after=mean(Vlocal,x),
                gravity_before=pot[l]*(m0+mean(H,block)),gravity_after=pot[l]*(m0+mean(H,x))))
        pure_out.append(np.stack(branches,axis=-1));energy.append(rows)
    assert abs(np.vdot(pure_out[0],pure_out[1]))<1e-12
    dims=[24,24,n,2];names=['sender_A','receiver_B','mode_F','path_P']
    local_h=[np.kron(m['H'],np.eye(4)),np.kron(m['H'],np.eye(4)),omega*N,np.zeros((2,2))]
    logzmat=math.log(float(np.sum(np.exp(-beta*np.linalg.eigvalsh(local_h[0])))))
    # The oscillator reference is the INFINITE normalized Gibbs state.
    logz=[logzmat,logzmat,-math.log1p(-math.exp(-beta*omega)),math.log(2)]
    all_rows=[];marginals={}
    for when,states in [('initial',pure_in),('final',pure_out)]:
        cond=[[reduce_pure(v,ax) for ax in range(4)] for v in states]
        avg=[(cond[0][ax]+cond[1][ax])/2 for ax in range(4)]
        marginals[when]=(cond,avg)
        ent=[entropy(r) for r in avg]
        energies=[float(np.trace(h@r).real) for h,r in zip(local_h,avg)]
        relatives=[beta*e+z-s for e,z,s in zip(energies,logz,ent)]
        holevo=[ent[i]-(entropy(cond[0][i])+entropy(cond[1][i]))/2 for i in range(4)]
        tc=sum(ent)-math.log(2)
        all_rows.append(dict(time=when,marginal_entropies=dict(zip(names,ent)),
            bare_energies=dict(zip(names,energies)),
            local_relative_entropies=dict(zip(names,relatives)),
            local_readiness=sum(relatives),total_correlation=tc,
            label_Holevo_information=dict(zip(names,holevo)),global_entropy=math.log(2)))
    ini,fin=all_rows
    dv=sum(r['local_interaction_after']-r['local_interaction_before']+
           r['gravity_after']-r['gravity_before'] for rs in energy for r in rs)/4
    budget_res=ini['local_readiness']-(fin['local_readiness']+fin['total_correlation']+beta*dv)
    assert abs(budget_res)<1e-10 and abs(ini['total_correlation'])<1e-11
    assert abs(fin['label_Holevo_information']['sender_A']-math.log(2))<1e-11
    cfinal=marginals['final'][0]
    label_probs=[float(np.trace(source_effect@cfinal[s][0]).real) for s in (0,1)]
    rec_probs=[float(np.trace(receiver_effect@cfinal[s][1]).real) for s in (0,1)]
    for s in (0,1):assert abs(rec_probs[s]-prior['rows'][s]['receiver_probability_unconditional'])<1e-9
    assert abs(label_probs[0]-1)<1e-11 and abs(label_probs[1])<1e-11
    error=(rec_probs[0]+1-rec_probs[1])/2
    mi_lower=math.log(2)-h2(error+eps)
    # Distinguish the produced phase-sensitive record from the conserved sector
    # pointer that protects the sender. The receiver's sector weight is identical
    # for both labels. Dephasing is a DIAGNOSTIC channel, not added dynamics.
    pb=source_effect; qb=np.eye(24)-pb
    sector_probs=[float(np.trace(pb@cfinal[s][1]).real) for s in (0,1)]
    dephased=[pb@cfinal[s][1]@pb+qb@cfinal[s][1]@qb for s in (0,1)]
    dephased_read=[float(np.trace(receiver_effect@r).real) for r in dephased]
    dephased_distance=float(np.sum(abs(np.linalg.eigvalsh(dephased[0]-dephased[1])))/2)
    assert max(abs(p-.5) for p in sector_probs)<1e-11
    assert dephased_distance+2*eps<.001
    # Energy bound from the SAME complete square, for exact and truncated states.
    largest_E=max(rs[0]['full_local_before'] for rs in energy)
    sqrtN_bound=math.sqrt((largest_E-np.linalg.eigvalsh(Hm)[0])/omega)+g*norm(S)/omega
    assert sqrtN_bound<1
    f24=eps*math.log(23)+h2(eps);f2=h2(eps)
    fosc=2*eps*thermal_entropy(1/eps)+h2(eps)
    entropy_error=2*f24+f2+fosc+1e-9
    c_norm_on_state=g*norm(S)*(1+math.sqrt(2))
    # H_L is conserved and initially exactly represented. Thus H0 differences
    # follow from the interaction: bounded Vcap,D plus energy-controlled C.
    v_error=2*eps*(norm(Vcap)+norm(D)+c_norm_on_state)+1e-10
    readiness_error=entropy_error+beta*v_error
    fin['total_correlation_interval']=[fin['total_correlation']-entropy_error,
                                      fin['total_correlation']+entropy_error]
    fin['local_readiness_interval']=[fin['local_readiness']-readiness_error,
                                    fin['local_readiness']+readiness_error]
    fin['receiver_label_information_error']=2*f24+1e-9
    assert mi_lower>.675 and fin['total_correlation']-entropy_error>mi_lower
    assert fin['local_readiness']+readiness_error<ini['local_readiness']
    # Same global quantum channel preserves all label information. This does
    # not copy arbitrary unknown quantum states; the source is a fixed sector bit.
    files=[Path(__file__),STAGE/'965/material_field_window.py',STAGE/'956/native_material_interface.py',
        STAGE/'974/joint_energy_source.py',STAGE/'974/joint_energy_source_results.json',
        STAGE/'972/native_thermal_record_results.json']
    return dict(round=975,date='2026-10-07',all_scientific_checks_passed=True,
        unchanged_physical_hamiltonian_preparation_and_endpoint=True,
        reference=dict(beta=beta,partitions=names,reference_is_product_of_bare_Gibbs_states=True,
            reference_is_actual_joint_equilibrium=False,physical_heat_bath_added=False),
        states=all_rows,energy_rows=energy,
        resource_balance=dict(initial_readiness=ini['local_readiness'],
            final_local_readiness=fin['local_readiness'],final_total_correlation=fin['total_correlation'],
            interaction_energy_change=dv,beta_interaction_change=beta*dv,
            identity_residual=budget_res,initial_correlations=0.),
        record=dict(source_label_probabilities=label_probs,receiver_probabilities=rec_probs,
            error_probability=error,error_upper=error+eps,
            certified_binary_information_lower=mi_lower,
            quantum_receiver_information=fin['label_Holevo_information']['receiver_B'],
            sender_retains_label_exactly=True,arbitrary_unknown_quantum_state_copied=False,
            all_time_record_retention_proven=False,
            receiver_conserved_sector_probabilities=sector_probs,
            sector_dephasing_is_diagnostic_only=True,
            after_sector_dephasing_same_effect=dephased_read,
            after_sector_dephasing_trace_distance=dephased_distance,
            after_sector_dephasing_infinite_distance_upper=dephased_distance+2*eps),
        error_bounds=dict(state_trace_distance=eps,oscillator_mean_number_bound=1.,
            derived_sqrt_mean_number_bound=sqrtN_bound,finite_material_entropy=f24,
            path_entropy=f2,oscillator_entropy=fosc,total_correlation=entropy_error,
            interaction_energy=v_error,local_readiness=readiness_error,
            infinite_oscillator_entropy_controlled=True),
        scope=dict(resource_identity_is_mature_theorem=True,
            snapshot_nondestructive_classical_record_certified=True,
            global_entropy_production=False,actual_heat_or_work_cost_identified=False,
            continuous_dissipation_required_for_storage=False,
            stable_macroscopic_record_dynamics_completed=False,
            full_cyclic_reset_completed=False,cosmic_low_entropy_origin_derived=False,
            full_goal_completed=False),
        adoption='local readiness is redistributed into classical records and other correlations; no new bath or arrow mechanism is claimed',
        references=['https://arxiv.org/abs/1404.2169','https://arxiv.org/html/1507.07775'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    assert a['source_hashes']==b['source_hashes']
    for key in ('initial_readiness','final_local_readiness','final_total_correlation','interaction_energy_change'):
        assert abs(a['resource_balance'][key]-b['resource_balance'][key])<1e-8,key
    assert abs(a['record']['quantum_receiver_information']-b['record']['quantum_receiver_information'])<1e-8

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    elif TARGET.exists():compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
