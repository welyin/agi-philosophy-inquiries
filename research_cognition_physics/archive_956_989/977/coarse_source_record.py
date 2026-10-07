"""977: a finite classical source label with quantum records inside each bin.
This is a task-specific representation, not microscopic physical dephasing.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'coarse_source_record_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def kron(*xs):
    out=np.ones((1,1))
    for x in xs:out=np.kron(out,x)
    return out
def norm(x):return float(np.linalg.norm(x,2))
def trace_norm(x):return float(sum(abs(np.linalg.eigvalsh((x+x.conj().T)/2))))

def run():
    native=load('native968',STAGE/'968/internal_relay.py')
    m,h,q,f,W,identities=native.material()
    phases=load('phase965',STAGE/'965/material_field_window.py').phases
    core976=load('core976',STAGE/'976/finite_internal_reference.py')
    previous=read(STAGE/'976/finite_internal_reference_results.json')
    T=previous['parameters']['T'];t0=previous['parameters']['t0'];J=m['J']
    phi=np.array(previous['parameters']['Phi']);lapse=1+phi
    mean_lapse=float(np.mean(lapse));half_difference=float((phi[0]-phi[1])/2)
    interval=np.array(previous['parameters']['coordinate_interval']);width=.01
    n=11;I=np.eye(4);ann=np.diag(np.sqrt(np.arange(1,n)),1)
    Hm=kron(h,I)+kron(I,h)+.2*kron(q,q);S=kron(q,I)+kron(I,q)
    H=kron(Hm,np.eye(n))+.5*kron(np.eye(16),np.diag(np.arange(n)))
    H+=.003*kron(S,ann+ann.T)+.000018*kron(S@S,np.eye(n))
    # Preserve the EXACT conserved sender P_A sectors in the finite model,
    # including degenerate eigenvalues. This makes the reference-uniform
    # output statement structural, rather than inferred from two test states.
    split=3*4*n
    assert np.array_equal(H[:split,split:],np.zeros((split,len(H)-split)))
    va,ua=np.linalg.eigh(H[:split,:split]);vb,ub=np.linalg.eigh(H[split:,split:])
    values=np.concatenate([va,vb]);V=np.zeros_like(H)
    V[:split,:split]=ua;V[split:,split:]=ub;dim=len(values)
    cert=core976.certificate(H,Hm,S,n,values,V,kron(W,W,np.eye(n)[:,:2]),float(lapse.max()*interval[1]))
    assert cert['total_state_vector_error']<3e-8
    # The finite candidate is defined by this real spectral approximation.
    # Its exact polar eigenbasis exists; norm defect is accounted below.
    gram=norm(V.T@V-np.eye(dim));assert gram<1e-12
    spectral_residual=norm(H@V-V*values)
    assert spectral_residual<1e-12
    energies=(values[:,None]+np.array([-J,0.])[None,:]).ravel()
    labels=np.rint(energies/width).astype(int);centres=labels.astype(float)*width
    same=labels[:,None]==labels[None,:]
    edge_distance=float(min(width/2-abs(energies-centres)))
    assert edge_distance>1e-7
    # Compress the actual B/R SWAP effect to R=span(g,t); no arbitrary
    # scalar receiver effect is substituted. Ordering: A, B, photon, R.
    E=.5*np.eye(2*dim)
    for a in range(2):
        for b in range(2):
            e=np.zeros((2,2));e[a,b]=1
            E+=.5*kron(I,np.outer(W[:,b],W[:,a]),np.eye(n),e)
    Vbig=kron(V,np.eye(2));E=Vbig.T@E@Vbig
    assert min(np.linalg.eigvalsh(E))>-1e-12 and max(np.linalg.eigvalsh(E))<1+1e-12
    Eb=E*same
    photon=np.zeros(n);photon[:2]=1/math.sqrt(2)
    plus=W@np.ones(2)/math.sqrt(2)
    initial=np.column_stack([np.kron(np.kron(W[:,s],plus),photon) for s in (0,1)])
    coefficients=np.repeat((V.T@initial)/math.sqrt(2),2,axis=0)
    assert np.count_nonzero(coefficients[2*split:,0])==0
    assert np.count_nonzero(coefficients[:2*split,1])==0
    assert np.count_nonzero(E[:2*split,2*split:])==0
    # Global source and interaction energy are retained in the bin quantum H.
    mass=200+energies;classical_mass=200+centres
    Hcq=np.column_stack([mean_lapse*mass+sign*half_difference*classical_mass for sign in (1,-1)])
    assert float(Hcq.min())>199.67
    generator_error=float(half_difference*max(abs(energies-centres)))
    source_derivative_error=float(max(abs(energies-centres))/2)
    centre_state_error=float(interval[1]*half_difference*width/2)
    numerical_budget=1e-5
    assert 4*cert['total_state_vector_error']+8*gram<1e-6
    def output(c,time,effect,coarse=False):
        frequencies=np.column_stack([mean_lapse*energies+sign*half_difference*centres for sign in (1,-1)]) if coarse else energies[:,None]*lapse
        states=c[:,None]*np.column_stack([phases(frequencies[:,j],float(time)) for j in (0,1)])
        # Known 200*path rest-mass phase is omitted on BOTH sides, a common
        # output unitary; trace distances and all joint comparisons are intact.
        total=(states.conj().T@states).T/2
        one=(states.conj().T@effect@states).T/2
        return [total-one,one]
    rows=[]
    for s in (0,1):
        c=coefficients[:,s]
        # Uniform in time, no oscillatory grid is used to prove this bound.
        dropped_weight=float(np.sum(abs(np.outer(c.conj(),c)*(E-Eb))))
        deletion_bound=2*dropped_weight
        total_bound=deletion_bound+centre_state_error+numerical_budget
        assert total_bound<.025
        data=[]
        for time in (float(interval[0]),t0,float(interval[1])):
            exact=output(c,time,E);pinched=output(c,time,Eb);cq=output(c,time,Eb,True)
            dephased=output(c,time,np.diag(np.diag(E)))
            diff=lambda a,b:sum(trace_norm(x-y) for x,y in zip(a,b))/2
            ddrop=diff(exact,pinched);dcentre=diff(pinched,cq);dtotal=diff(exact,cq)
            assert ddrop<=deletion_bound+1e-9
            assert dcentre<=centre_state_error+1e-9
            assert dtotal<total_bound
            for state in (exact,pinched,cq,dephased):
                assert abs(sum(np.trace(x).real for x in state)-1)<1e-11
                assert min(float(np.linalg.eigvalsh(x).min()) for x in state)>-1e-12
            data.append(dict(time=time,exact_record_probability=float(np.trace(exact[1]).real),
                cq_record_probability=float(np.trace(cq[1]).real),
                rank_one_energy_diagonal_record_probability=float(np.trace(dephased[1]).real),
                actual_pinching_joint_distance=ddrop,actual_centre_joint_distance=dcentre,
                actual_total_joint_distance=dtotal))
        weights=abs(c)**2
        bins=[]
        for k in np.unique(labels):
            prob=float(weights[labels==k].sum())
            if prob>1e-12:bins.append(dict(label=int(k),centre_mass=200+float(k)*width,probability=prob))
        mean_before=float(weights@mass);mean_coarse=float(weights@classical_mass)
        assert abs(mean_before-mean_coarse)<=width/2+1e-11
        # Removing inter-bin coherences keeps every function of total energy
        # exactly. Unlike pure classicalization, quantum blocks stay intact.
        rows.append(dict(label=s,dropped_effect_weight=dropped_weight,
            pinching_menu_uniform_distance_bound=deletion_bound,
            full_interval_joint_distance_bound=total_bound,
            exact_internal_mass_mean=mean_before,coarse_mass_mean=mean_coarse,
            same_energy_mean_and_distribution_under_pinching=True,
            selected_bin_probabilities=bins,endpoints=data))
    maximum=max(r['full_interval_joint_distance_bound'] for r in rows)
    certified_record=previous['interval_certificate']['uniform_infinite_contrast_lower']-2*maximum
    assert certified_record>.43
    # Sender P_A is conserved; its off-diagonal code matrix units disappear
    # on the A trace for this menu, before and after the bin map. Hence the
    # per-label maximum also covers an arbitrary sender and passive reference
    # on this specified OUTPUT, not preservation of all internal coherences.
    PA=kron(np.diag([1.,1.,1.,0.]),I,np.eye(n))
    pa=Vbig.T@kron(PA,np.eye(2))@Vbig
    commutator=norm(pa@E-E@pa)
    crossbin_pa=norm(pa*(~same))
    assert commutator<1e-11 and crossbin_pa<1e-9
    differences=[rows[1]['endpoints'][j]['cq_record_probability']-rows[0]['endpoints'][j]['cq_record_probability'] for j in range(3)]
    diagonal_differences=[abs(rows[1]['endpoints'][j]['rank_one_energy_diagonal_record_probability']-rows[0]['endpoints'][j]['rank_one_energy_diagonal_record_probability']) for j in range(3)]
    # Conditional distribution is not a single mean source: actual source
    # visibility remains roughly .7045 whereas a single mean gives unit modulus.
    characteristic=[]
    for s in (0,1):
        weights=abs(coefficients[:,s])**2
        chi=complex(weights@phases(energies,math.pi))
        reduced=complex(weights@phases(centres,math.pi))
        error=abs(chi-reduced)
        assert error<=math.pi*width/2+1e-10
        characteristic.append(dict(label=s,exact_visibility=abs(chi),coarse_visibility=abs(reduced),
            characteristic_error=error,mean_source_visibility=1.))
    files=[Path(__file__),HERE/'drafts/coarse_source_decision.md',
        STAGE/'956/native_material_interface.py',STAGE/'965/material_field_window.py',
        STAGE/'968/internal_relay.py',STAGE/'976/finite_internal_reference.py',
        STAGE/'976/finite_internal_reference_results.json',STAGE/'973/common_source_distribution_results.json']
    return dict(round=977,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(bin_width=width,time_interval=interval.tolist(),material_dimension=176,
            joint_internal_dimension=352,nonempty_spectral_bins=len(set(labels)),
            minimal_finite_spectrum_distance_to_bin_edge=edge_distance),
        algebra=dict(native_isometry_identities=identities,eigenbasis_gram_defect=gram,
            spectral_residual=spectral_residual,source_sector_effect_commutator=commutator,
            source_sector_off_bin_residual=crossbin_pa,minimum_cq_generator=float(Hcq.min()),
            source_derivative_error=source_derivative_error,
            full_energy_commutes_with_cq_generator_exactly=True,
            finite_spectral_polar_basis_defines_positive_candidate=True),
        uniform_budget=dict(centre_state_error=centre_state_error,
            finite_generator_difference=generator_error,numerical_and_infinite_budget=numerical_budget,
            infinite_occupation_certificate=cert,maximum_joint_output_distance=maximum,
            record_contrast_lower=certified_record,uniform_on_entire_interval=True),
        rows=rows,cq_endpoint_record_contrasts=differences,
        rank_one_energy_diagonal_endpoint_record_contrasts=diagonal_differences,
        characteristic_functions=characteristic,
        scope=dict(same_preparation_and_hamiltonian_compared=True,
            classical_bin_and_quantum_fibre_adopted=True,
            pinching_claimed_as_physical_irreversibility=False,
            all_internal_quantum_information_preserved=False,
            exact_infinite_H_spectral_projectors_approximated=False,
            arbitrary_sender_reference_on_declared_joint_output=True,
            repeat_measurement_or_future_control_covered=False,
            full_source_replaced_by_single_deterministic_mean=False,
            dynamical_classical_metric_derived=False,
            Einstein_Langevin_noise_kernel_matched=False,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-6,abs_tol=1e-10),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write:assert not OUT.exists()
    data=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(data,read(OUT))
    print(json.dumps({k:v for k,v in data.items() if k not in ('source_hashes','rows')},ensure_ascii=False,indent=2))
