"""976: a finite native reference and fixed relational record readout.
Finite-time, same-source certificate; no storage-lifetime or pulse search.
"""
from pathlib import Path
from decimal import Decimal, localcontext
from fractions import Fraction
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent; STAGE=HERE.parent; ROOT=HERE.parents[2]
OUT=HERE/'finite_internal_reference_results.json'
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
def read(path):return json.loads(path.read_text('utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def kron(*args):
    out=np.ones((1,1))
    for a in args:out=np.kron(out,a)
    return out
def gamma(n):
    u=np.finfo(float).eps/2;return n*u/(1-n*u)
def norm(a):return float(np.linalg.norm(a,2))
def fro(a):return float(np.linalg.norm(a))

def construction_error(H,Hm,S,n):
    """Exact rational/radical entries, including off-diagonal S and S^2."""
    size=len(H);errors=np.zeros(size);S2=S@S
    with localcontext() as ctx:
        ctx.prec=70;g=Decimal(3)/1000;w=Decimal(1)/2
        for i in range(16):
            for j in range(16):
                rational=Fraction(float(Hm[i,j])).limit_denominator(100)
                r=Decimal(rational.numerator)/Decimal(rational.denominator)
                assert abs(float(r)-Hm[i,j])<1e-15
                assert S[i,j]==round(S[i,j]) and S2[i,j]==round(S2[i,j])
                for p in range(n):
                    target=r+g*g/w*int(S2[i,j])+(w*p if i==j else 0)
                    errors[i*n+p]+=float(abs(Decimal.from_float(float(H[i*n+p,j*n+p]))-target))+1e-60
                for p in range(n-1):
                    target=g*int(S[i,j])*Decimal(p+1).sqrt()
                    for a,b in ((p,p+1),(p+1,p)):
                        errors[i*n+a]+=float(abs(Decimal.from_float(float(H[i*n+a,j*n+b]))-target))+1e-60
    return float(max(errors))*(1+1e-10)+1e-55

def certificate(H,Hm,S,n,values,V,Vin,time):
    dim=len(H);co=V.T@Vin;roundco=gamma(dim)*(abs(V.T)@abs(Vin))
    upper=np.linalg.norm(co,axis=1)+np.linalg.norm(roundco,axis=1)+1e-13
    errH=construction_error(H,Hm,S,n)
    residual=H@V-V*values
    u=np.finfo(float).eps/2
    rounding=gamma(dim)*(abs(H)@abs(V))+3*u*(abs(H@V)+abs(V*values))
    finite=np.linalg.norm(residual,axis=0)+np.linalg.norm(rounding,axis=0)+errH*np.linalg.norm(V,axis=0)
    boundary=.003*math.sqrt(n)*(S@V.reshape(16,n,dim)[:,-1,:])
    boundary_norm=np.linalg.norm(boundary,axis=0)*(1+gamma(100))+1e-15
    initial=fro(V@co-Vin)+fro(gamma(dim)*(abs(V)@abs(co)))+fro(abs(V)@roundco)+1e-11
    tail=time*float(boundary_norm@upper);spectral=time*float(finite@upper)
    total=initial+tail+spectral+1e-8
    return dict(maximum_proper_time=time,initial_defect_bound=initial,
        infinite_occupation_boundary=tail,finite_eigen_residual=spectral,
        defining_matrix_error=errH,phase_and_definition_guard=1e-8,
        total_state_vector_error=total,input_isometry_defect=fro(Vin.T@Vin-np.eye(Vin.shape[1])))

def run():
    old=load('native968',STAGE/'968/internal_relay.py')
    m,h,q,F,W,identities=old.material()
    phase=load('phase965',STAGE/'965/material_field_window.py').phases
    prior=read(STAGE/'974/joint_energy_source_results.json')
    T=read(STAGE/'958/capacitive_material_write_results.json')['exact_interface']['controlled_phase_time']
    J=m['J'];t0=T+math.pi;delta=math.pi/t0
    Phi=np.array([-delta,-2*delta]);lapse=1+Phi
    interval=np.array([.9*t0,1.1*t0]);n=11;I=np.eye(4)
    Hm=kron(h,I)+kron(I,h)+.2*kron(q,q)
    S=kron(q,I)+kron(I,q)
    ann=np.diag(np.sqrt(np.arange(1,n)),1);N=np.diag(np.arange(n,dtype=float))
    H=kron(Hm,np.eye(n))+.5*kron(np.eye(16),N)+.003*kron(S,ann+ann.T)+.000018*kron(S@S,np.eye(n))
    values,V=np.linalg.eigh(H)
    Vin=kron(W,W,np.eye(n)[:,:2])
    cert=certificate(H,Hm,S,n,values,V,Vin,float(lapse.max()*interval[1]))
    assert cert['total_state_vector_error']<5e-7
    probability_error=2e-6
    assert 2*cert['total_state_vector_error']+cert['total_state_vector_error']**2<probability_error
    photon=np.zeros(n);photon[:2]=1/math.sqrt(2)
    g=W[:,0];t=W[:,1];r0=(g+t)/math.sqrt(2)
    psi=np.column_stack([np.kron(np.kron(W[:,s],r0),photon) for s in (0,1)])
    c=V.T@psi
    O0=kron(I,(np.outer(g,g)+np.outer(t,t))/2,np.eye(n))
    Op=kron(I,np.outer(g,t)/2,np.eye(n))
    eo=V.T@O0@V;ep=V.T@Op@V
    dc=np.outer(c[:,1].conj(),c[:,1])-np.outer(c[:,0].conj(),c[:,0])
    freq=values[:,None]-values[None,:]
    L=float(np.sum(abs(dc*eo)*abs(freq))+2*np.sum(abs(dc*ep)*abs(freq+J)))
    # All-time derivative bound for the finite spectral formula, then transfer
    # to the exact infinite-occupation evolution with a uniform value error.
    Ljoint=float(sum(lapse)*L/4)*(1+1e-10)+1e-10
    intervals=max(256,math.ceil(Ljoint*float(interval[1]-interval[0])/(2*.002)))
    assert intervals<100000,'Fixed physical contract failed available certificate budget; do not change it.'
    grid=np.linspace(*interval,intervals+1)
    def branch_probs(tau,high_precision=False):
        tau=np.atleast_1d(tau)
        z=np.stack([phase(values,float(x)) for x in tau]) if high_precision else np.exp(-1j*tau[:,None]*values)
        states=np.einsum('ij,kj,js->kis',V,z,c,optimize=True).reshape(len(tau),4,4,n,2)
        rp=np.exp(1j*J*tau)
        ref=(rp[:,None]*g+t)/math.sqrt(2)
        overlap=np.einsum('kb,kabns->kans',ref.conj(),states,optimize=True)
        p=np.sum(abs(overlap)**2,axis=(1,2))
        return (1+p)/2
    pg=np.zeros((len(grid),2))
    for start in range(0,len(grid),512):
        sl=slice(start,min(start+512,len(grid)))
        pg[sl]=sum(branch_probs(l*grid[sl]) for l in lapse)/2
    contrast=pg[:,1]-pg[:,0]
    gap=Ljoint*float(grid[1]-grid[0])/2
    lower=float(min(contrast))-gap-2*probability_error-1e-8
    # Independent direct states at three fixed coordinates check spectral sums,
    # high-precision phase evaluation and the finite derivative implementation.
    endpoints=[];formula_errors=[];phase_errors=[]
    for x in (float(interval[0]),t0,float(interval[1])):
        rows=[]
        for l in lapse:
            tau=float(l*x);direct=branch_probs([tau],True)[0]
            fast=branch_probs([tau])[0]
            ef=np.exp(1j*freq*tau)
            f=(np.sum(dc*eo*ef)+2*np.real(np.sum(dc*ep*ef*np.exp(1j*J*tau)))).real
            formula_errors.append(abs(float(direct[1]-direct[0])-float(f/2)))
            phase_errors.append(float(max(abs(direct-fast))))
            rows.append(direct)
        mean=np.mean(rows,axis=0)
        endpoints.append(dict(coordinate_time=x,probabilities=mean.tolist(),contrast=float(mean[1]-mean[0])))
    assert max(formula_errors)<1e-9 and max(phase_errors)<1e-9
    assert lower>.4
    # E_sym is fixed, positive and invariant under common free material motion.
    swap=np.eye(16).reshape(4,4,16).transpose(1,0,2).reshape(16,16)
    Esym=(np.eye(16)+swap)/2;hBR=kron(h,I)+kron(I,h)
    invariant=norm(Esym@hBR-hBR@Esym)
    assert invariant<1e-14
    reference_mean=float(np.vdot(r0,h@r0).real)
    reference_variance=float(np.vdot(r0,h@h@r0).real-reference_mean**2)
    assert abs(reference_mean+J/2)<1e-14 and abs(reference_variance-J**2/4)<1e-14
    # The old matter-field marginal on each branch is unchanged; the path
    # marginal is not, since its characteristic function includes R.
    source_rows=[];deltatau=math.pi
    ref_chi=complex((1+np.exp(1j*J*deltatau))/2)
    for s in (0,1):
        energy=float(np.vdot(psi[:,s],H@psi[:,s]).real)
        variance=float(np.vdot(H@psi[:,s],H@psi[:,s]).real-energy**2)
        chi=complex(np.sum(abs(c[:,s])**2*phase(values,deltatau)))
        chi_joint=chi*ref_chi
        source_rows.append(dict(label=s,local_energy_mean_before_reference=energy,
            mass_mean_with_reference=200+energy+reference_mean,
            mass_variance_with_reference=variance+reference_variance,
            old_path_visibility=abs(chi),new_path_visibility=abs(chi_joint),
            new_internal_characteristic_real=chi_joint.real,new_internal_characteristic_imag=chi_joint.imag))
    # Check actual same-process marginal against 974's endpoint, using its fixed
    # receiver phase instead of the new relational effect.
    oldtarget=(np.exp(1j*J*T)*g+t)/math.sqrt(2)
    fixed=[]
    for l in lapse:
        states=(V@(phase(values,float(l*t0))[:,None]*c)).reshape(4,4,n,2)
        overlap=np.einsum('b,abns->ans',oldtarget.conj(),states)
        fixed.append(np.sum(abs(overlap)**2,axis=(0,1)))
    old_probs=np.mean(fixed,axis=0)
    # Prior schema checked explicitly, not inferred from its number of tests.
    old_expected=np.array([row['receiver_probability_unconditional'] for row in prior['rows']])
    old_error=float(max(abs(old_probs-old_expected)))
    assert old_error<1e-9
    sector_lower=200+float(min(np.linalg.eigvalsh(Hm)))+float(min(np.linalg.eigvalsh(h)))
    # Full 24-dimensional materials, not merely the populated 4D sectors:
    # H_L=H_m+omega(a+g S/omega)^dagger(a+g S/omega), h_i>=-J,
    # ||Q_A Q_B||<=1. This global analytic lower bound proves positivity.
    positivity_lower=200-3*J-.2
    assert sector_lower>199.88 and positivity_lower>199.68
    binary_error=(1-lower)/2
    binary_information=math.log(2)+binary_error*math.log(binary_error)+(1-binary_error)*math.log(1-binary_error)
    files=[Path(__file__),HERE/'drafts/finite_reference_decision.md',
        STAGE/'956/native_material_interface.py',STAGE/'958/capacitive_material_write_results.json',
        STAGE/'965/material_field_window.py',STAGE/'968/internal_relay.py',
        STAGE/'974/joint_energy_source_results.json',STAGE/'975/record_resource_audit_results.json']
    return dict(round=976,date='2026-10-07',all_scientific_checks_passed=True,
        material_reduction_identities=identities,parameters=dict(T=T,t0=t0,J=J,cutoff=10,
            Phi=Phi.tolist(),reference_baseline_mass=100,coordinate_interval=interval.tolist()),
        fixed_relational_readout=dict(effect='(I+SWAP_BR)/2',free_common_invariance_residual=invariant,
            extra_native_material_reference=True,external_time_dependent_effect=False,
            reference_preparation_is_resource=True,isolation_is_additional_input=True),
        interval_certificate=dict(finite_spectral_derivative_bound=L,joint_derivative_bound=Ljoint,
            grid_points=len(grid),grid_minimum_contrast=float(min(contrast)),mesh_gap_bound=gap,
            uniform_infinite_contrast_lower=lower,probability_error_per_input=probability_error,
            finite_spectral_formula_error=max(formula_errors),phase_check_error=max(phase_errors),
            equal_prior_binary_error_upper=binary_error,binary_information_lower=binary_information,
            whole_interval_not_only_grid_certified=True,endpoints=endpoints),
        infinite_occupation_certificate=cert,
        reference_account=dict(internal_energy_mean=reference_mean,internal_energy_variance=reference_variance,
            reference_characteristic_real=ref_chi.real,reference_characteristic_imag=ref_chi.imag,
            source_mass_positive_lower=positivity_lower,populated_sector_mass_lower=sector_lower,source_rows=source_rows,
            receiver_old_effect_probabilities=old_probs.tolist(),old_endpoint_check=old_error,
            source_includes_complete_reference_energy=True,old_path_predictions_unchanged=False),
        scope=dict(finite_relational_record_adopted=True,reference_mass_and_energy_included=True,
            phase_record_converted_to_protected_population=False,actual_joint_apparatus_derived=False,
            arbitrary_noise_stability_proven=False,repeatable_nondestructive_readout_proven=False,
            macroscopic_irreversible_memory_proven=False,full_microscopic_SM_matching_proven=False,
            moving_support_or_quantum_metric_completed=False,full_goal_completed=False),
        references=['https://arxiv.org/abs/quant-ph/0203016','https://arxiv.org/abs/1703.10434'],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-10),(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    if args.write:assert not OUT.exists()
    result=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as dest:json.dump(result,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(result,read(OUT))
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
