"""604: an explicit kinetic matching test for the original finite CAR matter.

A declared naive gauge-covariant nearest-neighbor spatial derivative is tested.
It was not fixed by 598. Exact frozen-background spectra/source identities do
not claim a nonperturbative chiral continuum or full Gauss thermal calculation.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_quantum_response_matching as response
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_chiral_source_matching_results.json'
SIG=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1,-1]).astype(complex)]
CHARGES=dict(Q=1,u=4,d=-2,L=-3,e=-6,nu=0)
LEFT={'Q','L'}

def kinetic_matrices():
    ks=[]
    for sig in SIG:
        k=np.zeros((32,32),complex)
        for name,sl in matter.SLICES.items():
            size=sl.stop-sl.start
            k[sl,sl]=(-1 if name in LEFT else 1)*np.kron(np.eye(size//2),sig)
        ks.append(k)
    zero=np.zeros((32,32),complex)
    # Hole block is -h(-p)^T, not -h(p)^T.
    return [np.block([[k,zero],[zero,k.T]]) for k in ks]

def mass(s=response.S0):
    return response.matrices(s)[0][0]

def original_kinetic_check():
    G=kinetic_matrices();B=mass();I=np.eye(64)
    err=max(float(np.max(abs(x@y+y@x-2*(i==j)*I))) for i,x in enumerate(G) for j,y in enumerate(G))
    anti=max(float(np.max(abs(x@B+B@x))) for x in G)
    assert max(err,anti)<1e-13
    phase=.271
    R=matter.representation(np.eye(3),np.eye(2),np.exp(1j*phase))
    Z=np.zeros_like(R);RR=np.block([[R,Z],[Z,R.conj()]])
    covariance=max(float(np.max(abs(RR@g-g@RR))) for g in G)
    assert covariance<1e-13
    a=.3;q=np.array([.07,-.09,.11])
    rows=[];spectra=[]
    m2=np.linalg.eigvalsh(B@B)
    for bits in itertools.product((0,1),repeat=3):
        p=np.array(bits)*np.pi/a+q
        d=np.sin(a*p)/a;H=B+sum(x*t for x,t in zip(G,d))
        square_error=float(np.max(abs(H@H-B@B-np.dot(d,d)*I)))
        eigen_error=float(np.max(abs(np.linalg.eigvalsh(H@H)-(m2+np.dot(d,d)))))
        assert max(square_error,eigen_error)<1e-12
        spectra.append(np.linalg.eigvalsh(H))
        rows.append(dict(corner=list(bits),orientation=(-1)**sum(bits),square_error=square_error,eigen_error=eigen_error))
    spread=float(np.max(abs(np.array(spectra)-spectra[0])))
    assert spread<1e-12 and sum(r['orientation'] for r in rows)==0
    return dict(CAR_modes=32,spatial_dimensions_input=3,Nambu_dimension=64,
                Clifford_error=err,mass_anticommutator_error=anti,gauge_kinetic_covariance_error=covariance,
                corners=rows,corner_spectral_spread=spread,
                original_positive_masses=np.linalg.eigvalsh(B)[32:].tolist(),
                hopping_is_newly_declared_branch=True)

def thermo(N,a,beta,s=response.S0,central=False):
    assert N%4==0
    indices=np.arange(-N//4,N//4) if central else np.arange(-N//2,N//2)
    phases=2*np.pi*indices/N
    xx,yy,zz=np.meshgrid(phases,phases,phases,indexing='ij')
    d2=(np.sin(xx)**2+np.sin(yy)**2+np.sin(zz)**2).ravel()/a**2
    masses=np.linalg.eigvalsh(mass(s))[32:]
    assert masses.min()>0
    E=np.sqrt(d2[:,None]+masses[None,:]**2)
    n=1/(1+np.exp(beta*E))
    free=-np.sum(np.logaddexp(0,-beta*E))/beta
    energy=np.sum(E*n);variance=np.sum(E**2*n*(1-n))
    return np.array([free,energy,variance],float)

def common_source_check():
    rows=[]
    for N,a,beta in ((8,.35,.7),(12,.25,1.7),(16,.2,3.)):
        full=thermo(N,a,beta);folded=thermo(N,a,beta,central=True)
        residual=float(np.max(abs(full-8*folded)))
        assert residual<2e-11*max(1,float(np.max(abs(full))))
        source=[]
        eps=2e-5
        for central in (False,True):
            scalar=(thermo(N,a,beta,response.S0+eps,central)[0]-
                    thermo(N,a,beta,response.S0-eps,central)[0])/(2*eps)
            dilation=(thermo(N,a*np.exp(eps),beta,central=central)[0]-
                      thermo(N,a*np.exp(-eps),beta,central=central)[0])/(2*eps)
            source.append(np.array([scalar,dilation]))
        error=float(np.max(abs(source[0]-8*source[1])))
        assert error<2e-7*max(1,float(np.max(abs(source[0]))))
        rows.append(dict(N=N,spacing=a,beta=beta,
            full_thermal_free_energy_energy_variance=full.tolist(),
            central_fold_values=folded.tolist(),eightfold_identity_error=residual,
            full_scalar_and_dilation_sources=source[0].tolist(),
            central_scalar_and_dilation_sources=source[1].tolist(),source_identity_error=error))
    return dict(rows=rows,finite_periodic_grid_folding_is_exact=True,
                thermal_excess_is_state_difference_not_gravitational_vacuum_subtraction=True,
                no_full_interacting_Gauss_trace_claim=True)

def repair_and_running_check():
    left={CHARGES[x] for x in LEFT};right={CHARGES[x] for x in CHARGES if x not in LEFT}
    assert not left&right
    charge=np.zeros(32)
    for name,sl in matter.SLICES.items():charge[sl]=CHARGES[name]
    phi=np.array([0.,response.HIGGS,0.,0.,response.S0])
    h,_=matter.mass_matrices(phi)
    comm=(charge[:,None]-charge[None,:])*h
    assert np.linalg.norm(comm)>1e-3
    rows=[]
    for r in (.5,1.,2.):
        theta=2*np.arctan(1/r)
        value=r*(1-np.cos(theta))-np.sin(theta)
        slope=r*np.sin(theta)-np.cos(theta)
        assert abs(value)<1e-14 and abs(slope-1)<1e-14
        rows.append(dict(Wilson_scalar_r=r,additional_zero_phase=float(theta),
                         lower_band_value=float(value),nonzero_normal_derivative=float(slope)))
    def counts(copies,N=3,ng=3):
        ns=5;nf=4*copies*ng*(N+1);nv=N*N+3
        c=Q(ns+3*nf+12*nv,120)
        b2=Q(-43+2*copies*ng*(N+1),6)
        invariant=c-Q(3,10)*b2
        return dict(copies=copies,physical_Weyl_count=nf,c_W=str(c),b_weak=str(b2),
                    old_spectral_defect=str(invariant))
    target=counts(1);naive=counts(8)
    assert target['old_spectral_defect']==naive['old_spectral_defect']==str(Q(407,120))
    # Exact coefficients of the general flavor product u=copies*ng*(N+1).
    assert Q(12,120)-Q(3,10)*Q(2,6)==0
    return dict(left_integer_hypercharges=sorted(left),right_integer_hypercharges=sorted(right),
        no_constant_left_right_gauge_intertwiner=True,
        fixed_vacuum_mass_gauge_commutator_norm=float(np.linalg.norm(comm)),
        identity_Wilson_term_gapless_surface_rows=rows,
        conditional_target_IR_counts=target,conditional_naive_IR_counts=naive,
        continuum_log_coefficients_require_declared_all_active_Lorentz_window=True,
        no_full_lattice_curved_effective_action_computed=True)

def run():
    a=original_kinetic_check();b=common_source_check();c=repair_and_running_check()
    deps=('research_note_374.md','research_note_531.md','research_note_553.md','research_note_598.md',
          'research_note_599.md','research_note_602.md','research_note_603.md',
          'joint_fermion_gauss_completion.py','joint_quantum_response_matching.py')
    return dict(round=604,tests_run=3,failures=0,errors=0,kinetic=a,common_sources=b,
        repair_and_running=c,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_mass_representation_retained=True,naive_hopping_added_as_test_branch=True,
                   explicit_branch_rejected_not_general_unification=True,
                   eightfold_scalar_thermal_and_scale_sources_proved=True,
                   old_single_combination_cannot_test_species_matching=True,
                   no_general_no_go_or_chiral_continuum_completion=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=604,tests=3,all_passed=True,
        counts=result['repair_and_running']['conditional_naive_IR_counts'])))
