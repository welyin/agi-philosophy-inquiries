"""1083: one flagged-qutrit source for comparison and source routing.

Fixed 18-dimensional comparator, retained fine records, and a passive reference.
Exact Fraction flags/cycle margins are separate from finite dense-matrix checks.
Default prints JSON; --check recomputes without writing files.
"""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import json
import platform
import numpy as np

I2=np.eye(2,dtype=complex)
I3=np.eye(3,dtype=complex)
PAULI=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
READY=np.diag([1,1,0]).astype(complex)
VAC=np.diag([0,0,1]).astype(complex)
ZERO=(F(0),F(0),F(0))
TOL=3e-12


def fv(*xs):return tuple(F(x) for x in xs)
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def sub(x,y):return tuple(a-b for a,b in zip(x,y))
def mul(a,x):return tuple(F(a)*b for b in x)
def dot(x,y):return sum((a*b for a,b in zip(x,y)),F(0))
def outvec(x):return [str(a) for a in x]

def kron(*matrices):
    result=np.array([[1]],complex)
    for matrix in matrices:result=np.kron(result,matrix)
    return result

def rho(p):
    return (I2+sum((float(x)*s for x,s in zip(p,PAULI)),np.zeros((2,2),complex)))/2

def omega(q,p):
    result=np.zeros((3,3),complex)
    result[:2,:2]=float(q)*rho(p)
    result[2,2]=float(1-q)
    return result

def conditional_source(flag,p):
    return omega(F(1),p) if flag else VAC.copy()

def trace_out(state,dims,keep):
    keep=list(keep)
    traced=[i for i in range(len(dims)) if i not in keep]
    n=len(dims)
    axes=keep+traced+[i+n for i in keep]+[i+n for i in traced]
    dk=int(np.prod([dims[i] for i in keep]))
    dt=int(np.prod([dims[i] for i in traced]))
    tensor=state.reshape(tuple(dims)*2).transpose(axes).reshape(dk,dt,dk,dt)
    return np.einsum('abcb->ac',tensor)

def reordered_product_qr_ab(qr,ab):
    # Q,R,A,B -> Q,A,B,R, for the fixed comparator tensored with I_R.
    return np.kron(qr,ab).reshape(2,2,3,3,2,2,3,3).transpose(0,2,3,1,4,6,7,5).reshape(36,36)

def swap_qa(which):
    s=np.zeros((18,18),complex)
    for q,a,b in product(range(2),range(3),range(3)):
        before=q*9+a*3+b
        target=[q,a,b]
        index=1 if which=='A' else 2
        if target[index]<2:
            target[0],target[index]=target[index],target[0]
        after=target[0]*9+target[1]*3+target[2]
        s[after,before]=1
    return s

def source_swap():
    s=np.zeros((9,9),complex)
    for a,b in product(range(3),repeat=2):s[b*3+a,a*3+b]=1
    return s

def blockdiag(a,b):
    out=np.zeros((len(a)+len(b),len(a)+len(b)),complex)
    out[:len(a),:len(a)]=a
    out[len(a):,len(a):]=b
    return out

def comparator():
    sa,sb=swap_qa('A'),swap_qa('B')
    identity=np.eye(18,dtype=complex)
    both=kron(I2,READY,READY)
    branches=[]
    # Store CP branch as (coarse outcome, weight, unscaled operator, label).
    # Its Kraus operator is sqrt(weight)*operator. Binary rational entries
    # keep this effect/completeness check exact in ordinary floating point.
    for label,swap in [('A',sa),('B',sb)]:
        for sign in (1,-1):
            projector=(identity+sign*swap)/2
            outcome=0 if (label=='A' and sign==1) or (label=='B' and sign==-1) else 1
            branches.append((outcome,F(1,2),projector@both,f'both:{label}:swap{sign:+d}'))
    for a,b in ((0,0),(0,1),(1,0)):
        projector=kron(I2,READY if a else VAC,READY if b else VAC)
        for outcome in (0,1):
            branches.append((outcome,F(1,2),projector,f'flags{a}{b}:coin{outcome}'))
    effects=[sum((float(w)*(k.conj().T@k) for o,w,k,_ in branches if o==outcome),np.zeros((18,18),complex)) for outcome in (0,1)]
    formula=identity/2+both@(sa-sb)@both/4
    return branches,effects,formula,sa,sb

def effective_probe(effect,program_ab):
    return np.einsum('iajb,ba->ij',effect.reshape(2,9,2,9),program_ab)

def formula_effect(qx,px,qy,py):
    return I2/2+sum((float(qx*qy*(x-y)/8)*s for x,y,s in zip(px,py,PAULI)),np.zeros((2,2),complex))

class Audit:
    def __init__(self):
        self.calls=0
        self.groups=[]
        self.max_residual=0.0
    def require(self,condition,label):
        self.calls+=1
        if not condition:raise AssertionError(label)
    def near(self,a,b,label):
        error=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
        self.max_residual=max(self.max_residual,error)
        self.require(error<TOL,label)
        return error
    def group(self,name,data):
        self.groups.append({'name':name,'passed':True,'certificate':data})


def flagged_distribution(correlated):
    if not correlated:
        return {bits:F(1,16) for bits in product((0,1),repeat=4)}
    return {(a,z,z,b):F(1,8) for a,z,b in product((0,1),repeat=3)}

def flag_statistics(distribution):
    q=[sum((p*bits[i] for bits,p in distribution.items()),F(0)) for i in range(4)]
    h=[[sum((p*bits[i]*bits[j] for bits,p in distribution.items()),F(0)) for j in range(4)] for i in range(4)]
    return q,h

def full_source(distribution,positions):
    state=np.zeros((81,81),complex)
    for flags,p in distribution.items():
        state+=float(p)*kron(*(conditional_source(bit,v) for bit,v in zip(flags,positions)))
    return state

def pair_source(distribution,positions,i,j):
    if i==j:raise ValueError('One unknown physical source is not duplicated; equal labels use the neutral branch.')
    state=np.zeros((9,9),complex)
    for flags,p in distribution.items():
        state+=float(p)*kron(conditional_source(flags[i],positions[i]),conditional_source(flags[j],positions[j]))
    return state

def mixed_vector(p,q,h,positions):
    result=ZERO
    for i,pi in enumerate(p):
        for j,qj in enumerate(q):
            result=add(result,mul(pi*qj*h[i][j],sub(positions[i],positions[j])))
    return result

def run():
    audit=Audit()
    branches,effects,global_formula,sa,sb=comparator()
    audit.require(len(branches)==10,'the physical comparator has ten retained fine branches')
    audit.near(sum(effects),np.eye(18),'complete comparator instrument is trace preserving')
    audit.near(effects[0],global_formula,'coarse report equals the fixed 18D operator formula')
    audit.near(sa@sa,np.eye(18),'ready-block SWAP QA is an involution')
    audit.near(sb@sb,np.eye(18),'ready-block SWAP QB is an involution')
    spectra=np.linalg.eigvalsh(effects[0])
    audit.require(float(spectra.min())>-TOL and float(spectra.max())<1+TOL,'global report effect is legal')
    cases=[(F(2,3),fv(F(1,3),F(1,4),F(1,5)),F(1,4),fv(F(-1,2),F(1,3),0)),
           (F(1),fv(1,0,0),F(1),fv(0,1,0)),
           (F(1,5),fv(0,0,0),F(3,4),fv(0,0,-1))]
    case_reports=[]
    for qx,px,qy,py in cases:
        ab=kron(omega(qx,px),omega(qy,py))
        actual=effective_probe(effects[0],ab)
        expected=formula_effect(qx,px,qy,py)
        err=audit.near(actual,expected,'source-independent fixed comparator gives qx*qy weighted difference')
        audit.require(np.linalg.eigvalsh(actual).min()>-TOL and np.linalg.eigvalsh(I2-actual).min()>-TOL,'induced probe effect is legal')
        case_reports.append({'qx':str(qx),'px':outvec(px),'qy':str(qy),'py':outvec(py),
                             'formula_residual':err,'actual_eigenvalues':[float(x) for x in np.linalg.eigvalsh(actual)]})
    audit.group('fixed_18_dimensional_comparator',{
        'source_Hilbert_space':'C^2_ready direct-sum C_vac; basis ready0, ready1, vac2',
        'dimension':18,'fine_branch_count':len(branches),
        'fine_branch_labels':[label for _,_,_,label in branches],
        'effect_formula':'T_plus=I18/2+P11*(S_QA-S_QB)*P11/4',
        'induced_effect':'E_xy=I2/2+q_x*q_y*(p_x-p_y).sigma/8',
        'global_effect_spectrum':[float(x) for x in spectra],
        'finite_program_cases':case_reports,
        'scope':'The matrix checks calibrate the operator identity; universality is established analytically.'})

    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    asym=np.array([1,0,0,2],complex)/np.sqrt(5)
    basis=np.array([0,1,0,0],complex)
    references=[np.outer(bell,bell.conj()),np.outer(asym,asym.conj()),
                np.outer(bell,bell.conj())/3+2*np.outer(basis,basis.conj())/3]
    ref_results=[]
    for ci,(qx,px,qy,py) in enumerate(cases):
        ab=kron(omega(qx,px),omega(qy,py))
        e=formula_effect(qx,px,qy,py)
        for ri,qr in enumerate(references):
            joint=reordered_product_qr_ab(qr,ab)
            outputs=[np.zeros((36,36),complex),np.zeros((36,36),complex)]
            total_trace=0.0
            for outcome,w,k,label in branches:
                kr=np.kron(k,I2)
                fine=float(w)*kr@joint@kr.conj().T
                outputs[outcome]+=fine
                total_trace+=float(np.trace(fine).real)
            actual=blockdiag(*(trace_out(o,[2,3,3,2],[3]) for o in outputs))
            expected=blockdiag(*(trace_out(np.kron(f,I2)@qr,[2,2],[1]) for f in (e,I2-e)))
            err=audit.near(actual,expected,'complete classical-report and passive-R map matches the effect')
            audit.near(total_trace,1,'all retained fine branches have total probability one')
            audit.near(trace_out(actual,[2,2],[1]),trace_out(qr,[2,2],[1]),'passive reference marginal is preserved')
            audit.require(np.linalg.eigvalsh((actual+actual.conj().T)/2).min()>-TOL,'report-reference state is positive')
            ref_results.append({'program_case':ci,'reference_case':ri,'CQ_report_residual':err})
    audit.group('passive_reference_and_full_report_map',{
        'reference_dimension':2,'program_probe_independence':'rho_QR tensor Omega_x tensor Omega_y',
        'finite_reference_cases':len(references),'program_cases':len(cases),'checks':ref_results,
        'unknown_probe_uses_per_trial':1,'unknown_probe_copied':False,
        'scope':'Full C_report tensor R is checked. Program poststates and fine source/branch records are retained, not asserted equivalent from q,p alone.'})

    s=source_swap()
    meeting_results=[]
    for qx,px,qy,py in cases:
        ox,oy=omega(qx,px),omega(qy,py)
        ab=kron(ox,oy)
        retained=blockdiag(ab/2,s@ab@s.conj().T/2) # coin, active port A, unused port B
        output=trace_out(retained,[2,3,3],[1])
        expected=(ox+oy)/2
        err=audit.near(output,expected,'actual fair controlled SWAP gives the averaged whole source')
        qm=(qx+qy)/2
        bm=mul(F(1,2),add(mul(qx,px),mul(qy,py)))
        pm=mul(1/qm,bm)
        audit.near(output,omega(qm,pm),'same q and B are affine source coordinates')
        audit.near(np.trace(output[:2,:2]),float(qm),'ready probability is the same comparison weight')
        b_measured=np.array([np.trace(output[:2,:2]@sig).real for sig in PAULI])
        audit.near(b_measured,[float(x) for x in bm],'ready Bloch moment is B=q*p')
        audit.near(trace_out(retained,[2,3,3],[0]),I2/2,'routing coin remains stored')
        recovery=blockdiag(np.eye(9),s)
        audit.near(recovery@retained@recovery.conj().T,blockdiag(ab/2,ab/2),'stored coin permits recovery of the original two complete source registers')
        meeting_results.append({'qx':str(qx),'qy':str(qy),'q_m':str(qm),'B_m':outvec(bm),'p_m':outvec(pm),
                                'whole_source_average_residual':err})
    # Retaining both routed outputs does not create two independent copies.
    z0=np.diag([1,0,0]).astype(complex)
    z1=np.diag([0,1,0]).astype(complex)
    joint_route=(kron(z0,z1)+kron(z1,z0))/2
    one=(z0+z1)/2
    independent=kron(one,one)
    clone_gap=float(np.sum(np.abs(np.linalg.eigvalsh(joint_route-independent)))/2)
    audit.near(clone_gap,F(1,2),'two equal routed marginals are correlated, not independent copies')
    # Source-reference CP control: A ready-Bell with R, B fixed and independent.
    qa=F(2,3)
    ar=np.zeros((6,6),complex)
    ar[:4,:4]=float(qa)*np.outer(bell,bell.conj())
    ar[4:,4:]=float(1-qa)*I2/2
    ob=omega(F(3,4),fv(F(1,3),F(-1,4),0))
    abr=np.kron(ar,ob).reshape(3,2,3,3,2,3).transpose(0,2,1,3,5,4).reshape(18,18)
    sr=np.kron(s,I2)
    retained_r=blockdiag(abr/2,sr@abr@sr.conj().T/2)
    audit.near(trace_out(retained_r,[2,3,3,2],[3]),trace_out(abr,[3,3,2],[2]),'source passive reference remains internally accounted for')
    audit.near(trace_out(retained_r,[2,3,3,2],[1]),
               (trace_out(abr,[3,3,2],[0])+trace_out(abr,[3,3,2],[1]))/2,
               'whole-source averaging is valid with a passive source reference')
    audit.group('whole_source_routing_and_same_weight',{
        'source_channel':'fair identity or source SWAP; keep coin label and the unselected source',
        'moment_rule':'q_m=(q_x+q_y)/2; B_m=(q_x*p_x+q_y*p_y)/2; p_m=B_m/q_m',
        'conditional_ready_participation':'q_x/(q_x+q_y) versus q_y/(q_x+q_y)',
        'finite_cases':meeting_results,'independent_clone_trace_distance':clone_gap,
        'controls_depend_on_unknown_q_or_p':False,
        'scope':'The active source marginal closes for this routing and fresh independent comparisons. Arbitrary reused partners or retained logs are not compressed into this marginal.'})

    positions=[ZERO,fv(1,0,0),fv(0,1,0),fv(0,0,1)]
    distributions=[flagged_distribution(False),flagged_distribution(True)]
    statistics=[flag_statistics(dist) for dist in distributions]
    full=[full_source(dist,positions) for dist in distributions]
    expected_local=[omega(F(1,2),p) for p in positions]
    for si,state in enumerate(full):
        audit.near(np.trace(state),1,'global flagged program preparation is normalized')
        audit.require(np.linalg.eigvalsh(state).min()>-TOL,'global program preparation is positive')
        for i in range(4):
            audit.near(trace_out(state,[3,3,3,3],[i]),expected_local[i],
                       'independent and correlated batches have identical single-source Omega')
        q,h=statistics[si]
        audit.require(q==[F(1,2)]*4,'all four readiness marginals are exactly one half')
        for i,j in combinations(range(4),2):
            expected_h=F(1,2) if si==1 and (i,j)==(1,2) else F(1,4)
            audit.require(h[i][j]==expected_h,'joint readiness probability has the specified exact value')
            pair=pair_source(distributions[si],positions,i,j)
            actual=effective_probe(effects[0],pair)
            expected=I2/2+sum((float(h[i][j]*(x-y)/8)*sig for x,y,sig in zip(positions[i],positions[j],PAULI)),np.zeros((2,2),complex))
            audit.near(actual,expected,'unchanged comparator responds to joint readiness, not just marginal q products')
    audit.require(np.max(np.abs(full[0]-full[1]))>1e-3,'joint sources differ despite identical local sources')
    audit.group('identical_local_sources_different_ready_correlations',{
        'positions':list(map(outvec,positions)),
        'independent_flags':[{'bits':list(k),'probability':str(v)} for k,v in distributions[0].items()],
        'correlated_flags':[{'bits':list(k),'probability':str(v)} for k,v in distributions[1].items()],
        'correlated_rule':'R0,Z,R3 independent fair bits, and R1=R2=Z',
        'readiness_marginals':[outvec(q) for q,h in statistics],
        'joint_readiness_matrices':[[outvec(row) for row in h] for q,h in statistics],
        'program_rule':'At flag 1 use rho(p_i), at flag 0 use vac; conditional ready program does not depend on other flags.',
        'actual_effect':'I/2+h_ij*(p_i-p_j).sigma/8 for distinct i,j',
        'same_label_branch':'If lottery selections i=j, the controller reports a neutral coin directly; it never copies one unknown source into two input ports.',
        'correlation_is_extra_spatial_dimension':False})

    cycle=[(F(0),F(1,6),F(2,3),F(1,6)),
           (F(0),F(1,3),F(1,6),F(1,2)),
           (F(1,6),F(1,6),F(0),F(2,3))]
    r=fv(F(1,4),F(1,2),F(3,4))
    all_cycles=[]
    for q,h in statistics:
        margins=[dot(r,mixed_vector(cycle[i],cycle[(i+1)%3],h,positions))/8 for i in range(3)]
        all_cycles.append(margins)
    audit.require(all_cycles[0]==[F(-1,768),F(0),F(1,768)],'independent flags have no strict cycle')
    audit.require(all_cycles[1]==[F(1,4608),F(1,4608),F(1,2304)],'shared-ready flag creates a strict cycle with the same comparator')
    integer_matrix=[[0,-1,-2,-3],[1,0,-2,-2],[2,2,0,-1],[3,2,1,0]]
    integer_lotteries=[(0,1,4,1),(0,2,1,3),(1,1,0,4)]
    bilinear=[sum(integer_lotteries[i][a]*integer_matrix[a][b]*integer_lotteries[(i+1)%3][b] for a in range(4) for b in range(4)) for i in range(3)]
    audit.require(bilinear==[1,1,2],'small integer matrix independently verifies the three cyclic margins')
    for i,j in product(range(4),repeat=2):
        exact=dot(r,mul(statistics[1][1][i][j],sub(positions[i],positions[j])))/8
        audit.require(exact==F(integer_matrix[i][j],128),'integer matrix entries equal the physical report margins')
    tieP=(F(0),F(1,3),F(0),F(2,3))
    tieQ=(F(0),F(0),F(1),F(0))
    tieR=(F(1),F(0),F(0),F(0))
    tr=fv(F(1,2),0,F(-1,2))
    tie_all=[]
    for q,h in statistics:
        tie_all.append([dot(tr,mixed_vector(p,q2,h,positions))/8 for p,q2 in ((tieP,tieQ),(tieQ,tieR),(tieP,tieR))])
    audit.require(tie_all[0]==[F(-1,192),F(0),F(-1,192)],'independent readiness preserves binary tie substitution')
    audit.require(tie_all[1]==[F(0),F(0),F(-1,192)],'correlated readiness breaks binary tie substitution')
    audit.group('finite_same_instrument_cycle_and_tie_witnesses',{
        'cycle_probe':outvec(r),'cycle_lotteries':list(map(outvec,cycle)),
        'independent_cycle_probability_margins':outvec(all_cycles[0]),
        'correlated_cycle_probability_margins':outvec(all_cycles[1]),
        'integer_report_margin_matrix_denominator':128,'integer_report_margin_matrix':integer_matrix,
        'integer_lotteries_denominator':6,'integer_lotteries':[list(row) for row in integer_lotteries],
        'integer_cyclic_products':bilinear,
        'binary_probe':outvec(tr),'binary_lotteries':list(map(outvec,[tieP,tieQ,tieR])),
        'binary_pair_order':['P,Q','Q,R','P,R'],
        'independent_binary_probability_margins':outvec(tie_all[0]),
        'correlated_binary_probability_margins':outvec(tie_all[1]),
        'source_marginals_and_instrument_held_fixed':True})

    epsilon=F(1,9216)
    lower=[margin-epsilon for margin in all_cycles[1]]
    audit.require(min(lower)==epsilon and all(x>0 for x in lower),'half-margin error budget preserves the strict three-cycle')
    audit.group('finite_probability_error_scope',{
        'per_coarse_probability_absolute_error_upper':str(epsilon),
        'strict_cycle_remaining_positive_margins':outvec(lower),
        'independent_fixed_sampling_finite_precision_claim':'Only probability error bounds are assumed; no sample cost is derived.',
        'exact_ties_claimed_robust_under_nonzero_noise':False,
        'finite_matrix_tolerance':TOL})

    return {'round':1083,'status':'matrix_and_exact_fraction_checks_passed',
            'scientific_calibration_groups':1,'validation_group_count':len(audit.groups),
            'require_call_count':audit.calls,'maximum_matrix_residual':audit.max_residual,
            'groups':audit.groups,
            'scope':{'actual_contact_or_position_identity_derived':False,
                     'actual_spatial_dimension_derived':False,
                     'source_independence_derived_from_cognitive_principles':False,
                     'same_q_is_comparison_weight_and_ready_participation_rate':True,
                     'only_single_source_marginals_determine_correlated_comparisons':False,
                     'all_natural_future_tasks_closed':False,
                     'unknown_probe_or_program_copied':False,
                     'ready_vac_coherence_or_arbitrary_source_probe_correlations_in_scalar_formula':False,
                     'retained_unselected_sources_or_fine_records_deleted':False},
            'runtime':{'python':platform.python_version(),'numpy':np.__version__},
            'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='read-only exact comparison with saved JSON on the same runtime')
    args=parser.parse_args()
    result=run()
    if args.check:
        path=Path(__file__).with_name('results.json')
        if result!=json.loads(path.read_text(encoding='utf-8')):
            raise AssertionError('saved results differ from recomputation')
        print(json.dumps({'round':1083,'status':'exact_saved_results_match',
                          'validation_group_count':result['validation_group_count'],
                          'require_call_count':result['require_call_count'],
                          'maximum_matrix_residual':result['maximum_matrix_residual'],
                          'checker_sha256':result['checker_sha256'],
                          'results_sha256':hashlib.sha256(path.read_bytes()).hexdigest()},ensure_ascii=False,indent=2))
    else:
        print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
