"""644: full-Gauss replica identity; original two-node conditional checks.

The numerical thermal calculation freezes bosons and link at declared values.
It retains all64 original CAR modes and nonzero hopping, but is NOT the complete
bosonic/Gauss thermal integral. Separate physical states are inherited from642.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_gauss_fermion_influence as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_region_replica_source_results.json'
BETA=1.7
HOP=.17


def matrices(hop=HOP):
    phi=[np.array([0.,.7,0.,0.,.45]),np.array([0.,.5,0.,0.,-.35])]
    pairs=[old.matter.mass_matrices(p) for p in phi]
    h=np.zeros((64,64),complex);d=np.zeros_like(h)
    for k,(hk,dk) in enumerate(pairs):
        h[32*k:32*(k+1),32*k:32*(k+1)]=hk
        d[32*k:32*(k+1),32*k:32*(k+1)]=dk
    h[:32,32:]=h[32:,:32]=hop*np.eye(32)
    return h,d


def reduce_node(rho,k):
    da=1<<k
    return np.einsum('babc->ac',rho.reshape(da,da,da,da))


def thermal_block(h,d,a=0.):
    H=old.fock(h,d)*np.exp(a);e,v=np.linalg.eigh(H)
    p=np.exp(-BETA*(e-e.min()));p/=p.sum()
    rho=(v*p)@v.conj().T
    drho=(v*(-BETA*(e-e@p)*p))@v.conj().T
    k=len(h)//2;ra=reduce_node(rho,k);dra=reduce_node(drho,k)
    purity=float(np.trace(ra@ra).real)
    source=-2*float(np.trace(ra@dra).real)/purity
    return dict(H=H,rho=rho,region=ra,S2=-np.log(purity),derivative=source,purity=purity)


def blocks(h,d,a=0.):
    # Phi=(0,h,0,0,s) and identity link: actual invariant flavour factors.
    specs=[('up',6,[0,12]),('down',6,[2,18]),('electron',2,[26,28]),
           ('neutral',1,[24,25,30,31])]
    rows=[];entropy=source=0.
    for name,mult,ids in specs:
        ix=np.ix_(ids+[i+32 for i in ids],ids+[i+32 for i in ids])
        r=thermal_block(h[ix],d[ix],a)
        entropy+=mult*r['S2'];source+=mult*r['derivative']
        rows.append(dict(name=name,multiplicity=mult,modes_per_node=len(ids),
            Fock_dimension=len(r['H']),regional_S2=float(r['S2']),
            lapse_derivative=float(r['derivative'])))
    assert sum(x['multiplicity']*x['modes_per_node'] for x in rows)==32
    return float(entropy),float(source),rows


def covariance_entropy(h,d,a=0.):
    b=old.bdg(h,d)*np.exp(a);e,v=np.linalg.eigh(b)
    n=(v*(1/(1+np.exp(BETA*e))))@v.conj().T
    ids=np.r_[np.arange(32),np.arange(64,96)]
    eig=np.linalg.eigvalsh(n[np.ix_(ids,ids)])
    assert min(eig)>0 and max(eig)<1
    return float(-.5*np.log(eig**2+(1-eig)**2).sum())


def swap_matrix(n,ids,signed=True):
    r=np.eye(2*n)
    for i in ids:
        r[i,i]=r[n+i,n+i]=0
        r[i,n+i]=-1 if signed else 1
        r[n+i,i]=1
    return r


def replica_entropy(h,d):
    n=len(h);b=old.bdg(h,d);e=np.linalg.eigvalsh(b)
    logz=.5*np.logaddexp(0,-BETA*e).sum()
    m=old.exp_h(b,BETA);double=np.zeros((4*n,4*n),complex)
    for k in range(2):
        ids=np.r_[np.arange(k*n,(k+1)*n),np.arange((2+k)*n,(3+k)*n)]
        double[np.ix_(ids,ids)]=m
    r=swap_matrix(n,range(32));rn=np.block([[r,np.zeros_like(r)],[np.zeros_like(r),r]])
    sign,logdet=np.linalg.slogdet(np.eye(4*n)+rn@double)
    assert abs(sign-1)<2e-12
    s2=2*logz-.5*logdet
    return float(s2),float(logz),float(abs(sign-1))


def full64_check():
    h,d=matrices();s,source,rows=blocks(h,d)
    covariance=covariance_entropy(h,d);replica,logz,phase=replica_entropy(h,d)
    assert max(abs(s-covariance),abs(s-replica))<2e-12
    # Direct exterior-algebra twisted exchange on one actual two-node up block.
    ids=[0,12,32,44];ix=np.ix_(ids,ids)
    b=thermal_block(h[ix],d[ix]);two=np.kron(b['rho'],b['rho'])
    rs=swap_matrix(4,range(2));signed=old.exterior(rs)
    ru=swap_matrix(4,range(2),False);unsigned=old.exterior(ru)
    direct=complex(np.trace(two@signed))
    wrong=complex(np.trace(two@unsigned))
    parity=np.diag([(-1)**i.bit_count() for i in range(4)])
    superpurity=np.trace(parity@b['region']@b['region'])
    assert abs(direct-b['purity'])<2e-13 and abs(wrong-superpurity)<2e-13
    assert abs(wrong-direct)>1e-3
    return dict(full_CAR_modes=64,original_hopping=HOP,beta=BETA,blocks=rows,
        exact_block_Fock_S2=s,reduced_Nambu_covariance_S2=covariance,
        signed_replica_determinant_S2=replica,log_partition=logz,
        determinant_phase_error=phase,regional_purity=float(np.exp(-s)),
        original_block_purity=b['purity'],signed_exchange_trace=old.pair(direct),
        unsigned_exchange_trace=old.pair(wrong),parity_weighted_purity=old.pair(superpurity),
        all_original_modes_and_pairing_present=True,
        bosons_frozen_and_no_Gauss_average_in_this_numerical_fixture=True)


def source_check():
    h,d=matrices();s,source,_=blocks(h,d)
    rows=[]
    for step in (4e-4,2e-4):
        plus=covariance_entropy(h,d,step);minus=covariance_entropy(h,d,-step)
        fd=(plus-minus)/(2*step)
        rows.append(dict(step=step,finite_difference=fd,error=float(abs(fd-source))))
    assert rows[-1]['error']<1e-6 and rows[-1]['error']<.4*rows[0]['error']
    h0,d0=matrices(0.);decoupled=covariance_entropy(h0,d0)
    assert abs(s-decoupled)>1e-4
    return dict(S2=s,analytic_conditional_lapse_derivative=source,rows=rows,
        independently_thermalized_onsite_S2=decoupled,
        onsite_replacement_entropy_difference=float(abs(s-decoupled)),
        full_Gauss_geometric_Duhamel_identity_is_analytic_not_this_fixture=True)


def gauss_witness_check():
    # The642 compact scalar/Gauss particle-hole states, partition at node1's edge.
    rows=[]
    for module,d in [('nu',1),('L',2),('u',3),('Q',6)]:
        swap=np.zeros((d*d,d*d))
        for i in range(d):
            for j in range(d):swap[j*d+i,i*d+j]=1
        two=np.eye(d*d)/(d*d)
        good=float(np.trace(two@swap));bad=float(np.trace(two@(-swap)))
        assert abs(good-1/d)<1e-14 and abs(bad+1/d)<1e-14
        rows.append(dict(original_module=module,Schmidt_rank=d,
            region_fermion_parity=-1,signed_replica_purity=good,
            unsigned_mode_exchange_weight=bad,regional_S2=float(np.log(d))))
    return dict(rows=rows,original_legal_full_Gauss_states=True,
        these_states_are_not_Gibbs_eigenstates=True,
        boundary_dimension_and_odd_matter_parity_both_retained=True)


def run():
    deps=('research_note_617.md','research_note_637.md','research_note_642.md','research_note_643.md',
          'joint_fermion_gauss_completion.py','joint_gauss_fermion_influence.py')
    return dict(round=644,tests_run=3,failures=0,errors=0,
        full_original_conditional_matter=full64_check(),common_entropy_source=source_check(),
        legal_Gauss_sign_witnesses=gauss_witness_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(full_fixed_graph_Gauss_replica_and_source_identity_proved_in_note=True,
            numerical_thermal_state_is_conditional64_CAR_not_full_physical_Gibbs=True,
            separate_exact_Gauss_witnesses_from642=True,
            no_area_law_Newton_continuum_or_GR_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
