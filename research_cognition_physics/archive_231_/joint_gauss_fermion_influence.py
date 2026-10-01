"""643: original full-CAR conditional influence, Gauss closure and sources.

The full nonlinear heat-kernel representation and its integrability are proved
in the note. Numerics evaluate actual 32-mode ONSITE conditional path factors,
not the bosonic bridge integral or the full interacting Gauss partition function.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_region_energy_gluing as groups

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_fermion_influence_results.json'
PHI=np.array([[.45,.31,-.17,.26,.62],[-.22,.53,.35,.12,-.41],
              [.37,-.19,.28,-.44,.55],[.13,.48,-.39,.21,.33]])
TIMES=np.array([.4,.3,.5,.25])


def exp_h(h,z):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-z*e))@v.conj().T


def word(state,ops):
    """Rightmost operator first; bit order fixes exterior-algebra signs."""
    sign=1
    for j,create in ops:
        occupied=(state>>j)&1
        if occupied==create:return None,0
        sign*=(-1)**((state&((1<<j)-1)).bit_count())
        state^=1<<j
    return state,sign


def fock(h,d):
    n=len(h);out=np.zeros((1<<n,1<<n),complex);pair=np.zeros_like(out)
    for i,j in zip(*np.nonzero(abs(h)>0)):
        for s in range(1<<n):
            t,a=word(s,[(j,False),(i,True)])
            if a:out[t,s]+=h[i,j]*a
    for i in range(n):
        for j in range(i+1,n):
            if not d[i,j]:continue
            for s in range(1<<n):
                t,a=word(s,[(j,True),(i,True)])
                if a:pair[t,s]+=d[i,j]*a
    out+=pair+pair.conj().T
    assert np.max(abs(out-out.conj().T))<1e-13
    return out


def exterior(r):
    n=len(r);out=np.zeros((1<<n,1<<n),complex)
    for k in range(n+1):
        subsets=list(itertools.combinations(range(n),k))
        for a in subsets:
            ia=sum(1<<i for i in a)
            for b in subsets:
                ib=sum(1<<i for i in b)
                out[ia,ib]=np.linalg.det(r[np.ix_(a,b)]) if k else 1.
    return out


def bdg(h,d):return np.block([[h,d],[d.conj().T,-h.T]])


def chain(mats,zs,a=0.):
    """U=E_last...E_first and derivative under B_j -> exp(a) B_j."""
    u=np.eye(len(mats[0]),dtype=complex);du=np.zeros_like(u)
    for b,z in zip(mats,zs):
        e=exp_h(b,z*np.exp(a));de=-z*np.exp(a)*b@e
        du=de@u+e@du;u=e@u
    return u,du


def data(phi=PHI):
    pairs=[matter.mass_matrices(p) for p in phi]
    return dict(pairs=pairs,quark=[h[:24,:24] for h,d in pairs],
                lepton=[fock(h[24:,24:],d[24:,24:]) for h,d in pairs],
                nambu=[bdg(h,d) for h,d in pairs])


def conditional(ds,r,zs=TIMES,a=0.,lift=None):
    rq=r[:24,:24];rl=r[24:,24:]
    if lift is None:lift=exterior(rl)
    uq,dq=chain(ds['quark'],zs,a);ul,dl=chain(ds['lepton'],zs,a)
    mq=rq@uq;dmq=rq@dq
    q=np.linalg.det(np.eye(24)+mq)
    qprime=q*np.trace(np.linalg.solve(np.eye(24)+mq,dmq))
    lepton=np.trace(lift@ul);lprime=np.trace(lift@dl)
    value=q*lepton/2**32;derivative=(qprime*lepton+q*lprime)/2**32
    u,du=chain(ds['nambu'],zs,a)
    rn=np.block([[r,np.zeros_like(r)],[np.zeros_like(r),r.conj()]])
    m=rn@u;dm=rn@du
    determinant=np.linalg.det((np.eye(64)+m)/2)
    square_error=abs(value*value-determinant)
    gaussian_derivative=.5*value*np.trace(np.linalg.solve(np.eye(64)+m,dm))
    return value,derivative,float(square_error),complex(gaussian_derivative)


def fixture():
    c=groups.gauge.group_exp(np.array([.17,-.11,.08,.05,.09,-.06,.04,.12]),3)
    w=groups.gauge.group_exp(np.array([.21,-.16,.28]),2)
    return c,w,np.exp(.13j)


def pair(x):return [float(x.real),float(x.imag)]


def conditional_trace_check():
    ds=data();r=matter.representation(*fixture());lift=exterior(r[24:,24:])
    rows=[]
    for name,zs in [('Euclidean',TIMES),('mixed_contour',np.array([.7,.4j,-.6j,.25j]))]:
        value,derivative,err,direct=conditional(ds,r,zs,lift=lift)
        assert err<2e-12 and abs(derivative-direct)<2e-12
        rows.append(dict(contour=name,z=[pair(z) for z in zs],
            full32_normalized_trace=pair(value),trace_square_determinant_error=err,
            conditional_lapse_insertion=pair(derivative),insertion_error=float(abs(derivative-direct))))
    # A distinct explicit Fock check includes paired and unpaired original modes.
    ids=[24,25,30,31]
    h,d=ds['pairs'][0];hf,df=h[np.ix_(ids,ids)],d[np.ix_(ids,ids)]
    hfock=fock(hf,df);tr=np.trace(exp_h(hfock,.63))
    exact_error=float(abs(tr*tr-np.linalg.det(np.eye(8)+exp_h(bdg(hf,df),.63))))
    comm=max(float(np.linalg.norm(ds['nambu'][i]@ds['nambu'][j]-
          ds['nambu'][j]@ds['nambu'][i],2)) for i in range(4) for j in range(i))
    assert exact_error<2e-11 and comm>1e-3
    return dict(rows=rows,independent_4_mode_Fock_error=exact_error,
        original_BdG_path_noncommutativity=comm,
        normalized_by_Fock_dimension=2**32,explicit_paired_Fock_dimension=256,
        all_original32_modes_retained=True,
        no_full_bosonic_path_integration_or_Gauss_average_performed=True)


def gauge_and_order_check():
    ds=data();g=fixture();r=matter.representation(*g)
    value,*_=conditional(ds,r)
    rng=np.random.default_rng(643)
    k=groups.sample(rng)
    transformed=[]
    for p in PHI:
        X=k[2]**3*k[1]@(p[:2]+1j*p[2:4]);transformed.append(np.r_[X.real,X.imag,p[4]])
    gp=groups.product(groups.product(k,g),groups.inverse(k))
    same,*_=conditional(data(np.array(transformed)),matter.representation(*gp))
    without,*_=conditional(ds,np.eye(32))
    center=(np.exp(2j*np.pi/3)*np.eye(3),-np.eye(2),np.exp(1j*np.pi/3))
    rc=matter.representation(*groups.product(center,g))
    center_error=float(np.max(abs(rc-r)))
    # Replace the ordered product by the exponential of its integrated original B.
    meanh=sum(t*p[0] for t,p in zip(TIMES,ds['pairs']))
    meand=sum(t*p[1] for t,p in zip(TIMES,ds['pairs']))
    naive=dict(quark=[meanh[:24,:24]],lepton=[fock(meanh[24:,24:],meand[24:,24:])],
               nambu=[bdg(meanh,meand)])
    averaged,*_=conditional(naive,r,zs=[1.])
    omitted_pair=dict(quark=ds['quark'],
        lepton=[fock(h[24:,24:],np.zeros((8,8))) for h,d in ds['pairs']],
        nambu=[bdg(h,np.zeros_like(d)) for h,d in ds['pairs']])
    no_majorana,*_=conditional(omitted_pair,r)
    assert abs(same-value)<2e-12 and center_error<3e-14
    assert abs(without-value)>1e-4 and abs(averaged-value)>1e-6 and abs(no_majorana-value)>1e-6
    return dict(original_simultaneous_gauge_covariance_error=float(abs(same-value)),
        quotient_lift_error=center_error,full32_gauge_determinant_error=float(abs(np.linalg.det(r)-1)),
        original_weight=pair(value),omitting_Gauss_closure_weight=pair(without),
        replacing_ordered_history_by_average_weight=pair(averaged),
        deleting_original_Majorana_weight=pair(no_majorana),
        closure_omission_difference=float(abs(without-value)),
        ordering_omission_difference=float(abs(averaged-value)),
        Majorana_omission_difference=float(abs(no_majorana-value)),
        these_are_conditional_factors_not_physical_partition_functions=True)


def source_check():
    ds=data();r=matter.representation(*fixture());lift=exterior(r[24:,24:])
    value,analytic,err,gaussian=conditional(ds,r,lift=lift)
    rows=[]
    for step in (2e-4,1e-4):
        plus=conditional(ds,r,a=step,lift=lift)[0]
        minus=conditional(ds,r,a=-step,lift=lift)[0]
        fd=(plus-minus)/(2*step)
        rows.append(dict(step=step,finite_difference=pair(fd),error=float(abs(fd-analytic))))
    assert rows[-1]['error']<1e-8 and rows[-1]['error']<.35*rows[0]['error']
    return dict(conditional_weight=pair(value),original_mass_lapse_insertion=pair(analytic),
        logarithmic_conditional_derivative=pair(analytic/value),rows=rows,
        original_full_H_source_also_requires_bosonic_kinetic_and_potential_terms=True,
        differentiation_of_the_full_source_is_operator_Duhamel_not_a_Brownian_score=True)


def run():
    deps=('research_note_598.md','research_note_603.md','research_note_623.md','research_note_642.md',
          'joint_fermion_gauss_completion.py','joint_region_energy_gluing.py')
    return dict(round=643,tests_run=3,failures=0,errors=0,
        conditional_full_matter=conditional_trace_check(),gauge_order_and_pairing=gauge_and_order_check(),
        common_source=source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(full_original_nonlinear_Gauss_heat_representation_proved_in_note=True,
            finite_graph_and_given_positive_geometry=True,
            fermion_elimination_keeps_order_and_Gauss_closure=True,
            numerical_results_are_original32_mode_conditional_factors=True,
            no_continuum_chiral_GR_or_autonomous_design_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
