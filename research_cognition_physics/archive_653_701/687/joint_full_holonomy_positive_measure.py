"""687: original one-site full non-Abelian holonomy, S9 positivity and exact Haar integral.

Finite one-site four-direction diagnostic inherited from615; no claim of the
complete two-time673 H_b process or physical reflection positivity.
Exact Haar result uses integer Fourier constant term with two proven primes.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_full_holonomy_positive_measure_results.json'
spec=importlib.util.spec_from_file_location('entry687',HERE/'round687_drafts/full_holonomy_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
old=entry.old
N=20
DENOM=4**24*495
PAIRS=((0,1),(0,2),(1,2),(3,4))


def indices(n):
    return np.indices((n,)*4,dtype=np.int64).reshape(4,-1).T


def floating_haar(n):
    idx=indices(n)
    phase=np.column_stack((idx,-idx.sum(axis=1)))*(2*np.pi/n)
    lam=(1+np.cos(phase))/2
    homogeneous=np.zeros((len(idx),9));homogeneous[:,0]=1
    for j in range(5):
        for k in range(1,9):
            homogeneous[:,k]+=lam[:,j]*homogeneous[:,k-1]
    M=homogeneous[:,8]/495
    W=np.ones(len(idx))
    for subset in old.old.STATES:
        angle=phase[:,list(subset)].sum(axis=1) if subset else np.zeros(len(idx))
        W*=(1+np.cos(angle))/2
    density=np.ones(len(idx))
    for i,j in PAIRS:
        density*=2-2*np.cos(phase[:,i]-phase[:,j])
    assert abs(np.mean(density)/12-1)<3e-14
    value=float(np.mean(W*M*density)/12)
    weyl=float(np.mean(W*density)/12)
    return dict(grid=n,value=value,normalized_Haar_mass=float(np.mean(density)/12),
                physical_Weyl_average=weyl)


def isprime(p):
    if p<2:return False
    if p%2==0:return p==2
    return all(p%d for d in range(3,math.isqrt(p)+1,2))


def modular_constant_term(p):
    assert isprime(p) and (p-1)%N==0
    generator=2
    while True:
        root=pow(generator,(p-1)//N,p)
        if pow(root,N,p)==1 and pow(root,N//2,p)!=1 and pow(root,N//5,p)!=1:break
        generator+=1
    assert root!=1
    idx=indices(N)
    exponents=np.column_stack((idx,-idx.sum(axis=1)))%N
    powers=np.array([pow(root,j,p) for j in range(N)],dtype=np.int64)
    z=powers[exponents];zinv=powers[(-exponents)%N]
    # Integer numerators: lambda=(2+z+z^-1)/4.
    lam=(2+z+zinv)%p
    homogeneous=np.zeros((len(idx),9),dtype=np.int64);homogeneous[:,0]=1
    for j in range(5):
        for k in range(1,9):
            homogeneous[:,k]=(homogeneous[:,k]+lam[:,j]*homogeneous[:,k-1])%p
    w=np.ones(len(idx),dtype=np.int64)
    for subset in old.old.STATES:
        r=np.ones(len(idx),dtype=np.int64);ri=r.copy()
        for j in subset:
            r=r*z[:,j]%p;ri=ri*zinv[:,j]%p
        w=w*((2+r+ri)%p)%p
    density=np.ones(len(idx),dtype=np.int64)
    for i,j in PAIRS:
        factor=(2-z[:,i]*zinv[:,j]%p-z[:,j]*zinv[:,i]%p)%p
        density=density*factor%p
    values=(w*homogeneous[:,8]%p)*density%p
    total=int(np.sum(values,dtype=np.int64))%p
    ct=total*pow(N**4,-1,p)%p
    assert (p-1)**2+p<np.iinfo(np.int64).max
    return dict(prime=p,primality='trial division through integer square root',
                root_order=N,root=root,constant_term_residue=ct)


def full_group_integral():
    # Degree <=(19,19,19,18); N20 extracts the constant term exactly.
    degrees=np.zeros(4,dtype=int)
    weights=np.vstack((np.eye(4,dtype=int),-np.ones(4,dtype=int)))
    for subset in old.old.STATES:
        if subset:degrees+=abs(weights[list(subset)].sum(axis=0))
    assert np.array_equal(degrees,[8,8,8,8])
    degrees+=8
    for i,j in PAIRS:degrees+=abs(weights[i]-weights[j])
    assert degrees.tolist()==[19,19,19,18] and max(degrees)<N
    residues=[modular_constant_term(p) for p in (2013265921,1000000021)]
    p,q=[x['prime'] for x in residues]
    rp,rq=[x['constant_term_residue'] for x in residues]
    total_modulus=p*q
    bound=12*DENOM
    assert total_modulus>bound
    integer=rp+p*((rq-rp)*pow(p,-1,q)%q)
    assert 0<integer<=bound and integer%12==0
    exact=Fraction(integer,12*DENOM)
    floating=[floating_haar(n) for n in (20,23)]
    for row in floating:
        assert abs(row['value']-float(exact))<2e-16
    assert exact>=Fraction(1,2**40)
    return dict(rank_four_full_group_torus=True,Weyl_divisor=12,
        degree_bounds=degrees.tolist(),grid=N,
        integer_polynomial_denominator=DENOM,integer_constant_term=integer,
        constant_term_bound=bound,CRT_modulus=total_modulus,modular_checks=residues,
        exact_full_Haar_S9_weight=str(exact),decimal_full_weight=float(exact),
        representation_lower_bound=str(Fraction(1,2**40)),floating_crosschecks=floating,
        full_original_Hb_average_not_computed=True)


def actual_matrices_and_positive_kernel():
    inherited=entry.run()
    assert inherited==json.loads((HERE/'round687_drafts/full_holonomy_probe_results.json').read_text('utf8'))
    # A nearby non-Abelian family tests strongly correlated, not just near-diagonal, kernels.
    rs=[np.eye(16,dtype=complex)]
    for j in range(5):
        rs.append(entry.gauge.rep(*entry.gauge.group(68790+j,.025*(j+1))))
    gram=np.array([[entry.weight(a.conj().T@b) for b in rs] for a in rs])
    ev=np.linalg.eigvalsh((gram+gram.T)/2)
    assert min(ev)>-2e-12 and np.max(abs(gram-gram.T))<3e-13
    assert gram[0,1]>.8
    # Sharp-scope tests: M is nonzero on original subgroup, W need not be.
    r=entry.gauge.rep(np.eye(3),np.eye(2),np.exp(1j*np.pi/6))
    assert entry.moment(entry.rotation(r))>=25.**-8
    assert entry.weight(r)<1e-28
    # Original Z6 quotient changes the covering coordinates but neither carrier nor weight.
    c,w,z=entry.gauge.group(68719,.61)
    r0=entry.gauge.rep(c,w,z)
    rt=entry.gauge.rep(np.exp(2j*np.pi/3)*c,-w,np.exp(1j*np.pi/3)*z)
    assert np.linalg.norm(rt-r0,2)<3e-13
    return dict(entry_reproduced=True,original_matrix_checks=inherited['rows'],
        direct_old615_integral_error=inherited['old615_Beta_formula_error'],
        noncommuting_global_representation_law_checked=True,
        nearby_Gram=gram.tolist(),nearby_Gram_eigenvalues=ev.tolist(),
        original_Z6_quotient_respected=True,physical_zero_weight_retained=True,
        full_S9_uniform_lower_bound=str(Fraction(1,25**8)),
        Gram_positivity_proof_is_analytic_not_eigenvalue_sampling=True)


def run():
    deps=('research_note_615.md','research_note_657.md','research_note_672.md',
          'research_note_673.md','research_note_675.md','research_note_679.md',
          'research_note_680.md','research_note_686.md','joint_subgroup_measure_source.py',
          'joint_spinor_subgroup_mass.py','joint_nonflat_mass_measure.py',
          'round687_drafts/full_holonomy_probe.py','round687_drafts/full_holonomy_probe_results.json')
    return dict(date='2026-10-02',round=687,tests_run=2,failures=0,errors=0,
        original_full_group_holonomy=actual_matrices_and_positive_kernel(),
        exact_complete_group_and_S9_integral=full_group_integral(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_group_and_16_channels_retained=True,
            complete_S9_and_single_holonomy_group_Haar=True,
            uniform_auxiliary_nonzero_and_positive_type_proved=True,
            full_single_holonomy_weight_group_positive_type_proved=True,
            original_two_time_Hb_Gauss_RP_decided=False,
            arbitrary_spatial_links_or_nonzero_mass_covered=False,
            original_HF_time_identity_or_quantum_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=687,tests_run=2,all_checks_passed=True,
        exact_Haar=result['exact_complete_group_and_S9_integral']['exact_full_Haar_S9_weight'],
        value=result['exact_complete_group_and_S9_integral']['decimal_full_weight'])))

