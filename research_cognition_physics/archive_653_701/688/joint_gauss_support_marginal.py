"""688: exact Gauss multiplicities and a CAR-only marginal witness.

The spherical kernel, its spectrum and temporal transfer are inherited from653,
not new results. Current new calculations use the same actual 32-CAR character.
No claim identifies this conditional branch with original full H_F or physical Y.
"""
import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import joint_full_holonomy_positive_measure as old
import joint_chiral_character_state as original
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_support_marginal_results.json'
spec=importlib.util.spec_from_file_location('entry688',HERE/'round688_drafts/temporal_sphere_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)
GRID=20


def modular_characters(p):
    assert old.isprime(p) and (p-1)%GRID==0
    gen=2
    while True:
        root=pow(gen,(p-1)//GRID,p)
        if pow(root,GRID,p)==1 and pow(root,10,p)!=1 and pow(root,4,p)!=1:break
        gen+=1
    idx=old.indices(GRID)
    exps=np.column_stack((idx,-idx.sum(axis=1)))%GRID
    power=np.array([pow(root,k,p) for k in range(GRID)],dtype=np.int64)
    z=power[exps];zi=power[-exps%GRID]
    sym=np.zeros((len(idx),9),dtype=np.int64);sym[:,0]=1
    for eigenvalue in np.concatenate((z,zi),axis=1).T:
        for ell in range(1,9):
            sym[:,ell]=(sym[:,ell]+eigenvalue*sym[:,ell-1])%p
    harmonic=sym.copy();harmonic[:,2:]=(sym[:,2:]-sym[:,:-2])%p
    F=np.ones(len(idx),dtype=np.int64)
    for subset in old.old.old.STATES:
        q=np.ones(len(idx),dtype=np.int64);qi=q.copy()
        for j in subset:q=q*z[:,j]%p;qi=qi*zi[:,j]%p
        F=F*((2+q+qi)%p)%p
    density=np.ones(len(idx),dtype=np.int64)
    for i,j in old.PAIRS:
        density=density*((2-z[:,i]*zi[:,j]%p-z[:,j]*zi[:,i]%p)%p)%p
    normalizer=pow(GRID**4,-1,p)
    fd=F*density%p
    total=[];alone=[]
    for ell in range(9):
        total.append(int(np.sum(fd*harmonic[:,ell]%p,dtype=np.int64))%p*normalizer%p)
        alone.append(int(np.sum(density*harmonic[:,ell]%p,dtype=np.int64))%p*normalizer%p)
    assert (p-1)**2+p<np.iinfo(np.int64).max
    return dict(prime=p,root=root,CT_full=total,CT_aux_only=alone)


def exact_gauss_support():
    # Explicitly reuse old653 integer spectrum and dimension, rather than deriving it again.
    saved=json.loads((HERE/'joint_chiral_character_state_results.json').read_text('utf8'))
    previous=saved['temporal_gluing']['spherical_harmonic_spectrum']
    modes=entry.eigenvalues()
    assert [(r['l'],r['dimension'],Fraction(r['eigenvalue'])) for r in previous]==modes
    residue=[modular_characters(p) for p in (2013265921,1000000021)]
    p,q=[r['prime'] for r in residue];modulus=p*q
    rows=[]
    for ell,dim,mu in modes:
        def crt(a,b):return a+p*((b-a)*pow(p,-1,q)%q)
        C=crt(residue[0]['CT_full'][ell],residue[1]['CT_full'][ell])
        S=crt(residue[0]['CT_aux_only'][ell],residue[1]['CT_aux_only'][ell])
        bound=12*2**32*dim
        assert modulus>bound and 0<=C<=bound and C%12==0
        assert 0<=S<=12*dim and S%12==0
        m,s=C//12,S//12
        assert s==int(ell%2==0)
        rows.append(dict(ell=ell,dimension=dim,mu=str(mu),
            Gauss_multiplicity=m,auxiliary_singlets=s,
            constant_term=C,CT_bound=bound))
    n0=rows[0]['Gauss_multiplicity']
    assert n0>0 and rows[1]['Gauss_multiplicity']>0
    assert all(r['Gauss_multiplicity']>=n0*r['auxiliary_singlets'] for r in rows)
    Z1=sum(r['Gauss_multiplicity']*Fraction(r['mu']) for r in rows)/2**32
    reference=Fraction(json.loads((HERE/'joint_full_holonomy_positive_measure_results.json').read_text('utf8'))['exact_complete_group_and_S9_integral']['exact_full_Haar_S9_weight'])
    assert Z1==reference
    return dict(original_653_spectrum_reused_exactly=True,grid=GRID,
        degree_bounds=[19,19,19,18],modular_residues=residue,CRT_modulus=modulus,
        rows=rows,original_CAR_singlet_dimension=n0,
        exact_round687_weight_recovered=str(Z1))


def marginal_certificate(support):
    rows=support['rows'];m=[r['Gauss_multiplicity'] for r in rows]
    singlets=[r['auxiliary_singlets'] for r in rows]
    mu=[Fraction(r['mu']) for r in rows];n0=m[0]
    results=[]
    for N in (1,2,3,4,8,16):
        total=sum(count*x**N for count,x in zip(m,mu))
        already_singlet=n0*sum(count*x**N for count,x in zip(singlets,mu))
        p=already_singlet/total;distance=1-p
        assert 0<p<1
        residual=sum((m[i]-n0*singlets[i])*(mu[i]/mu[0])**N for i in range(1,9))
        assert residual>0
        assert distance==residual/sum(m[i]*(mu[i]/mu[0])**N for i in range(9))
        upper=sum(m[i]-n0*singlets[i] for i in range(1,9))/n0*float(Fraction(8,17)**N)
        assert float(distance)<=upper+1e-16
        results.append(dict(N=N,raw_projected_partition=str(total/Fraction(2**32)**N),
            CAR_singlet_probability=str(p),probability_float=float(p),
            trace_distance_to_CAR_only_Gauss_state=str(distance),distance_float=float(distance),
            analytic_upper_bound=float(upper)))
    assert results[-1]['distance_float']<results[0]['distance_float']
    # One physical CAR character factor: singlet multiplicity always a multiple of4,
    # since the original neutral two spin modes are gauge singlets.
    assert all(count%4==0 for count in m)
    return dict(rows=results,neutral_two_mode_factor_four_retained=True,
        witness='Pi_F = Haar projection of the original 32-CAR representation',
        original_E_measurement_not_needed_for_this_reconstructed_CAR_witness=True,
        original_Weyl_Y_equal_time_insertion_dictionary_not_proved=True,
        comparison_state='Pi_F / dim(F^G), conditional CAR-only zero-H reference',
        full_original_HF_state_not_substituted=True)


def run():
    s=exact_gauss_support()
    # Entry is preserved as repeated653 validation, not a new scientific group.
    inherited=entry.run()
    assert inherited==json.loads((HERE/'round688_drafts/temporal_sphere_probe_results.json').read_text('utf8'))
    deps=('research_note_653.md','research_note_654.md','research_note_661.md','research_note_687.md',
          'joint_chiral_character_state.py','joint_chiral_character_state_results.json',
          'joint_full_holonomy_positive_measure.py','joint_full_holonomy_positive_measure_results.json',
          'round688_drafts/temporal_sphere_probe.py','round688_drafts/temporal_sphere_probe_results.json',
          'round688_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=688,tests_run=2,failures=0,errors=0,
        exact_projected_support=s,CAR_marginal=marginal_certificate(s),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(new_spherical_transfer_theorem_claimed=False,
            previous687_novelty_scope_corrected=True,
            original_Gauss_multiplicities_exact=True,
            reconstructed_CAR_marginal_witness_exact=True,
            entire_original_dynamic_Hb_Gauss_RP_decided=False,
            original_Weyl_Y_source_identity_proved=False,
            original_HF_or_quantum_GR_completed=False,old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=688,all_checks_passed=True,
        multiplicities=[r['Gauss_multiplicity'] for r in result['exact_projected_support']['rows']],
        distances=[r['distance_float'] for r in result['CAR_marginal']['rows']])))
