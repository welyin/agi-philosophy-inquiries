"""691: actual neutral phase symmetry, graded RP reduction, and odd Gauss probe.

Symmetry implementation on the reflection quotient is not identified with a
thermal left-charge observable, nor an autonomous internal instrument.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'joint_neutral_symmetry_reduction_results.json'
def module(name,path):
    s=importlib.util.spec_from_file_location(name,HERE/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=module('spatial690','round690_drafts/spatial_car_source_probe.py')
cut=module('cut691','round691_drafts/projected_cut_current_probe.py')
old=previous.old;pf=old.internal.pfaffian


def full_covariance(C):
    z=np.zeros_like(C);return np.block([[z,C],[-C.T,z]])


def neutral_os():
    f=previous.physical(np.broadcast_to(np.eye(16,dtype=complex),(2,4,16,16)).copy())
    P=np.zeros((8,128),complex)
    for t,k,s in itertools.product(range(2),repeat=3):
        for x in range(2):P[4*t+2*k+s,32*(2*t+x)+16*s]=(-1)**(k*x)/np.sqrt(2)
    C=P@f['C']@P.T;G=full_covariance(C)
    pos=[4,5,6,7,12,13,14,15];ref=[8,9,10,11,0,1,2,3]
    gram1=G[np.ix_(ref,pos)];a=np.sqrt(5)-2
    weights=np.array([1.,1.,a,a,1.,1.,a,a]);assert np.max(abs(gram1-np.diag(weights)))<3e-12
    monomials=[tuple(i for i in range(8) if mask>>i&1) for mask in range(256)]
    charges=np.array([sum(1 if i<4 else -1 for i in m) for m in monomials])
    norms=np.array([math.prod(weights[list(m)]) for m in monomials])
    worst=0.
    def wick(i,j):
        ids=[ref[k] for k in reversed(monomials[i])]+[pos[k] for k in monomials[j]]
        if len(ids)%2:return 0.
        return pf(G[np.ix_(ids,ids)]) if ids else 1.
    for i in range(256):worst=max(worst,float(abs(wick(i,i)-norms[i])))
    rng=np.random.default_rng(69111)
    for _ in range(120):
        i,j=rng.choice(256,2,replace=False)
        worst=max(worst,float(abs(wick(i,j))))
    assert worst<3e-12
    unit=0.
    for theta in (.2,.71,2*np.pi):
        U=np.exp(1j*theta*charges)
        unit=max(unit,float(np.max(abs(U.conj()*norms*U-norms))))
    dims={str(q):int(sum(charges==q)) for q in range(-4,5)}
    assert [dims[str(q)] for q in range(-4,5)]==[1,8,28,56,70,56,28,8,1]
    assert np.min(norms)>0 and sum(dims.values())==256
    # Local parity is the degree8 monomial: charge0, orthogonal to vacuum.
    assert charges[-1]==0 and abs(wick(0,255))<3e-13
    return dict(actual_neutral_one_field_Gram=weights.tolist(),monomial_dimension=256,
        full_diagonal_and_sampled_offdiagonal_Wick_error=worst,phase_isometry_error=unit,
        charge_sector_dimensions=dims,smallest_monomial_norm=float(min(norms)),
        local_parity_vector_norm_squared=float(norms[-1]),local_parity_vacuum_overlap=0.,
        symmetry_implementer_fixes_vacuum=True,
        thermal_symmetry_implementer_not_identified_with_left_charge_observable=True),f


def factor_and_odd_probe(free):
    rows=[];nu=[32*site+16*s for site in range(4) for s in range(2)]
    rest=[j for j in range(128) if j not in nu]
    for seed in (None,69123,69124):
        links=np.broadcast_to(np.eye(16,dtype=complex),(2,4,16,16)).copy()
        if seed is not None:
            for mu in range(2):
                for i in range(4):links[mu,i]=old.rep(*old.group(seed+10*mu+i,.025))
        f=free if seed is None else previous.physical(links);C=f['C'];G=full_covariance(C)
        off=float(max(np.max(abs(C[np.ix_(nu,rest)])),np.max(abs(C[np.ix_(rest,nu)]))))
        assert off<3e-12
        worst=0.
        for dn,dr in ((2,2),(4,2),(2,4),(4,4),(6,6)):
            # Balanced actual cross-reflection insertions give nonzero values;
            # random monomials would mostly vanish by charge selection.
            kn,kr=dn//2,dr//2
            I=nu[4:4+kn]+[128+i for i in nu[:kn]]
            J=[64+i for i in range(1,kr+1)]+[128+i for i in range(1,kr+1)]
            combined=pf(G[np.ix_(I+J,I+J)])
            product=pf(G[np.ix_(I,I)])*pf(G[np.ix_(J,J)])
            assert abs(product)>1e-7
            worst=max(worst,float(abs(combined-product)))
        assert worst<3e-12
        rows.append(dict(seed=seed,neutral_rest_covariance_cross_error=off,mixed_Wick_factor_error=worst))
    dictionary=old.mass.dictionary;J=dictionary.dictionary()
    B=[((6,10,11),2.),((7,9,11),-2.),((8,9,10),2.)]
    triples=list(itertools.combinations(range(16),3));target={indices:c for indices,c in B}
    error=0.
    for seed in (69151,69152,69153):
        R=J.conj().T@old.rep(*old.group(seed,.6))@J
        for cols in triples:
            value=sum(c*np.linalg.det(R[np.ix_(indices,cols)]) for indices,c in B)
            error=max(error,float(abs(value-target.get(cols,0))))
    assert error<3e-12
    # Actual source rows for B=u^c d^c d^c at positive time, x0, spin0.
    G=full_covariance(free['C']);theta=np.zeros((256,256))
    for site in range(4):
        t,x=divmod(site,2);rs=2*(1-t)+x
        for comp in range(32):
            theta[32*site+comp,128+32*rs+comp]=1
            theta[128+32*site+comp,32*rs+comp]=1
    fields=[]
    for i in range(16):
        z=np.zeros(256,complex);z[64:80]=J.conj().T[i,:];fields.append(z)
    neutral=np.zeros(256,complex);neutral[64]=neutral[96]=1/np.sqrt(2)
    bare=0j;even=0j
    for ia,ca in B:
        for ib,cb in B:
            left=[fields[k] for k in ia];right=[fields[k] for k in ib]
            z=np.array([r.conj()@theta for r in reversed(left)]+right)
            bare+=ca*cb*pf(z@G@z.T)
            z=np.array([r.conj()@theta for r in reversed([neutral]+left)]+[neutral]+right)
            even+=ca*cb*pf(z@G@z.T)
    a=np.sqrt(5)-2;expected=12*((1+a)/2)**3
    assert abs(bare-expected)<3e-12 and abs(even-bare)<3e-12
    return dict(actual_factor_rows=rows,odd_Gauss_polynomial='epsilon^(abc) u^c_a d^c_b d^c_c',
        local_actual_G_invariance_error=error,odd_polynomial_nonzero_terms=len(B),
        original_free_odd_baryon_reflected_norm=float(bare.real),
        original_free_even_neutral_baryon_reflected_norm=float(even.real),
        graded_factorization_error=float(abs(even-bare)),
        original_dynamic_rest_RP_not_inferred_from_free_sample=True)


def run():
    entry=cut.run();assert entry==json.loads(cut.TARGET.read_text('utf8'))
    os,f=neutral_os();factor=factor_and_odd_probe(f)
    deps=('research_note_661.md','research_note_670.md','research_note_673.md','research_note_679.md','research_note_680.md',
          'research_note_689.md','research_note_690.md','joint_spatial_charge_dictionary.py','joint_spatial_charge_dictionary_results.json',
          'round691_drafts/projected_cut_current_probe.py','round691_drafts/projected_cut_current_probe_results.json','round691_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=691,tests_run=2,failures=0,errors=0,neutral_OS_symmetry=os,
        original_graded_factorization=factor,
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_dynamic_massless_RP_equivalent_to_rest_functional_RP=True,
            original_auxiliary16_channels_measure_and_normalization_retained=True,
            even_full_algebra_still_requires_odd_rest_tests=True,
            neutral_phase_unitary_on_reflection_quotient=True,
            neutral_phase_implementer_not_proved_as_same_local_charge_or_instrument=True,
            general_dynamic_Q0_RP_decided=False,full_original_HF_or_common_continuum_or_GR_completed=False,
            old_space_and_full_goal_unchanged=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();out=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    else:assert out==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=691,all_checks_passed=True,
        neutral_quotient_dimension=out['neutral_OS_symmetry']['monomial_dimension'],
        graded_factorization_error=out['original_graded_factorization']['graded_factorization_error'])))
