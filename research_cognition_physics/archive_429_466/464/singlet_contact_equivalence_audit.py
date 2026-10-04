"""Round 464: complete contact-weight ambiguity inside a shared singlet sector.

This is an operational-interface boundary, not different physical universes
with the same complete physics. An entirely internal extension distinguishes
the previously equivalent generators. No probability threshold defines edges.
"""
import argparse
from fractions import Fraction
import io
import itertools
import json
import math
from pathlib import Path
import platform
import unittest
import numpy as np
import exchange_relation_audit as old
from local_descriptor_support_audit import exact_rank
from recursive_exchange_interface_audit import permutation

TARGET=Path(__file__).with_name('singlet_contact_equivalence_audit_results.json')
OBS={}
CACHE={}


def singlet_projector(n):
    if n not in CACHE:
        pairs=list(itertools.combinations(range(n),2))
        swaps=[old.swap(n,*e).real.astype(np.int64) for e in pairs]
        eye=np.eye(2**n,dtype=np.int64)
        casimir=sum(swaps)-n*(n-4)//4*eye
        numerator=eye.copy()
        denominator=1
        for j in range(1,n//2+1):
            factor=j*(j+1)
            numerator=numerator@(factor*eye-casimir)
            denominator*=factor
        CACHE[n]=(pairs,swaps,numerator,denominator,casimir)
    return CACHE[n]


def gram_data(n):
    pairs,swaps,q,den,_=singlet_projector(n)
    dim=math.comb(n,n//2)-math.comb(n,n//2-1)
    mean=Fraction(n-4,2*(n-1))
    perms=[permutation(n,*e) for e in pairs]
    gram=[]
    for pe in perms:
        row=[]
        for pf in perms:
            tr=int(np.sum(q[np.arange(2**n),pe[pf]]))
            row.append(Fraction(tr,den*dim)-mean**2)
        gram.append(row)
    return pairs,gram,mean,dim


def encoding(n):
    _,_,q,den,_=singlet_projector(n)
    values,vectors=np.linalg.eigh(q/den)
    return vectors[:,values>.5]


class Audit(unittest.TestCase):
    def close(self,a,b,tolerance=2e-11):
        self.assertLess(float(np.linalg.norm(a-b)),tolerance)

    def test_01_exact_complete_singlet_projectors(self):
        rows=[]
        for n in (4,6,8):
            _,_,q,den,c=singlet_projector(n)
            self.assertTrue(np.array_equal(q@q,den*q))
            self.assertFalse(np.any(c@q))
            dim=math.comb(n,n//2)-math.comb(n,n//2-1)
            self.assertEqual(int(np.trace(q)),den*dim)
            rows.append(dict(n=n,sector_dimension=dim,projector_denominator=den))
        OBS['sector_projectors']=rows

    def test_02_full_exact_gram_and_kernel_rank(self):
        rows=[]
        for n in (4,6,8):
            pairs,gram,s,dim=gram_data(n)
            a=1-s*s
            b=-a/(n-2)
            c=2*a/((n-2)*(n-3))
            for i,e in enumerate(pairs):
                for j,f in enumerate(pairs):
                    expected=a if e==f else (b if set(e)&set(f) else c)
                    self.assertEqual(gram[i][j],expected)
            scale=math.lcm(*(v.denominator for row in gram for v in row))
            integer=[[int(v*scale) for v in row] for row in gram]
            rank=exact_rank(integer)
            self.assertEqual(rank,n*(n-3)//2)
            incidence=np.array([[int(v in e) for e in pairs] for v in range(n)])
            self.assertEqual(exact_rank(incidence),n)
            self.assertFalse(np.any(np.array(integer,dtype=np.int64)@incidence.T))
            eigenvalue=a*(n-1)/(n-3)
            self.assertGreater(eigenvalue,0)
            rows.append(dict(n=n,edges=len(pairs),sector_dimension=dim,
                diagonal=str(a),shared_endpoint=str(b),disjoint=str(c),
                positive_eigenvalue=str(eigenvalue),rank=rank,kernel_dimension=n,
                exact_integer_gram_scale=scale))
        OBS['complete_weight_kernel']=dict(
            criterion='Jprime_ij-J_ij = a_i+a_j',
            scope='all even N>=4; scalar equality on the whole total-j=0 sector',
            certificates=rows,not_just_one_stationary_state=True)

    def test_03_each_local_star_is_scalar_on_the_whole_sector(self):
        for n in (4,6,8):
            pairs,swaps,q,den,_=singlet_projector(n)
            for i in range(n):
                star=sum(s for e,s in zip(pairs,swaps) if i in e)
                self.assertTrue(np.array_equal(star@q,(n-4)//2*q))
            weights=np.arange(1,n+1)
            shift=sum(int(weights[i]+weights[j])*s for (i,j),s in zip(pairs,swaps))
            scalar=(n-4)//2*int(weights.sum())
            self.assertTrue(np.array_equal(shift@q,scalar*q))
        OBS['star_identity']=dict(
            local_identity='sum_(j!=i) P0 Sij P0 = (N-4) P0/2',
            total_scalar_shift='(N-4) sum_i a_i / 2',
            independent_raw_star_directions=True,
            uniform_all_pair_shift_is_only_one_special_case=True)

    def test_04_positive_sparse_dense_weights_and_complete_internal_process(self):
        n=6
        pairs,swaps,_,_,_=singlet_projector(n)
        v=encoding(n)
        j=np.array([int(b==a+1) for a,b in pairs])
        add=np.array([a+b+2 for a,b in pairs])
        h=sum(int(w)*s for w,s in zip(j,swaps))
        hp=sum(int(w)*s for w,s in zip(j+add,swaps))
        a,b=v.T@h@v,v.T@hp@v
        scalar=21
        self.close(b-a,scalar*np.eye(5))
        rng=np.random.default_rng(46404)
        rho=old.density(10,rng)
        raw_pair=v.T@swaps[0]@v
        projectors=[(np.eye(5)+raw_pair)/2,(np.eye(5)-raw_pair)/2]
        # An actual instrument with retained outcome labels, followed by more
        # evolution: both subnormalized reference-carrying branches agree.
        first=old.evolve(a,.17)
        second=old.evolve(b,.17)
        x=np.kron(first,np.eye(2))@rho@np.kron(first.conj().T,np.eye(2))
        y=np.kron(second,np.eye(2))@rho@np.kron(second.conj().T,np.eye(2))
        probabilities=[]
        for p in projectors:
            k=np.kron(p,np.eye(2))
            ux=np.kron(old.evolve(a,.29),np.eye(2))
            uy=np.kron(old.evolve(b,.29),np.eye(2))
            xx=ux@k@x@k@ux.conj().T
            yy=uy@k@y@k@uy.conj().T
            self.close(xx,yy)
            probabilities.append(float(f'{np.trace(xx).real:.12g}'))
        self.assertAlmostEqual(sum(probabilities),1)
        OBS['internal_process_equivalence']=dict(n=6,old_contacts=int(np.count_nonzero(j)),
            new_contacts=int(np.count_nonzero(j+add)),old_support='P6',new_support='K6',
            all_present_weights_positive=True,scalar_shift=scalar,
            retained_instrument_probabilities=probabilities,
            arbitrary_unknown_sector_state_and_reference_covered_analytically=True,
            arbitrary_sector_preserving_adaptive_instruments_covered_analytically=True,
            old_and_new_physical_universes_claimed_distinct=False)

    def test_05_new_internal_relation_distinguishes_old_equivalence(self):
        n=6
        s=np.array([0,1,-1,0],dtype=np.int64)
        psi=np.kron(np.kron(s,s),s)  # norm squared 8
        self.assertEqual(int(psi@psi),8)
        _,_,q,den,casimir=singlet_projector(n)
        self.assertFalse(np.any(casimir@psi))
        raw=lambda a,b: old.swap(n,a,b).real.astype(np.int64)
        h0=raw(0,4)
        star=raw(0,1)+raw(0,2)+raw(0,3)
        h1=star+h0
        em= np.eye(64,dtype=np.int64)-raw(0,4)  # 2E
        k=h1-np.eye(64,dtype=np.int64)
        self.assertFalse(np.any(star@psi))
        self.assertTrue(np.array_equal(k@k@psi,psi))
        self.assertEqual(int(psi@em@psi),4)
        self.assertEqual(int((k@psi)@em@(k@psi)),16)
        self.assertFalse(np.any(h0@em-em@h0))
        self.assertFalse(np.any(casimir@h1-h1@casimir))
        rows=[]
        for t in (0.,.31,math.pi/2,1.8):
            v0=old.evolve(h0,t)@psi/math.sqrt(8)
            v1=old.evolve(h1,t)@psi/math.sqrt(8)
            p0=np.vdot(v0,em@v0).real/2
            p1=np.vdot(v1,em@v1).real/2
            self.assertAlmostEqual(p0,.25)
            self.assertAlmostEqual(p1,.25+.75*math.sin(t)**2)
            self.close(casimir@v1,np.zeros(64))
            rows.append(dict(time=float(f'{t:.12g}'),p0=float(f'{p0:.12g}'),p1=float(f'{p1:.12g}')))
        OBS['internal_extension']=dict(
            old_subjects='0,1,2,3 with all total-j=0 inputs',
            old_H0='0',old_H1='S01+S02+S03',added_common_contact='S04',
            new_internal_pair='singlet 45',initial='singlet 01 tensor singlet 23 tensor singlet 45',
            effect='singlet 04 projector',p_H0='1/4',p_H1='1/4 + (3/4) sin(t)^2',
            center_difference='3/4',center='pi/2',window_half_width='1/4',
            difference_on_window_lower='45/64',entire_six_body_total_j_zero_preserved=True,
            old_four_body_total_j_zero_not_preserved=True,external_axis_or_reference=False,
            numerical_reproduction=rows)

    def test_06_canonical_representative_does_not_select_positive_locality(self):
        extension_rows=[]
        for old_n in (4,6):
            pairs,gram,_,_=gram_data(old_n+2)
            indices=[i for i,(a,b) in enumerate(pairs) if b<old_n]
            restricted=[[gram[i][j] for j in indices] for i in indices]
            scale=math.lcm(*(x.denominator for row in restricted for x in row))
            integer=[[int(x*scale) for x in row] for row in restricted]
            self.assertEqual(exact_rank(integer),math.comb(old_n,2))
            constant=Fraction(3*(old_n+4),2*(old_n+1)**2*(old_n-1))
            standard=constant*(old_n+1)
            cycle=constant*old_n*(old_n+1)/2
            expected=[float(constant)]+[float(standard)]*(old_n-1)+[float(cycle)]*(old_n*(old_n-3)//2)
            actual=np.linalg.eigvalsh(np.array(integer,dtype=float)/scale)
            self.close(actual,np.array(expected))
            extension_rows.append(dict(old_n=old_n,new_n=old_n+2,
                old_weight_count=len(indices),restricted_gram_rank=len(indices),
                smallest_eigenvalue=str(constant),
                other_eigenvalues=[str(standard),str(cycle)]))
        OBS['complete_internal_reference_repair']=dict(
            internal_qubits_added=2,new_total_j_zero_interface=True,
            newly_added_contact_weights_equal_between_candidates=True,
            complete_old_weight_kernel_removed=True,
            minimum_gram_eigenvalue='3(N+4)/(2(N+1)^2(N-1))',
            inverse_bound='Euclidean norm(delta J) <= sector RMS error / sqrt(lambda_min)',
            full_new_sector_access_is_extra_contract=True,
            spontaneous_preparation_or_control_not_derived=True,
            certificates=extension_rows)
        n=6
        pairs,_,_,_,_=singlet_projector(n)
        incidence=np.array([[int(i in e) for e in pairs] for i in range(n)],dtype=float)
        j=np.array([int(b==a+1) for a,b in pairs],dtype=float)
        q=np.eye(len(pairs))-incidence.T@np.linalg.solve(incidence@incidence.T,incidence)
        canonical=q@j
        self.close(incidence@canonical,np.zeros(n))
        self.close(q@q,q)
        self.close(q@(j+incidence.T@np.arange(n)),canonical)
        self.assertLess(canonical.min(),0)
        self.assertGreater(canonical.max(),0)
        OBS['quotient_is_not_space']=dict(
            unique_zero_row_sum_representative=True,
            minimum_Euclidean_weight_norm_only=True,
            canonical_positive_or_sparse_not_guaranteed=True,
            example_canonical_min=float(f'{canonical.min():.12g}'),
            example_canonical_max=float(f'{canonical.max():.12g}'),
            sign_or_minimum_support_not_derived_cognitive_axioms=True,
            quotient_weight_dimension_not_spatial_dimension=True,
            actual_preparations_and_access_permissions_remain_inputs=True)


def run():
    OBS.clear()
    out=io.StringIO()
    result=unittest.TextTestRunner(stream=out).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not result.wasSuccessful():
        raise RuntimeError(out.getvalue())
    return dict(round=464,baseline_round=463,tests_run=result.testsRun,
        failures=len(result.failures),errors=len(result.errors),
        python=platform.python_version(),numpy=np.__version__,observations=OBS,
        scope=dict(complete_shared_singlet_process_weight_kernel_classified=True,
            singlet_only_interface_is_extra_assumption=True,
            global_correlations_derive_contact_graph=False,
            internal_extension_distinguishes_previous_equivalence=True,
            full_cognitive_principles_refuted=False,spatial_dimension_derived=False,
            full_GR_goal_completed=False,phase_closure_triggered=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dry-run',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.check:
        assert result==json.loads(TARGET.read_text(encoding='utf-8'))
    elif not args.dry_run:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
