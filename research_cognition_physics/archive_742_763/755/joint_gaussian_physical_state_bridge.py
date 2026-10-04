"""755: fixed original CAR/Gauss state identity versus source-only matching.

The global centre rule is inherited from the round710 entry. New checks concern
the proposed Gaussian state bridge, the round754 pair contract, and a full
32-mode physical alternative with the same bilinears but different records.
No Fock cutoff or new Hamiltonian is introduced. Small Fock matrices only
independently check the variance identity; the original 32-mode states are sparse.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_reference_constraint_strata as background

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_gaussian_physical_state_bridge_results.json'


def residue_weights(probabilities):
    weights = [F(1), F(0), F(0)]
    for p in probabilities:
        weights = [(1-p)*weights[r]+p*weights[(r-1)%3] for r in range(3)]
    return weights


def residue_check():
    checked = 0
    for m in range(1, 6):
        for probs in itertools.product((F(0), F(1), F(1,3)), repeat=m):
            w = residue_weights(probs)
            pure_neutral = all(p in (0,1) for p in probs) and sum(probs)%3 == 0
            assert (w[0] == 1) == pure_neutral and sum(w) == 1
            checked += 1
    m = 24
    w = residue_weights([F(1,2)]*m)
    assert w[0] == (1+F(1,2)**(m-1))/3
    omega = np.exp(2j*np.pi/3)
    formula_error = 0.
    for probs in ([F(1,2)]*24, [F(0),F(1),F(2,7),F(3,8)], [F(1)]*6):
        product = np.prod([1-float(p)+float(p)*omega for p in probs])
        formula_error = max(formula_error, abs((1+2*product.real)/3-float(residue_weights(probs)[0])))
    assert formula_error < 2e-15
    return dict(exact_small_spectral_cases=checked,
                Gaussian_quark_modes=m, half_occupied_physical_center_weight=str(w[0]),
                half_occupied_trace_distance_lower_bound=str(1-w[0]),
                lower_bound_float=float(1-w[0]), characteristic_formula_error=formula_error,
                full_Gauss_support_is_stronger_than_center_support=True,
                no_actual_Gibbs_covariance_computed=True)


def pair_check():
    rows=[]
    for label,probs in (
            ('vacuum_rest',[F(0)]*22),
            ('half_filled_rest',[F(1,2)]*22),
            ('sharp_lower_bound',[F(1,2)]+[F(0)]*21)):
        w=residue_weights(probs)
        minus=residue_weights([F(0),F(0)]+probs)[0]
        plus=residue_weights([F(1),F(1)]+probs)[0]
        assert minus==w[0] and plus==w[1] and minus+plus<=1
        lower=max(1-minus,1-plus)
        assert lower>=F(1,2)
        rows.append(dict(rest=label,minus_center_weight=str(minus),
                         plus_center_weight=str(plus),at_least_one_distance=str(lower)))
    assert rows[-1]['at_least_one_distance']=='1/2'
    return dict(rows=rows, same_complement_pair_contract_excluded=True,
                bosonic_only_dressing_cannot_change_global_center=True,
                finite_original_CAR_mapping_assumption_explicit=True,
                no_infinite_volume_number_operator_assumed=True)


def annihilate(state,j):
    out={}
    for mask,a in state.items():
        if (mask>>j)&1:
            out[mask^(1<<j)]=a*((-1)**((mask&((1<<j)-1)).bit_count()))
    return out


def create(state,j):
    out={}
    for mask,a in state.items():
        if not ((mask>>j)&1):
            out[mask|(1<<j)]=a*((-1)**((mask&((1<<j)-1)).bit_count()))
    return out


def inner(a,b):
    return sum(complex(x).conjugate()*b.get(mask,0) for mask,x in a.items())


def second_quantize(t,state):
    out={}
    for i,j in zip(*np.nonzero(t)):
        term=create(annihilate(state,int(j)),int(i))
        for mask,a in term.items():out[mask]=out.get(mask,0)+t[i,j]*a
    return out


def color_generators32():
    generators=[]
    for t in background.T:
        g=np.zeros((32,32),complex)
        g[:12,:12]=np.kron(t,np.eye(4))
        g[12:18,12:18]=np.kron(t,np.eye(2))
        g[18:24,18:24]=np.kron(t,np.eye(2))
        generators.append(g)
    return generators


def sparse_covariance(state):
    C=np.zeros((32,32),complex);A=np.zeros_like(C)
    for i in range(32):
        for j in range(32):
            C[i,j]=inner(state,create(annihilate(state,j),i))
            A[i,j]=inner(state,annihilate(annihilate(state,j),i))
    return C,A


def small_fock_variance_check():
    c=[]
    for j in range(3):
        a=np.zeros((8,8),complex)
        for m in range(8):
            for out,x in annihilate({m:1},j).items():a[out,m]=x
        c.append(a)
    total=0.
    for t in background.T:
        q=sum((t[i,j]*c[i].conj().T@c[j] for i in range(3) for j in range(3)),np.zeros((8,8),complex))
        total+=np.trace(q@q).real/8
    assert abs(total-1)<2e-15  # p(1-p)*sum_a tr(T_a^2)=1 for p=1/2.
    return float(total)


def same_bilinears_physical_record_check():
    empty={0:1.};full={(1<<32)-1:1.}
    C0,A0=sparse_covariance(empty);C1,A1=sparse_covariance(full)
    assert np.array_equal(C0,np.zeros((32,32)))
    assert np.array_equal(C1,np.eye(32)) and not np.any(A0) and not np.any(A1)
    generators=color_generators32()
    y=np.diag([1]*12+[4]*6+[-2]*6+[-3]*4+[-6]*2+[0]*2)
    # Full local invariance follows from det R=1. Independently check color
    # and hypercharge generators on both sparse original Fock vectors.
    invariant_error=max(abs(inner(v,second_quantize(g,v))) for v in (empty,full) for g in [*generators,y])
    norm_error=max(float(sum(abs(x)**2 for x in second_quantize(g,v).values()))
                   for v in (empty,full) for g in [*generators,y])
    assert invariant_error<2e-13 and norm_error<2e-25 and np.trace(y)==0
    rng=np.random.default_rng(755)
    determinant_error=0.;source_error=0.;center_error=0.
    phase=np.exp(2j*np.pi/3)
    R=matter.representation(phase*np.eye(3),np.eye(2),1)
    center_error=float(np.max(abs(R-np.diag([phase]*24+[1]*8))))
    for _ in range(5):
        R=matter.representation(matter.gauge.group_exp(rng.normal(size=8),3),
                               matter.gauge.group_exp(rng.normal(size=3),2),np.exp(1j*rng.normal()))
        determinant_error=max(determinant_error,float(abs(np.linalg.det(R)-1)))
        phi=rng.normal(size=5)*.2;h,delta=matter.mass_matrices(phi)
        for p in (.2,.5,.8):
            # Compare sparse Fock evaluation to the Gaussian two-point rule
            # on the original matrix mass and independent full source symbols.
            for t in [h,y,*generators]:
                lhs=(1-p)*inner(empty,second_quantize(t,empty))+p*inner(full,second_quantize(t,full))
                source_error=max(source_error,float(abs(lhs-p*np.trace(t))))
            assert np.max(abs((1-p)*A0+p*A1))==0
            assert np.max(abs((1-p)*C0+p*C1-p*np.eye(32)))==0
            # All pair expectations vanish exactly, including original Majorana.
            assert abs(np.sum(delta*((1-p)*A0+p*A1)))==0
    assert max(determinant_error,source_error,center_error)<1e-12
    traces=np.array([np.trace(g@g).real for g in generators])
    assert np.max(abs(traces-4))<2e-15
    rows=[]
    for p in (F(1,5),F(1,2),F(4,5)):
        actual=[1-p,F(0),F(0),p]
        gaussian=[(1-p)**2,p*(1-p),p*(1-p),p*p]
        tv=sum(abs(a-b) for a,b in zip(actual,gaussian))/2
        assert tv==2*p*(1-p)
        rows.append(dict(p=str(p),physical_two_sterile_record=[str(x) for x in actual],
                         Gaussian_two_sterile_record=[str(x) for x in gaussian],
                         record_total_variation=str(tv),four_point_Wick_defect=str(p*(1-p)),
                         physical_color_Gauss_square='0',
                         Gaussian_color_Gauss_square=str(32*p*(1-p))))
    return dict(original_modes=32,sparse_Fock_support_per_vector=1,
                full_local_group_determinant_error=determinant_error,
                generator_annihilation_norm_square_error=norm_error,
                center_representation_error=center_error,instantaneous_bilinear_source_error=source_error,
                small_Fock_independent_variance=small_fock_variance_check(),rows=rows,
                full_original_H_preserves_Gauss_not_the_two_vector_code=True,
                no_753_geometry_or_continuum_completion_claim=True,
                no_autonomous_instrument_implementation_claim=True)


def run():
    a=residue_check();b=pair_check();c=same_bilinears_physical_record_check()
    deps=('research_note_709.md','round710_drafts/residual_phase_entry.md','research_note_728.md',
          'research_note_730.md','research_note_742.md','research_note_751.md','research_note_753.md',
          'research_note_754.md','joint_fermion_gauss_completion.py','joint_reference_constraint_strata.py')
    return dict(round=755,tests_run=3,failures=0,errors=0,
                Gaussian_center_support=a,unchanged_complement_pair=b,
                physical_source_matching_and_records=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Necessary qualification of a finite original closed Gauss CAR bridge, not a new centre theorem. Number-conserving quasifree quark marginals require a neutral pure Slater symbol; fixed-complement two-mode filling cannot both lift faithfully. A full original 32-mode local singlet family matches all instantaneous bilinears but changes actual joint sterile records. No continuum number operator, new action, quantum Einstein constraint solution, autonomous preparation, or exact feedback is asserted.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
