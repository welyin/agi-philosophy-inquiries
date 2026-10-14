"""Independent exact review for 1182; author code is never imported.

Three groups: center/fusion, normalized Peter--Weyl fundamental transport,
and admissibility/deletion controls. --probe-only does not bind draft files;
--write creates the frozen review once; default only recomputes and compares.
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction as Q
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def algebra():
    path = HERE/'reviewer_algebra.py'
    spec = importlib.util.spec_from_file_location('review1182_exact_algebra', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fuse(spins2):
    """Exact irrep multiplicities by SU(2) tensor product, twice-spin labels."""
    multiplicities = {0: 1}
    for spin in spins2:
        out = defaultdict(int)
        for previous, count in multiplicities.items():
            for target in range(abs(previous-spin), previous+spin+1, 2):
                out[target] += count
        multiplicities = dict(out)
    assert sum((j+1)*n for j,n in multiplicities.items()) == __import__('math').prod(j+1 for j in spins2)
    return multiplicities


def allowed(k2, gamma):
    return all(((1+s*gamma)*Q(k2,2)).denominator == 1 for s in (-1,1))


def core_math():
    a = algebra()
    dimension_cases = {
        'four_original_legs': (3,3,3,3),
        'plus_one_spinor': (3,3,3,3,1),
        'plus_spinor_and_half_spin_reference': (3,3,3,3,1,1),
        'one_leg_to_spin_1_plus_spinor': (2,3,3,3,1),
        'one_leg_to_spin_2_plus_spinor': (4,3,3,3,1),
    }
    fusion = {label:fuse(spins) for label,spins in dimension_cases.items()}
    dimensions = {label:multiplicities.get(0,0) for label,multiplicities in fusion.items()}
    assert dimensions['four_original_legs'] == 4
    assert dimensions['plus_one_spinor'] == 0
    assert dimensions['plus_spinor_and_half_spin_reference'] == 13
    assert dimensions['one_leg_to_spin_1_plus_spinor'] == 6
    assert dimensions['one_leg_to_spin_2_plus_spinor'] == 7
    assert fusion['four_original_legs'][2] == 9
    assert dimensions['plus_spinor_and_half_spin_reference'] == 4+9
    assert (-1)**sum(dimension_cases['plus_one_spinor']) == -1
    # Directly recheck that the added matter/reference singlet has norm one
    # and total raising/lowering annihilate it, rather than treating the
    # existence of a neutral pair as an unknown-qubit encoding certificate.
    singlet = {(m,n):a.cg(1,m,1,n,0,0) for m in (-1,1) for n in (-1,1)}
    assert sum((c*c for c in singlet.values()),a.R()) == a.R(1)
    assert singlet[(1,-1)] + singlet[(-1,1)] == a.R()

    # f^j_mn = sqrt(2j+1) D^j_mn, Haar volume one. Multiplication by the
    # fundamental SU(2) matrix acts on the retained two-component spinor.
    # (j=3/2 Peter--Weyl sector) x spinor: 32 columns; j'=1,2: 68 rows.
    ins = [(m,n,b) for m in range(-3,4,2) for n in range(-3,4,2) for b in (-1,1)]
    outs = [(j,m,n,aa) for j in (2,4) for m in range(-j,j+1,2)
            for n in range(-j,j+1,2) for aa in (-1,1)]
    rowindex = {key:i for i,key in enumerate(outs)}
    u = [[a.R() for _ in ins] for _ in outs]
    for c,(m,n,b) in enumerate(ins):
        for aa in (-1,1):
            for j in (2,4):
                mp,np = m+aa,n+b
                if abs(mp)>j or abs(np)>j:
                    continue
                coef = (a.sqrtq(Q(4,j+1))*a.cg(1,aa,3,m,j,mp)
                        *a.cg(1,b,3,n,j,np))
                u[rowindex[(j,mp,np,aa)]][c] = coef
    ut = a.transpose(u)
    gram = a.multiply(ut,u)
    assert gram == a.identity(32)
    assert len(ins) == 32 and len(outs) == 68
    gamma = Q(1,3)
    assert allowed(3,gamma)
    assert not allowed(2,gamma) and not allowed(4,gamma)
    compressed = [[x if allowed(out[0],gamma) else a.R() for x in row]
                  for out,row in zip(outs,u)]
    assert a.zero(compressed)
    assert a.zero(a.multiply(a.transpose(compressed),compressed))
    # A fixed normal input supplies finite norm/probability data, not a delta
    # function on SU(2). Both output components and all spin labels are kept.
    col = ins.index((3,3,1))
    witness = [row[col] for row in u]
    witness_norm = sum((x*x for x in witness),a.R())
    assert witness_norm == a.R(1)
    spinor_weights = {aa:sum((w*w for out,w in zip(outs,witness) if out[3]==aa),a.R())
                       for aa in (-1,1)}
    assert spinor_weights[1] == a.R(Q(4,5))
    assert spinor_weights[-1] == a.R(Q(1,5))
    # Retaining the intermediate non-EPRL sector gives an exact inverse on
    # every original input and any reference; inserting Pi between steps
    # replaces this whole return by zero.
    assert gram == a.identity(32)
    assert not any(allowed(j,gamma) for j in (2,4))
    # Coefficient matrix is saved sparsely in exact radicals for independent
    # inspection. No grid quadrature or fitted matrix is used.
    sparse = [[r,c,x.text()] for r,row in enumerate(u) for c,x in enumerate(row) if x]
    sparse_radicals = [[r,c,{str(d):str(v) for d,v in sorted(x.a.items())}]
                       for r,row in enumerate(u) for c,x in enumerate(row) if x]
    return {
        'group_1_local_Gauss': {
            'fusion_multiplicities_twice_spin': {name:{str(j):n for j,n in sorted(v.items())}
                                                 for name,v in fusion.items()},
            'invariant_dimensions': dimensions,
            'fixed_four_legs_plus_one_spinor_center': -1,
            'full_fixed_leg_tensor_space_not_only_original_invariants': True,
            'two_Dirac_spinor_copies_remain_zero': 2*dimensions['plus_one_spinor'],
            'reference_singlet_exact': True,
            'old_geometry_times_reference_singlet_dimension': 4,
            'full_reference_invariant_dimension': 13,
            'reference_is_not_automatic_unknown_matter_encoding': True,
        },
        'group_2_full_Peter_Weyl_transport': {
            'input_basis_m2_n2_spinor2': ins,
            'output_basis_j2_m2_n2_spinor2': outs,
            'dimensions': [68,32],
            'sparse_exact_unitary_multiplication_restriction': sparse,
            'sparse_exact_radical_coefficients': sparse_radicals,
            'Udagger_U_exact_identity32': True,
            'Pi_gamma_U_Pi_gamma_exact_zero': True,
            'normal_highest_matrix_element_input_norm_squared': '1',
            'fixed_input_output_spinor_component_norms_squared': {str(k):v.text() for k,v in spinor_weights.items()},
            'full_intermediate_return_exact_identity32': True,
            'stepwise_compressed_return_exact_zero': True,
            'unknown_reference_scope': 'All 32 input columns; exact identities extend by tensoring any reference. No geometry-to-matter process implementation is inferred.',
        },
        'group_3_scope_and_deletion_controls': {
            'gamma': str(gamma),
            'fixed_k_admissible': True,
            'spin_1_admissible': False,
            'spin_2_admissible': False,
            'general_strict_gamma_proof': 'For |gamma|<1, the change of 2j_plus or 2j_minus under k -> k+/-1/2 is +/-(1+/-gamma)/2, nonzero and strictly smaller than one; it cannot be the integer difference of two admissible labels.',
            'boundary_gamma_plus1_neighbor_admissible': allowed(3,Q(1)) and allowed(2,Q(1)) and allowed(4,Q(1)),
            'boundary_gamma_minus1_neighbor_admissible': allowed(3,Q(-1)) and allowed(2,Q(-1)) and allowed(4,Q(-1)),
            'half_step_legs_restore_some_Gauss_states_but_leave_EPRL_domain': True,
            'not_all_fermions_or_EPRL_matter_couplings_excluded': True,
            'no_new_group_or_new_instrument_constructed': True,
        },
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--probe-only',action='store_true')
    ap.add_argument('--write',action='store_true')
    args=ap.parse_args()
    assert not (args.probe_only and args.write)
    checks=core_math()
    if args.probe_only:
        print('INDEPENDENT_EXACT_PROBE_PASS; groups=3; blockers=0')
        print(checks['group_1_local_Gauss']['invariant_dimensions'])
        return
    core=[STAGE/'research_note_1182.md']+[HERE/n for n in
           ('admission.md','proof.md','check.py','results.json','NEXT.md')]
    ar=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    sources=ar['source_sha256']
    assert sources and ar['groups']==3 and ar['status']=='PASS'
    for relative,expected in sources.items():
        assert sha(STAGE/relative)==expected
    # Compare every exact author matrix entry after the independently chosen
    # ascending, twice-spin basis is reordered to the author's basis. We do
    # not import or execute author functions to perform this comparison.
    mine=checks['group_2_full_Peter_Weyl_transport']
    theirs=ar['g2_fundamental_pure_label_compression']
    ii={tuple(v):i for i,v in enumerate(mine['input_basis_m2_n2_spinor2'])}
    oi={tuple(v):i for i,v in enumerate(mine['output_basis_j2_m2_n2_spinor2'])}
    radicals={(r,c):{int(d):Q(v) for d,v in value.items()}
              for r,c,value in mine['sparse_exact_radical_coefficients']}
    columns=[ii[tuple(int(2*Q(v)) for v in label)] for label in theirs['input_order_m_n_b']]
    rows=[oi[tuple(int(2*Q(v)) for v in label)] for label in theirs['output_order_J_M_N_a']]
    assert len(columns)==32 and len(rows)==68
    for r,row in enumerate(theirs['full_transport_U_exact']):
        assert len(row)==32
        for c,value in enumerate(row):
            assert {int(d):Q(v) for d,v in value.items()}==radicals.get((rows[r],columns[c]),{})
    assert len(theirs['full_transport_U_exact'])==68
    matches={'all_2176_exact_Peter_Weyl_entries_after_basis_reordering':True}
    result={
        'status':'PASS_WITH_SCOPE','groups':3,'blockers':[],
        'independence':'Prior scope advice disclosed; exact fusion and full 68x32 CG construction implemented before reading author source. Author source is not imported.',
        'reviewer_code_sha256':sha(Path(__file__)),
        'reviewer_algebra_sha256':sha(HERE/'reviewer_algebra.py'),
        'author_core_sha256':{p.relative_to(STAGE).as_posix():sha(p) for p in core},
        'author_sources_sha256_verified':sources,
        'author_result_crosscheck':matches,
        'checks':checks,
        'scope':'Only the declared fixed pure Euclidean EPRL code, common local SU(2) Gauss condition and stepwise projected fundamental transport. No failure of all fermion couplings, EPRL, or GR is asserted.',
    }
    result=json.loads(json.dumps(result))
    out=HERE/'review_checks.json'
    if args.write:
        with out.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('WRITE_ONCE_PASS; groups=3; blockers=0')
    else:
        assert result==json.loads(out.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH; groups=3; blockers=0')


if __name__=='__main__':
    main()
