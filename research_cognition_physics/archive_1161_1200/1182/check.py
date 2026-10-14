"""1182: fixed pure EPRL support and fundamental-spinor transport.

Exact CG arithmetic reused from the frozen author helper; no reviewer import.
Default compares the saved JSON without writing; --write creates it once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product, permutations
import json
import hashlib
from math import factorial, gcd, sqrt
from pathlib import Path
import os
for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[name] = "1"


class R:
    """Exact real sum q[d] sqrt(d), d positive and squarefree."""
    __slots__ = ("q",)
    def __init__(self, q=0):
        self.q = {d: Q(v) for d, v in q.items() if v} if isinstance(q, dict) else ({1: Q(q)} if q else {})
    def __add__(self, other):
        other = other if isinstance(other, R) else R(other)
        out = dict(self.q)
        for d, v in other.q.items():
            out[d] = out.get(d, Q(0)) + v
        return R(out)
    __radd__ = __add__
    def __neg__(self):
        return R({d: -v for d,v in self.q.items()})
    def __sub__(self, other):
        return self + (-other if isinstance(other,R) else -R(other))
    def __rsub__(self, other):
        return R(other) - self
    def __mul__(self, other):
        other = other if isinstance(other, R) else R(other)
        out = {}
        for d,a in self.q.items():
            for e,b in other.q.items():
                g = gcd(d,e)
                f = d*e//(g*g)
                out[f] = out.get(f,Q(0)) + a*b*g
        return R(out)
    __rmul__ = __mul__
    def __eq__(self, other):
        other = other if isinstance(other,R) else R(other)
        return self.q == other.q
    def __bool__(self):
        return bool(self.q)
    def value(self):
        return sum(float(v)*sqrt(d) for d,v in self.q.items())
    def data(self):
        return {str(d): f"{v.numerator}/{v.denominator}" for d,v in sorted(self.q.items())}
    def __repr__(self):
        return " + ".join(f"({v})" + (f"*sqrt({d})" if d!=1 else "") for d,v in sorted(self.q.items())) or "0"


def root(q):
    q = Q(q)
    assert q >= 0
    if not q:
        return R()
    n, outside, inside, p = q.numerator*q.denominator, 1, 1, 2
    while p*p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        outside *= p**(e//2)
        if e%2:
            inside *= p
        p += 1
    inside *= n
    ans = R({inside: Q(outside,q.denominator)})
    assert ans*ans == R(q)
    return ans


def integer(q):
    q = Q(q)
    assert q.denominator == 1
    return q.numerator


@lru_cache(None)
def cg(j1,m1,j2,m2,j,m):
    if m1+m2 != m or abs(m)>j or j>j1+j2 or j<abs(j1-j2):
        return R()
    facts = [j+j1-j2, j-j1+j2, j1+j2-j, j1+j2+j+1,
             j+m,j-m,j1-m1,j1+m1,j2-m2,j2+m2]
    if any(x<0 or Q(x).denominator!=1 for x in facts):
        return R()
    f = [factorial(integer(x)) for x in facts]
    rad = Q((2*j+1)*f[0]*f[1]*f[2],f[3])
    for v in f[4:]:
        rad *= v
    sm = Q(0)
    for z in range(integer(j1+j2+j)+2):
        args = [z,j1+j2-j-z,j1-m1-z,j2+m2-z,
                j-j2+m1+z,j-j1-m2+z]
        if any(x<0 or Q(x).denominator!=1 for x in args):
            continue
        den = 1
        for x in args:
            den *= factorial(integer(x))
        sm += Q((-1)**z,den)
    return root(rad)*sm


def magnetic(j):
    return [j-i for i in range(integer(2*j)+1)]


@lru_cache(None)
def invariant(j,a,ms):
    if sum(ms)!=0:
        return R()
    m=ms[0]+ms[1]
    if abs(m)>a:
        return R()
    return cg(j,ms[0],j,ms[1],a,m)*cg(j,ms[2],j,ms[3],a,-m)*((-1)**integer(a-m))*root(Q(1,2*a+1))


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),R())


def mm(a,b):
    return [[dot(row,col) for col in zip(*b)] for row in a]


def transpose(a):
    return [list(x) for x in zip(*a)]


def eye(n):
    return [[R(int(i==j)) for j in range(n)] for i in range(n)]


def subtract(a,b):
    return [[x-y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]


def determinant(a):
    n=len(a)
    out=R()
    for perm in permutations(range(n)):
        inv=sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=R((-1)**inv)
        for i,j in enumerate(perm):
            term=term*a[i][j]
        out=out+term
    return out


def encode(a):
    return [[x.data() for x in row] for row in a]


def fusion(twice_spins):
    """Integer SU(2) representation multiplicities, keyed by twice total spin."""
    counts={0:1}
    for twice_j in twice_spins:
        new={}
        for twice_old,mult in counts.items():
            for twice_total in range(abs(twice_old-twice_j),twice_old+twice_j+1,2):
                new[twice_total]=new.get(twice_total,0)+mult
        counts=new
    dimension=1
    for j in twice_spins:
        dimension*=j+1
    assert sum((j+1)*v for j,v in counts.items())==dimension
    return counts


def rational(x):
    assert not x.q or x.q.keys()=={1}
    return x.q.get(1,Q(0))


def zero(rows,cols):
    return [[R() for _ in range(cols)] for _ in range(rows)]


def mv(a,v):
    return [dot(row,v) for row in a]


def kron(a,b):
    return [[x*y for x in ra for y in rb] for ra in a for rb in b]


def add(a,b):
    return [[x+y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]


def admissible(gamma,k):
    assert -1<=gamma<=1
    return all(x>=0 and (2*x).denominator==1
               for x in ((1+gamma)*k/2,(1-gamma)*k/2))


def label_data(labels):
    return [[str(x) for x in row] for row in labels]


def source_hashes():
    stage=Path(__file__).resolve().parents[1]
    sources=[
        "_shared/goal.md",
        "../archive_585_628/research_note_590.md",
        "../archive_629_652/research_note_649.md",
        "research_note_1173.md",
        "research_note_1178.md",
        "research_note_1179.md",
        "research_note_1181.md",
        "1181/proof.md",
        "1181/check.py",
        "1181/results.json",
    ]
    return {p:hashlib.sha256((stage/p).read_bytes()).hexdigest() for p in sources}


def calculate():
    geometry=[3,3,3,3]
    pure=fusion(geometry)
    with_spinor=fusion(geometry+[1])
    with_reference=fusion(geometry+[1,1])
    changed_low=fusion([2,3,3,3,1])
    changed_high=fusion([4,3,3,3,1])
    assert pure[0]==4 and pure[2]==9
    assert with_spinor.get(0,0)==0
    assert (-1)**sum(geometry)==1 and (-1)**sum(geometry+[1])==-1
    assert with_reference[0]==13==pure[0]+pure[2]
    assert changed_low[0]==6 and changed_high[0]==7
    pair_count=[(a,b) for a in range(4) for b in range(4) if abs(a-b)<=1<=a+b]
    assert len(pair_count)==9

    # The all-gamma proof is analytic, not a parameter scan:
    # 2 Delta(j_plus)=(1+gamma)/2 is strictly between integers 0 and 1.
    endpoint_two_delta=[(1+Q(g))/2 for g in (-1,1)]
    assert endpoint_two_delta==[Q(0),Q(1)]
    k,gamma=Q(3,2),Q(1,3)
    assert admissible(gamma,k)
    assert not admissible(gamma,Q(1)) and not admissible(gamma,Q(2))
    assert all(admissible(Q(g),q) for g in (-1,1) for q in (Q(1),k,Q(2)))

    # Full pointwise fundamental SU(2) matrix multiplication on a normalized
    # Peter-Weyl link wavefunction, with spinor input and output indices retained.
    inputs=list(product(magnetic(k),magnetic(k),magnetic(Q(1,2))))
    outputs=[(j,m,n,a) for j in (Q(1),Q(2))
             for m,n,a in product(magnetic(j),magnetic(j),magnetic(Q(1,2)))]
    assert len(inputs)==32 and len(outputs)==68
    u=zero(68,32)
    for row,(j,mout,nout,a) in enumerate(outputs):
        factor=root(Q(4)/(2*j+1))
        for col,(m,n,b) in enumerate(inputs):
            u[row][col]=factor*cg(k,m,Q(1,2),a,j,mout)*cg(k,n,Q(1,2),b,j,nout)
    ut=transpose(u)
    gram=mm(ut,u)
    assert gram==eye(32)
    keep=[admissible(gamma,row[0]) for row in outputs]
    assert not any(keep)
    compressed=[[x if keep[row] else R() for x in line] for row,line in enumerate(u)]
    assert compressed==zero(68,32)
    assert mm(ut,compressed)==zero(32,32)
    outside=subtract(u,compressed)
    assert mm(transpose(outside),outside)==eye(32)
    range_projector=mm(u,ut)
    # The 68-dimensional carrier includes 36 unused orthogonal output directions.
    # The reverse U^dagger recovers on im(U), not as a 68->32 isometry everywhere.
    assert sum((range_projector[i][i] for i in range(68)),R())==R(32)
    assert mm(range_projector,u)==u

    witnesses=[]
    for b,expected in ((Q(1,2),{Q(1):Q(0),Q(2):Q(1)}),
                       (Q(-1,2),{Q(1):Q(3,4),Q(2):Q(1,4)})):
        col=inputs.index((k,k,b))
        coeff=[row[col] for row in u]
        weight={j:rational(sum((x*x for x,label in zip(coeff,outputs) if label[0]==j),R()))
                for j in (Q(1),Q(2))}
        assert weight==expected
        witnesses.append({
            "input_m_n_b":[str(k),str(k),str(b)],
            "J_weights_exact":{str(j):str(p) for j,p in weight.items()},
            "nonzero_outputs_exact":[{"J_M_N_a":[str(z) for z in label],"amplitude":x.data()}
                                     for x,label in zip(coeff,outputs) if x],
        })
    top_col=inputs.index((k,k,Q(1,2)))
    assert u[outputs.index((Q(2),Q(2),Q(2),Q(1,2)))][top_col]==root(Q(4,5))
    assert u[outputs.index((Q(2),Q(1),Q(2),Q(-1,2)))][top_col]==root(Q(1,5))

    # A real, charged reference is an extra resource. A known singlet is not
    # an encoding of an arbitrary unknown matter spinor.
    singlet=[R(),root(Q(1,2)),-root(Q(1,2)),R()]
    assert dot(singlet,singlet)==R(1)
    raising=[[R(),R(1)],[R(),R()]]
    lowering=transpose(raising)
    z=[[R(Q(1,2)),R()],[R(),R(Q(-1,2))]]
    for generator in (raising,lowering,z):
        total=add(kron(generator,eye(2)),kron(eye(2),generator))
        assert mv(total,singlet)==[R()]*4
    ps=[[x*y for y in singlet] for x in singlet]
    assert mm(ps,ps)==ps
    fixed_down_embedding=[[R(),R()],[R(1),R()],[R(),R()],[R(),R(1)]]
    projected=mm(ps,fixed_down_embedding)
    projected_gram=mm(transpose(projected),projected)
    assert projected_gram==[[R(Q(1,2)),R()],[R(),R()]]
    assert 4<4*2
    return {
        "round":1182,"status":"PASS","groups":3,
        "arithmetic":"integer SU2 fusion and exact rational square-root coefficients; no floating certificate",
        "source_sha256":source_hashes(),
        "primary_sources":{
            "common_fermion_Gauss_intertwiner":"https://arxiv.org/html/2112.08865",
            "declared_Dirac_spinfoam_coupling":"https://arxiv.org/html/1101.3264",
            "Euclidean_EPRL_map":"https://arxiv.org/html/0912.0540v2",
        },
        "g1_fixed_node_support":{
            "twice_geometry_spins":geometry,
            "pure_geometry_total_dimension":256,
            "with_one_spinor_total_dimension":512,
            "pure_geometry_spin_multiplicities":{str(Q(j,2)):v for j,v in sorted(pure.items())},
            "with_one_spinor_multiplicities":{str(Q(j,2)):v for j,v in sorted(with_spinor.items())},
            "pure_geometry_invariant_dimension":4,
            "same_SU2_one_spinor_invariant_dimension":0,
            "center_sign_geometry":1,
            "center_sign_geometry_plus_spinor":-1,
            "arbitrary_spectator_reference_does_not_change_zero_support":True,
            "requires_closed_internal_node_total_Gauss":True,
        },
        "g2_fundamental_pure_label_compression":{
            "general_gamma_range":"-1 < gamma < 1",
            "analytic_two_delta_j_plus":"(1+gamma)/2 lies strictly in (0,1), containing no integer",
            "general_theorem_is_not_a_scan":True,
            "endpoints_plus_minus_one_excluded":True,
            "fixed_gamma":"1/3",
            "fixed_admissible_labels":"k=3*n/2, n nonnegative integer",
            "link_input_spin":"3/2",
            "input_dimension":32,
            "output_carrier_dimension":68,
            "output_spins":["1","2"],
            "input_order_m_n_b":label_data(inputs),
            "output_order_J_M_N_a":label_data(outputs),
            "radical_encoding":"entry dictionary d -> rational coefficient of sqrt(d)",
            "full_transport_U_exact":encode(u),
            "Udagger_U_identity32_exact":True,
            "all_output_labels_outside_pure_set_exact":True,
            "projected_transport_zero_exact":True,
            "outside_transport_Gram_identity32_exact":True,
            "retained_probability_all_unknown_and_R_exact":"0",
            "outside_probability_all_unknown_and_R_exact":"1",
            "pure_label_effect_gap_against_any_normalized_inside_output_exact":"1",
            "distance_to_unnormalized_zero_exact":"1",
            "distance_to_any_normalized_inside_pure_vector_exact":"sqrt(2)",
            "fixed_nonzero_witnesses":witnesses,
            "single_matrix_entry_is_not_the_full_isometry":True,
        },
        "g3_keep_intermediate_and_support_escapes":{
            "round_trip_keep_all_intermediate_equals_identity32_exact":True,
            "round_trip_with_middle_pure_projection_equals_zero_exact":True,
            "reverse_recovery_only_on_32dim_image_within_68dim_carrier":True,
            "spin1_multiplicity_of_original_four_legs":9,
            "spin1_pair_count_a_b":[list(x) for x in pair_count],
            "with_two_spinors_invariant_dimension":13,
            "original_geometry_times_known_singlet_subcode_dimension":4,
            "original_geometry_times_unknown_spinor_input_dimension":8,
            "known_singlet_subcode_does_not_embed_full_unknown_code":True,
            "matter_reference_singlet_exact":[x.data() for x in singlet],
            "total_raising_lowering_z_annihilate_singlet_exact":True,
            "fixed_down_reference_then_singlet_projection_Gram_exact":encode(projected_gram),
            "fixed_reference_projection_rank":1,
            "changed_one_leg_to_spin1_invariant_dimension":6,
            "changed_one_leg_to_spin2_invariant_dimension":7,
            "changed_legs_not_EPRL_admissible_at_gamma_one_third":True,
        },
        "scope":{
            "same_local_SU2_matter_identity_is_added_condition":True,
            "not_pseudospin_or_neutral_composite_clock_no_go":True,
            "kinematic_link_code_not_original_node_physical_code":True,
            "fundamental_SU2_multiplier_is_selected_primitive_not_all_Dirac_dynamics":True,
            "does_not_exclude_open_charged_boundary_or_noninvariant_subregion":True,
            "does_not_exclude_other_kernels_or_retained_virtual_intermediate_states":True,
            "no_free_unknown_spinor_encoding_via_added_reference":True,
            "no_hardware_Hamiltonian_energy_clock_or_GR_claim":True,
            "no_dynamical_coefficient_selection_claim":True,
            "no_dark_sector_claim":True,
        },
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    ap.add_argument("--calculate-only",action="store_true",help="Recompute exact algebra without writing or requiring results.")
    args=ap.parse_args()
    assert not (args.write and args.calculate_only)
    result=calculate()
    target=Path(__file__).with_name("results.json")
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
        label="WRITE_ONCE_PASS"
    elif args.calculate_only:
        label="CALCULATE_ONLY_PASS"
    else:
        assert result==json.loads(target.read_text(encoding="utf-8"))
        label="READ_ONLY_MATCH"
    print(label,result["groups"],"groups")
    print("fixed-node singlet dimension:",result["g1_fixed_node_support"]["same_SU2_one_spinor_invariant_dimension"])
    print("transport: 68x32, Udagger U=I, pure-label compression=0")
    print("reference/changed-leg invariant dimensions: 13 / 6 / 7")


if __name__=="__main__":
    main()

