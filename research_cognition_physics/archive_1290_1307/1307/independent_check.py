"""Independent 1307 checks: weighted Hilbert momentum and neutral source kernel.
No author imports. Exact Gaussian rationals, not numerical PDE certification.
Default and --check read only; --write creates independent_results.json once.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent

def z(v=0, im=0):
    if isinstance(v, tuple):
        return v
    if isinstance(v, complex):
        return (Q(v.real), Q(v.imag))
    return (Q(v), Q(im))

ZERO = z()
ONE = z(1)

def add(a,b): return (a[0]+b[0],a[1]+b[1])
def neg(a): return (-a[0],-a[1])
def sub(a,b): return add(a,neg(b))
def times(a,b): return (a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def conj(a): return (a[0],-a[1])
def scale(a,c): return (a[0]*c,a[1]*c)
def mat(rows): return [[z(v) for v in row] for row in rows]
def zero(n,m=None): return [[ZERO for _ in range(n if m is None else m)] for _ in range(n)]
def eye(n): return [[ONE if i==j else ZERO for j in range(n)] for i in range(n)]
def adj(a): return [[conj(v) for v in row] for row in zip(*a)]
def trans(a): return [list(row) for row in zip(*a)]
def conjugate(a): return [[conj(v) for v in row] for row in a]
def plus(a,b): return [[add(x,y) for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]
def minus(a,b): return plus(a,times_scalar(b,-1))
def times_scalar(a,c): return [[scale(v,Q(c)) for v in row] for row in a]
def mm(a,b):
    bt=trans(b)
    return [[sum_z(times(x,y) for x,y in zip(row,col)) for col in bt] for row in a]
def sum_z(xs):
    out=ZERO
    for v in xs: out=add(out,v)
    return out
def kron(a,b):
    return [[times(a[i][j],b[k][l]) for j in range(len(a[0])) for l in range(len(b[0]))]
            for i in range(len(a)) for k in range(len(b))]
def blocks(a,b,c,d):
    return [ra+rb for ra,rb in zip(a,b)]+[rc+rd for rc,rd in zip(c,d)]
def comm(a,b): return minus(mm(a,b),mm(b,a))
def iszero(a): return all(v==ZERO for row in a for v in row)
def max_abs_component(a): return max(abs(c) for row in a for v in row for c in v)
def real_rank(a):
    a=[list(map(Q,row)) for row in a]
    if not a: return 0
    r=0
    for j in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][j]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        pivot=a[r][j]
        a[r]=[v/pivot for v in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:
                t=a[i][j]; a[i]=[u-t*v for u,v in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r
def c_rank(a):
    # Rank over C equals half the rank of its realification.
    r=[[v[0] for v in row]+[-v[1] for v in row] for row in a]
    r += [[v[1] for v in row]+[v[0] for v in row] for row in a]
    rank=real_rank(r)
    assert rank%2==0
    return rank//2
def hermitian_basis(n):
    out=[]
    for j in range(n):
        b=zero(n); b[j][j]=ONE;out.append(b)
    for j in range(n):
        for k in range(j+1,n):
            b=zero(n);b[j][k]=b[k][j]=ONE;out.append(b)
            b=zero(n);b[j][k]=z(0,1);b[k][j]=z(0,-1);out.append(b)
    return out
def source_nullity(y):
    cols=[]
    for b in hermitian_basis(len(y[0])):
        cols.append([c for row in mm(y,b) for v in row for c in v])
    return len(cols)-real_rank(list(zip(*cols)))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def checks():
    sx=mat([[0,1],[1,0]])
    sy=mat([[0,-1j],[1j,0]])
    sz=mat([[1,0],[0,-1]])
    paulis=[sx,sy,sz]
    p=[Q(1),Q(2),Q(3)]
    spin=zero(2)
    for pj,sigma in zip(p,paulis):spin=plus(spin,times_scalar(sigma,pj))
    y=mat([[1,1j,0],[0,2,1+1j],[1,0,3]])
    a=mat([[1,1+2j,-1j],[1-2j,-2,3],[1j,3,2]])
    assert a==adj(a) and c_rank(y)==3
    # All-left Weyl A becomes A^T in the right-handed conjugate dictionary.
    # M = v Y*; choose v=1, sufficient to check the exact intertwining.
    ar=trans(a);mass=conjugate(y)
    assert mm(mass,ar)==conjugate(mm(y,a))
    wrong_dictionary=minus(mm(mass,a),conjugate(mm(y,a)))
    assert not iszero(wrong_dictionary)
    I3,I2=eye(3),eye(2);Z6=zero(6)
    h=blocks(times_scalar(kron(I3,spin),-1),kron(mass,I2),
             kron(adj(mass),I2),kron(I3,spin))
    p_a=blocks(Z6,Z6,Z6,kron(ar,I2))
    assert h==adj(h) and p_a==adj(p_a)
    identities=[]
    nonzero_obstruction=[]
    for i,sigma in enumerate(paulis):
        alpha=blocks(times_scalar(kron(I3,sigma),-1),Z6,Z6,kron(I3,sigma))
        ap=mm(alpha,p_a)
        # T^{0i} from the symmetric Hilbert kinetic tensor.
        k_i=plus(times_scalar(p_a,p[i]/2),
                 times_scalar(plus(mm(h,ap),mm(ap,h)),Q(1,4)))
        actual=comm(h,k_i)
        predicted=plus(times_scalar(comm(h,p_a),p[i]/2),
                       times_scalar(comm(mm(h,h),ap),Q(1,4)))
        assert actual==predicted and k_i==adj(k_i)
        upper=[row[6:] for row in actual[:6]]
        expected=times_scalar(kron(mm(mass,ar),I2),p[i]/2)
        assert upper==expected and not iszero(upper)
        identities.append(True)
        nonzero_obstruction.append(str(max_abs_component(upper)))
    # Deliberate incomplete check: on a frozen Higgs, energy alone can be
    # supplemented by half a mass term. This is not a conserved full tensor.
    I=eye(2);Z=zero(2)
    one_h=blocks(times_scalar(spin,-1),I,I,spin)
    energy_loophole=blocks(Z,times_scalar(I,Q(1,2)),times_scalar(I,Q(1,2)),spin)
    assert iszero(comm(one_h,energy_loophole))
    one_pa=blocks(Z,Z,Z,I)
    alpha=blocks(times_scalar(sx,-1),Z,Z,sx)
    one_k=plus(times_scalar(one_pa,p[0]/2),
               times_scalar(plus(mm(one_h,mm(alpha,one_pa)),mm(mm(alpha,one_pa),one_h)),Q(1,4)))
    assert not iszero(comm(one_h,one_k))
    ys=[
        zero(3),
        mat([[1,1j,1],[2,2j,2],[-1,-1j,-1]]),
        mat([[1,1j,1],[0,1,1j],[1,1+1j,1+1j]]),
        y,
    ]
    kernel_table=[]
    for yy in ys:
        rank=c_rank(yy);nullity=source_nullity(yy)
        assert nullity==(3-rank)**2
        kernel_table.append({'complex_Y_rank':rank,'Hermitian_extra_source_dimension':nullity})
    assert [x['complex_Y_rank'] for x in kernel_table]==[0,1,2,3]
    y2=ys[2]
    v=mat([[-2],[-1j],[1]])
    projector=times_scalar(mm(v,adj(v)),Q(1,6))
    assert projector==adj(projector) and mm(projector,projector)==projector
    assert iszero(mm(y2,projector)) and not iszero(projector)
    # A generic perturbation outside the kernel is excluded.
    bad=plus(projector,times_scalar(eye(3),Q(1,7)))
    assert not iszero(mm(y2,bad))
    return {
      'status':'PASS_INDEPENDENT_LEADING_SOURCE_CLASS',
      'scope':'same full local parent, stated linear weighted kinetic source class; no quantum all-source uniqueness',
      'exact_complex_dictionary_M_AR_equals_conjugate_YA':True,
      'wrong_untransposed_right_weight_detected':True,
      'Hilbert_momentum_commutator_identity_each_spatial_axis':identities,
      'upper_right_obstruction_max_real_imag_coefficients':nonzero_obstruction,
      'frozen_Higgs_energy_only_countercheck_commutes':True,
      'same_countercheck_fails_full_spatial_momentum':True,
      'exact_kernel_dimension_table':kernel_table,
      'non_diagonal_complex_rank_two_kernel_projector_exact':True,
      'general_dimension_formula':'(n-rank(Y))^2; analytic proof, samples only check',
      'full_parent_linearized_Grassmann_degree_two_required':True,
      'all_local_shell_data_not_a_single_background_solution':True,
      'background_Ward_exchange_terms_not_deleted':True,
      'new_axioms':0,
    }

def payload():
    out=checks()
    sources=['proof.md','check.py','results.json','NEXT.md','../research_note_1307.md',
             '../_shared/round1307_admission.md',
             '../../archive_1009_1043/research_note_1027.md',
             '../../archive_1009_1043/research_note_1035.md',
             '../../archive_1046_1062/research_note_1060.md',
             '../1305/proof.md']
    out['source_sha256']={p:sha(HERE/p) for p in sources}
    out['independent_code_sha256']=sha(Path(__file__).resolve())
    return out

def main():
    parser=argparse.ArgumentParser()
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--write',action='store_true')
    modes.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=payload();target=HERE/'independent_results.json'
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
        print('CREATED',result['status'])
    else:
        saved=json.loads(target.read_text(encoding='utf-8'))
        assert result==saved,'Saved independent result or bound files changed'
        print('READ_ONLY_MATCH',result['status'])
    print('source_kernel_dimensions',[r['Hermitian_extra_source_dimension'] for r in result['exact_kernel_dimension_table']])

if __name__=='__main__':main()

