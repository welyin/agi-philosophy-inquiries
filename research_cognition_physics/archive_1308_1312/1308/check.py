"""1308 exact algebra and conditional joint selection; no new QFT simulation.

Default/--check is read-only. --write exclusively creates the first result.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import permutations
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'results.json'


class W:
    """Exact Q(omega), omega^2+omega+1=0, Im(omega)>0."""
    def __init__(self, a=0, b=0):
        if isinstance(a, W): self.a, self.b = a.a, a.b
        else: self.a, self.b = F(a), F(b)
    def __add__(self, other):
        z=W(other); return W(self.a+z.a,self.b+z.b)
    __radd__=__add__
    def __neg__(self): return W(-self.a,-self.b)
    def __sub__(self, other): return self+-W(other)
    def __rsub__(self, other): return W(other)+-self
    def __mul__(self, other):
        z=W(other)
        return W(self.a*z.a-self.b*z.b,self.a*z.b+self.b*z.a-self.b*z.b)
    __rmul__=__mul__
    def conj(self): return W(self.a-self.b,-self.b)
    def __eq__(self, other):
        z=W(other); return self.a==z.a and self.b==z.b
    def __bool__(self): return bool(self.a or self.b)
    def dump(self): return [str(self.a),str(self.b)]


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),W())
             for j in range(len(b[0]))] for i in range(len(a))]
def dag(a): return [[a[j][i].conj() for j in range(len(a))] for i in range(len(a[0]))]
def eye(n): return [[W(i==j) for j in range(n)] for i in range(n)]
def eq(a,b): return all(a[i][j]==b[i][j] for i in range(len(a)) for j in range(len(a[0])))


def field_table(N):
    # (Y, B, L, fermionic, weak multiplicity)
    return {
        'Q':(F(1,2*N),F(1,N),F(0),True,N),
        'uc':(F(-N-1,2*N),F(-1,N),F(0),True,0),
        'dc':(F(N-1,2*N),F(-1,N),F(0),True,0),
        'L':(F(-1,2),F(0),F(1),True,1),
        'ec':(F(1),F(0),F(-1),True,0),
        'N':(F(0),F(0),F(-1),True,0),
        'H':(F(1,2),F(0),F(0),False,1)}


def beta(N,n): return F(-43+2*n*(N+1),6)


def permutation_sign(p):
    return -1 if sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))%2 else 1


def wedge(word):
    """Coefficient of the sorted exterior monomial; repeated generators vanish."""
    if len(set(word))!=len(word): return 0,()
    return permutation_sign(word),tuple(sorted(word))


def projected_exterior_coefficient(N):
    # Project full epsilon-color/spin/weak polynomial onto the specific M_k.
    # Each color has one selected component, independently Grassmann valued.
    k=(N-1)//2
    components={a:((0,0) if a<k+1 else (1,1)) for a in range(N)}
    eps={(0,1):1,(1,0):-1}
    total=0; terms=0
    for p in permutations(range(N)):
        coefficient=permutation_sign(p)
        for j in range(k):
            a,b=components[p[2*j]],components[p[2*j+1]]
            coefficient*=eps.get((a[0],b[0]),0)*eps.get((a[1],b[1]),0)
        a=components[p[-1]]
        coefficient*=eps.get((a[0],1),0)*eps.get((a[1],1),0)
        if not coefficient: continue
        sign,word=wedge((*p,N))
        assert word==tuple(range(N+1))
        total+=coefficient*sign; terms+=1
    expected=2**k*factorial(k)*factorial(k+1)
    assert total==expected and total>0
    return dict(Nc=N,dimension=str(F(3*(N+1),2)),coefficient=total,
                nonzero_exterior_terms=terms,analytic_coefficient=expected,
                exact_r_allowed=(N+1)%4==0)


def flavor_witness(n):
    # A positive Hermitian mass matrix = rational orthogonal conjugate of a
    # 3-mode Fourier block plus positive, distinct spectators.
    omega=W(0,1)
    powers=[W(1),omega,omega*omega]
    V=[[powers[(i*j)%3] for j in range(3)] for i in range(3)]
    assert eq(mm(V,dag(V)),[[W(3*(i==j)) for j in range(3)] for i in range(3)])
    d=[[W((i+1) if i==j else 0) for j in range(3)] for i in range(3)]
    core=[[z*F(1,3) for z in row] for row in mm(mm(V,d),dag(V))]
    M=[[W() for _ in range(n)] for _ in range(n)]
    for i in range(3):
        for j in range(3): M[i][j]=core[i][j]
    for i in range(3,n): M[i][i]=W(i+1)
    c,s=F(99,101),F(20,101)
    R=eye(n)
    for i in range(2,n-1):
        r=eye(n);r[i][i]=r[i+1][i+1]=W(c);r[i][i+1]=W(s);r[i+1][i]=W(-s)
        R=mm(r,R)
    assert eq(mm(R,dag(R)),eye(n))
    M=mm(mm(R,M),dag(R))
    assert eq(M,dag(M))
    visited={0}
    while True:
        more=visited|{j for i in visited for j in range(n) if M[i][j]}
        if more==visited: break
        visited=more
    assert len(visited)==n
    triangle=M[0][1]*M[1][2]*M[2][0]
    assert triangle.b!=0
    return dict(generations=n,positive_simple_eigenvalues=list(range(1,n+1)),
                exact_unitary_conjugation=True,connected_flavor_graph=True,
                triangle_Qomega=triangle.dump(),imaginary_part_omega_coefficient=str(triangle.b),
                beta_w=str(beta(3,n)),phase_probability_certificate_recomputed=False)


def compute():
    charges=[]; beta_checks=[]
    for N in range(3,32):
        fields=field_table(N)
        xs={name:(N+2)*b-l-4*y for name,(y,b,l,_,_) in fields.items()}
        assert list(xs.values())==[1,1,-3,1,-3,1,-2]
        assert N*fields['Q'][0]+fields['L'][0]==0
        # Conjugate fields and r act with the opposite X charge.
        assert all((-x)%4 == (3 if fields[name][3] else 2) for name,x in xs.items())
        charges.append(dict(Nc=N,X=[str(x) for x in xs.values()],unit_exchange_r_residue=(N+1)%4))
        for n in range(1,9):
            # Actual diagonal fundamental generator eigenvalues, not closed form.
            trace_weyl=sum((v*v for _ in range(n) for name,info in fields.items()
                           if info[3] for _ in range(info[4]) for v in (F(1,2),F(-1,2))),F(0))
            trace_scalar=F(1,2)
            actual=F(-11,3)*2+F(2,3)*trace_weyl+F(1,3)*trace_scalar
            assert actual==beta(N,n)
            beta_checks.append((N,n,str(actual)))
    # Analytic inequalities bound N<=6 and n<=5; enumeration is a diagnostic.
    survivors=[(N,n) for N in range(3,64) for n in range(3,64)
               if (N+1)%4==0 and beta(N,n)<=0]
    assert survivors==[(3,3),(3,4),(3,5)]
    deletion={
        'without_E':dict(Nc=5,ng=3,unit_exchange=False,r_classical=True,r_quantum=False,phase=True,weak=True),
        'without_R':dict(Nc=5,ng=3,unit_exchange=True,r=False,phase=True,weak=True),
        'without_P':dict(Nc=7,ng=2,unit_exchange=True,r=True,phase=False,weak=True),
        'without_W':dict(Nc=7,ng=3,unit_exchange=True,r=True,phase=True,weak=False)}
    for row in deletion.values():
        row['beta_w']=str(beta(row['Nc'],row['ng']))
        assert (beta(row['Nc'],row['ng'])<=0)==row['weak']
        assert (row['ng']>=3)==row['phase']
    assert beta(3,3)==F(-19,6) and beta(3,4)==F(-11,6) and beta(3,5)==F(-1,2)
    bindings=[HERE/'proof.md',Path(__file__).resolve(),HERE.parent/'_admission/round1308_admission.md',
              BASE/'archive_531_553/research_note_531.md',BASE/'archive_531_553/research_note_553.md',
              BASE/'archive_585_628/research_note_604.md',BASE/'archive_1290_1307/1299/proof.md',
              BASE/'archive_1290_1307/1303/proof.md']
    return dict(schema='round1308_joint_selection_v1',status='PASS_CONDITIONAL_LEADING_ORDER',
                new_axioms=0,unit_operator_dimension_cap_assumed=False,
                analytic_scope='Nc>=3 and ng>=1; all-polynomial r rule; specified P; dimension-four one-loop W',
                field_X_checks=charges,exact_weak_trace_checks=len(beta_checks),
                exterior_witnesses=[projected_exterior_coefficient(N) for N in (3,5,7)],
                selected_Nc_ng=[list(v) for v in survivors],deletion_witnesses=deletion,
                structure_flavor_witnesses=[flavor_witness(n) for n in (3,4,5)],
                same_interacting_instrument_claim=False,all_global_backgrounds_for_general_Nc=False,
                unchanged_source_sha256={p.relative_to(BASE.parent).as_posix():sha256(p.read_bytes()).hexdigest() for p in bindings})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=compute()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        assert json.loads(OUT.read_text(encoding='utf-8'))==result,'Saved result or source binding differs'
    print(json.dumps(dict(status=result['status'],selected_Nc_ng=result['selected_Nc_ng'],
                         exact_weak_trace_checks=result['exact_weak_trace_checks'],
                         exterior_coefficients=[v['coefficient'] for v in result['exterior_witnesses']],
                         source_bindings=len(result['unchanged_source_sha256'])),ensure_ascii=False))


if __name__=='__main__': main()
