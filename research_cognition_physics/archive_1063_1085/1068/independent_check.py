"""Independent exact finite-algebra check of archive 1068.
Python standard library only. No imports from the author's implementation.
Pure states are represented by rational unnormalized B and their normalized
density matrices, so no floating-point square roots are used.

SCOPE: Exact finite examples and constructive maximum-rank witnesses for
2 <= n <= 12. These tests do not by themselves prove the universal CPTP
theorem or continuity step. For the charge rank upper bounds, each support
has the stated block shape: PP <= r, QQ <= q, cross <= 2*min(r,q).
The constructed matrices attain all three upper bounds.
"""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

class C:
    __slots__ = ("r", "i")
    def __init__(self, r=0, i=0):
        if isinstance(r, C):
            self.r, self.i = r.r, r.i
        else:
            self.r, self.i = F(r), F(i)
    def __add__(self, other):
        other = C(other)
        return C(self.r + other.r, self.i + other.i)
    __radd__ = __add__
    def __neg__(self):
        return C(-self.r, -self.i)
    def __sub__(self, other):
        return self + (-C(other))
    def __rsub__(self, other):
        return C(other) - self
    def __mul__(self, other):
        other = C(other)
        return C(self.r*other.r-self.i*other.i,
                 self.r*other.i+self.i*other.r)
    __rmul__ = __mul__
    def __truediv__(self, other):
        other = C(other)
        d = other.abs2()
        if d == 0:
            raise ZeroDivisionError
        p = self * other.conj()
        return C(p.r/d, p.i/d)
    def __rtruediv__(self, other):
        return C(other)/self
    def __eq__(self, other):
        if not isinstance(other, (C, int, F)):
            return False
        other = C(other)
        return self.r == other.r and self.i == other.i
    def __bool__(self):
        return bool(self.r or self.i)
    def conj(self):
        return C(self.r, -self.i)
    def abs2(self):
        return self.r*self.r+self.i*self.i
    def __repr__(self):
        if not self.i:
            return str(self.r)
        return f"({self.r}{'+' if self.i >= 0 else ''}{self.i}i)"

def zero(m, n=None):
    return [[C() for _ in range(m if n is None else n)] for _ in range(m)]

def eye(n):
    a = zero(n)
    for i in range(n):
        a[i][i] = C(1)
    return a

def diag(values):
    a = zero(len(values))
    for i, value in enumerate(values):
        a[i][i] = C(value)
    return a

def tr(a):
    return sum((a[i][i] for i in range(len(a))), C())

def transpose(a):
    return [list(row) for row in zip(*a)]

def conj(a):
    return [[x.conj() for x in row] for row in a]

def adj(a):
    return transpose(conj(a))

def scale(a, s):
    return [[x*s for x in row] for row in a]

def add(*matrices):
    a = zero(len(matrices[0]), len(matrices[0][0]))
    for b in matrices:
        a = [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]
    return a

def mul(a,b):
    assert len(a[0]) == len(b)
    c = zero(len(a), len(b[0]))
    for i, row in enumerate(a):
        for k, x in enumerate(row):
            if x:
                for j, y in enumerate(b[k]):
                    if y:
                        c[i][j] = c[i][j] + x*y
    return c

def kron(a,b):
    return [[x*y for x in ar for y in br] for ar in a for br in b]

def rank(a):
    a = [[C(x) for x in row] for row in a]
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row,len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        p = a[row][col]
        a[row] = [v/p for v in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                f = a[i][col]
                a[i] = [v-f*w for v,w in zip(a[i],a[row])]
        row += 1
        if row == len(a):
            break
    return row

def pure(v):
    norm = sum((x.abs2() for x in v), F())
    assert norm > 0
    return [[x*y.conj()/norm for y in v] for x in v]

def rho_b(b):
    return pure([x for row in b for x in row])

def conjugate(u, rho):
    return mul(mul(u,rho),adj(u))

def b_action(u,b,v=None):
    return mul(mul(u,b),transpose(u if v is None else v))

def partial_b(rho,n):
    return [[sum((rho[i*n+k][j*n+k] for k in range(n)),C())
             for j in range(n)] for i in range(n)]

def purity(rho):
    return tr(mul(rho,rho))

def exact_rotation(n, a, b):
    u = eye(n)
    u[a][a] = u[b][b] = C(F(3,5))
    u[a][b],u[b][a] = C(F(4,5)),C(F(-4,5))
    return u

def dephase(rho):
    return diag([rho[i][i] for i in range(len(rho))])

def average_conjugations(us,rho):
    return scale(add(*(conjugate(u,rho) for u in us)),F(1,len(us)))

results = []
assertions = 0
def check(condition,message):
    global assertions
    assertions += 1
    if not condition:
        raise AssertionError(message)

def report(name, **evidence):
    results.append({"name":name,"status":"PASS",**evidence})

z = C(F(3,5),F(4,5))
imag = C(0,1)
check(z.abs2()==1,"exact phase lies on S1")
check(z*z != z and z != 1 and z*z != 1,"three charges distinct")

# For diagonal P, matrix units E_ab form an eigenbasis for B -> V B V^T.
# Every eigenspace is therefore exactly one of the PP, cross, QQ supports.
partitions = 0
rank_witnesses = 0
for n in range(2,13):
    for r in range(1,n):
        q = n-r
        v = diag([z]*r+[1]*q)
        charges = [z*z,z,C(1)]
        counts = [0,0,0]
        for i in range(n):
            for j in range(n):
                c = v[i][i]*v[j][j]
                matches = [k for k,x in enumerate(charges) if c==x]
                check(len(matches)==1,"one unique charge per matrix unit")
                counts[matches[0]] += 1
        check(counts==[r*r,2*r*q,q*q],"exact charge-space dimensions")
        pp, cross, qq = zero(n),zero(n),zero(n)
        for i in range(r):
            pp[i][i]=C(1)
        for i in range(r,n):
            qq[i][i]=C(1)
        for i in range(min(r,q)):
            cross[i][r+i]=C(1)
            cross[r+i][i]=C(1)
        for b,charge,bound in zip([pp,cross,qq],charges,[r,2*min(r,q),q]):
            check(b_action(v,b)==scale(b,charge),"maximum-rank witness has correct charge")
            check(rank(b)==bound,"maximum-rank witness attains block-shape upper bound")
            rank_witnesses += 1
        check((rank(cross)==n)==(n==2*r),"full rank cross iff balanced dimensions")
        check(rank(pp)<n and rank(qq)<n,"same-sector charges never full rank")
        partitions += 1
report("three charge supports and maximum ranks", n_range="2..12",
       nontrivial_partitions=partitions, maximum_rank_witnesses=rank_witnesses,
       exact_ranks="r, 2*min(r,n-r), n-r",
       invertible_support="cross, exactly when n=2*r")

# n=2 unequal Schmidt weights: phase resource but not common SU(2)-invariant.
b = [[C(),C(F(3,5))],[C(F(4,5)),C()]]
rho = rho_b(b)
u = diag([z,1])
check(rank(b)==2,"n2 full Schmidt rank")
check(tr(rho)==1 and purity(rho)==1,"n2 normalized pure relation")
check(b_action(u,b)==scale(b,z),"n2 common phase invariant ray")
check(conjugate(kron(u,u),rho)==rho,"n2 common phase invariant density")
check(partial_b(rho,2)==diag([F(9,25),F(16,25)]),"n2 unequal exact Schmidt weights")
su2 = exact_rotation(2,0,1)
check(mul(su2,adj(su2))==eye(2),"chosen SU2 matrix unitary")
check(su2[0][0]*su2[1][1]-su2[0][1]*su2[1][0]==1,"chosen SU2 determinant one")
check(conjugate(kron(su2,su2),rho)!=rho,"n2 resource not common SU2 invariant")
report("n2 partial entanglement need not be SU2 invariant",
       marginal_eigenvalues=["9/25","16/25"], Schmidt_rank=2)

# n=4, rank(P)=2 retains a full-rank, generally non-maximal, pure relation.
b = zero(4)
b[0][2],b[1][3],b[2][0],b[3][1] = map(C,[1,2,3,4])
u = diag([z,z,1,1])
check(rank(b)==4,"n4 rank2 sectors full Schmidt rank")
check(b_action(u,b)==scale(b,z),"n4 balanced resource phase covariance")
check(conjugate(kron(u,u),rho_b(b))==rho_b(b),"n4 balanced density invariant")
report("n4 binary sectors of rank two", rank_P=2, Schmidt_rank=4)

# n=3 embedding only constrains the participating two-dimensional support.
b = zero(3)
b[0][1],b[1][0] = C(1),C(-1)
u = diag([z,1,1])
rho = rho_b(b)
check(rank(b)==2,"embedded singlet not full hardware Schmidt rank")
check(b_action(u,b)==scale(b,z),"embedded singlet phase covariance")
check(partial_b(rho,3)==diag([F(1,2),F(1,2),0]),"unused third level")
report("n3 embedded singlet", Schmidt_rank=2, hardware_dimension=3,
       marginal_diagonal=["1/2","1/2","0"])

# Mixed antisymmetric relation with full-rank marginals in n=3.
n=3
swap = zero(n*n)
for i in range(n):
    for j in range(n):
        swap[j*n+i][i*n+j]=C(1)
rho = scale(add(eye(n*n),scale(swap,-1)),F(1,n*(n-1)))
check(tr(rho)==1,"antisymmetric state trace")
check(rank(rho)==3 and purity(rho)==F(1,3),"antisymmetric state mixed rank three")
check(partial_b(rho,n)==scale(eye(n),F(1,n)),"antisymmetric marginal maximally mixed")
rotation01 = exact_rotation(3,0,1)
rotation12 = exact_rotation(3,1,2)
phase = diag([z,1,1])
unitaries = [eye(3),phase,rotation01,rotation12,mul(phase,mul(rotation01,rotation12))]
for u in unitaries:
    check(mul(u,adj(u))==eye(n),"antisymmetric test setting unitary")
    uu = kron(u,u)
    check(mul(swap,uu)==mul(uu,swap),"swap commutes with common unitary")
    check(conjugate(uu,rho)==rho,"mixed antisymmetric state invariant")
# The commutation identity is general by S(A x A)|i,j>=(A x A)S|i,j>.
report("n3 mixed antisymmetric relation", state_rank=3,
       marginal_rank=3, purity="1/3", exact_unitary_settings=len(unitaries))

# Classical correlated mixed relation survives independent complete dephasing.
classical = zero(9)
for i in range(3):
    classical[i*3+i][i*3+i]=C(F(1,3))
local_projectors = []
for i in range(3):
    p=zero(3)
    p[i][i]=C(1)
    local_projectors.append(p)
independent = add(*(conjugate(kron(p,q),classical)
                    for p in local_projectors for q in local_projectors))
check(independent==classical,"classical mixture independently dephasing invariant")
plus3 = pure([C(1),C(1),C()])
check(rank(dephase(plus3))==2,"local dephasing not unitary")
report("n3 classical mixed relation and independent dephasing",
       state_rank=rank(classical), local_pure_output_rank=2)

# Disconnected endpoint menu works in odd dimension; an i setting fails.
b = eye(3)
rho = rho_b(b)
minus_endpoint = diag([-1,1,1])
middle_i = diag([imag,1,1])
check(b_action(minus_endpoint,b)==b,"n3 minus endpoint invariant")
check(conjugate(kron(minus_endpoint,minus_endpoint),rho)==rho,"endpoint density invariant")
middle_rho = conjugate(kron(middle_i,middle_i),rho)
check(middle_rho!=rho,"n3 i intermediate does not preserve relation")
check(tr(mul(rho,middle_rho))==F(1,9),"intermediate exact fidelity")
report("n3 discrete minus phase and failed i intermediate",
       endpoint="invariant", i_intermediate="changed", fidelity="1/9")

# Statistical preservation is weaker than pointwise preservation within Q.
j = [[C(),C(1)],[C(-1),C()]]
b = zero(4)
for offset in (0,2):
    for i in range(2):
        for k in range(2):
            b[offset+i][offset+k]=j[i][k]
u = diag([z,z.conj(),1,1])
p = diag([1,0,0,0])
q_superposition = pure([C(),C(1),C(1),C()])
q_after = conjugate(u,q_superposition)
check(rank(b)==4,"J2+J2 full Schmidt rank")
check(tr(rho_b(b))==1,"J2+J2 relation normalization")
check(b_action(u,b)==b,"J2+J2 exact common invariant vector")
check(conjugate(kron(u,u),rho_b(b))==rho_b(b),"J2+J2 density invariant")
check(mul(mul(adj(u),p),u)==p,"P statistics preserved for every input")
check(q_after!=q_superposition,"Q superposition changed")
check(q_superposition[1][2]==F(1,2),"initial Q coherence")
check(q_after[1][2]==C(F(3,10),F(-2,5)),"changed Q coherence")
report("n4 only P statistics with J2 direct sum J2", rank_P=1,
       Schmidt_rank=4, Q_coherence_before="1/2",
       Q_coherence_after=str(q_after[1][2]))

# The second conjugate representation changes the physical role.
role_dimensions=[]
for n in (2,3,4):
    b=eye(n)
    u=mul(diag([z]+[1]*(n-1)),exact_rotation(n,0,1))
    check(mul(u,adj(u))==eye(n),"role test unitary")
    check(b_action(u,b,conj(u))==b,"U tensor conjugate U invariant vector")
    check(conjugate(kron(u,conj(u)),rho_b(b))==rho_b(b),"conjugate role invariant state")
    check(conjugate(kron(u,u),rho_b(b))!=rho_b(b),"same roles need not preserve maximally entangled state")
    role_dimensions.append(n)
report("U tensor conjugate U versus same role", dimensions=role_dimensions)

# Shared ZZ noise and tensor product of two dephasing channels are different.
z_pauli = diag([1,-1])
id2 = eye(2)
bell = rho_b(eye(2))
zz = kron(z_pauli,z_pauli)
shared = average_conjugations([eye(4),zz],bell)
independent = average_conjugations(
    [kron(a,b) for a in [id2,z_pauli] for b in [id2,z_pauli]],bell)
expected = diag([F(1,2),0,0,F(1,2)])
check(shared==bell,"correlated ZZ channel preserves Bell pure state")
check(independent==expected and independent!=bell,"independent dephasing loses Bell coherence")
check(purity(shared)==1 and rank(shared)==1,"shared output pure")
check(purity(independent)==F(1,2) and rank(independent)==2,"independent output mixed")
single_plus = pure([C(1),C(1)])
single_after = average_conjugations([id2,z_pauli],single_plus)
check(single_after==scale(eye(2),F(1,2)),"single local channel is dephasing")
check(rank(single_after)==2,"single local channel not unitary")
check(partial_b(shared,2)==partial_b(independent,2),"marginals alone cannot distinguish joint channels")
report("correlated ZZ versus independent local dephasing",
       shared_rank=1, shared_purity="1", independent_rank=2,
       independent_purity="1/2", equal_Bell_output_marginals=True)

print(json.dumps({
    "checker":"independent Fraction complex finite-algebra verifier",
    "arithmetic":"exact fractions; no numpy, sympy, floating point or author imports",
    "assertions":assertions,
    "test_groups":len(results),
    "results":results,
    "scope":"Exact finite cases and constructive charge-rank bounds for n=2..12. "
            "Not a formal proof of universal CPTP classification or the continuity argument."
},ensure_ascii=False,indent=2))
print("SCRIPT_SHA256", hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

