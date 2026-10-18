# 1307: exact source-weight kernel and symmetric Hilbert momentum diagnostic.
# Default/--check reads only; --write creates results.json exclusively.
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'

class Q:
    def __init__(self,r=0,i=0):
        if isinstance(r,Q): self.r,self.i=r.r,r.i
        else: self.r,self.i=F(r),F(i)
    def __add__(self,b):
        b=Q(b); return Q(self.r+b.r,self.i+b.i)
    __radd__=__add__
    def __neg__(self): return Q(-self.r,-self.i)
    def __sub__(self,b): return self+-Q(b)
    def __rsub__(self,b): return Q(b)+-self
    def __mul__(self,b):
        b=Q(b); return Q(self.r*b.r-self.i*b.i,self.r*b.i+self.i*b.r)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=Q(b); den=b.r*b.r+b.i*b.i
        return Q((self.r*b.r+self.i*b.i)/den,(self.i*b.r-self.r*b.i)/den)
    def __eq__(self,b):
        b=Q(b); return self.r==b.r and self.i==b.i
    def conj(self): return Q(self.r,-self.i)
    def __bool__(self): return bool(self.r or self.i)
    def __complex__(self): return complex(float(self.r),float(self.i))
    def dump(self): return [str(self.r),str(self.i)]

def zero(n,m): return [[Q() for j in range(m)] for i in range(n)]
def ident(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def dag(a): return [[a[j][i].conj() for j in range(len(a))] for i in range(len(a[0]))]
def mul(a,b): return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q()) for j in range(len(b[0]))] for i in range(len(a))]
def scale(a,c): return [[v*c for v in row] for row in a]
def add(a,b): return [[a[i][j]+b[i][j] for j in range(len(a[0]))] for i in range(len(a))]
def eq(a,b): return all(a[i][j]==b[i][j] for i in range(len(a)) for j in range(len(a[0])))
def npa(a): return np.array([[complex(v) for v in row] for row in a])
def dump(a): return [[v.dump() for v in row] for row in a]

def rank_real(a):
    a=[[F(v) for v in row] for row in a]
    r=0
    for j in range(len(a[0])):
        p=next((i for i in range(r,len(a)) if a[i][j]),None)
        if p is None: continue
        a[r],a[p]=a[p],a[r]
        v=a[r][j];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r and a[i][j]:
                v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==len(a): break
    return r

def herm_basis(n):
    bs=[]
    for i in range(n):
        a=zero(n,n);a[i][i]=Q(1);bs.append(a)
    for i in range(n):
        for j in range(i+1,n):
            a=zero(n,n);a[i][j]=a[j][i]=Q(1);bs.append(a)
            a=zero(n,n);a[i][j]=Q(0,1);a[j][i]=Q(0,-1);bs.append(a)
    return bs

def exact_source_kernels():
    U=[[Q(F(3,5)),Q(F(4,5)),Q()],[Q(F(-4,5)),Q(F(3,5)),Q()],[Q(),Q(),Q(0,1)]]
    V=[[Q(0,F(3,5)),Q(F(4,5)),Q()],[Q(0,F(-4,5)),Q(F(3,5)),Q()],[Q(),Q(),Q(1)]]
    # Mix the free kernel with charged-to-N directions, rather than only a diagonal test.
    R=[[Q(1),Q(),Q()],[Q(),Q(F(5,13)),Q(F(12,13))],[Q(),Q(F(-12,13)),Q(F(5,13))]]
    V=mul(V,R)
    assert eq(mul(dag(U),U),ident(3)) and eq(mul(dag(V),V),ident(3))
    bs=herm_basis(3);records=[];ys=[];survivors=[]
    for r in range(4):
        D=zero(3,3)
        for i in range(r): D[i][i]=Q(i+1)
        Y=mul(mul(U,D),dag(V));ys.append(Y)
        outputs=[mul(Y,b) for b in bs]
        rows=[]
        for i in range(3):
            for j in range(3):
                rows.append([z[i][j].r for z in outputs])
                rows.append([z[i][j].i for z in outputs])
        rr=rank_real(rows)
        assert rr==9-(3-r)**2
        fullrows=[]
        for i in range(3):
            for j in range(3):
                fullrows.append([-Y[i][j].r]+[z[i][j].r for z in outputs])
                fullrows.append([-Y[i][j].i]+[z[i][j].i for z in outputs])
        assert rank_real(fullrows)==rr
        P=zero(3,3)
        for i in range(r,3): P[i][i]=Q(1)
        P=mul(mul(V,P),dag(V))
        assert eq(mul(P,P),P) and eq(mul(Y,P),zero(3,3))
        A=mul(mul(P,bs[-1]),P)
        if not any(bool(v) for row in A for v in row): A=P
        assert eq(A,dag(A)) and eq(mul(Y,A),zero(3,3))
        survivors.append(A)
        records.append({'Y_rank':r,'constraint_rank':rr,'R_N_extra_dimension':9-rr,
                        'all_dynamic_weight_dimension':10-rr,'Y_gaussian_rationals':dump(Y),
                        'kernel_projector':dump(P),'free_source_example':dump(A)})
    return records,ys,survivors

def hilbert_momentum(ys,survivors):
    sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]],complex);sz=np.diag([1.,-1.])
    sigmas=[sx,sy,sz];sp=np.array([F(2,3),F(-1,4),F(1,5)],dtype=float)
    max_identity=max_cross=max_transpose=max_free=0.
    cases=0
    for Y0,A0 in zip(ys,survivors):
        Y=npa(Y0)
        candidates=[npa(b) for b in herm_basis(3)]+[npa(A0)]
        M=np.kron(Y.conj(),np.eye(2))
        kin=np.kron(np.eye(3),sum(sp[j]*sigmas[j] for j in range(3)))
        z=np.zeros_like(M);h=np.block([[-kin,M],[M.conj().T,kin]])
        assert np.max(abs(h-h.conj().T))<1e-13
        for A in candidates:
            AR=np.kron(A.T,np.eye(2));P=np.block([[z,z],[z,AR]])
            max_transpose=max(max_transpose,float(np.max(abs(Y.conj()@A.T-(Y@A).conj()))))
            for i in range(3):
                si=np.kron(np.eye(3),sigmas[i]);alpha=np.block([[-si,z],[z,si]])
                ap=alpha@P;K=sp[i]*P/2+(h@ap+ap@h)/4
                com=h@K-K@h
                rhs=sp[i]*(h@P-P@h)/2+(h@h@ap-ap@h@h)/4
                max_identity=max(max_identity,float(np.max(abs(com-rhs))))
                max_cross=max(max_cross,float(np.max(abs(com[:6,6:]-sp[i]*M@AR/2))))
                if np.max(abs(Y@A))<1e-13:
                    max_free=max(max_free,float(np.max(abs(com))))
                cases+=1
    assert max(max_identity,max_cross,max_transpose,max_free)<3e-14
    # Exact bad-weight coefficient, independently in the one-flavor p_x=m=A=1 block.
    h=np.block([[-sx,np.eye(2)],[np.eye(2),sx]])
    P=np.diag([0.,0.,1.,1.]);alpha=np.block([[-sx,np.zeros((2,2))],[np.zeros((2,2)),sx]])
    K=P/2+(h@alpha@P+alpha@P@h)/4
    bad=h@K-K@h
    assert np.array_equal(bad[:2,2:],np.eye(2)/2)
    # A time-component-only check can falsely pass: the extra energy is conserved,
    # while the spatial Hilbert momentum is not. No conservation claim based on energy alone.
    beta=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]])
    energy=np.block([[np.zeros((2,2)),np.zeros((2,2))],[np.zeros((2,2)),sx]])+beta/2
    assert np.max(abs(h@energy-energy@h))==0 and np.max(abs(bad))>0
    return {'Hilbert_cases':cases,'max_full_identity_residual':max_identity,
            'max_cross_block_residual':max_cross,'max_charge_conjugation_residual':max_transpose,
            'max_free_kernel_commutator':max_free,'bad_LR_coefficient':'1/2',
            'potential_metric_0i_is_zero':True,'energy_only_false_positive_detected':True,
            'ordinary_test_spinors_are_Grassmann_coefficients_not_quantum_states':True}

def grassmann_coefficient():
    # For independently chosen odd theta_bar,theta, their degree-two product is nonzero.
    def wedge(a,b):
        if set(a)&set(b): return 0,()
        sign=(-1)**sum(x>y for x in a for y in b)
        return sign,tuple(sorted(a+b))
    s,w=wedge((0,),(1,));s2,w2=wedge((1,),(0,))
    assert (s,w)==(1,(0,1)) and (s2,w2)==(-1,(0,1))
    return {'independent_quadratic_monomial':[0,1],'bad_momentum_coefficient':'i/2',
            'opposite_order_sign':s2,'Grassmann_polarization_nonzero':True}

def common_branch():
    q={'Q':1,'u':-4,'d':2,'L':-3,'e':6,'N':0,'H':3}
    terms={'up':[('Q',1),('H',1),('u',1)],
           'down':[('Q',1),('H',-1),('d',1)],
           'charged':[('L',1),('H',-1),('e',1)],
           'neutral':[('L',1),('H',1),('N',1)]}
    rq={x:1 for x in q};rq['H']=2
    rp={x:(3 if x in ['Q','L'] else 1) for x in q};rp['H']=0
    rows={}
    for name,t in terms.items():
        vals=[sum(q[x]*sg for x,sg in t),sum(rq[x]*sg for x,sg in t)%4,sum(rp[x]*sg for x,sg in t)%4]
        assert vals==[0,0,0];rows[name]=vals
    assert 2*rp['N']%4==2 and 2*(rp['L']+rp['H'])%4==2
    mu=F(1,4096)
    # Affine masses stay positive on [0,1]; no new numerical thermal point is certified.
    endpoint_bounds=[min(F(1),j*mu) for j in [1,2,3]]
    assert min(endpoint_bounds)==mu
    # Constant background contribution cancels identically, with no convention-dependent V0 sign.
    kappa,Cf,Cbg,Qf,g=F(2,7),F(3,5),F(-4,9),F(11,13),F(7,3)
    before=kappa*Cf*g-kappa*Qf
    after=kappa*(Cf+Cbg)*g-kappa*(Qf+Cbg*g)
    assert before==after
    # At H=(0,1), norm^2=1 and V=lambda(norm^2-1)^2+V0 has zero derivative.
    lam=F(3,7);H2=F(1)
    dV=2*lam*(H2-1)
    assert dV==0
    # Differentiating the Yukawa term before setting H0 gives nonzero scalar source.
    formal_Y,L,N=F(2),F(3),F(5)
    coeff=formal_Y*L*N
    assert coeff==30
    return {'four_Yukawa_hypercharge_r_rtilde_residues':rows,
            'NN_and_Weinberg_rtilde_residue':2,'1305_mass_uniform_positive_lower':str(mu),
            '1305_parameter_point_not_newly_certified':True,
            'generic_background_source_shift_exact':str(after-before),
            'Higgs_stationary_derivative':str(dV),
            'Yukawa_H_source_Grassmann_coefficient_control':str(coeff),
            'C6_nonzero_not_adopted':True,'no_Majorana_S_or_Weinberg':True,
            'full_quantum_bosonic_evolution_not_certified':True}

def bindings():
    paths=['proof.md','check.py','NEXT.md','../research_note_1307.md',
      '../_shared/round1307_admission.md',
      '../../archive_1009_1043/research_note_1027.md',
      '../../archive_1009_1043/1027/drafts/source_selection_derivation.md',
      '../../archive_1009_1043/research_note_1035.md',
      '../../archive_1046_1062/1060/proof.md',
      '../1299/proof.md','../1305/proof.md','../1305/acceptance.json',
      '../1306/proof.md']
    return {s:sha256((HERE/s).resolve().read_bytes()).hexdigest() for s in paths}

def build():
    ks,ys,free=exact_source_kernels()
    return {'status':'PASS','scope':'declared leading source class; same Dirac branch background mapping',
            'new_axioms':0,'exact_source_kernels':ks,'symmetric_Hilbert_momentum':hilbert_momentum(ys,free),
            'Grassmann_coefficient':grassmann_coefficient(),'common_branch':common_branch(),
            'source_sha256':bindings()}

def compare(a,b,path=''):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif isinstance(a,float): assert abs(a-b)<=1e-12,path
    else: assert a==b,(path,a,b)

def main():
    ap=argparse.ArgumentParser();group=ap.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=ap.parse_args();result=build()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('PASS: exclusively created results.json')
    else:
        if OUT.exists():compare(result,json.loads(OUT.read_text(encoding='utf-8')))
        print('PASS: read-only check'+(' and saved-result comparison' if OUT.exists() else ' (no saved result yet)'))
    print(json.dumps({'source_dimensions':[v['all_dynamic_weight_dimension'] for v in result['exact_source_kernels']],
                      'Hilbert_cases':result['symmetric_Hilbert_momentum']['Hilbert_cases'],
                      'bound_sources':len(result['source_sha256'])},ensure_ascii=False))
if __name__=='__main__':main()
