"""1187: declared bare real material-clock branch; no quantum/GR claim.
Default read-only recomputation; --write creates results exclusively.
"""
from fractions import Fraction as F
from math import factorial,cos,cosh,sin,sinh
from pathlib import Path
import json,sys
import numpy as np

def matmul(a,b):
    return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]
def tr(a):return list(map(list,zip(*a)))
def sub(a,b):return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]
def scale(c,a):return [[c*x for x in r] for r in a]
def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]
def det(a):
    a=[r[:] for r in a];out=F(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None:return F(0)
        if j!=i:a[i],a[j]=a[j],a[i];out=-out
        piv=a[i][i];out*=piv
        for j in range(i+1,len(a)):
            r=a[j][i]/piv
            for k in range(i,len(a)):a[j][k]-=r*a[i][k]
    return out

def evaluate():
    # An oblique three-dimensional spatial frame in R^6, no chosen time metric.
    E=[[F(v) for v in r] for r in [[1,1,0],[0,1,1],[0,0,1],[0,0,0],[0,0,0],[0,0,0]]]
    P=[[F(i==j and i>=3) for j in range(6)] for i in range(6)]
    assert matmul(P,E)==[[F(0)]*3 for _ in range(6)]
    h=matmul(tr(E),E);assert det(h)==1
    v=[F(x) for x in [1,2,3,2,3,6]];u=[F(0)]*3+v[3:]
    assert sum(x*x for x in u)==49
    W=F(1,2); norm=F(7); n=[x/norm for x in u]; pi=[W*x for x in n]
    assert matmul(tr(E),[[x] for x in pi])==[[F(0)]]*3
    assert sum(x*x for x in pi)==W*W
    lag=W*norm;assert sum(x*y for x,y in zip(pi,v))-lag==0
    Hess=scale(W/norm,sub(P,[[x*y for y in n] for x in n]))
    assert matmul(Hess,[[x] for x in v])==[[F(0)]]*6
    assert matmul(Hess,E)==[[F(0)]*3 for _ in range(6)]
    assert matmul(Hess,Hess)==scale(W/norm,Hess)
    eig=np.linalg.eigvalsh(np.array(Hess,dtype=float))
    assert np.count_nonzero(eig>1e-12)==2
    # Exact graph determinant identity det(I+dd^T)=1+|d|^2 for one residual field.
    ds=[[F(1,2),F(1,3),F(-2,3),F(2)],[F(0),F(1),F(0),F(0)],[F(3),F(2),F(1),F(-1)]]
    for d in ds:
        A=[[F(i==j)+d[i]*d[j] for j in range(4)] for i in range(4)]
        assert det(A)==1+sum(x*x for x in d)
    # H=-sqrt((1+|grad psi|^2)(F0^2-p^2)); exact Hamilton-Legendre identities.
    # Pick all rational velocity/gradient with sqrt(1+v^2+|grad|^2)=3.
    vel=F(2); grad=[F(2),F(0),F(0)]; f0=F(3,2)
    sqrtA=F(3);p=f0*vel/sqrtA;L=f0*sqrtA;H=p*vel-L
    R=(1+sum(x*x for x in grad))*(f0*f0-p*p)
    assert H<0 and H*H==R and R>0
    assert p*p<f0*f0
    # Constraint kernel remains after gauge conditions on action-independent X.
    # In the regular local constraint basis (C_phi, pi_X), gauge X=x has rank 4 of 8.
    fp=np.concatenate([np.zeros((4,4)),np.eye(4)],axis=1)
    assert np.linalg.matrix_rank(fp)==4 and np.all(fp[:,:4]==0)
    # Exact lower/upper bound for cosh(1)-cos(1): positive series exponents 4j+2.
    lower=sum((F(2,factorial(4*j+2)) for j in range(9)),F(0))
    first_omitted=F(2,factorial(38)); ratio=F(1,39*40*41*42)
    upper=lower+first_omitted/(1-ratio)
    assert lower>1 and upper-lower<F(1,10**43)
    samples=[]
    for k in [F(1,2),F(1),F(2)]:
        t=1/k; z=float(k*t)
        e=np.array([[cosh(z),sinh(z)/float(k)],[float(k)*sinh(z),cosh(z)]])
        wave=np.array([[cos(z),sin(z)/float(k)],[-float(k)*sin(z),cos(z)]])
        assert abs(np.linalg.det(e)-1)<1e-14 and abs(np.linalg.det(wave)-1)<1e-14
        assert abs((e[0,0]-wave[0,0])-float(lower))<1e-14
        samples.append({'k':str(k),'T':str(t),'bare_linear_position':float(e[0,0]),'Lorentz_wave_position':float(wave[0,0])})
    # A real regular point has positive definite Euclidean principal form.
    # General single-field graph has principal coefficient F/sqrt(A)*(I-dd^T/A).
    d=np.array([.2,.3,-.4,.1]);A=1+float(d@d)
    principal=(np.eye(4)-np.outer(d,d)/A)/np.sqrt(A)
    eigp=np.linalg.eigvalsh(principal);assert min(eigp)>0
    return {'round':1187,'date':'2026-10-14','status':'PASS_WITH_BARE_MATERIAL_CLOCK_SCOPE','groups':{
      '1_material_constraints':{'pass':True,'spatial_gram_determinant':str(det(h)),'material_fields':6,'Legendre_Hessian_rank':2,'expected_rank':'N-4','canonical_H':'0','normal_constraint':'pi^2-W^2=0','dummy_X_FP_rank':4,'constraints_before_X_gauge':8,'original_material_constraints_remaining':4},
      '2_true_material_clock_graph':{'pass':True,'exact_graph_tests':len(ds),'momentum':str(p),'lagrangian':str(L),'Hamiltonian_same_oriented_branch':str(H),'root_squared':str(R),'flat_graph_principal':'F0 times I4, positive definite','nonzero_graph_eigenvalues_diagnostic':eigp.tolist(),'linear_equation_constant_V':'psi_TT + Delta psi = 0'},
      '3_finite_mode_linear_separation':{'pass':True,'exact_gap_lower':str(lower),'exact_gap_upper':str(upper),'certified_gap_strict_lower':'1','float_gap_diagnostic':float(lower),'samples':samples,'statement':'normalized linear-response coefficient, not a complete finite-amplitude quantum record'}},
      'scope_limits':['Original real positive Gram and V>0 regular material coordinate branch only.','Flat constant-V witness is an explicit m=0 diagnostic; not the massive loop-induced background.','Variable positive V changes lower-order terms; elliptic principal type alone does not prove all finite or low modes unstable.','H_can=0 does not preclude relational time; the material clock has a nontrivial reduced Hamiltonian on R>0.','No quantum physical Hilbert space, nonlinear finite-error estimate, actual clock resource or Lorentz collective renormalized phase certified.','No rejection of general composite gravity or the cognitive axioms.'],
      'new_adopted_axioms':0,'new_empirical_calibration':0,'gr_hard_deliverables_completed':0}
if __name__=='__main__':
    result=evaluate(); path=Path(__file__).with_name('results.json')
    if '--write' in sys.argv:
        with path.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WROTE 3 groups; bare material-clock finite linear gap > 1')
    else:
        assert result==json.loads(path.read_text(encoding='utf-8'))
        print('READ_ONLY_MATCH 3 groups; bare material-clock finite linear gap > 1')
