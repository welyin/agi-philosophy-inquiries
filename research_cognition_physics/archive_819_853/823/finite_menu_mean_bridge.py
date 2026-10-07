"""823: exact finite calibration of the local mean-matching construction.

This is a symplectic/Gaussian coefficient calculation, not a finite matrix
representation of the CCR and not a computation of the original curved W.
All algebraic checks use Fraction; no new dependencies or plots are needed.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json

HERE=Path(__file__).resolve().parent
TARGET=HERE/'finite_menu_mean_bridge_results.json'

def mat(rows):return [[F(x) for x in row] for row in rows]
def zeros(n,m):return mat([[0]*m for _ in range(n)])
def eye(n):return mat([[int(i==j) for j in range(n)] for i in range(n)])
def tr(a):return list(map(list,zip(*a)))
def mm(a,b):return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]
def add(a,b):return [[x+y for x,y in zip(r,s)] for r,s in zip(a,b)]
def neg(a):return [[-x for x in r] for r in a]
def sub(a,b):return add(a,neg(b))
def inv(a):
    n=len(a);m=[r.copy()+i for r,i in zip(a,eye(n))]
    for j in range(n):
        k=next(k for k in range(j,n) if m[k][j]);m[j],m[k]=m[k],m[j]
        pivot=m[j][j];m[j]=[x/pivot for x in m[j]]
        for i in range(n):
            if i!=j:
                factor=m[i][j];m[i]=[x-factor*y for x,y in zip(m[i],m[j])]
    return [r[n:] for r in m]
def pivots(a):
    a=[r.copy() for r in a];ans=[]
    for j in range(len(a)):
        v=a[j][j];assert v>0;ans.append(v)
        for i in range(j+1,len(a)):
            for k in range(j+1,len(a)):
                a[i][k]-=a[i][j]*a[j][k]/v
    return ans
def maxabs(a):return max(abs(x) for r in a for x in r)
def serial(a):return [[str(x) for x in r] for r in a]

def local_menu():
    n=6;lift=eye(n)+zeros(2,n)
    top=mat([[1,0],[0,1],[1,1],[0,-1],[1,0],[0,1]])
    K=top+eye(2)
    quotient=[r+[-x for x in t] for r,t in zip(eye(n),top)]
    skew=mat([[0,'1/3',0,0,'1/5',0],['1/4',0,0,0,0,'1/7']])
    L=add([([F(0)]*n)+r for r in eye(2)],mm(skew,quotient))
    assert mm(L,K)==eye(2)
    Pi=sub(eye(8),mm(K,L));assert mm(Pi,K)==zeros(8,2)
    J=zeros(n,n)
    for i in range(3):J[i][i+3]=1;J[i+3][i]=-1
    E=mm(mm(mm(mm(Pi,lift),J),tr(lift)),tr(Pi))
    assert tr(E)==neg(E)
    j=zeros(8,3);j[0][0]=1;j[1][1]=1;j[2][2]=1;j[3][2]=F(1,2)
    w=mm(E,j);g=mm(tr(Pi),w);S=mm(tr(w),w)
    positive=pivots(S)
    assert mm(mm(tr(g),E),j)==S
    # The four columns are the common term and the three Bloch coefficients.
    actual=mat([['1/3',2,-1,'1/2'],[0,'2/3',1,-2],[1,-1,'1/4',3]])
    target=mat([[0,0,0,0],[1,'1/2',0,0],[-1,0,'2/3',0]])
    d=sub(target,actual)
    controls=neg(mm(mm(g,inv(S)),d))
    shift=mm(E,controls)
    assert mm(tr(j),shift)==d
    assert mm(L,shift)==zeros(2,4)
    assert mm(tr(K),controls)==zeros(2,4)
    # A genuine three-mode Gaussian covariance V+iJ/2, checked modewise.
    variances=[F(1,2),F(3,4),F(1),F(1,2),F(3,4),F(1)]
    assert all(variances[i]*variances[i+3]>=F(1,4) for i in range(3))
    V=zeros(n,n)
    for i,v in enumerate(variances):V[i][i]=v
    Cov=mm(mm(mm(mm(Pi,lift),V),tr(lift)),tr(Pi))
    Q=mm(mm(tr(j),Cov),j);pivots(Q)
    costs=mm(mm(tr(controls),Cov),controls)
    lower=mm(mm(tr(d),inv(Q)),d)
    lower=[[x/4 for x in r] for r in lower]
    assert all(costs[i][i]>=lower[i][i] for i in range(4))
    # Redundant leading tests need not have redundant quadratic contacts.
    dep=mat([[1,0,0],[0,1,0],[0,0,1],[1,2,0],[-1,0,1]])
    contact=mat([[0,1,0,0],[0,0,0,1],[1,0,0,0],[0,0,1,0],[0,1,0,1]])
    true_actual=add(mm(dep,actual),contact)
    true_target=add(mm(dep,target),contact)
    corrected=add(true_actual,mm(dep,mm(tr(j),shift)))
    assert corrected==true_target
    relation=mat([[-1,-2,0,1,0],[1,0,-1,0,1]])
    assert mm(relation,dep)==zeros(2,3)
    good=sub(true_target,true_actual)
    assert mm(relation,good)==zeros(2,4)
    wrong=sub(mm(dep,target),true_actual)
    defect=mm(relation,wrong);assert maxabs(defect)>0
    # The same controls act affinely on any logical state; check two endpoints.
    bloch_a=mat([[1],['1/3'],['-1/4'],['1/5']])
    bloch_b=mat([[1],['-1/5'],['1/7'],['-1/3']])
    mix=[[F(2,5)*a[0]+F(3,5)*b[0]] for a,b in zip(bloch_a,bloch_b)]
    assert mm(controls,mix)==add([[F(2,5)*r[0]] for r in mm(controls,bloch_a)],
                                  [[F(3,5)*r[0]] for r in mm(controls,bloch_b)])
    return dict(exact_Gram_LDL_pivots=[str(x) for x in positive],
        matching_residual='0',gauge_residual='0',all_four_affine_components_checked=True,
        actual_menu_covariance=serial(Q),
        preparation_variances=[str(costs[i][i]) for i in range(4)],
        Robertson_lower_bounds=[str(lower[i][i]) for i in range(4)],
        logical_variance_sum=str(sum(costs[i][i] for i in range(1,4))),
        logical_variance_lower_sum=str(sum(lower[i][i] for i in range(1,4))),
        dropping_target_contact_relation_defect=serial(defect),
        finite_matrix_CCR_claimed=False)

def canonical_dictionary():
    constraint=mat([[1,2,1,0],[0,1,0,1]])
    right=mat([[0,0],[0,0],[-1,0],[0,-1]])
    assert mm(constraint,right)==neg(eye(2))
    q=mat([[1,'1/2',-1],[0,2,'1/3']])
    T=mat([[1,0,0,0],[0,1,0,0],[1,0,1,0],[0,1,'1/2',1]])
    k=mat([[0,1,0],[1,0,0],[0,'1/3',1],[1,0,'1/2']])
    ztarget=mm(right,q)
    u=mm(inv(T),sub(ztarget,k))
    qlag=add(q,mm(constraint,k))
    assert add(mm(mm(constraint,T),u),qlag)==zeros(2,3)
    assert add(mm(T,u),k)==ztarget
    wrong=mm(inv(T),ztarget)
    defect=add(mm(mm(constraint,T),wrong),qlag)
    assert maxabs(defect)>0
    assert ztarget[:2]==zeros(2,3)
    return dict(correct_constraint_residual='0',true_canonical_mean_residual='0',
        ignoring_canonical_contraction_constraint_defect=serial(defect),
        canonical_fixed_components=[0,1],
        diagnostic_not_original_canonical_contact_evaluation=True)

def run():
    return dict(round=823,all_checks_passed=True,finite_menu=local_menu(),
        canonical_dictionary=canonical_dictionary(),
        original_local_continuum_embedding_and_interpolation_are_analytic=True,
        original_Green_and_Hadamard_matrices_numerically_evaluated=False,
        original_actual_mean_coefficients_numerically_evaluated=False,
        common_preparation_for_each_fixed_finite_menu_proven_analytically=True,
        one_preparation_matching_every_menu_or_whole_initial_slice_proven=False,
        exact_nonlinear_canonical_contact_of_original_model_computed=False,
        formal_test_groups_added=1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite frozen results.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
