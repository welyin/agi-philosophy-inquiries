"""Round1195: exact grading and original eight-variable DV kernel. Readonly by default."""
from pathlib import Path
from math import comb
from fractions import Fraction
import argparse,json
import numpy as np

def add(*polys):
    out={}
    for p in polys:
        for key,v in p.items(): out[key]=out.get(key,0)+v
    return {k:v for k,v in out.items() if v!=0}
def scale(p,c,charge=0,l1=0,l2=0):
    return {(m,q+charge,a+l1,b+l2):v*c for (m,q,a,b),v in p.items() if v*c!=0}
def sign_product(a,b):
    return -1 if sum((a>>(j+1)).bit_count() for j in range(8) if b>>j&1)%2 else 1
def mul(p,q):
    out={}
    for (m,c,a,b),v in p.items():
        for (n,d,e,f),w in q.items():
            if m&n: continue
            key=(m|n,c+d,a+e,b+f)
            out[key]=out.get(key,0)+v*w*sign_product(m,n)
    return {k:v for k,v in out.items() if v!=0}
def var(j):return {(1<<j,0,0,0):1}
def bilinear(sitebar,sitepsi,matrix,charge):
    out={}
    for r in range(2):
        for c in range(2):
            if matrix[r][c]!=0:
                term=mul(var(sitebar*4+2+r),var(sitepsi*4+c))
                out=add(out,scale(term,matrix[r][c],charge=charge(r,c)))
    return out
def original_kernel():
    s=(1,-1)
    sigmas=(((0,1),(1,0)),((0,-1j),(1j,0)))
    action={}
    for sigma in sigmas:
        eij=scale(add(bilinear(1,1,sigma,lambda r,c:s[c]-s[r]),
                       scale(bilinear(0,0,sigma,lambda r,c:0),-1)),1j)
        eji=scale(add(bilinear(0,0,sigma,lambda r,c:s[r]-s[c]),
                       scale(bilinear(1,1,sigma,lambda r,c:0),-1)),1j)
        fij=add(bilinear(0,1,sigma,lambda r,c:s[c]),
                 scale(bilinear(1,0,sigma,lambda r,c:-s[r]),-1))
        action=add(action,scale(add(mul(eij,eij),mul(eji,eji)),1/8,l1=1),
                   scale(mul(fij,fij),1/4,l2=1))
    assert all(m.bit_count()==4 for m,q,a,b in action)
    expS=add({(0,0,0,0):1},action,scale(mul(action,action),1/2))
    assert not mul(mul(action,action),action)
    return {key:v for key,v in expS.items() if key[1]==0}
def symbolic_matrix(kernel):
    basis=[0,3,5,10,12,15]
    matrix=[[{} for _ in basis] for _ in basis]
    for j,mask in enumerate(basis):
        contracted=scale(mul(kernel,{(mask,0,0,0):1}),-1)
        for (m,q,a,b),v in contracted.items():
            if m&15 !=15:continue
            target=m>>4
            assert target in basis
            i=basis.index(target)
            matrix[i][j][(a,b)]=matrix[i][j].get((a,b),0)+v
    return matrix
def expected_symbolic():
    E=[[{} for _ in range(6)] for _ in range(6)]
    E[0][0]={(1,0):1};E[5][5]={(1,0):1}
    E[0][5]={(0,0):-1};E[5][0]={(2,0):-2,(0,2):-2}
    for i,sg in [(1,-1),(2,1),(3,1),(4,-1)]:E[i][i]={(0,1):sg}
    return E
def evaluate(M,a,b):
    return np.array([[sum(v*a**i*b**j for (i,j),v in cell.items()) for cell in row] for row in M],complex)
def run():
    D=[sum(comb(16,j) for j in range(17) if j%12==r) for r in range(12)]
    shift=(-16-96)%12
    rbound=sum(min(D[r],D[(r+shift)%12]) for r in range(12))
    evenbound=sum(min(D[r],D[(r+shift)%12]) for r in range(0,12,2))
    assert rbound==24991 and evenbound==13951
    kernel=original_kernel();M=symbolic_matrix(kernel);E=expected_symbolic()
    assert M==E,(M,E)
    T=evaluate(M,1,.5).real;q=2.5
    W=np.diag([q,1,1,1,1,1]);sqrtW=np.diag(np.sqrt(np.diag(W)))
    A=sqrtW@T@np.linalg.inv(sqrtW)
    eigen=np.linalg.eigvalsh(A)
    expected=np.sort([1-np.sqrt(q),1+np.sqrt(q),-.5,-.5,.5,.5])
    H2e=-.5*np.log(eigen**2)
    signed=np.sign(eigen)
    reconstruction=float(np.max(abs(signed*np.exp(-H2e)-eigen)))
    trace_exact=[Fraction(8),Fraction(17),Fraction(179,4)]
    traces=[float(np.trace(np.linalg.matrix_power(T,n))) for n in (2,3,4)]
    assert all(abs(v-float(w))<1e-12 for v,w in zip(traces,trace_exact))
    C=np.zeros((6,6));basis=[0,3,5,10,12,15]
    for (mask,charge,a,b),v in kernel.items():
        C[basis.index(mask&15),basis.index(mask>>4)]+=float((v*1**a*.5**b).real)
    Q=np.zeros((6,6))
    for i,x in enumerate(basis):
        for j,y in enumerate(basis):
            if not x&y and x|y==15:Q[i,j]=sign_product(x,y)
    assert np.array_equal(T,-C.T@Q)
    verification={
        'symbolic_original_action_matches_matrix':True,
        'WT_minus_TadjW':float(np.linalg.norm(W@T-T.T@W)),
        'symmetrized_asymmetry':float(np.linalg.norm(A-A.T)),
        'spectrum_residual':float(np.max(abs(eigen-expected))),
        'signed_block_reconstruction':reconstruction,
        'berezin_pairing_eigenvalues':np.linalg.eigvalsh(Q).tolist()}
    assert max(verification[k] for k in ('WT_minus_TadjW','symmetrized_asymmetry','spectrum_residual','signed_block_reconstruction'))<1e-12
    return {
      'status':'PASS_WITH_SCOPE','check_groups':3,
      'wetterich_cell':{'boundary_variables_each':16,'integrated_variables':96,'action_degree':12,
        'grade_shift':shift,'grade_dimensions':D,'rank_upper_full_from_grade':rbound,
        'even_boundary_dimension':32768,'rank_upper_even':evenbound,'nullity_lower_even':32768-evenbound,
        'grade8_kernel_dimension_lower':comb(16,8)-comb(16,4)-comb(16,16),
        'scope':'unmodified single-cell boundary amplitude; not full physical quotient or multi-layer transfer'},
      'dv_original_kernel':{'couplings':{'lambda1':1,'lambda2':.5,'mu':0},
        'basis_masks':basis,'matrix':T.tolist(),'eigenvalues':eigen.tolist(),'negative_eigenvalues':int(sum(eigen<0)),
        'partition_traces_2_3_4':traces,'exact_traces':[str(x) for x in trace_exact],
        'minimum_modulus':float(min(abs(eigen))),'maximum_modulus':float(max(abs(eigen)))},
      'blocking':{'positive_metric_diagonal':np.diag(W).tolist(),'T2_eigenvalues':(eigen**2).tolist(),
        'H2_eigenvalues_at_a1':H2e.tolist(),'sign_involution_eigenvalues':signed.tolist(),
        'scope':'finite algebraic two-step positive reconstruction only; original reflection/observable positivity and real-time process not certified'},
      'verification':verification,
      'not_proved':['full physical spinor Hilbert space','same original quantum instruments or arbitrary-reference implementation','GR or dark sector','failure of all spinor routes']}
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    obj=run()
    if args.write_results:
        with Path(__file__).with_name('results.json').open('x',encoding='utf-8') as f:
            json.dump(obj,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(obj,ensure_ascii=False,indent=2))
