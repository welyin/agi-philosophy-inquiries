"""Independent exact checks for 1306. No result claims a detector construction.
Default / --check read only. --write creates its result once.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse, hashlib, json
import numpy as np
HERE=Path(__file__).resolve().parent

def mul(a,b):
    return [[sum((x*y for x,y in zip(row,col)),Q(0)) for col in zip(*b)] for row in a]
def eye(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def transpose(a): return [list(x) for x in zip(*a)]
def power(a,n):
    r=eye(len(a))
    for _ in range(n): r=mul(r,a)
    return r
def determinant(a):
    b=[list(map(Q,row)) for row in a];d=Q(1);n=len(b)
    for j in range(n):
        k=next((i for i in range(j,n) if b[i][j]),None)
        if k is None:return Q(0)
        if k!=j:b[j],b[k]=b[k],b[j];d=-d
        pivot=b[j][j];d*=pivot
        for i in range(j+1,n):
            ratio=b[i][j]/pivot
            b[i]=[x-ratio*y for x,y in zip(b[i],b[j])]
    return d
def inverse(a):
    n=len(a);b=[list(map(Q,row))+unit for row,unit in zip(a,eye(n))]
    for j in range(n):
        k=next(i for i in range(j,n) if b[i][j]);b[j],b[k]=b[k],b[j]
        pivot=b[j][j];b[j]=[x/pivot for x in b[j]]
        for i in range(n):
            if i!=j:
                r=b[i][j];b[i]=[x-r*y for x,y in zip(b[i],b[j])]
    return [row[n:] for row in b]
def krylov(a):
    k=mul(a,a)
    return [[power(k,j)[i][0] for j in range(len(a))] for i in range(len(a))]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run_checks():
    MB=[[0,1,0,0],[1,1,1,0],[0,1,2,1],[0,0,1,3]]
    MA=[[0,1,0,0],[1,1,0,0],[0,0,2,0],[0,0,0,3]]
    CB,CA=krylov(MB),krylov(MA)
    assert determinant(MB)==-5 and determinant(CB)==27
    assert determinant(CA)==0
    assert all(CA[i][j]==0 for i in (2,3) for j in range(4))
    inv=inverse(CB);assert mul(CB,inv)==eye(4)
    K=mul(MB,MB)
    # Each j-th species row is recovered from the same first left field.
    for j in range(4):
        recovered=[sum(inv[k][j]*power(K,k)[0][i] for k in range(4)) for i in range(4)]
        assert recovered==eye(4)[j]
    mu=Q(1,2**30);a=Q(1,32)
    norm_M=Q(5)
    sigma_M=abs(determinant(MB))/norm_M**3
    assert sigma_M==Q(1,25)
    lead=sigma_M**3/8
    assert lead==Q(1,125000) # coefficient of W mu^3 / pi
    ratio=100*mu/lead
    assert ratio==Q(390625,33554432) and ratio<Q(1,2)
    x=(5*mu)**2/Q(16)
    assert x<Q(1,2) # exp(x) <= 1/(1-x)<2
    # Cauchy propagator interior: d_t S-alpha.grad S-i beta M S.
    # Each remainder series is bounded by its first term times exp(x).
    derivative_coefficient=Q(625)*Q(13,8)/64*2
    mass_coefficient=Q(3125,128)*mu
    assert derivative_coefficient+mass_coefficient<100
    w_lower_over_pi=a**3/12
    singular_lower=w_lower_over_pi*mu**3/Q(250000)
    assert singular_lower>0
    t0,eta,deadline=Q(15,32),Q(1,128),Q(3,4)
    tube_radius=a+eta
    assert t0+eta+tube_radius<deadline
    assert Q(7,16)<t0-eta<t0+eta<Q(1,2)
    # Independent spectral/Pauli diagnostic: 4 Weyl species x 2 spin x Nambu.
    M=np.asarray(MB,dtype=float)
    sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1.,-1.]);momentum=np.array([.21,-.13,.17])
    spin=momentum[0]*sx+momentum[1]*sy+momentum[2]*sz
    h=np.block([[-np.kron(M*0+np.eye(4),spin),np.kron(M,np.eye(2))],
                [np.kron(M,np.eye(2)),np.kron(M*0+np.eye(4),spin)]])
    expected=np.kron(np.eye(2),np.kron(M@M,np.eye(2))+np.dot(momentum,momentum)*np.eye(8))
    square_error=float(np.max(np.abs(h@h-expected)))
    assert square_error<1e-12
    vals=np.linalg.eigvalsh(M)
    # Common effect bound: ||Xi_l||,||Xi_m||<=C implies ||i[Xi_l,Xi_m]/(2C²)||<=1.
    # Use nonorthogonal Hermitian examples: never assume the A-branch modes canonical.
    common_C=Q(3)
    X=np.kron(sx,np.eye(2))+0.25*np.kron(sz,sx)
    Y=0.4*np.kron(sy,sz)+np.kron(np.eye(2),sy)
    effect=(np.eye(4)+1j*(X@Y-Y@X)/(2*float(common_C)**2))/2
    erange=np.linalg.eigvalsh(effect)
    assert np.linalg.norm(X,2)<float(common_C) and np.linalg.norm(Y,2)<float(common_C)
    assert erange[0]>=0 and erange[-1]<=1
    evidence={
       'status':'PASS_INDEPENDENT_WITH_EXPLICIT_ACCESS_SCOPE',
       'dirac_rank_both':1,
       'active_mass_det_B':str(determinant(MB)),
       'mass_squared_krylov_det_B':str(determinant(CB)),
       'mass_squared_krylov_det_A':str(determinant(CA)),
       'B_left_field_polynomial_coefficients_columns':[[str(x) for x in row] for row in inv],
       'A_N2_N3_all_left_dark':True,
       'strict_tail_to_leading_ratio_bound':str(ratio),
       'strict_singular_lower_bound':str(singular_lower),
       'strict_singular_lower_float':float(singular_lower),
       'uniform_kernel_remainder_over_W_mu4_pi_bound':str(derivative_coefficient+mass_coefficient),
       'deadline_margin':str(deadline-t0-eta-tube_radius),
       'Nambu_square_diagnostic_error':square_error,
       'signed_mass_diagnostic':vals.tolist(),
       'common_effect_diagnostic_range':erange.tolist(),
       'bounded_effects_do_not_certify_detector_implementation':True,
       'background_quadratic_only':True,
       'same_fixed_spacetime_modes_required_in_both_branches':True,
       'conditional_O_not_full_CC_or_M4':True,
    }
    evidence['code_sha256']=sha(Path(__file__).resolve())
    evidence['files_sha256']={
       'proof.md':sha(HERE/'proof.md'),
       '../_shared/round1306_admission.md':sha(HERE.parent/'_shared/round1306_admission.md'),
    }
    return evidence

def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group()
    g.add_argument('--write',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();result=run_checks();target=HERE/'independent_results.json'
    if args.write:
       with target.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
       print('CREATED',result['status'])
    else:
       saved=json.loads(target.read_text(encoding='utf-8'))
       assert result==saved,'Saved independent result mismatch'
       print('READ_ONLY_MATCH',result['status'])
if __name__=='__main__':main()
