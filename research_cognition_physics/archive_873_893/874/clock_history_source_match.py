"""874: original source-label clock reduction, exact measure identity and records.
The field statement is a local classical/formal-tree identity. Finite CAR
history calculations calibrate a declared spectral ordering, not full QG.
"""
from pathlib import Path
from fractions import Fraction as F
from itertools import permutations
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'873'))
import clock_mass_source_bridge as algebra
import constraint_clock_kernel_probe as probe
TARGET=HERE/'clock_history_source_match_results.json'
add,scale,mul,const=algebra.add,algebra.scale,algebra.mul,algebra.const
def ja(*xs):return [add(*(x[i] for x in xs)) for i in range(3)]
def js(x,c):return [scale(a,c) for a in x]
def jm(x,y):return [add(*(mul(x[j],y[i-j]) for j in range(i+1))) for i in range(3)]
def jc(x):return [x,{},{}]
def power(x,n):
    z=jc(const(1))
    for _ in range(n):z=jm(z,x)
    return z
def exact_measure():
    a,b,p0=F(2),F(1,3),F(2,5);v=2*a*p0+b
    ar={(0,0,5):F(1,13)}
    br={(0,0,3):F(2,9),(0,0,10):F(-1,8)}
    cr={(0,0,12):F(4,9),(0,0,6):F(1,6),(0,0,15):F(5,17)}
    mass={(1,0,3):F(3,7),(0,1,12):F(5,11)}
    # z=C0(pi): pi(z)=p0+z/v-a z^2/v^3.
    pi=[const(p0),const(1/v),const(-a/v**3)]
    inverse_speed=[const(1/v),const(-2*a/v**3),const(6*a*a/v**5)]
    q=ja(jm(jc(ar),jm(pi,pi)),jm(jc(br),pi),jc(cr),jc(mass))
    qpi=ja(js(jm(jc(ar),pi),2),jc(br))
    jac=ja(jc(const(1)),jm(qpi,inverse_speed))
    # Independent original implicit root, retaining all nilpotent terms.
    dp={}
    for _ in range(6):
        pp=add(const(p0),dp)
        nonlinear=add(scale(mul(dp,dp),a),mul(ar,mul(pp,pp)),mul(br,pp),cr,mass)
        nd=scale(nonlinear,-1/v)
        if nd==dp:break
        dp=nd
    else:raise AssertionError('implicit root')
    root=add(const(p0),dp)
    rows=[]
    for degree in range(5):
        test=power(pi,degree)
        numerator=jm(jac,test)
        # delta(z+Q)=delta(z)+Q delta'(z)+Q^2 delta''(z)/2;
        # distribution pairing gives coefficient0 - coefficient1 + coefficient2.
        aq=jm(numerator,q);aqq=jm(aq,q)
        integral=add(numerator[0],scale(aq[1],-1),aqq[2])
        expected=const(1)
        for _ in range(degree):expected=mul(expected,root)
        assert integral==expected
        rows.append(dict(test_clock_momentum_power=degree,exact_residual_zero=True,
            mixed_source_four_leg_coefficient=str(integral.get((1,1,15),F(0)))))
    # Missing full C_pi leaves a density, already visible for F=1.
    A=inverse_speed
    wrong=add(A[0],scale(jm(A,q)[1],-1),jm(jm(A,q),q)[2])
    assert wrong!=const(1)
    # Remove only common background v: it does not fix the source-dependent density.
    normalized_wrong=scale(wrong,v)
    assert normalized_wrong!=const(1)
    assert normalized_wrong.get((1,1,15),F(0))!=0
    # A nontrivial finite spatial gauge block with nilpotent clock factor.
    speed=add(const(v),scale(dp,2*a),scale(mul(ar,root),2),br)
    J=[[F(2),F(1,3),F(0)],[F(-1,4),F(3),F(1,5)],[F(1,7),F(0),F(4)]]
    def det(M):
        total={}
        for perm in permutations(range(len(M))):
            inv=sum(perm[i]>perm[j] for i in range(len(M)) for j in range(i+1,len(M)))
            term=const((-1)**inv)
            for i,j in enumerate(perm):term=mul(term,M[i][j])
            total=add(total,term)
        return total
    detj=det([[const(x) for x in row] for row in J])[(0,0,0)]
    full=[[speed,{}, {},{}]]
    for i in range(3):
        full.append([add(const(F(i+1,6)),cr)]+[const(x) for x in J[i]])
    assert det(full)==scale(speed,detj)
    return dict(exact_test_moment_rows=rows,spatial_J_determinant=str(detj),
        spatial_and_clock_determinants_factor_exactly=True,
        nilpotent_clock_and_spin_like_momentum_terms_retained=True,
        background_normalized_missing_FP_mixed_coefficient=str(normalized_wrong[(1,1,15)]),
        all_full_measure_residuals_zero=True)
def finite_records():
    bg=json.loads((HERE.parent/'861/magnetic_reduced_hamiltonian_results.json').read_text('utf-8'))['original_859_background_clock_reduction'][1]
    a,v=bg['a'],bg['actual_clock_speed'];qstar=v*v/(4*a);hstar=v/(2*a);time=1/hstar
    nodes,weights=np.polynomial.legendre.leggauss(64);weights*=.75*(1-nodes*nodes)
    Y=np.array([[0,-1j],[1j,0]],complex);rot=np.cos(.61)*np.eye(2)-1j*np.sin(.61)*Y
    base=[qstar*np.diag([.08,.12]),qstar*rot@np.diag([.04,.13])@rot.conj().T]
    dirs=[qstar*np.diag([1.,0.]),qstar*np.diag([0.,1.])]
    c=[probe.annihilation(i) for i in range(2)]
    n=[x.conj().T@x for x in c];I=np.eye(4)
    px=.5*(n[0]+n[1]+c[0].conj().T@c[1]+c[1].conj().T@c[0])
    assert np.linalg.norm(px@px-px)<1e-14
    rho=np.zeros((4,4),complex);rho[0,0]=rho[3,3]=.1
    rho[1:3,1:3]=np.array([[.8,.5],[.5,.8]])/2
    def kernel(B,width=0,FP=True):
        ev,Q=np.linalg.eigh(B)
        if width:
            z=ev[:,None]-width*qstar*nodes
            speed=np.sqrt(v*v-4*a*z);f=2*z/(v+speed)
            vals=np.exp(-1j*time*f)
            if not FP:vals*=v/speed
            vals=vals@weights
        else:
            speed=np.sqrt(v*v-4*a*ev);f=2*ev/(v+speed)
            vals=np.exp(-1j*time*f)
            if not FP:vals*=v/speed
        return (Q*vals)@Q.conj().T
    def histories(s=0.,t=0.,width=0.,FP=True):
        B1=probe.lift(base[0]+s*dirs[0]+t*dirs[1])
        B2=probe.lift(base[1])
        U1,U2=kernel(B1,width,FP),kernel(B2,width,FP)
        values=[]
        for P in (I-px,px):
            for Q in (I-n[0],n[0]):
                C=Q@U2@P@U1
                values.append(float(np.trace(C@rho@C.conj().T).real))
        return np.array(values)
    exact=histories();assert exact.min()>0 and abs(sum(exact)-1)<1e-14
    wrong=histories(FP=False);assert abs(sum(wrong)-1)>.05
    rows=[]
    for width in (.04,.02,.01):
        p=histories(width=width)
        rows.append(dict(delta_width_fraction=width,probability_error=float(max(abs(p-exact))),
                         normalization_defect=float(abs(sum(p)-1))))
    assert rows[-1]['probability_error']<rows[0]['probability_error']/12
    eps=1e-4
    def mixed(width=0.):
        return (histories(eps,eps,width)-histories(eps,-eps,width)-histories(-eps,eps,width)+histories(-eps,-eps,width))/(4*eps*eps)
    desired=mixed();assert max(abs(desired))>1e-3
    errors=[]
    for width in (.04,.02,.01):errors.append(float(max(abs(mixed(width)-desired))))
    assert errors[-1]<errors[0]/10
    return dict(exact_joint_history_probabilities=exact.tolist(),
        original870_preparation_m=.5,original870_preparation_r=.8,
        middle_record='occupation in plus-x mode',final_record='occupation in mode1',
        exact_probability_sum=float(sum(exact)),regularized_history_rows=rows,
        mixed_mass_record_source=desired.tolist(),mixed_source_errors=errors,
        missing_FP_probability_sum=float(sum(wrong)),
        scope='declared finite CAR spectral clock history; not original field probabilities')
def run():
    previous=probe.run()
    assert previous==json.loads(probe.TARGET.read_text('utf-8'))
    return dict(round=874,date='2026-10-06',formal_reports=874,
        fresh_numbered_groups=1,cumulative_numbered_groups=3659,all_checks_passed=True,
        previous_working_probe_reproduced=True,exact_constraint_measure=exact_measure(),
        finite_common_record_and_source=finite_records(),
        analytic_scope='Original regular h/Y chart, canonical measure and mass source family, classical formal coefficients; the finite CAR ordered history identity is separately exact for its declared spectral ordering.',
        original873_mixed_source_recovered=True,
        fermion_second_class_and_spin_mass_independent_factors_retained=True,
        all_original_field_quantum_measures_matched=False,
        original_finite_coupling_or_graph_continuum_bridge_completed=False,
        full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
