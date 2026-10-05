"""799: exact constrained Cauchy and mixed Wick representation diagnostics.

The finite recurrences calibrate the identities of the continuous proof.
They are not a discretization of the original interacting quantum gravity model.
"""
from pathlib import Path
from fractions import Fraction as Q
import importlib.util
import json
import sys
sys.dont_write_bytecode=True
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'798'))
import internal_register_probe as base
C,cmat,cscale,dag,csame=base.C,base.matrix,base.scale,base.dag,base.same
def mat(a): return np.array([[Q(x) for x in row] for row in a],dtype=object)
def eye(n): return mat([[i==j for j in range(n)] for i in range(n)])
def zero(m,n): return mat([[0]*n for _ in range(m)])
def same(a,b): return all(x==0 for x in (a-b).flat)
def cq(a): return cmat(a.tolist())
def defect(a):
    i,j=next((i,j) for i in range(a.shape[0]) for j in range(a.shape[1]) if a[i,j])
    return dict(row=i,column=j,value=str(a[i,j]))
def wick4(w):
    return w[0,1]*w[2,3]+w[0,2]*w[1,3]+w[0,3]*w[1,2]

def constrained_cauchy():
    steps=[mat([[1,Q(1,3)],[0,1]]),mat([[1,0],[Q(-2,5),1]]),
           mat([[1,Q(-3,7)],[0,1]]),mat([[1,0],[Q(4,9),1]])]
    ns=5
    J=mat([[0,1],[-1,0]])
    us=[eye(2)]
    for m in steps:
        assert same(m.T@J@m,J)
        us.append(m@us[-1])
    U=np.vstack(us)
    D=zero(8,10)
    for t,m in enumerate(steps):
        D[2*t:2*t+2,2*t:2*t+2]=-m
        D[2*t:2*t+2,2*t+2:2*t+4]=eye(2)
    assert same(D@U,zero(8,2))
    diff=zero(ns,ns)
    for t in range(ns-1):
        diff[t,t]=-Q(t+2,3);diff[t,t+1]=Q(t+2,3)
    kp=zero(10,5)
    mix=zero(5,10)
    for t in range(ns):
        kp[2*t,:]=diff[t,:]
        kp[2*t+1,:]=diff[t,:]*Q(t+1,5)
        mix[t,0::2]=diff[t,:]*Q(1,4)
        mix[t,2*t+1]+=Q(1,7)
    K=np.vstack((kp,eye(5)))
    F=np.hstack((eye(10),-kp))
    E=np.vstack((eye(10),zero(5,10)))
    L=np.hstack((zero(5,10),eye(5)))+mix@F
    PI=eye(15)-K@L
    P=F.T@D.T@D@F
    assert same(L@K,eye(5)) and same(PI@PI,PI)
    assert same(P@PI,P) and same(L@PI,zero(5,15))
    Ufull=PI@E@U
    assert same(P@Ufull,zero(15,2))
    E0=zero(10,2);E0[:2,:]=eye(2)
    localizer=F.T@E0@U.T@E.T@PI.T
    assert same(K.T@localizer,zero(5,15))
    assert same(localizer@localizer,localizer)
    assert same(Ufull.T@localizer,Ufull.T)
    # Construct every equation primitive, rather than only testing covariance.
    for j in range(15):
        source=eye(15)[:,j:j+1]
        f=E.T@PI.T@source
        rem=f-E0@U.T@f
        g=zero(8,1);g[6:8]=rem[8:10]
        for t in range(3,0,-1):
            g[2*(t-1):2*t]=rem[2*t:2*t+2]+steps[t].T@g[2*t:2*t+2]
        assert same(D.T@g,rem)
        y=zero(10,1)
        for t,m in enumerate(steps):
            y[2*t+2:2*t+4]=m@y[2*t:2*t+2]+g[2*t:2*t+2]
        primitive=PI@E@y
        assert same(L@primitive,zero(5,1))
        assert same(PI.T@source-localizer@source,P@primitive)
    # Same fixed quantum covariance; no new state at the early slice.
    covariance=cmat([[Q(3,4),C(0,Q(1,2))],[C(0,Q(-1,2)),1]])
    W=cq(Ufull)@covariance@cq(Ufull.T)
    T=cq(localizer)
    assert csame(dag(T)@W@T,W)
    # Mixed finite kernels: localization applied in each Wick slot.
    probes=[]
    for n in range(6):
        probes.append(mat([[Q(((i+1)*(n+2))%9-4,7)] for i in range(15)]))
    wcount=0
    for shift in range(3):
        tests=np.hstack(probes[shift:shift+4])
        wp=cq(tests.T)@W@cq(tests)
        local=localizer@tests
        wl=cq(local.T)@W@cq(local)
        assert wick4(wp)==wick4(wl)
        wcount+=1
    chi=zero(15,15)
    values=[0,0,Q(1,3),1,1]
    for t,x in enumerate(values):
        chi[2*t,2*t]=x;chi[2*t+1,2*t+1]=x;chi[10+t,10+t]=x
    bad=L@chi@Ufull
    assert any(bad.flat)
    assert same(L@PI@chi@Ufull,zero(5,2))
    # Omitting dual projection leaves a gauge-dependent source.
    wrong=K.T@(eye(15)-PI.T+localizer)
    assert any(wrong.flat)
    support=[i for i in range(15) if any(localizer[i,:])]
    assert support==[0,1,10,11],support
    return dict(nonstationary_steps=4,field_dimension=15,physical_initial_data_dimension=2,
                source_equation_primitives_checked=15,
                all_constrained_localizer_identities_exact=True,
                localizer_support_rows=support,
                same_quantum_two_point_kernel_preserved=True,
                Wick_four_point_checks=wcount,
                cut_field_without_projection_defect=defect(bad),
                omit_dual_source_projection_defect=defect(wrong),
                incoming_algebra_equals_finite_record_M2=False)

def mixed_records():
    a,b=base.car(2)
    I=base.eye(4)
    n=dag(a)@a
    pair=dag(a)@dag(b)
    even=base.proj(4,0)+base.proj(4,3)
    # Quadratic parity-even finite evolutions, not the original full H.
    v1=I-even+cscale(even,Q(3,5))+cscale(pair+dag(pair),C(0,Q(4,5)))
    v2=I+cscale(n,C(Q(-2,5),Q(4,5)))
    v3=I-even+cscale(even,Q(5,13))+cscale(pair+dag(pair),C(0,Q(12,13)))
    steps=[v1,v2,v3]
    us=[I]
    for v in steps:
        assert csame(dag(v)@v,I)
        us.append(v@us[-1])
    operators=[dag(us[1])@a@us[1],dag(us[2])@b@us[2],
               dag(us[3])@dag(b)@us[3],dag(us[3])@dag(a)@us[3]]
    rho=cmat([[Q(21,50),0,0,0],[0,Q(14,50),0,0],
              [0,0,Q(9,50),0],[0,0,0,Q(6,50)]])
    def avg(x): return np.trace(rho@x)
    ws=[[avg(x@y) for y in operators] for x in operators]
    val=avg(operators[0]@operators[1]@operators[2]@operators[3])
    wick=ws[0][1]*ws[2][3]-ws[0][2]*ws[1][3]+ws[0][3]*ws[1][2]
    assert val==wick
    wrong=ws[0][1]*ws[2][3]+ws[0][2]*ws[1][3]+ws[0][3]*ws[1][2]
    # A second ordering detects the fermionic sign if this covariance entry vanishes.
    op=[operators[0],operators[3],operators[1],operators[2]]
    ww=[[avg(x@y) for y in op] for x in op]
    fv=avg(op[0]@op[1]@op[2]@op[3])
    ff=ww[0][1]*ww[2][3]-ww[0][2]*ww[1][3]+ww[0][3]*ww[1][2]
    assert fv==ff
    sign_defect=ww[0][1]*ww[2][3]+ww[0][2]*ww[1][3]+ww[0][3]*ww[1][2]-fv
    assert sign_defect
    record=dag(us[-1])@n@us[-1]
    assert csame(record@record,record) and csame(dag(record),record)
    p=avg(record)
    assert not p.im and 0<=p.re<=1
    assert p!=avg(n)
    # The algebra of n alone cannot express this actual pulled-back record.
    offdiag=record@ n-n@record
    assert any(offdiag.flat)
    return dict(time_dependent_quadratic_CAR_steps=3,
                same_fixed_preparation_used=True,ordered_CAR_Wick_checks=2,
                ordered_four_point_value=str(val),
                omitted_fermion_sign_defect=str(sign_defect),
                original_occupation_probability=str(avg(n)),
                pulled_back_record_probability=str(p),
                finite_initial_occupation_menu_insufficient=True,
                extra_CAR_phase_space_used_not_extra_particle=True,
                continuous_native_H_simulated=False)

def run():
    return dict(round=799,all_checks_passed=True,
                constrained_cauchy=constrained_cauchy(),mixed_record=mixed_records(),
                finite_tests_establish_continuum_theorem=False,
                original_interacting_early_local_algebra_time_slice_proven=False,
                autonomous_preparation_or_readout_proven=False)
if __name__=='__main__':
    out=run()
    (HERE/'constrained_cauchy_records_results.json').write_text(
        json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(out,ensure_ascii=False,indent=2))

