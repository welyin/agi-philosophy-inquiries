"""860 working calibration: one mixed process and the full probe pullback.
The rational oscillator is not the original gravitational PDE. It checks
adjoints and source terms in R G_R R*, including a derivative reference.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'852'))
import native_outgoing_probe as old
TARGET=HERE/'compensated_probe_kernel_probe_results.json'

def A(rows):return np.array([[F(x) for x in row] for row in rows],dtype=object)
def zeros(a,b):return np.full((a,b),F(0),dtype=object)

def kernel(epsilon):
    nt=14;dt=F(1,4)
    # A positive, mutually coupled reference/probe process.
    K=A([[F(3,2),epsilon*F(2,3)],[epsilon*F(2,3),F(5,4)]])
    assert K[0,0]*K[1,1]-K[0,1]**2>0
    G=zeros(2*nt,2*nt)
    for s in range(1,nt-1):
        for field in range(2):
            z=zeros(nt,2)
            for t in range(1,nt-1):
                source=np.array([F(t==s and a==field) for a in range(2)],dtype=object)
                z[t+1]=2*z[t]-z[t-1]-dt**2*(K@z[t])+dt**2*source
            G[:,2*s+field]=z.reshape(-1)
    R=zeros(nt,2*nt);bare=zeros(nt,2*nt)
    d,e=F(3,7),F(2,5)
    for t in range(1,nt-1):
        R[t,2*t+1]=bare[t,2*t+1]=1
        R[t,2*t]=-epsilon*d
        R[t,2*(t+1)]-=epsilon*e/(2*dt)
        R[t,2*(t-1)]+=epsilon*e/(2*dt)
    Q=R@G@R.T
    # Source/response adjoint identity, with derivative stencils retained.
    z=np.array([F((i%7)-3,11) for i in range(2*nt)],dtype=object)
    f=np.array([F((i%5)-2,13) if 2<=i<=nt-3 else F(0) for i in range(nt)],dtype=object)
    assert f@(R@z)==(R.T@f)@z
    assert np.array_equal(Q.T,R@G.T@R.T)
    source=np.array([F(3,7) if i==3 else F(4,7) if i==4 else F(0) for i in range(nt)],dtype=object)
    later=np.array([F(i==9) for i in range(nt)],dtype=object)
    transfer=later@Q@source
    omitted=later@(R@G@bare.T)@source
    frozen=later@(bare@G@bare.T)@source
    # The finite derivative stencil has a one-step halo; only disjoint supports
    # outside that halo are used as a causality check here.
    assert all(Q[t,s]==0 for t in range(nt) for s in range(nt) if t<s-1)
    return transfer,omitted,frozen

def reference_identity():
    J=A([[2,1,0,1],[0,3,1,0],[1,0,2,1],[0,1,0,2]])
    Ji=np.array(old.inv(J.tolist()),dtype=object)
    dp=np.array([F(0),F(1,50),F(0),F(0)],dtype=object)
    grad=dp@Ji
    xi=np.array([F(2,3),F(-3,7),F(5,11),F(1,13)],dtype=object)
    assert dp@xi-grad@(J@xi)==0
    # If p itself is a full coordinate, calibrating f(X_new)=p fixes f=X_new^2.
    Jnew=J.copy();Jnew[2]=dp
    inverse=np.array(old.inv(Jnew.tolist()),dtype=object)
    newgrad=dp@inverse
    assert list(newgrad)==[0,0,1,0]
    return dict(old_reference_profile_gradient=[str(x) for x in grad],linear_gauge_residual='0',old_reference_off_shell_probe_derivative='1',new_coordinate_only_profile_gradient=[str(x) for x in newgrad],new_coordinate_only_subtraction_identically_zero=True)

def run():
    base=kernel(F(0))[0];assert base!=0
    rows=[]
    for ep in (F(1,50),F(1,25),F(1,10)):
        value,omitted,frozen=kernel(ep)
        assert kernel(-ep)==(value,omitted,frozen)
        assert value!=0 and value!=omitted and value!=frozen
        rows.append(dict(epsilon=str(ep),joint_transfer=str(value),joint_transfer_float=float(value),base_transfer=str(base),change_over_epsilon_squared=float((value-base)/ep**2),omitted_reference_source=str(omitted),omission_error_float=float(value-omitted),frozen_reference_transfer=str(frozen),exact_evenness=True))
    return dict(kind='round_860_working_calibration',formal_reports=859,new_numbered_scientific_groups=0,
        all_checks_passed=True,reference_identity=reference_identity(),rows=rows,
        derivative_reference_adjoint_checked=True,mixed_process_retained=True,
        original_curved_PDE_solved=False,new_Hadamard_family_constructed=False,
        full_quantum_source_and_receiver_transport_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
