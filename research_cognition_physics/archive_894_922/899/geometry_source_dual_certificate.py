"""899: dual geometry/source error and joint parameter jets.
Exact rational calibration reuses the full mixed 860 oscillator and all four
872 Euler source coefficients. It is NOT the original curved PDE computation.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, json, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'860'))
import compensated_probe_kernel_probe as old
TARGET=HERE/'geometry_source_dual_certificate_results.json'
def zero(m,n):return np.full((m,n),F(0),dtype=object)
def eye(n):return np.array([[F(i==j) for j in range(n)] for i in range(n)],dtype=object)
def maxabs(a):return max(map(abs,a.flat),default=F(0))
def as_scalar(v):return dict(exact=str(v),value=float(v))
def build(ep):
    nt=14;dt=F(1,4);n=2*(nt-2)
    k=np.array([[F(3,2),ep*F(2,3)],[ep*F(2,3),F(5,4)]],dtype=object)
    dk=np.array([[F(0),F(2,3)],[F(2,3),F(0)]],dtype=object)
    # Unknown slots y[2],...,y[13]; equations at t=1,...,12.
    A=zero(n,n);dA=zero(n,n)
    for t in range(1,nt-1):
        row=2*(t-1);A[row:row+2,row:row+2]=eye(2)/dt**2
        if t>=2:
            col=2*(t-2);A[row:row+2,col:col+2]=k-2*eye(2)/dt**2
            dA[row:row+2,col:col+2]=dk
        if t>=3:A[row:row+2,2*(t-3):2*(t-3)+2]=eye(2)/dt**2
    # The retarded inverse is calculated by the same recurrence as frozen860.
    G=zero(n,n)
    for source_index in range(n):
        yy=zero(nt,2)
        for t in range(1,nt-1):
            force=np.array([F(source_index==2*(t-1)+j) for j in range(2)],dtype=object)
            yy[t+1]=2*yy[t]-yy[t-1]-dt**2*(k@yy[t])+dt**2*force
        G[:,source_index]=yy[2:].reshape(-1)
    assert np.array_equal(A@G,eye(n)) and np.array_equal(G@A,eye(n))
    Rfull=zero(nt,2*nt);Rpfull=zero(nt,2*nt)
    for t in range(1,nt-1):
        Rfull[t,2*t+1]=1
        Rpfull[t,2*t]=-F(3,7)
        Rpfull[t,2*(t+1)]=-F(2,5)/(2*dt)
        Rpfull[t,2*(t-1)]=F(2,5)/(2*dt)
    Rfull+=ep*Rpfull
    # Source endpoint: R^T acts at force slots t=1,...,12.
    R=Rfull[:,2:2*(nt-1)];Rp=Rpfull[:,2:2*(nt-1)]
    pulse=np.array([F(3,7) if t==3 else F(4,7) if t==4 else F(0) for t in range(nt)],dtype=object)
    assert np.all((Rfull.T@pulse)[:2]==0) and np.all((Rfull.T@pulse)[-2:]==0)
    # Frozen860's probe readout and its original transfer are unchanged.
    original_read=Rfull[9,4:]
    original=original_read@G@R.T@pulse
    assert original==old.kernel(ep)[0]
    # An unrenormalized field-0 readout, not the original spacetime metric.
    ell=np.array([F(t==9 and c==0)+F(t==10 and c==0)/3 for t in range(2,nt) for c in range(2)],dtype=object)
    return A,dA,G,R,Rp,pulse,ell,original

def run():
    saved=json.loads((HERE.parent/'872/read_vertex_source_closure_results.json').read_text('utf-8'))
    parts=[F(p['15']) for p in saved['four_leg_source_calibration']['source_parts']]
    total=sum(parts,F(0));assert total==F(83,2250)
    ep=F(1,50);A,dA,G,R,Rp,pulse,ell,original=build(ep)
    z=G.T@ell
    zp=-G.T@dA.T@z
    zpp=-2*G.T@dA.T@zp
    qs=[R@z,Rp@z+R@zp,2*Rp@zp+R@zpp]
    J=total*pulse;s=R.T@J;sp=Rp.T@J
    y=G@s;yp=G@(sp-dA@y);ypp=G@(-2*dA@yp)
    direct=[ell@v for v in (y,yp,ypp)]
    dual=[q@J for q in qs]
    assert direct==dual and all(v!=0 for v in direct)
    individual=[parts[i]*(qs[0]@pulse) for i in range(4)]
    assert sum(individual,F(0))==direct[0]
    omissions=[direct[0]-individual[i] for i in range(4)]
    assert all(v!=direct[0] for v in omissions)
    # Freeze geometry/propagation pullback while differentiating parameters:
    # this loses all q jets in this fixed-source calibration.
    rows=[]
    for denominator in (100,1000,10000):
        for k,q in enumerate(qs):
            qhat=np.array([F(round(x*denominator),denominator) for x in q],dtype=object)
            query=qhat@J;error=abs(query-direct[k]);eta=maxabs(q-qhat)
            bound=sum(map(abs,J),F(0))*eta
            assert error<=bound
            rows.append(dict(jet_order=k,rounding_denominator=denominator,
                predicted=as_scalar(query),error=as_scalar(error),test_error_sup=as_scalar(eta),
                source_l1=as_scalar(sum(map(abs,J),F(0))),bound=as_scalar(bound)))
    # Full equation residual certifies the adjoint; it is not the projected
    # residual inside a numerical subspace.
    zhat=np.array([F(round(x*1000),1000) for x in z],dtype=object)
    res=A.T@zhat-ell
    dz=G.T@res
    assert np.array_equal(dz,zhat-z)
    assert ell@y-zhat@s==-res@y
    residual_bound=sum(map(abs,J),F(0))*maxabs(R@dz)
    actual=abs(ell@y-(R@zhat)@J)
    assert actual<=residual_bound
    # An independent original homogeneous solution detects missing initial data.
    nt=14;dt=F(1,4);K=np.array([[F(3,2),ep*F(2,3)],[ep*F(2,3),F(5,4)]],dtype=object)
    hom=zero(nt,2);hom[0]=[F(1,3),F(-1,5)];hom[1]=[F(2,7),F(1,11)]
    for t in range(1,nt-1):hom[t+1]=2*hom[t]-hom[t-1]-dt**2*K@hom[t]
    hobs=ell@hom[2:].reshape(-1);assert hobs!=0
    return dict(round=899,date='2026-10-06',formal_reports=899,cumulative_numbered_groups=3684,
        fresh_numbered_groups=1,all_checks_passed=True,
        argument_scope='Conditional finite geometry/source error transfer in the current860/869/872 causal branch; exact rational diagnostic only, not the original curved source or quantum gravity solver.',
        original860_mixed_transfer_reproduced=as_scalar(original),
        complete872_source_coefficients=[str(x) for x in parts],complete_source=str(total),
        fixed_background_parameter=str(ep),
        field_readout_parameter_jets=[as_scalar(x) for x in direct],
        four_contributions=[as_scalar(x) for x in individual],
        omitted_source_prediction=[as_scalar(x) for x in omissions],
        frozen_adjoint_first_derivative_defect=as_scalar(abs(direct[1])),
        frozen_adjoint_second_derivative_defect=as_scalar(abs(direct[2])),
        dual_test_budget_rows=rows,
        full_adjoint_residual_error=as_scalar(actual),full_adjoint_residual_bound=as_scalar(residual_bound),
        missing_initial_data_defect=as_scalar(abs(hobs)),
        exact_rational_forward_adjoint_and_jet_checks=True,
        source_ward_not_tested_by_oscillator=True,
        original_curved_source_and_adjoint_numerically_evaluated=False,
        actual_Qeff_physical_parameter_budget_completed=False,
        full_local_constraint_or_nonlinear_feedback_from_finite_menu=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else: assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k!='dual_test_budget_rows'},ensure_ascii=False,indent=2))
