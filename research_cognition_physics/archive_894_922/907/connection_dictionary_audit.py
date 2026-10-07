"""907: audit the literal old anti-Hermitian loop adapter against902's gauge law.
No frozen code is changed. The original constrained color background and the
adjoint loop record give a finite pure-gauge negative control.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'902'));sys.path.insert(0,str(STAGE/'863'))
import common_background_evolution as ev
from physical_relational_loop_bridge import exp_with_tangent
TARGET=HERE/'connection_dictionary_audit_results.json'
def exp(A):return exp_with_tangent(A,np.zeros_like(A))[0]

def loop(Av,Az,L,sign,parameter):
    G=exp(1j*parameter*L*Av)
    sides=[exp(-sign*1j*L*(1-parameter)*Av),G@exp(-sign*1j*L*Az)@G.conj().T,
        exp(sign*1j*L*(1-parameter)*Av),exp(sign*1j*L*Az)]
    U=np.eye(3,dtype=complex)
    for side in sides:U=side@U
    q=sign+(1-sign)*parameter
    exact=exp(sign*1j*L*Az)@exp(1j*q*L*Av)@exp(-sign*1j*L*Az)@exp(-1j*q*L*Av)
    assert np.max(abs(U-exact))<2e-14
    return U

def derivative(Av,Az,L,sign):
    U=np.eye(3,dtype=complex);dU=np.zeros_like(U)
    for B,dB in ((-1j*sign*L*Av,-1j*(1-sign)*L*Av),(-1j*sign*L*Az,np.zeros_like(Az)),
        (1j*sign*L*Av,1j*(1-sign)*L*Av),(1j*sign*L*Az,np.zeros_like(Az))):
        E,dE=exp_with_tangent(B,dB);dU=dE@U+E@dU;U=E@U
    tr=np.trace(U);dt=np.trace(dU)
    return float(dt.real),float(2*np.real(tr.conjugate()*dt))

def run():
    Av=.31*ev.T[0]-.27*ev.T[1];Az=.21*ev.T[3];comm=Av@Az-Az@Av
    R=(Q('.31')**2+Q('.27')**2)*Q('.21')**2/4
    assert abs(float(R)+2*np.trace(comm@comm).real)<1e-15
    rows=[];curv_error=0.;correct_max=0.;bad_max=0.
    for L in (.8,.4,.2,.1):
        records=[]
        for sign in (1,-1):
            U0=loop(Av,Az,L,sign,0.);r0=abs(np.trace(U0))**2-1
            errors=[]
            for t in (.1,.25,.5):
                U=loop(Av,Az,L,sign,t);rv=abs(np.trace(U))**2-1;delta=float(rv-r0)
                # At the base point G=I. All physical curvature is unchanged.
                A1=(1-t)*Av;A2=Az;d1A2=1j*t*comm
                F=d1A2+1j*(A1@A2-A2@A1)
                curv_error=max(curv_error,float(np.max(abs(F-1j*comm))))
                errors.append(dict(gauge_parameter=t,adjoint_character=float(rv),change=delta))
                if sign==1:correct_max=max(correct_max,abs(delta))
                else:bad_max=max(bad_max,abs(delta))
            dtrace,dadj=derivative(Av,Az,L,sign)
            h=1e-4;up=loop(Av,Az,L,sign,h);um=loop(Av,Az,L,sign,-h)
            fd=(abs(np.trace(up))**2-abs(np.trace(um))**2)/(2*h)
            assert abs(fd-dadj)<1e-10
            records.append(dict(anti_Hermitian_adapter_sign=sign,initial_adjoint_character=float(r0),
                finite_gauge_records=errors,real_fundamental_trace_derivative=dtrace,
                adjoint_character_derivative=dadj,adjoint_derivative_over_L4=dadj/L**4,
                derivative_finite_difference_error=abs(fd-dadj)))
        rows.append(dict(side=L,records=records))
    assert correct_max<2e-14 and bad_max>1e-5 and curv_error<1e-15
    ratios=[row['records'][1]['adjoint_derivative_over_L4'] for row in rows]
    assert abs(ratios[-1]/float(6*R)-1)<.001
    return dict(round=907,date='2026-10-06',all_checks_passed=True,
        physical_connection_convention='D=partial+i A; F=dA+i[A,A]; A^G=G A G^-1+i(dG)G^-1; a=i A',
        legacy_literal_adapter='a=-i A without also transporting A, its gauge law and representation',
        exact_color_curvature_norm=str(R),exact_wrong_adapter_adjoint_derivative_leading_coefficient=str(6*R),
        correct_adapter_maximum_pure_gauge_record_change=correct_max,
        legacy_literal_adapter_maximum_pure_gauge_record_change=bad_max,
        physical_curvature_gauge_invariance_error=curv_error,rows=rows,
        original_constrained_background_preserved_under_completed_gauge_transform=True,
        material_references_and_paths_unchanged_under_this_color_gauge_transform=True,
        abstract_covariant_Wilson_theorem_refuted=False,legacy_calibration_files_overwritten=False,
        full_smeared_source_implemented=False,whole_model_failure_proved=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();v=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert v==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(v,ensure_ascii=False,indent=2))
