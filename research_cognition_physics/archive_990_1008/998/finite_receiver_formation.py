"""Finite internal receiver screen in an explicitly classical thermodynamic model.

Reuses virial/heat-capacity mechanisms; does not derive gravity or build memory.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'finite_receiver_formation_results.json'


def read(p):
    return json.loads(p.read_text('utf-8-sig'))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def n(x):
    return dict(exact=str(x),value=float(x))


def log_interval(x,terms=10):
    """Exact positive artanh series with a geometric upper bound on its tail."""
    assert x>1
    z=(x-1)/(x+1)
    lo=2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
    remainder=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return lo,lo+remainder


def bisect_increasing(fun,lo,hi,steps=52):
    assert fun(lo)<0<fun(hi)
    for _ in range(steps):
        mid=(lo+hi)/2
        if fun(mid)<0:lo=mid
        else:hi=mid
    assert fun(lo)<=0<=fun(hi)
    return lo,hi


def run():
    c,A=F(3),F(3,4)
    Tg0,Tr0,Tf=F(497,1024),F(1,4),F(1,2)
    E0=-c*Tg0+A*Tr0**4
    Ef=-c*Tf+A*Tf**4
    Q=c*(Tf-Tg0)
    assert E0==Ef==F(-93,64)
    assert Q==A*(Tf**4-Tr0**4)==F(45,1024)>0
    # Virial sequence: K=c*T, U=-2c*T, E=K+U.
    dK=c*(Tf-Tg0);dU=-2*dK;dEg=dK+dU
    assert dK==Q and dU==-2*Q and dEg+Q==0
    Cr=4*A*Tf**3
    drift_per_conductance=1/c-1/Cr
    entropy_second_derivative=drift_per_conductance/Tf**2
    assert Cr==F(3,8)<c
    assert drift_per_conductance==F(-7,3)
    assert entropy_second_derivative==F(-28,3)
    # On this entire interval radiation heats faster per unit received energy.
    assert 4*A*Tr0**3>0 and 4*A*Tf**3<c
    loglo,loghi=log_interval(Tf/Tg0)
    dSr=F(4,3)*A*(Tf**3-Tr0**3)
    dSlo=dSr-c*loghi;dShi=dSr-c*loglo
    assert dSr==F(7,64) and 0<dSlo<dShi
    assert dShi-dSlo<F(1,10**35)
    # Same total energy also has a higher-temperature, thermally unstable root.
    Tc=F(1);Emin=-F(3,4)*c*Tc
    assert 4*A*Tc**3==c and Emin==F(-9,4)<E0<0
    polynomial=lambda t:A*t**4-c*t-E0
    hi_lo,hi_hi=bisect_increasing(polynomial,Tc,F(2))
    assert hi_lo>1 and 4*A*hi_lo**3>c
    high_lambda_lo=1/c-1/(4*A*hi_lo**3)
    high_lambda_hi=1/c-1/(4*A*hi_hi**3)
    assert 0<high_lambda_lo<high_lambda_hi
    # Independent polynomial roots and finite energy-transfer curve diagnostics.
    roots=np.roots([float(A),0,0,-float(c),-float(E0)])
    positive=sorted(float(z.real) for z in roots if abs(z.imag)<1e-10 and z.real>0)
    assert len(positive)==2 and abs(positive[0]-float(Tf))<1e-12
    assert abs(positive[1]-float((hi_lo+hi_hi)/2))<1e-12
    fractions=[F(k,16) for k in range(17)]
    curve=[]
    for f in fractions:
        q=Q*f
        tg=Tg0+q/c
        tr=(float(Tr0**4+q/A))**.25
        drift=float(tg)-tr
        st=-float(c)*math.log(float(tg/Tg0))+float(F(4,3)*A)*(tr**3-float(Tr0**3))
        energy=-float(c*tg)+float(A)*tr**4
        assert abs(energy-float(E0))<1e-14
        assert drift>-1e-14
        curve.append(dict(Q_fraction=str(f),T_gravity=float(tg),T_radiation=tr,
                          temperature_difference=drift,entropy_change=st))
    assert all(curve[i+1]['entropy_change']>curve[i]['entropy_change'] for i in range(16))
    floating_entropy=dSr-c*F(str(math.log(float(Tf/Tg0))))
    assert abs(float(floating_entropy)-float(dSlo))<1e-14
    # A cold infinite reservoir is a distinct equation, with unstable response.
    fixed_bath_drift_per_conductance=1/c
    # A lower total energy has no equal-temperature solution on this assumed branch.
    noeq_Tg,noeq_Tr=F(1),F(1,2)
    noeq_E=-c*noeq_Tg+A*noeq_Tr**4
    assert noeq_E==F(-189,64)<Emin
    sources=[HERE/'drafts/formation_selection.md',HERE/'drafts/NEXT.md',
        HERE/'drafts/boundary_adoption_audit_checks.json',
        STAGE/'research_note_975.md',STAGE/'research_note_992.md',
        STAGE.parent/'archive_301_341/research_note_333.md',Path(__file__)]
    return dict(round=998,date='2026-10-07',all_scientific_checks_passed=True,
        kind='classical_finite_thermodynamic_mechanism_screen',
        inputs=dict(c_gravity=n(c),radiation_coefficient=n(A),
                    initial_T_gravity=n(Tg0),initial_T_radiation=n(Tr0),
                    normalization='T/T_unit; E/(C_unit*T_unit); S/C_unit'),
        stable_case=dict(total_energy=n(E0),final_T=n(Tf),released_energy=n(Q),
            kinetic_change=n(dK),potential_change=n(dU),source_energy_change=n(dEg),
            receiver_energy_change=n(Q),radius_ratio=n(Tg0/Tf),
            radiation_heat_response=n(Cr),linear_drift_over_conductance=n(drift_per_conductance),
            entropy_curvature_wrt_transferred_energy=n(entropy_second_derivative),
            receiver_entropy_change=n(dSr),
            total_entropy_change_lower=n(dSlo),total_entropy_change_upper=n(dShi),
            curve=curve),
        same_energy_unstable_case=dict(T_lower=n(hi_lo),T_upper=n(hi_hi),
            linear_drift_over_conductance_lower=n(high_lambda_lo),
            linear_drift_over_conductance_upper=n(high_lambda_hi)),
        domain_boundary=dict(critical_T=n(Tc),minimum_equilibrium_energy=n(Emin),
            illustrative_no_equilibrium_energy=n(noeq_E),
            infinite_fixed_bath_linear_drift_over_conductance=n(fixed_bath_drift_per_conductance)),
        thermal_stability_only=True,mechanical_stability_certified=False,
        unbound_to_bound_formation_simulated=False,actual_heat_transfer_rate_calculated=False,
        actual_work_extracted=False,record_device_constructed=False,
        cosmological_history_generated=False,quantum_or_gravity_derived=False,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sources})


def compare(actual,saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for k in actual:compare(actual[k],saved[k])
    elif isinstance(actual,list):
        assert len(actual)==len(saved)
        for x,y in zip(actual,saved):compare(x,y)
    elif isinstance(actual,float):assert math.isclose(actual,saved,abs_tol=1e-14,rel_tol=1e-13)
    else:assert actual==saved,(actual,saved)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps(dict(round=998,all_scientific_checks_passed=True,
        Q=out['stable_case']['released_energy']['value'],
        entropy_increase=out['stable_case']['total_entropy_change_lower']['value'],
        second_root=out['same_energy_unstable_case']['T_lower']['value'],
        full_goal_completed=False),ensure_ascii=False,indent=2))
