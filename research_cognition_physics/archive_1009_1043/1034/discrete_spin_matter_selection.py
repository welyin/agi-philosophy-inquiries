"""1034: audit an ADOPTED Spin-Z4 anomaly theorem and its mass compatibility.

The bordism theorem is a cited input, not numerically proved. Exact arithmetic
checks its specialization; matrix checks concern a stated tree mass sector only.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
OUT=HERE/'discrete_spin_matter_selection_results.json'
HISTORY=(
    'archive_585_628/research_note_616.md',
    'archive_585_628/research_note_628.md',
    'archive_531_553/research_note_539.md',
    'archive_956_989/981/drafts/common_parent_contract_v1.md',
    'archive_1009_/research_note_1010.md',
    'archive_1009_/research_note_1012.md',
    'archive_1009_/1009/input_dependency_ledger_v0_1.md',
    'archive_1009_/1033/NEXT.md',
)
# Left-handed Weyl fields: multiplicity, q=6Y, b=3(B-L), color index,
# weak index. Dynkin indices here omit the common factor 1/2.
FERMIONS={
    'q':(6,1,1,2,3), 'uc':(3,-4,-1,1,0),
    'dc':(3,2,-1,1,0), 'l':(2,-3,-3,0,1),
    'ec':(1,6,3,0,0), 'N':(1,0,3,0,0),
}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def xcharge(q,b):
    x=F(5*b-2*q,3)
    assert x.denominator==1
    return int(x)


def eta_twisted(charges):
    s1=sum(charges);s3=sum(q**3 for q in charges)
    return F(11*s3-5*s1,96)%1,F(2*s3+s1,4)%1


def eta_untwisted(charges):
    # Hsieh (2.45), n=4, ordinary Spin x Z4; a DIFFERENT background class.
    s1=sum(charges);s3=sum(q**3 for q in charges)
    return F(30*s3,24)%1,F(2*s1,4)%1


def charge_table():
    rows=[]
    for name,(d,q,b,t3,t2) in FERMIONS.items():
        x=xcharge(q,b)
        assert x%4==1
        rows.append(dict(field=name,multiplicity=d,q_6Y=q,b_3BL=b,X=x,
                         Z4=x%4,twice_SU3_index=t3,twice_SU2_index=t2))
    return rows


def anomaly_checks():
    periodic=[]
    for q in range(-65,66,2):
        a,b=eta_twisted([q]);sign=1 if q%4==1 else -1
        assert a==F(sign,16)%1 and b==F(3*sign,4)%1
        assert eta_twisted([q+4])==(a,b)
        periodic.append(dict(charge=q,eta_X=str(a),eta_Y=str(b),sign=sign))
    sm=[xcharge(q,b) for name,(d,q,b,_,__) in FERMIONS.items()
        if name!='N' for _ in range(d)]
    cases=[]
    for ng in range(1,7):
        for plus in range(0,21):
            for minus in range(0,21):
                # ±1 are representatives of the two odd Z4 charge classes.
                a,b=eta_twisted(sm*ng+[1]*plus+[-1]*minus)
                nu=(15*ng+plus-minus)%16
                assert a==F(nu,16) and b==F(3*nu,4)%1
                assert (a==b==0)==(nu==0)
                cases.append(dict(generations=ng,plus=plus,minus=minus,nu=nu,
                                  audited_obstruction_cancels=nu==0))
    matches=[r for r in cases if r['generations']==3 and r['nu']==0]
    minimum=min(r['plus']+r['minus'] for r in matches)
    minimal=[(r['plus'],r['minus']) for r in matches if r['plus']+r['minus']==minimum]
    assert minimum==3 and minimal==[(3,0)]
    positive_only=[n for n in range(65) if eta_twisted(sm*3+[1]*n)==(0,0)]
    assert positive_only==[3,19,35,51]
    assert eta_untwisted([1]*4)==(0,0)
    assert eta_twisted([1]*4)==(F(1,4),F(0))
    # Vectorlike invariant masses do not change the audited obstruction.
    pair=[eta_twisted([q,-q]) for q in range(-15,16,2)]
    assert all(v==(0,0) for v in pair)
    return dict(odd_charge_periodicity=periodic,menu_cases=cases,
                three_generation_minimum=minimum,minimal_menus=minimal,
                positive_only_solutions_to_64=positive_only,
                ordinary_spin_four_fermions=[str(v) for v in eta_untwisted([1]*4)],
                twisted_spin_four_fermions=[str(v) for v in eta_twisted([1]*4)],
                vectorlike_pairs_checked=len(pair))


def continuous_lift_checks():
    rows=[]
    for ng in (1,2,3,4):
        for count in range(0,7):
            tr={k:F(0) for k in ('X','X3','Y2X','YX2','SU3_2X','SU2_2X')}
            for name,(d,q,b,t3,t2) in FERMIONS.items():
                copies=count if name=='N' else ng
                x=xcharge(q,b);y=F(q,6)
                tr['X']+=copies*d*x;tr['X3']+=copies*d*x**3
                tr['Y2X']+=copies*d*y*y*x;tr['YX2']+=copies*d*y*x*x
                tr['SU3_2X']+=copies*t3*x;tr['SU2_2X']+=copies*t2*x
            assert tr['X']==5*(count-ng) and tr['X3']==125*(count-ng)
            assert all(tr[k]==0 for k in ('Y2X','YX2','SU3_2X','SU2_2X'))
            rows.append(dict(generations=ng,N_with_X5=count,
                             anomaly_traces={k:str(v) for k,v in tr.items()}))
    return rows


def operator_checks():
    x={r['field']:r['X'] for r in charge_table()};x.update(H=-2,Hbar=2,S=2,Nminus=-1)
    words={
        'up_Yukawa':['q','H','uc'], 'down_Yukawa':['q','Hbar','dc'],
        'charged_lepton_Yukawa':['l','Hbar','ec'],
        'neutrino_Dirac':['l','H','N'], 'bare_Majorana':['N','N'],
        'Weinberg':['l','H','l','H'], 'scalar_Majorana':['S','N','N'],
        'vectorlike_singlet_mass':['N','Nminus'],
        'Higgs_norm':['Hbar','H'],
    }
    rows=[]
    for name,word in words.items():
        q=sum(x[f] for f in word)%4
        rows.append(dict(operator=name,fields=word,Z4=q,invariant=q==0))
    assert {r['operator'] for r in rows if not r['invariant']}=={'bare_Majorana','Weinberg'}
    assert all((2+4*k)%4==2 for k in range(33))
    # Composite LH has charge -i; sterile N has charge +i.
    U=np.diag([-1j]*3+[1j]*3)
    rng=np.random.default_rng(1034)
    checks=[]
    for n in range(6):
        y=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/10
        r=rng.normal(size=(3,3))/10
        mass=r+r.T+np.diag([2.,3.,4.])
        zero=np.zeros((3,3))
        k=np.block([[zero,y],[y.T,mass]])
        k_opposite=np.block([[zero,y],[y.T,-mass]])
        cov=float(np.max(np.abs(U.T@k_opposite@U-k)))
        broken=float(np.max(np.abs(U.T@k@U-k)))
        c=-y@np.linalg.solve(mass,y.T)
        elimination=np.vstack([np.eye(3),-np.linalg.solve(mass,y.T)])
        schur=float(np.max(np.abs(elimination.T@k@elimination-c)))
        assert cov<1e-13 and schur<1e-13 and broken>1
        assert np.max(np.abs(-y@np.linalg.solve(-mass,y.T)+c))<1e-13
        checks.append(dict(seed_case=n,parent_covariance_residual=cov,
                           frozen_mass_symmetry_defect=broken,Schur_residual=schur))
    return dict(words=rows,Higgs_norm_insertions_checked=33,matrix_cases=checks)


def surviving_freedom():
    # Exact rational tree-level matching family, not numerical pole masses.
    rows=[]
    for t in (F(1,2),F(1),F(2),F(3)):
        y=[t*F(1,10),t*F(1,5),t*F(3,10)]
        mass=[t*t*F(2),t*t*F(3),t*t*F(5)]
        c=[-yi*yi/mi for yi,mi in zip(y,mass)]
        c6=[yi*yi/(mi*mi) for yi,mi in zip(y,mass)]
        assert c==[-F(1,200),-F(1,75),-F(9,500)]
        rows.append(dict(scale=str(t),Y_diagonal=[str(v) for v in y],
                         M_diagonal=[str(v) for v in mass],
                         C5_diagonal=[str(v) for v in c],
                         C6_diagonal=[str(v) for v in c6],
                         C6_trace=str(sum(c6))))
    assert len({r['C6_trace'] for r in rows})==4
    return rows


def run():
    return dict(round=1034,date='2026-10-08',new_calibration_groups=1,
                cumulative_research_groups=3811,new_adopted_cognitive_axioms=0,
                goal_completed=False,external_independent_agent_review_completed=False,
                bordism_classification_adopted_from_literature=True,
                full_mixed_global_anomaly_classification_proved_here=False,
                existing_P981_has_new_exact_Z4=False,
                neutrino_masses_or_generation_number_predicted=False,
                physical_implementation_or_error_certified=False,
                code_sha256=sha(Path(__file__)),
                historical_source_sha256={p:sha(BASE/p) for p in HISTORY},
                charge_table=charge_table(),anomalies=anomaly_checks(),
                continuous_lift=continuous_lift_checks(),
                mass_operators=operator_checks(),remaining_matching_freedom=surviving_freedom())


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,abs_tol=1e-12,rel_tol=1e-9),(a,b)
    else:assert a==b,(a,b)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true')
    args=p.parse_args();r=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(r,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:compare(r,json.loads(OUT.read_text('utf8')))
    print(json.dumps(dict(round=1034,menu_cases=len(r['anomalies']['menu_cases']),
        minimum=r['anomalies']['minimal_menus'],
        positive_only=r['anomalies']['positive_only_solutions_to_64'],
        matrix_cases=len(r['mass_operators']['matrix_cases']),
        surviving_C6_traces=[v['C6_trace'] for v in r['remaining_matching_freedom']]),ensure_ascii=False))
