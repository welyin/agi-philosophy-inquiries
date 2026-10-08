"""1032: algebraic electric completion, not physical endpoint preparation.

Integer word construction is checked separately against center characters and
tensor highest-weight matrices. No finite grid proves the universal theorem.
Default is read-only; --write creates the scientific result once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / 'electric_completion_selection_results.json'
PS = (1, 2, 3, 6)
HISTORY = (
    'archive_531_553/research_note_531.md',
    'archive_1009_/research_note_1021.md',
    'archive_1009_/1031/next_selection_audit.md',
    'archive_1009_/1009/input_dependency_ledger_v0_1.md',
    'archive_956_989/957/drafts/unified_operation_hypotheses_v0_2.md',
    'archive_990_1008/1008/overall_operation_hypothesis_v2_2.md',
)
SM = {'q': (1, 0, 1, 1), 'uc': (0, 1, 0, -4),
      'dc': (0, 1, 0, 2), 'l': (0, 0, 1, -3),
      'ec': (0, 0, 0, 6), 'H': (0, 0, 1, 3)}


def center_weight(rep):
    a, b, n, q = rep
    return 2*(a+2*b)+3*n+q


def allowed(rep, p):
    return center_weight(rep) % p == 0


def word(rep, singlet_charge=6):
    a, b, n, q = rep
    q0 = -2*a+2*b+3*n
    if (q-q0) % singlet_charge:
        return None
    k = (q-q0)//singlet_charge
    return dict(D=a, Dbar=b, H=n, positive_singlet=max(k, 0),
                negative_singlet=max(-k, 0), singlet_charge=singlet_charge,
                length=a+b+n+abs(k), q0=q0)


def rebuild(w):
    return (w['D'], w['Dbar'], w['H'],
            -2*w['D']+2*w['Dbar']+3*w['H']
            +w['singlet_charge']*(w['positive_singlet']-w['negative_singlet']))


def center_certificate():
    kernel = [k for k in range(6) if all(k*center_weight(r) % 6 == 0 for r in SM.values())]
    assert kernel == list(range(6))
    cases = []
    for p in PS:
        gamma = list(range(0, 6, 6//p))
        witness = (0, 0, 0, p)
        assert allowed(witness, p)
        completion_kernel = [k for k in kernel if k*p % 6 == 0]
        assert completion_kernel == gamma
        angle = Fraction(p, 6)
        cases.append(dict(p=p, quotient_subgroup_exponents=gamma,
            residual_electric_center_order=6//p,
            pure_singlet_witness=list(witness), sm_word_exists=word(witness) is not None,
            witness_z_phase_turns=str(angle),
            completed_kernel_exponents=completion_kernel,
            new_left_Weyl_charges=[p, -p], anomaly_cubic=p**3+(-p)**3,
            anomaly_gravitational=p-p, invariant_mass_charge=p-p))
    return dict(existing_SM_center_weights={k:center_weight(v) for k,v in SM.items()},
                common_kernel_exponents=kernel, cases=cases)


def finite_word_grid():
    rows = []
    for p in PS:
        allowed_count = original_count = completed_count = 0
        longest = 0
        for a in range(5):
            for b in range(5):
                for n in range(5):
                    for q in range(-24, 25):
                        r = (a,b,n,q)
                        if not allowed(r,p):
                            continue
                        allowed_count += 1
                        old = word(r)
                        assert (old is not None) == allowed(r,6)
                        if old is not None:
                            assert rebuild(old) == r
                            original_count += 1
                        new = word(r,p)
                        assert new is not None and rebuild(new) == r
                        completed_count += 1
                        longest = max(longest,new['length'])
        assert completed_count == allowed_count
        assert (original_count == allowed_count) == (p == 6)
        rows.append(dict(p=p, allowed=allowed_count, SM_generated=original_count,
                         completed_generated=completed_count, maximum_constructed_length=longest))
    return dict(grid={'a':[0,4],'b':[0,4],'n':[0,4],'q':[-24,24]}, cases=rows)


def coproduct(matrices):
    if not matrices:
        return np.zeros((1,1),complex)
    dims=[m.shape[0] for m in matrices]
    out=np.zeros((math.prod(dims),)*2,complex)
    for j,mat in enumerate(matrices):
        term=np.ones((1,1),complex)
        for k,d in enumerate(dims):
            term=np.kron(term,mat if k==j else np.eye(d))
        out+=term
    return out


def vector_product(vectors):
    out=np.ones(1,complex)
    for v in vectors:
        out=np.kron(out,v)
    return out


def highest_weight_certificate():
    e12=np.zeros((3,3));e12[0,1]=1
    e23=np.zeros((3,3));e23[1,2]=1
    h1=np.diag([1,-1,0]);h2=np.diag([0,1,-1])
    e=np.eye(3)
    rows=[]
    for a in range(4):
        for b in range(4):
            if a+b>4:
                continue
            v=vector_product([e[0]]*a+[e[2]]*b)
            matrices=[coproduct([m]*a+[-m.T]*b) for m in (e12,e23,h1,h2)]
            residual=max(np.linalg.norm(matrices[0]@v),np.linalg.norm(matrices[1]@v),
                         np.linalg.norm(matrices[2]@v-a*v),np.linalg.norm(matrices[3]@v-b*v))
            assert residual<1e-12 and abs(np.vdot(v,v)-1)<1e-12
            rows.append(dict(a=a,b=b,tensor_dimension=len(v),residual=float(residual)))
    weak=[]
    raising=np.array([[0,1],[0,0]])
    h=np.diag([1,-1])
    for n in range(7):
        v=vector_product([np.array([1,0])]*n)
        residual=max(np.linalg.norm(coproduct([raising]*n)@v),
                     np.linalg.norm(coproduct([h]*n)@v-n*v))
        assert residual<1e-12
        weak.append(dict(n=n,tensor_dimension=len(v),residual=float(residual)))
    return dict(color_cases=rows,weak_cases=weak,
                tensor_factors_are_distinguishable_slots=True,
                not_a_bound_state_or_local_same_mode_fermion_construction=True)


def signed_closure_certificate():
    # Independent finite abelian character-group closure, not the word algorithm.
    modulus=12
    def plus(x,y):return ((x[0]+y[0])%3,(x[1]+y[1])%2,(x[2]+y[2])%modulus)
    generators=[(1,0,-2),(0,1,3),(0,0,6)]
    rows=[]
    for p in PS:
        for extended in (False,True):
            gs=generators+([(0,0,p)] if extended else [])
            gs=gs+[(-x,-y,-z) for x,y,z in gs]
            reached={(0,0,0)};todo=[(0,0,0)]
            while todo:
                x=todo.pop()
                for g in gs:
                    z=plus(x,g)
                    if z not in reached:
                        reached.add(z);todo.append(z)
            target={(t,n,q) for t in range(3) for n in range(2) for q in range(modulus)
                    if (2*t+3*n+q)%p==0}
            expected=target if extended else {x for x in target if (2*x[0]+3*x[1]+x[2])%6==0}
            assert reached==expected
            rows.append(dict(p=p,extended=extended,reached=len(reached),allowed=len(target),
                             generated_all=reached==target))
    return dict(charge_modulus=modulus,cases=rows,
                finite_character_quotient_only_not_all_irreps=True)


def run():
    return dict(round=1032,status='local_calculations_verified_independent_agent_review_incomplete',
        new_calibration_groups=1,cumulative_research_groups=3809,
        new_adopted_cognitive_axioms=0,proposed_selection_condition='E_alg',
        complete_existing_endpoint_menu_is_additional_input=True,
        external_independent_agent_review_completed=False,goal_completed=False,
        center=center_certificate(),words=finite_word_grid(),
        highest_weights=highest_weight_certificate(),closure=signed_closure_certificate(),
        physical_endpoints_prepared=False,resource_uniform_bound_proved=False,
        magnetic_completeness_proved=False,full_quantum_gravity_model_proved=False,
        no_actual_low_energy_matching_error_certified=True,
        code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        historical_source_sha256={s:hashlib.sha256((BASE/s).read_bytes()).hexdigest() for s in HISTORY})


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12)
    else:assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        with OUT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps({'round':1032,'verified':True,'word_grid':result['words']['cases'],
                      'color_tensor_cases':len(result['highest_weights']['color_cases']),
                      'weak_tensor_cases':len(result['highest_weights']['weak_cases']),
                      'independent_agent_review_completed':False},ensure_ascii=False))
