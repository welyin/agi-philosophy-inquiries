"""Bounded initial causal screening; no future integration or quantum model.
Reuses and binds the saved full 5^3 Hamiltonian. Default read-only comparison
uses explicit tolerances for floating diagnostics; --write creates results once.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import math
import numpy as np

PARENT_HASH='e9bbcda917d9b9cebcffbc7c474320c418f952e50fd2e48034c3da3f97a4a759'
DIAGNOSTIC_TOL=1e-8
COMPARE_ABS=5e-10
COMPARE_REL=1e-9


def load_parent():
    p=Path(__file__).parent.parent/'finite_grid_candidate_after1256'/'diagnostic.py'
    actual=sha256(p.read_bytes()).hexdigest()
    assert actual==PARENT_HASH, 'Saved parent Hamiltonian changed; review is required'
    spec=importlib.util.spec_from_file_location('saved_grid_parent',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m


def fs(x): return f'{x.numerator}/{x.denominator}'


def exact_initial_algebra():
    delta=Q(2);v=delta**3;alpha=Q(1,10)
    A0=Q(187,200);U=Q(167,200)
    A=A0-alpha**2/2;J=A-U
    assert J==Q(19,200)>0
    assert Q(7,40)+alpha**2/2+A==Q(111,100)
    R=A;D0=Q(0)
    assert R>0 and R*R+D0*D0==A*A and -A+R==0
    hxx=-v/(4*delta**2)
    force=-hxx/v
    t3=force/6
    assert hxx==Q(-1,2) and force==Q(1,16) and t3==Q(1,96)
    return dict(
      status='PASS_EXACT_INITIAL_ALGEBRA',
      alpha=fs(alpha),J_source_cell=fs(J),J_other_cells='1/10',
      J_change_source=fs(-alpha**2/2),source_Tnn='111/100',
      source_extra_spatial_pressure=fs(alpha**2/2),
      momentum_constraint='0',hamiltonian_constraint='0: positive R=A, Ec=A',
      cell_weight=fs(v),H_X_source_X_remote=fs(hxx),H_PX_source_PX_source=fs(v),
      remote_PX_force_coefficient=fs(force),remote_X_t_cubed_coefficient=fs(t3),
      physical_input_scope='Finite local PX preparation also changes the same clock J; no new free source.',
      meaning='Exact initial source and explicitly derived jet coefficients, not a finite-time remainder bound.')


def full_tangent_diagnostic(m):
    y=m.initial(0.0);source=(0,0,0);remote=(2,0,0);eps=1e-5
    impulse=np.zeros(m.GRID);impulse[source]=1
    lap=sum(m.diff(m.diff(impulse,i),i) for i in range(3))-impulse
    rows={}
    for slot,name in [(15,'A_e_PX'),(14,'A_e_X')]:
        plus=y.copy();minus=y.copy()
        plus[(slot,)+source]+=eps;minus[(slot,)+source]-=eps
        direction=(m.flow_complex_step(plus)-m.flow_complex_step(minus))/(2*eps)
        expected=np.zeros_like(y)
        if slot==15:
            expected[(14,)+source]=1
        else:
            expected[15]=lap
            for a in (6,7,8):expected[(a,)+source]=-1/20
        err=float(np.max(np.abs(direction-expected)))
        assert err<DIAGNOSTIC_TOL
        rows[name]=dict(
          full_vector_max_error=err,
          X_source=float(direction[(14,)+source]),X_remote=float(direction[(14,)+remote]),
          PX_remote=float(direction[(15,)+remote]),
          metric_h_max=float(np.max(np.abs(direction[:6]))),
          geometry_pi_at_source=[float(direction[(a,)+source]) for a in (6,7,8)],
          geometry_pi_at_remote=[float(direction[(a,)+remote]) for a in (6,7,8)],
          first_clock_internal_max=float(np.max(np.abs(direction[22:28]))),
          second_flow_max=float(np.max(np.abs(direction[28:36]))))
    assert abs(rows['A_e_X']['PX_remote']-1/16)<DIAGNOSTIC_TOL
    return dict(
      status='PASS_FLOATING_DIAGNOSTIC_NOT_A_PROOF_OF_FUTURE_CC',
      state_dimension=4500,source=list(source),remote=list(remote),
      outer_central_difference_step=eps,inner_complex_step='1e-28',
      absolute_assertion_tolerance=DIAGNOSTIC_TOL,
      method='Four full-H vector fields, each with all 4500 canonical partial derivatives; no reduced-wave substitution.',
      columns=rows,
      A_cubed_e_PX_remote_X_via_exact_initial_X_row=rows['A_e_X']['PX_remote'],
      note='At t=0 A_e_PX=e_X. A^2_e_PX is remote in PX, not X. The remote X entry occurs in A^3.')


def local_hessian_diagnostic(m):
    y=m.initial(0.0);source=(0,0,0);remote=(2,0,0);a=.1;base=m.mother(y)
    def h(x,z):
        q=y.copy();q[(14,)+source]+=x;q[(14,)+remote]+=z
        return float(m.mother(q)['H'])
    hxx=(h(a,a)-h(a,-a)-h(-a,a)+h(-a,-a))/(4*a*a)
    plus=y.copy();minus=y.copy();plus[(15,)+source]+=a;minus[(15,)+source]-=a
    hpp=(m.mother(plus)['H']-2*base['H']+m.mother(minus)['H'])/(a*a)
    jchange=m.mother(plus)['J'][source]-base['J'][source]
    assert abs(hxx+.5)<DIAGNOSTIC_TOL and abs(hpp-8)<DIAGNOSTIC_TOL
    assert abs(jchange+a*a/2)<DIAGNOSTIC_TOL
    return dict(H_XX_remote=float(hxx),H_PX_PX_source=float(hpp),
                finite_preparation_J_change=float(jchange),step=a,
                status='PASS_FLOATING_DIAGNOSTIC')


def calculate():
    m=load_parent()
    return dict(
      schema='finite_grid_causal_screen_after1256.v1',
      status='FIXED_GRID_CLASSICAL_RESPONSE_REMAINS_CC_TASK_NOT_CERTIFIED',
      check_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
      parent_diagnostic_sha256=PARENT_HASH,
      exact_algebra=exact_initial_algebra(),
      full_tangent=full_tangent_diagnostic(m),
      local_hessian=local_hessian_diagnostic(m),
      saved_float_comparison=dict(absolute_tolerance=COMPARE_ABS,relative_tolerance=COMPARE_REL),
      not_claimed=['actual local kick permission','actual early X readout permission',
                   'nonzero probability difference for the original K three-outcome task',
                   'fixed-time numerical remainder certificate','general HDA failure theorem',
                   'failure of all discretizations or of the research program'],
      decision='Do not continue by shrinking epsilon alone. A finite control/effect menu and an independent spatial/bandwidth causality budget are required. No research number.')


def compare(a,b,path='result'):
    assert type(a)==type(b), (path,type(a),type(b))
    if isinstance(a,dict):
        assert set(a)==set(b),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f'{path}[{i}]')
    elif isinstance(a,float):
        assert math.isfinite(a) and math.isfinite(b),path
        assert math.isclose(a,b,abs_tol=COMPARE_ABS,rel_tol=COMPARE_REL),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    target=Path(__file__).with_name('results.json')
    if args.write and target.exists():raise FileExistsError('Refusing to overwrite results.json')
    result=calculate()
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('CREATED_ONCE',target)
    else:
        compare(json.loads(target.read_text(encoding='utf-8')),result)
        print('READ_ONLY_RECOMPUTATION_MATCH_WITH_DECLARED_FLOAT_TOLERANCES')
    print(result['status'])
    for name,row in result['full_tangent']['columns'].items():
        print(name,'max_error',row['full_vector_max_error'],'remote_PX',row['PX_remote'])
    print('EXACT_REMOTE_X_T3',result['exact_algebra']['remote_X_t_cubed_coefficient'])

if __name__=='__main__':main()
