"""1271: verify the physical baseline and complete-type join, not old heat data.

Only old model constructors are imported. The new all-history assertion is
proved by typed structural induction in proof.md, not by these finite tests.
Default recomputes without writes; --write-once creates results.json exclusively.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
SOURCES=[
 'archive_1106_1119/1111/check.py',
 'archive_1106_1119/1112/check.py',
 'archive_1120_1158/1158/proof.md',
 'archive_1161_1200/1166/proof.md',
 'archive_1201_1218/_shared/thermal_roles_adoption_after1218.md',
 'archive_1219_1233/1219/candidate_proof.md',
 'archive_1219_1233/1219/results.json',
 'archive_1219_1233/1219/independent_review.md',
 'archive_1234_/1234/proof.md',
 'archive_1234_/1234/NEXT.md',
]
TOL=3e-9

def load(name,rel):
    spec=importlib.util.spec_from_file_location(name,BASE/rel)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod

def err(a):
    x=float(np.linalg.norm(a));assert x<TOL,x
    return x

def calculate():
    old=load('old1111_for1271','archive_1106_1119/1111/check.py')
    mat=load('old1112_for1271','archive_1106_1119/1112/check.py')
    m0=old.model();m=mat.Material(1)
    I=np.eye(288,dtype=complex)
    shape=(4,8,3,3,1,1,288)
    hplus=m.h_action(I.reshape(shape)).reshape(288,288)
    common={
      'same_FK_matrix':err(m.hfk-m0['hfk']),
      'same_Gamma':err(m.gamma-m0['gamma']),
      'full_physical_mass_shift':err(hplus-m0['hout']-2*I),
      'direct_full_generator':err(hplus-m.baseline*I-m.hfk),
    }
    bare=np.kron(np.eye(32),np.kron(np.diag([0,0,1]),np.eye(3)))
    C=m.gamma@bare@m.gamma.conj().T
    common['same_dressed_source_control']=err(hplus@C-C@hplus)
    phase_checks=[]
    for t in [1/8,math.pi/14,math.pi/2]:
        actual=m.wait(I.reshape(shape),t).reshape(288,288)
        prior=old.evolve(m0['hout'],t,I)
        residual=err(actual-np.exp(-2j*t)*prior)
        # This nonzero deletion defect is a full amplitude dictionary test;
        # no new relative-path interferometer permission is inferred.
        deletion=float(np.linalg.norm(actual-prior)/math.sqrt(288))
        predicted=2*abs(math.sin(t))
        err(deletion-predicted)
        phase_checks.append({'proper_time':t,'full_unitary_residual':residual,
          'omitting_work_mass_amplitude_defect':deletion,
          'exact_defect_formula':'2*abs(sin(proper_time))'})
    assert phase_checks[-1]['omitting_work_mass_amplitude_defect']>1.99
    # Real work mass, not an ignorable constant: same fourP and proper clock.
    eta=np.diag([1.,-1.,-1.,-1.])
    u=np.array([5/4,3/4,0.,0.]);sigma=F(7,10)
    full_phase=sum(eta[mu,mu]*u[mu]*hplus*float(sigma)*u[mu] for mu in range(4))
    common['proper_fourP_phase']=err(full_phase-float(sigma)*hplus)
    shift=[float(2*v) for v in u]
    assert shift==[2.5,1.5,0.,0.]
    common['mass_shell_after_shift']=err(sum(eta[mu,mu]*(u[mu]*hplus)@(u[mu]*hplus) for mu in range(4))-hplus@hplus)
    # Full energy certificate for n independent complete copies, no clock reset.
    upper=m.baseline+3+float(F(92,3))
    exact_upper='152/3 + 2*sqrt(19)'
    err(upper-(float(F(152,3))+2*math.sqrt(19)))
    inventory=[{'fresh_complete_copies':n,'energy_upper':n*upper,
                'symbolic_upper':f'{n}*({exact_upper})'} for n in [1,2,3,8]]
    # The old clock's energy-shell weights remain unchanged under new gates
    # commuting with H_D. This does NOT purport to measure a new shell effect.
    vals,vecs=np.linalg.eigh(m.hc)
    weights=np.abs(vecs[0,:])**2
    err(vals-np.arange(4));err(weights-np.array([1,3,3,1])/8)
    assert F(3,4)-F(1,2)==F(1,4)
    samples=[]
    for u0 in [F(1,5),F(1,2),F(1),F(2),F(7)]:
        middle=u0/(1+u0*u0)
        assert middle<=F(1,2)
        samples.append({'exp_minus_beta':str(u0),'gibbs_middle_weight':str(middle),
                        'gap_to_ready_middle_weight':str(F(3,4)-middle)})
    # Stored old output is inheritance evidence only. We do not run calculate().
    old_result=json.loads((BASE/'archive_1219_1233/1219/results.json').read_text('utf-8-sig'))
    return {'round':1271,'status':'PASS_AUTHOR_WITH_SCOPE','author_groups':1,
      'tolerance':TOL,'new_join':common,'physical_phase_checks':phase_checks,
      'boosted_added_work_fourP':shift,'full_single_task_energy_upper':upper,
      'full_single_task_energy_upper_exact':exact_upper,'fresh_copy_inventory':inventory,
      'ready_middle_energy_weight':'3/4','all_finite_beta_gibbs_middle_upper':'1/2',
      'original_D_gibbs_distance_lower':'1/4','finite_beta_diagnostics':samples,
      'inherited_thermal_result_status':old_result.get('status'),
      'old_thermal_groups_recomputed':0,
      'source_sha256':{rel:hashlib.sha256((BASE/rel).read_bytes()).hexdigest() for rel in SOURCES},
      'scope':[
        'One new complete-material join group; old heat exchange groups not repeated.',
        'All-history closure is a typed analytic induction, not a random numerical test.',
        'Two-dimensional fixed-motion internal thermal role; full four-dimensional carrier retained.',
        'The work mass +2 changes full fourP and phases; it is not a removable physical zero.',
        'Fresh complete copies include D, clock, stocks, seeds and records; no free reuse/reset.',
        'Original D Gibbs preparation and arbitrary actual TR menus are not inferred.',
        'No CT adoption, full-cognitive-programme refutation or GR goal completion.'
      ]}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write-once',action='store_true');args=parser.parse_args()
    result=calculate();dest=HERE/'results.json'
    if args.write_once:
        with dest.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        action='CREATED_ONCE'
    else:
        stored=json.loads(dest.read_text('utf-8-sig'))
        # Reproducibility within the declared numerical tolerance and exact metadata.
        def match(a,b):
            if isinstance(a,float) and isinstance(b,(float,int)):return abs(a-b)<=TOL
            if isinstance(a,dict):return a.keys()==b.keys() and all(match(a[k],b[k]) for k in a)
            if isinstance(a,list):return len(a)==len(b) and all(match(x,y) for x,y in zip(a,b))
            return a==b
        assert match(result,stored),'saved result mismatch'
        action='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'action':action,'author_groups':1,
      'full_task_energy_upper':result['full_single_task_energy_upper'],
      'D_Gibbs_gap':'1/4','new_adopted_axioms':0}))

if __name__=='__main__':main()
