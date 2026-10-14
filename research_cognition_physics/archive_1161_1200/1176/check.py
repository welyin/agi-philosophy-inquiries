"""1176 default: read-only audit recomputation, fixed-certificate bindings.

--write creates results.json once. --full-certificate additionally recomputes the
2048-step outward ODE certificate (~4 min); it still does not overwrite sources.
Default does not claim to rerun that integration or simulate the full K-body H.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
SCREEN=STAGE/'_admission'/'after1175_selection'
CAND=SCREEN/'tensor_candidate'
OLD=STAGE/'_admission'/'after1174_direction'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def calculate(full=False):
    pv=module('preview',CAND/'preview.py')
    sys.modules['preview']=pv
    sy=module('round1176_symmetry',CAND/'symmetry_check.py')
    nb=module('round1176_numerical_bound',CAND/'numerical_bound.py')
    first=json.loads((CAND/'results_preview.json').read_text(encoding='utf-8'))
    num=json.loads((CAND/'numerical_bound_checks.json').read_text(encoding='utf-8'))
    guard=json.loads((CAND/'polar_endpoint_guard.json').read_text(encoding='utf-8'))
    old_numeric=OLD/'mean_field_numerical_bound.py'
    for name,path in [('numerical_bound.py',CAND/'numerical_bound.py'),
                      ('preview.py',CAND/'preview.py'),
                      ('mean_field_numerical_bound.py',old_numeric)]:
        assert num['source_sha256'][name]==sha(path),(name,'changed source')
    for name,value in guard['source_sha256'].items():
        assert value==sha(CAND/name),(name,'changed endpoint source')
    assert guard['status']=='POSTHOC_POLAR_ENDPOINT_GUARD_PASS'
    assert guard['normalized_menu_bitwise_matched_float_leaves']==36
    assert guard['total_probability_allowance']<num['terminal_probability_arithmetic_padding']
    exact_entry=F(guard['exact_stored_minus_raw_entry_l1_rational'])
    assert F(guard['stored_minus_raw_operator_upper'])>exact_entry
    if full:
        rerun=nb.run(steps=2048)
        pv.compare(rerun,num)
        endpoint=module('round1176_endpoint_replay',CAND/'polar_endpoint_guard.py')
        pv.compare(endpoint.run(),guard)
    symmetry=sy.compute()
    pv.compare(symmetry,json.loads((CAND/'results_symmetry.json').read_text(encoding='utf-8')))
    _,cnorm,_=nb.interpolation_matrix()
    assert cnorm==num['interpolation_row_norm']
    assert num['source_modes_represented']==65**3 and num['compressed_modes']==33**3
    assert num['steps']==2048 and num['physics_unchanged']

    ell=F(16,65); V=65**3
    eps,nu,g,omega=F(3,20),F(7,20),F(1,5),F(1,2)
    bs=2*eps*2+6*nu/ell+4*g/ell
    bp=2*eps*F(4,5)+18*nu/ell+12*g/ell
    w=2*V*bs; lb=V*bp
    ETAK=1e-5; INST=2e-5; NUM=9e-5; KEXP=70881071
    log_required=2*(math.log10(4*(float(lb)+1e-8)*2/ETAK)
                   +12*(float(w)+1e-8)/math.log(10))
    assert log_required<KEXP
    q=2*np.sum(pv.WEIGHTS*pv.KABS[:,0])
    em0=4/65
    zssecond=-4*float(omega)*(float(eps)*em0+float(nu)*q)
    ztsecond=-4*float(omega)*float(g)*q
    base=(float(eps)-float(nu))*zssecond
    asecond=np.array([base-float(g)*ztsecond,base+float(g)*ztsecond,base])
    c2=4*(float(g)*ztsecond)**2
    assert c2>0.68
    assert max(abs(asecond-np.array(first['background']['active']['initial_a_second'])))<1e-12
    gram=np.diag([float(eps**2*4+3*nu**2/(2*ell**2)),float(g**2/ell**2)])
    fop=math.sqrt(float(max(np.diag(gram))))
    assert fop<num['Fop_upper']

    # Conservative outward decimal constants. They exceed the certified values.
    p_up=0.392699081698725
    d_up=0.000611563643
    J_up=0.045662192326534
    V_up=1.868888431192047
    assert J_up>num['J_integral_abs_gzt_upper']
    assert V_up>num['V_integral_vbar_plus_abs_gzt_upper']
    assert num['numerical_error_pair']<NUM
    lattice=float(np.nextafter(2*d_up*J_up+4*p_up*d_up*J_up*V_up,np.inf))
    budget=float(np.nextafter(lattice+2*ETAK+2*INST+NUM,np.inf))
    menu=num['normalized_original_menu']
    y=menu['y_clock_minusY_joint']; ml=menu['x_massless_position']
    lower_y=float(np.nextafter(abs(y['signed_gap'])-budget,-np.inf))
    lower_ml=float(np.nextafter(abs(ml['signed_gap'])-budget,-np.inf))
    old_scalar_upper=2*ETAK+2*INST
    excess_y=float(np.nextafter(lower_y-old_scalar_upper,-np.inf))
    assert lower_y>8.7e-4 and lower_ml>3e-3 and excess_y>8.1e-4
    scalar=first['scalar_original_menu']
    assert max(abs(x['signed_gap']) for x in scalar.values())<1e-14
    assert max(x['complete_48code_effect_distance'] for x in scalar.values())<1e-14
    assert max(symmetry['errors'].values())<1e-13
    # This is a conditional declared numerical arithmetic model, not a new
    # empirical instrument certification or a formal verification of all libm.
    paths=[
        STAGE/'1175'/'proof.md',OLD/'mean_field_bound.md',old_numeric,
        SCREEN/'tensor_covariance_bridge.md',SCREEN/'tensor_covariance_check.py',
        SCREEN/'tensor_covariance_checks.json',SCREEN/'correlated_lattice_bound.md',
        CAND/'preview.py',CAND/'results_preview.json',
        CAND/'symmetry_check.py',CAND/'results_symmetry.json',
        CAND/'numerical_bound.py',CAND/'numerical_bound_checks.json',CAND/'numerical_bound.md',
        CAND/'polar_endpoint_guard.py',CAND/'polar_endpoint_guard.json',CAND/'polar_endpoint_guard.md',
    ]
    return {
      'round':1176,'date':'2026-10-13','status':'PASS_WITH_EXPLICIT_EFFECTIVE_SCOPE',
      'sources':{str(p.relative_to(STAGE)).replace('\\','/'):sha(p) for p in paths},
      'group1_parent_source_kernel':{
        'source_CAR_modes_per_site':4,'probe_CAR_modes_per_site':12,
        'relation_spins_per_site_per_replica':2,'replica_diagonal_excluded':True,
        'bS_exact':str(bs),'bP_exact':str(bp),'w':float(w),'L':float(lb),
        'K_as_power_of_ten':KEXP,'sufficient_log10_K_outward':log_required,
        'etaK_per_branch':ETAK,'exact_BG_each_branch_retained':True,
        'same_full_BG_energy_spectral_distribution':True,
        'same_spectrum_of_full_task_with_fixed_probe_claimed':False,
        'full_Kbody_parent_simulated':False,
      },
      'group2_autonomous_common_transport':{
        'parameters':first['parameters'],'unknown_prepared_dimension':48,
        'all_mass_frequency_and_reference_coherence_kept':True,
        'Em_initial_exact':str(F(4,65)),'Etx_initial':float(q),
        'initial_zs_second':zssecond,'initial_zt_second':ztsecond,
        'initial_a_second':asecond.tolist(),'initial_comparator_Weyl_squared':c2,
        'Weyl_is_comparator_property_not_curvature_instrument_tomography':True,
        'source_Gram':gram.tolist(),'Fop':fop,
        'preview_background':first['background']['active'],
        'numeric_certificate_status':num['status'],
      },
      'group3_fixed_original_menu_and_budget':{
        'all_original_menus':menu,
        'certified_numeric_error_pair':num['numerical_error_pair'],
        'polar_endpoint_arithmetic_evidence':guard,
        'adopted_numeric_pair_budget':NUM,
        'dstar_up':d_up,'p_up':p_up,'J_up':J_up,'V_up':V_up,
        'correlated_lattice_difference_bound':lattice,
        'finite_K_pair_bound':2*ETAK,'instrument_per_branch_INPUT':INST,
        'total_pair_probability_budget':budget,
        'fixed_y_clock_site_lower':lower_y,'fixed_massless_x_lower':lower_ml,
        'old_scalar_finite_parent_and_instrument_difference_upper':old_scalar_upper,
        'new_vs_old_scalar_menu_separation_margin_y':excess_y,
        'single_branch_lattice_error_not_reduced_by_difference_bound':True,
        'actual_instruments_built':False,
      },
      'group4_reflection_and_deletion':{
        'full48_parity_check':symmetry,'scalar_rule_saved_menus':scalar,
        'scalar_kernel_exact_in_Hartree_probe_menu':True,
        'scalar_kernel_exact_in_full_finite_BG_claimed':False,
        'source_A_B_physical_change_keeps_original_probe_and_records_fixed':True,
        'common_passive_transport_changes_all_inputs_H_and_records':True,
      },
      'execution_scope':{
        'default':'recompute finite algebra, constants, full budget and source bindings; read stored independent ODE certificate',
        'full_recompute_command':'python -B -X utf8 check.py --full-certificate',
        'full_option_recomputes_certificate_without_rewriting_it':True,
      },
      'not_claimed':[
        'Einstein, healthy massless spin-2, complete HDA or nonlinear constraints',
        'unique generation of g, Omega, relation dynamics or all tensor directions',
        'unrestricted local-scalar no-go or arbitrary-background classicalization',
        'all continuum position PVMs, high-energy histories or curvature tomography',
        'practical K preparation, manufactured readout or dark-sector evidence',
      ],
    }

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--write',action='store_true')
    ap.add_argument('--full-certificate',action='store_true')
    args=ap.parse_args()
    result=calculate(full=args.full_certificate)
    target=HERE/'results.json'
    if args.write:
        with target.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
        mode='WRITE_ONCE'
    else:
        pv=module('round1176_readonly_compare',CAND/'preview.py')
        pv.compare(result,json.loads(target.read_text(encoding='utf-8')));mode='READ_ONLY_MATCH'
    out=result['group3_fixed_original_menu_and_budget']
    print(json.dumps({'status':mode,'groups':4,'full_certificate_recomputed':args.full_certificate,
        'y_lower':out['fixed_y_clock_site_lower'],'massless_lower':out['fixed_massless_x_lower'],
        'separation_from_scalar_upper_y':out['new_vs_old_scalar_menu_separation_margin_y'],
        'total_budget':out['total_pair_probability_budget']},indent=2))
