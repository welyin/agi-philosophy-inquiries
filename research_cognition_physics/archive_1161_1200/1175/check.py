"""Round1175: default read-only recomputation; --write creates results once.

Runs the fixed source/probe comparator and the separate outward residual bound.
Never constructs or claims to simulate the astronomical full K-body parent.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
SCREEN=HERE.parent/'_admission'/'after1174_direction'
CAND=SCREEN/'mean_field_candidate'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def calc():
    preview=module('round1175_frozen_preview',CAND/'preview.py')
    ext=module('round1175_fixed_direct_sum',CAND/'unified_probe_extension.py')
    nb=module('round1175_numerical_certificate',SCREEN/'mean_field_numerical_bound.py')
    first=preview.calculate()
    preview.compare(first,json.loads((CAND/'results_preview.json').read_text(encoding='utf-8')))
    extended=ext.compute()
    preview.compare(extended,json.loads((CAND/'results_unified_probe.json').read_text(encoding='utf-8')))
    numerical=nb.run(steps=2048)
    # Strictly larger than ell^2 (2 pi / 16)^3 / 6 using pi < 3.141592653589794.
    dstar=0.000611563643
    tau_error=numerical['clock_tau_error']
    xi_error=numerical['clock_xi_error']
    tau_upper=np.nextafter(numerical['stored_tau']+tau_error,np.inf)
    xi_upper=np.nextafter(numerical['stored_xi']+xi_error,np.inf)
    eta_lattice=np.nextafter(dstar*(tau_upper+xi_upper),np.inf)
    # M=0 block is analytic, so its numerical bound uses the same integrals.
    p_upper=2*3.141592653589794/16
    zero_active_error=p_upper*xi_error
    zero_flat_error=p_upper*tau_error
    uniform_actual=max(numerical['probe_total_numerical_error'],zero_active_error)
    uniform_flat=max(numerical['flat_clock_calibration_operator_error'],zero_flat_error)
    # Includes reported Gram/normalization/effect arithmetic defects conservatively.
    arithmetic_allowance=1e-9
    numerical_sum_bound=uniform_actual+uniform_flat+arithmetic_allowance
    NUM=9e-5
    assert numerical_sum_bound < NUM,(numerical_sum_bound,NUM)
    ETAK=1e-4; INST=2e-5
    budget=float(np.nextafter(eta_lattice+2*ETAK+2*INST+NUM,np.inf))
    y=first['active_mobile_clock_and_site']['y_joint_clock_minusY']
    massless=first['massless_site_and_control']['x_position']
    lower_y=abs(y['signed_gap'])-budget
    lower_m0=abs(massless['signed_gap'])-budget
    assert lower_y>0 and lower_m0>0,(lower_y,lower_m0)
    assert extended['full_48_column_gram_error']<1e-10
    assert first['pure_lapse_proper_time_control']['full_prepared_unitary_error']<1e-10
    assert first['vacuum_flat_control']['full_prepared_unitary_error']<1e-10
    for branch in first['background'].values():
        assert branch['diagnostic_residuals']['energy_per_site']<1e-7
    # Passive relabelling transports the state AND the complete effect.
    # This is a zero-change control, not a different physical source.
    active_for_reference=preview.evolve(preview.NU)
    e=preview.effect_matrix(active_for_reference['u'],1,True)
    phases=np.repeat(np.exp(2j*np.pi*preview.MODES[:,1]/preview.GRID),8)
    translated=phases*preview.WITNESS
    e_translated=phases[:,None]*e*phases.conj()[None,:]
    ref_error=float(abs(np.vdot(translated,e_translated@translated)
                        -np.vdot(preview.WITNESS,e@preview.WITNESS)))
    assert ref_error<1e-12
    source_files=[
        SCREEN/'mean_field_bound.md',SCREEN/'mean_field_numerical_bound.md',
        SCREEN/'mean_field_numerical_bound.py',CAND/'preview.py',CAND/'results_preview.json',
        CAND/'unified_probe_extension.py',CAND/'results_unified_probe.json',
    ]
    sources={str(p.relative_to(HERE.parent)).replace('\\','/'):
             hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    return {
      'round':1175,'date':'2026-10-13','status':'PASS_WITH_EXPLICIT_EFFECTIVE_SCOPE',
      'sources':sources,
      'group1_parent_and_full_background_bridge':{
        'full_CAR_modes_source_per_site':4,'full_CAR_modes_probe_per_site':12,
        'replica_diagonal_terms_included':False,'exact_BG_output_kept':True,
        'unknown_probe_reference_kept':True,'resource_constants':extended['resources'],
        'K_defined_as_power_of_10':52275242,'finite_K_budget_per_branch':ETAK,
        'full_Kbody_parent_numerically_simulated':False,
      },
      'group2_autonomous_common_transport':{
        'parameters':first['parameters'],'unified_mass_labels':[0,.3,.5],
        'unknown_prepared_dimension':48,'background':first['background']['active'],
        'full_48_column_gram_error':extended['full_48_column_gram_error'],
        'cross_mass_coherent_norm_error':extended['cross_mass_coherent_norm_error'],
        'old_menu_reuse_error':extended['maximum_old_menu_reuse_error'],
        'initial_comparator_tidal':first['analytic_initial_tidal'],
        'all_three_mobile_clock_effects':first['active_mobile_clock_and_site'],
      },
      'group3_fixed_menu_and_error_certificate':{
        'numerical_enclosure':numerical,'uniform_48_code_numeric_sum_bound':numerical_sum_bound,
        'numeric_budget':NUM,'arithmetic_and_polar_allowance':arithmetic_allowance,
        'lattice_dstar_upper':dstar,'tau_upper':float(tau_upper),'xi_upper':float(xi_upper),
        'two_branch_lattice_error':float(eta_lattice),'two_branch_finite_K_error':2*ETAK,
        'instrument_contract_error_per_branch_INPUT':INST,
        'total_two_branch_probability_budget':budget,
        'fixed_y_clock_site_probabilities':y,'fixed_massless_x_probabilities':massless,
        'certified_y_joint_gap_lower':lower_y,'certified_massless_x_gap_lower':lower_m0,
        'x_joint_gap_is_NOT_claimed_above_budget':True,
        'instruments_built_in_this_round':False,
      },
      'group4_delimited_controls':{
        'pure_lapse_comparator':first['pure_lapse_proper_time_control'],
        'vacuum_comparator':first['vacuum_flat_control'],
        'massless_all_axes':first['massless_site_and_control'],
        'finite_K_controls_exactly_flat':False,
        'passive_reference_joint_probability_error':ref_error,
        'passive_reference_change_requires_joint_transport_of_state_H_instruments':True,
      },
      'not_claimed':[
        'Einstein dynamics or completed GR hard deliverables',
        'dimension generation or arbitrary background classicalization',
        'all continuum position PVMs and their arbitrary intervening histories',
        'practical preparation of K replicas or manufactured instruments',
        'dark energy or observational cosmology evidence',
      ],
    }

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=calc();target=HERE/'results.json'
    if args.write:
        with target.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False)
        mode='WRITE_ONCE'
    else:
        p=module('round1175_readonly_comparison',CAND/'preview.py')
        p.compare(result,json.loads(target.read_text(encoding='utf-8')));mode='READ_ONLY_MATCH'
    g=result['group3_fixed_menu_and_error_certificate']
    print(json.dumps({'status':mode,'groups':4,'y_lower':g['certified_y_joint_gap_lower'],
                      'massless_lower':g['certified_massless_x_gap_lower'],
                      'total_budget':g['total_two_branch_probability_budget']},indent=2))
