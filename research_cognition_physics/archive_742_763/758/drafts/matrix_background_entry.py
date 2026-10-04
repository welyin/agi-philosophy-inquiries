"""Unnumbered758 entry: original background and finite spectral source checks.

These are entry diagnostics, not a proof of full semiclassical dynamics.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
import joint_correlated_hadamard_noise as old
TARGET=HERE/'matrix_background_entry_results.json'

def run():
    q,psi,tensor,info,color=old.background.completed(12,1.)
    phi=q['phi'];h=np.linalg.norm(phi[...,:4],axis=-1);s=phi[...,4]
    F=old.matter_source.matter.original.F(phi)
    R=3*F/4*np.cos(2*s)-5*s/48*np.sin(2*s)
    # Uniform rectangle suffices, independent of sampling once field bounds hold.
    lower=Q(3,4)*Q(1477,800)*Q(7,25)-Q(3,40)
    assert lower==Q(25017,80000) and lower>0
    assert h.max()<.75 and np.max(abs(s))<.6 and R.min()>float(lower)
    point=phi[0,3,2];hh,dd=old.matter_source.matter.mass_matrices(point)
    sl=np.ix_(np.arange(24,32),np.arange(24,32));c=old.car(8)
    B=old.fock_quadratic(hh[sl],dd[sl],c)
    vac=np.zeros((64,64));vac[0,0]=1
    rho=np.kron(old.tau(.5,.5),vac);ref=np.kron(old.tau(.5,.25),vac)
    difference=rho-ref;Ys=old.matter_source.matter.Y['s'];fp=old.matter_source.matter.original.F(point)
    expected_second=.5*abs(Ys)**2*point[4]**2/fp
    actual_second=float(np.trace(difference@B@B).real)
    assert abs(actual_second-expected_second)<2e-14
    ev,U=np.linalg.eigh(B);weights=np.diag(U.conj().T@difference@U)
    rows=[]
    for u in (.1,.05,.025):
        gap=np.sum(weights*np.exp(1j*u*ev))
        rows.append(dict(spectral_test_parameter=u,gap=[float(gap.real),float(gap.imag)],
                         divided_by_u_squared=float(gap.real/u**2)))
    assert abs(rows[-1]['divided_by_u_squared']+actual_second/2)<2e-6
    deps=('research_note_525.md','research_note_574.md','research_note_726.md','research_note_727.md',
          'research_note_753.md','research_note_756.md','research_note_757.md','joint_reference_constraint_strata.py')
    return dict(type='unnumbered_entry_diagnostic',latest_scientific_round=757,
        actual753_field_rectangle=dict(h_min=float(h.min()),h_max=float(h.max()),s_min=float(s.min()),s_max=float(s.max()),
                                      R_min=float(R.min()),R_max=float(R.max()),conditional_uniform_R_lower=str(lower)),
        full_original_lepton_Fock=dict(dimension=256,initial_mean_gap=float(np.trace(difference@B).real),
            second_moment_gap=actual_second,analytic_second_gap=expected_second,characteristic_rows=rows,
            expected_small_parameter_ratio=-actual_second/2),
        dependency_hashes={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},
        single_point_mass_calibration_not_full_graph_spectrum=True,
        spectral_parameter_not_physical_time=True,
        continuous_or_semiclassical_dynamic_limit_proven=False,
        no_new_numbered_scientific_tests=True,all_entry_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('latest_scientific_round','all_entry_checks_passed','no_new_numbered_scientific_tests')}))
