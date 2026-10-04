"""Unfrozen 569 entry: inherited-scale hypercharge units and neutral Gauss mode.
Classical coefficient matching prescription, not lattice-to-MS quantum matching.
"""
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_scalar_propagation_matching as inherited
from protected_pair_rg import inputs,gauge_squared

def run():
    old,bound=inputs()
    row=json.loads((HERE.parent/'joint_singlet_common_mass_rg_results.json').read_text('utf8'))['examples'][2]
    gy2,gw2,gc2=gauge_squared(-row['u'],old['matched_inverse_couplings'])
    _,u,_=inherited.parameters()
    h2=float(u[0]);N=3;charge_h=N;factor=2*N
    # c=1, epsilon=v=1, k_h=k_s=1: explicit classical matching choice.
    bw=gw2/2; bu=gy2/(2*factor**2)
    beta_w=2/bw;beta_u=1/(2*bu)
    mu=h2/4
    A=np.diag([2*bw,2*bu])
    root=np.diag(np.sqrt(np.diag(A)))
    K=mu*np.outer([1.,-factor],[1.,-factor])
    physical=root@K@root
    values,basis=np.linalg.eigh(physical)
    target=h2*(gw2+gy2)/4
    assert abs(values[0])<1e-14 and abs(values[1]-target)<1e-14
    photon=np.array([np.sqrt(gy2),np.sqrt(gw2)])
    photon/=np.linalg.norm(photon)
    assert np.linalg.norm(physical@photon)<1e-14
    slopes=np.diag([bw*beta_w/2,2*bu*beta_u])
    assert np.max(abs(slopes-np.eye(2)))<1e-14
    rows=[]
    for lam in (0.,.2,3.):
        spectrum=np.linalg.eigvalsh(physical+lam*slopes)
        assert np.max(abs(spectrum-[lam,target+lam]))<1e-14
        rows.append(dict(graph_Laplacian=lam,neutral_transverse_frequencies_squared=spectrum.tolist()))
    # Nonzero longitudinal mode: one unbroken constraint remains.
    lam=3.;piw=.23;piu=-factor*piw;pphi=-np.sqrt(lam)*piw
    G=np.array([pphi+np.sqrt(lam)*piw,-factor*pphi+np.sqrt(lam)*piu])
    assert np.max(abs(G))<1e-14
    reduced_gap=(2*(bw+factor**2*bu)+4*lam/h2)*mu
    assert abs(reduced_gap-(target+lam))<1e-14
    charges=np.array([1,-4,2,-3,6,0,3])
    phase_error=float(np.max(abs(np.exp(2j*np.pi*charges)-1)))
    assert phase_error<2e-14
    wrong=h2*(gw2+factor**2*gy2)/4
    return dict(status='unfrozen_entry_not_completed_round',
        same_frozen_matter_scale_GeV=row['mu_GeV'],same_rg_u=row['u'],
        inherited_couplings_squared=dict(Y=float(gy2),weak=float(gw2),colour=float(gc2)),
        Higgs_integer_charge=charge_h,physical_Y_to_integer_Q=factor,
        b_weak=float(bw),b_circle=float(bu),beta_weak=float(beta_w),beta_circle=float(beta_u),
        squared_Higgs_amplitude=h2,neutral_mass_squared_matrix=physical.tolist(),
        neutral_zero_momentum_eigenvalues=values.tolist(),target_Z_gap_squared=float(target),
        original_SU2_only_W_gap_squared=float(gw2*h2/4),
        photon_direction_in_normalized_fields=photon.tolist(),rows=rows,
        longitudinal_Gauss_residual=float(np.max(abs(G))),
        longitudinal_Z_frequency_squared=float(reduced_gap),photon_longitudinal_mode_not_retained=True,
        integer_circle_identity_phase_error=phase_error,
        forgetting_charge_unit_Z_gap_squared=float(wrong),
        compact_full_colour_centre_quotient_not_yet_constructed=True,
        classical_matching_to_MS_coefficients_is_an_explicit_prescription=True,
        no_new_completed_scientific_round=True)

if __name__=='__main__':
    result=run()
    with (HERE/'electroweak_entry_probe_results.json').open('x',encoding='utf8',newline='\n') as f:
        f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
