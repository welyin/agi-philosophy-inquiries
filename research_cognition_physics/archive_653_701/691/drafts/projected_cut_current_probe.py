"""691 entry: projected chiral cut variation and physicalY support.

Mature projected local chiral variation is inherited from Kikukawa-Yamada1998.
The actual neutral cut/source/same-time dictionary comparison is tested here.
Formal Jacobian-one identities do not establish a physical charge operator.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_nonflat_mass_measure as old
TARGET=HERE/'projected_cut_current_probe_results.json'


def expi(h,theta):
    e,u=np.linalg.eigh((h+h.conj().T)/2)
    return (u*np.exp(1j*theta*e))@u.conj().T


def run():
    links=np.broadcast_to(np.eye(16,dtype=complex),(2,4,16,16)).copy()
    u,v,d,hw=old.kernel(links)
    ev,vs=np.linalg.eigh(old.internal.spin.GAMMA[3]);wp=vs[:,ev>.5];wm=old.internal.spin.G5@wp
    jp=np.kron(np.kron(np.eye(4),(wp+wm)/np.sqrt(2)),np.eye(16))
    jm=np.kron(np.kron(np.eye(4),(wp-wm)/np.sqrt(2)),np.eye(16))
    K=jp.conj().T@d@v;A=jm.conj().T@v;pv=v@v.conj().T
    alpha=np.zeros((256,256),complex)
    for x in range(2):
        for spin in range(4):alpha[64*x+16*spin,64*x+16*spin]=1
    Tv=v.conj().T@alpha@v;Tb=jp.conj().T@alpha@jp
    assert abs(np.trace(Tv)-np.trace(Tb))<2e-12
    variation=K@Tv-Tb@K
    current=jp.conj().T@(d@alpha-alpha@d)@v
    identity=float(np.max(np.abs(variation-current)))
    assert identity<2e-13
    defect=A@Tv-Tb@A
    right=jm.conj().T@(pv-np.eye(256))@alpha@v
    assert np.max(abs(defect-right))<2e-13
    assert np.linalg.norm(defect,2)>.1
    rows=[];reference=np.linalg.det(K)
    neutral_t0=[32*x+16*s for x in range(2) for s in range(2)]
    neutral_t1=[64+i for i in neutral_t0]
    for theta in (.2,.7,1.2):
        changed=expi(Tb,-theta)@K@expi(Tv,theta)
        ratio=np.linalg.det(changed)/reference
        assert abs(ratio-1)<2e-12
        B=(changed-K)@np.linalg.inv(A)
        # B is the actualbarw B w source, computed fromoriginal chiral matrix.
        positive_negative=float(np.linalg.norm(B[np.ix_(neutral_t0,neutral_t1)],2))
        negative_positive=float(np.linalg.norm(B[np.ix_(neutral_t1,neutral_t0)],2))
        assert min(positive_negative,negative_positive)>.01
        # A product of normalized independent half-polynomials with no linear
        # terms cannot contain this nonzero mixed quadratic coefficient.
        covariance=-A@np.linalg.inv(K)
        via_Y=np.linalg.det(np.eye(128)-B@covariance)
        assert abs(via_Y-ratio)<2e-12
        naive=(1-(1-(np.sqrt(5)-2)**2)*np.sin(theta/2)**2)**2
        assert naive<1
        rows.append(dict(theta=theta,raw_projected_cut_ratio_real=float(ratio.real),
            unit_Ward_error=float(abs(ratio-1)),original_Y_source_error=float(abs(via_Y-ratio)),
            cross_time_source_norms=[positive_negative,negative_positive],
            unchanged689_two_local_charge_symbols=float(naive)))
    deps=('research_note_690.md','joint_spatial_charge_dictionary_results.json',
          'research_note_670.md','research_note_679.md','joint_nonflat_mass_measure.py')
    return dict(date='2026-10-02',entry_round=691,latest_formal_round=690,not_formal_round=True,
        maturity_source='https://arxiv.org/pdf/hep-lat/9808026',
        maturity_equations=['6','18-26','35-37'],
        projected_current_identity_error=identity,
        original_physical_local_phase_defect_norm=float(np.linalg.norm(defect,2)),
        neutral_measure_trace_difference=float(abs(np.trace(Tv)-np.trace(Tb))),
        rows=rows,all_checks_passed=True,
        source_only_variable_change_not_claimed_as_physical_charge=True,
        single_product_of_separate_halves_blocked_by_mixed_quadratic_coefficient=True,
        sum_of_products_or_extended_boundary_realization_not_excluded=True,
        general_chiral_measure_anomaly_not_assumed_absent=True,
        old_space_full_goal_and_four_branches_unchanged=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=691,all_checks_passed=True,
        Ward_error=max(r['unit_Ward_error'] for r in result['rows']),
        smallest_cross_time_source_norm=min(min(r['cross_time_source_norms']) for r in result['rows']))))
