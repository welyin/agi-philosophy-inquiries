"""690 entry: same originalY symbols with spatial propagation and a neutral mode.

Use670 actual256-dimensional Wilson matrix, original16-channel representation,
and689 spin convention. No auxiliary field is frozen in the analytic factorization.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_nonflat_mass_measure as old
TARGET=HERE/'spatial_car_source_probe_results.json'


def physical(links):
    u,v,d,h=old.kernel(links)
    ev,basis=np.linalg.eigh(old.internal.spin.GAMMA[3]);wp=basis[:,ev>.5];wm=old.internal.spin.G5@wp
    jp=np.kron(np.kron(np.eye(4),(wp+wm)/np.sqrt(2)),np.eye(16))
    jm=np.kron(np.kron(np.eye(4),(wp-wm)/np.sqrt(2)),np.eye(16))
    K=jp.conj().T@d@v;A=jm.conj().T@v;L=K+A
    C=-A@np.linalg.inv(K)
    # First internal even-exterior basis vector is the actual neutral mode.
    neutral=np.eye(16)[:,0]
    for g in links.reshape(-1,16,16):assert np.max(np.abs(g@neutral-neutral))<2e-13
    rows=np.zeros((2,128),complex)
    for t in range(2):
        for x in range(2):rows[t,32*(2*t+x)]=(-1)**x/np.sqrt(2)
    neutral_C=rows@C@rows.conj().T
    q=(A-K)@np.linalg.inv(L)
    return dict(K=K,A=A,L=L,C=C,q=q,neutral_C=neutral_C,neutral_rows=rows,
        gap=float(np.linalg.svd(h,compute_uv=False)[-1]))


def run():
    rows=[];first=None
    exact_amp=np.sqrt(5)-2
    for seed in (None,69001,69002,69003):
        links=np.zeros((2,4,16,16),complex)
        for mu in range(2):
            for i in range(4):links[mu,i]=np.eye(16) if seed is None else old.rep(*old.group(seed+10*mu+i,.03))
        f=physical(links)
        if first is None:first=f
        C=f['neutral_C'];expected=exact_amp*np.array([[0.,1.],[-1.,0.]])
        # Ordered sign follows the original AP shift; never replaceC by abs(C).
        error=float(np.max(np.abs(C-expected)))
        assert error<4e-12,(seed,C,expected)
        n=.5+.5*np.diag(C)
        n01=(1+C[0,0]+C[1,1]+C[0,0]*C[1,1]-C[0,1]*C[1,0])/4
        exact_n01=(1+exact_amp**2)/4
        assert max(abs(n-.5))<2e-12 and abs(n01-exact_n01)<3e-12
        # A rank-one H source: original689 symbol is1+(h-1)n_t.
        hs=[np.exp(.37j),np.exp(-.29j)]
        H=np.eye(128,dtype=complex)
        for t in range(2):
            p=f['neutral_rows'][t]
            H+=(hs[t]-1)*np.outer(p.conj(),p)
        D=(np.eye(128)+H)/2;B=(H-np.eye(128))/2
        actual=np.linalg.det(D@f['K']-B@f['A'])/np.linalg.det(f['K'])
        predicted=1+(hs[0]-1)*n[0]+(hs[1]-1)*n[1]+(hs[0]-1)*(hs[1]-1)*n01
        assert abs(actual-predicted)<3e-12
        # General chart exact factorization, not proof of a nearest-time transfer.
        chart=2**(-128)*np.linalg.det(f['L'])*np.linalg.det(np.eye(128)-H@f['q'])/np.linalg.det(f['K'])
        assert abs(chart-actual)<4e-12
        charge=[]
        for theta in (.31,.83,np.pi):
            Hg=np.eye(128,dtype=complex)
            for t in range(2):
                phase=np.exp((1 if t==0 else -1)*1j*theta)
                for x in range(2):
                    for spin in range(2):Hg[32*(2*t+x)+16*spin,32*(2*t+x)+16*spin]=phase
            Dg=(np.eye(128)+Hg)/2;Bg=(Hg-np.eye(128))/2
            two_insertions=np.linalg.det(Dg@f['K']-Bg@f['A'])/np.linalg.det(f['K'])
            neutral_prediction=(1-(1-exact_amp**2)*np.sin(theta/2)**2)**2
            assert abs(two_insertions-neutral_prediction)<3e-12
            charge.append(dict(theta=float(theta),original_inverse_charge_pair_real=float(two_insertions.real),
                expected_same_conserved_CAR_charge_pair=1.,
                original_formula_error=float(abs(two_insertions-neutral_prediction))))
        rows.append(dict(seed=seed,neutral_covariance_real=C.real.tolist(),
            neutral_covariance_imag_max=float(np.max(abs(C.imag))),
            neutral_free_covariance_error=error,mean_occupation_real=float(n[0].real),
            two_time_occupation_real=float(n01.real),same_conserved_projector_defect=float(abs(n01-n[0])),
            actual_source_error=float(abs(actual-predicted)),general_chart_error=float(abs(chart-actual)),
            original_Wilson_gap=f['gap'],neutral_charge_pairs=charge))
    neutral_q=first['neutral_rows']@first['q']@first['neutral_rows'].conj().T
    assert np.max(np.abs(np.diag(neutral_q)+2/np.sqrt(5)))<3e-12
    assert np.max(np.abs(first['C'][0:64,0:64]))<4e-12 # all first-time modes, free2-site spatial grid
    deps=('research_note_612.md','research_note_646.md','research_note_653.md','research_note_670.md',
          'research_note_673.md','research_note_679.md','research_note_680.md','research_note_689.md',
          'joint_nonflat_mass_measure.py','joint_physical_car_projector.py')
    return dict(date='2026-10-02',entry_round=690,latest_formal_round=689,not_formal_round=True,
        rows=rows,exact_neutral_amplitude='sqrt(5)-2',exact_two_time_occupation='(5-2*sqrt(5))/2',
        exact_gap_to_conserved_occupation='sqrt(5)-2',
        exact_neutral_parity_pair='(sqrt(5)-2)^4 = 161-72*sqrt(5)',
        neutral_effective_Q_real=neutral_q.real.tolist(),free_equal_time_all_physical_covariance_zero=True,
        gauge_neutral_factor_survives_entire_average_analytically=True,
        original_Q0_RP_failure_not_claimed=True,original_HF_not_identified=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=690,all_checks_passed=True,
        exact_two_time_occupation=result['exact_two_time_occupation'],
        max_source_error=max(r['actual_source_error'] for r in result['rows']))))
