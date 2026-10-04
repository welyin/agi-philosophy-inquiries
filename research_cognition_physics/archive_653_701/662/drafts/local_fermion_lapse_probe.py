"""662 entry: full original CAR contributes to the two-local-time bracket.

Conditional finite matrix calculation, not a full Gauss/path integral or ADM
closure proof. Reuses655 actual mass/hopping without redesigning a Hamiltonian.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_full_graph_transfer_sources as old


def blocks(data):
    h=data['h'];d=data['d'];d1=np.zeros_like(d);d2=np.zeros_like(d)
    d1[:32,:32]=d[:32,:32];d2[32:,32:]=d[32:,32:]
    a=old.nambu(h/2,d1);b=old.nambu(h/2,d2)
    assert old.norm(a+b-data['target_BdG'])<2e-13
    return a,b


def run():
    data=old.graph_data();a,b=blocks(data)
    current=(a@b-b@a)/1j
    n=np.array([1.2,.7]);m=np.array([.6,1.3]);factor=n[0]*m[1]-n[1]*m[0]
    an=n[0]*a+n[1]*b;am=m[0]*a+m[1]*b
    identity=old.norm((an@am-am@an)/1j-factor*current)
    assert identity<2e-13 and old.norm(current)>.05
    # Existing exact neutral invariant sector, four CAR modes, actual Fock check.
    neutral=old.graph_data(neutral=True);a,b=blocks(neutral)
    modes=[30,31,62,63];ids=modes+[x+64 for x in modes]
    aa=a[np.ix_(ids,ids)];bb=b[np.ix_(ids,ids)]
    def lift(matrix):
        h=matrix[:4,:4];d=matrix[:4,4:]
        return old.car.fock(h,d)-.5*np.trace(h)*np.eye(16)
    af=lift(aa);bf=lift(bb);jf=(af@bf-bf@af)/1j
    predicted=lift((aa@bb-bb@aa)/1j)
    fock_error=old.norm(jf-predicted);assert fock_error<2e-13
    h=af+bf;rho=old.car.exp_h(h,1.3);rho/=np.trace(rho)
    mean=np.trace(rho@jf);variance=np.trace(rho@jf@jf)-mean*mean
    assert abs(mean)<2e-13 and variance.real>1e-5 and abs(variance.imag)<1e-13
    # Honest real-time phase preparation under a selected local contribution.
    # This is a diagnostic unitary, not the original scalar s record instrument.
    ev,v=np.linalg.eigh(af);u=(v*np.exp(-.4j*ev))@v.conj().T
    prepared=u@rho@u.conj().T
    prepared_mean=np.trace(prepared@jf)
    assert abs(prepared_mean.imag)<1e-13 and abs(prepared_mean.real)>1e-6
    return dict(date='2026-10-02',status='662 actual entry, not complete round',
        original_two_node_CAR_modes=64,original_hopping_strength=.23,
        selected_allocation='onsite mass at each endpoint, hopping half at each endpoint',
        full_CAR_local_bracket_norm=old.norm(current),positive_lapses=[n.tolist(),m.tolist()],
        bilinear_lapse_bracket_identity_error=identity,
        neutral_Fock_identity_error=fock_error,
        neutral_stationary_current_mean=[float(mean.real),float(mean.imag)],
        neutral_stationary_current_variance=float(variance.real),
        neutral_locally_prepared_current=float(prepared_mean.real),
        fixed_background_not_full_bosonic_Gauss_state=True,
        additional_matrix_current_not_an_ADM_shift_identification=True,
        original_allocation_not_claimed_unique=True,
        dependencies={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
            ('joint_full_graph_transfer_sources.py','joint_fermion_gauss_completion.py',
             'research_note_588.md','research_note_623.md')})


if __name__=='__main__':
    result=run()
    with (HERE/'local_fermion_lapse_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(result,ensure_ascii=False))
