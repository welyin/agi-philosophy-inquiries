"""Finite complete-H exchange diagnostic, applying rounds 421 and 1215.

No geometric readout or autonomous local apparatus is certified. Replay is
read-only; --write creates a result once. All battery boundary shells remain.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

SOURCE = Path(__file__).resolve()
RESULT = SOURCE.with_suffix('.json')


def norm(a):
    return float(np.linalg.norm(a, ord=2))


def run():
    eye = np.eye(2)
    x = np.array([[0., 1.], [1., 0.]])
    z = np.diag([1., -1.])
    n = np.diag([0., 1.])
    ns = np.kron(n, eye) + np.kron(eye, n)
    nd = np.kron(n, eye) - np.kron(eye, n)
    swap = np.array([[1.,0.,0.,0.], [0.,0.,1.,0.],
                     [0.,1.,0.,0.], [0.,0.,0.,1.]])
    V = np.kron(eye, swap)
    hsym = 2*np.eye(8) + .6*np.kron(x, np.eye(4)) + np.kron(eye, ns)
    h = hsym + .8*np.kron(z, nd)
    delta = V.T @ h @ V - h
    assert norm(delta + 1.6*np.kron(z, nd)) < 2e-14
    assert norm(V.T @ hsym @ V - hsym) < 2e-14
    assert norm(.6*(x@z-z@x)) > 1
    energies, modes = np.linalg.eigh(h)
    ticks = np.rint(5*energies).astype(int)
    assert np.max(np.abs(energies - ticks/5)) < 3e-15
    assert ticks.tolist() == [7,10,10,13,17,20,20,23]
    ue = modes.T @ V @ modes
    assert norm(ue.T @ ue - np.eye(8)) < 3e-14
    assert norm(ue.T - ue) < 3e-14
    margin = int(ticks.max()-ticks.min())
    diameter = float(energies.max()-energies.min())
    assert abs(diameter-16/5) < 2e-14
    assert abs(norm(delta)-8/5) < 2e-14
    cases = []
    for width in [65,257,1025]:
        size = width + 2*margin
        balances = np.arange(size)/5
        b = np.zeros(size)
        b[margin:margin+width] = 1/math.sqrt(width)
        # A[i,c,j] is W(|j> tensor |b>); only safe columns are populated.
        a = np.zeros((8,size,8))
        for j in range(8):
            for i in range(8):
                shift = int(ticks[j]-ticks[i])
                a[i,margin+shift:margin+width+shift,j] = ue[i,j]/math.sqrt(width)
        iso = a.reshape(8*size,8)
        ideal_iso = (ue[:,None,:] * b[None,:,None]).reshape(8*size,8)
        assert norm(iso.T @ iso-np.eye(8)) < 3e-13
        iso_error = norm(iso-ideal_iso)
        # Bound applies to all unknown inputs and arbitrary reference R.
        coarse_bound = math.sqrt(2*8*margin/width)
        assert iso_error <= coarse_bound+2e-13
        pull_h = np.einsum('icj,i,ick->jk', a, energies, a)
        pull_c = np.einsum('icj,c,ick->jk', a, balances, a)
        initial_c = float(np.dot(b*b,balances))
        energy_balance_gap = norm(pull_h+pull_c-np.diag(energies)-initial_c*np.eye(8))
        assert energy_balance_gap < 8e-11
        actual_delta = pull_h-np.diag(energies)
        ideal_delta = ue.T @ np.diag(energies) @ ue-np.diag(energies)
        energy_error = norm(actual_delta-ideal_delta)
        assert energy_error <= diameter*iso_error+1e-12
        worst_work = norm(pull_c-initial_c*np.eye(8))
        assert worst_work+diameter*iso_error >= 8/5-1e-12

        # Independent channel construction from overlaps of translated packets.
        choi = np.zeros((64,64))
        for i in range(8):
            for j in range(8):
                for l in range(8):
                    for k in range(8):
                        shift = int(ticks[j]-ticks[i]-ticks[k]+ticks[l])
                        overlap = max(0.,1-abs(shift)/width)
                        choi[8*i+j,8*l+k] = ue[i,j]*ue[l,k]*overlap/8
        entangled_columns = a.transpose(0,2,1).reshape(64,size)/math.sqrt(8)
        actual_choi = entangled_columns @ entangled_columns.T
        assert norm(choi-actual_choi) < 2e-14
        target = ue.reshape(64)/math.sqrt(8)
        choi_error = float(np.sum(np.abs(np.linalg.eigvalsh(
            actual_choi-np.outer(target,target))))/2)
        assert choi_error <= iso_error+2e-13

        # All finite total-energy shells are explicitly checked, including edges.
        full, boundary = 0, 0
        total_dimension = 0
        for energy_tick in range(int(ticks.min()), int(ticks.max())+size):
            members = [j for j in range(8) if 0 <= energy_tick-ticks[j] < size]
            if not members:
                continue
            block = ue if len(members)==8 else np.eye(len(members))
            total_dimension += len(members)
            assert norm(block.T @ block-np.eye(len(members))) < 3e-14
            if len(members)==8: full += 1
            else: boundary += 1
        assert total_dimension == 8*size
        # W is Hermitian involutive. G=(I-W)/20 >=0, ||G||<=1/10.
        # For tau=10*pi and Q ticks/5, exp(-i*tau*(Q+G))=W shellwise.
        g = (np.eye(8)-ue)/20
        eg, vg = np.linalg.eigh(g)
        assert eg.min() > -2e-14 and eg.max() <= .1+2e-14
        restored = (vg*np.exp(-1j*10*math.pi*eg)) @ vg.T
        assert norm(restored-ue) < 3e-14
        # One fixed physical input saturates the ideal material-energy change.
        physical = np.zeros(8); physical[1] = 1 # Zs=+, nA=0, nB=1
        initial = modes.T @ physical
        material_gain = float(initial @ actual_delta @ initial)
        battery_gain = float(initial @ (pull_c-initial_c*np.eye(8)) @ initial)
        assert abs(material_gain+battery_gain) < 8e-11
        cases.append(dict(
            packet_width=width, battery_levels=size,
            initial_battery_mean=initial_c, maximum_battery_energy=float(balances[-1]),
            full_shells=full, boundary_shells=boundary,
            all_reference_half_diamond_upper=iso_error,
            analytic_coarse_upper=coarse_bound,
            maximally_entangled_trace_error=choi_error,
            energy_balance_gap=energy_balance_gap,
            worst_absolute_battery_mean_change=worst_work,
            material_energy_operator_error=energy_error,
            witness_material_gain=material_gain, witness_battery_gain=battery_gain))

    assert cases[-1]['all_reference_half_diamond_upper'] < cases[0]['all_reference_half_diamond_upper']
    return dict(
        schema='resource_role_exchange_non_numbered_v1', status='PASS',
        source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        numbered_science_round=False,
        model=dict(omega='3/5', chi='4/5', shift=2, energies_ticks=ticks.tolist(),
                   energy_tick='1/5', source_defect_norm='8/5', energy_diameter='16/5',
                   noncommuting_internal_generator=True),
        analytic_contracts=dict(
            finite_arbitrary_H='1215 lifting in complete energy basis; finite Z^r packet',
            approximate_full_input='half diamond <= sqrt(2*d*r*M/N)',
            additive_energy_budget='delta + epsilon*diam(H) >= norm(Vdag H V-H)',
            exact_full_input='421: exact unitary and additive endpoint conservation force commutation',
            strong_budget_for_example='delta + (16/5)*epsilon >= 8/5'),
        cases=cases,
        autonomous_matrix=dict(tau='10*pi', persistent_G='(I-W)/20',
                               G_norm_le='1/10', G_energy_change=0,
                               status='finite mathematical realization, not local apparatus'),
        conclusion='Growing finite work resources do not select chi=0 in this permission test.',
        not_claimed=['A full cognitive-axiom countermodel', 'A certified geometry response',
                     'Spatial locality, causal apparatus or actual preparation',
                     'Einstein dynamics', 'Uniform resource budget or exact finite implementation'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        with RESULT.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        status='WRITTEN'
    else:
        assert json.loads(RESULT.read_text('utf-8')) == result
        status='READ_ONLY_MATCH'
    print(json.dumps(dict(status=status,numbered_science_round=False,
                          last_case=result['cases'][-1]),ensure_ascii=False))
