"""1051 finite diagnostics; the infinite-dimensional theorem is in proof.md.

Default: recompute and compare the stored result, read-only.
--write: create the result once. No plotting, no external dependencies beyond NumPy.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULT = HERE / 'effective_lie_interface_results.json'
HISTORY = (
    'archive_370_428/research_note_425.md',
    'archive_370_428/research_note_427.md',
    'archive_370_428/research_note_417.md',
    'archive_819_853/research_note_831.md',
    'archive_1009_1043/research_note_1021.md',
    'archive_1044_/research_note_1045.md',
    'archive_370_428/research_note_386.md',
    'archive_370_428/research_note_485.md',
)


def td(a, b):
    return float(np.sum(np.abs(np.linalg.eigvalsh((a-b + (a-b).conj().T)/2)))/2)


def pure(v):
    return np.outer(v, v.conj())


def twirl(rho, d, r, m):
    # Exact finite restriction of the infinite tail torus, not Haar sampling.
    mask = np.eye(d, dtype=bool)
    mask[:m+1, :m+1] = True
    tensor = rho.reshape(d, r, d, r)
    return (tensor * mask[:, None, :, None]).reshape(d*r, d*r)


def twirl_kraus(rho, d, r, m):
    projections = []
    head = np.zeros(d); head[:m+1] = 1
    projections.append(np.kron(np.diag(head), np.eye(r)))
    for j in range(m+1, d):
        p = np.zeros(d); p[j] = 1
        projections.append(np.kron(np.diag(p), np.eye(r)))
    answer = sum(p @ rho @ p for p in projections)
    completeness = sum(p @ p for p in projections)
    return answer, float(np.linalg.norm(completeness-np.eye(d*r)))


def rotate(d, j, angle):
    u = np.eye(d, dtype=complex)
    c, s = math.cos(angle), math.sin(angle)
    u[0, 0] = u[j, j] = c
    u[j, 0] = s; u[0, j] = -s
    return u


def apply_branches(branches, stage, d, r):
    answer = {}
    for history, rho in branches.items():
        if stage in (0, 2):
            u = rotate(d, 3, (.002 if stage == 0 else -.0015))
            if stage == 2 and history[-1] == 1:
                u = rotate(d, 3, .0025)
            a = np.kron(u, np.eye(r))
            answer[history] = a @ rho @ a.conj().T
        else:
            probabilities = .35 + .2*np.arange(d)/(d-1)
            for outcome in (0, 1):
                p = probabilities if outcome == 0 else 1-probabilities
                k = np.diag(np.sqrt(p))
                # An actual outcome-dependent poststate update, not just effects.
                u = rotate(d, 1, (.001 if outcome == 0 else -.001))
                a = np.kron(u @ k, np.eye(r))
                answer[history + (outcome,)] = a @ rho @ a.conj().T
    return answer


def branch_distance(a, b):
    assert a.keys() == b.keys()
    return sum(td(a[k], b[k]) for k in a)


def build():
    rng = np.random.default_rng(1051)
    d, r = 10, 3
    energies = np.array([0.] + [2.**j for j in range(1, d)])
    kr = np.kron(np.diag(energies), np.eye(r))
    maximums = dict(kraus=0., trace=0., idempotence=0., covariance=0.,
                    representative=0., k_moment=0., nested=0.)
    budget_rows = []
    for m in range(1, 7):
        max_distance = 0.
        for sample in range(7):
            z = rng.normal(size=(d-1, r)) + 1j*rng.normal(size=(d-1, r))
            z /= np.sqrt(energies[1:, None])
            z /= np.linalg.norm(z)
            tail_energy = float(np.sum(energies[1:, None] * np.abs(z)**2))
            p = .5/tail_energy
            ground = rng.normal(size=r) + 1j*rng.normal(size=r)
            ground /= np.linalg.norm(ground)
            psi = np.concatenate((math.sqrt(1-p)*ground, math.sqrt(p)*z.flatten()))
            rho = pure(psi)
            out = twirl(rho, d, r, m)
            independent, complete = twirl_kraus(rho, d, r, m)
            maximums['kraus'] = max(maximums['kraus'], float(np.linalg.norm(out-independent)), complete)
            maximums['trace'] = max(maximums['trace'], abs(float(np.trace(out).real)-1))
            maximums['idempotence'] = max(maximums['idempotence'], float(np.linalg.norm(twirl(out,d,r,m)-out)))
            assert np.linalg.eigvalsh(out).min() > -1e-12
            energy = float(np.trace(kr @ rho).real)
            assert abs(energy-.5) < 1e-12
            for power in (1, 2, 3):
                diag = np.repeat(energies**power, r)
                error = abs(float(np.dot(diag, np.diag(out-rho).real)))
                maximums['k_moment'] = max(maximums['k_moment'], error)
            phases = rng.uniform(-np.pi,np.pi,size=d); phases[0] = 0
            ug = np.kron(np.diag(np.exp(1j*phases)), np.eye(r))
            conjugate = ug @ rho @ ug.conj().T
            maximums['covariance'] = max(maximums['covariance'],
                float(np.linalg.norm(twirl(conjugate,d,r,m)-ug@out@ug.conj().T)))
            same_rep = phases.copy(); same_rep[m+1:] += rng.uniform(-np.pi,np.pi,size=d-m-1)
            uh = np.kron(np.diag(np.exp(1j*same_rep)), np.eye(r))
            maximums['representative'] = max(maximums['representative'],
                float(np.linalg.norm(ug@out@ug.conj().T-uh@out@uh.conj().T)))
            fine = twirl(rho,d,r,min(m+1,d-1))
            maximums['nested'] = max(maximums['nested'],float(np.linalg.norm(twirl(fine,d,r,m)-out)))
            max_distance = max(max_distance,td(out,rho))
        bound = min(1.,2*math.sqrt(.5/2**(m+1)))
        assert max_distance <= bound+1e-12
        budget_rows.append(dict(head_last_index=m, budget=.5, samples=7,
                               max_joint_distance=max_distance, proven_upper_bound=bound))

    reference_rows = []
    for m in range(1, 7):
        p = .5/2**(m+1)
        psi = np.zeros(d*2,complex); psi[0] = math.sqrt(1-p); psi[(m+1)*2+1] = math.sqrt(p)
        rho = pure(psi); out = twirl(rho,d,2,m)
        target = math.sqrt(p*(1-p))
        actual = td(rho,out)
        reduced = np.trace((rho-out).reshape(d,2,d,2),axis1=1,axis2=3)
        assert abs(actual-target) < 1e-12 and np.linalg.norm(reduced) < 1e-12
        reference_rows.append(dict(m=m,p=p,joint_distance=actual,formula=target,
                                   system_marginal_distance=td(reduced,np.zeros_like(reduced))))

    # Four true operations, two measurements, four retained classical branches.
    hd, hr, hm = 6, 2, 2
    hk = np.kron(np.diag([0]+[2**j for j in range(1,hd)]),np.eye(hr))
    start = np.zeros(hd*hr,complex);start[0]=math.sqrt(.999);start[3]=math.sqrt(.001)
    true = {():pure(start)}
    effective = {():twirl(true[()],hd,hr,hm)}
    true_energies = []; local_losses = []
    for j in range(5):
        true_energies.append(float(sum(np.trace(hk@rho).real for rho in true.values())))
        averaged = {h:twirl(rho,hd,hr,hm) for h,rho in true.items()}
        local_losses.append(branch_distance(true,averaged))
        if j < 4:
            true = apply_branches(true,j,hd,hr)
            effective = apply_branches(effective,j,hd,hr)
            effective = {h:twirl(rho,hd,hr,hm) for h,rho in effective.items()}
    actual_history = branch_distance(true,effective)
    hybrid_bound = sum(local_losses)
    budget_bound = sum(min(1.,2*math.sqrt(e/2**(hm+1))) for e in true_energies)
    assert len(true) == 4 and actual_history <= hybrid_bound+1e-12 <= budget_bound+1e-12
    assert abs(sum(np.trace(x).real for x in true.values())-1) < 1e-12
    assert max(true_energies) < .003

    unbudgeted = []
    for length in (2,4,8,16,32):
        psi = np.ones(length)/math.sqrt(length)
        distance = td(pure(psi),np.eye(length)/length)
        assert abs(distance-(1-1/length)) < 1e-12
        unbudgeted.append(dict(tail_levels=length,distance=distance,formula=1-1/length))

    source_rows = []
    delta = .01
    for length in (4,8,16,32):
        weights = 2.**(-np.arange(1,length+1,dtype=float)/2)
        weights /= np.linalg.norm(weights)
        v = np.concatenate(([delta], math.sqrt(1-delta**2)*weights))
        p = pure(v); q = np.eye(length+1)-p
        rho = np.zeros((length+1,length+1));rho[0,0]=1
        out = p@rho@p + q@rho@q
        kvals = np.array([0.]+[2.**j for j in range(1,length+1)])
        cost = float(np.dot(kvals,np.diag(out)))
        analytic = 2*delta**2*(1-delta**2)*length/(1-2.**(-length))
        assert abs(cost-analytic) < 1e-12
        assert abs(td(rho,out)-delta*math.sqrt(1-delta**2)) < 1e-12
        source_rows.append(dict(tail_levels=length,trace_distance=td(rho,out),
                                output_budget=cost,exact_expression_value=analytic))

    plus = pure(np.array([1.,1.])/math.sqrt(2))
    hadamard = np.array([[1.,1.],[1.,-1.]])/math.sqrt(2)
    left = np.diag(np.diag(hadamard@plus@hadamard))
    right = np.diag(np.diag(hadamard@np.diag(np.diag(plus))@hadamard))
    incompatible = td(left,right)
    assert abs(incompatible-.5)<1e-12
    assert max(maximums.values()) < 1e-11
    history = {}
    for name in HISTORY:
        path = ROOT/name
        # 485 was archived separately; resolve its known name without copying it.
        if not path.exists():
            found = list(ROOT.glob('archive_*/'+Path(name).name))
            assert len(found)==1,(name,found)
            path = found[0]
        history[path.relative_to(ROOT).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return dict(round=1051,date='2026-10-08',all_scientific_checks_passed=True,
        new_science_groups=1,new_cognitive_axioms=0,roadmap_complete=False,
        diagnostics_are_finite_restrictions_not_infinite_dimensional_proof=True,
        haar_twirl_is_effective_dictionary_not_free_physical_control=True,
        effective_quantum_state_space_claimed_finite=False,
        historical_sha256=history,max_residuals=maximums,budget_samples=budget_rows,
        reference_witnesses=reference_rows,
        finite_history=dict(operations=4,retained_record_branches=4,
            true_prefix_budgets=true_energies,local_joint_losses=local_losses,
            full_final_joint_distance=actual_history,exact_hybrid_upper_bound=hybrid_bound,
            budget_upper_bound=budget_bound,effective_prefix_budget_required=False),
        no_budget_witnesses=unbudgeted,source_failure_truncations=source_rows,
        cross_scale_dynamics_nonintertwining_distance=incompatible)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert np.isclose(a,b,rtol=2e-9,atol=2e-11),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=build()
    if args.write:
        with RESULT.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(result,json.loads(RESULT.read_text(encoding='utf8')))
    print(json.dumps(dict(round=1051,all_scientific_checks_passed=True,
        mode='exclusive_write' if args.write else 'read_only_compare',
        max_residual=max(result['max_residuals'].values()),
        finite_history=result['finite_history'],new_science_groups=1),ensure_ascii=False))


if __name__=='__main__':main()
