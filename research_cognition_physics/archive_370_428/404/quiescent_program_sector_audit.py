"""Round 404: normal quiescent-sector states versus an infinite classical program.

Direct implementation of Schaeffer (2015), section 6, for small finite cases.
The state-space and all-time claims are analytic; simulations do not prove them.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

SPEEDS=(-3,-2,-1,1,2,3)
TARGET=Path(__file__).with_name('quiescent_program_sector_audit_results.json')


def bits(j):
    return [(j >> (5-k)) & 1 for k in range(6)]


def integer(a):
    return sum(int(x) << (5-k) for k,x in enumerate(a))


def local_column(j):
    a=bits(j)
    b=a.copy()
    if all(a[:3]) or all(a[3:]):
        if all(a[:3]):
            b[3],b[4],b[5]=a[5],a[3],a[4]
        if all(a[3:]):
            b[0],b[1],b[2]=a[2],a[0],a[1]
        return [(integer(b),1.)]
    if a[4]+a[5] == 1 and a[1]+a[0] == 1:
        if (a[4] and a[1]) or (a[5] and a[0]):
            b[2],b[3]=a[3],a[2]
        elif a[4] and a[0]:
            b[2] ^= a[3]
        else:
            b[3] ^= a[2]
        return [(integer(b),1.)]
    h_target=None
    if a[4] and a[5] and a[0] and not a[1]:
        h_target=2
    elif a[5] and a[1] and a[0] and not a[4]:
        h_target=3
    if h_target is not None:
        out=[]
        for value in (0,1):
            b=a.copy()
            b[h_target]=value
            out.append((integer(b),(-1.)**(a[h_target]*value)/np.sqrt(2)))
        return out
    if a[4] and a[5] and a[1] and not a[0]:
        return [(j,np.exp(1j*np.pi*a[2]/4))]
    if a[4] and a[1] and a[0] and not a[5]:
        return [(j,np.exp(1j*np.pi*a[3]/4))]
    return [(j,1.)]


TABLE=[local_column(j) for j in range(64)]
REVERSE_TABLE=[[(j,np.conj(a)) for j in range(64) for i,a in TABLE[j] if i==out]
               for out in range(64)]


def local_matrix():
    u=np.zeros((64,64),complex)
    for j,entries in enumerate(TABLE):
        for i,amplitude in entries:
            u[i,j]=amplitude
    return u


def config_at(position,j):
    return tuple((position,k) for k,b in enumerate(bits(j)) if b)


def shift(config):
    return tuple(sorted((p+SPEEDS[k],k) for p,k in config))


def step(state,reverse=False):
    if reverse:
        state={tuple(sorted((p-SPEEDS[k],k) for p,k in c)):v for c,v in state.items()}
    table=REVERSE_TABLE if reverse else TABLE
    answer={}
    for config,vector in state.items():
        branches=[(config,1.+0j)]
        for position in sorted({p for p,k in config}):
            new=[]
            for current,amplitude in branches:
                occupied={k for p,k in current if p==position}
                j=integer([int(k in occupied) for k in range(6)])
                other=tuple((p,k) for p,k in current if p!=position)
                for out,coefficient in table[j]:
                    updated=tuple(sorted(other+config_at(position,out)))
                    new.append((updated,amplitude*coefficient))
            branches=new
        for current,amplitude in branches:
            moved=current if reverse else shift(current)
            answer[moved]=answer.get(moved,np.zeros_like(vector))+amplitude*vector
    return {c:v for c,v in answer.items() if np.linalg.norm(v)>1e-14}


def evolve(state,t):
    for _ in range(abs(t)):
        state=step(state,reverse=t<0)
    return state


def norm2(state):
    return float(sum(np.vdot(v,v).real for v in state.values()))


def reference_state(state):
    return sum(np.outer(v,v.conj()) for v in state.values())


def nonvacuum_probability(state,radius=0):
    return float(sum(np.vdot(v,v).real for c,v in state.items()
                     if any(abs(p)<=radius for p,k in c)))


def local_reference_distance_to_vacuum(state,radius=0):
    local_keys={()}
    blocks={}
    for config,vector in state.items():
        here=tuple((p,k) for p,k in config if abs(p)<=radius)
        outside=tuple((p,k) for p,k in config if abs(p)>radius)
        local_keys.add(here)
        blocks.setdefault(outside,{})[here]=vector
    keys=sorted(local_keys)
    index={key:j for j,key in enumerate(keys)}
    refdim=len(next(iter(state.values())))
    reduced=np.zeros((len(keys)*refdim,)*2,complex)
    for block in blocks.values():
        v=np.zeros(len(keys)*refdim,complex)
        for key,vector in block.items():
            start=index[key]*refdim
            v[start:start+refdim]=vector
        reduced+=np.outer(v,v.conj())
    target=np.zeros_like(reduced)
    start=index[()]*refdim
    target[start:start+refdim,start:start+refdim]=reference_state(state)
    return float(np.abs(np.linalg.eigvalsh(reduced-target)).sum()/2)


def choi_input(program=None):
    if program is None:
        program={():1.}
    state={}
    for j in range(64):
        v=np.zeros(64,complex)
        v[j]=1/8
        for config,amplitude in program.items():
            assert all(p!=0 for p,k in config)
            key=tuple(sorted(config+config_at(0,j)))
            state[key]=amplitude*v
    return state


def delayed_program(delay):
    # At t=delay three controls meet at cell 0; before that no two meet.
    return tuple(sorted((-SPEEDS[k]*delay,k) for k in (0,1,5)))


@lru_cache(None)
def report():
    u=local_matrix()
    boundary=[]
    for half in (range(3),range(3,6)):
        selected=[j for j in range(64) if all(bits(j)[k]==0 for k in half)]
        outside=[j for j in range(64) if j not in selected]
        boundary.append(float(np.linalg.norm(u[np.ix_(outside,selected)])))
    full=[]
    for t in (-8,-3,-1,1,3,8):
        state=evolve(choi_input(),t)
        full.append(dict(time=t,nonvacuum=nonvacuum_probability(state),
                         reference_error=float(np.linalg.norm(reference_state(state)-np.eye(64)/64)),
                         choi_receiver_distance_to_reset=local_reference_distance_to_vacuum(state),
                         norm_error=abs(norm2(state)-1)))
    rng=np.random.default_rng(404)
    front=[]
    for radius in (1,2,3):
        initial={}
        for j in range(7):
            occupied=tuple((p,k) for p in range(-radius,radius+1) for k in range(6)
                           if rng.random()<.24)
            initial[occupied]=rng.normal(size=2)+1j*rng.normal(size=2)
        scale=np.sqrt(norm2(initial))
        state={c:v/scale for c,v in initial.items()}
        r0=reference_state(state)
        for t in range(1,radius+3):
            state=step(state)
            violations=0
            for config in state:
                violations+=sum((SPEEDS[k]>0 and p<t-radius) or
                                (SPEEDS[k]<0 and p>radius-t) for p,k in config)
            front.append(dict(initial_radius=radius,time=t,
                              configuration_count=len(state),front_violations=int(violations),
                              central_nonvacuum=nonvacuum_probability(state),
                              norm_error=abs(norm2(state)-1),
                              reference_error=float(np.linalg.norm(reference_state(state)-r0))))
    tails=[]
    for q,t in ((.01,5),(.09,5),(.01,-5),(.09,-5)):
        initial=choi_input({():np.sqrt(1-q),delayed_program(t):1j*np.sqrt(q)})
        actual=evolve(initial,t)
        tails.append(dict(program_tail_probability=q,cutoff_radius=0,time=t,
                          complete_reference_choi_distance=local_reference_distance_to_vacuum(actual),
                          half_diamond_upper_bound=np.sqrt(q),
                          central_nonvacuum=nonvacuum_probability(actual),
                          reference_error=float(np.linalg.norm(reference_state(actual)-np.eye(64)/64)),
                          norm_error=abs(norm2(actual)-1)))
    delayed=[]
    for t in (2,5,11,23):
        state=evolve({delayed_program(t):np.ones(1,complex)},t)
        delayed.append(dict(delay=t,initial_radius=3*t,
                            nonvacuum_at_delay=nonvacuum_probability(state),
                            norm_error=abs(norm2(state)-1)))
    # Finite all-one prefixes on eight positive-speed sites; global vectors are
    # orthogonal while every fixed finite subsystem stabilizes.
    prefix=[]
    vectors=[]
    for n in range(1,9):
        v=np.zeros(256,complex)
        v[sum(1 << (7-k) for k in range(n))]=1.
        vectors.append(v)
    for n in (2,4,6):
        v,w=vectors[n-1],vectors[n]
        vr=v.reshape(4,64)
        wr=w.reshape(4,64)
        prefix.append(dict(prefix=n,next_prefix=n+1,
                           global_vector_distance=float(np.linalg.norm(v-w)),
                           global_pure_trace_distance=float(np.sqrt(1-abs(np.vdot(v,w))**2)),
                           first_two_cells_marginal_distance=float(np.linalg.norm(vr@vr.conj().T-wr@wr.conj().T))))
    return dict(round=404,
                scope='Schaeffer quiescent-sector LQCA: normal fixed programs give local vacuum asymptotics, incompatible with arbitrary-accuracy repeated unitary return except finite exact hits. Infinite classical program limits require a different state-space completion. No general non-three-dimensional cognition or GR countermodel is claimed.',
                local_unitarity_error=float(np.linalg.norm(u.conj().T@u-np.eye(64))),
                quiescent_error=float(np.linalg.norm(u[:,0]-np.eye(64)[:,0])),
                one_sign_vacuum_invariance_errors=boundary,
                complete_input_checks=full,finite_front_checks=front,
                coherent_program_tail_checks=tails,nonuniform_delay_checks=delayed,
                infinite_program_prefix_checks=prefix,
                ordinary_finite_task_universality_rejected=False,
                normal_quiescent_sector_explicit=True,
                all_fixed_normal_programs_have_local_reset_limit=True,
                local_reset_implies_global_information_destruction=False,
                convergence_uniform_over_all_program_states=False,
                infinite_classical_program_belongs_to_original_hilbert_sector=False,
                strong_universality_corollary_requires_state_space_clarification=True,
                all_infinite_sector_extensions_excluded=False,
                autonomous_continuous_hamiltonian_derived=False,
                full_cognition_to_gr_refuted=False,spatial_dimension_generated=False)


class Audit(unittest.TestCase):
    def test_local_quantum_rule_and_both_vacuum_fronts(self):
        r=report()
        self.assertLess(r['local_unitarity_error'],1e-12)
        self.assertLess(r['quiescent_error'],1e-12)
        self.assertLess(max(r['one_sign_vacuum_invariance_errors']),1e-12)

    def test_full_six_qubit_input_is_transferred_out_not_destroyed(self):
        for row in report()['complete_input_checks']:
            self.assertLess(row['nonvacuum'],1e-12)
            self.assertLess(row['reference_error'],1e-12)
            self.assertLess(row['choi_receiver_distance_to_reset'],1e-12)
            self.assertLess(row['norm_error'],1e-12)

    def test_quantum_front_bound_for_finite_support_with_reference(self):
        for row in report()['finite_front_checks']:
            self.assertEqual(row['front_violations'],0)
            self.assertLess(row['norm_error'],1e-11)
            self.assertLess(row['reference_error'],1e-11)
            if row['time']>row['initial_radius']:
                self.assertLess(row['central_nonvacuum'],1e-12)

    def test_coherent_tail_bound_with_complete_input_reference(self):
        for row in report()['coherent_program_tail_checks']:
            self.assertLessEqual(row['complete_reference_choi_distance'],row['half_diamond_upper_bound']+1e-11)
            self.assertAlmostEqual(row['complete_reference_choi_distance'],np.sqrt(row['program_tail_probability']))
            self.assertAlmostEqual(row['central_nonvacuum'],row['program_tail_probability'])
            self.assertLess(row['reference_error'],1e-11)
            self.assertLess(row['norm_error'],1e-11)

    def test_no_uniform_reset_time_over_all_finite_programs(self):
        for row in report()['nonuniform_delay_checks']:
            self.assertAlmostEqual(row['nonvacuum_at_delay'],1.)
            self.assertEqual(row['initial_radius'],3*row['delay'])
            self.assertLess(row['norm_error'],1e-12)

    def test_infinite_program_prefixes_are_locally_consistent_but_not_norm_cauchy(self):
        for row in report()['infinite_program_prefix_checks']:
            self.assertAlmostEqual(row['global_vector_distance'],np.sqrt(2))
            self.assertAlmostEqual(row['global_pure_trace_distance'],1.)
            self.assertLess(row['first_two_cells_marginal_distance'],1e-12)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    checks=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not checks.wasSuccessful():
        raise SystemExit(1)
    result=dict(report())
    result['checks']=dict(run=checks.testsRun,failures=len(checks.failures),errors=len(checks.errors))
    result['runtime']=dict(python=platform.python_version(),numpy=np.__version__)
    if args.write_results:
        with TARGET.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
