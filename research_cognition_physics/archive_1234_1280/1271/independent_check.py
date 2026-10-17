"""1271 independent finite identity checks; no author-code imports.

One independent group: baseline, transported resources/phases, and the
four-layer distinction. Full process-domain quantifiers belong to proof.md.
Default only compares; --write-once exclusively creates independent_results.json.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import argparse
import json

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir())

@dataclass(frozen=True)
class Q19:
    a: F
    b: F = F(0)
    def __add__(self, other):
        other = other if isinstance(other, Q19) else Q19(F(other))
        return Q19(self.a + other.a, self.b + other.b)
    __radd__ = __add__
    def __neg__(self): return Q19(-self.a, -self.b)
    def __sub__(self, other): return self + (-other if isinstance(other,Q19) else -F(other))
    def __mul__(self, other):
        other = other if isinstance(other, Q19) else Q19(F(other))
        return Q19(self.a*other.a+19*self.b*other.b,
                   self.a*other.b+self.b*other.a)
    __rmul__ = __mul__
    def as_json(self): return {'rational':str(self.a), 'sqrt19':str(self.b)}

ETA = [F(1), F(-1), F(-1), F(-1)]
IDENTITY = [[F(i == j) for j in range(4)] for i in range(4)]

def matmul(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def apply(a,v): return [sum(a[i][j]*v[j] for j in range(4)) for i in range(4)]

def minkowski(a,b): return sum(ETA[i]*a[i]*b[i] for i in range(4))

def boost(axis,v,gamma):
    assert gamma*gamma*(1-v*v)==1
    b=[row[:] for row in IDENTITY]
    b[0][0]=b[axis+1][axis+1]=gamma
    b[0][axis+1]=b[axis+1][0]=gamma*v
    return b

def calculate():
    bx=boost(0,F(3,5),F(5,4))
    by=boost(1,F(-4,5),F(5,3))
    bz=boost(2,F(12,13),F(13,5))
    frames=[IDENTITY,bx,by,bz,matmul(by,bx),matmul(bz,matmul(by,bx))]
    lorentz_checks=0
    shell_checks=0
    phase_checks=0
    tau=[F(1,7),F(2,9),F(5,11)]
    assert matmul(bx,by)!=matmul(by,bx)
    for frame in frames:
        for i in range(4):
            for j in range(4):
                assert sum(ETA[k]*frame[k][i]*frame[k][j] for k in range(4)) == (ETA[i] if i==j else 0)
                lorentz_checks+=1
        u=apply(frame,[F(1),F(0),F(0),F(0)])
        assert minkowski(u,u)==1 and u[0]>0
        for level in range(4):
            old=Q19(F(15+level),F(2))
            new=Q19(F(17+level),F(2))
            pold=[old*x for x in u]
            pnew=[new*x for x in u]
            assert minkowski(pnew,pnew)==new*new
            assert minkowski(pold,pold)==old*old
            assert all(pnew[i]-pold[i]==Q19(2*u[i]) for i in range(4))
            shell_checks+=1
            old_action=Q19(F(0)); new_action=Q19(F(0))
            for duration in tau:
                delta=[duration*x for x in u]
                assert minkowski(pnew,delta)==new*duration
                assert minkowski(pold,delta)==old*duration
                old_action+=minkowski(pold,delta)
                new_action+=minkowski(pnew,delta)
                phase_checks+=1
            assert new_action-old_action==Q19(2*sum(tau))
        # Every full supported preparation includes the added inventory energy.
        old_total=Q19(F(15)+F(92,3),F(2))
        new_total=old_total+2
        assert all(apply(frame,[new_total,Q19(0),Q19(0),Q19(0)])[i]
                   -apply(frame,[old_total,Q19(0),Q19(0),Q19(0)])[i]
                   ==Q19(2*u[i]) for i in range(4))
    # Finite consecutive legs with distinct velocities, then one common boost.
    velocities=[apply(bx,[F(1),F(0),F(0),F(0)]),
                apply(by,[F(1),F(0),F(0),F(0)]),
                apply(bz,[F(1),F(0),F(0),F(0)])]
    baseline_action=F(0)
    for u,duration in zip(velocities,tau):
        delta=[duration*x for x in u]
        p=[2*x for x in u]
        assert minkowski(apply(frames[-1],p),apply(frames[-1],delta))==2*duration
        baseline_action+=minkowski(p,delta)
    assert baseline_action==2*sum(tau)
    # Algebraic branch-phase deletion only, not a licensed interferometer.
    path_tau_pi=(F(1,8),F(3,8))
    relative_phase_pi=2*(path_tau_pi[1]-path_tau_pi[0])
    assert relative_phase_pi==F(1,2)  # exp(-i*pi/2) = -i
    phase_re,phase_im=F(0),F(-1)
    overlap_squared=((1+phase_re)**2+phase_im**2)/4
    assert overlap_squared==F(1,2)
    # Exact original four-layer obstruction: q_mid = u/(1+u^2).
    weights=[F(1,8),F(3,8),F(3,8),F(1,8)]
    assert sum(weights)==1 and sum(weights[1:3])==F(3,4)
    # (1+u^2)-2u = (u-1)^2; positive denominator for every u>0.
    gap_polynomial=[F(1),F(0)-2,F(1)]
    linear=[F(-1),F(1)]
    square_polynomial=[sum(linear[i]*linear[j] for i in range(2) for j in range(2) if i+j==degree) for degree in range(3)]
    assert gap_polynomial==square_polynomial
    sample_diagnostics=[]
    for u in [F(1,10),F(1,3),F(1),F(3),F(10)]:
        layer=[u**j for j in range(4)]
        middle=(layer[1]+layer[2])/sum(layer)
        assert middle==u/(1+u*u)
        assert F(1,2)-middle==(u-1)**2/(2*(1+u*u))
        assert middle<=F(1,2)
        sample_diagnostics.append({'u':str(u),'gibbs_middle':str(middle)})
    inputs={
      '1219_candidate':'research_cognition_physics/archive_1219_1233/1219/candidate_proof.md',
      '1234_proof':'research_cognition_physics/archive_1234_/1234/proof.md',
      '1166_proof':'research_cognition_physics/archive_1161_1200/1166/proof.md',
    }
    bindings={name:{'path':path,'sha256':sha256((ROOT/path).read_bytes()).hexdigest()}
              for name,path in inputs.items()}
    return {
      'verdict':'PASS_INDEPENDENT_1271_FINITE_IDENTITIES_ONLY',
      'independent_groups':1,
      'source_bindings':bindings,
      'independent_code_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
      'exact_algebra':'Fraction and Q(sqrt(19)); no author-code imports',
      'lorentz_component_checks':lorentz_checks,
      'full_mass_shell_checks':shell_checks,
      'leg_action_checks':phase_checks,
      'new_complete_D_base':Q19(F(17),F(2)).as_json(),
      'new_active_initial_energy_upper':Q19(F(152,3),F(2)).as_json(),
      'finite_path_added_action':str(baseline_action),
      'algebraic_unequal_path_extra_phase':'exp(-i*pi/2)',
      'algebraic_omitted_phase_overlap_squared':str(overlap_squared),
      'original_middle_weight':'3/4',
      'all_beta_gibbs_middle_bound':'1/2',
      'resulting_trace_distance_lower':'1/4',
      'gibbs_samples_diagnostics_only':sample_diagnostics,
      'source_weights_preservation_premise':'[U, H_D tensor I]=0 on full type; verified by inherited analytic gate identities, not this finite spectrum check',
      'scope':{
        'full_domain_claims':'analytic final proof only',
        'actual_h_reference':'fixed-motion internal department; never full D Gibbs',
        'new_1219_heat_groups':0,
        'new_phase_instrument_permission':False,
        'CT_or_entire_programme_counterexample':False,
      },
    }

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args()
    result=calculate()
    target=HERE/'independent_results.json'
    if args.write_once:
        with target.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
        action='WRITE_ONCE_PASS'
    else:
        assert json.loads(target.read_text(encoding='utf-8'))==result
        action='READ_ONLY_MATCH'
    print(json.dumps({'verdict':result['verdict'],'action':action,
                      'independent_groups':1,'new_active_upper':result['new_active_initial_energy_upper']},ensure_ascii=False))

if __name__=='__main__': main()
