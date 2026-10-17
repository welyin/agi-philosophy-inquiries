"""1252: a fixed physical kappa interval, retained K, and inverse stability.

Exact outward intervals and time jets are reused read-only from 1242/1247.
Default recomputes and matches the saved result; --write creates it once.
No physical interval, pulse, or readout precision is optimized.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent
EXPECTED = {
    '1242/check.py': '2a25be6a9b371c0689f03458c0151dcd534339d983989ca8c914056749aabc16',
    '1246/check.py': 'eadd74be3258a9adbab20f5dfc7d060f50166562b29b96a5617e459e4dabe43b',
    '1247/check.py': 'cb34c985dac837b7cbd450ae7a218cab2f4b38596127037c021453ec654ea1b8',
    '1247/results.json': 'dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06',
}
for name, digest in EXPECTED.items():
    assert hashlib.sha256((BASE/name).read_bytes()).hexdigest() == digest, name
spec = importlib.util.spec_from_file_location('r1247_reuse_1252', BASE/'1247/check.py')
old = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(old)
I, Jet, S, hull = old.I, old.Jet, old.S, old.hull
NAMES = old.NAMES
KAPPA0, KAPPA1 = F(3), F(3001, 1000)
RHO0 = F(111, 100)
PARAMETER_CELLS = 1
STEPS_PER_PHASE, ORDER = 200, 8
SENSITIVITY_MIN = F(1, 20000)


class Dual:
    """Exact first parameter derivative, over either I or time-Jet values."""
    def __init__(self, v, d=0): self.v, self.d = v, d
    def coerce(self, b): return b if isinstance(b, Dual) else Dual(b, 0)
    def __add__(self, b):
        b = self.coerce(b); return Dual(self.v+b.v, self.d+b.d)
    __radd__ = __add__
    def __neg__(self): return Dual(-self.v, -self.d)
    def __sub__(self, b): return self+-self.coerce(b)
    def __rsub__(self, b): return self.coerce(b)+-self
    def __mul__(self, b):
        b = self.coerce(b)
        return Dual(self.v*b.v, self.d*b.v+self.v*b.d)
    __rmul__ = __mul__
    def __pow__(self, n):
        z = Dual(1, 0)
        for _ in range(n): z = self*z
        return z


def rhs(z, kappa, phase):
    # Each y and y_kappa follows its own actual history. Nothing is frozen.
    vals = [Dual(v, d) for v, d in zip(z[:15], z[15:])]
    out = old.rhs(vals, 1, phase)
    out[1] = out[1]*Dual(kappa, I(1))*F(1, 3)
    return ([x.v if isinstance(x, Dual) else x for x in out]
            + [x.d if isinstance(x, Dual) else 0 for x in out])


def coefficients(y, kappa, phase, order):
    c = [[I(v)] for v in y]
    for k in range(order):
        f = rhs([Jet(v) for v in c], kappa, phase)
        for j in range(len(y)):
            v = f[j].c[k] if isinstance(f[j], Jet) else I(f[j] if k == 0 else 0)
            c[j].append(v.div(k+1))
    return c


def initial_interval(klo, khi):
    y = old.initial_interval()
    lo, hi = old.sqrt_interval(klo*RHO0/3), old.sqrt_interval(khi*RHO0/3)
    y[1] = I(lo.lo, hi.hi, raw=True)
    dy = [I(0) for _ in y]
    # H_kappa = rho0/(6 H), positive reciprocal with exact rational bounds.
    dy[1] = I(RHO0*F(S, 6*y[1].hi), RHO0*F(S, 6*y[1].lo))
    return y+dy


def integrate_cell(klo, khi):
    kap = I(klo, khi)
    dt = F(1, STEPS_PER_PHASE)
    h, time = I(dt), I(0, dt)
    y = initial_interval(klo, khi)
    radii = [2*dt,3*dt,F(1,5)*dt,F(1,2)*dt,F(1,30)*dt,
             F(11,10)*dt,0,F(1,4)*dt,F(3,5)*dt,F(3,5)*dt,F(1,20)*dt,
             F(1,5)*dt,F(1,5)*dt,F(1,10)*dt,F(1,10)*dt]+[F(1,10)]*15
    radii[21] = 0
    domain = [I('.2','1.02'),I('.1','1.09'),I('-.15','.03'),I('-.15','.05'),
              I('-.03','.002'),I('-.01','2.02'),I(0,2),I('.04','.17'),
              I('-.3','.52'),I('-.52','.03'),I('-.02','.01'),
              I('-.1','.12'),I('-.15','.03'),I('-.1','.03'),I('-.1','.03')]
    hp = [I(1)]
    for _ in range(ORDER+1): hp.append(hp[-1]*h)
    all_tube = positive = None
    checkpoints = {}
    max_state_width = max_sensitivity_width = max_remainder = F(0)
    for phase in (0, 1):
        assert y[6].lo == y[6].hi == phase*S
        for step in range(STEPS_PER_PHASE):
            box = [a+I(-r,r) for a,r in zip(y,radii)]
            box[6] = I(phase+step*dt, phase+(step+1)*dt)
            assert all(a.subset(b) for a,b in zip(box[:15],domain)), ('domain',klo,khi,phase,step)
            proposed = [a+time*b for a,b in zip(y,rhs(box,kap,phase))]
            assert all(a.subset(b) for a,b in zip(proposed,box)), ('Picard',klo,khi,phase,step)
            c = coefficients(y,kap,phase,ORDER)
            tail = coefficients(box,kap,phase,ORDER+1)
            out = []
            for j in range(30):
                rem = tail[j][-1]*hp[ORDER+1]
                out.append(sum((c[j][k]*hp[k] for k in range(ORDER+1)),I(0))+rem)
                max_remainder = max(max_remainder,abs(F(rem.lo,S)),abs(F(rem.hi,S)))
            assert all(a.subset(b) for a,b in zip(out,box)), ('endpoint',klo,khi,phase,step)
            y = out
            Mc,D,rod,dust,rho,P,memory,osc = old.energy(box[:15],1,phase)
            pos = [box[0],box[1],box[7],Mc,1+F(1,5)*box[2],
                   F(1,5)+F(1,10)*box[2],memory,D+F(1,100)*box[2]]
            assert all(a.lo>0 for a in pos), ('positive',klo,khi,phase,step)
            positive = pos if positive is None else [hull(a,b) for a,b in zip(positive,pos)]
            rho_end = old.energy(y[:15],1,phase)[4]
            assert (y[1]**2-rho_end*kap*F(1,3)).contains_zero()
            erho = old.energy([Dual(v,d) for v,d in zip(y[:15],y[15:])],1,phase)[4]
            dc = 2*y[1]*y[16]-erho.v*F(1,3)-erho.d*kap*F(1,3)
            assert dc.contains_zero()
            all_tube = box if all_tube is None else [hull(a,b) for a,b in zip(all_tube,box)]
            max_state_width = max(max_state_width,*(a.width() for a in y[:15]))
            max_sensitivity_width = max(max_sensitivity_width,*(a.width() for a in y[15:]))
        checkpoints[str(phase+1)] = {
            'state': {n:a.out() for n,a in zip(NAMES,y[:15])},
            'sensitivity': {n:a.out() for n,a in zip(NAMES,y[15:])},
        }
    assert y[29].hi < -I(SENSITIVITY_MIN).hi, ('uniform derivative',klo,khi,y[29].out())
    off = rhs(y,kap,2)
    assert off[14].lo == off[14].hi == 0
    assert off[29].lo == off[29].hi == 0
    assert off[7].lo == off[7].hi == 0
    return dict(kappa=[str(klo),str(khi)],checkpoints=checkpoints,
        positive_all_time={n:a.out() for n,a in zip(
            ['u','H','J','M_c','F','M_r','memory_complete','dust_mass'],positive)},
        all_time_state_tube={n:a.out() for n,a in zip(NAMES,all_tube[:15])},
        all_time_sensitivity_tube={n:a.out() for n,a in zip(NAMES,all_tube[15:])},
        max_state_width=str(max_state_width),max_sensitivity_width=str(max_sensitivity_width),
        max_step_remainder=str(max_remainder),K_derivative=y[29].out(),
        pulse_junctions_exact=True,constraint_and_derivative_enclose_zero=True,
        constraint_proof_is_analytic_continuity_not_zero_enclosure=True,
        retained_K_and_derivative_off_pulse='0')


def source_checks():
    checks = []
    for kap in (KAPPA0,(KAPPA0+KAPPA1)/2,KAPPA1):
        for phase in (0,1,2):
            y=[F(4,5),F(7,10),F(-1,50),F(-1,40),F(-1,1000),F(1,3),
               phase+F(1,4),F(1,10),F(2,5),F(-1,5),F(0),
               F(2,25),F(-1,50),F(-1,1000),F(-1,500)]
            dy=old.rhs(y,1,phase);dy[1]*=kap/3
            u,H,p,v,z,light,theta,J,q,r,psi,X,V,Y,K=y
            up,Hp,pp,vp,zp,Lp,tp,Jp,qp,rp,psip,Xp,Vp,Yp,Kp=dy
            w,wp,f,fp=old.pulses(theta,phase)
            Mc,D,rod,dust,rho,P,memory,osc=old.energy(y,1,phase)
            Mdot=Jp+q*qp+r*rp+K*Kp+F(1,4)*(wp*q*p+w*qp*p+w*q*pp)+F(1,10)*(fp*X*Y+f*Xp*Y+f*X*Yp)
            assert Mdot==F(1,4)*w*q*v+F(1,10)*f*Y*V
            rod_dot=F(3,100)*v*u*u-2*H*rod
            dust_dot=(Mdot+F(1,100)*v)*u**3-3*H*dust
            rhodot=v*vp+p*v+V*Vp+X*V+rod_dot+dust_dot
            assert rhodot+3*H*(rho+P)==0
            assert 2*H*Hp-kap*rhodot/3==0
            checks.append(dict(kappa=str(kap),phase=phase,continuity='0',constraint_derivative='0'))
    return checks


def inference_checks():
    eta = F(1,10**10)
    pred_error = F(1,10**19)
    prior = json.loads((BASE/'1247/results.json').read_text(encoding='utf-8'))
    baseline = I(*prior['group_finite_memory']['alpha1']['checkpoints']['2']['state']['K'])
    assert baseline.width()/2 <= pred_error
    center = F(baseline.lo+baseline.hi,2*S)
    kappa_bound = (2*eta+pred_error)/SENSITIVITY_MIN
    field_bound = RHO0*kappa_bound
    endpoint_gap = SENSITIVITY_MIN*(KAPPA1-KAPPA0)
    assert endpoint_gap>2*eta and kappa_bound<KAPPA1-KAPPA0
    # Off-pulse extension: positive sources + continuity supply rho <= rho0.
    assert KAPPA1*RHO0/3 < F(11,10)**2
    assert RHO0 < F(9,8)
    assert 1-F(11,10)*F(1,4)==F(29,40)
    # Two deliberate omissions do not lie in the admitted parameter family.
    wrong_initial_constraint = -(KAPPA1-KAPPA0)*RHO0/3
    assert wrong_initial_constraint == -F(37,100000)
    wrong_initial_derivative = -RHO0/3
    omitted_memory_constraint = KAPPA1*F(1,100)/3
    assert wrong_initial_derivative != 0 and omitted_memory_constraint > 0
    return dict(status='PASS',eta_K=str(eta),reference_prediction_error=str(pred_error),
        baseline_kappa='3',baseline_K_interval=baseline.out(),baseline_midpoint=str(center),
        acceptance_rule='abs(observed_K-baseline_midpoint)<=eta_K',
        parameter_error_bound=str(kappa_bound),orthonormal_field_residual_bound=str(field_bound),
        uniform_source_component_bound=str(RHO0),
        general_inverse_rule='abs(kappa-kappa_hat) <= (read_error+prediction_error+distance_to_prediction)/m',
        m=str(SENSITIVITY_MIN),endpoint_prediction_gap_at_least=str(endpoint_gap),
        endpoint_separation_after_two_read_errors_at_least=str(endpoint_gap-2*eta),
        positive_hold_window=['2','9/4'],H_upper='11/10',u_hold_ratio_lower='29/40',
        retention_does_not_remove_initial_clock_or_gate_conditions=True,
        omitted_initial_geometric_response_constraint=str(wrong_initial_constraint),
        omitted_initial_sensitivity_constraint=str(wrong_initial_derivative),
        omitted_memory_baseline_initial_constraint=str(omitted_memory_constraint),
        only_model_inputs_fixed_not_claimed_empirical_calibration=True)


def run():
    source = source_checks()
    records=[]
    width=(KAPPA1-KAPPA0)/PARAMETER_CELLS
    for j in range(PARAMETER_CELLS):
        records.append(integrate_cell(KAPPA0+j*width,KAPPA0+(j+1)*width))
    assert records[0]['kappa'][0]==str(KAPPA0) and records[-1]['kappa'][1]==str(KAPPA1)
    for a,b in zip(records,records[1:]): assert a['kappa'][1]==b['kappa'][0]
    derivative=I(*records[0]['K_derivative'])
    for record in records[1:]: derivative=hull(derivative,I(*record['K_derivative']))
    assert derivative.hi < -I(SENSITIVITY_MIN).hi
    return dict(round=1252,status='PASS',science_groups=2,
        arithmetic='10^-36 integer outward I, exact Fraction, polynomial Dual over time Jet; no sampled parameter proof',
        group_parameter_family=dict(status='PASS',alpha=1,kappa_interval=[str(KAPPA0),str(KAPPA1)],
            parameter_cells=PARAMETER_CELLS,exact_closed_interval_union=True,
            steps_per_phase=STEPS_PER_PHASE,time_order=ORDER,remainder_order=ORDER+1,
            uniform_K_derivative=derivative.out(),m=str(SENSITIVITY_MIN),
            source_identity_diagnostics=source,cells=records),
        group_inverse_and_sources=inference_checks(),
        inputs_sha256=EXPECTED,
        author_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        new_adopted_roles=0,CT_adopted=False,
        scope=['Known full initial material state; expanding constrained initial geometry varies with kappa.',
               'Only the declared CT-selected one-parameter classical effective family, Lambda=0.',
               'Existing local write/hold and declared terminal K error, not 1249/1251 local quantum implementation.',
               'No arbitrary residual control, no derivation of the measured Newton value, no Jacobson reference claim.'])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as file:json.dump(result,file,ensure_ascii=False,indent=2);file.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(path.read_text(encoding='utf-8'))==result,'saved result mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status='PASS',mode=mode,groups=2,
        K_derivative=result['group_parameter_family']['uniform_K_derivative'],
        parameter_bound=result['group_inverse_and_sources']['parameter_error_bound'],
        field_bound=result['group_inverse_and_sources']['orthonormal_field_residual_bound']),ensure_ascii=False))
