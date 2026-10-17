"""Independent 1252 checks.

Exact rational chain rules test the complete source and initial-constraint
tangent.  A separately written constraint-reduced Decimal evolution uses a,
not the author's independent H evolution.  Step convergence is diagnostic,
never a substitute for the author's interval certificate.
Default execution only compares the saved JSON; --write exclusively creates it.
"""
from __future__ import annotations

import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULT = HERE / 'independent_results.json'


class Dual:
    def __init__(self, value, tangent=0):
        self.v, self.d = value, tangent

    @staticmethod
    def wrap(x):
        return x if isinstance(x, Dual) else Dual(x)

    def __add__(self, x):
        x = self.wrap(x)
        return Dual(self.v + x.v, self.d + x.d)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __sub__(self, x):
        return self + (-self.wrap(x))

    def __rsub__(self, x):
        return self.wrap(x) - self

    def __mul__(self, x):
        x = self.wrap(x)
        return Dual(self.v*x.v, self.d*x.v+self.v*x.d)

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = self.wrap(x)
        return Dual(self.v/x.v, (self.d*x.v-self.v*x.d)/(x.v*x.v))

    def __rtruediv__(self, x):
        return self.wrap(x)/self

    def __pow__(self, n):
        if n < 0:
            return 1/(self**(-n))
        if n == 0:
            return Dual(self.v*0+1, self.d*0)
        return Dual(self.v**n, n*self.v**(n-1)*self.d)

    def sqrt(self):
        root = self.v.sqrt()
        return Dual(root, self.d/(2*root))


def pulse(t):
    return 64*t**3*(1-t)**3, 192*t**2*(1-t)**2*(1-2*t)


def exact_source_row(phase, s):
    # Independent local state and its arbitrary tangent, not an ODE trajectory.
    raw = ('4/5', '7/10', '-1/50', '-1/40', '1/10', '2/5',
           '-1/5', '1/10', '-1/30', '1/200', '-1/300')
    states = [Dual(F(x), F(i-4, 170)) for i, x in enumerate(raw)]
    u,H,p,v,J,q,r,X,V,Y,K = states
    kap = Dual(F(6001,2000), F(1))
    w,wp = pulse(s) if phase == 0 else (F(0),F(0))
    f,fp = pulse(s) if phase == 1 else (F(0),F(0))
    mass = F(71,100)+J+(q*q+r*r+K*K)/2+w*q*p/4+f*X*Y/10
    rod = F(3,20)*(1+p/5)*u*u
    dust = (mass+F(1,50)+p/100)*u**3
    rho = (v*v+p*p+V*V+X*X)/2+rod+dust
    pressure = (v*v-p*p+V*V-X*X)/2-rod/3
    total = rho+pressure
    up,Hp = -H*u, -kap*total/2
    vp = -3*H*v-p-F(3,100)*u*u-(F(1,100)+w*q/4)*u**3
    Vp = -3*H*V-X-f*Y*u**3/10
    Jp = -wp*q*p/4-fp*X*Y/10
    qp,rp,Yp,Kp = r,-q-w*p/4,K,-f*X/10
    massp = Jp+q*qp+r*rp+K*Kp+(wp*q*p+w*qp*p+w*q*v)/4
    massp += (fp*X*Y+f*V*Y+f*X*Yp)/10
    exchange = massp-w*q*v/4-f*Y*V/10
    rodp = F(3,100)*v*u*u-2*H*rod
    dustp = (massp+v/100)*u**3-3*H*dust
    rhop = v*vp+p*v+V*Vp+X*V+rodp+dustp
    continuity = rhop+3*H*total
    constraint_rate = 2*H*Hp-kap*rhop/3
    for item in (exchange,continuity,constraint_rate):
        assert item.v == 0 and item.d == 0
    assert Hp.d == -total.v/2-kap.v*total.d/2
    # Deleting the program's actual energy reaction breaks the same-source law.
    omit_control = (wp*q*p/4+fp*X*Y/10)*u**3
    assert omit_control.v != 0
    # Deleting only the declared positive memory baseline changes actual source.
    omitted_rho = F(1,100)*u**3
    assert omitted_rho.v > 0
    return dict(phase=phase,local_pulse_time=str(s),
                complete_exchange='0',continuity_and_tangent='0',
                constraint_rate_and_tangent='0',
                omitted_control_continuity=str(omit_control.v),
                omitted_memory_source=str(omitted_rho.v),
                geometric_residual_if_memory_source_omitted=str(kap.v*omitted_rho.v))


def exact_group():
    rows=[exact_source_row(phase,F(1,4)) for phase in (0,1)]
    rho0=F(111,100)
    # 2 H0 H0_kappa=rho0/3; no irrational arithmetic is required.
    correct_constraint_tangent=rho0/3-rho0/3
    wrong_fixed_H0_tangent=-rho0/3
    assert correct_constraint_tangent==0 and wrong_fixed_H0_tangent==F(-37,100)
    m,eta=F(1,20000),F(1,10**10)
    dk=F(1,1000)
    assert m*dk-2*eta==F(249,5_000_000_000)>0
    component=rho0*eta/m
    assert component==F(111,50_000_000)
    pred_error=F(1,10**19)
    accepted_baseline_parameter_bound=(2*eta+pred_error)/m
    old_result=json.loads((HERE.parent/'1247'/'results.json').read_text(encoding='utf-8'))
    old_K=old_result['group_finite_memory']['alpha1']['checkpoints']['2']['state']['K']
    assert (F(old_K[1])-F(old_K[0]))/2<=pred_error
    # Positive homogeneous rods/dust and canonical scalar kinetic/potential
    # imply |P|<=rho; this uses the author's uniform positive-domain proof.
    return dict(status='PASS',method='exact rational independent chain rules',
                source_rows=rows,initial_constraint_tangent='0',
                fixed_H0_deletion_tangent=str(wrong_fixed_H0_tangent),
                inverse_budget=dict(required_sensitivity_lower=str(m),
                    terminal_read_error=str(eta),parameter_interval_width=str(dk),
                    predicted_endpoint_separation_lower=str(m*dk),
                    endpoint_two_read_margin=str(m*dk-2*eta),
                    single_inverse_error=str(eta/m),
                    feasible_parameter_set_diameter=str(2*eta/m),
                    orthonormal_T_component_bound=str(rho0),
                    single_inverse_residual_component_bound=str(component),
                    feasible_set_residual_component_diameter=str(2*component),
                    baseline_prediction_error=str(pred_error),
                    old_K_interval=old_K,
                    baseline_acceptance_rule='abs(observed_K-old_K_midpoint)<=eta_K',
                    accepted_baseline_parameter_bound=str(accepted_baseline_parameter_bound),
                    accepted_baseline_field_bound=str(rho0*accepted_baseline_parameter_bound)),
                residual_scope='Same actual Einstein-family solution; normalized L1 tensor tests or orthonormal components. Not arbitrary off-family equations.',
                assumptions='Complete known preparation, fixed units and other parameters, Lambda=0, actual kappa in the declared interval, inherited terminal numeric read contract.')


def decimal_wave(t,phase):
    zero=D(0)
    if phase==0:
        w,wp=pulse(t)
        return w,wp,zero,zero
    f,fp=pulse(t-1)
    return zero,zero,f,fp


def decimal_energy(t,y,phase):
    a,p,v,J,q,r,z,L,psi,X,V,Y,K=y
    w,wp,f,fp=decimal_wave(t,phase)
    mass=D('.71')+J+(q*q+r*r+K*K)/2+w*q*p/4+f*X*Y/10
    rods=D('.15')*(1+p/5)/a**2
    dust=(mass+D('.02')+p/100)/a**3
    rho=(v*v+p*p+V*V+X*X)/2+rods+dust
    pressure=(v*v-p*p+V*V-X*X)/2-rods/3
    return mass,rods,dust,rho,pressure


def decimal_rhs(t,y,kap,phase):
    a,p,v,J,q,r,z,L,psi,X,V,Y,K=y
    w,wp,f,fp=decimal_wave(t,phase)
    mass,rods,dust,rho,pressure=decimal_energy(t,y,phase)
    H=(kap*rho/3).sqrt()
    return [a*H,v,-3*H*v-p-D('.03')/a**2-(D('.01')+w*q/4)/a**3,
            -wp*q*p/4-fp*X*Y/10,r,-q-w*p/4,
            p/10,1/a,w*q*p/4,V,-3*H*V-X-f*Y/(10*a**3),K,-f*X/10]


def integrate(kappa,steps):
    vals=('1','0','0','.1','.5','0','0','0','0','.1','0','0','0')
    y=[Dual(D(v),D(0)) for v in vals]
    kap=Dual(D(kappa),D(1))
    dt=D(1)/steps
    minima={name:D(10) for name in ('J','mass','F','Mr')}
    max_rho=D(0)
    for phase in (0,1):
        for n in range(steps):
            t=D(phase)+D(n)/steps
            k1=decimal_rhs(t,y,kap,phase)
            k2=decimal_rhs(t+dt/2,[a+dt*b/2 for a,b in zip(y,k1)],kap,phase)
            k3=decimal_rhs(t+dt/2,[a+dt*b/2 for a,b in zip(y,k2)],kap,phase)
            k4=decimal_rhs(t+dt,[a+dt*b for a,b in zip(y,k3)],kap,phase)
            mass,rod,dust,rho,pressure=decimal_energy(t,y,phase)
            assert rho.v > 0
            assert abs(pressure.v)<=rho.v
            max_rho=max(max_rho,rho.v)
            for key,value in {'J':y[3].v,'mass':mass.v,'F':1+y[1].v/5,
                              'Mr':D('.2')+y[1].v/10}.items():
                minima[key]=min(minima[key],value)
            y=[a+dt*(b+2*c+2*d+e)/6 for a,b,c,d,e in zip(y,k1,k2,k3,k4)]
    assert min(minima.values())>0
    assert max_rho<=D('1.11')
    mass,rod,dust,rho,pressure=decimal_energy(D(2),y,1)
    return dict(kappa=str(kap.v),steps_per_pulse=steps,K=str(y[-1].v),
                K_kappa=str(y[-1].d),a=str(y[0].v),rho=str(rho.v),
                H=str((kap*rho/3).sqrt().v),minima={k:str(v) for k,v in minima.items()},
                max_sampled_rho=str(max_rho))


def diagnostic_group():
    rows=[]
    with localcontext() as ctx:
        ctx.prec=45
        for kap in ('3','3.0005','3.001'):
            coarse=integrate(kap,256)
            fine=integrate(kap,512)
            dK=abs(D(coarse['K'])-D(fine['K']))
            ds=abs(D(coarse['K_kappa'])-D(fine['K_kappa']))
            assert dK<D('1e-11') and ds<D('1e-11')
            assert -D(fine['K_kappa'])>D('0.00005')
            rows.append(dict(coarse=coarse,fine=fine,
                             K_step_difference=str(dK),tangent_step_difference=str(ds)))
        endpoint_gap=D(rows[0]['fine']['K'])-D(rows[-1]['fine']['K'])
        assert endpoint_gap>D('5e-8')
    return dict(status='PASS',method='45-digit Decimal dual-number RK4 in a, with H reconstructed from the constraint',
                rows=rows,endpoint_gap_diagnostic=str(endpoint_gap),
                scope='Finite diagnostic points and step agreement only. Uniform parameter monotonicity and rigorous truncation bounds require the independently audited author interval proof.')


def bindings():
    paths=[HERE/'admission.md',HERE/'proof.md',HERE/'check.py',HERE/'results.json',
           HERE/'NEXT.md',HERE.parent/'research_note_1252.md',HERE.parent/'1247'/'proof.md',
           HERE.parent/'1247'/'results.json']
    return {str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths if p.exists()}


def author_certificate_metadata():
    path=HERE/'results.json'
    if not path.exists():
        return {'status':'PENDING_AUTHOR_CERTIFICATE'}
    data=json.loads(path.read_text(encoding='utf-8'))
    group=data['group_parameter_family']
    cells=group['cells']
    assert group['kappa_interval']==['3','3001/1000']
    assert F(cells[0]['kappa'][0])==3
    assert F(cells[-1]['kappa'][1])==F(3001,1000)
    for left,right in zip(cells,cells[1:]):
        assert F(left['kappa'][1])==F(right['kappa'][0])
    for cell in cells:
        assert F(cell['K_derivative'][1]) < -F(1,20000)
        assert cell['pulse_junctions_exact']
        assert cell['constraint_proof_is_analytic_continuity_not_zero_enclosure']
        for interval in cell['positive_all_time'].values():
            assert F(interval[0])>0
    assert data['author_code_sha256']==hashlib.sha256((HERE/'check.py').read_bytes()).hexdigest()
    inverse=data['group_inverse_and_sources']
    assert F(inverse['reference_prediction_error'])==F(1,10**19)
    assert F(inverse['parameter_error_bound'])==F(2000000001,500000000000000)
    assert F(inverse['orthonormal_field_residual_bound'])==F(222000000111,50000000000000000)
    return dict(status='BOUND_AUTHOR_CERTIFICATE_METADATA_CHECKED',
                cells=len(cells),uniform_K_derivative=group['uniform_K_derivative'],
                scope='Exact stored-bound/coverage check; full interval-code correctness is separately reviewed and its default run replayed. This is not an independent second interval implementation.')


def calculate():
    return dict(round=1252,scientific_groups=2,status='PASS_WITH_SCOPE',
                exact_source_and_inverse=exact_group(),
                independent_constraint_reduced_evolution=diagnostic_group(),
                author_certificate=author_certificate_metadata(),
                source_bindings=bindings())


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--calculate-only',action='store_true')
    args=parser.parse_args()
    out=calculate()
    if args.calculate_only:
        print(json.dumps(out,ensure_ascii=False,indent=2))
    elif args.write:
        assert out['author_certificate']['status']=='BOUND_AUTHOR_CERTIFICATE_METADATA_CHECKED'
        assert len(out['source_bindings'])==8, 'Final author/report files must exist before saving.'
        with RESULT.open('x',encoding='utf-8') as f:
            json.dump(out,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('INITIAL_SAVE_PASS_WITH_SCOPE')
    else:
        assert json.loads(RESULT.read_text(encoding='utf-8'))==out
        print('READ_ONLY_MATCH')


if __name__=='__main__':
    main()
