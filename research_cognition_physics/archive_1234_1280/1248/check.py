"""1248: finite common-history samples -> a weak acceleration residual.

The scientific certificate is rational: frozen 1247 endpoint inclusions,
holding-window derivative bounds, exact linear-hat weights and error budgets.
The short RK4 trajectory is only a calibration reference, never measured data.
No integration over the already certified interval [0, 2] is repeated.
Default execution compares with the saved result without modifying any file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
OLD_RESULT = PARENT / '1247' / 'results.json'
EXPECTED_OLD_SHA = 'dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06'
EXPECTED_OLD_CODE_SHA = 'cb34c985dac837b7cbd450ae7a218cab2f4b38596127037c021453ec654ea1b8'
L, N = F(1, 4), 128
STEP = L / N
ETA_A = ETA_RHO = ETA_P = F(1, 10**8)
ETA_K = F(1, 10**10)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact(x):
    """JSON-friendly exact number, with a convenience decimal display."""
    return dict(rational=str(x), decimal=float(x))


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def derivative(a):
    return [i*a[i] for i in range(1, len(a))]


def integral(a, left, right):
    return sum((c*(right**(i+1)-left**(i+1))/(i+1)
                for i, c in enumerate(a)), F(0))


def weights():
    # In s coordinates phi dt = f(s) ds, phi'' dt = f''(s) ds/L**2.
    f = [F(0), F(0), F(30), F(-60), F(30)]
    f2 = [c/L**2 for c in derivative(derivative(f))]
    wa, wb = [F(0)]*(N+1), [F(0)]*(N+1)
    for j in range(N):
        lo, hi = F(j, N), F(j+1, N)
        left_hat, right_hat = [N*hi, -F(N)], [-N*lo, F(N)]
        for index, hat in [(j, left_hat), (j+1, right_hat)]:
            wa[index] += integral(multiply(hat, f2), lo, hi)
            wb[index] += integral(multiply(hat, f), lo, hi)
    times = [2+j*STEP for j in range(N+1)]
    assert sum(wa) == 0 and sum(wb) == 1
    assert sum(x*t for x, t in zip(wa, times)) == 0
    assert sum(x*t for x, t in zip(wb, times)) == F(17, 8)
    assert all(x >= 0 for x in wb)
    assert wa == wa[::-1] and wb == wb[::-1]
    assert sum(map(abs, wa)) < 384
    # ||phi''||_1=640/sqrt(3); this exact rational square test proves <384.
    assert F(640)**2 < 3*F(384)**2
    return wa, wb, dict(status='PASS', window=['2', '9/4'], N=N,
        cells=N, nodes=N+1, h=str(STEP),
        phi_norm='1', phi_second_L1_exact='640/sqrt(3)',
        phi_second_L1_upper='384',
        a_weights=[str(x) for x in wa], b_weights=[str(x) for x in wb],
        a_weight_sum='0', b_weight_sum='1',
        a_linear_moment='0', b_linear_moment='17/8',
        a_weight_absolute_sum=str(sum(map(abs, wa))),
        weights_are_exact_not_quadrature_of_sampled_phi=True)


def load_endpoints():
    assert sha(OLD_RESULT) == EXPECTED_OLD_SHA, 'frozen 1247 results changed'
    assert sha(PARENT/'1247'/'check.py') == EXPECTED_OLD_CODE_SHA
    saved = json.loads(OLD_RESULT.read_text('utf-8'))
    assert saved['status'] == 'PASS' and saved['science_groups'] == 2
    out = []
    for alpha in [0, 1]:
        item = saved['group_finite_memory'][f'alpha{alpha}']['checkpoints']['2']
        state = {k: tuple(map(F, v)) for k, v in item['state'].items()}
        energy = {k: tuple(map(F, v)) for k, v in item['energy'].items()}
        assert F(187, 500) < state['u'][0] <= state['u'][1] < F(3, 8)
        assert 0 < state['H'][0] <= state['H'][1] < F(27, 100)
        assert 0 < energy['rho'][0] <= energy['rho'][1] < F(3, 40)
        assert F(9, 10) < energy['M_c'][0] <= energy['M_c'][1] < F(117, 125)
        assert 0 < state['J'][0]
        assert max(map(abs, state['K'])) < F(3, 1000)
        assert F(1, 100) < energy['memory_complete'][0]
        assert energy['memory_complete'][1] < F(11, 1000)
        out.append((state, energy))
    return out


def derivative_and_error_certificate():
    # DEC and the already proved post-pulse continuity/positive-domain bootstrap
    # transport these endpoint bounds; finite sampling does not prove them.
    RHO, H, U, A, Z = F(3,40), F(27,100), F(3,8), F(3), F(2,5)
    HP, UP, UPP = F(9,40), F(11,100), F(3,25)
    VP, VPP = F(3,4), F(13,10)
    lower_u = F(187,500)*(1-H*L)
    assert lower_u > F(1,3) and 2*RHO < Z**2
    assert H*U < UP
    assert (HP+H**2)*U < UPP
    assert 3*H*Z+Z+F(3,100)*U**2+F(1,100)*U**3 < VP
    vpp_raw = 3*HP*Z+3*H*VP+Z+F(6,100)*U*UP+F(3,100)*U**2*UP
    assert vpp_raw < VPP

    # A second chain certifies the intermediate bounds displayed in proof.md.
    # All terms use absolute values, not cancellations of the field equation.
    FF = F(27,25)
    rod1 = F(3,100)*Z*U**2+F(3,10)*FF*U*UP
    rod2 = (F(3,100)*VP*U**2+F(3,25)*Z*U*UP
            +F(3,10)*FF*(UP**2+U*UPP))
    assert rod1 < F(2,125) and rod2 < F(3,100)
    pressure1 = 2*Z*VP+2*Z**2+F(2,125)/3
    pressure2 = 2*VP**2+2*Z*VPP+2*Z**2+2*Z*VP+F(3,100)/3
    assert pressure1 < 1 and pressure2 < F(16,5)
    rho1 = 3*H*2*RHO
    rho2 = 3*HP*2*RHO+3*H*(F(1,8)+1)
    assert rho1 < F(1,8) and rho2 < F(11,10)

    # S=rho+3P=2(v^2+V^2)-(p^2+X^2)+E*u^3.  The rod term
    # cancels algebraically here; the full rod source remains in the ODE.
    E, EP, EPP = F(24,25), F(1,100)*Z, F(1,100)*VP
    assert F(117,125)+F(1,50)+F(1,100)*Z == E
    d1 = EP*U**3+3*E*U**2*UP
    d2 = EPP*U**3+6*EP*U**2*UP+6*E*U*UP**2+3*E*U**2*UPP
    S0, S1, S2 = F(3,10), F(25,8), F(11)
    s1_raw = 8*Z*VP+4*Z**2+d1
    s2_raw = 8*(VP**2+Z*VPP)+4*(Z**2+Z*VP)+d2
    assert 4*RHO == S0 and s1_raw < S1 and s2_raw < S2
    AP, BA, BB, BB_MINUS = F(81,100), F(1), F(20), F(21)
    a2_raw = A*(HP+H**2)
    b2_raw = (BA*S0+2*AP*S1+A*S2)/2
    assert A*H == AP and a2_raw < BA and b2_raw < BB
    # The removed mass is constant after t=2.  Its interpolation requires its
    # own derivative budget, instead of reusing BB for a different function.
    mk_upper = F(11,1000)
    removed_b2 = mk_upper*(UP**2+U*UPP)
    assert removed_b2 < F(1,1000) and BB+removed_b2 < BB_MINUS

    eta_s = ETA_RHO+3*ETA_P
    eta_b = (A*eta_s+S0*ETA_A+ETA_A*eta_s)/2
    # Subtraction d(a,K)=m_K/(2a^2) is evaluated on the same recorded a and
    # retained K.  The whole connecting segment has a>5/2 and |K|<.003+eta_K.
    measured_a_lower = F(5,2)
    measured_k_upper = F(3,1000)+ETA_K
    measured_m_upper = F(1,100)+measured_k_upper**2/2
    assert F(8,3)-ETA_A > measured_a_lower
    assert measured_m_upper < mk_upper
    subtraction_error = (mk_upper/measured_a_lower**3*ETA_A
                         +measured_k_upper/(2*measured_a_lower**2)*ETA_K)
    assert subtraction_error < F(1,10**10)
    eta_b_minus = eta_b+F(1,10**10)
    interp = STEP**2*(BA*384+BB)/8
    interp_minus = STEP**2*(BA*384+BB_MINUS)/8
    budget = interp+384*ETA_A+eta_b
    budget_minus = interp_minus+384*ETA_A+eta_b_minus
    deletion_lower = F(1,1800)
    margin = deletion_lower-budget-budget_minus
    assert margin > F(1619,10**7) > 0
    return dict(status='PASS', bounds={k: exact(v) for k, v in dict(
        rho=RHO,H=H,u_upper=U,u_lower=lower_u,a_upper=A,
        scalar_absolute=Z,H_prime=HP,u_prime=UP,u_second=UPP,
        scalar_velocity_prime=VP,scalar_velocity_second=VPP,
        raw_scalar_velocity_second=vpp_raw,
        raw_rod_prime=rod1,raw_rod_second=rod2,
        raw_pressure_prime=pressure1,raw_pressure_second=pressure2,
        raw_rho_prime=rho1,raw_rho_second=rho2,
        dust_prime=d1,dust_second=d2,
        S=S0,S_prime=S1,S_second=S2,raw_S_prime=s1_raw,raw_S_second=s2_raw,
        a_prime=AP,raw_a_second=a2_raw,raw_b_second=b2_raw,
        B_a=BA,B_b=BB,B_b_omitted=BB_MINUS,removed_b_second=removed_b2).items()},
        errors={k: exact(v) for k, v in dict(eta_a=ETA_A,eta_rho=ETA_RHO,
        eta_pressure=ETA_P,eta_K=ETA_K,eta_b=eta_b,
        subtraction_error_bound=subtraction_error,eta_b_omitted=eta_b_minus,
        interpolation=interp,interpolation_omitted=interp_minus,
        total_full=budget,total_omitted=budget_minus,
        memory_deletion_lower=deletion_lower,separation_margin=margin).items()},
        full_and_omitted_have_distinct_error_budgets=True,
        derivative_bounds_are_analytic_not_grid_maxima=True,
        W_zero_not_used_to_compute_discrete_estimator=True,
        one_test_does_not_establish_all_field_components=True)


def diagnostic(endpoints, wa, wb):
    """Short holding-window calibration, entirely separate from certification."""
    import numpy as np
    result = {}
    for alpha, (state, energy) in enumerate(endpoints):
        mid = lambda x: float((x[0]+x[1])/2)
        y = np.array([mid(state[k]) for k in ['u','H','p','v','X','V']],dtype=float)
        mass, kval = mid(energy['M_c']), mid(state['K'])
        memory = .01+kval*kval/2
        def source(z):
            u,H,p,v,X,V = z
            rod=.15*(1+.2*p)*u*u
            dust=(mass+.02+.01*p)*u**3
            rho=.5*(v*v+p*p+V*V+X*X)+rod+dust
            pressure=.5*(v*v-p*p+V*V-X*X)-rod/3
            return rho,pressure,rod,dust
        def rhs(z):
            u,H,p,v,X,V=z
            rho,pressure,rod,dust=source(z)
            return np.array([-H*u,-1.5*(v*v+V*V+2*rod/3+dust),v,
                -3*H*v-p-.03*u*u-.01*u**3,V,-3*H*V-X])
        aa=[];bb=[];bm=[];rr=[];pp=[];states=[]
        dt=float(STEP)/8
        for j in range(N+1):
            rho,pressure,rod,dust=source(y)
            a=1/y[0];b=.5*a*(rho+3*pressure)
            aa.append(a);bb.append(b);bm.append(b-.5*memory*y[0]**2)
            rr.append(rho);pp.append(pressure);states.append(y.copy())
            if j==N: break
            for _ in range(8):
                k1=rhs(y);k2=rhs(y+dt*k1/2);k3=rhs(y+dt*k2/2);k4=rhs(y+dt*k3)
                y=y+dt*(k1+2*k2+2*k3+k4)/6
        # Convert finite binary calibration numbers to exact Fractions before
        # applying the exact weights: this is not a proof of RK integration.
        wf=sum((w*F(float(v)) for w,v in zip(wa,aa)),F(0))
        full=wf+sum((w*F(float(v)) for w,v in zip(wb,bb)),F(0))
        omitted=wf+sum((w*F(float(v)) for w,v in zip(wb,bm)),F(0))
        defects=[abs(z[1]**2-rho) for z,rho in zip(states,rr)]
        result[f'alpha{alpha}']=dict(W_full=float(full),W_omitted=float(omitted),
            discrete_difference=float(full-omitted),max_constraint_defect=float(max(defects)),
            endpoint=dict(zip(['u','H','p','v','X','V'],map(float,y))),
            midpoint_M_c=mass,midpoint_K=kval,
            menu=dict(a=aa,rho=rr,pressure=pp))
    return dict(method='binary64 RK4, 8 substeps per fixed cell on [2,9/4] only',
        data_are_model_calibration_not_actual_measurements=True,
        not_an_interval_integration_certificate=True,branches=result)


def run():
    endpoints=load_endpoints()
    wa,wb,weight_group=weights()
    cert=derivative_and_error_certificate()
    source_paths={'1247/results.json':OLD_RESULT,'1247/check.py':PARENT/'1247'/'check.py',
                  '1247/proof.md':PARENT/'1247'/'proof.md','1248/check.py':Path(__file__)}
    return dict(round=1248,status='PASS',science_groups=2,
        arithmetic='Exact Fraction rational weights and budgets; frozen outward endpoint inclusions',
        group_weights_and_error=weight_group,
        group_common_domain_and_separation=cert,
        endpoint_evidence={f'alpha{a}':dict(
            u=[str(x) for x in s['u']],H=[str(x) for x in s['H']],
            rho=[str(x) for x in e['rho']],M_c=[str(x) for x in e['M_c']],
            K=[str(x) for x in s['K']],memory_mass=[str(x) for x in e['memory_complete']])
            for a,(s,e) in enumerate(endpoints)},
        diagnostic=diagnostic(endpoints,wa,wb),
        source_sha256={name:sha(path) for name,path in source_paths.items()},
        existing_0_to_2_not_reintegrated=True,new_adopted_axioms=0,whole_GR_proved=False,
        scope='Conditional classical same-history weak residual certificate. Full finite sampling and retained K precision are explicit contracts. Omitted-source accounting is not a second Einstein solution. No actual instrument precision or independent gravity generation is proved.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    output=run();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:
            json.dump(output,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(path.read_text('utf-8'))==output,'saved result mismatch'
        mode='READ_ONLY_MATCH'
    errors=output['group_common_domain_and_separation']['errors']
    print(json.dumps(dict(status='PASS',mode=mode,groups=2,
        full_budget=errors['total_full'],omitted_budget=errors['total_omitted'],
        separation=errors['separation_margin']),ensure_ascii=False))
