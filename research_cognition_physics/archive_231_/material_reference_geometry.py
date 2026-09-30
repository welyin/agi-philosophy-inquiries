"""Reproduce the classical material-reference/geometry interface; NumPy only.

All gravitational, scalar-field, dimensional and detector assumptions are inputs.
No claim of a quantum measurement bound, a microscopic limit, or dimension selection.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / "material_reference_geometry_results.json"


def source(d, a, v, b):
    rho = 0.5 * (v*v + d*b*b/(a*a))
    p = 0.5 * (v*v + (2-d)*b*b/(a*a))
    return rho, p


def rhs(y, d, b, kappa=1.0):
    a, h, v, clock, conformal, work = y
    rho, p = source(d, a, v, b)
    return np.array([a*h, -kappa*(rho+p)/(d-1), -d*h*v,
                     v, 1/a, d*h*p*a**d])


def trajectory(d, q=0.6, b=0.4, end=10.0, steps=1024):
    rho, _ = source(d, 1.0, q, b)
    y = np.array([1., math.sqrt(2*rho/(d*(d-1))), q, 0., 0., 0.])
    out = [y.copy()]
    dt = end/steps
    for _ in range(steps):
        k1 = rhs(y, d, b)
        k2 = rhs(y+dt*k1/2, d, b)
        k3 = rhs(y+dt*k2/2, d, b)
        k4 = rhs(y+dt*k3, d, b)
        y = y + dt*(k1+2*k2+2*k3+k4)/6
        out.append(y.copy())
    return np.linspace(0., end, steps+1), np.asarray(out)


def clock_parameters(d, q, b, kappa=1.0):
    assert d >= 2 and q > 0 and b > 0 and kappa > 0
    c = math.sqrt(d/(kappa*(d-1)))
    r = q/(math.sqrt(d)*b)
    return c, r, c*math.asinh(r)  # a(0)=1, T(0)=0


def clock_of_a(a, d, q, b):
    c, r, limit = clock_parameters(d, q, b)
    return limit - c*np.arcsinh(r/np.asarray(a)**(d-1))


def a_of_clock(clock, d, q, b):
    c, r, limit = clock_parameters(d, q, b)
    if np.any(np.asarray(clock) >= limit):
        raise ValueError("The proposed clock reading lies outside this expanding chart.")
    return (r/np.sinh((limit-np.asarray(clock))/c))**(1/(d-1))


def adot(a, d, q, b):
    return np.sqrt((q*q*np.asarray(a)**(2-2*d)+d*b*b)/(d*(d-1)))


def integral(f, lo, hi):
    # Independent Gauss quadrature, never uses the ODE trajectory.
    x, w = np.polynomial.legendre.leggauss(96)
    nodes = (lo+hi)/2 + (hi-lo)*x/2
    return float((hi-lo)*np.dot(w, f(nodes))/2)


def metric(a, d):
    return np.diag([-1.]+[a*a]*d)


def gradients(a, d, q, b):
    return np.diag([q/a**d]+[b]*d)


def stress_from_action(g, f):
    kinetic = np.einsum("Am,mn,An->", f, np.linalg.inv(g), f)
    return f.T @ f - g*kinetic/2


def geometry_checks():
    worst = 0.
    for d in (2,3,4,7):
        for a in (0.75,1.,2.,4.):
            q, b = 0.6, 0.4
            g = metric(a,d)
            f = gradients(a,d,q,b)
            tensor = stress_from_action(g,f)
            rho, p = source(d,a,q/a**d,b)
            expected = np.diag([rho]+[a*a*p]*d)
            worst = max(worst,float(np.max(abs(tensor-expected))))
            h2 = 2*rho/(d*(d-1))
            dh = -(rho+p)/(d-1)
            # Independently check spatial Einstein equation G_ii/a^2=p.
            assert abs(-(d-1)*dh-d*(d-1)*h2/2-p) < 1e-12
            # The scalar continuity equation is independently differentiated.
            drho = -d*math.sqrt(h2)*(q*q/a**(2*d)+b*b/a**2)
            assert abs(drho+d*math.sqrt(h2)*(rho+p)) < 1e-12
    assert worst < 1e-12
    return {"dimensions_tested":[2,3,4,7],"stress_component_max_residual":worst}


def evolution_checks():
    rows=[]
    for d in (2,3,4):
        q,b=0.6,0.4
        coarse=trajectory(d,q,b,steps=512)[1]
        time, fine=trajectory(d,q,b,steps=1024)
        a,h,v,tau,chi,work=fine.T
        rho,p=source(d,a,v,b)
        errors=dict(
            constraint=float(np.max(abs(d*(d-1)*h*h/2-rho))),
            clock_charge=float(np.max(abs(a**d*v-q))),
            pressure_work=float(np.max(abs(a**d*rho-rho[0]+work))),
            exact_clock=float(np.max(abs(tau-clock_of_a(a,d,q,b)))),
            step_comparison=float(np.max(abs(fine[::2]-coarse))))
        aa=float(a[-1])
        qt=integral(lambda x:1/adot(x,d,q,b),1,aa)
        qc=integral(lambda x:1/(x*adot(x,d,q,b)),1,aa)
        errors["independent_time_quadrature"]=abs(qt-time[-1])
        errors["independent_null_quadrature"]=abs(qc-chi[-1])
        assert max(errors.values()) < 2e-7, (d,errors)
        rows.append(dict(d=d,end_time=float(time[-1]),end_a=aa,
                         end_clock=float(tau[-1]),errors=errors))
    return rows


def relational_invariance_checks():
    rng=np.random.default_rng(521)  # seed is not a numbered-round claim
    worst_metric=worst_stress=worst_interval=0.
    cases=0
    for d in (2,3,4):
        for a in (1.,1.5,3.):
            g=metric(a,d)
            f=gradients(a,d,0.6,0.4)
            invf=np.linalg.inv(f)
            expected=invf.T@g@invf
            # y coordinates -> x coordinates; arbitrary well-conditioned local Jacobian.
            # Every such matrix is itself an invertible affine chart map.
            jac=np.eye(d+1)+0.14*rng.normal(size=(d+1,d+1))
            assert np.linalg.det(jac)>0.1
            gy=jac.T@g@jac
            fy=f@jac
            gi=np.linalg.inv(fy).T@gy@np.linalg.inv(fy)
            worst_metric=max(worst_metric,float(np.max(abs(gi-expected))))
            ty=stress_from_action(gy,fy)
            tx=stress_from_action(g,f)
            worst_stress=max(worst_stress,float(np.max(abs(ty-jac.T@tx@jac))))
            # A change of reference reporting convention is tensorial, not
            # equality of component tables. It is not a new canonical field solution.
            chart=np.eye(d+1)
            chart[0,0]=1.2
            chart[1:,0]=rng.normal(size=d)/5
            chart[1:,1:]+=np.diag(np.linspace(0.1,0.3,d))
            gin=np.linalg.inv(chart).T@expected@np.linalg.inv(chart)
            dx=rng.normal(size=d+1)
            worst_interval=max(worst_interval,abs(float(dx@expected@dx-
                                          (chart@dx)@gin@(chart@dx))))
            cases+=1
    assert max(worst_metric,worst_stress,worst_interval)<1e-9
    return dict(cases=cases,metric_residual=worst_metric,
                stress_covariance_residual=worst_stress,
                reporting_chart_interval_residual=worst_interval)


def clock_domain_checks():
    rows=[]
    for d in (2,3,4,7):
        q,b=0.6,0.4
        c,r,limit=clock_parameters(d,q,b)
        av=np.array([1.,1.2,2.,5.])
        tv=clock_of_a(av,d,q,b)
        inverse_error=float(np.max(abs(a_of_clock(tv,d,q,b)-av)))
        assert inverse_error<1e-10
        # Tail is evaluated without subtracting two almost equal finite values.
        tail=lambda a:c*np.arcsinh(r/a**(d-1))
        ratio=float(tail(20.)/tail(10.))
        assert abs(ratio/2**(-(d-1))-1)<0.002
        failed=False
        try:
            a_of_clock(limit,d,q,b)
        except ValueError:
            failed=True
        assert failed
        rows.append(dict(d=d,finite_clock_limit=limit,inverse_error=inverse_error,
                         tail_ratio_20_to_10=ratio,asymptotic_ratio=2**(-(d-1))))
    return rows


def finite_readout_checks():
    # The uncertainty sets are deliberately assumed, not generated by a detector.
    # Positive full uncertainty radii are used, not a mere upper bound on
    # the noise that happened on a particular shot.
    d,q,b,a=3,0.6,0.4,2.
    eps_t,eps_s=0.01,0.02
    clock=float(clock_of_a(a,d,q,b))
    a_plus=float(a_of_clock(clock+eps_t,d,q,b))
    dt=integral(lambda x:1/adot(x,d,q,b),a,a_plus)
    differential=eps_t*a**d/q
    assert dt>differential
    target_dt=0.20
    target_l=0.12
    necessary_density=eps_t**2/(2*target_dt**2)+d*eps_s**2/(2*target_l**2)
    rho,_=source(d,a,q/a**d,b)
    actual_l=a*eps_s/b
    assert actual_l <= target_l and dt <= target_dt
    assert rho>=necessary_density
    # These readings cannot correspond to any future event in this chart.
    limit=clock_parameters(d,q,b)[2]
    near_end_clock=limit-eps_t/2
    assert near_end_clock+eps_t>limit
    # Smooth value relabelling must rescale its uncertainty as well.
    # If \tilde T = 5T, both eps and the gradient multiply by 5.
    assert abs((5*eps_t)/(5*q/a**d)-differential)<1e-14
    return dict(a=a,clock_error_radius=eps_t,rod_vector_error_radius=eps_s,
                exact_future_time_ambiguity=dt,
                local_linear_time_ambiguity=differential,
                exact_same_slice_distance_ambiguity=actual_l,
                stipulated_time_tolerance=target_dt,stipulated_distance_tolerance=target_l,
                conditional_necessary_reference_density=necessary_density,
                actual_reference_density=rho,
                near_clock_endpoint_full_error_interval_invalid=True,
                bound_is_a_quantum_measurement_lower_bound=False)


def propagation_and_failure_checks():
    d,a,q,b=3,1.7,0.6,0.4
    g=metric(a,d)
    f=gradients(a,d,q,b)
    gx=np.linalg.inv(f).T@g@np.linalg.inv(f)
    tangent=np.zeros(d+1)
    tangent[0]=1
    tangent[1]=1/a
    tx=f@tangent
    null1=float(tangent@g@tangent)
    null2=float(tx@gx@tx)
    assert abs(null1)<1e-12 and abs(null2)<1e-12
    speed=tx[1]/tx[0]
    assert abs(speed-b*a**(d-1)/q)<1e-12
    lost_clock=int(np.linalg.matrix_rank(gradients(a,d,0,b)))
    lost_rods=int(np.linalg.matrix_rank(gradients(a,d,q,0)))
    assert lost_clock==d and lost_rods==1
    # Holding these nonzero physical sources on a Minkowski metric fails Einstein.
    tensor=stress_from_action(metric(1,d),gradients(1,d,q,b))
    minkowski_residual=float(np.max(abs(tensor)))
    assert minkowski_residual>0
    return dict(null_in_spacetime_chart=null1,null_in_reference_chart=null2,
                reference_coordinate_null_speed=float(speed),
                rank_when_clock_stops=lost_clock,rank_when_rods_vanish=lost_rods,
                deleting_reference_backreaction_residual=minkowski_residual)


def run():
    checks={
        "action_stress_and_equations":geometry_checks(),
        "independent_coupled_evolution":evolution_checks(),
        "simultaneous_coordinate_changes":relational_invariance_checks(),
        "finite_clock_range":clock_domain_checks(),
        "finite_readout_contract":finite_readout_checks(),
        "null_propagation_and_failure_controls":propagation_and_failure_checks()}
    return dict(date="2026-09-30",diagnostic_tests=len(checks),failures=0,errors=0,
                numbered_round_created=False,numbered_test_increment=0,
                scope=dict(classical_effective_model=True,einstein_action_input=True,
                           spacetime_dimension_input=True,reference_stress_retained=True,
                           actual_detector_implementation=False,
                           original_graph_continuum_map=False,
                           spatial_dimension_or_GR_derived=False,
                           unified_framework_complete=False),
                diagnostics=checks)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    answer=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(answer,ensure_ascii=False,indent=2)+"\n")
    if args.check:
        assert json.loads(TARGET.read_text("utf8"))==json.loads(json.dumps(answer))
    print(json.dumps(answer,ensure_ascii=False,indent=2))
