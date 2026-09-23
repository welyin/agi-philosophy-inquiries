"""Round 294: invariant interval densities and simultaneous clock coordinates.
A smooth conditional criterion; existence of the densities is not inferred
from causality or from a finite collection of numerical samples.
"""
import unittest
import numpy as np
from growing_stream_audit import main
from local_clock_gauge_audit import clocks, rechart
from temporal_distance_audit import path_arrival


def density_fixture():
    phi,inverse=clocks()
    scales=[1.2,.8,1.5]
    rates=[.3,.4,.25]
    offsets=[2.,-3.,4.]
    prime=[lambda t,a=a,s=s:a*s*np.cosh(s*t)
           for a,s in zip(scales,rates)]
    density=[lambda x,a=a,s=s,b=b:1/(a*s*np.sqrt(1+((x-b)/a)**2))
             for a,s,b in zip(scales,rates,offsets)]
    delays={(0,1):.6,(1,0):.4,(0,2):.8,(2,0):float(np.sqrt(2)-.8)}
    arcs={(i,j):(lambda x,i=i,j=j,d=d:phi[j](inverse[i](x)+d))
          for (i,j),d in delays.items()}
    jacobian={(i,j):(lambda x,i=i,j=j,d=d:
              prime[j](inverse[i](x)+d)/prime[i](inverse[i](x)))
              for (i,j),d in delays.items()}
    return phi,inverse,density,delays,arcs,jacobian


def integrate(function,a,b):
    nodes,weights=np.polynomial.legendre.leggauss(80)
    x=(a+b)/2+(b-a)*nodes/2
    return float((b-a)/2*np.dot(weights,function(x)))


def crossing_pair(amplitude=.1):
    f=lambda x:np.asarray(x)+1.
    g=lambda x:np.asarray(x)+1.+amplitude*np.sin(2*np.pi*np.asarray(x))
    derivative=lambda x:1.+2*np.pi*amplitude*np.cos(2*np.pi*np.asarray(x))
    return f,g,derivative


def periodic_clock(x,amplitude=.05):
    return np.asarray(x)+amplitude*np.sin(2*np.pi*np.asarray(x))


def periodic_shift_error(shift,amplitude=.05):
    grid=np.linspace(0,1,20001)
    residual=periodic_clock(grid+shift,amplitude)-periodic_clock(grid,amplitude)-shift
    return float(np.max(abs(residual)))


def report():
    phi,inverse,density,delays,arcs,jacobian=density_fixture()
    residuals=[]
    recovered=[]
    for edge in arcs:
        i,j=edge
        errors=[]
        shifts=[]
        for t in np.linspace(-3,3,25):
            x=phi[i](t)
            errors.append(abs(density[j](arcs[edge](x))*jacobian[edge](x)-density[i](x)))
            shifts.append(float(inverse[j](arcs[edge](x))-inverse[i](x)))
        residuals.extend(errors)
        recovered.append({'edge':list(edge),'known_input_delay':delays[edge],
                          'recovered_shift_min':min(shifts),'recovered_shift_max':max(shifts)})
    f,g,gprime=crossing_pair()
    grid=np.linspace(-3,3,121)
    commute=float(np.max(abs(f(g(grid))-g(f(grid)))))
    return {'round':294,
            'scope':'Smooth iff criterion for constant-delay coordinates via positive complete interval densities; causal order or commutation alone is insufficient, and no physical time standard or spacetime law is derived.',
            'density_transport_positive_example':{
                'maximum_transport_residual':float(max(residuals)),
                'edges':recovered,
                'root_loop_shifts':[1.,float(np.sqrt(2))],
                'clock_functions_are_supplied_for_this_constructive_example':True},
            'commuting_but_not_simultaneously_translatable':{
                'maps':'F(x)=x+1; G(x)=x+1+0.1*sin(2*pi*x)',
                'sampled_commutation_residual':commute,
                'values_at_zero':[float(f(0.)),float(g(0.))],
                'values_at_quarter':[float(f(.25)),float(g(.25))],
                'values_at_three_quarters':[float(f(.75)),float(g(.75))],
                'return_order_reverses_between_quarter_and_three_quarters':True,
                'G_derivative_at_zero':float(gprime(0.)),
                'positive_density_ratio_required_by_F':1.,
                'positive_density_ratio_required_by_G':float(1/gprime(0.)),
                'translation_obstruction':'Maps coincide at one point but differ elsewhere; distinct translations cannot do that.',
                'minimum_G_derivative':float(1-.2*np.pi),
                'minimum_G_forward_increment':.9},
            'clock_uniqueness_boundary':{
                'nonlinear_clock':'v(u)=u+0.05*sin(2*pi*u)',
                'unit_shift_error':periodic_shift_error(1.),
                'double_shift_error':periodic_shift_error(2.),
                'sqrt2_shift_error':periodic_shift_error(float(np.sqrt(2))),
                'sqrt2_analytic_max_error':float(.1*abs(np.sin(np.pi*np.sqrt(2)))),
                'near_unit_shift_1_001_error':periodic_shift_error(1.001),
                'interpretation':'Commensurate periods leave periodic clock freedom; two exact incommensurate periods remove it under continuity.'},
            'resource_and_inference_limits':['Density transport is a condition to establish, not an automatic consequence of cognition or a fitted derivative.',
                                            'Finite samples cannot certify an identity for all times, exact irrationality, or global completeness.',
                                            'Edge constants depend on endpoint zero choices; their closed-loop sums do not.',
                                            'Formal inverse coordinate maps are not backward physical transmissions.']}


class Audit(unittest.TestCase):
    def test_01_positive_density_transport(self):
        phi,_,density,_,arcs,jac=density_fixture()
        for (i,j),h in arcs.items():
            for t in np.linspace(-5,5,31):
                x=phi[i](t)
                self.assertGreater(density[i](x),0.)
                self.assertAlmostEqual(density[j](h(x))*jac[i,j](x),density[i](x),places=12)

    def test_02_integrated_density_gives_constant_edge_shift(self):
        phi,inverse,density,delays,arcs,_=density_fixture()
        for (i,j),h in arcs.items():
            for t in (-2.,0.,3.):
                x=phi[i](t)
                # Integrals from each node's own coordinate at reference t=0.
                ui=integrate(density[i],phi[i](0.),x)
                uj=integrate(density[j],phi[j](0.),h(x))
                self.assertAlmostEqual(ui,inverse[i](x),places=11)
                self.assertAlmostEqual(uj-ui,delays[i,j],places=11)

    def test_03_analytic_derivatives_match_finite_differences(self):
        phi,_,_,_,arcs,jac=density_fixture()
        for (i,j),h in arcs.items():
            for t in (-1.,0.,2.):
                x=phi[i](t)
                step=1e-5
                numerical=(h(x+step)-h(x-step))/(2*step)
                self.assertAlmostEqual(numerical,jac[i,j](x),places=8)

    def test_04_transformed_loop_intervals_add(self):
        phi,inverse,_,_,arcs,_=density_fixture()
        for t in (-2.,0.,3.):
            for path,shift in (([0,1,0],1.),([0,2,0],np.sqrt(2)),
                               ([0,1,0,2,0],1+np.sqrt(2))):
                end=path_arrival(arcs,path,phi[0](t))
                self.assertAlmostEqual(inverse[0](end)-t,shift,places=12)

    def test_05_commuting_maps_can_cross(self):
        f,g,gprime=crossing_pair()
        grid=np.linspace(-5,5,201)
        np.testing.assert_allclose(f(g(grid)),g(f(grid)),atol=1e-12)
        self.assertTrue(np.all(gprime(grid)>0))
        self.assertTrue(np.all(g(grid)>grid))
        self.assertEqual(float(f(0.)),float(g(0.)))
        self.assertGreater(float(g(.25)),float(f(.25)))
        self.assertLess(float(g(.75)),float(f(.75)))

    def test_06_crossing_derivatives_contradict_positive_density(self):
        _,_,gprime=crossing_pair()
        # F transport at 0 forces w(1)=w(0), G transport forces another ratio.
        ratio_f=1.
        ratio_g=1/gprime(0.)
        self.assertGreater(abs(ratio_f-ratio_g),.3)
        self.assertAlmostEqual(float(gprime(0.)),1+.2*np.pi)

    def test_07_commensurate_loops_do_not_fix_subtick_scale(self):
        for shift in (1.,2.,3.):
            self.assertLess(periodic_shift_error(shift),1e-12)
        self.assertGreater(abs(float(periodic_clock(.25))-.25),.04)
        self.assertGreater(1-.1*np.pi,0.)

    def test_08_incommensurate_and_near_resonant_shift_residual(self):
        shift=float(np.sqrt(2))
        self.assertAlmostEqual(periodic_shift_error(shift),
                               .1*abs(np.sin(np.pi*shift)),places=7)
        self.assertGreater(periodic_shift_error(shift),.09)
        self.assertLess(periodic_shift_error(1.001),.000315)
        # Numerical sqrt(2) is illustrative, not a numerical proof of irrationality.

    def test_09_density_criterion_is_coordinate_covariant(self):
        phi,_,density,_,arcs,jac=density_fixture()
        psi=[lambda x:np.sinh(.2*x) for _ in range(3)]
        psi_inv=[lambda y:np.arcsinh(y)/.2 for _ in range(3)]
        changed=rechart(arcs,psi,psi_inv)
        new_density=[lambda y,i=i:density[i](psi_inv[i](y))/(.2*np.cosh(.2*psi_inv[i](y)))
                     for i in range(3)]
        for (i,j),h in arcs.items():
            for t in (-1.,0.,1.):
                x=phi[i](t)
                y=psi[i](x)
                new_jac=(.2*np.cosh(.2*h(x)))*jac[i,j](x)/(.2*np.cosh(.2*x))
                self.assertAlmostEqual(new_density[j](changed[i,j](y))*new_jac,
                                       new_density[i](y),places=11)


if __name__ == '__main__':
    main(__name__,'interval_density_audit',report)
