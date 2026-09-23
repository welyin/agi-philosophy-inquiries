"""Round 279: screening, finite field energy and variational attainment.

The gap gamma is an added quadratic response coefficient, not derived mass.
Finite Dirichlet boxes check exact identities; asymptotic claims use proofs.
"""
import math
import unittest
import numpy as np
from growing_stream_audit import main
from defect_capacity_audit import (
    apply_laplacian, center_moments, point_source, solve_box)


def line_screened(gap):
    if gap <= 0:
        raise ValueError('Positive screening coefficient required.')
    scale = math.sqrt(gap*(gap+4))
    ratio = 2/(gap+2+scale)
    return {'response':1/scale,'ratio':ratio,
            'potential_sum':(1+ratio)/((1-ratio)*scale),
            'l2_squared':(1+ratio**2)/((1-ratio**2)*scale**2)}


def boundary_flux(u):
    return sum(float(np.take(u,side,axis=axis).sum())
               for axis in range(u.ndim) for side in (0,-1))


def energy_budget(dimension, radius, gap):
    b = point_source(dimension,radius)
    u = solve_box(b,gap)
    gradient = float(np.sum(u*apply_laplacian(u)))
    onsite = gap*float(np.sum(u*u))
    pairing = float(np.sum(b*u))
    outward = boundary_flux(u)
    induced = gap*float(u.sum())
    return {'dimension':dimension,'radius':radius,'gap':gap,
            'gradient_energy':gradient,'onsite_energy':onsite,
            'source_pairing':pairing,
            'functional_at_minimum':(gradient+onsite)/2-pairing,
            'boundary_outward_flux':outward,'screening_response_total':induced}


def report():
    rows = []
    for d in (1,2,3,4):
        for gap in (0.0,0.01,0.1,1.0):
            row = center_moments(d,32,gap)
            rows.append({'dimension':d,'radius':32,'gap':gap,
                         **row,'onsite_energy':gap*row['potential_l2_squared'],
                         'functional_infimum':-row['response']/2})
    norms = [{'dimension':d,'radius':r,**center_moments(d,r)}
             for d in (3,4,5) for r in (4,8,16)]
    return {'round':279,
            'scope':'Specified scalar response on supplied lattices; gapless capacity versus screened response and l2 attainment, not a mass or dimension derivation.',
            'screening_cases':rows,'gapless_energy_and_norm':norms,
            'closed_budget_examples':[energy_budget(d,8,0.1) for d in (1,2,3)],
            'infinite_screened_line':[{'gap':g,**line_screened(g)} for g in (0.01,0.1,1.0)],
            'analytic_criteria':{
                'gapless_nonzero_charge_finite_gradient_energy':'d>2',
                'gapless_nonzero_charge_l2_potential':'d>4',
                'positive_gap_finite_energy_and_l2':'all integer d>=1',
                'gapless_variational_infimum':'-infinity for d<=2; finite for d>=3',
                'gapless_l2_minimizer':'absent for d=3,4 despite finite infimum; present for d>=5'}}


class ScreeningTests(unittest.TestCase):
    def test_01_equation_and_complete_energy_ledger(self):
        for d in (1,2,3):
            for gap in (0.01,0.1,1.0):
                b = point_source(d,3)
                u = solve_box(b,gap)
                np.testing.assert_allclose(apply_laplacian(u)+gap*u,b,atol=1e-14)
                row = energy_budget(d,3,gap)
                self.assertAlmostEqual(row['gradient_energy']+row['onsite_energy'],
                                       row['source_pairing'],places=12)
                self.assertAlmostEqual(row['functional_at_minimum'],
                                       -row['source_pairing']/2,places=12)

    def test_02_positive_gap_coercivity_bound(self):
        for d in (1,2,3,4):
            for gap in (0.01,0.1,1.0):
                value = center_moments(d,12,gap)
                self.assertLessEqual(value['response'],1/gap)
                self.assertLessEqual(value['potential_l2_squared'],1/gap**2)

    def test_03_exact_infinite_line_profile(self):
        for gap in (0.01,0.1,1.0):
            model = line_screened(gap)
            x = np.arange(-100,101)
            u = model['response']*model['ratio']**abs(x)
            residual = (2+gap)*u[1:-1]-u[:-2]-u[2:]
            expected = np.zeros_like(residual)
            expected[99] = 1
            np.testing.assert_allclose(residual,expected,atol=4e-15)
            self.assertAlmostEqual(gap*model['potential_sum'],1,places=12)

    def test_04_box_convergence_to_exact_screened_line(self):
        for gap in (0.01,0.1,1.0):
            exact = line_screened(gap)
            finite = center_moments(1,256,gap)
            self.assertAlmostEqual(finite['response'],exact['response'],places=11)
            self.assertAlmostEqual(finite['potential_l2_squared'],exact['l2_squared'],places=9)

    def test_05_flux_accounting_requires_screening_response(self):
        for d in (1,2,3):
            row = energy_budget(d,6,0.4)
            self.assertAlmostEqual(row['boundary_outward_flux']+
                                   row['screening_response_total'],1,places=12)
            self.assertGreater(row['screening_response_total'],0)
            self.assertLess(row['boundary_outward_flux'],1)

    def test_06_resolvent_derivative_matches_norm(self):
        for d in (1,2,3,4):
            h, gap = 1e-5,0.3
            derivative = (center_moments(d,8,gap+h)['response']-
                          center_moments(d,8,gap-h)['response'])/(2*h)
            self.assertAlmostEqual(-derivative,
                                   center_moments(d,8,gap)['potential_l2_squared'],places=7)

    def test_07_finite_minimum_by_independent_perturbations(self):
        rng = np.random.default_rng(279)
        b = point_source(2,3)
        gap = 0.2
        u = solve_box(b,gap)
        def functional(x):
            return float(np.sum(x*(apply_laplacian(x)+gap*x))/2-np.sum(b*x))
        for _ in range(10):
            v = rng.normal(size=b.shape)
            cost = float(np.sum(v*(apply_laplacian(v)+gap*v))/2)
            self.assertAlmostEqual(functional(u+v)-functional(u),cost,places=11)

    def test_08_finite_energy_is_distinct_from_potential_norm(self):
        # Finite diagnostics of the analytic Fourier thresholds, not a fit proof.
        for d in (3,4):
            a,b = center_moments(d,4),center_moments(d,32)
            self.assertLess(b['response'],0.3)
            self.assertGreater(b['potential_l2_squared'],a['potential_l2_squared'])
        self.assertGreater(center_moments(3,32)['potential_l2_squared'],1.0)

    def test_09_screening_increases_response_cost_when_removed(self):
        for d in (1,2,3,4):
            values = [center_moments(d,8,g)['response'] for g in (1,0.1,0.01,0)]
            self.assertTrue(all(a<b for a,b in zip(values,values[1:])))


if __name__ == '__main__':
    main(__name__,'screened_defect_audit',report)
