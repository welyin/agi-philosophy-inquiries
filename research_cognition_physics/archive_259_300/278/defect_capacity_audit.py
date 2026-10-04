"""Round 278: charged versus neutral finite-energy defects.

Uses the classical Thomson/Nash-Williams tools on explicitly supplied lattices.
The quadratic response, conductances, source and grounded boundary are inputs.
"""
from fractions import Fraction
import itertools
import math
import unittest
import numpy as np
from growing_stream_audit import main


def apply_laplacian(u):
    out = 2*u.ndim*u.copy()
    for axis in range(u.ndim):
        left = [slice(None)]*u.ndim
        right = left.copy()
        left[axis], right[axis] = slice(1,None), slice(None,-1)
        out[tuple(left)] -= u[tuple(right)]
        out[tuple(right)] -= u[tuple(left)]
    return out


def solve_box(source, gap=0.0):
    """Dirichlet outside an odd cube: every outside nearest neighbor is ground."""
    n = source.shape[0]
    assert gap >= 0 and all(size == n for size in source.shape)
    k = np.arange(1,n+1)
    transform = np.sqrt(2/(n+1))*np.sin(np.pi*np.outer(k,k)/(n+1))
    modes = source.copy()
    eigen = np.zeros(source.shape)
    one = 4*np.sin(np.pi*k/(2*(n+1)))**2
    for axis in range(source.ndim):
        modes = np.moveaxis(np.tensordot(transform.T,modes,axes=(1,axis)),0,axis)
        shape = [1]*source.ndim
        shape[axis] = n
        eigen += one.reshape(shape)
    modes /= eigen+gap
    for axis in range(source.ndim):
        modes = np.moveaxis(np.tensordot(transform,modes,axes=(1,axis)),0,axis)
    return modes


def point_source(dimension, radius, dipole=False):
    source = np.zeros((2*radius+1,)*dimension)
    source[(radius,)*dimension] = 1.0
    if dipole:
        source[(radius+1,)+(radius,)*(dimension-1)] = -1.0
    return source


def center_moments(dimension, radius, gap=0.0):
    """Exact separable finite-box eigen-sum; slices bound peak memory."""
    k = np.arange(1,2*radius+2,2)
    one = 4*np.sin(np.pi*k/(4*radius+4))**2
    others = np.zeros((radius+1,)*(dimension-1))
    for axis in range(dimension-1):
        shape = [1]*(dimension-1)
        shape[axis] = radius+1
        others += one.reshape(shape)
    response = norm = 0.0
    for value in one:
        inverse = 1/(value+others+gap)
        response += float(np.mean(inverse))
        norm += float(np.mean(inverse**2))
    response /= radius+1
    norm /= radius+1
    return {'response':response,'potential_l2_squared':norm,
            'gradient_energy':response-gap*norm}


def cut_lower_bound(dimension, radius):
    return sum(1/(2*dimension*(2*r+1)**(dimension-1)) for r in range(radius+1))


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for rest in compositions(total-first,parts-1):
                yield (first,)+rest


def orthant_flow_level(dimension, level):
    count = math.comb(level+dimension-1,dimension-1)
    edges = {}
    for x in compositions(level,dimension):
        for axis in range(dimension):
            y = list(x)
            y[axis] += 1
            edges[(x,tuple(y))] = Fraction(x[axis]+1,count*(level+dimension))
    return edges


def dense_laplacian(dimension, radius):
    n = 2*radius+1
    matrix = np.eye(n**dimension)*2*dimension
    for x in itertools.product(range(n),repeat=dimension):
        u = np.ravel_multi_index(x,(n,)*dimension)
        for axis in range(dimension):
            if x[axis]+1 < n:
                y = list(x)
                y[axis] += 1
                v = np.ravel_multi_index(tuple(y),(n,)*dimension)
                matrix[u,v] = matrix[v,u] = -1
    return matrix


def report():
    rows = []
    for d in (1,2,3,4):
        for radius in (4,8,16,32):
            row = center_moments(d,radius)
            rows.append({'dimension':d,'radius':radius,'interior_nodes':(2*radius+1)**d,
                         'unit_charge_energy':row['response'],
                         'nash_williams_lower_bound':cut_lower_bound(d,radius)})
    dipoles = []
    for d in (1,2,3,4):
        radius = 3
        b = point_source(d,radius,True)
        u = solve_box(b)
        dipoles.append({'dimension':d,'radius':radius,
                        'energy':float(np.sum(b*u)),
                        'net_source':float(b.sum()),
                        'single_edge_feasible_energy':1.0})
    return {'round':278,
            'scope':'Classical capacity criterion under specified unscreened quadratic response; lattice calibration, not derived physical space or charge.',
            'grounded_monopoles':rows,'neutral_dipoles':dipoles,
            'orthant_flow_energy_upper_bounds':{
                str(d):float(Fraction(d-1,d-2)) for d in (3,4,5)},
            'weighted_ray':{'edge_conductance':'c(n,n+1)=(n+1)^2',
                            'first_64_edges_energy':sum(1/n**2 for n in range(1,65)),
                            'infinite_energy_upper_bound':2.0,
                            'uniform_conductance_upper_bound':False},
            'analytic_conclusions':[
                'Unit-conductance Z^d supports finite-energy net flux to infinity iff d>=3.',
                'Any finite neutral source on a connected positive-conductance network has a finite-support feasible flow.',
                'A closed finite graph cannot support nonzero total divergence without a declared return path.',
                'The criterion does not single out dimension three.']}


class CapacityTests(unittest.TestCase):
    def test_01_spectral_solver_against_independent_dense_equations(self):
        for d in (1,2,3,4):
            for dipole in (False,True):
                b = point_source(d,1,dipole)
                u = solve_box(b)
                exact = np.linalg.solve(dense_laplacian(d,1),b.ravel())
                np.testing.assert_allclose(u.ravel(),exact,atol=2e-14)

    def test_02_current_divergence_and_quadratic_energy(self):
        rng = np.random.default_rng(278)
        for d in (1,2,3):
            b = rng.normal(size=(5,)*d)
            u = solve_box(b)
            np.testing.assert_allclose(apply_laplacian(u),b,atol=3e-14)
            # Sum undirected internal and boundary-edge energies explicitly.
            energy = 0.0
            for axis in range(d):
                energy += float(np.sum(np.diff(u,axis=axis)**2))
                energy += float(np.sum(np.take(u,0,axis=axis)**2))
                energy += float(np.sum(np.take(u,-1,axis=axis)**2))
            self.assertAlmostEqual(energy,float(np.sum(b*u)),places=11)

    def test_03_line_exact_monopole_and_dipole(self):
        for radius in (1,3,8,20):
            self.assertAlmostEqual(center_moments(1,radius)['response'],(radius+1)/2)
            b = point_source(1,radius,True)
            self.assertAlmostEqual(float(b@solve_box(b)),1-1/(2*radius+2))

    def test_04_shell_lower_bounds_and_ground_recession(self):
        for d in (1,2,3,4):
            values = [center_moments(d,r)['response'] for r in (2,4,8,16)]
            self.assertTrue(all(a<b for a,b in zip(values,values[1:])))
            for r,value in zip((2,4,8,16),values):
                self.assertGreaterEqual(value+1e-13,cut_lower_bound(d,r))

    def test_05_constructed_orthant_flow_conserves_current(self):
        for d in (3,4,5):
            previous = {}
            for n in range(7):
                edges = orthant_flow_level(d,n)
                incoming,outgoing = {},{}
                for (x,y),value in previous.items():
                    incoming[y] = incoming.get(y,0)+value
                for (x,y),value in edges.items():
                    outgoing[x] = outgoing.get(x,0)+value
                self.assertEqual(sum(edges.values()),1)
                if n:
                    self.assertEqual(incoming,outgoing)
                self.assertLessEqual(sum(value**2 for value in edges.values()),
                                     Fraction(1,math.comb(n+d-1,d-1)))
                previous = edges

    def test_06_neutral_path_flow_is_finite_without_infinity_sink(self):
        for length in (1,2,7,31):
            incidence = np.zeros((length+1,length))
            for edge in range(length):
                incidence[edge,edge] = 1
                incidence[edge+1,edge] = -1
            flow = np.ones(length)
            charge = incidence@flow
            np.testing.assert_array_equal(charge,np.r_[1,np.zeros(length-1),-1])
            self.assertEqual(float(flow@flow),length)
            self.assertEqual(float(charge.sum()),0)

    def test_07_center_eigen_sum_against_solved_field(self):
        for d in (1,2,3,4):
            b = point_source(d,2)
            u = solve_box(b)
            m = center_moments(d,2)
            self.assertAlmostEqual(m['response'],float(np.sum(b*u)),places=12)
            self.assertAlmostEqual(m['potential_l2_squared'],float(np.sum(u*u)),places=12)

    def test_08_closed_network_requires_neutrality(self):
        n = 7
        lap = 2*np.eye(n)-np.roll(np.eye(n),1,axis=1)-np.roll(np.eye(n),-1,axis=1)
        np.testing.assert_array_equal(np.ones(n)@lap,np.zeros(n))
        b = np.eye(n)[0]-np.eye(n)[1]
        np.testing.assert_allclose(lap@np.linalg.pinv(lap)@b,b,atol=1e-14)
        monopole = np.eye(n)[0]
        self.assertGreater(np.linalg.norm(lap@np.linalg.pinv(lap)@monopole-monopole),0.3)

    def test_09_higher_dimension_survives_and_dipole_bound(self):
        for d in (3,4):
            self.assertLess(center_moments(d,16)['response'],(d-1)/(d-2))
        for row in report()['neutral_dipoles']:
            self.assertGreater(row['energy'],0)
            self.assertLessEqual(row['energy'],1)

    def test_10_conductance_assumption_is_material(self):
        count = 32
        c = np.arange(1,count+1,dtype=float)**2
        incidence = np.zeros((count,count))
        for n in range(count):
            incidence[n,n] = 1
            if n+1<count:
                incidence[n+1,n] = -1
        lap = (incidence*c)@incidence.T
        b = np.eye(count)[0]
        u = np.linalg.solve(lap,b)
        current = c*(incidence.T@u)
        np.testing.assert_allclose(current,np.ones(count),atol=3e-12)
        self.assertAlmostEqual(u[0],float(np.sum(1/c)),places=11)
        self.assertLess(u[0]+1/count,2.0)


if __name__ == '__main__':
    main(__name__,'defect_capacity_audit',report)
