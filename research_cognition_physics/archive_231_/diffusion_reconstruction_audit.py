"""Round 283: reconstruct a supplied diffusion network from operational records.

All endpoints are distinguishable; the generator is time-independent, symmetric
and Markovian. Clock units and experimental sampling costs are stated explicitly.
"""
import math
import unittest
import numpy as np
from growing_stream_audit import main
from conserved_screening_audit import laplacian
from operational_dimension_audit import decorated_torus, distances
from fisher_hanoi_scaling_audit import hanoi


def graph_laplacian(rows):
    n = len(rows)
    return laplacian(n, [(a, b, 1.) for a, row in enumerate(rows) for b in row if a < b])


def path_laplacian(n):
    return laplacian(n, [(j, j+1, 1.) for j in range(n-1)])


def heat(matrix, time):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors*np.exp(-time*values)) @ vectors.T


def heat_uniformized(matrix, time):
    """Positive Poisson expansion avoids cancellation of rare-event probabilities."""
    rate = float(np.max(np.diag(matrix)))
    if rate == 0:
        return np.eye(len(matrix))
    transition = np.eye(len(matrix))-matrix/rate
    assert transition.min() >= -1e-13
    mu = rate*time
    power = np.eye(len(matrix))
    weight = math.exp(-mu)
    out = weight*power
    for k in range(1, max(100, int(4*mu+len(matrix)+40))):
        power = power @ transition
        weight *= mu/k
        out += weight*power
    return out


def squared_diffusion_distances(probability):
    # Counting measure convention: 1/N times the uniform-stationary convention.
    return np.sum((probability[:, None, :]-probability[None, :, :])**2, axis=2)


def symmetric_log(matrix):
    values, vectors = np.linalg.eigh((matrix+matrix.T)/2)
    if values[0] <= 0:
        raise ValueError('Reconstruction requires a positive-definite kernel.')
    return (vectors*np.log(values)) @ vectors.T


def reconstruct_from_distances(squared, time):
    n = len(squared)
    projector = np.ones((n, n))/n
    centered = np.eye(n)-projector
    kernel = -centered @ squared @ centered/2+projector
    return -symmetric_log(kernel)/(2*time), kernel


def trials_for_detection(probability, failure=0.05):
    return math.ceil(math.log(failure)/math.log1p(-probability))


def poisson_tail(mu, threshold):
    weight = math.exp(-mu)*mu**threshold/math.factorial(threshold)
    total = weight
    for k in range(threshold+1, threshold+100):
        weight *= mu/k
        total += weight
    return total


def supplied_cases():
    rows, _ = hanoi(2)
    torus, _ = decorated_torus(2, 3)
    return [('weighted_path',laplacian(6,[(j,j+1,0.6+0.2*j) for j in range(5)])),
            ('hanoi_9',graph_laplacian(rows)),('three_port_torus_36',graph_laplacian(torus))]


def report():
    recovered = []
    for name, matrix in supplied_cases():
        time = 0.4
        probability = heat(matrix, time)
        reconstructed, kernel = reconstruct_from_distances(squared_diffusion_distances(probability),time)
        recovered.append({'name':name,'endpoints':len(matrix),'observation_time':time,
                          'maximum_generator_error':float(np.max(abs(reconstructed-matrix))),
                          'minimum_centered_nonconstant_eigenvalue':float(np.exp(-2*time*np.linalg.eigvalsh(matrix)[-1])),
                          'exact_diffusion_coordinate_rank':len(matrix)-1,
                          'pair_distance_mean':float(np.mean(np.sqrt(squared_diffusion_distances(probability))))})
    rare = []
    matrix = path_laplacian(5)
    for time in (0.1,0.01,0.001):
        p = float(heat_uniformized(matrix,time)[0,4])
        p2 = float(heat_uniformized(matrix,2*time)[0,4])
        rare.append({'hop_distance':4,'time':time,'endpoint_probability':p,
                     'measured_two_time_log_slope':math.log(p2/p)/math.log(2),
                     'independent_probes_for_95_percent_one_or_more_detections':trials_for_detection(p)})
    return {'round':283,
            'scope':'Exact identification within supplied finite symmetric diffusion models with fully resolved endpoints; not a derivation of a natural generator, physical clock, three dimensions or a strict light cone.',
            'reconstruction_cases':recovered,'rare_event_costs':rare,
            'calibration':'Changing clock units by a multiplies numerical time by a and divides the reconstructed generator by a.',
            'next_gap':'Fully resolved endpoint identity is essential; a hidden covering network can have the same entire visible Markov process.'}


class Audit(unittest.TestCase):
    def test_01_positive_uniformization_matches_spectral_heat(self):
        for _, matrix in supplied_cases():
            for time in (0.03,0.4,1.):
                probability = heat_uniformized(matrix,time)
                np.testing.assert_allclose(probability,heat(matrix,time),atol=2e-14)
                np.testing.assert_allclose(probability.sum(axis=1),1.,atol=2e-14)
                self.assertGreater(float(probability.min()),0.)

    def test_02_first_nonzero_matrix_power_matches_graph_distance(self):
        rows, _ = hanoi(2)
        matrix = graph_laplacian(rows)
        adjacency = np.diag(np.diag(matrix))-matrix
        for source in range(len(rows)):
            hops = distances(rows,source)
            for target in range(len(rows)):
                h = int(hops[target])
                for k in range(h):
                    self.assertAlmostEqual(np.linalg.matrix_power(-matrix,k)[source,target],0.)
                self.assertAlmostEqual(np.linalg.matrix_power(-matrix,h)[source,target],
                                       np.linalg.matrix_power(adjacency,h)[source,target])

    def test_03_short_time_statistics_recover_hop_count(self):
        matrix = path_laplacian(5)
        small = heat_uniformized(matrix,0.001)
        double = heat_uniformized(matrix,0.002)
        for h in range(1,5):
            slope = math.log(double[0,h]/small[0,h])/math.log(2)
            self.assertLess(abs(slope-h),0.003)
            leading = 0.001**h/math.factorial(h)
            self.assertLess(abs(small[0,h]/leading-1),0.003)

    def test_04_diffusion_distances_recover_centered_heat_kernel(self):
        for _, matrix in supplied_cases():
            n = len(matrix)
            probability = heat(matrix,0.4)
            squared = squared_diffusion_distances(probability)
            center = np.eye(n)-np.ones((n,n))/n
            expected = heat(matrix,0.8)-np.ones((n,n))/n
            np.testing.assert_allclose(-center @ squared @ center/2,expected,atol=2e-14)

    def test_05_generator_and_weighted_adjacency_reconstruct(self):
        for _, matrix in supplied_cases():
            for time in (0.1,0.4,0.7):
                probability = heat(matrix,time)
                reconstructed,_ = reconstruct_from_distances(squared_diffusion_distances(probability),time)
                np.testing.assert_allclose(reconstructed,matrix,atol=2e-11)
                np.testing.assert_array_equal(reconstructed < -1e-8,matrix < -1e-8)

    def test_06_permutations_and_clock_changes_calibrate(self):
        matrix = supplied_cases()[1][1]
        permutation = np.random.default_rng(283).permutation(len(matrix))
        time, scale = 0.4,3.7
        probability = heat(matrix,time)
        np.testing.assert_allclose(heat(matrix/scale,scale*time),probability,atol=2e-14)
        squared = squared_diffusion_distances(probability)[np.ix_(permutation,permutation)]
        reconstructed,_ = reconstruct_from_distances(squared,scale*time)
        np.testing.assert_allclose(reconstructed,matrix[np.ix_(permutation,permutation)]/scale,atol=1e-12)

    def test_07_finite_diffusion_embedding_rank_is_not_spatial_dimension(self):
        for _, matrix in supplied_cases():
            for time in (0.1,0.4):
                n = len(matrix)
                centered = heat(matrix,2*time)-np.ones((n,n))/n
                self.assertEqual(np.linalg.matrix_rank(centered,tol=1e-10),n-1)

    def test_08_inverse_log_has_quantified_conditioning(self):
        matrix = supplied_cases()[1][1]
        perturbation = laplacian(len(matrix),[(0,1,1e-4)])
        time = 0.7
        first,second = heat(matrix,2*time),heat(matrix+perturbation,2*time)
        floor = min(np.linalg.eigvalsh(first)[0],np.linalg.eigvalsh(second)[0])
        bound = np.linalg.norm(first-second,2)/(2*time*floor)
        self.assertLessEqual(np.linalg.norm(perturbation,2),bound*(1+1e-10))
        reconstructed,_ = reconstruct_from_distances(squared_diffusion_distances(heat(matrix+perturbation,time)),time)
        np.testing.assert_allclose(reconstructed,matrix+perturbation,atol=2e-11)

    def test_09_rare_event_sample_budget_is_recorded(self):
        probabilities = [float(heat_uniformized(path_laplacian(5),t)[0,4]) for t in (0.1,0.01,0.001)]
        counts = [trials_for_detection(p) for p in probabilities]
        self.assertTrue(all(a < b for a,b in zip(counts,counts[1:])))
        self.assertGreater(counts[-1],10**13)
        for p,count in zip(probabilities,counts):
            self.assertLessEqual(math.exp(count*math.log1p(-p)),0.05+1e-14)

    def test_10_diffusion_has_tails_not_a_strict_tick_cone(self):
        matrix = path_laplacian(9)
        time = 0.05
        probability = heat_uniformized(matrix,time)[0]
        for hops in range(1,9):
            self.assertGreater(probability[hops],0.)
            self.assertLessEqual(sum(probability[hops:]),poisson_tail(2*time,hops)*(1+1e-12))


if __name__ == '__main__':
    main(__name__,'diffusion_reconstruction_audit',report)
