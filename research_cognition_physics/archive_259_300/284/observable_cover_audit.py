"""Round 284: visible geometry, graph covers and an identity-memory witness.

The readout restriction is explicit. Equality of label-only stochastic records
does not imply equivalence of all quantum instruments or of marker experiments.
"""
from collections import Counter
import unittest
import numpy as np
from growing_stream_audit import main
from diffusion_reconstruction_audit import graph_laplacian, heat
from operational_dimension_audit import decorated_torus, distances


def complete_four():
    return [tuple(j for j in range(4) if j != i) for i in range(4)]


def cube():
    return [tuple(i ^ (1 << k) for k in range(3)) for i in range(8)]


def cube_label(v):
    return min(v,7-v)


def projection_matrix(labels, size):
    out = np.zeros((size,len(labels)))
    out[labels,np.arange(len(labels))] = 1.
    return out


def path_histories(rows, source, steps, label=lambda x:x):
    current = [(source,(label(source),))]
    for _ in range(steps):
        current = [(w,record+(label(w),)) for v,record in current for w in rows[v]]
    return Counter(record for _,record in current)


def lift_path(rows, source, visible_path, label):
    assert label(source) == visible_path[0]
    actual = [source]
    for target in visible_path[1:]:
        candidates = [w for w in rows[actual[-1]] if label(w) == target]
        assert len(candidates) == 1
        actual.append(candidates[0])
    return actual


def cover_ball(base, depth):
    """Truncated nonbacktracking-path tree; boundary leaves are not full degrees."""
    assert depth >= 0 and all(len(set(row)) == 3 for row in base)
    rows,labels,levels,parents = [[]],[0],[0],[-1]
    for v in range(3*2**depth-2):
        if v >= len(rows):
            break
        if levels[v] == depth:
            continue
        previous_label = labels[parents[v]] if parents[v] >= 0 else None
        for target in base[labels[v]]:
            if target == previous_label:
                continue
            child = len(rows)
            rows.append([v])
            rows[v].append(child)
            labels.append(target)
            levels.append(levels[v]+1)
            parents.append(v)
    return rows,labels,levels


def transition(rows):
    out = np.zeros((len(rows),len(rows)))
    for v,row in enumerate(rows):
        for w in row:
            out[w,v] = 1/len(row)
    return out


def report():
    small,large = complete_four(),cube()
    c = projection_matrix([cube_label(v) for v in range(8)],4)
    small_l,large_l = graph_laplacian(small),graph_laplacian(large)
    cases = []
    for time in (0.1,0.5,1.,2.):
        p_small,p_large = heat(small_l,time),heat(large_l,time)
        cases.append({'time':time,'base_true_return':float(p_small[0,0]),
                      'cover_true_return':float(p_large[0,0]),
                      'cover_visible_return':float((c @ p_large)[0,0]),
                      'maximum_full_visible_kernel_error':float(np.max(abs(c @ p_large-p_small @ c)))})
    route = [0,1,3,0]
    lifted = lift_path(large,0,route,cube_label)
    base, _ = decorated_torus(3,7)
    _,labels,levels = cover_ball(base,8)
    return {'round':284,
            'scope':'Explicitly restricted classical label readout on graph covers; memory marking expands the interface. No equivalence of full quantum theories and no derivation or uniqueness of spatial dimension.',
            'finite_cover':{'base_nodes':4,'cover_nodes':8,'degree_in_both':3,
                            'fibers':[[i for i in range(8) if cube_label(i)==j] for j in range(4)]},
            'continuous_time_cases':cases,
            'marker_witness':{'visible_route':route,'actual_cover_route':lifted,
                              'base_marker_result':1,'cover_marker_result':0,
                              'traversed_edges':3,'base_initialized_flag_bits':4,
                              'cover_initialized_flag_bits':8,
                              'uncontrolled_three_step_marked_return_base':2/9,
                              'uncontrolled_three_step_marked_return_cover':0},
            'truncated_universal_cover_example':{'base_nodes':len(base),'depth':8,
                                                 'tree_nodes':len(labels),
                                                 'distinct_visible_labels':len(set(labels)),
                                                 'underlying_tree_ball_formula':3*2**8-2,
                                                 'boundary_nodes_are_truncated':True},
            'conclusion':'Visible transition statistics can define a unique effective network without identifying the full underlying graph; persistent local identity can distinguish the cover.'}


class Audit(unittest.TestCase):
    def test_01_cube_is_a_three_regular_cover_of_complete_four(self):
        rows = cube()
        for v in range(8):
            self.assertEqual(len(set(rows[v])),3)
            self.assertEqual(set(cube_label(w) for w in rows[v]),set(complete_four()[cube_label(v)]))

    def test_02_all_hidden_initial_states_project_to_same_heat_process(self):
        c = projection_matrix([cube_label(v) for v in range(8)],4)
        small,large = graph_laplacian(complete_four()),graph_laplacian(cube())
        np.testing.assert_allclose(c @ large,small @ c,atol=1e-14)
        for time in (0.,0.1,0.5,1.,3.):
            np.testing.assert_allclose(c @ heat(large,time),heat(small,time) @ c,atol=2e-14)

    def test_03_full_finite_label_histories_have_equal_probabilities(self):
        for source in (0,3,7):
            for steps in range(1,7):
                small = path_histories(complete_four(),cube_label(source),steps)
                large = path_histories(cube(),source,steps,cube_label)
                self.assertEqual(small,large)
                self.assertEqual(sum(large.values()),3**steps)

    def test_04_exact_return_formulas_and_label_aliasing(self):
        small,large = graph_laplacian(complete_four()),graph_laplacian(cube())
        for time in (0.1,0.5,1.,3.):
            ps,pl = heat(small,time),heat(large,time)
            self.assertAlmostEqual(ps[0,0],(1+3*np.exp(-4*time))/4)
            self.assertAlmostEqual(pl[0,0],((1+np.exp(-2*time))/2)**3)
            self.assertAlmostEqual(pl[0,0]+pl[7,0],ps[0,0])
            self.assertGreater(ps[0,0]-pl[0,0],0.)

    def test_05_persistent_mark_distinguishes_controlled_triangle(self):
        route = [0,1,3,0]
        small_path = lift_path(complete_four(),0,route,lambda x:x)
        large_path = lift_path(cube(),0,route,cube_label)
        small_marks,large_marks = np.zeros(4,dtype=int),np.zeros(8,dtype=int)
        small_marks[0] = large_marks[0] = 1
        self.assertEqual(small_path,[0,1,3,0])
        self.assertEqual(large_path,[0,1,3,7])
        self.assertEqual(small_marks[small_path[-1]],1)
        self.assertEqual(large_marks[large_path[-1]],0)

    def test_06_marker_witness_needs_no_postselection(self):
        small = np.linalg.matrix_power(transition(complete_four()),3)
        large = np.linalg.matrix_power(transition(cube()),3)
        self.assertAlmostEqual(small[0,0],2/9)
        self.assertEqual(large[0,0],0.)
        self.assertAlmostEqual(large[0,0]+large[7,0],2/9)

    def test_07_universal_cover_interior_has_three_ports_and_tree_growth(self):
        base,_ = decorated_torus(3,3)
        for depth in range(1,9):
            rows,labels,levels = cover_ball(base,depth)
            self.assertEqual(len(rows),3*2**depth-2)
            self.assertEqual(sum(len(r) for r in rows)//2,len(rows)-1)
            for v,row in enumerate(rows):
                if levels[v] < depth:
                    self.assertEqual(len(row),3)
                    self.assertEqual(set(labels[w] for w in row),set(base[labels[v]]))
                else:
                    self.assertEqual(len(row),1)

    def test_08_truncated_cover_walk_is_exact_before_boundary_transitions(self):
        base,_ = decorated_torus(2,3)
        depth = 6
        rows,labels,levels = cover_ball(base,depth)
        c = projection_matrix(labels,len(base))
        fine = np.eye(len(rows))[:,0]
        coarse = np.eye(len(base))[:,0]
        base_transition = transition(base)
        for _ in range(depth):
            moved = np.zeros(len(rows))
            for v,p in enumerate(fine):
                if p:
                    self.assertLess(levels[v],depth)
                    for w in rows[v]:
                        moved[w] += p/3
            fine = moved
            coarse = base_transition @ coarse
            np.testing.assert_allclose(c @ fine,coarse,atol=1e-14)

    def test_09_visible_ball_is_base_ball_despite_larger_cover(self):
        base,_ = decorated_torus(3,7)
        rows,labels,levels = cover_ball(base,8)
        base_distances = distances(base,0)
        for radius in range(9):
            visible = {label for label,level in zip(labels,levels) if level <= radius}
            expected = set(np.flatnonzero(base_distances <= radius))
            self.assertEqual(visible,expected)
        self.assertGreater(len(rows),len(set(labels)))


if __name__ == '__main__':
    main(__name__,'observable_cover_audit',report)
