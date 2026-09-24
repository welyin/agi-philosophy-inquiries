"""Round 413: coordinates from calibrated displacement effects, without budget queries.

All geometry and effect access are declared model inputs. Reconstruction uses
only endpoint identities and four Born probabilities per directional query.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import platform
import unittest
import numpy as np

TARGET=Path(__file__).with_name("direction_only_coordinate_audit_results.json")
SIGMA=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=complex)
ETA=.7
PREPS=[np.eye(2)/2]+[(np.eye(2)+ETA*s)/2 for s in SIGMA]


def unit(v):
    v=np.asarray(v,dtype=float)
    norm=np.linalg.norm(v)
    if norm==0:
        raise ValueError("zero directional contrast")
    return v/norm


def rotation():
    axis=unit([1.,-2.,3.])
    x,y,z=axis
    k=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
    angle=.63
    return np.eye(3)+np.sin(angle)*k+(1-np.cos(angle))*k@k


class EffectOracle:
    """Stationary endpoints and repeatable probe effects are explicit inputs."""
    def __init__(self, scale=1., frame=None, blind=False):
        self.scale=scale
        self.frame=np.eye(3) if frame is None else frame
        self.blind=blind
        self._positions=[]
        self.query_count=0

    def add(self, x):
        self._positions.append(self.scale*self.frame@np.asarray(x,dtype=float))
        return len(self._positions)-1

    def truth(self, token):
        return self._positions[token].copy()

    def effect(self, a, b):
        # The model scale changes the contrast profile as well as all endpoints.
        y=(self._positions[b]-self._positions[a])/self.scale
        r=np.linalg.norm(y)
        if self.blind:
            c=.5
            w=.5*y/(1+r*r)*((r-2)**2/(1+(r-2)**2))
        else:
            c=.45+.03*r*r/(1+r*r)
            w=.8*y/(1+r*r)
        return c*np.eye(2)+.5*np.einsum('i,ijk->jk',w,SIGMA)

    def probabilities(self, a, b):
        self.query_count+=1
        effect=self.effect(a,b)
        return np.array([np.trace(effect@rho).real for rho in PREPS])


def direction(probabilities, probability_error=0., contrast_floor=1e-12):
    p=np.asarray(probabilities,dtype=float)
    w=2*(p[1:]-p[0])/ETA
    size=float(np.linalg.norm(w))
    error=4*np.sqrt(3)*probability_error/ETA
    if size<=max(error,contrast_floor):
        raise ValueError("direction not certified above probability error")
    return w/size, min(2.,2*error/size), size


def intersect(a, b, u, v):
    m=np.column_stack((u,-v))
    singular=np.linalg.svd(m,compute_uv=False)
    if singular[-1]<=1e-10:
        raise ValueError("collinear directional queries")
    scales=np.linalg.lstsq(m,b-a,rcond=None)[0]
    return a+scales[0]*u,scales,float(singular[-1])


def certified_intersection(a,b,u,v,ea,eb,eu,ev):
    x,t,sigma=intersect(a,b,u,v)
    matrix_error=np.hypot(eu,ev)
    if sigma<=matrix_error:
        raise ValueError("no separation from angular degeneracy")
    bound_t=(ea+eb+matrix_error*np.linalg.norm(t))/(sigma-matrix_error)
    bound_x=ea+(1+eu)*bound_t+abs(t[0])*eu
    return x,float(bound_x),dict(sigma=sigma,matrix_error=float(matrix_error),scale_error_bound=float(bound_t))


class DirectionChart:
    def __init__(self, probabilities, anchor_tokens):
        self.probabilities=probabilities
        self.anchors=tuple(anchor_tokens)
        o,a,c=self.anchors
        baseline=direction(probabilities(o,a))[0]
        u=direction(probabilities(o,c))[0]
        v=direction(probabilities(a,c))[0]
        third,_,_=intersect(np.zeros(3),baseline,u,v)
        self.locations=[np.zeros(3),baseline,third]

    def locate(self, token):
        if token in self.anchors:
            return self.locations[self.anchors.index(token)].copy()
        bearings=[direction(self.probabilities(a,token))[0] for a in self.anchors]
        pair=max(((i,j) for i in range(3) for j in range(i+1,3)),
                 key=lambda pair:1-abs(float(bearings[pair[0]]@bearings[pair[1]])))
        i,j=pair
        return intersect(self.locations[i],self.locations[j],bearings[i],bearings[j])[0]


@lru_cache(None)
def report():
    oracle=EffectOracle()
    anchor_values=[[0.,0.,0.],[1.,0.,0.],[.3,.9,.2]]
    anchors=[oracle.add(v) for v in anchor_values]
    values=[[.7,-.4,.8],[-.5,.3,.4],[2.,0.,0.],[.5,0.,0.],[0.,0.,1.4],[.2,-.8,1.1]]
    points=[oracle.add(v) for v in values]
    chart=DirectionChart(oracle.probabilities,anchors)
    reconstruction=max(np.linalg.norm(chart.locate(p)-oracle.truth(p)) for p in points+anchors)
    probability_formula_error=0.
    min_eigenvalue=1.
    max_eigenvalue=0.
    angular_error=0.
    for a in anchors:
        for p in points:
            effect=oracle.effect(a,p)
            eigenvalues=np.linalg.eigvalsh(effect)
            min_eigenvalue=min(min_eigenvalue,float(eigenvalues[0]))
            max_eigenvalue=max(max_eigenvalue,float(eigenvalues[-1]))
            w=np.array([np.trace(effect@s).real for s in SIGMA])
            c=np.trace(effect).real/2
            direct=np.concatenate(([c],c+ETA*w/2))
            probability_formula_error=max(probability_formula_error,np.linalg.norm(direct-oracle.probabilities(a,p)))
            recovered,_,_=direction(oracle.probabilities(a,p))
            angular_error=max(angular_error,np.linalg.norm(recovered-unit(oracle.truth(p)-oracle.truth(a))))

    # A full displacement can be composed; the simulator does not expose vectors to the chart.
    p,q=points[:2]
    combined=oracle.add(oracle.truth(p)+oracle.truth(q))
    addition_error=np.linalg.norm(chart.locate(combined)-chart.locate(p)-chart.locate(q))

    r=rotation()
    changed=EffectOracle(scale=7.,frame=r)
    changed_anchors=[changed.add(v) for v in anchor_values]
    changed_points=[changed.add(v) for v in values]
    changed_chart=DirectionChart(changed.probabilities,changed_anchors)
    frame_error=max(np.linalg.norm(changed_chart.locate(q)-r@chart.locate(p)) for p,q in zip(points,changed_points))

    # Raw probability errors, including the measured baseline direction, enter the bound.
    eps=2e-6
    patterns=[np.array([1.,-.3,.7,-.8]),np.array([-.5,1.,-.7,.4]),np.array([.8,-.6,.2,-1.])]
    bearing_pairs=[(anchors[0],anchors[1]),(anchors[0],points[0]),(anchors[1],points[0])]
    observed=[direction(oracle.probabilities(a,b)+eps*z,eps) for (a,b),z in zip(bearing_pairs,patterns)]
    (baseline,eb,_),(u,eu,_),(v,ev,_)=observed
    estimated,bound,diagnostic=certified_intersection(np.zeros(3),baseline,u,v,0.,eb,eu,ev)
    actual_error=float(np.linalg.norm(estimated-oracle.truth(points[0])))

    # All normalized directions from two collinear anchors coincide for these endpoints.
    line_a=oracle.add([2.,0.,0.]);line_b=oracle.add([3.,0.,0.])
    line_difference=max(np.linalg.norm(direction(oracle.probabilities(a,line_a))[0]-direction(oracle.probabilities(a,line_b))[0]) for a in anchors[:2])
    resolved=float(np.linalg.norm(chart.locate(line_a)-chart.locate(line_b)))

    # The old unit-shell direction contract permits a blind shell away from it.
    blind=EffectOracle(blind=True)
    ba=[blind.add(v) for v in ([0,0,0],[1,0,0],[0,1,0])]
    z=np.sqrt(3.5)
    bp=[blind.add(v) for v in ([.5,.5,z],[.5,.5,-z])]
    blind_difference=max(np.linalg.norm(blind.probabilities(a,bp[0])-blind.probabilities(a,bp[1])) for a in ba)
    unit_probe=blind.add([0,0,1])
    _,_,unit_contrast=direction(blind.probabilities(ba[0],unit_probe))
    rejection=0
    for a in ba:
        for p in bp:
            try:direction(blind.probabilities(a,p))
            except ValueError:rejection+=1

    return dict(round=413,
        scope="Conditional direction-only coordinate readout replacing numerical budget access; endpoint structure and nonzero covariant displacement effects remain inputs",
        probes=dict(minimum_effect_eigenvalue=min_eigenvalue,maximum_effect_eigenvalue=max_eigenvalue,independent_born_probability_error=float(probability_formula_error),direction_reconstruction_error=float(angular_error),preparations_per_probability_family=4),
        coordinates=dict(maximum_reconstruction_error=float(reconstruction),addition_error=float(addition_error),common_frame_and_unit_change_error=float(frame_error),budget_queries=0,known_length_unit="one actual baseline chosen as unit; no absolute length derived"),
        raw_probability_certificate=dict(probability_error=eps,position_error=actual_error,position_error_upper_bound=bound,**diagnostic),
        degeneracy=dict(two_anchor_normalized_direction_difference=float(line_difference),third_anchor_separates_endpoints=resolved),
        blind_shell=dict(unit_shell_contrast=unit_contrast,distinct_target_separation=2*z,all_anchor_probability_difference=float(blind_difference),rejected_queries=rejection),
        full_displacement_effect_access_is_additional_input=True,
        budget_numerical_oracle_used=False,budget_structure_removed=False,
        global_nonzero_contrast_inferred_from_unit_shell=False,
        unknown_single_state_tomography_or_cloning_claimed=False,
        position_topology_generated_from_cognitive_axioms=False,
        exact_group_contract_certified_by_finite_tests=False,new_cognitive_axiom_adopted=False,
        phase_closure_triggered=False)


class Audit(unittest.TestCase):
    def test_01_legal_effects_and_independent_born_probabilities(self):
        p=report()["probes"]
        self.assertGreater(p["minimum_effect_eigenvalue"],0)
        self.assertLess(p["maximum_effect_eigenvalue"],1)
        self.assertLess(p["independent_born_probability_error"],1e-14)

    def test_02_unknown_scalar_bias_and_radial_contrast_removed(self):
        self.assertLess(report()["probes"]["direction_reconstruction_error"],1e-13)

    def test_03_coordinates_from_only_probability_queries(self):
        r=report()["coordinates"]
        self.assertLess(r["maximum_reconstruction_error"],1e-12)
        self.assertEqual(r["budget_queries"],0)

    def test_04_addition_and_common_reference_change(self):
        r=report()["coordinates"]
        self.assertLess(r["addition_error"],1e-12)
        self.assertLess(r["common_frame_and_unit_change_error"],1e-12)

    def test_05_raw_probability_errors_have_a_posteriori_bound(self):
        r=report()["raw_probability_certificate"]
        self.assertLess(r["matrix_error"],r["sigma"])
        self.assertLessEqual(r["position_error"],r["position_error_upper_bound"])

    def test_06_two_anchor_collinearity_and_actual_extra_query(self):
        r=report()["degeneracy"]
        self.assertLess(r["two_anchor_normalized_direction_difference"],1e-13)
        self.assertAlmostEqual(r["third_anchor_separates_endpoints"],1.,places=12)
        with self.assertRaises(ValueError):
            intersect(np.zeros(3),np.array([1.,0,0]),np.array([1.,0,0]),np.array([1.,0,0]))

    def test_07_blind_shell_cannot_be_certified_from_unit_shell(self):
        r=report()["blind_shell"]
        self.assertGreater(r["unit_shell_contrast"],.1)
        self.assertGreater(r["distinct_target_separation"],3)
        self.assertLess(r["all_anchor_probability_difference"],1e-13)
        self.assertEqual(r["rejected_queries"],6)

    def test_08_weak_contrast_and_angular_uncertainty_are_rejected(self):
        with self.assertRaises(ValueError):direction([.5,.5,.5,.5],1e-4)
        with self.assertRaises(ValueError):
            certified_intersection(np.zeros(3),np.array([1.,0,0]),unit([1.,1e-5,0]),unit([1.,-1e-5,0]),0.,0.,1e-3,1e-3)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    tests=unittest.TextTestRunner().run(unittest.defaultTestLoader.loadTestsFromTestCase(Audit))
    if not tests.wasSuccessful():raise SystemExit(1)
    result=dict(report(),checks=dict(run=tests.testsRun,failures=len(tests.failures),errors=len(tests.errors)),runtime=dict(python=platform.python_version(),numpy=np.__version__))
    if args.write_results:
        with TARGET.open("x",encoding="utf-8",newline="\n") as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps(result,ensure_ascii=False,indent=2))
