"""Round 1015: conditional Yang-Mills interaction selection.

Finite first-jet checks supplement, rather than replace, the local deformation
classification and the analytic identities in research_note_1015.md.
"""
from pathlib import Path
from itertools import permutations
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
OUT = HERE / "vector_consistency_selection_results.json"
ETA = np.diag([-1., 1., 1., 1.])
METRIC_SIGNS = np.diag(ETA)


def three_form(n, seeds):
    out = np.zeros((n,n,n), dtype=np.int64)
    for inds, value in seeds:
        for p in permutations(range(3)):
            inversions = sum(p[i] > p[j] for i in range(3) for j in range(i+1,3))
            out[tuple(inds[i] for i in p)] = value * (-1)**inversions
    return out


def bracket(f, x, y):
    return np.einsum("abc,...b,...c->...a", f, x, y)


def jacobi(f, x, y, z):
    return bracket(f,x,bracket(f,y,z)) + bracket(f,y,bracket(f,z,x)) + bracket(f,z,bracket(f,x,y))


def jacobi_tensor(f):
    return (np.einsum("abe,ecd->abcd",f,f) + np.einsum("ace,edb->abcd",f,f)
            + np.einsum("ade,ebc->abcd",f,f))


def contract(x, y):
    return float(np.einsum("m,n,mna,mna->", METRIC_SIGNS, METRIC_SIGNS, x, y))


def field_parts(f, a, da):
    f0 = da - np.swapaxes(da,0,1)
    b = np.stack([np.stack([bracket(f,x,y) for y in a]) for x in a])
    return f0,b


def action_parts(f, a, da):
    f0,b = field_parts(f,a,da)
    return np.array([-.25*contract(f0,f0), -.5*contract(f0,b), -.25*contract(b,b)])


def density(f,a,da,g,quartic=1.):
    return float(action_parts(f,a,da) @ np.array([1.,g,quartic*g*g]))


def gauge_tangent(f,a,da,e,de,dde,g):
    va = de + g*bracket(f,a,e)
    vda = np.empty_like(da, dtype=float)
    for mu in range(4):
        for nu in range(4):
            vda[mu,nu] = dde[mu,nu] + g*(bracket(f,da[mu,nu],e) + bracket(f,a[nu],de[mu]))
    return va,vda


def variation(f,a,da,e,de,dde,g,quartic=1.):
    va,vda = gauge_tangent(f,a,da,e,de,dde,g)
    f0,b = field_parts(f,a,da)
    vf0 = vda - np.swapaxes(vda,0,1)
    vb = np.empty_like(b)
    for mu in range(4):
        for nu in range(4):
            vb[mu,nu] = bracket(f,va[mu],a[nu])+bracket(f,a[mu],va[nu])
    variation_parts = np.array([-.5*contract(f0,vf0),
        -.5*(contract(vf0,b)+contract(f0,vb)), -.5*contract(b,vb)])
    df = vf0+g*vb
    full = f0+g*b
    cov = df-g*bracket(f,full,e)
    j = np.stack([np.stack([jacobi(f,x,y,e) for y in a]) for x in a])
    return dict(tangent=(va,vda),action=float(variation_parts@[1.,g,quartic*g*g]),
        covariance_defect=cov,jacobi_covariance_prediction=g*g*j,
        expanded_action_parts=action_parts(f,a,da))


def closure(f,a,e,de,z,dz,g):
    va = de+g*bracket(f,a,e)
    wa = dz+g*bracket(f,a,z)
    comm = g*(bracket(f,va,z)-bracket(f,wa,e))
    beta = g*bracket(f,e,z)
    dbeta = g*(bracket(f,de,z)+bracket(f,e,dz))
    residual = comm-(dbeta+g*bracket(f,a,beta))
    predicted = -g*g*jacobi(f,a,e,z)
    return residual,predicted


def maxabs(a):
    return float(np.max(np.abs(a)))


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def run():
    su2 = three_form(3,[((0,1,2),1)])
    su2u1 = three_form(4,[((0,1,2),1)])
    product = three_form(6,[((0,1,2),1),((3,4,5),2)])
    bad = three_form(5,[((0,1,2),1),((0,3,4),1)])
    abelian = np.zeros((3,3,3),dtype=np.int64)
    families = [("abelian_three_vectors",abelian), ("su2",su2),
                ("su2_plus_abelian",su2u1), ("two_su2_with_relative_strength_two",product)]
    exact = []
    for name,f in families + [("antisymmetric_non_jacobi",bad)]:
        assert np.array_equal(f,-np.swapaxes(f,0,1))
        assert np.array_equal(f,-np.swapaxes(f,1,2))
        j = jacobi_tensor(f)
        nonzero = int(np.count_nonzero(j))
        assert (nonzero == 0) == (name != "antisymmetric_non_jacobi")
        exact.append(dict(family=name,vectors=len(f),jacobi_nonzero_entries=nonzero,
                          largest_integer_jacobi_component=int(np.max(np.abs(j)))))

    rng = np.random.default_rng(1015)
    checks = []
    largest_action = largest_covariance = largest_closure = largest_difference = 0.
    for name,f in families:
        n = len(f)
        local_action = local_cov = local_close = local_diff = 0.
        for _ in range(12):
            a = rng.integers(-3,4,(4,n))/5
            da = rng.integers(-3,4,(4,4,n))/7
            e = rng.integers(-3,4,n)/3
            de = rng.integers(-3,4,(4,n))/11
            dde0 = rng.integers(-3,4,(4,4,n))/13
            dde = (dde0+np.swapaxes(dde0,0,1))/2
            z = rng.integers(-3,4,n)/4
            dz = rng.integers(-3,4,(4,n))/9
            g = .375
            v = variation(f,a,da,e,de,dde,g)
            c, cp = closure(f,a,e,de,z,dz,g)
            local_action = max(local_action,abs(v["action"]))
            local_cov = max(local_cov,maxabs(v["covariance_defect"]))
            local_close = max(local_close,maxabs(c))
            assert maxabs(c-cp) < 1e-12
            assert maxabs(v["covariance_defect"]-v["jacobi_covariance_prediction"]) < 1e-12
            va,vda = v["tangent"]
            step = 1e-5
            central = (density(f,a+step*va,da+step*vda,g)-density(f,a-step*va,da-step*vda,g))/(2*step)
            local_diff = max(local_diff,abs(central-v["action"]))
            f0,b = field_parts(f,a,da)
            assert abs(density(f,a,da,g)+.25*contract(f0+g*b,f0+g*b)) < 1e-12
        assert max(local_action,local_cov,local_close) < 1e-11
        assert local_diff < 1e-7
        largest_action=max(largest_action,local_action)
        largest_covariance=max(largest_covariance,local_cov)
        largest_closure=max(largest_closure,local_close)
        largest_difference=max(largest_difference,local_diff)
        checks.append(dict(family=name,first_jets=12,action_residual=local_action,
            covariance_residual=local_cov,closure_residual=local_close,
            independent_directional_difference_error=local_diff))

    # A fully antisymmetric cubic coefficient need not extend to a Lie algebra.
    a=np.zeros((4,5)); a[0,1]=1; a[1,2]=1
    da=np.zeros((4,4,5)); da[0,1,4]=1
    e=np.eye(5)[3]; de=np.zeros((4,5)); dde=np.zeros((4,4,5))
    basis=np.eye(5,dtype=np.int64)
    witness=jacobi(bad,basis[1],basis[2],basis[3])
    assert np.array_equal(witness,np.array([0,0,0,0,-1]))
    failure=[]
    for gq in [F(1,4),F(1,2),F(1)]:
        g=float(gq)
        v=variation(bad,a,da,e,de,dde,g)
        assert maxabs(v["covariance_defect"]-v["jacobi_covariance_prediction"]) == 0
        assert v["action"] == -g*g
        step=1e-5
        va,vda=v["tangent"]
        fd=(density(bad,a+step*va,da+step*vda,g)-density(bad,a-step*va,da-step*vda,g))/(2*step)
        assert abs(fd-v["action"]) < 1e-7
        failure.append(dict(g=str(gq),action_variation=str(-gq*gq),
            covariance_defect_max=str(gq*gq),independent_difference_error=abs(fd-v["action"])))
    # Separate gauge-commutator witness, including its sign.
    ac=np.zeros((4,5)); ac[0]=basis[1]
    cr,cp=closure(bad,ac,basis[2],np.zeros_like(ac),basis[3],np.zeros_like(ac),.5)
    assert maxabs(cr-cp)==0 and maxabs(cr)==.25

    # Independently vary the quartic coefficient with the cubic vertex fixed.
    a=np.zeros((4,3)); a[0,0]=1; a[1,1]=1
    da=np.zeros((4,4,3)); e=np.zeros(3); de=np.zeros((4,3)); de[0,0]=1
    dde=np.zeros((4,4,3))
    quartic=[]
    for alpha in [F(0),F(1,2),F(1),F(3,2)]:
        v=variation(su2,a,da,e,de,dde,.5,float(alpha))
        expected=(alpha-1)/4
        assert v["action"]==float(expected)
        quartic.append(dict(quartic_multiplier=str(alpha),g="1/2",gauge_variation=str(expected)))

    # Jacobi alone is insufficient for a positive invariant internal metric.
    affine=np.zeros((2,2,2),dtype=np.int64)
    affine[1,0,1]=1; affine[1,1,0]=-1
    assert not np.count_nonzero(jacobi_tensor(affine))
    # ad(e1)=diag(0,1); invariance forces k12=k22=0, hence no positive k.
    ad=np.diag([0,1])
    kbases=[np.array([[1,0],[0,0]]),np.array([[0,1],[1,0]]),np.array([[0,0],[0,1]])]
    metric_obstructions=[(ad.T@k+k@ad).tolist() for k in kbases]
    assert metric_obstructions==[[[0,0],[0,0]],[[0,1],[1,0]],[[0,0],[0,2]]]
    a=np.zeros((4,2)); da=np.zeros((4,4,2)); da[0,1,1]=1
    va=variation(affine,a,da,np.array([1.,0.]),np.zeros((4,2)),np.zeros((4,4,2)),.5)
    assert va["action"] == -.5

    # A dimension-eight Abelian operator is an allowed effective freedom, not
    # a counterexample to the dimension-four classification.
    one=np.zeros((1,1,1)); a=np.zeros((4,1)); da=np.zeros((4,4,1)); da[0,1,0]=1
    f0,_=field_parts(one,a,da)
    invariant=contract(f0,f0)
    highdim=F(3,10)*F(int(invariant))**2/F(2)**4
    tangent=gauge_tangent(one,a,da,np.ones(1),np.ones((4,1)),np.ones((4,4,1)),1.)
    changed_f0,_=field_parts(one,a+tangent[0],da+tangent[1])
    assert maxabs(changed_f0-f0)==0

    sources=["archive_301_341/research_note_326.md", "archive_342_369/research_note_358.md",
             "archive_1009_/1009/input_dependency_ledger_v0_1.md",
             "archive_956_989/981/drafts/common_parent_contract_v1.md"]
    return dict(round=1015,date="2026-10-08",new_calibration_groups=1,
        cumulative_research_groups=3793,cumulative_count_pending_primary_confirmation=False,
        all_scientific_calibrations_passed=True,exact_structure_tests=exact,
        first_jet_calibrations=checks,first_jet_count=48,
        max_action_residual=largest_action,max_covariance_residual=largest_covariance,
        max_closure_residual=largest_closure,max_independent_directional_error=largest_difference,
        non_jacobi_witness=dict(seeds={"f123":1,"f145":1},
            jacobi_e2_e3_e4=witness.tolist(),failures=failure,closure_defect_at_g_half=.25),
        independent_quartic_test=quartic,
        noncompact_affine_boundary=dict(jacobi_passed=True,invariant_metric_obstructions=metric_obstructions,
            positive_invariant_metric_exists=False,delta_metric_first_order_action_variation=-.5),
        higher_dimension_abelian_freedom=dict(operator="(F_munu F^munu)^2",operator_dimension=8,
            coefficient="3/10",cutoff="2",sample_density_shift=str(highdim),gauge_invariance_residual=0.),
        claim_boundaries=dict(SM_group_selected=False,vector_field_number_selected=False,
            physical_couplings_predicted=False,spacetime_dimension_derived=False,
            massless_vector_modes_generated=False,matter_representation_generated=False,
            nonzero_interactions_forced=False,full_quantum_existence_proved=False,
            all_EFT_operators_fixed=False,new_cognitive_axiom=False,goal_completed=False),
        historical_source_sha256={p:sha(BASE/p) for p in sources},code_sha256=sha(Path(__file__)))


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("--write",action="store_true")
    args=parser.parse_args(); result=run()
    if args.write:
        with OUT.open("x",encoding="utf8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False); stream.write("\n")
    else:
        assert result==json.loads(OUT.read_text("utf8"))
    print(json.dumps({k:v for k,v in result.items() if k in ["round","new_calibration_groups",
        "cumulative_research_groups","all_scientific_calibrations_passed","first_jet_count",
        "max_action_residual","max_independent_directional_error"]},ensure_ascii=False))
