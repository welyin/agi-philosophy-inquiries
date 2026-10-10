"""Count-zero exact calibration: strict-between meetings and visible image.

No scientific round is added.  The general geometric arguments are in
screen.md; this file checks finite rational certificates, not spatial dimension.
Default prints JSON only.  --check compares saved results without writing.
"""
from fractions import Fraction as F
from itertools import product
from math import isqrt
from pathlib import Path
import argparse
import hashlib
import json


def vec(*xs):
    return tuple(F(x) for x in xs)


def sub(a,b):
    return tuple(x-y for x,y in zip(a,b))


def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),F(0))


def meet(x,y,t=F(1,2)):
    return tuple((1-t)*a+t*b for a,b in zip(x,y))


def visible(x):
    return x[:3]


def det3(rows):
    a,b,c=rows
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def run():
    vertices=[vec(*x) for x in product((-1,1),repeat=4)]
    probes=[vec(1,0,0),vec(0,1,0),vec(F(1,3),F(-1,4),F(1,5))]
    assert all(dot(r,r)<=1 for r in probes)
    pairs=0
    strict_tests=0
    tie_tests=0
    max_effect_square=F(0)
    for x in vertices:
        for y in vertices:
            m=meet(x,y)
            assert all(-1<=q<=1 for q in m)
            assert visible(m)==meet(visible(x),visible(y))
            for r in probes:
                qx,qy,qm=dot(r,visible(x)),dot(r,visible(y)),dot(r,visible(m))
                assert qm==(qx+qy)/2
                if qx==qy:
                    assert qm==qx
                    tie_tests+=1
                else:
                    assert min(qx,qy)<qm<max(qx,qy)
                    strict_tests+=1
            d=sub(visible(x),visible(y))
            effect_square=dot(d,d)/64
            max_effect_square=max(max_effect_square,effect_square)
            assert effect_square<=F(3,16)<F(1,4)
            assert sub(visible(y),visible(x))==tuple(-q for q in d)
            pairs+=1
    assert pairs==256 and max_effect_square==F(3,16)
    # For every, not just vertex, pair in [-1,1]^4 each visible difference
    # coordinate is bounded by 2.  This elementary bound gives 12/64 globally.
    component_bound=F(2)
    global_effect_square=3*component_bound**2/64
    assert global_effect_square==F(3,16)

    x=vec(-1,0,1)
    y=vec(1,1,-1)
    t=F(1,3)
    m=meet(x,y,t)
    assert m==vec(F(-1,3),F(1,3),F(1,3))
    assert m!=meet(x,y)
    nonhalf=[]
    for r in probes:
        qx,qy,qm=dot(r,x),dot(r,y),dot(r,m)
        assert qm==(1-t)*qx+t*qy
        assert qx==qy==qm or min(qx,qy)<qm<max(qx,qy)
        nonhalf.append({'probe':[str(q) for q in r],
                        'left':str(qx),'right':str(qy),'meeting':str(qm)})

    hidden_x=vec(0,0,0,1)
    hidden_y=vec(0,0,0,-1)
    hidden_m=meet(hidden_x,hidden_y)
    assert hidden_x!=hidden_y and hidden_m==vec(0,0,0,0)
    assert visible(hidden_x)==visible(hidden_y)==visible(hidden_m)==vec(0,0,0)
    assert sub(visible(hidden_x),visible(hidden_y))==vec(0,0,0)
    # Thus their complete comparison effect is I/2, not merely equal on
    # the finite probe set above.  The fourth coordinate remains distinct.

    anchors=[vec(0,0,0,0),vec(1,0,0,0),vec(0,1,0,0),vec(0,0,1,0)]
    anchor_rows=[list(sub(visible(x),visible(anchors[0]))) for x in anchors[1:]]
    anchor_det=det3(anchor_rows)
    assert anchor_det==1

    # Rational points are closed under rational non-half meetings.  These
    # exact finite examples calibrate the field formula, not its universality.
    rational_pairs=[(vec(F(1,3),F(-2,7),F(5,11)),vec(F(-4,5),F(1,13),F(-3,8))),
                    (vec(-1,0,1),vec(1,F(3,10),-1))]
    rational_examples=[]
    for a,b in rational_pairs:
        for weight in (F(1,2),F(1,3)):
            c=meet(a,b,weight)
            assert all(isinstance(q,F) and -1<=q<=1 for q in c)
            assert all(c[k]==(1-weight)*a[k]+weight*b[k] for k in range(3))
            rational_examples.append({'left':[str(q) for q in a],
                                      'right':[str(q) for q in b],
                                      'weight':str(weight),'meeting':[str(q) for q in c]})
    # q_n < sqrt(2)/2 < u_n are exact integer-square brackets, shrinking
    # inside [0,1].  The omitted limit is irrational by the standard
    # parity proof.  No floating-point sqrt or dimensional inference is used.
    brackets=[]
    previous=None
    for n in range(1,7):
        denominator=10**n
        numerator=isqrt(2*denominator**2)
        lower=F(numerator,2*denominator)
        upper=F(numerator+1,2*denominator)
        assert lower**2<F(1,2)<upper**2
        assert upper-lower==F(1,2*denominator)
        assert 0<lower<upper<1
        if previous is not None:
            assert previous[0]<=lower and upper<=previous[1]
        previous=(lower,upper)
        brackets.append({'n':n,'lower':str(lower),'upper':str(upper),
                         'width':str(upper-lower),
                         'lower_square':str(lower**2),'upper_square':str(upper**2)})

    return {
        'status':'exact_fraction_screen_passed',
        'formal_round_increment':0,
        'scientific_calibration_increment':0,
        'validation_group_count':5,
        'arithmetic':'fractions.Fraction and integer square roots; no scientific floating point',
        'four_cube_meeting':{
            'domain':'X=[-1,1]^3 x [-1,1]',
            'visible_map':'v(x1,x2,x3,h)=(x1,x2,x3)',
            'actual_model_operation':'m(x,y)=(x+y)/2 on all four coordinates',
            'vertex_count':len(vertices),'ordered_vertex_pairs':pairs,
            'finite_physical_probe_count':len(probes),
            'strict_probe_tests':strict_tests,'tie_probe_tests':tie_tests,
            'global_all_probe_identity':'r.v(m)=(r.v(x)+r.v(y))/2',
            'effect':'E_xy=I/2+(v(x)-v(y)).sigma/8',
            'global_component_difference_bound':str(component_bound),
            'global_effect_Bloch_norm_squared_upper':str(global_effect_square),
            'positivity_threshold':'1/4',
            'effect_positivity_margin':str(F(1,4)-global_effect_square),
            'interpretation':'An abstract admissible countermodel, not an implementation by the old autonomous Hamiltonian.'},
        'nonhalf_meeting':{'weight':str(t),'left':[str(q) for q in x],
                           'right':[str(q) for q in y],'meeting':[str(q) for q in m],
                           'finite_checks':nonhalf},
        'hidden_antipodal_pair':{'left':[str(q) for q in hidden_x],
                                 'right':[str(q) for q in hidden_y],
                                 'meeting':[str(q) for q in hidden_m],
                                 'visible_difference':['0','0','0'],
                                 'complete_comparison_effect':'I/2',
                                 'full_endpoints_equal':False},
        'four_anchor_affine_rank':{'anchors':[[str(q) for q in a] for a in anchors],
                                  'visible_difference_rows':[[str(q) for q in row] for row in anchor_rows],
                                  'determinant':str(anchor_det),'affine_rank':3},
        'rational_closure_boundary':{
            'domain':'S=Q^3 intersect [-1,1]^3, inherited Euclidean topology',
            'meetings':'arithmetic average; t=1/3 also preserves rationality',
            'closure_formula':'(1-a/b)*(m/n)+(a/b)*(p/q)=((b-a)*m*q+a*p*n)/(b*n*q)',
            'finite_rational_examples':rational_examples,
            'omitted_limit':'(sqrt(2)/2,0,0)',
            'integer_square_brackets':brackets,
            'analytic_boundary':'Rational density and irrationality imply S is not closed and has empty R3 interior; these are analytic facts, not conclusions of finite sampling.'},
        'scope':{'general_convexity_proof_done_by_this_code':False,
                 'random_mixture_records_claimed_to_be_actual_meetings':False,
                 'actual_qubit_process_for_meeting_constructed':False,
                 'actual_spatial_dimension_derived':False,
                 'finite_samples_establish_topological_dimension':False,
                 'full_endpoint_injectivity_imposed_as_new_universal_principle':False},
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check',action='store_true',help='read-only exact comparison with results.json')
    args=parser.parse_args()
    result=run()
    if args.check:
        path=Path(__file__).with_name('results.json')
        if result!=json.loads(path.read_text(encoding='utf-8')):
            raise AssertionError('saved result differs from exact recomputation')
        print(json.dumps({'status':'exact_saved_results_match','formal_round_increment':0,
                          'validation_group_count':result['validation_group_count'],
                          'checker_sha256':result['checker_sha256'],
                          'results_sha256':hashlib.sha256(path.read_bytes()).hexdigest()},
                         ensure_ascii=False,indent=2))
    else:
        print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
