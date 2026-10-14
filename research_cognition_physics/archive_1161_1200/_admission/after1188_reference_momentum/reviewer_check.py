"""Independent finite-jet and partial velocity checks; default is read-only.

No author imports. This uses an ambient flat metric and the pullback chain
rule, rather than the author's chart-Christoffel calculation. Not an on-shell
solution or a complete Hamiltonian/quantum constraint calculation.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import json

ZERO = Q(0)
def matrix(rows):
    return [[Q(v) for v in row] for row in rows]
def plus(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]
def times(c, a):
    return [[c*a[i][j] for j in range(3)] for i in range(3)]
def pair(a, b):
    return sum((a[i][j]*b[i][j] for i in range(3) for j in range(3)), ZERO)
def trace(a, hi):
    return pair(a, hi)
def raise_two(a, hi):
    return [[sum((hi[i][k]*hi[j][l]*a[k][l] for k in range(3) for l in range(3)), ZERO)
             for j in range(3)] for i in range(3)]
def norm(a, hi):
    return pair(a, raise_two(a, hi))
def fref(k, c, hi, alpha, beta):
    b = plus(k, times(-1, c))
    return alpha*norm(b, hi) + beta*trace(b, hi)**2
def response(b, hi, alpha, beta):
    return plus(times(alpha, raise_two(b, hi)), times(beta*trace(b, hi), hi))
def delta(values):
    return values[0]+values[2]-2*values[1]

def run():
    # det h=1; this nonorthogonal basis differs from the author fixtures.
    h = matrix([[1, 1, 0], [1, 2, 0], [0, 0, 1]])
    hi = matrix([[2, -1, 0], [-1, 1, 0], [0, 0, 1]])
    for i in range(3):
        for j in range(3):
            assert sum((h[i][a]*hi[a][j] for a in range(3)), ZERO) == Q(i==j)
    shapes = {
        'rank_one': matrix([[1, 1, 0], [1, 1, 0], [0, 0, 0]]),
        'trace_free': matrix([[1, 1, 0], [1, 0, 0], [0, 0, 0]])}
    assert (trace(shapes['rank_one'],hi), norm(shapes['rank_one'],hi)) == (1,1)
    assert (trace(shapes['trace_free'],hi), norm(shapes['trace_free'],hi)) == (0,2)
    k0 = matrix([[Q(2,7), Q(-1,9), Q(1,5)], [Q(-1,9), Q(3,8), 0], [Q(1,5),0,Q(-2,11)]])
    d = matrix([[Q(1,3), Q(1,7), 0], [Q(1,7), Q(-2,5), Q(1,11)], [0,Q(1,11),Q(4,9)]])
    q, alpha, beta, s, lapse = Q(5,4), Q(3,8), Q(-1,10), Q(2,5), Q(3,2)
    # Same physical T=q*t+q*d_ij*x^i*x^j/2 in all embeddings t=f(x).
    # At the event df=0. Embedding chain rule adds q*D2f to pullback Hessian.
    # Ambient normalized gradient has spatial derivative du_j/dx_i=-d_ij.
    ambient_b = times(-1, d)
    basis = []
    for i in range(3):
        for j in range(i,3):
            e = matrix([[0]*3 for _ in range(3)])
            e[i][j] = e[j][i] = Q(1)
            basis.append(e)
    cases = []
    velocity_checks = 0
    for name, shape in shapes.items():
        complete, bare, mixed, momenta, ks, cs = [], [], [], [], [], []
        for sign in (-1,0,1):
            k = plus(k0, times(sign*s,shape))
            pulled_hessian = times(q, plus(d,k))
            c = times(1/q,pulled_hessian)
            b = plus(k,times(-1,c))
            assert b == ambient_b
            # Projected ambient Hessian is recovered by subtracting the
            # embedding curvature q*K; normalization then gives B=-D.
            assert plus(pulled_hessian,times(-q,k)) == times(q,d)
            complete.append(fref(k,c,hi,alpha,beta))
            bare.append(fref(k,times(0,c),hi,alpha,beta))
            mixed.append(complete[-1]-bare[-1])
            p = response(b,hi,alpha,beta)
            momenta.append(p); ks.append(k); cs.append(c)
            velocity = times(2*lapse,k)
            step=Q(1,101)
            for e in basis:
                # Hold h,N,shift,T and independent clock jet fixed; vary hdot.
                kp=times(1/(2*lapse),plus(velocity,times(step,e)))
                km=times(1/(2*lapse),plus(velocity,times(-step,e)))
                derivative=lapse*(fref(kp,c,hi,alpha,beta)-fref(km,c,hi,alpha,beta))/(2*step)
                assert derivative==pair(p,e)
                velocity_checks+=1
        assert momenta[0]==momenta[1]==momenta[2]
        assert complete[0]==complete[1]==complete[2]
        expected=2*s*s*(alpha*norm(shape,hi)+beta*trace(shape,hi)**2)
        assert delta(bare)==expected and delta(mixed)==-expected and delta(complete)==0
        assert expected=={'rank_one':Q(11,125),'trace_free':Q(6,25)}[name]
        # The C contribution to the partial metric momentum cancels the
        # spurious change inferred by transporting the clock-adapted K term.
        naked_change=plus(response(ks[2],hi,alpha,beta),times(-1,response(ks[1],hi,alpha,beta)))
        omitted_change=plus(times(-1,response(cs[2],hi,alpha,beta)),response(cs[1],hi,alpha,beta))
        assert plus(naked_change,omitted_change)==matrix([[0]*3 for _ in range(3)])
        cases.append({'shape':name,'bare_density_second_difference':str(delta(bare)),
                      'mixed_density_second_difference':str(delta(mixed)),
                      'full_density_second_difference':str(delta(complete)),
                      'reference_partial_momentum_change':'0'})
    # Repeat the published orthonormal coefficient certificate independently.
    ap,bp,sp=Q(2,5),Q(-1,7),Q(1,2)
    assert 2*sp*sp*(ap+bp)==Q(9,70)
    assert 2*sp*sp*(2*ap)==Q(2,5)
    apparent_a=1/(1+ap)
    apparent_lagrange=(1-bp)/(1+ap)
    apparent_hamilton=apparent_lagrange/(3*apparent_lagrange-1)
    assert (apparent_a,apparent_lagrange,apparent_hamilton)==(Q(5,7),Q(40,49),Q(40,71))
    return {'status':'PASS_INTERNAL_NONBLIND_RANGE_REVIEW','numbered_round':None,
            'method':'Independent ambient embedding/pullback and hdot variations; exact Fraction arithmetic; no author imports',
            'groups':{
                'fixed_reference_projection':{'pass':True,'metric_determinant':'1','clock_gradient':str(q),'shape_cases':cases},
                'partial_metric_velocity':{'pass':True,'independent_hdot_direction_checks':velocity_checks,'lapse':str(lapse),'same_reference_partial_momentum':True},
                'published_coefficient_crosscheck':{'pass':True,'bare_second_differences':['9/70','2/5'],
                    'clock_adapted_a':str(apparent_a),'clock_adapted_lambda_L':str(apparent_lagrange),'clock_adapted_lambda_H':str(apparent_hamilton)}},
            'not_certified':['On-shell EOM or Noether conservation for these local jets',
                'Full canonical reduction or Hamiltonian constraint cancellation',
                'Original quantum-process implementation, stability or empirical calibration',
                'Generated EH/GR, new axiom, dark matter or dark energy identification']}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=run()
    target=Path(__file__).with_name('review_checks.json')
    if args.write:
        with target.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,indent=2,ensure_ascii=False)
            stream.write('\n')
    else:
        assert json.loads(target.read_text(encoding='utf-8-sig'))==result
    print('PASS: independent local reference jet, action and partial metric velocity audit; unnumbered')