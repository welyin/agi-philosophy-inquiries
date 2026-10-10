"""Independent 1054 review: rational bounds and moment-recursion dynamics.

No author module is imported. Composite Simpson radial moments and an angular
power series replace the author's Gauss quadrature and star diagonalization.
Numerical agreement is diagnostic, not the continuous theorem or QED error.
Default compares this review's saved result; --write-results creates it once.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'review_independent_results.json'


def run():
    t, r, g = Q(1, 100), Q(1, 4), Q(1, 10)
    z = 2*Q(15, 7)*t
    fp = Q(7, 12)*(2*Q(23, 90)*r+Q(23, 180)*r)
    dremainder = Q(2, 15)*z**4/(6*(1-z))
    slope = g*g*t*t*fp-dremainder
    premainder = z**4/(24*(1-z))
    plower = g*g*t*t*Q(49, 96)-premainder
    hessian = 2*(2*t/7)**2+2*(t/7)*(4*t/7+4*t*t/49)
    jvariation = 6*hessian/Q(100000)
    scalars = dict(minus_F_prime_lower=fp, derivative_remainder_upper=dremainder,
                   exact_slope_lower=slope, probability_remainder_upper=premainder,
                   exact_probability_lower=plower, hessian_component_bound=hessian,
                   jacobian_perturbation_upper=jvariation,
                   direction_t1_gap_lower_in_ball=Q(6889,4000000)-Q(8,700000))
    saved = json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    for k, v in scalars.items():
        assert str(v) == saved['rational_certificate']['exact'][k], k
    assert slope-jvariation > Q(1,200000000)
    assert plower > Q(1,3000000)

    # Integrate only physical radial moments. No normalized quadrature repair.
    ngrid, order = 1024, 18
    k = np.linspace(1., 2., ngrid+1)
    sw = np.ones(ngrid+1)
    sw[1:-1:2], sw[2:-1:2] = 4., 2.
    w = sw/(3*ngrid)*(8/3)*np.sin(np.pi*(k-1))**4
    bt = np.zeros_like(k)
    bl = np.zeros_like(k)
    dt = np.zeros_like(k)
    dl = np.zeros_like(k)
    for m in range(13):
        common = (-1)**m*(k*float(r))**(2*m)/math.factorial(2*m)
        at = .5*(1/(2*m+1)+1/(2*m+3))
        al = 1/(2*m+1)-1/(2*m+3)
        bt += at*common
        bl += al*common
        if m:
            dt += at*common*(2*m/float(r))
            dl += al*common*(2*m/float(r))
    moment = np.array([np.dot(w, k**n) for n in range(order)])
    injection = np.array([[np.dot(w*k**n, b) for b in (bt,bl)] for n in range(order)])
    injection_d = np.array([[np.dot(w*k**n, b) for b in (dt,dl)] for n in range(order)])
    cn, dn = np.zeros((order+1,2)), np.zeros((order+1,2))
    for n in range(order):
        cn[n+1] = 1.5*cn[n]+float(g)*injection[n]
        dn[n+1] = 1.5*dn[n]+float(g)*injection_d[n]
        for j in range(n):
            cn[n+1] += 2*float(g*g)*moment[n-1-j]*cn[j]
            dn[n+1] += 2*float(g*g)*moment[n-1-j]*dn[j]
    timeweights = np.array([(-1j*float(t))**n/math.factorial(n) for n in range(order+1)])
    c, dc = timeweights@cn, timeweights@dn
    p = float((2*abs(c[0])**2+abs(c[1])**2)/2)
    dp = float(2*(c[0].conjugate()*dc[0]).real+(c[1].conjugate()*dc[1]).real)
    post_z = float((abs(c[1])**2-2*abs(c[0])**2)/(2*p))
    sx = np.array([[0,1],[1,0]],complex)
    sy = np.array([[0,-1j],[1j,0]],complex)
    sz = np.diag([1.,-1.]).astype(complex)
    ke = np.vstack([c[0]*sx,c[0]*sy,c[1]*sz])/np.sqrt(2)
    bell = np.array([1,0,0,1],complex)/np.sqrt(2)
    v = (np.kron(ke,np.eye(2))@bell).reshape(3,2,2)
    reference = np.einsum('jar,jas->rs',v,v.conj())
    pure = np.outer(v.ravel(),v.ravel().conj())/p
    dephased = np.zeros_like(pure)
    for j in range(3):
        dephased[4*j:4*j+4,4*j:4*j+4] = pure[4*j:4*j+4,4*j:4*j+4]
    coherence_distance = float(np.sum(abs(np.linalg.eigvalsh(pure-dephased)))/2)
    comparisons = {
        'p': (p,saved['probability_and_slope']['p']),
        'minus_p_prime': (-dp,saved['probability_and_slope']['minus_p_prime']),
        'conditional_excited_spin_z': (post_z,saved['complete_instrument_diagnostic']['conditional_excited_spin_z'])}
    differences = {name: abs(a-b) for name,(a,b) in comparisons.items()}
    assert differences['p'] < 1e-18 and differences['minus_p_prime'] < 1e-18
    assert differences['conditional_excited_spin_z'] < 1e-12
    assert np.linalg.norm(reference-p*np.eye(2)/2) < 1e-20
    assert np.linalg.norm(ke.conj().T@ke-p*np.eye(2)) < 1e-20
    assert coherence_distance > .6
    return {
        'round':1054, 'independent_check_passed':True,
        'method':'angular moment series + unrenormalized composite Simpson + Hamiltonian block-moment recurrence',
        'rational_checks':{k:str(v) for k,v in scalars.items()},
        'radial_normalization':float(w.sum()), 'p':p, 'minus_p_prime':-dp,
        'conditional_excited_spin_z':post_z,
        'source_coherent_vs_dephased_Bell_output_distance':coherence_distance,
        'Bell_reference_residual':float(np.linalg.norm(reference-p*np.eye(2)/2)),
        'differences_from_author':differences,
        'source_tag_coherence_is_retained_not_new_access_permission':True,
        'continuous_theorem_or_full_QED_error_claimed_by_numerics':False,
        'reviewed_author_receipt_sha256':hashlib.sha256((HERE/'research_round_1054_checks.json').read_bytes()).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with OUT.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
    else:
        previous = json.loads(OUT.read_text(encoding='utf-8'))
        for key in ('rational_checks','method','reviewed_author_receipt_sha256'):
            assert result[key] == previous[key]
        for key in ('p','minus_p_prime','conditional_excited_spin_z','source_coherent_vs_dephased_Bell_output_distance'):
            assert np.isclose(result[key],previous[key],rtol=1e-11,atol=1e-18)
    print(json.dumps(result,ensure_ascii=False))
