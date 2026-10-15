"""Independent 1221 audit using closed thermal spectra, not the author module.

The original flux compression is two dimensional; its original sterile CAR
part has two copies of the one-particle energies (+t,-t). This permits a
separate calculation without diagonalising the author's 32-by-32 matrix.
This checks a physical compression only, not a full-model Gibbs approximation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]


def h2(p):
    return -p * math.log(p) - (1-p) * math.log1p(-p)


def coefficients(s=.12, t=.8):
    source = BASE/'archive_554_584/569/joint_quotient_gauge_completion_results.json'
    old = json.loads(source.read_text('utf-8'))['evidence']['units_check']
    S = np.array([[.24,.31,-.12],[.31,-.07,.17],[-.12,.17,-.17]])
    # Direct exponential series, independent of the author's eigh exponential.
    A = t*S
    E = np.eye(3)
    term = np.eye(3)
    for n in range(1,31):
        term = term@A/n
        E += term
    pref = math.exp(-2*s)/.73
    b = old['electric_coefficients_colour_weak_integer_circle'][1]
    a = old['residual_positive_coefficients'][1]
    k = b*pref*(2*E[0,0]+2*E[1,1]-2*E[0,1])
    g = a*pref*E[2,2]
    dE = S@E
    dg = a*pref*dE[2,2]
    return k,g,dg


def closed(s=.12, t=.8):
    beta, hopping = 1.7,.17
    k,g,dg = coefficients(s,t)
    center = k-g/2
    radius = math.hypot(center,g)
    th = math.tanh(beta*radius)
    p = .5+.5*center/radius*th
    off = .5*g/radius*th
    flux_entropy = h2((1+th)/2)
    coherence = h2(p)-flux_entropy
    logz_flux = -beta*center+math.log(2*math.cosh(beta*radius))
    reverse_flux = -h2(p)+beta*(1-p)*2*center+logz_flux
    q = 1/(1+math.exp(beta*hopping))
    car_entropy = 4*h2(q)
    mi = 4*math.log(2)-car_entropy
    car_energy = -2*hopping*math.tanh(beta*hopping/2)
    reverse_car = 4*math.log(math.cosh(beta*hopping/2))
    gap = 2*g*off-car_energy
    c = coherence+mi
    reverse = reverse_flux+reverse_car
    source_s = -4*g*off
    source_t = 2*dg*off
    result = dict(correlation_relative_entropy=c, sector_coherence=coherence,
        conditional_multiplicity_mutual_information=mi,
        reverse_relative_entropy=reverse, boundary_energy_gap=gap,
        sector_probabilities=[p,1-p], center_Shannon=h2(p),
        one_side_forced_matching_entropy=(1-p)*math.log(9),
        extended_mutual_information=h2(p)+2*(1-p)*math.log(9)+c,
        naive_extended_product_outside_matching_weight=1-p*p-(1-p)**2/81,
        fixed_state_shear_source_gap=source_t,
        fixed_state_conformal_source_gap=source_s)
    assert abs(c+reverse-beta*gap)<3e-14
    return result


def parity_and_cut():
    # Independently build the four original neutral modes in bit order A,A,B,B.
    dim=16
    H=np.zeros((dim,dim))
    for state in range(dim):
        for spin in (0,1):
            for src,dst in ((spin,spin+2),(spin+2,spin)):
                i,j=3-src,3-dst
                if state&(1<<i) and not state&(1<<j):
                    sign=(-1)**((state>>(i+1)).bit_count())
                    removed=state^(1<<i)
                    sign*=(-1)**((removed>>(j+1)).bit_count())
                    H[removed^(1<<j),state]+=.17*sign
    ev,V=np.linalg.eigh(H)
    w=np.exp(-1.7*ev);w/=sum(w)
    rho=(V*w)@V.T
    P=np.diag([(-1)**j.bit_count() for j in range(dim)])
    assert np.linalg.norm(P@rho-rho@P)<1e-14
    a=np.zeros((4,4));b=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            a[i,j]=sum(rho[4*i+k,4*j+k] for k in range(4))
            b[i,j]=sum(rho[4*k+i,4*k+j] for k in range(4))
    margin=max(float(np.max(abs(a-np.eye(4)/4))),float(np.max(abs(b-np.eye(4)/4))))
    assert margin<1e-14
    # The matching vectors have ranks 1 and 9. For each, the product of
    # maximally mixed edge marginals has overlap 1/rank**2 with the singlet.
    errors=[]
    for rank in (1,9):
        singlet=np.eye(rank).reshape(-1)/math.sqrt(rank)
        errors.append(abs(float(singlet@singlet)-1))
        product=np.eye(rank*rank)/(rank*rank)
        errors.append(abs(float(singlet@product@singlet)-1/(rank*rank)))
    assert max(errors)<1e-14
    return dict(original_CAR_dimension=16, local_even_marginal_error=margin,
                matching_overlap_error=max(errors), full_thermal_state_claim=False)


def run():
    values=closed()
    saved=json.loads((HERE/'results.json').read_text('utf-8'))
    actual=saved['groups']['matched_entropy_energy']
    actual=dict(actual,**saved['groups']['original_boundary_sources'])
    error=max(float(np.max(abs(np.asarray(v)-np.asarray(actual[k])))) for k,v in values.items())
    assert error<2e-12, error
    h=2e-5
    derivative=(closed(.12+h)['boundary_energy_gap']-closed(.12-h)['boundary_energy_gap'])/(2*h)
    state=derivative-values['fixed_state_conformal_source_gap']
    assert abs(state)>.3
    assert abs(derivative-actual['moving_pair_total_energy_gap_derivative'])<2e-10
    assert abs(state-actual['indispensable_state_derivative'])<2e-10
    assert min(values['sector_coherence'],values['conditional_multiplicity_mutual_information'])>.01
    return dict(status='PASS',check_groups=3,closed_spectrum_values=values,
        author_comparison_maximum_error=error,parity_and_matching=parity_and_cut(),
        moving_pair_derivative=derivative,required_state_derivative=state,
        numerical_scope='Independent exact-compression diagnostic, not full Gauss Gibbs evaluation',
        author_code_imported=False,
        source_sha256={'archive_1219_/1221/check.py':hashlib.sha256((HERE/'check.py').read_bytes()).hexdigest(),
                       'archive_1219_/1221/results.json':hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest()},
        independent_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=run()
    dest=HERE/'independent_results.json'
    if args.write:
        dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif dest.exists():
        old=json.loads(dest.read_text('utf-8'))
        assert old['source_sha256']==result['source_sha256']
        assert old['independent_code_sha256']==result['independent_code_sha256']
        assert old['status']==result['status']=='PASS'
        for k,v in old['closed_spectrum_values'].items():
            assert np.max(abs(np.asarray(v)-np.asarray(result['closed_spectrum_values'][k])))<1e-12
    else:
        raise SystemExit('Missing independent_results.json; explicit --write required')
    print(json.dumps(result,ensure_ascii=False,indent=2))
