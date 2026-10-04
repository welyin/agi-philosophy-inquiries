"""715: reference-weighted readout and the energy-domain obstruction.

The full interacting Gibbs asymptotic is analytic. Numerical integrals use
the explicit original Haar Gauss wavefunction, not a surrogate Gibbs state.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_holonomy_readout_scale as old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_sharp_record_domain_results.json'


@lru_cache(None)
def rule(n,left,right):
    x,w=np.polynomial.legendre.leggauss(n)
    return left+(right-left)*(x+1)/2,(right-left)*w/2


def haar_weighted_density(x,n=192):
    a=np.abs(np.asarray(x))[...,None]
    theta,w=rule(n,0.,float(np.pi/2))
    s=np.sin(theta);c=np.cos(theta)
    integral=(1-a*a)**2*np.sum(w*s*s*c*c/np.sqrt(a*a+(1-a*a)*s*s),axis=-1,keepdims=True)
    return np.squeeze(4*old.CASIMIR/np.pi**2*integral,axis=-1)


def density_check():
    x,w=rule(192,0.,1.)
    g=haar_weighted_density(x)
    moments=[float(2*np.sum(w*g*x**k)) for k in (0,2)]
    expected=[old.CASIMIR/8,old.CASIMIR/64]
    error=max(abs(a-b) for a,b in zip(moments,expected))
    g0=4*old.CASIMIR/(3*np.pi**2)
    assert error<1e-10 and abs(float(haar_weighted_density(0.))-g0)<1e-13
    return dict(weighted_Casimir=old.CASIMIR,haar_weighted_density_at_zero=g0,
                equivalent_original_couplings_expression=(old.BW+12*old.B0)/np.pi**2,
                computed_moments=moments,exact_moments=expected,moment_error=error)


def kernel(x,eta):
    return 1/np.cosh(x)**4/(1-eta*eta*np.tanh(x)**2)


def mean_injection(z,n=160):
    eta=.6;kxN=4/.73*np.exp(-.24)
    x,w=rule(n,0.,float(min(z,18.)))
    return float(kxN*eta**2*z/2*np.sum(w*kernel(x,eta)*haar_weighted_density(x/z,n)))


def direct_haar(z,n=192):
    theta=np.pi*np.arange(1,n+1)/(n+1)
    u=np.cos(theta)[:,None]
    weights=(2/(n+1)*np.sin(theta)**2)[:,None]
    phase=2*np.pi*(np.arange(4*n)+.5)/(4*n)
    f=u*np.cos(phase)
    gam=old.BW*(1-u*u)*np.cos(phase)**2/4+9*old.B0*u*u*np.sin(phase)**2
    return float(4/.73*np.exp(-.24)*.6**2*z*z/4*np.sum(weights*gam*kernel(z*f,.6))/(4*n))


def asymptotic_check():
    eta=.6
    I=2/eta**2-2*(1-eta**2)*np.arctanh(eta)/eta**3
    x,w=rule(192,0.,18.)
    kernel_error=abs(2*float(np.sum(w*kernel(x,eta)))-I)
    limit=4/.73*np.exp(-.24)*eta**2/4*(4*old.CASIMIR/(3*np.pi**2))*I
    rows=[]
    for z in (1.,2.,4.,8.,16.,32.,64.,128.):
        lo=mean_injection(z,128);hi=mean_injection(z,256)
        rows.append(dict(gain=z,mean_injection_in_actual_Haar_state=hi,
                         divided_by_gain=hi/z,relative_to_analytic_slope=hi/z/limit,
                         two_quadratures_difference=abs(lo-hi),
                         conformal_mean_source_change=-2*hi))
    direct_errors=[abs(direct_haar(z)-mean_injection(z,256)) for z in (1.,2.,4.)]
    assert max(direct_errors)<1e-9 and kernel_error<1e-11
    assert max(row['two_quadratures_difference'] for row in rows)<1e-8
    assert abs(rows[-1]['relative_to_analytic_slope']-1)<.002
    return dict(kernel_integral=I,kernel_quadrature_error=kernel_error,
                exact_Haar_state_linear_slope=limit,rows=rows,
                direct_original_group_integral_errors=direct_errors,
                original_interacting_Gibbs_slope_not_numerically_computed=True)


def branch_l2(z,n):
    theta=np.pi*np.arange(1,n+1)/(n+1)
    u=np.cos(theta)[:,None]
    weights=(2/(n+1)*np.sin(theta)**2)[:,None]
    phase=2*np.pi*(np.arange(4*n)+.5)/(4*n)
    f=u*np.cos(phase)
    smooth=np.sqrt((1+.6*np.tanh(z*f))/2)
    limit=np.sqrt((1+.6*np.sign(f))/2)
    return float(np.sum(weights*(smooth-limit)**2)/(4*n))


def state_limit_check():
    rows=[]
    for z in (2.,8.,32.):
        a=branch_l2(z,192);b=branch_l2(z,384)
        rows.append(dict(gain=z,original_Haar_branch_L2_error_squared=b,
                         quadrature_order_difference=abs(a-b),
                         cq_trace_distance_upper_from_vector_bound=2*np.sqrt(2*b)))
    assert all(rows[i+1]['original_Haar_branch_L2_error_squared']<rows[i]['original_Haar_branch_L2_error_squared'] for i in range(2))
    eta=.6
    jump=np.sqrt((1+eta)/2)-np.sqrt((1-eta)/2)
    assert jump>0
    return dict(eta=eta,nonzero_sharp_Kraus_jump=jump,rows=rows,
                numerical_state_is_normal_original_Gauss_Haar_wavefunction=True,
                infinite_energy_of_full_Gibbs_sharp_state_is_analytic=True,
                finite_gain_instrument_not_an_autonomous_implementation=True)


def run():
    deps=('research_note_589.md','research_note_598.md','research_note_603.md','research_note_623.md',
          'research_note_633.md','research_note_638.md','research_note_643.md','research_note_713.md',
          'research_note_714.md','round715_drafts/observable_gain_entry.md')
    return dict(round=715,tests_run=3,failures=0,errors=0,
                weighted_density=density_check(),reference_mean=asymptotic_check(),
                state_limit=state_limit_check(),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Fixed original finite graph and positive geometry. Full interacting physical Gibbs '
                      'linear sharpening cost and infinite-energy limit are analytic; numerical coefficients '
                      'use an explicitly identified normal Haar Gauss state. No joint spatial continuum '
                      'no-go, autonomous detector or GR derivation.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
