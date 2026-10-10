"""Count-0 algebra diagnostics for adopted normal-contact transport.

No experimental fitting. e=h=k_B=1 except the explicitly labelled SI values.
Default is read-only; --write-results exclusively creates the first result file.
"""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
import argparse
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent


def quadrature(order=24):
    x, w = np.polynomial.legendre.leggauss(order)
    centers = np.arange(-59.0, 60.0, 2.0)
    return (centers[:, None] + x).ravel(), np.tile(w, len(centers))


def tau(energy):
    # One spin-resolved channel, smooth and everywhere in [0.1,0.9].
    return 0.1 + 0.8 / (1.0 + ((energy-0.4)/0.7)**2)


def fermi(energy, mu, temperature):
    return 1.0/(1.0+np.exp(np.clip((energy-mu)/temperature, -700, 700)))


def currents(mu_l, t_l, mu_r, t_r, order=24, zero_shift=0.0):
    grid, weights = quadrature(order)
    energy = grid + zero_shift
    trans = tau(grid)
    fl = fermi(energy, mu_l+zero_shift, t_l)
    fr = fermi(energy, mu_r+zero_shift, t_r)
    density = trans*(fl-fr)
    jn = float(weights @ density)
    je = float(weights @ (energy*density))
    pl = float(weights @ ((mu_l+zero_shift-energy)*density))
    pr = float(weights @ ((energy-mu_r-zero_shift)*density))
    noise = float(2*weights @ (trans*(fl*(1-fl)+fr*(1-fr))
                               +trans*(1-trans)*(fl-fr)**2))
    entropy = pl/t_l+pr/t_r
    affinity = (energy-mu_r-zero_shift)/t_r-(energy-mu_l-zero_shift)/t_l
    assert np.min(density*affinity) >= -1e-15
    assert abs(pl+pr-(mu_l-mu_r)*jn) < 1e-12
    assert noise >= 0 and entropy >= -1e-12
    return dict(jn=jn, je=je, p_left=pl, p_right=pr,
                noise=noise, entropy=entropy)


def equilibrium(mu, temperature, order=24):
    energy, weights = quadrature(order)
    f = fermi(energy, mu, temperature)
    kernel = tau(energy)*f*(1-f)/temperature
    moments = [float(weights @ ((energy-mu)**j*kernel)) for j in range(3)]
    l0, l1, l2 = moments
    conductance = l0
    k_open = (l2-l1*l1/l0)/temperature
    k_fixed_v = l2/temperature
    noise = float(4*weights @ (tau(energy)*f*(1-f)))
    assert k_open >= 0 and k_fixed_v >= k_open
    assert abs(noise-4*temperature*conductance) < 1e-12
    return dict(moments=moments, conductance=conductance,
                k_open=k_open, k_fixed_v=k_fixed_v, noise=noise)


def calculate():
    # Exact zero-temperature window [-1/2,1/2], tau=1/2+E/4.
    jn, je = F(1,2), F(1,48)
    pl, pr = jn/2-je, je+jn/2
    noise = 2*(F(1,4)-F(1,16)*F(1,12))
    assert (pl, pr, noise)==(F(11,48), F(13,48), F(47,96))
    assert pl+pr==jn and pl!=pr
    # Same conductance and constant-energy thermal moment; different noise.
    channels = [(F(1), F(0)), (F(1,2), F(1,2))]
    assert [sum(v) for v in channels]==[F(1),F(1)]
    partitions=[sum(t*(1-t) for t in v) for v in channels]
    assert partitions==[F(0), F(1,2)]

    eq=[equilibrium(mu,t) for mu,t in [(0.,.7),(.2,1.0)]]
    finite=[]
    max_refinement=0.0
    for pars in [(.4,1.0,-.3,1.0),(.2,1.2,-.1,.8),(-.2,.8,.3,1.1)]:
        a=currents(*pars)
        b=currents(*pars,order=32)
        shifted=currents(*pars,zero_shift=7.0)
        for key in a:
            max_refinement=max(max_refinement,abs(a[key]-b[key]))
        assert max_refinement < 1e-11
        for key in ('jn','p_left','p_right','noise','entropy'):
            assert abs(a[key]-shifted[key])<1e-12
        assert abs(shifted['je']-a['je']-7*a['jn'])<1e-12
        finite.append({'input_muL_TL_muR_TR':list(pars),'output':a})

    # Constant-transmission Sommerfeld moments, numerical normalization check.
    energy, weights=quadrature()
    f=fermi(energy,0.0,1.0)
    m=[float(weights @ (energy**j*f*(1-f))) for j in range(3)]
    assert abs(m[0]-1)<1e-12 and abs(m[1])<1e-12
    assert abs(m[2]-math.pi**2/3)<1e-11
    with localcontext() as ctx:
        ctx.prec=45
        d=Decimal
        e,h,k=d('1.602176634e-19'),d('6.62607015e-34'),d('1.380649e-23')
        pi=d('3.141592653589793238462643383279502884197169399')
        si={'T_kelvin':'305','G_spin_pair_siemens':str(2*e*e/h),
            'K_spin_pair_watt_per_kelvin':str(2*pi*pi*k*k*d(305)/(3*h)),
            'Lorenz_SI':str(pi*pi*k*k/(3*e*e))}
    return dict(count_new_science=0, all_checks_passed=True,
                exact_zero_T={'JN':str(jn),'JE':str(je),'PL':str(pl),
                              'PR':str(pr),'SII':str(noise)},
                equal_G_partition_sums=[str(x) for x in partitions],
                equilibrium=eq, finite_bias=finite,
                constant_tau_fermi_moments=m,
                quadrature_refinement_max_difference=max_refinement,
                SI_normalization_only=si,
                scope='Adopted-formula diagnostics, not experiment or rigorous continuum error certificate.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-results', action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=HERE/'results.json'
    if args.write_results:
        with path.open('x',encoding='utf-8') as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write('\n')
    else:
        saved=json.loads(path.read_text(encoding='utf-8'))
        def compare(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for key in a: compare(a[key],b[key])
            elif isinstance(a,list):
                assert len(a)==len(b)
                for x,y in zip(a,b): compare(x,y)
            elif isinstance(a,float):
                assert abs(a-b)<1e-11
            else:
                assert a==b
        compare(result,saved)
    print(json.dumps({'passed':True,'new_science':0,
                      'quadrature_refinement':result['quadrature_refinement_max_difference']}))
