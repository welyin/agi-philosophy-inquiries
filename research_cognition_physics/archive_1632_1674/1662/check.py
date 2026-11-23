"""Join a family-universal overlap to the threshold of the same heavy Dirac matter."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_bytes())


def ends(row, widths):
    return row[0]-widths*row[1], row[0]+widths*row[-1]


def gauge_comparison(cfg, gain):
    g = read(cfg['gauge'])
    b = np.array([41/10, -19/6, -7])
    slope = (b[2]-b[1])/(b[0]-b[1])
    coefficient = (1-.4*slope)/math.pi
    data = []
    scales = []
    invariant_weights = np.array([5, -3, -2])
    triplet = np.array([4/15, 0, 2/3])
    doublet = np.array([2/5, 2/3, 0])
    np.testing.assert_allclose(np.array([triplet, doublet, b]) @ invariant_weights,
                               [0, 0, 44], atol=1e-14)
    for row in itertools.product(*(ends(g[k], cfg['widths']) for k in ('g1_GUT','g2','g3'))):
        a1,a2,a3 = 4*math.pi/np.square(row)
        mismatch = a3-a2-slope*(a1-a2)
        data.append(mismatch)
        running_length = 2*math.pi*(5*a1-3*a2-2*a3)/44
        scales.append(g['mu_ref_GeV']*math.exp(running_length))
    u = [min(data),max(data)]
    predicted = [-coefficient*math.log(k) for k in gain[::-1]]
    required = [math.exp(-point/coefficient) for point in u[::-1]]
    inverse_g3 = []
    for g1,g2,k in itertools.product(ends(g['g1_GUT'],cfg['widths']),ends(g['g2'],cfg['widths']),gain):
        a1,a2 = 4*math.pi/np.square([g1,g2])
        inverse_g3.append(a2+slope*(a1-a2)-coefficient*math.log(k))
    g3 = sorted([math.sqrt(4*math.pi/x) for x in (min(inverse_g3),max(inverse_g3))])
    measured_g3 = list(ends(g['g3'],cfg['widths']))
    assert max(u)<min(predicted) and g3[1]<measured_g3[0]
    assert max(scales)<float(cfg['matching'])
    return {'U_from_gauge_data':u,'U_from_flavour_bounds':predicted,
            'unaccounted_projection_at_least':predicted[0]-u[1],
            'gain_required_by_gauge_sector':required, 'g3_holdout_envelope':g3,
            'g3_empirical_box':measured_g3,
            'scale_GeV_from_gauge_alone':[min(scales),max(scales)]}


def threshold_test(k):
    # A nonsymmetric bare block including two spectators. Right row mixing is
    # chosen by an isometry so its norm changes every ordinary direction equally.
    rng = np.random.default_rng(1662)
    M = np.diag([1.1,2.2,3.4,4.8,6.1])+0.06*rng.normal(size=(5,5))
    raw = rng.normal(size=(5,3))+1j*rng.normal(size=(5,3))
    embedding,_ = np.linalg.qr(raw)
    def spectrum(bare,tensor):
        mix = bare@tensor
        stack = np.column_stack([mix,bare])
        heavy = np.linalg.svd(stack,compute_uv=False)
        light,_ = np.linalg.qr(np.vstack([np.eye(3),-tensor]))
        np.testing.assert_allclose(stack@light,0,atol=2e-14)
        return heavy, light[:3]@light[:3].conj().T

    ell = .9
    T_L = math.sqrt(ell**-2-1)*embedding
    T_D = math.sqrt((k/ell)**2-1)*embedding.conj()
    d, Pd = spectrum(M,T_D)
    l, Pl = spectrum(M,T_L)
    np.testing.assert_allclose(Pl,k*k*Pd,atol=2e-14)
    ratio = d.prod()/l.prod()
    np.testing.assert_allclose(ratio,k**3)
    deform = np.diag([1,1,1,1,2.5])
    altered, Pnew = spectrum(deform@M,T_D)
    np.testing.assert_allclose(Pnew,Pd,atol=2e-14)
    np.testing.assert_allclose(altered.prod()/l.prod(),ratio*2.5)
    return {'k':k,'squared_D_overlap':float(Pd[0,0].real),
            'squared_L_overlap':float(Pl[0,0].real),'heavy_product_ratio':ratio,
            'independent_bare_splitting_example':altered.prod()/l.prod()}


def compute():
    cfg = read('inputs.json')
    rows = read(cfg['spectra'])['scales'][cfg['matching']]
    d,s,b,e,mu,tau = (ends(x, cfg['widths']) for x in rows)
    theta = read(cfg['ckm'])['high_scale']['angles_theta12_theta13_theta23_delta']
    t13, t23 = (ends(theta[j],cfg['widths']) for j in (1,2))
    # Exact monotonic endpoints on the positive-angle intervals; the inequalities
    # discard only positive row-norm contributions.
    vtb_min = math.cos(t13[1])*math.cos(t23[1])
    vub_min = math.sin(t13[0])
    ca,cb,_=cfg['clebsch']
    k_min = math.sqrt(e[0]*mu[0]*vtb_min/(abs(ca*cb)*d[1]*s[1]))
    k_max = tau[1]/(abs(cb)*b[0]*vub_min)
    assert k_min>1 and k_max>k_min
    witness=read(cfg['charged_witness'])['extended_box_GJ_test']
    k=witness['k_times_yb']/witness['marginal_coordinates'][2]
    assert k_min<k<k_max
    return {'necessary_k_range':[k_min,k_max],'generations':3,
            'concrete_shared_overlap':threshold_test(k),
            'gauge_comparison':gauge_comparison(cfg,[k_min,k_max])}


def tidy(x):
    if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):return [tidy(v) for v in x]
    return float(f'{x:.10g}') if isinstance(x,(float,np.floating)) else x


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read('results.json')
        print('PASS: mass and gauge windows, universal overlaps, heavy spectrum, competitor')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as stream:
            stream.write(json.dumps(result,indent=2)+'\n')
