"""Gauge threshold test for the same scalar representation used in the type-II mass term."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SM = np.array([41/10, -19/6, -7])
TRIPLET_FERMION = np.array([4/15, 0, 2/3])
DOUBLET_FERMION = np.array([2/5, 2/3, 0])
SCALARS = np.array([[3/5, 2/3, 0], [1/30, 1/2, 1/3], [8/15, 0, 5/6]])
READOUT = np.array([[-1, 1, 0], [0, -1, 1]])


def read(name):
    return json.loads((HERE / name).read_bytes())


def scalar_ratios(gauge, k, length):
    # Unknowns X=sum log(m_Delta/m_B), Y=sum log(m_S/m_B).
    # The common scalar log relative to B drops out of the two differences.
    a = 4*math.pi/np.square(gauge)
    residual = READOUT @ (2*math.pi*a - SM*length + 3*TRIPLET_FERMION*math.log(k))
    return np.linalg.solve(-READOUT @ SCALARS[[0,2]].T, residual)


def explicit_model(cfg, measured, k, ratios, lower):
    upper = cfg['matching_GeV']
    witness = cfg['spectrum_witness']
    mean = witness['doublet_geometric_mean_GeV']
    small = witness['three_identical_doublet_masses_GeV']
    doublets = [small]*3 + [mean*(mean/small)**3]
    delta, sextet = mean*np.exp(ratios/4)
    masses = np.column_stack([np.full(4,delta), doublets, np.full(4,sextet)])
    assert masses.min()>=lower and masses.max()<=upper

    # Three paired Dirac fermions realizing the same measured gain per family.
    m_l = witness['ordinary_heavy_lepton_GeV']
    m_d = k*m_l
    inverse = 4*math.pi/np.square(measured)
    corrections = (SM*math.log(upper/lower)
                   + 3*TRIPLET_FERMION*math.log(upper/m_d)
                   + 3*DOUBLET_FERMION*math.log(upper/m_l)
                   + np.log(upper/masses).sum(axis=0) @ SCALARS)/(2*math.pi)
    unified = inverse-corrections
    np.testing.assert_allclose(unified,np.full(3,unified.mean()),rtol=0,atol=2e-12)
    assert unified.min()>0
    return {'doublet_masses_GeV':doublets,'common_triplet_mass_GeV':delta,
            'common_sextet_mass_GeV':sextet,'fermion_doublet_GeV':m_l,
            'fermion_triplet_GeV':m_d,'common_inverse_alpha':unified.mean(),
            'maximum_inverse_coupling_difference':float(np.ptp(unified)),
            'alternative_neutrino_active_doublet_GeV':[doublets[0],doublets[-1]]}


def compute():
    cfg = read('inputs.json')
    g = read(cfg['gauge'])
    k_bounds = read(cfg['flavour_bounds'])['two_width_flavour_gain']
    # Deliberately round bounds outwards before the full-box proof.
    k_box = [math.floor(k_bounds[0]*1e5)/1e5, math.ceil(k_bounds[1]*1e3)/1e3]
    k_point = read(cfg['overlap'])['concrete_shared_overlap']['k']
    gauge_rows = [g[key] for key in ('g1_GUT','g2','g3')]
    gauge_box = [[v-cfg['widths']*e,v+cfg['widths']*e] for v,e in gauge_rows]
    upper, lower = cfg['matching_GeV'], g['mu_ref_GeV']
    length = math.log(upper/lower)
    extrema = np.array([scalar_ratios(gp,k,length) for gp in itertools.product(*gauge_box)
                        for k in k_box])
    lo, hi = extrema.min(axis=0), extrema.max(axis=0)
    assert 0 < lo[1] < hi[1] < lo[0] < hi[0] < 4*length
    assert lo[0]>3*length
    bounds = [lower*math.exp((4*length-lo[0])/(5-index)) for index in range(1,5)]
    centre = [v for v,_ in gauge_rows]
    ratios = scalar_ratios(centre,k_point,length)
    witness = explicit_model(cfg,centre,k_point,ratios,lower)
    return {'effective_gain_box':k_box,'log_mass_range':length,
            'sum_log_triplet_over_doublet':[lo[0],hi[0]],
            'sum_log_sextet_over_doublet':[lo[1],hi[1]],
            'least_number_of_scalar_15s':4,
            'four_copies_upper_GeV_on_ordered_doublets':bounds,
            'four_copies_geometric_mean_upper_GeV':upper*math.exp(-lo[0]/4),
            'central_gain_and_log_ratios':[k_point,*ratios],
            'central_four_copy_geometric_mean_upper_GeV':upper*math.exp(-ratios[0]/4),
            'central_witness':witness}


def tidy(data):
    if isinstance(data,dict):
        return {key:tidy(value) for key,value in data.items()}
    if isinstance(data,list):
        return [tidy(value) for value in data]
    if isinstance(data,(float,np.floating)):
        return 0.0 if abs(data)<1e-13 else float(f'{data:.11g}')
    return data


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=tidy(compute())
    if args.check:
        assert result==read('results.json')
        print('PASS: both gauge readouts, exact count, mass bounds, positive common coupling and explicit spectrum')
    else:
        with (HERE/'results.json').open('x',encoding='utf-8') as outfile:
            outfile.write(json.dumps(result,indent=2)+'\n')
