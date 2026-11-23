"""Independent verification with two CMS collision energies and both incoming partons."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# CMS AK7, central |y|. Signed absolute shifts in pb/GeV, in the two
# non-overlapping pT bins [737,846] and [1588,1784] GeV respectively.
CMS13 = (0.27377, 0.0012515)
ERROR13 = (
    ("all uncorr. unc. 1", 0.00090705, 0.000024274),
    ("AbsoluteMPFBias", 0.0041948, 0.000024465),
    ("AbsoluteScale", 0.0025374, 0.000015001),
    ("fake", 2.5336e-7, -1.1151e-8),
    ("FlavorQCD", 0.0014427, -0.000010599),
    ("Fragmentation", -0.00088046, -0.0000066991),
    ("JER", -0.0010376, -0.0000052552),
    ("miss", 0.0003273, 0.0000022216),
    ("PileUpDataMC", 0.00034058, 4.344e-7),
    ("PileUpPtBB", -0.00010479, -0.0000028333),
    ("PileUpPtEC1", -1.7069e-7, 8.8224e-9),
    ("PileUpPtEC2", -1.6157e-8, 1.8062e-10),
    ("PileUpPtRef", -0.00017244, 0.0000011427),
    ("Pref", -0.0010761, -0.00000304),
    ("PU", 0.00011803, 5.1263e-7),
    ("RelativeBal", 0.0031018, 0.000019695),
    ("RelativeFSR", 0.0018661, 0.000027223),
    ("RelativeJEREC1", 1.3346e-7, 3.0922e-8),
    ("RelativeJEREC2", -2.7034e-7, 1.0338e-8),
    ("RelativePtBB", 0.0015245, 0.000012022),
    ("RelativePtEC1", 1.5248e-7, 4.7922e-9),
    ("RelativePtEC2", -2.3319e-8, -1.5522e-10),
    ("RelativeSample", 0.000070329, -0.0000090765),
    ("RelativeStatEC", 3.2671e-7, 2.7969e-8),
    ("RelativeStatFSR", 0.00029801, 0.0000017596),
    ("SinglePionECAL", 0.00094717, 0.0000074182),
    ("SinglePionHCAL", 0.00050305, 0.0000048855),
    ("TimePtEta", 0.0076905, 0.000044665),
    ("Lumi", 0.0032852, 0.000015018),
    ("model", 0.00079342, 0.0000039183),
)
# (lo, hi, value, absolute statistical error, systematic down/up magnitude)
CMS8 = (
    (737, 790, 0.0659701, 0.0002567, 0.00249017, 0.00245152),
    (790, 846, 0.0398936, 0.0001939, 0.00153886, 0.00151361),
    (1588, 1784, 0.0000720384, 0.000004876, 0.00000464370, 0.00000488690),
)


def shape_double_ratio(h13, l13, h8, l8):
    return (h13/l13)/(h8/l8)


def cross_sections():
    (low13, high13) = CMS13
    pieces = np.asarray(CMS8)
    widths = pieces[:,1]-pieces[:,0]
    weights = widths[:2]/widths[:2].sum()
    centres = pieces[:,2]
    low8, high8 = weights@centres[:2], centres[2]
    syst13 = np.asarray([v[1:] for v in ERROR13])
    dx13 = np.linalg.norm(syst13,axis=0)
    dx8 = np.hypot(pieces[:,3:4],pieces[:,4:6])
    low8_errors = weights@dx8[:2]
    high8_errors = dx8[2]
    point = shape_double_ratio(high13,low13,high8,low8)
    box = [
        shape_double_ratio(high13-dx13[1],low13+dx13[0], high8+high8_errors[1],low8-low8_errors[0]),
        shape_double_ratio(high13+dx13[1],low13-dx13[0], high8-high8_errors[0],low8+low8_errors[1])]
    # Check the monotone corner calculation with every independent input sign.
    import itertools
    vertices=[]
    for signs in itertools.product((-1,1),repeat=5):
        l13,h13=CMS13+np.asarray(signs[:2])*dx13
        bins=centres+[s*dx8[i,int(s>0)] for i,s in enumerate(signs[2:])]
        vertices.append(shape_double_ratio(h13,l13,bins[2],weights@bins[:2]))
    np.testing.assert_allclose(box,[min(vertices),max(vertices)],rtol=1e-15)
    return {'13TeV_high_over_low':high13/low13,'8TeV_low_pb_per_GeV':low8,
            '8TeV_high_over_low':high8/low8,'double_ratio':point,
            'marginal_box_for_double_ratio':box}


def first_log():
    # Per weak doublet of complex fundamental scalars. For the leading
    # collinear shape z(1-z) the colour trace normalizes one tagged daughter.
    z,w = np.polynomial.legendre.leggauss(6)
    z=(z+1)/2
    w=w/2
    tr=0.5
    real_momentum = 4*float(w@(tr*z**2*(1-z)))
    delta_b=2*tr/3
    np.testing.assert_allclose(real_momentum,delta_b/2,rtol=1e-12)
    # Test p-n_g by differentiating the hard factor and individual beams.
    eps=1e-6
    slopes={}
    for power in [2,3]:
        row=[]
        for incoming_gluons in [2,1,0]:
            def rate(t):
                return (1+t)**power*(1-t)**incoming_gluons
            slope=(rate(eps)-rate(-eps))/(2*eps)
            assert abs(slope-(power-incoming_gluons)) < 1e-8
            row.append(power-incoming_gluons)
        slopes[str(power)] = row
    return {'additional_real_scalar_momentum_kernel':real_momentum,
            'virtual_gluon_delta_kernel':-delta_b/2,
            'hard_minus_incoming_gluon_coefficients_gg_qg_qq':slopes}


def clean(value):
    if isinstance(value,dict):
        return {k:clean(v) for k,v in value.items()}
    if isinstance(value,list):
        return [clean(v) for v in value]
    return float(f'{value:.12g}') if isinstance(value,float) else value


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=clean({'cross_energy_shape':cross_sections(),'threshold_coefficient':first_log()})
    p=HERE/'independent_results.json'
    if args.check:
        assert result==json.loads(p.read_bytes())
        print('PASS: matched 8/13-TeV widths, 32 uncertainty corners, momentum sum and two beams')
    else:
        with p.open('x',encoding='utf-8') as stream:
            stream.write(json.dumps(result,indent=2)+'\n')

