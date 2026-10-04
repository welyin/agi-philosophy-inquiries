"""736 entry: reuse630/631 spectra to test the response's logarithmic order.

No new masses, loop coefficients, physical cutoff or actual reference reset.
This is the previously declared constant-background comparison branch.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_continuum_source_spectrum as scalar
import joint_tensor_stress_spectrum as tensor
TARGET=HERE/'response_regularity_entry_results.json'


def slope(Q,masses,weights,channel,n=256):
    # Differentiate the subtracted Euclidean ratio under its absolutely
    # convergent integral, then use u=log(omega/(2m)) for each original mass.
    nodes,quad=scalar.quadrature(n)
    total=0.
    for m,d in zip(masses,weights):
        extent=np.log(Q/(2*m))+18
        omega=2*m*np.exp(extent*nodes)
        if channel=='scalar':
            rho=scalar.density(omega,np.array([m]),np.array([d]))
            integrand=2*Q*Q/np.pi*rho/(omega*omega+Q*Q)**2
        else:
            rho=tensor.spin2(omega,np.array([m]),np.array([d]))
            integrand=2*Q*Q/np.pi*rho/(omega*omega*(omega*omega+Q*Q)**2)
        total+=extent*np.dot(quad,integrand)
    return float(total)


def ratio(Q,masses,weights,channel):
    # These are the old analytic dispersive remainders, evaluated on the
    # imaginary axis (no threshold pole); not newly fitted response kernels.
    if channel=='scalar':
        return float(scalar.finite_subtracted(1j*Q,masses,weights,384).real/Q**2)
    return float(-tensor.tensor_dispersion(1j*Q,masses,weights,384).real/Q**4)


def run():
    _,masses,weights,_,_=scalar.original_data()
    expected={'scalar':float(np.dot(weights,masses**2)/(4*np.pi**2)),
              'tensor':float(weights.sum()/(80*np.pi**2))}
    rows=[]
    for channel in ('scalar','tensor'):
        Q=12*max(masses)
        h=2e-4
        fd=(ratio(Q*np.exp(h),masses,weights,channel)-ratio(Q*np.exp(-h),masses,weights,channel))/(2*h)
        direct=slope(Q,masses,weights,channel)
        assert abs(fd-direct)<2e-10
        points=[]
        for factor in (32,128):
            q=factor*max(masses)
            a=slope(q,masses,weights,channel,192)
            b=slope(q,masses,weights,channel,384)
            assert abs(a-b)<2e-10
            points.append(dict(Q=float(q),slope=b,ratio_to_limit=b/expected[channel],
                              quadrature_change=abs(a-b)))
        assert points[1]['ratio_to_limit']>points[0]['ratio_to_limit']>.97
        assert points[1]['ratio_to_limit']<1.000001
        rows.append(dict(channel=channel,logarithmic_coefficient=expected[channel],
                         inherited_coefficient_not_a_new_loop_calculation=True,
                         analytic_integrand_slope=direct,independent_remainder_difference=fd,
                         derivative_error=abs(fd-direct),points=points))
    deps=('research_note_630.md','research_note_631.md','research_note_601.md',
          'research_note_734.md','research_note_735.md',
          'joint_continuum_source_spectrum.py','joint_tensor_stress_spectrum.py')
    return dict(entry_round=736,new_formal_round=False,latest_completed_round=735,
                formal_tests_unchanged=3409,rows=rows,
                dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in deps},
                scope='Previously declared flat constant-mass whole-generation vacuum branch only. Its nonlocal responses have second/fourth order times a logarithm; finite local polynomials cannot remove that logarithm. This rejects a proposed no-log Sobolev estimate, not the unified model or existence of semiclassical solutions.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as handle:
            handle.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
