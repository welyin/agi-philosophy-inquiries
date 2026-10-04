"""709 entry: local phase-conjugation calibration in the inherited H5 model.
No Gibbs spectrum, quantum recovery map or independent research round.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_geodesic_spatial_block as old
TARGET=HERE/'conditional_reference_entry_results.json'


def run():
    points=np.array([[.12,-.14,.08,.03,.35],[-.07,.18,-.05,.11,-.16],[.08,.02,.12,-.09,.21]])
    XX=np.array([old.hyper(p) for p in points]);z=XX[:,5];Z=z.sum();n=len(points)
    gamma=n+z@z/old.R**2
    dg=np.cos(Z)**2*gamma
    lapg=-np.sin(Z)*gamma+np.cos(Z)*5/old.R**2*Z
    alpha=.43;hbar=.7;w=.8;g=np.sin(Z);step=2e-4
    # Hsc e^(i alpha g) / e^(i alpha g) at a constant test function.
    # Uses every orthonormal geodesic direction in the original product target.
    lap=0j
    for p,X,zi in zip(points,XX,z):
        E=old.frame(p)
        for a in range(5):
            plus=np.cosh(step/old.R)*X+old.R*np.sinh(step/old.R)*E[:,a]
            minus=np.cosh(step/old.R)*X-old.R*np.sinh(step/old.R)*E[:,a]
            lap+=(np.exp(1j*alpha*np.sin(Z-zi+plus[5]))
                  +np.exp(1j*alpha*np.sin(Z-zi+minus[5]))-2*np.exp(1j*alpha*g))/step**2
    actual=-hbar*hbar/(2*w)*lap*np.exp(-1j*alpha*g)
    expected=hbar*hbar/(2*w)*(alpha*alpha*dg-1j*alpha*lapg)
    error=float(abs(actual-expected));assert error<2e-7
    deps=('research_note_577.md','research_note_598.md','research_note_641.md',
          'research_note_682.md','research_note_708.md','joint_block_fluctuation_contract_results.json')
    return dict(date='2026-10-03',entry_round=709,latest_formal_round=708,new_formal_round=False,
        product_target_dimensions=15,local_gradient_coefficient=float(dg),
        phase_kinetic_actual=[float(actual.real),float(actual.imag)],
        phase_kinetic_expected=[float(expected.real),float(expected.imag)],
        absolute_error=error,full_Gibbs_coefficient_numerically_computed=False,
        conditional_expectation_identities_analytic=True,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r))
