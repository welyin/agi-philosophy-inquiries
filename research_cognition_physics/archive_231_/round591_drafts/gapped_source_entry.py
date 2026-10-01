"""591 entry: reuse 579 coercivity for the original full-metric low-energy sector.

No ground-state eigenvalues or induced masses are numerically estimated.
The script checks explicit confinement certificates for a declared compact
geometry family; spectral compactness and gap are proved in the entry note.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_full_spatial_metric as original
TARGET=HERE/'gapped_source_entry_results.json'


def run():
    path=HERE.parent/'round579_drafts/residual_energy_entry_results.json'
    old=json.loads(path.read_text('utf8'));c=old['constants']
    a,beta,rho0=c['potential_prefactor'],c['beta'],c['rho0']
    q=original.original.shared_source(8);x,y,z=np.moveaxis(q['grid'],-1,0)
    psi=1.1+.05*np.cos(x)+.02*np.sin(y)
    # Six Frobenius-orthonormal symmetric directions, with three unit
    # diagonal directions. |Tr(sum sin(theta_a) T_a)| <= 3.
    radius=.2;wmin=float(q['eps']**3*np.min(psi**6)*np.exp(-1.5*radius))
    budget=1000.;rows=[]
    for probability in (.1,.01,.001):
        R=max(rho0,np.log(budget/(probability*wmin*a))/beta)
        tail=budget/(wmin*a*np.exp(beta*R))
        Fmax=original.original.M/np.cosh(R/np.sqrt(6))**2
        assert tail<=probability*(1+2e-14)
        rows.append(dict(requested_escape_bound=probability,target_geodesic_radius=float(R),
                         certified_escape_bound=float(tail),F_boundary_threshold=float(Fmax)))
    # All needed inverse-gap weights decrease on Delta>0. This samples
    # the scalar inequality only; Delta is not the original H spectrum.
    gap=.3;ground_plus_one=2.4;d=np.geomspace(gap,1e6,701);ratios=[]
    for power in (1,2,3,5):
        values=(d+ground_plus_one)/d**power
        bound=gap**(1-power)+ground_plus_one/gap**power
        assert np.max(values)<=bound*(1+1e-14)
        ratios.append(dict(power=power,max_scalar_ratio=float(np.max(values)/bound)))
    names=('research_note_579.md','round579_drafts/residual_energy_entry_results.json',
           'research_note_589.md','research_note_590.md','joint_full_spatial_metric.py')
    return dict(status='591 gapped-source entry; not a completed numbered round',checks_passed=True,
                inherited_potential_constants=dict(a=a,beta=beta,rho0=rho0),
                declared_geometry_family=dict(N=8,log_matrix_amplitude=radius,uniform_node_volume_lower=wmin),
                illustrative_energy_budget=budget,nonempty_budget_sector_not_numerically_certified=True,
                confinement_rows=rows,scalar_inverse_gap_inequality=ratios,
                actual_ground_energy_not_computed=True,actual_gap_not_computed=True,
                adiabatic_process_error_not_yet_proved=True,
                dependency_hashes={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
