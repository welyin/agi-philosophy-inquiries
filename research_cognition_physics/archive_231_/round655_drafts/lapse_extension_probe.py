"""655 entry: same flat transfer can have inequivalent lapse extensions.

Uses the actual654 mass and653 auxiliary spectrum. The two extensions are
declared choices, not two predictions of one already specified geometry theory.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_mass_state_time_limit as original


def run():
    data=json.loads((ROOT/'joint_chiral_character_state_results.json').read_text('utf8'))
    spec=data['temporal_gluing']['spherical_harmonic_spectrum']
    logs=np.log([float(Fraction(r['eigenvalue'])) for r in spec])
    dims=np.array([r['dimension'] for r in spec]);relative=logs-logs[0]
    sigmas=np.linalg.svd(original.delta(original.PHIS[0]),compute_uv=False)
    beta=2.;rows=[];k=32*np.log(2)-logs[0]
    for n in (2,4,8,16,32,64):
        a=beta/n;energy=2/a*np.arcsinh(a*sigmas/2)
        def logz(eta,family,normalized=False):
            if family=='rebuild_mass_at_fixed_grid':
                x=beta/a*np.arcsinh(a*eta*sigmas/2)
                aux=n*logs[0]+np.log(np.sum(dims*np.exp(n*relative)))
                value=-32*n*np.log(2)+np.sum(np.logaddexp(x,-x))+aux
                return float(value+n*k if normalized else value)
            x=beta*eta*energy/2
            aux=n*eta*logs[0]+np.log(np.sum(dims*np.exp(n*eta*relative)))
            value=-32*n*eta*np.log(2)+np.sum(np.logaddexp(x,-x))+aux
            return float(value+n*eta*k if normalized else value)
        assert abs(logz(1.,'rebuild_mass_at_fixed_grid')-logz(1.,'lapse_of_same_generator'))<1e-11
        p=np.sum(beta*sigmas/(2*np.sqrt(1+(a*sigmas/2)**2))*np.tanh(beta*energy/2))
        weights=dims*np.exp(n*relative);weights/=sum(weights)
        hn=np.sum(beta*energy/2*np.tanh(beta*energy/2))+n*(weights@relative)
        h=hn-n*k;step=1e-5
        fd=[]
        for family,expected,normalized in [('rebuild_mass_at_fixed_grid',p,False),
            ('lapse_of_same_generator',h,False),('lapse_of_same_generator',hn,True)]:
            actual=(logz(1+step,family,normalized)-logz(1-step,family,normalized))/(2*step)
            assert abs(actual-expected)<2e-7
            fd.append(float(abs(actual-expected)))
        rows.append(dict(time_sites=n,time_step=a,same_base_log_partition=logz(1.,'rebuild_mass_at_fixed_grid'),
            fixed_grid_raw_lapse_source=float(p),same_generator_raw_lapse_source=float(h),
            same_generator_vacuum_subtracted_source=float(hn),vacuum_source_difference=float(n*k),
            normalized_source_difference=float(hn-p),finite_difference_errors=fd))
    return dict(date='2026-10-02',status='655 entry; not a completed round',rows=rows,
        dependency_hashes={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in
        ('joint_mass_state_time_limit.py','joint_chiral_character_state_results.json','research_note_643.md')},
        conclusion='Agreement of the flat transfer at lapse=1 does not identify its geometric source; an off-background prescription is required. No cosmological constant inference.')


if __name__=='__main__':
    result=run()
    with (HERE/'lapse_extension_probe_results.json').open('x',encoding='utf8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps(dict(status=result['status'],cases=len(result['rows']))))
