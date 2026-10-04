"""593 entry: correlated quantum geometry, actual reads, and source balance.

The geometry marginal is a genuine Gram matrix of a smooth joint wave packet.
It is not a solved Einstein-constraint state. The constraint mismatch is a
conditional analytic statement; no gravity constraint is fabricated here.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_geometry_work_noise as prior
spec=importlib.util.spec_from_file_location('read592',HERE.parent/'round592_drafts/record_band_entry.py')
read=importlib.util.module_from_spec(spec);spec.loader.exec_module(read)
TARGET=HERE/'joint_constraint_entry_results.json'


def run():
    rows=[];radius=.12;correlation=.25
    for nodes in (32,48,64):
        h,s,p,kss=prior.packet_data(nodes)
        z,weight=np.polynomial.legendre.leggauss(nodes)
        theta=.3+.4*z
        pg=weight*np.exp(-2/(1-z*z));pg/=pg.sum()
        density=pg[:,None,None]*p[None,:,:]*(1+correlation*np.sin(theta)[:,None,None]*np.sin(s)[None,:,:])
        density/=density.sum()
        # Quadrature-embedded amplitudes include the positive curved measure.
        amplitude=np.sqrt(density).reshape(nodes,-1)
        L,dL,ddL=read.instrument(s)
        marginal=amplitude@amplitude.T
        after=sum((amplitude*L[...,r].reshape(1,-1))@(amplitude*L[...,r].reshape(1,-1)).T for r in range(2))
        error=float(np.linalg.norm(after-marginal,ord='fro'))
        purity=float(np.trace(marginal@marginal))
        assert error<2e-14 and purity<1-1e-10
        A=np.sum(dL*dL,axis=-1)
        def injection(th):
            volume=prior.W*np.exp(6*radius*np.sin(th))
            return prior.HBAR**2/(2*volume[:,None,None])*kss[None,:,:]*A[None,:,:]
        q=injection(theta);total=float(np.sum(density*q))
        force_jump=float(np.sum(density*(-6*radius*np.cos(theta))[:,None,None]*q))
        errors=[]
        for step in (.02,.01,.005):
            finite=float(np.sum(density*(injection(theta+step)-injection(theta-step)))/(2*step))
            errors.append(abs(finite-force_jump))
        assert total>0 and force_jump<0 and errors[-1]<errors[0]/12
        rows.append(dict(nodes=nodes,geometry_reduced_state_error=error,
                         geometry_reduced_purity=purity,total_energy_injection=total,
                         geometry_force_source_jump=force_jump,source_derivative_errors=errors))
    assert abs(rows[-1]['total_energy_injection']-rows[-2]['total_energy_injection'])<2e-12
    c=read.constants();wmin=prior.W*np.exp(-6*radius)
    joint_c1=c['commutator_energy_coefficient']*prior.W/wmin
    joint_c0=c['commutator_norm_coefficient']*(prior.W/wmin)**2
    gain=prior.HBAR**2*read.M/(32*wmin)
    assert rows[-1]['total_energy_injection']<gain
    names=('research_note_590.md','research_note_592.md','joint_geometry_work_noise.py',
           'round592_drafts/record_band_entry.py')
    return dict(status='593 joint-constraint entry; not a completed numbered round',checks_passed=True,
                packet_rows=rows,compact_geometry=dict(local_log_scale_radius=radius,w_min=wmin),
                joint_moment_constants=dict(c1=joint_c1,c0=joint_c0,single_read_gain_bound=gain),
                joint_packet_is_entangled_and_Gauss_compatible=True,
                constraint_state_not_constructed=True,Einstein_constraint_not_implemented=True,
                if_initial_total_constraint_mean_zero_then_after_equals_injection=True,
                this_last_statement_is_conditional_analytic_not_numerical=True,
                full_operator_claims_require_declared_finite_domains=True,
                dependency_hashes={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
