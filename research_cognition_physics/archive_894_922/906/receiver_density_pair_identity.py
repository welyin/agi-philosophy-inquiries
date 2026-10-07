"""906: on the872 lambda^2 four-leg coefficient, cancel the density pair together.
Tests use the actual neutral doublet, curved jets, and forced Dirac equation.
They do not evaluate the original physical source j or the full mean response.
"""
from pathlib import Path
import argparse,json
import numpy as np
import compact_receiver_propagation as r
from compact_receiver_checks import point_action,point_geometry
TARGET=Path(__file__).with_name('receiver_density_pair_identity_results.json')

def run():
    rng=np.random.default_rng(909);maximum=0.;forced_error=0.;omission=0.
    sy=r.SIG[1];c=.75
    for _ in range(20):
        a=.03*rng.normal(size=(4,4));g=np.diag([-1.,1.,1.,1.])+a+a.T
        d=.04*rng.normal(size=(4,4,4));d=(d+d.transpose(0,2,1))/2;y=dict(g=g);z=point_geometry(g,d)
        z.update(vol=np.sqrt(np.linalg.det(g[1:,1:])),invgamma=np.linalg.inv(g[1:,1:]));f=r.shared.frame_data(y,z)
        u=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2));v=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
        dux=rng.normal(size=(3,4,2))+1j*rng.normal(size=(3,4,2));dvx=rng.normal(size=(3,4,2))+1j*rng.normal(size=(3,4,2))
        ell=float(rng.normal());Su=u@sy.T
        um=[];vm=[];mixed_L=0.;mixed_full=np.zeros((4,4));mixed_short=np.zeros((4,4))
        for a in range(2):
            du=np.r_[(-1j*point_action(y,z,f,u[:,a],dux[:,:,a]))[None,:],dux[:,:,a]]
            dv=np.r_[(-1j*point_action(y,z,f,v[:,a],dvx[:,:,a])-1j*z['alpha']*c*ell*r.BETA@Su[:,a])[None,:],dvx[:,:,a]]
            orig=r.bilinear_jet(y,z,f,u[:,a],du);forced=r.bilinear_jet(y,z,f,v[:,a],dv)
            forced_error=max(forced_error,float(np.max(abs(forced['dirac_residual']-c*ell*r.BETA@Su[:,a]/np.sqrt(z['vol'])))))
            plus=r.bilinear_jet(y,z,f,u[:,a]+v[:,a],du+dv);minus=r.bilinear_jet(y,z,f,u[:,a]-v[:,a],du-dv)
            mixed_L+=(plus['lagrangian']-minus['lagrangian'])/4
            mixed_full+=(plus['stress']-minus['stress'])/4
            mixed_short+=(plus['onshell_stress']-minus['onshell_stress'])/4
        O=float(np.einsum('af,ab,bf->',u.conj(),r.BETA,Su).real/z['vol'])
        metric_density_term=-c*ell*O*g
        error=max(abs(2*mixed_L-c*ell*O),np.max(abs(2*(mixed_full-mixed_short)+metric_density_term)))
        maximum=max(maximum,float(error));omission=max(omission,float(np.max(abs(metric_density_term))))
    assert maximum<5e-12 and forced_error<5e-13 and omission>.01,(maximum,forced_error,omission)
    return dict(round=906,date='2026-10-06',all_checks_passed=True,curved_forced_doublet_jets=20,
        forced_Dirac_equation_error=forced_error,paired_density_cancellation_error=maximum,
        largest_uncancelled_piece_if_only_one_is_omitted=omission,
        statement='On the872 lambda^2 four-receiver-leg coefficient, 2 L_R(psi1,psi32)=c ell(b21) O(psi1). The corresponding minimal-Dirac Hilbert density source cancels exactly the -c nu_star ell(b21) O term of F1_prime, before and after the same Pi_adjoint.',
        remaining_b_nu_term_required=True,other_old_F0_vertices_retained=True,full_source_not_zero_by_this_identity=True,
        physical_source_or_mean_numerically_evaluated=False,extra_model_principle=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();v=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert v==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(v,ensure_ascii=False,indent=2))
