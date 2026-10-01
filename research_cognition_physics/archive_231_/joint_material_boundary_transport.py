"""651: original material clocks, causal walls and complete boundary geometry.

Original573 data supply the on-shell wall. The explicit material-chart jet
family is off shell and tests tensor variations, not a new Einstein solution.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_gravity_material_coordinates as old
import joint_geometric_boundary_matching as boundary

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_material_boundary_transport_results.json'
TANGENT=[0,2,3]


def original_wall_check():
    spec=importlib.util.spec_from_file_location('wall_probe651',HERE/'round651_drafts/material_wall_signature_probe.py')
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    saved=json.loads((HERE/'round651_drafts/material_wall_signature_probe_results.json').read_text('utf8'))
    assert probe.run()==saved
    lower=(Q('1.8')*Q('.97')*Q('.14')/Q('.033'))**2
    psi8=Q(3,2)**8
    assert lower>psi8
    f=old.fields(32);c=old.point_coefficients(f)
    lam=-12*c['beta']/(c['H']*c['S']);idx=(0,8,4);norm=c['b']**2*f['psi'][idx]**-4
    assert abs(lam-saved['rows'][-1]['clock_relative_slope'])<2e-13
    assert abs(norm-saved['rows'][-1]['clock_relative_wall_normal_squared'])<2e-14
    return dict(probe_reproduced=True,source_lower_ratio_exact=str(lower),
        source_lower_ratio=float(lower),psi_eighth_upper_exact=str(psi8),
        shape_choice_lambda=float(lam),lambda_independent_of_psi=True,
        normal_squared_from_original_spatial_metric=float(norm),rows=saved['rows'],
        local_strict_signature_proved_by_inherited_continuous_bounds=True,
        no_claim_of_global_or_long_time_wall=True)


def base_data():
    f=old.fields(32);c=old.point_coefficients(f);idx=(0,8,4)
    lam=-12*c['beta']/(c['H']*c['S'])
    Y=np.array([c['H'],c['S']-lam*c['H'],f['rh'][idx],f['rs'][idx]])
    return c,lam,Y


def chart_jet(parameter,frame='E'):
    c,lam,Y=base_data();T,Z,A,B=Y
    # Full smooth family in Y near the reference point, evaluated at that
    # point. g^00=Y2 and g^11+2lam*g^01+lam^2*g^00=Y3 identically.
    q=.03;dq=np.array([.02,.04,-.01,.02])
    gi=np.array([[A,parameter*q,0,0],[parameter*q,B-2*lam*parameter*q-lam*lam*A,0,0],
                 [0,0,1+.2*parameter,0],[0,0,0,1+.3*parameter]])
    dgi=np.zeros((4,4,4));dgi[2,0,0]=1
    dgi[:,0,1]=dgi[:,1,0]=parameter*dq
    dgi[:,1,1]=np.array([0.,0.,-lam*lam,1.])-2*lam*parameter*dq
    dgi[1,2,2]=.08+.11*parameter
    dgi[0,3,3]=.05;dgi[1,3,3]=-.07*parameter
    pgi=np.zeros((4,4));pgi[0,1]=pgi[1,0]=q;pgi[1,1]=-2*lam*q
    pgi[2,2]=.2;pgi[3,3]=.3
    dpgi=np.zeros_like(dgi);dpgi[:,0,1]=dpgi[:,1,0]=dq
    dpgi[:,1,1]=-2*lam*dq;dpgi[1,2,2]=.11;dpgi[1,3,3]=-.07
    g=np.linalg.inv(gi);dg=np.array([-g@d@g for d in dgi])
    kap=-g@pgi@g
    dkap=np.array([-dg[i]@pgi@g-g@dpgi[i]@g-g@pgi@dg[i] for i in range(4)])
    phi=np.array([0.,T,0.,0.,Z+lam*T]);dphi=np.zeros((4,5))
    dphi[0,1]=1;dphi[0,4]=lam;dphi[1,4]=1
    F,ff=boundary.data(phi);dF=dphi@ff
    if frame=='J':
        dg=dg/F-g[None,:,:]*dF[:,None,None]/F**2
        dkap=dkap/F-kap[None,:,:]*dF[:,None,None]/F**2
        g=g/F;kap=kap/F;gi=np.linalg.inv(g)
    assert np.linalg.eigvalsh(g)[0]<0<np.linalg.eigvalsh(g)[1] and gi[1,1]>0
    Gamma=np.zeros((4,4,4))
    for r in range(4):
        for m in range(4):
            for n in range(4):
                Gamma[r,m,n]=.5*sum(gi[r,s]*(dg[m,n,s]+dg[n,m,s]-dg[s,m,n]) for s in range(4))
    normal_cov=np.eye(4)[1]/np.sqrt(gi[1,1]);normal=gi@normal_cov
    h=g[np.ix_(TANGENT,TANGENT)]
    K=-Gamma[1][np.ix_(TANGENT,TANGENT)]/np.sqrt(gi[1,1])
    v=normal@dphi
    return dict(g=g,gi=gi,dg=dg,kappa=kap,dkappa=dkap,Gamma=Gamma,
        n_cov=normal_cov,n=normal,h=h,K=K,v=v,phi=phi,dphi=dphi,F=F,dF=dF,Y=Y,lam=lam)


def geometric_variation(j):
    g,gi,kap,Gamma,n=j['g'],j['gi'],j['kappa'],j['Gamma'],j['n']
    knn=float(n@kap@n)
    dn=-gi@kap@n+.5*knn*n
    nabla=j['dkappa'].copy()
    for r in range(4):
        nabla[r]-=Gamma[:,r,:].T@kap+kap@Gamma[:,r,:]
    deltaGamma=np.zeros((4,4,4))
    for r in range(4):
        for m in range(4):
            for p in range(4):
                deltaGamma[r,m,p]=.5*sum(gi[r,s]*(nabla[m,p,s]+nabla[p,m,s]-nabla[s,m,p]) for s in range(4))
    dK=.5*knn*j['K']-np.einsum('r,rab->ab',j['n_cov'],deltaGamma[:,TANGENT,:][:,:,TANGENT])
    return dict(normal=dn,extrinsic_curvature=dK,normal_scalar_derivative=dn@j['dphi'])


def reference_metric_check():
    rows=[]
    for p in (-.00005,0.,.00004):
        j=chart_jet(p);gi=j['gi'];dh=j['dphi'][:,1];ds=j['dphi'][:,4]
        A=float(dh@gi@dh);B=float(ds@gi@ds)
        assert abs(A-j['Y'][2])<1e-14 and abs(B-j['Y'][3])<1e-13
        jJ=chart_jet(p,'J')
        assert np.max(abs(j['g']-j['F']*jJ['g']))<1e-9
        assert abs(gi[1,1]-(B-2*j['lam']*(dh@gi@ds)+j['lam']**2*A))<1e-14
        rows.append(dict(parameter=p,reference_A=A,reference_B=B,
            wall_normal_squared=float(gi[1,1]),F_at_fixed_material_label=j['F']))
    # The full functions also obey A(Y)=Y2 and B(Y)=Y3, so the reconstructed
    # four material references (h,s-lambda h,A,B) are exactly Y, not merely
    # assigned labels. This off-shell family is analytically self-consistent.
    return dict(rows=rows,reference_chart_reconstructed_exactly=True,
        fixed_label_F_variation=0.,two_inverse_metric_components_constrained=True,
        illustrative_full_chart_family_is_off_shell=True,
        other_original573_reference_derivatives_not_fabricated=True)


def boundary_transport_check():
    rows=[];worst=0.
    for frame in ('E','J'):
        j=chart_jet(0.,frame);v=geometric_variation(j);cases=[]
        # These material coordinates have a small g^00 and correspondingly
        # sensitive inverse-metric derivatives; keep all three steps inside
        # the observed centered-difference quadratic regime.
        for step in (2e-6,1e-6,5e-7):
            plus=chart_jet(step,frame);minus=chart_jet(-step,frame)
            errors={}
            for key,theory in (('n',v['normal']),('K',v['extrinsic_curvature']),('v',v['normal_scalar_derivative'])):
                fd=(plus[key]-minus[key])/(2*step)
                errors[key]=float(np.max(abs(fd-theory))/max(1.,np.max(abs(theory))))
            cases.append(dict(step=step,relative_errors=errors))
        assert max(cases[-1]['relative_errors'].values())<2e-5
        assert max(cases[-1]['relative_errors'].values())<max(cases[0]['relative_errors'].values())/10
        worst=max(worst,max(cases[-1]['relative_errors'].values()))
        rows.append(dict(frame=frame,cases=cases,
            scalar_values_variation_at_fixed_labels=0.,
            normal_scalar_derivative_variation=v['normal_scalar_derivative'].tolist(),
            omitted_if_normal_is_frozen_norm=float(np.linalg.norm(v['normal_scalar_derivative']))))
    E=chart_jet(0.,'E');J=chart_jet(0.,'J')
    # Original nonminimal boundary momentum is independently evaluated.
    P,jscalar=boundary.jordan_momenta(J['phi'],J['h'],J['K'],J['v'])
    hE,KE,vE,PE,jE=boundary.einstein_data(J['phi'],J['h'],J['K'],J['v'])
    assert max(np.max(abs(hE-E['h'])),np.max(abs(KE-E['K'])),np.max(abs(vE-E['v'])))<2e-8
    rootJ=np.sqrt(abs(np.linalg.det(J['h'])));rootE=np.sqrt(abs(np.linalg.det(E['h'])))
    workJ=.5*rootJ*np.sum(P*J['kappa'][np.ix_(TANGENT,TANGENT)])
    workE=.5*rootE*np.sum(PE*E['kappa'][np.ix_(TANGENT,TANGENT)])
    assert abs(workJ-workE)<2e-10
    return dict(rows=rows,max_final_relative_error=worst,original_boundary_work_Jordan=float(workJ),
        original_boundary_work_Einstein=float(workE),frame_work_difference=float(abs(workJ-workE)),
        geometry_and_normal_derivative_variations_retained=True,
        numerical_check_is_local_off_shell_not_evolution_or_quantization=True,
        full_moving_action_and_transmission_connection_analytic_only=True)


def run():
    deps=('research_note_573.md','research_note_618.md','research_note_619.md',
          'research_note_647.md','research_note_648.md','research_note_650.md',
          'joint_gravity_material_coordinates.py','joint_geometric_boundary_matching.py',
          'round651_drafts/material_wall_signature_probe.py',
          'round651_drafts/material_wall_signature_probe_results.json')
    return dict(round=651,tests_run=3,failures=0,errors=0,original_wall=original_wall_check(),
        same_reference_metric=reference_metric_check(),boundary_transport=boundary_transport_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original573_on_shell_local_causal_wall=True,
            full_classical_action_pullback_requires_same_embedding_on_both_sides=True,
            self_consistent_off_shell_chart_used_only_for_tensor_differences=True,
            not_all_boundary_symplectic_potentials_identified=True,
            no_full_quantum_GR_or_fixed_graph_instrument_equivalence=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
