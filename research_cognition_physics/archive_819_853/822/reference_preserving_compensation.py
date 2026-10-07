"""822 working result: gauge-only momentum repair preserves scalar clock jets.
Uses the original 753 slice and actual EW connection. Sources in the numerical
calibration remain the explicitly prescribed 754 diagnostic.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'820'))
import color_all_mode_certificate as fixed
previous=fixed.previous;geo=previous.geo
TARGET=HERE/'reference_preserving_compensation_results.json'

def exact_diagonal():
    A=Q('.41');B=Q('.29');C=Q('.12');u=Q('.38');v=Q('.2')
    mean4=lambda u:1+3*u*u+3*u**4/8
    derivative2=lambda u:4*u*u+u**4
    coeff=[2*A**4*B**4*derivative2(u)*mean4(v),
           2*A**4*B**4*mean4(u)*derivative2(v),
           2*A**4*C**4*mean4(u)]
    assert min(coeff)>0
    return coeff

def run():
    d=previous.setup();q=d['q'];psi=d['psi']
    A,Ec,Fc,Cc=fixed.probe.data()
    dp,dE,dE0,info=previous.old.inverse_gauss(q,d['k'],d['sigma'])
    fw,_=previous.old.curvature(q)
    xy=np.sum(fw[...,0,1,:]**2,axis=-1);xz=np.sum(fw[...,0,2,:]**2,axis=-1)
    menus=[];gauss=[]
    zero_p=np.zeros_like(dp);zero_e0=np.zeros_like(dE0)
    for pair,direction,density in (((0,1),0,xy),((0,1),1,xy),((0,2),2,xz)):
        i,j=pair;weight=geo.derivative(density,direction);weighted=weight[...,None]*fw[...,i,j,:]
        em=np.zeros_like(dE);em[...,i,:]=previous.old.covariant(weighted,j,q)
        em[...,j,:]=-previous.old.covariant(weighted,i,q)
        menus.append(em)
        gauss.append(float(np.max(abs(previous.old.gauss(q,d['k'],zero_p,em,zero_e0)))))
    matrix=np.column_stack([q['dx']**3*previous.old.momentum(q,zero_p,m,zero_e0).sum(axis=(0,1,2)) for m in menus])
    fractions=exact_diagonal();expected=np.diag([float(f)*np.pi**3 for f in fractions])
    matrix_error=float(np.max(abs(matrix-expected)))
    assert matrix_error<1e-13 and max(gauss)<1e-12 and min(np.diag(matrix))>0
    g=-d['sigma_c']
    target_integral=-q['dx']**3*(d['J']+previous.old.momentum(q,dp,dE,dE0)).sum(axis=(0,1,2))
    target_integral-=A@(q['dx']**3*g.sum(axis=(0,1,2)))
    coeff=np.linalg.solve(matrix,target_integral)
    for value,menu in zip(coeff,menus):dE+=value*menu
    bc,bw,b0=geo.old.PAR['b']
    pe=(psi**-12*np.sum(previous.old.kinverse(q,q['p'])*dp,axis=-1)
        +2*psi**-8*(bw*np.sum(q['f']['E']*dE,axis=(-1,-2))
                     +b0*np.sum(q['f']['E0']*dE0,axis=-1)))
    mom=previous.old.momentum(q,dp,dE,dE0)
    dec=fixed.solve_color(g,-d['J']-mom,-(d['rho']+pe)*psi**8/(2*bc))
    residuals=dict(EW_Gauss=float(np.max(abs(previous.old.gauss(q,d['k'],dp,dE,dE0)+d['sigma']))),
        color_Gauss=float(np.max(abs(sum(previous.D(dec[...,i,:],i,Cc) for i in range(3))+d['sigma_c']))),
        momentum=float(np.max(abs(mom+np.einsum('...ja,ija->...i',dec,Fc)+d['J']))),
        energy=float(np.max(abs(pe+2*bc*psi**-8*np.sum(Ec*dec,axis=(-1,-2))+d['rho']))))
    assert max(residuals.values())<2e-10
    velocity=psi[...,None]**-6*previous.old.kinverse(q,q['p'])
    change=psi[...,None]**-6*previous.old.kinverse(q,dp)
    clocks=np.stack((-2*velocity[...,1]*change[...,1],-2*velocity[...,4]*change[...,4]),axis=-1)
    old_change=psi[...,None]**-6*previous.old.kinverse(q,d['dp'])
    old_clocks=np.stack((-2*velocity[...,1]*old_change[...,1],-2*velocity[...,4]*old_change[...,4]),axis=-1)
    tangency=float(np.max(abs(np.sum(q['phi']*dp,axis=-1))))
    assert tangency<1e-13 and np.max(abs(clocks))<1e-13
    assert np.max(abs(old_clocks))>1e-4
    return dict(round=822,formal_round_completed=False,all_working_checks_passed=True,
        gauge_only_momentum_matrix=matrix.tolist(),matrix_formula_error=matrix_error,
        positive_diagonal_divided_by_pi_cubed=[str(v) for v in fractions],
        menu_Gauss_errors=gauss,correction_coefficients=coeff.tolist(),
        full_linear_constraints=residuals,
        scalar_momentum_tangency_error=tangency,
        old_compensation_scalar_reference_change_max=float(np.max(abs(old_clocks))),
        new_compensation_scalar_reference_change_max=float(np.max(abs(clocks))),
        scalar_h_s_values_and_first_jets_preserved=True,
        four_reference_values_and_spatial_jets_preserved=True,
        original_special_point_reference_determinant_preserved_analytically=True,
        all_reference_normal_derivatives_or_full_relational_metric_preserved=False,
        source_is_declared754_diagnostic_not_actual819=True,
        actual_quantum_absolute_means_matched=False,
        autonomous_implementation_proven=False,formal_test_groups_added=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
