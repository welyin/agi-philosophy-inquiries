"""900: rational continuum barriers, reference lower bounds and C1 residual transfer.
The PDE arguments are in the note. Exact arithmetic checks their constants;
it does not convert old collocation residuals into continuum residual bounds.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,math
HERE=Path(__file__).resolve().parent;RESEARCH=HERE.parents[1]
TARGET=HERE/'continuous_reference_bounds_results.json'
def val(x):return dict(exact=str(x),value=float(x))
def run():
    raw=json.loads((RESEARCH/'archive_531_553/544/joint_singlet_common_mass_rg_results.json').read_text('utf-8'),parse_float=str)['examples'][2]['state']
    lh,ls,p=map(F,(raw['lambda_H'],raw['lambda_s'],raw['p']))
    cx,cy=F(raw['x'])/4,F(raw['y'])/4;det=lh*ls-p*p
    assert det>0 and lh>0 and ls>0
    uh=(ls*cx-p*cy)/det;us=(lh*cy-p*cx)/det
    assert F('.65')**2<uh<F('.67')**2 and F('.54')**2<us<F('.55')**2
    assert max(lh+abs(p),ls+abs(p))<1
    dH=F('.67')**2*(F('1.05')**2-1)
    dS=F('.55')**2*(F('1.06')**2-1)
    Uplus=(dH*dH+dS*dS)/(4*F('1.8')**2)
    assert Uplus<F('.001')
    old=json.loads((RESEARCH/'archive_554_584/573/joint_gravity_material_coordinates_results.json').read_text('utf-8'))['evidence']['coefficient_bounds_check']['rational_bounds']
    A=F(old['A_upper']['exact']);B=F(old['B_upper']['exact'])+F(1,400)
    yc_low=F('.68')*(F('.23')**2+F('.19')**2+F('.17')**2)
    yc_high=F('1716229971')/F('20000000000')
    Y=F(old['Y_upper']['exact'])+yc_high
    Cmin=2-F(1,400);Cmax=F('2.002')
    lo=F(18,25);hi=F(9,8)
    lower_margin=2*yc_low-Cmax*lo**8
    upper_margin=Cmin*hi**5-B*hi-A*hi**-7-2*Y*hi**-3
    assert lower_margin>0 and upper_margin>0
    a_min=5*Cmin*lo**4-B+6*yc_low*hi**-4
    a_max=5*Cmax*hi**4+7*A*lo**-8+6*Y*lo**-4
    assert a_min>F('2.8')
    hmin=F('.95')*F('.65')
    bz_min=F('.12')**2*hmin**2/8
    yz_min=F('.12')**2*((F('.41')*F('.62'))**2+(F('.29')*F('.8'))**2)/(2*F('.421'))
    force_min=bz_min*lo+2*yz_min*hi**-3
    sine_amplitude=force_min/(32+a_max)
    assert sine_amplitude*(32+a_max)==force_min
    # sin(2z) >= sqrt(1/2) > 7/10 on [pi/8,3pi/8].
    dz_lower=F(7,10)*sine_amplitude
    magnetic=F('.31')**2*F('.27')**2+(F('.31')**2*F('.21')**2+F('.27')**2*F('.21')**2)/4
    magz_lower=8*magnetic*hi**-9*dz_lower
    brmin=F('.06')*F('.54');brmax=F('.06')*F('.55')
    Ch_min=F('1.8')*(F('1.05')*F('.65')*F('.54')/12)*(F('.0056')/brmax-F('.025'))
    vh_min=Ch_min*hi**-6;ep=F(1,50)
    spatial_det=brmin*ep*magz_lower;spacetime_det=vh_min*spatial_det
    # A simple all-frequency torus gradient bound. exp(1)>271/100 follows
    # from sum_{n=0}^5 1/n!, pi>3 and 2/sqrt(3)<7/6.
    assert sum((F(1,math.factorial(n)) for n in range(6)),F(0))>F(271,100)
    r=F(100,271)
    heat_integral=F(7,6)+(2*r/(1-r)**2)*(1+2*r/(1-r))**2
    assert heat_integral<10
    c0=1/a_min;c1=F(10,8)*(1+a_max/a_min)
    laplace_bound=(Cmax*hi**5+B*hi+A*lo**-7+2*Y*lo**-3)/8
    gradient_bound=10*laplace_bound
    mxy_bound=8*magnetic*lo**-9*gradient_bound
    inverse_infinity_bound=max(1/ep,1/brmin,(1+mxy_bound/ep+mxy_bound/brmin)/magz_lower)
    # At the special line, delta grad M: use gradient of the approximant,
    # not an unproved numerical value of the exact solution.
    Mgrad_res_constant_without_candidate_gradient=8*magnetic*lo**-9*c1
    Mgrad_res_gradient_multiplier=72*magnetic*lo**-10*c0
    working=json.loads((HERE/'current_reference_initial_budget_results.json').read_text('utf-8'))
    diagnostics=[]
    for row in working['rows']:
        assert float(lo)<row['psi']<float(hi)
        assert row['actual_spatial_material_J'][2][2]>float(magz_lower)
        assert row['spacetime_determinant']>float(spacetime_det)
        diagnostics.append(dict(N=row['N'],observed_magz_to_analytic_lower=row['actual_spatial_material_J'][2][2]/float(magz_lower),
            observed_Jacobian_to_analytic_lower=row['spacetime_determinant']/float(spacetime_det)))
    return dict(round=900,date='2026-10-06',formal_reports=900,cumulative_numbered_groups=3685,fresh_numbered_groups=1,
        all_checks_passed=True,argument_scope='Original859 continuous Cauchy branch: rational barriers, quantitative z-gradient and material rank on a compact initial line, and a conditional C0/C1 certificate from a true continuum residual. No spacetime patch or quantum-source certificate yet.',
        original_frozen_decimal_potential_parameters_used=True,potential_Uplus_upper=val(Uplus),
        coefficients={k:val(v) for k,v in dict(A_max=A,B_max=B,Y_max=Y,Y_min=yc_low,C_min=Cmin,C_max=Cmax).items()},
        psi_lower=val(lo),psi_upper=val(hi),lower_barrier_margin=val(lower_margin),upper_barrier_margin=val(upper_margin),
        reaction_derivative_lower=val(a_min),reaction_derivative_upper=val(a_max),
        b_z_lower=val(bz_min),y_z_lower=val(yz_min),forcing_lower=val(force_min),sine_subsolution_amplitude=val(sine_amplitude),
        interval='z in [pi/8,3pi/8], x/y periodic for psi_z; material rank on x=0,y=pi/2 only',
        negative_psi_z_lower=val(dz_lower),magnetic_z_lower=val(magz_lower),clock_velocity_lower=val(vh_min),
        fixed_probe_amplitude=str(ep),spatial_Jacobian_lower=val(spatial_det),spacetime_Jacobian_lower=val(spacetime_det),
        heat_gradient_integral_upper=val(heat_integral),C0_error_per_continuum_residual=val(c0),each_C1_error_per_continuum_residual=val(c1),
        exact_solution_each_gradient_bound=val(gradient_bound),spatial_inverse_infinity_upper=val(inverse_infinity_bound),
        magnetic_gradient_error_residual_constant=val(Mgrad_res_constant_without_candidate_gradient),
        magnetic_gradient_error_residual_times_candidate_gradient=val(Mgrad_res_gradient_multiplier),
        existing_initial_data_diagnostics=diagnostics,
        collocation_residual_reinterpreted_as_continuum=False,full_continuum_residual_of_numeric_psi_bounded=False,
        full_spacetime_reference_patch_certified=False,original_quantum_modes_or_source_M_T_computed=False,
        no_new_fields_or_physical_thresholds=True,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in r.items() if k not in ('coefficients','existing_initial_data_diagnostics')},ensure_ascii=False,indent=2))
