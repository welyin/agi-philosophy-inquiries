"""873: mass-source mixing induced by the original material-clock constraint.
Exact Grassmann implicit solving plus a finite-CAR calibration at stored861
background coefficients. Not a quantization of the full reduced field theory.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
TARGET=HERE/'clock_mass_source_bridge_results.json'
# Keys are (first mass power, second mass power, exterior-algebra mask).
def add(*xs):
    z={}
    for x in xs:
        for k,v in x.items():z[k]=z.get(k,F(0))+v
    return {k:v for k,v in z.items() if v}
def scale(x,a):return {k:v*a for k,v in x.items() if v*a}
def mul(x,y):
    z={}
    for (i,j,m),a in x.items():
        for (k,l,n),b in y.items():
            if m&n:continue
            sign=(-1)**sum((n&((1<<r)-1)).bit_count() for r in range(4) if (m>>r)&1)
            key=(i+k,j+l,m|n);z[key]=z.get(key,F(0))+sign*a*b
    return {k:v for k,v in z.items() if v}
def const(a):return {(0,0,0):F(a)}
def exact_implicit():
    a=F(2);b=F(1,3);p0=F(2,5);v=2*a*p0+b
    u1=F(3,7);u2=F(5,11)
    mass={(1,0,3):u1,(0,1,12):u2}
    # Unrelated, mass-independent lower fermion terms are retained. They
    # model the general analytic remainder after the same spatial reduction.
    ar={(0,0,5):F(1,13)}
    br={(0,0,3):F(2,9),(0,0,10):F(-1,8)}
    cr={(0,0,12):F(4,9),(0,0,6):F(1,6),(0,0,15):F(5,17)}
    def other(dp,keep_curvature=True):
        pi=add(const(p0),dp)
        terms=[mul(ar,mul(pi,pi)),mul(br,pi),cr,mass]
        if keep_curvature:terms.append(scale(mul(dp,dp),a))
        return add(*terms)
    def solve(keep):
        dp={}
        for _ in range(6):
            nxt=scale(other(dp,keep),-1/v)
            if nxt==dp:break
            dp=nxt
        else:raise AssertionError('nilpotent iteration did not terminate')
        assert not add(scale(dp,v),other(dp,keep))
        return dp
    dp=solve(True);naive=solve(False);key=(1,1,15)
    got=-dp.get(key,F(0));wanted=2*a*u1*u2/v**3
    assert got==wanted and got!=0
    assert naive.get(key,F(0))==0
    bad_res=add(scale(naive,v),other(naive,True))
    assert bad_res[key]!=0
    return dict(alpha=str(a),beta=str(b),background_clock_momentum=str(p0),
        clock_speed=str(v),mass_density_amplitudes=[str(u1),str(u2)],
        mixed_mass_four_leg_coefficient=str(got),analytic_coefficient=str(wanted),
        exact_full_constraint_residual_zero=True,
        mass_independent_fermion_terms_retained=True,
        frozen_clock_mixed_coefficient='0',
        frozen_clock_original_constraint_residual=str(bad_res[key]))
def original_coefficients():
    saved=json.loads((HERE.parent/'861/magnetic_reduced_hamiltonian_results.json').read_text('utf-8'))
    rows=[]
    for bg in saved['original_859_background_clock_reduction']:
        a=bg['a'];v=bg['actual_clock_speed']
        kappa=2*a/v**3;u_star=v*v/(4*a);h_star=v/(2*a)
        x1,x2=.08,.12
        s1,s2=x1*u_star,x2*u_star
        def shift(s):return 2*s/(v+np.sqrt(v*v-4*a*s))
        values=np.array([0.,shift(s2),shift(s1),shift(s1+s2)])
        defect=values[3]-values[2]-values[1]
        assert defect>0
        # This calibrated finite CAR root is spectral functional calculus.
        # No original continuum mass-mode normalization is inferred from it.
        time=1/h_star;phases=np.exp(-1j*time*values)
        gamma_pred=phases[1]*phases[2]
        lift_defect=abs(phases[3]-gamma_pred)
        determinant=(1+phases[1])*(1+phases[2])
        fock_trace=sum(phases)
        assert abs(abs(fock_trace-determinant)-lift_defect)<3e-15
        assert lift_defect>1e-3
        # Central mixed derivative independently checks the clock coefficient.
        checks=[]
        for e in (.004,.002,.001):
            ss=e*u_star
            mixed=(shift(2*ss)-2*shift(0.)+shift(-2*ss))/(4*ss*ss)
            rel=abs(mixed/kappa-1)
            checks.append(dict(source_step_fraction=e,mixed_derivative=float(mixed),relative_error=float(rel)))
        assert checks[-1]['relative_error']<1e-5
        assert checks[-1]['relative_error']<checks[0]['relative_error']/10
        assert abs(kappa*u_star**2-h_star/4)<1e-16
        # The ZZ coefficient cannot belong to any diagonal quadratic CAR generator.
        zz=defect/4
        residuals=[]
        # The root shift lowers pi_h by e; use the background root identity.
        for source,e in zip((0.,s2,s1,s1+s2),values):
            residuals.append(float(abs(a*e*e-v*e+source)))
        assert max(residuals)<1e-19
        rows.append(dict(grid_N=bg['grid_N'],probe_amplitude=bg['probe_amplitude'],
            alpha=a,clock_speed=v,mixed_mass_density_coefficient=kappa,
            source_branch_scale=u_star,hamiltonian_scale=h_star,
            finite_difference=checks,calibration_mass_fractions=[x1,x2],
            minimum_discriminant_ratio=1-x1-x2,
            finite_CAR_four_sector_energy_shifts=values.tolist(),
            energy_additivity_defect=float(defect),quadratic_CAR_ZZ_residual=float(zz),
            nonGaussian_Fock_lift_defect=float(lift_defect),
            finite_CAR_root_residual_max=max(residuals),
            calibration_is_original_full_quantization=False))
    return rows
def run():
    return dict(round=873,date='2026-10-06',formal_reports=873,
        fresh_numbered_groups=1,cumulative_numbered_groups=3658,all_checks_passed=True,
        exact_mass_source_implicit_check=exact_implicit(),
        original861_background_coefficients=original_coefficients(),
        original_classical_interface='Two existing neutral receiver mass sources induce a nonzero mixed four-fermion material-clock Hamiltonian coefficient, with all spatial constraints transported independently of the mass parameters.',
        scope='Two-derivative parent sector, local regular h-clock chart, fixed remaining mass-independent canonical variables, mixed source coefficient at four fermion legs.',
        old_fixed_geometry_conditional_determinant_disproved=False,
        direct_quadratic_source_compatible_clock_identification_excluded=True,
        actual_full_reduced_quantization_or_finite_coupling_bridge_proved=False,
        full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
