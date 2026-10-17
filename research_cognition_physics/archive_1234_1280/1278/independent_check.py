"""1278 independent audit: original (unnormalized) pair variables, exact Taylor bounds.
No author code is imported. Matrix checks are finite algebra diagnostics only.
Default is read-only comparison; --write creates the result exclusively.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in HERE.parents if (p / 'research_cognition_physics').is_dir())
T = F(1, 8)
ORDER = 64


def need(ok, text):
    if not ok:
        raise AssertionError(text)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rat(x):
    x = F(x)
    return {'numerator': x.numerator, 'denominator': x.denominator}


def outward(lo, hi, digits=24):
    scale = 10 ** digits
    a = lo.numerator * scale // lo.denominator
    b = -((-hi.numerator * scale) // hi.denominator)
    return [rat(F(a, scale)), rat(F(b, scale))]


def clean(x):
    return float(format(float(x), '.12g'))


def conv(c, i, j, n):
    return sum((c[i][m] * c[j][n-m] for m in range(n+1)), F(0))


def series(program):
    # Variables are rx,ry,rz,gx,gy,gz, rather than the author's normalized u.
    c = [[F(v)] for v in (0, 0, F(-3,20), 1, 0, 0)]
    # Sparse monomials: coefficient, first variable, optional second variable.
    rhs = [
        [(F(-2),5,1)],
        [(F(2),5,0),(F(-2*program),2,None)],
        [(F(2*program),1,None)],
        [(F(-2),2,4)],
        [(F(2),2,3),(F(2),5,None)],
        [(F(-2),4,None)],
    ]
    for n in range(ORDER):
        next_coeff = []
        for equation in rhs:
            q = F(0)
            for weight, i, j in equation:
                q += weight * (c[i][n] if j is None else conv(c,i,j,n))
            next_coeff.append(q/F(n+1))
        for i in range(6):
            c[i].append(next_coeff[i])
    # The positive majorant y'=4y^2,y(0)=1 has coefficients 4^n.
    # Every row has coefficient absolute sum <=4; linear terms use the
    # constant series 1 <= y. This inductively controls ALL orders.
    need(all(sum(abs(q[0]) for q in eq) <= 4 for eq in rhs), 'analytic majorant rows')
    need(all(abs(c[i][n]) <= 4**n for i in range(6) for n in range(ORDER+1)),
         'independent coefficient recurrence')
    values = [sum((c[i][n]*T**n for n in range(ORDER+1)),F(0)) for i in range(6)]
    integral = sum((c[5][n]*T**(n+1)/F(n+1) for n in range(ORDER+1)),F(0))
    need(4*T == F(1,2), 'fixed radius')
    tail = (4*T)**(ORDER+1)/(1-4*T)
    need(tail == F(1,2**64), 'uniform coordinate tail')
    invariants = 0
    for n in range(ORDER+1):
        r2 = sum((conv(c,j,j,n) for j in range(3)),F(0))
        g2 = sum((conv(c,j,j,n) for j in range(3,6)),F(0))
        energy = program*c[0][n] - c[3][n] + conv(c,2,5,n)
        need(r2 == (F(9,400) if n == 0 else 0), 'pair Bloch norm')
        need(g2 == (F(1) if n == 0 else 0), 'relation Bloch norm')
        need(energy == (F(-1) if n == 0 else 0), 'interaction plus relation energy')
        invariants += 3
    if program == 0:
        need(c[2] == [F(-3,20)]+[F(0)]*ORDER, 'no exchange pair distribution')
    return c, values, integral, tail, invariants


def finite_bounds():
    zero, one = series(0), series(1)
    coordinate_tail = zero[3]
    integral_tail = T*coordinate_tail
    di = one[2]-zero[2]
    dr = F(1,5)*di
    dr_tail = F(2,5)*integral_tail
    need(-dr-dr_tail > F(1,10**7), 'finite non-rescaling ideal gap')
    need(one[0][5][4]-zero[0][5][4] == F(-1,10), 'first differential coefficient')
    heat = (one[1][2]+F(3,20))/2
    heat_tail = coordinate_tail/2
    need(heat-heat_tail > F(1,500), 'positive A energy transfer')
    # The full probe Fock norm is inherited as an operator bound, not an
    # old numerical Hartree/probe error for a different source preparation.
    volume, ell = 65**3, F(16,65)
    eps, nu = F(3,20), F(7,20)
    lp = volume*(2*eps*(F(3,10)+F(1,2))+18*nu/ell)
    w = 2*volume
    exponent = 6*w*T
    ceiling = -((-exponent.numerator)//exponent.denominator)
    need(lp < 2**23 and 4*lp*T < 2**22, 'full probe Fock prefactor')
    need(ceiling == 411938, 'mean-field exponential ceiling')
    # e<4, 10>2^3, hence eta <2^(22+2*ceil(6wT)-900000).
    power = 22+2*ceiling-3*300000
    need(power == -76102 and 2**40 > 10**12 and -power > 40,
         'finite-copy mean-field <1e-12')
    # Exact lattice k=sin(2pi/65)/ell; elementary outward pi interval and
    # sin(x)>=x-x^3/6, sin(x)<=x on this positive interval.
    xlo, xhi = 2*F(157,50)/65, 2*F(22,7)/65
    klo = (xlo-xhi**3/6)/ell
    khi = xhi/ell
    need(F(3,8) < klo < khi < F(2,5), 'exact lattice wave number interval')
    tau_max, xi_max = (1+eps)*T, (1+nu)*T
    p_error = F(1,10**10)+F(1,10**12)
    clock_arg = tau_max/5+2*p_error
    position_arg = khi*xi_max+4*p_error
    need(clock_arg < F(1,10) and position_arg < F(1,10),
         'error-expanded principal arcsin branches')
    # sqrt(1-y^2)>99/100 for |y|<1/10.
    need(F(99,100)**2 < 1-F(1,10)**2, 'inverse square-root bound')
    need(F(10)/F(99,100) < 11, 'clock inverse Lipschitz')
    need(F(4)/(klo*F(99,100)) < 14, 'signal inverse Lipschitz')
    pair_error = 2*(11+14)*p_error
    need(pair_error == F(101,20_000_000_000), 'two-program inverse error')
    robust_coarse = F(1,10**7)-pair_error
    need(robust_coarse > F(9,10**8), 'conditioned operational gap')
    # A same-t difference of xi-tau alone does not reject an arbitrary common
    # reparametrization. Match tau1(T)=tau0(t*) instead. The old curve obeys
    # dxi/dtau=v/N<=27/17; t* is an analysis parameter, not a new stop protocol.
    slope = F(27,17)
    curve_coefficient = nu-slope*eps
    need(curve_coefficient == F(19,170), 'equal-clock curve coefficient')
    curve_value = curve_coefficient*di
    curve_tail = 2*curve_coefficient*integral_tail
    curve_ideal_lower = -curve_value-curve_tail
    curve_error = 2*(14+slope*11)*p_error
    need(curve_error == F(1070,17)*p_error, 'curve inverse error')
    need(curve_ideal_lower > F(67,10**9), 'ideal equal-clock curve gap')
    need(curve_ideal_lower-curve_error > F(6,10**8), 'robust equal-clock curve gap')
    need(T+eps*(one[2]-integral_tail)>0 and di+2*integral_tail<0,
         'unique old-curve comparison t* in (0,T)')
    # Positivity domains from exact Bloch invariants, independent of samples.
    need(1-eps == F(17,20) and 1-nu == F(13,20), 'positive N and v')
    programs = []
    for k, entry in enumerate((zero,one)):
        programs.append({
            'program': k,
            'rz_interval': outward(entry[1][2]-coordinate_tail, entry[1][2]+coordinate_tail),
            'gz_interval': outward(entry[1][5]-coordinate_tail, entry[1][5]+coordinate_tail),
            'integral_gz_interval': outward(entry[2]-integral_tail, entry[2]+integral_tail),
            'tau_interval': outward(T+eps*(entry[2]-integral_tail),T+eps*(entry[2]+integral_tail)),
            'xi_interval': outward(T+nu*(entry[2]-integral_tail),T+nu*(entry[2]+integral_tail)),
        })
    return {
        'original_pair_variables_used': True,
        'order': ORDER,
        'T': rat(T),
        'coefficient_comparisons': 2*6*(ORDER+1),
        'exact_invariant_coefficient_checks': zero[4]+one[4],
        'coordinate_uniform_tail': rat(coordinate_tail),
        'integral_uniform_tail': rat(integral_tail),
        'programs': programs,
        'delta_integral_gz_interval': outward(di-2*integral_tail,di+2*integral_tail),
        'delta_xi_minus_tau_interval': outward(dr-dr_tail,dr+dr_tail),
        'positive_A_energy_transfer_interval': outward(heat-heat_tail,heat+heat_tail),
        'ideal_gap_strict_lower':rat(F(1,10**7)),
        'full_probe_L': rat(lp), 'cross_cluster_w':w,
        'copies':'10^600000','sqrt_copies_power10':300000,
        'mean_field_strict_upper_power2':power,
        'mean_field_used_upper':rat(F(1,10**12)),
        'instrument_preparation_contract':rat(F(1,10**10)),
        'exact_lattice_k_interval':[rat(klo),rat(khi)],
        'inverse_lipschitz':{'clock':11,'position':14},
        'error_expanded_inverse_arguments':[rat(clock_arg),rat(position_arg)],
        'two_program_inverse_error_upper':rat(pair_error),
        'coarse_robust_gap_lower':rat(robust_coarse),
        'reported_robust_gap_strict_lower':rat(F(9,10**8)),
        'non_reparametrization_certificate':{
          'old_curve_slope_upper':rat(slope),
          'coefficient_nu_minus_slope_epsilon':rat(curve_coefficient),
          'signed_witness_interval':outward(curve_value-curve_tail,curve_value+curve_tail),
          'inverse_error_upper':rat(curve_error),
          'ideal_gap_strict_lower':rat(F(67,10**9)),
          'robust_gap_strict_lower':rat(F(6,10**8)),
          'comparison_time_is_proof_parameter_not_new_stop_permission':True,
        },
        'N_domain':[rat(F(17,20)),rat(F(23,20))],
        'v_domain':[rat(F(13,20)),rat(F(27,20))],
    }


def matrix_diagnostics():
    eye2=np.eye(2,dtype=complex)
    x=np.array([[0,1],[1,0]],dtype=complex)
    z=np.diag([1,-1]).astype(complex)
    n=np.diag([0,1]).astype(complex)
    na=np.kron(n,eye2); nb=np.kron(eye2,n)
    nt=na+nb; d=na-nb
    ex=np.zeros((4,4),complex); ex[1,2]=ex[2,1]=1
    ey=(d@ex-ex@d)/(2j)
    p1=n.copy()
    expand=lambda ab, gg, pp: np.kron(np.kron(ab,gg),pp)
    dc=expand(d,eye2,eye2); zc=expand(np.eye(4),z,eye2)
    nc=expand(nt,eye2,eye2); pc=expand(np.eye(4),eye2,p1)
    h=expand(nt,eye2,eye2)+expand(ex,eye2,p1)-expand(np.eye(4),x,eye2)
    W=np.kron(dc,zc)+np.kron(zc,dc)
    comm=lambda a,b:a@b-b@a
    residual=max(np.max(np.abs(comm(h,nc))),np.max(np.abs(comm(h,pc))))
    need(residual==0, 'exact bare-energy and program conservation')
    need(np.max(np.abs(comm(d,ex)))==2, 'distribution is not conserved during exchange')
    need(np.max(np.abs(W-np.diag(np.diag(W))))==0, 'cross interaction diagonal')
    need(np.max(np.abs(np.diag(W)))==2, 'cross interaction norm 2')
    need(np.linalg.eigvalsh(h).min()>-1.00000000000001 and
         np.linalg.eigvalsh(h).max()<3.00000000000001,'local energy bound')
    # Whole finite-parent positive shift can conservatively include
    # 4 V K + V(K-1) + ||H_P0|| + L +1; no term is removed from the ledger.
    # Numerical commutator reconstruction of all six Hartree equations.
    ode_res=0.0
    samples=[(np.array([0.,0.,-.15]),np.array([1.,0.,0.])),
             (np.array([.03,.04,-.10]),np.array([.5,-.2,.3])),
             (np.array([-.04,.02,.07]),np.array([-.2,.6,-.1]))]
    observables=[np.kron(ex,eye2),np.kron(ey,eye2),np.kron(d,eye2)]
    observables += [np.kron(np.eye(4),a) for a in (x,np.array([[0,-1j],[1j,0]]),z)]
    for k in (0,1):
      for rr,gg in samples:
        rhoab=np.zeros((4,4),complex); rhoab[0,0]=.45;rhoab[3,3]=.10
        rhoab[1,1]=.225-rr[2]/2;rhoab[2,2]=.225+rr[2]/2
        rhoab += (rr[0]*ex+rr[1]*ey)/2
        y=np.array([[0,-1j],[1j,0]])
        rhog=(eye2+gg[0]*x+gg[1]*y+gg[2]*z)/2
        rho=np.kron(rhoab,rhog)
        hm=np.kron(nt+k*ex+gg[2]*d,eye2)+np.kron(np.eye(4),-x+rr[2]*z)
        dot=-1j*comm(hm,rho)
        actual=np.array([np.trace(o@dot).real for o in observables])
        expected=np.r_[2*np.cross([k,0,gg[2]],rr),2*np.cross([-1,0,rr[2]],gg)]
        ode_res=max(ode_res,float(np.max(np.abs(actual-expected))))
    need(ode_res<1e-14,'Hamiltonian-to-ODE diagnostic')
    # Genuine local mixed preparation, retained purification, no Gibbs-of-h claim.
    weights=np.array([.45,.30,.15,.10])
    purified=np.zeros((4,4),complex)
    for j in range(4):purified[j,j]=np.sqrt(weights[j])
    marginal=purified@purified.conj().T
    preparation_res=float(np.max(np.abs(marginal-np.diag(weights))))
    need(preparation_res<1e-14,'retained thermal purification')
    need(abs(np.trace(d@marginal).real+.15)<1e-14,'initial D=-3/20')
    need(abs(np.trace(nt@marginal).real-.65)<1e-14,'total excitation 13/20')
    return {
      'cluster_matrix_dimension':16,'two_cluster_interaction_dimension':256,
      'commutator_residual':clean(residual),'cross_interaction_norm':2,
      'local_h_eigenvalue_interval':[-1,3],
      'hartree_matrix_cases':6,'hartree_ode_residual':clean(ode_res),
      'thermal_purification_residual':clean(preparation_res),
      'role':'finite algebra diagnostics; no huge-copy matrix simulation',
    }


def output():
    files=[
      HERE/'admission.md', HERE/'fixed_spec.md', HERE/'source_baseline.json',
      ROOT/'research_cognition_physics/archive_1161_1200/1175/proof.md',
      ROOT/'research_cognition_physics/archive_1161_1200/1179/proof.md',
      ROOT/'research_cognition_physics/archive_1161_1200/_admission/after1174_direction/mean_field_bound.md',
      Path(__file__).resolve(),
    ]
    return {
      'round':1278,'status':'PASS_INDEPENDENT_WITH_SCOPE','independent_groups':1,
      'source_bindings':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files},
      'rational_certificate':finite_bounds(),'matrix_diagnostics':matrix_diagnostics(),
      'scope':{
        'same_rule_two_known_program_inputs':True,
        'unknown_probe_and_reference_retained_by_preparation_isometry':True,
        'background_and_purification_retained':True,
        'new_1e_minus_10_instrument_contract_is_conditional_input':True,
        'no_old_source_numeric_error_inherited':True,
        'no_continuum_lattice_error_needed_for_exact_endpoint_inversion':True,
        'finite_frequency_sampling_error_not_included':True,
        'huge_copy_count_is_not_engineering_feasibility':True,
        'no_universal_stress_GR_geometry_entropy_or_observer_temperature_claim':True,
      },
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=output()
    destination=HERE/'independent_results.json'
    if args.write:
      with destination.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(result,f,ensure_ascii=False,indent=2,sort_keys=True)
        f.write('\n')
      mode='WRITE_ONCE'
    else:
      expected=json.loads(destination.read_text(encoding='utf-8'))
      need(expected==result,'read-only result mismatch')
      mode='READ_ONLY_MATCH'
    print(json.dumps({'status':result['status'],'mode':mode,
                     'xi_minus_tau_gap_strict_lower':'9e-8',
                     'non_reparametrization_gap_strict_lower':'6e-8',
                     'results_sha256':sha(destination)},ensure_ascii=False))

if __name__=='__main__':
    main()
