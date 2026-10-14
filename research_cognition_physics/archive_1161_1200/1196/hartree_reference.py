"""Independent Hartree reference for the fixed 1178 parent; not a GR certificate.

Default: calculate and print JSON without writing. --write: exclusively create
hartree_reference_result.json; an existing result is never overwritten.
Only Python/NumPy are used. Exact fractions certify the global contraction and
strict positivity bounds. Spectral/finite-sum decimal values are diagnostics.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np


def calculate():
    nlat=65
    ell=F(16,65)
    eps,nu,omega,m=F(3,20),F(7,20),F(1,2),F(2)
    qcert=F(6,25)
    # |b| >= lower_b; |b'| <= bp_without_sqrt3/sqrt(3).
    lower_b=2*eps*(1-nu)*m/(1+nu)+eps*omega
    bp_without_sqrt3=4*m*(eps-nu)**2/(3*(1-nu)**3)
    qsq=omega**4*bp_without_sqrt3**2/(3*(omega**2+lower_b**2)**3)
    assert qsq < qcert*qcert
    # Two transverse momenta +/-16, k_y >= (999/1000)/ell,
    # since cos(pi/130)>999/1000 follows from pi^2<10.
    a2=(m*(1+eps))**2+((1+nu)/ell)**2
    b2=((1-nu)*F(999,1000)/ell)**2
    pair_floor=b2/(2*(a2+b2))
    corr_floor=pair_floor/(nlat*nlat)
    fock_floor=corr_floor**8
    assert F(10,2*130*130)<F(1,1000)

    p=np.arange(-32,33)
    kval=np.sin(2*np.pi*p/nlat)/float(ell)
    ksq=kval[:,None,None]**2+kval[None,:,None]**2+kval[None,None,:]**2
    ep,nn,om,mm=map(float,(eps,nu,omega,m))

    def data(z):
        lapse=1+ep*z
        speed=1+nn*z
        energy=np.sqrt((mm*lapse)**2+speed*speed*ksq)
        em=-2*mm*mm*lapse*np.mean(1/energy)
        et=-2*speed*np.mean(ksq/energy)
        b=ep*em+nn*et-ep*om
        bp=-2*mm*mm*(ep-nn)**2*np.mean(ksq/energy**3)
        mapped=-b/np.sqrt(om*om+b*b)
        fp=-om*om*bp/(om*om+b*b)**1.5
        functional=-2*np.mean(energy)-om*np.sqrt(1-z*z)-om-ep*om*z
        return lapse,speed,energy,em,et,b,bp,mapped,fp,functional

    z=0.0
    for iteration in range(100):
        dz=data(z)
        updated=dz[7]
        if abs(updated-z)<2e-15:
            z=float(updated)
            break
        z=float(updated)
    else:
        raise RuntimeError('Fixed-point iteration failed')
    lapse,speed,energy,em,et,b,bp,mapped,fp,functional=data(z)
    assert abs(mapped-z)<1e-13
    xs=np.sqrt(1-z*z)
    # Hartree full energy variation including reciprocal -epsilon*Omega*z_s*x_t.
    first_variation=b+om*z/xs
    assert abs(first_variation)<2e-12
    d_source=ep*em+nn*et
    assert abs(d_source-(-2*np.mean((ep*mm*mm*lapse+nn*speed*ksq)/energy)))<2e-13

    eye=np.eye(2)
    sx=np.array([[0,1],[1,0]],complex)
    sz=np.diag([1,-1])
    beta=np.kron(sz,eye)
    ax=np.kron(sx,sx)
    # Convention C_(i,a),(j,b)=<c^dagger_(j,b)c_(i,a)>.
    phase=np.exp(-2j*np.pi*p/nlat)[:,None,None]
    c0=np.eye(4)/2-mm*lapse*np.mean(1/energy)*beta/2
    c01=(np.mean(phase)*np.eye(4)/2
          -mm*lapse*np.mean(phase/energy)*beta/2
          -speed*np.mean(phase*kval[:,None,None]/energy)*ax/2)
    cov=np.block([[c0,c01],[c01.conj().T,c0]])
    eig,vec=np.linalg.eigh(cov)
    assert np.max(np.abs(cov-cov.conj().T))<1e-14
    assert 0<eig[0]<=eig[-1]<1
    modular=(vec*np.log((1-eig)/eig))@vec.conj().T
    # This is the source CAR region only: eight complex modes / 256 Fock states.
    mu_numeric=float(np.prod(np.minimum(eig,1-eig)))
    mod_eig,mod_vec=np.linalg.eigh(modular)
    rebuilt=(mod_vec*(1/(1+np.exp(mod_eig))))@mod_vec.conj().T
    assert np.linalg.norm(rebuilt-cov,2)<2e-14

    volume=nlat**3
    bsource=2*eps*m+6*nu/ell+4*F(1,5)/ell
    bprobe=2*eps*(F(3,10)+F(1,2))+18*nu/ell+12*F(1,5)/ell
    w=2*volume*bsource+2*volume*eps*omega
    lprobe=volume*bprobe
    return {
      'status':'PASS_REFERENCE_SCOPE_ONLY',
      'parent':'1178, E0=0, full source CAR occupations retained',
      'parameters':{'Nlat':nlat,'ell':str(ell),'epsilon':str(eps),'nu':str(nu),
                    'Omega':str(omega),'mass_source':str(m),'g':str(F(1,5))},
      'analytic_global_certificate':{
        'interval':'[-1,1]','b_absolute_lower_exact':str(lower_b),
        'bprime_bound_coefficient_over_sqrt3_exact':str(bp_without_sqrt3),
        'map_derivative_squared_upper_exact':str(qsq),
        'strict_contraction_upper_exact':str(qcert),
        'comparison_margin_exact':str(qcert*qcert-qsq),
        'unique_isotropic_fixed_point':True,
        'claim_global_minimum_all_nonuniform_states':False},
      'fixed_point_floating_diagnostics':{
        'iterations':iteration+1,'z_s':z,'x_s':float(xs),'z_t':0.0,'x_t':1.0,
        'N':float(lapse),'v':float(speed),'b':float(b),'bprime':float(bp),
        'map_derivative':float(fp),'fixed_point_residual':float(mapped-z),
        'energy_first_variation':float(first_variation),
        'source_mass_energy_per_site':float(em),'source_hopping_energy_per_site':float(et),
        'complete_Hartree_energy_per_site_excluding_global_baseline':float(functional),
        'occupied_rank_each_momentum':2,'source_particles_each_copy':2*volume},
      'two_site_source_reference':{
        'sites':[[0,0,0],[1,0,0]],'complex_CAR_modes':8,'Fock_dimension':256,
        'covariance_eigenvalues':eig.tolist(),
        'crossblock_operator_norm':float(np.linalg.norm(c01,2)),
        'crossblock_Frobenius_norm':float(np.linalg.norm(c01)),
        'modular_reconstruction_residual':float(np.linalg.norm(rebuilt-cov,2)),
        'Fock_min_eigenvalue_numeric':mu_numeric,
        'analytic_covariance_lower_and_complement_lower_exact':str(corr_floor),
        'analytic_Fock_lower_exact':str(fock_floor),
        'analytic_Fock_lower_float':float(fock_floor),
        'is_full_parent_region':False},
      'finite_K_reuse':{
        'whole_cluster_W_norm_bound_exact':str(w),
        'probe_coupling_norm_bound_exact':str(lprobe),
        'alpha_bound':'(exp(12*w*T)-1)/K',
        'source_region_trace_distance_bound':'sqrt(alpha)',
        'probe_isometry_error_bound':'4*L*T*exp(6*w*T)/sqrt(K)',
        'old_1178_site_independence_bound_reused':False,
        'log_bound_automatic_from_trace_distance':False},
      'scope':{
        'target_boost_used_to_choose_reference':False,
        'full_finite_K_ground_or_Gibbs_state_certified':False,
        'all_parent_region_faithfulness_certified':False,
        'original_probe_reference_implementation_certified':False,
        'full_source_boost_identity_certified':False,
        'new_axiom_adopted':False,'GR_generated':False}
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate()
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    output=json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if args.write:
        path=Path(__file__).with_name('hartree_reference_result.json')
        with path.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(output)
    print(output,end='')

if __name__=='__main__':
    main()
