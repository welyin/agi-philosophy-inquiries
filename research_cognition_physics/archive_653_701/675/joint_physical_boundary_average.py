"""675: physical scalar marginal, original nonflat sphere contraction and ground limit.

Exact sphere identity inherited from657. Numeric minors validate its pullback,
not the complete 128-mode integral or positivity of the full Gauss kernel.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as base
import joint_auxiliary_reflection_gluing as sphere
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_physical_boundary_average_results.json'
old=base.old

def pf(a):
    return 1.+0j if len(a)==0 else base.pf(a)

def background(seed,scale):
    links=np.array([[base.prior.rep(*base.prior.group(seed+8*mu+i,scale)) for i in range(4)] for mu in range(2)])
    return links,base.kernel(links)

def physical_factor(u,v,mat,lam=.37):
    jm,jp,m,mb,pair=mat;r=len(m)//2;rv=v.shape[1]
    w=jm.conj().T@v;kl=jp.conj().T@v
    l=np.zeros((2*r,rv+r),complex)
    l[:r,:rv]=w;l[r:,rv:]=np.eye(r)
    n0=np.block([[np.zeros((rv,rv)),-kl.T],[kl,np.zeros((r,r))]])
    nw=n0+lam*l.T@pair@l
    phase=np.linalg.det(np.column_stack((u,v)))
    return pf(nw)/phase,nw

def original_factor_check():
    rows=[];rng=np.random.default_rng(67531)
    phis=base.mass.car.PHI.copy()
    for seed,scale in ((67531,.28),(67351,.9)):
        _,(u,v,d,h,gap)=background(seed,scale)
        values=[];constant=None
        for _ in range(2):
            e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
            mat=base.fixed_matrices(e,phis)
            physical,nw=physical_factor(u,v,mat)
            aux=pf(u.T@mat[2]@u)
            full=base.regular(u,v,d,mat)['weight']
            denom=max(abs(full),abs(aux*physical),1e-280)
            rel=float(abs(full-aux*physical)/denom)
            assert rel<2e-9 and abs(aux)<1+2e-12
            if constant is None:constant=physical
            assert abs(physical-constant)<1e-13*max(abs(constant),1e-200)
            pairing_norm=float(np.linalg.norm(mat[-1],2))
            bound_log=(len(nw)/2)*np.log1p(.37*pairing_norm)
            assert np.log(max(abs(full),1e-280))<=bound_log+1e-10
            values.append(dict(full_weight=old.cpair(full),auxiliary=old.cpair(aux),
                E_independent_physical_factor=old.cpair(physical),factor_relative=rel,
                original_pairing_norm=pairing_norm,log_scalar_bound=float(bound_log)))
        rows.append(dict(seed=seed,scale=scale,Wilson_gap=gap,ru=u.shape[1],rv=v.shape[1],values=values))
    return dict(rows=rows,no_inverse_mass_or_Dirac_or_weight=True)

def creation_quadratic(z,vec):
    return sum(z[i,j]*sphere.create_pair(vec,i,j) for i in range(len(z)) for j in range(i+1,len(z)))

def exact_minor(coeff):
    """Pull657 polynomial to8 modes across four original sites."""
    n=len(coeff[0][0]);vac=np.zeros(1<<n,complex);vac[0]=1
    def cx(x,vec):
        return sum(creation_quadratic(z,creation_quadratic(z,vec)) for z in coeff[x])
    first=[cx(x,vac) for x in range(len(coeff))]
    fourth=sum(first)/20
    eighth=sum(cx(x,first[x]) for x in range(len(coeff)))/960
    for x in range(len(coeff)):
        for y in range(x):eighth+=cx(x,first[y])/400
    return vac+fourth+eighth

def independent_cubature(coeff):
    """Integrate the degree4 Pf8 exactly, without global product quadrature.

    Pure one-site degree4 uses the inherited1044-point sphere rule. Two-site
    degree2+2 uses coordinate axes; sign pairing removes degree1+3 terms.
    All other degree partitions have an odd site and integrate to zero.
    """
    d=10;sites=len(coeff);value=0j
    for x in range(sites):
        for a in range(d):value+=2*pf(coeff[x][a])/(d*(d+2))
        for signs in itertools.product((-1,1),repeat=d):
            a=np.einsum('a,aij->ij',np.array(signs)/np.sqrt(d),coeff[x])
            value+=pf(a)*d/((d+2)*2**d)
    for x in range(sites):
        for y in range(x):
            for a in coeff[x]:
                for b in coeff[y]:
                    cross=(pf(a+b)+pf(-a+b))/2-pf(a)-pf(b)
                    value+=cross/d**2
    return value

def contraction_check():
    rows=[];rng=np.random.default_rng(67544)
    for seed,scale in ((67531,.28),(67351,.9)):
        links,(u,v,d,h,gap)=background(seed,scale)
        # Original spectral subspace; choose8 orthonormal linear combinations.
        z=rng.normal(size=(u.shape[1],8))+1j*rng.normal(size=(u.shape[1],8))
        mixing=np.linalg.qr(z)[0];small=u@mixing
        coeff=[]
        for x in range(4):
            ux=small[64*x:64*(x+1)]
            coeff.append(np.array([ux.T@np.kron(base.internal.B,t)@ux for t in base.internal.T]))
        coeff=np.array(coeff)
        polynomial=exact_minor(coeff)
        cubature=independent_cubature(coeff)
        err=float(abs(polynomial[-1]-cubature))
        assert err<2e-15
        # A4 coefficient from direct second moments, no quartic quadrature.
        a4=sum(pf(a[:4,:4]) for block in coeff for a in block)/10
        err4=float(abs(polynomial[15]-a4));assert err4<2e-15
        # Exterior rephasing is retained by exact integration.
        raw=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));w=np.linalg.qr(raw)[0]
        rotated=np.array([[w.T@a@w for a in block] for block in coeff])
        rotated_value=exact_minor(rotated)[-1]
        phase_error=float(abs(rotated_value-np.linalg.det(w)*polynomial[-1]))
        assert phase_error<2e-15
        # Degree3 Pf6 cancels under simultaneous E -> -E; no sector claim.
        e=rng.normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
        a6=np.einsum('xa,xaij->ij',e,coeff)[:6,:6]
        odd_average=(pf(a6)+pf(-a6))/2
        assert abs(odd_average)<2e-15
        rows.append(dict(seed=seed,scale=scale,Wilson_gap=gap,original_chiral_dimension=u.shape[1],
            checked_minor_modes=8,retained_internal_channels=16,original_sites=4,
            exact_eight_mode_coefficient=old.cpair(polynomial[-1]),
            independent_cubature_coefficient=old.cpair(cubature),eight_mode_error=err,
            four_mode_error=err4,frame_phase_error=phase_error,
            degree_six_antithetic_cancellation=old.cpair(odd_average)))
    return dict(rows=rows,original_128_mode_integral_not_evaluated=True,
        full_S9_moment_formula_reused_from657=True,
        no_minor_sign_extrapolated_to_full_physical_kernel=True)

def run():
    pilot=json.loads((HERE/'round675_drafts/ground_haar_probe_results.json').read_text('utf8'))
    deps=('joint_gauss_boundary_functional.py','joint_auxiliary_reflection_gluing.py',
          'research_note_591.md','research_note_657.md','research_note_673.md','research_note_674.md',
          'round675_drafts/ground_haar_entry.md','round675_drafts/ground_haar_probe_results.json')
    return dict(date='2026-10-02',round=675,tests_run=2,failures=0,errors=0,
        full_original_E_factorization=original_factor_check(),nonflat_exact_sphere_pullback=contraction_check(),
        prior_pilot=dict(samples=pilot['iid_base_samples'],
            largest_absolute_sample_share=pilot['largest_absolute_sample_share'],
            numerical_sign_unresolved=True),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Complete physical scalar marginal keeps auxiliary S9 and both original Haar boundaries. Exact657 local polynomial pulled into original nonflat chiral subspace, with selection rule and inverse-free physical factor.591 spectral gap yields controlled fixed-lambda ground limit and a necessary averaged-kernel positivity test. Only4/8-mode contractions evaluated exactly; no full integral sign, RP, normalization, original HF time, continuum or quantum GR identification.',
        all_checks_passed=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=675,tests_run=2,all_checks_passed=True)))
