"""884: exact original electron block and time-smeared source obstruction.
All early working evidence is preserved. The limit proof is in research_note884;
quadrature only calibrates the explicit finite integral, not the proof.
"""
from pathlib import Path
import argparse,json
import numpy as np
import cut_source_variance_probe as probe
import cut_time_smeared_lower_bound as working
old=probe.old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'cut_source_variance_certificate_results.json'
ELECTRON=[26,27,28,29]

def low_data(N,w,m2):
    u,wp=probe.solve_u(N,w)
    d=N/2-.25-u/2
    A=d+w;R=np.sqrt(A*A+4*m2)
    delta=2*m2/(R+A)
    mix=delta/R
    energy=w+delta
    slope=(1-mix)*wp+mix/2
    return dict(u=u,wp=wp,d=d,energy=energy,slope=slope,delta=delta,mix=mix)

def occupation_variance(E):
    r=np.exp(-old.BETA*abs(E))
    return r/(1+r)**2

def integrate(N,m2,order=48):
    x,weights=np.polynomial.legendre.leggauss(order)
    ws=.75*x+1.25;weights=.75*weights
    total=0.
    max_w_slope_error=0.;max_energy_error=0.;max_slope_ratio_error=0.
    for w,weight in zip(ws,weights):
        r=low_data(N,w,m2)
        total+=weight*12*r['slope']**2/r['wp']*occupation_variance(r['energy'])
        max_w_slope_error=max(max_w_slope_error,abs(r['wp']/(w/old.ETA*np.log(N)**2)-1))
        max_energy_error=max(max_energy_error,abs(r['energy']-w))
        max_slope_ratio_error=max(max_slope_ratio_error,abs(r['slope']/r['wp']-1))
    l=low_data(N,.5,m2);h=low_data(N,2.,m2)
    return dict(N=N,electron_time_smeared_variance_lower_bound=total,
        bound_over_logN_squared=total/np.log(N)**2,holonomy_strip_width=(h['u']-l['u'])/6,
        max_original_mass_energy_correction=max_energy_error,
        max_original_mass_slope_ratio_correction=max_slope_ratio_error,
        max_cut_slope_asymptotic_relative_error=max_w_slope_error)

def matrix_audit(m,q,mass):
    M=np.exp(old.XI)*m.MASS;m2=abs(mass)**2
    C=np.zeros((64,4),complex)
    vplus=np.array([1.,1.])/np.sqrt(2);vminus=np.array([1.,-1.])/np.sqrt(2)
    C[26:28,0]=vplus;C[28:30,1]=vplus
    C[26:28,2]=vminus;C[28:30,3]=vminus
    outside=np.eye(64)-C@C.conj().T
    assert np.max(abs(outside@M@C))<1e-13
    rows=[]
    for N in (17,65,257):
        for w in (.5,1.,2.):
            r=low_data(N,w,m2)
            alpha=-(.5-r['u'])/6;k=N//2
            h=old.H(m,q,N,[k,0,0],alpha)
            g=m.GAMMA[0]@np.diag(q*probe.source(N,k+alpha*q))
            blocks=np.zeros((4,4),complex);sources=np.zeros((4,4),complex)
            for b,s in enumerate((1,-1)):
                j=2*b
                blocks[j:j+2,j:j+2]=np.array([[-s*r['d'],mass],[mass.conjugate(),s*w]])
                sources[j:j+2,j:j+2]=np.diag([3*s,6*s*r['wp']])
            err=float(np.max(abs(C.conj().T@h@C-blocks)))
            gerr=float(np.max(abs(C.conj().T@g@C-sources)))
            leakage=float(np.max(abs(outside@h@C)))
            eigen_error=0.;slope_error=0.
            for s in (1,-1):
                B=np.array([[-s*r['d'],mass],[mass.conjugate(),s*w]])
                G=np.diag([3*s,6*s*r['wp']])
                e,v=np.linalg.eigh(B);j=int(np.argmin(abs(e-s*r['energy'])))
                eigen_error=max(eigen_error,float(abs(e[j]-s*r['energy'])))
                slope_error=max(slope_error,float(abs(np.vdot(v[:,j],G@v[:,j]).real-s*6*r['slope'])))
            all_lower=0.
            for kk in (k,-k):
                hh=old.H(m,q,N,[kk,0,0],alpha)
                gg=m.GAMMA[0]@np.diag(q*probe.source(N,kk+alpha*q))
                all_lower+=working.diagonal_lower(hh,gg)/2
            selected=72*r['slope']**2*occupation_variance(r['energy'])
            assert max(err,gerr,leakage,eigen_error,slope_error)<3e-8
            assert all_lower>=selected-1e-6
            rows.append(dict(N=N,w=w,original_block_error=err,source_block_error=gerr,
                invariant_block_leakage=leakage,exact_eigen_error=eigen_error,
                Hellmann_Feynman_error=slope_error,full64_diagonal_lower=all_lower,
                selected_electron_lower=selected))
    return rows

def run():
    m=old.load();q,_=old.charges(m)
    M=np.exp(old.XI)*m.MASS
    mass=complex(M[26,28]);m2=abs(mass)**2
    assert abs(M[27,29]-mass)<1e-14
    audit=matrix_audit(m,q,mass)
    rows=[integrate(N,m2) for N in (17,33,65,129,257,1025,4097,65537,1048577)]
    checks=[]
    for N in (17,257,4097):
        a=integrate(N,m2,32);b=integrate(N,m2,64)
        err=abs(a['electron_time_smeared_variance_lower_bound']-b['electron_time_smeared_variance_lower_bound'])
        checks.append(dict(N=N,quadrature32_vs64_absolute_error=err))
        assert err<2e-9
    coeffs=[]
    for order in (32,64):
        x,w=np.polynomial.legendre.leggauss(order);v=.75*x+1.25
        coeffs.append(float(12/old.ETA*np.dot(.75*w,v*occupation_variance(v))))
    analytic_positive_lower=12/old.ETA*((2**2-.5**2)/2)*occupation_variance(2.)
    assert abs(coeffs[0]-coeffs[1])<1e-11 and coeffs[1]>analytic_positive_lower>0
    return dict(round=884,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3669,
        argument_scope='For883 smooth-cut hopping with a fixed flat hypercharge average of positive density at -1/12 and its declared conditional free thermal preparations, a fixed smooth time-smeared current variance diverges at least as log(N)^2. The original continuous free theory has a finite corresponding averaged variance. This excludes that joint source limit, not the full dynamical Gauss model or all spectral regulators.',
        original_electron_particle_indices=ELECTRON,
        original_electron_mass_real=mass.real,original_electron_mass_imag=mass.imag,
        original_electron_mass_squared=m2,integer_hypercharges=dict(left=-3,right=-6),
        Nambu_counting='1/2 for complete opposite-momentum pair; selected four low BdG directions give72 slope^2 n(1-n).',
        exact_original_matrix_audit=audit,asymptotic_lower_bound_rows=rows,
        numerical_limit_coefficient=coeffs[1],
        conservative_analytic_positive_coefficient_lower=analytic_positive_lower,
        coefficient_formula='(12/eta) integral_[0.5,2] w exp(-beta w)/(1+exp(-beta w))^2 dw',
        quadrature_consistency=checks,
        fixed_real_compact_smooth_time_test=True,time_integral=1,
        whole_continuous_time_smeared_variance_finite_proved=True,
        logarithmic_source_obstruction_proved=True,
        full_dynamic_Gauss_measure_derived=False,
        all_spectral_regularizations_excluded=False,
        source_Hessian_can_cancel_positive_observable_variance=False,
        full_interacting_Q_E_bridge_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
