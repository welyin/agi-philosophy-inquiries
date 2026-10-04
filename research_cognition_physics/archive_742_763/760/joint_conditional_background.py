"""760: conditional amplitudes of the original mixed Gauss process.
Numerical checks calibrate identities; they are not full graph propagation.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_record_mass_feedback as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_conditional_background_results.json'

def jets(x):
    rng=np.random.default_rng(760)
    A=[rng.normal(size=(4,3))+1j*rng.normal(size=(4,3)) for _ in range(3)]
    e=np.exp(-x*x/2)
    p=A[0]+x*A[1]+x*x*A[2]
    dp=A[1]+2*x*A[2]
    return e*p,e*(dp-x*p),e*(2*A[2]-2*x*dp+(x*x-1)*p)

def factors(w,dw,ddw):
    n=float(np.vdot(w,w).real)
    dn=2*float(np.vdot(w,dw).real)
    ddn=2*float((np.vdot(dw,dw)+np.vdot(w,ddw)).real)
    s=np.sqrt(n); u=dn/(2*n); v=ddn/(2*n)-dn*dn/(4*n*n)
    phi=w/s
    dphi=dw/s-u*phi
    ddphi=ddw/s-2*u*dw/s+(2*u*u-v)*phi
    return n,dn,ddn,s,u,v,phi,dphi,ddphi

def curved_equations_check():
    rng=np.random.default_rng(7601)
    q=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)); b0=(q+q.conj().T)/2
    q=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)); m=(q+q.conj().T)/2
    errors=dict(continuity=0.,conditional_equation=0.,kinetic_decomposition=0.,source_current=0.)
    for hb in (.13,.8):
        for x in (-.71,-.2,.17,.63):
            w,dw,ddw=jets(x)
            n,dn,ddn,s,u,v,ph,dph,ddph=factors(w,dw,ddw)
            # Original H5 x5 subalgebra: Delta f=(1+x^2/6)f''+5x f'/6.
            a=1+x*x/6; drift=5*x/6; V=.4+x*x+x**4
            B=b0+x*m
            Hw=-hb*hb/2*(a*ddw+drift*dw)+V*w+B@w
            wt=-1j*Hw/hb
            nt=2*float(np.vdot(w,wt).real)
            divj=hb*(a*np.vdot(w,ddw).imag+drift*np.vdot(w,dw).imag)
            errors['continuity']=max(errors['continuity'],abs(nt+divj))
            pht=wt/s-ph*nt/(2*n)
            rhs=-hb*hb/2*(a*ddph+drift*dph+2*a*u*dph+(a*v+drift*u)*ph)
            rhs+=V*ph+B@ph+1j*hb*divj/(2*n)*ph
            errors['conditional_equation']=max(errors['conditional_equation'],float(np.max(abs(1j*hb*pht-rhs))))
            ar=hb*float(np.vdot(ph,dph).imag)
            g=float(np.vdot(dph,dph).real-abs(np.vdot(ph,dph))**2)
            lhs=hb*hb/2*a*float(np.vdot(dw,dw).real)
            rhs=hb*hb*a*dn*dn/(8*n)+n*a*ar*ar/2+n*hb*hb*a*g/2
            errors['kinetic_decomposition']=max(errors['kinetic_decomposition'],abs(lhs-rhs))
            c=float(np.vdot(ph,m@ph).real)
            eta=hb*float(np.vdot(ph,(m-c*np.eye(4))@dph).imag)
            direct=hb*float(np.vdot(w,m@dw).imag)
            errors['source_current']=max(errors['source_current'],abs(direct-n*(ar*c+eta)))
    assert max(errors.values())<2e-12,errors
    return dict(samples=8,maximum_absolute_errors=errors,
                original_curved_kinetic_coefficients=True,
                matrix_jet_calibration_not_full_Gauss_evolution=True)

def original_phase_check():
    # Same inherited ten-mode coefficient fixture as759, not a new trace theorem.
    h,d=old.mass_x(np.array([0.,.51,0.,0.,.39]))
    ids=[0,4,8,12,14,16,24,25,30,31];ix=np.ix_(ids,ids);h,d=h[ix],d[ix]
    states=[b for b in range(1<<10) if (b&63).bit_count()%3==0]
    pos={b:i for i,b in enumerate(states)};dim=len(states)
    M=np.zeros((dim,dim),complex)
    for j,b in enumerate(states):
        for a,c in old.car.quadratic({b:1.+0j},h,d).items():
            assert a in pos
            M[pos[a],j]=c
    assert np.max(abs(M-M.conj().T))<1e-14
    kappa=float(np.vdot(M,M).real/dim)
    evals,evecs=np.linalg.eigh(M);base=np.eye(dim)/np.sqrt(dim)
    rows=[]
    for hb in (1.,.1,.01):
        phase=np.exp(-1j*.173*evals/hb)
        ph=(evecs*phase)@evecs.conj().T/np.sqrt(dim)
        dp=-1j*(M@ph)/hb
        sigma=ph@ph.conj().T
        current=hb*float(np.vdot(ph,dp).imag)
        source_current=hb*float(np.vdot(ph,M@dp).imag)
        weighted_metric=hb*hb*float(np.vdot(dp,dp).real-abs(np.vdot(ph,dp))**2)
        sigma_error=float(np.max(abs(sigma-np.eye(dim)/dim)))
        # The cross-configuration kernel changes although its diagonal does not.
        kernel_difference=float(np.linalg.norm(ph@base.conj().T-np.eye(dim)/dim))
        assert sigma_error<1e-14 and abs(current)<2e-13
        assert abs(source_current+kappa)<2e-13 and abs(weighted_metric-kappa)<2e-13
        assert kernel_difference>1e-3
        rows.append(dict(hbar=hb,conditional_density_error=sigma_error,
                         total_current=current,source_current=source_current,
                         hbar_squared_conditional_metric=weighted_metric,
                         off_diagonal_kernel_difference=kernel_difference))
    return dict(dimension=dim,original_CAR_mass_variance=kappa,rows=rows,
                same_conditional_density_and_current_different_source_current=True,
                finite_original_coefficient_calibration_not_full_periodic_simulation=True)

def original_record_check():
    errors=dict(probability_partition=0.,conditional_amplitude=0.,
                source_moment=0.,kinetic_injection=0.)
    rng=np.random.default_rng(7603)
    q=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4));C=(q+q.conj().T)/2
    for x in (-.51,-.13,.37,.81):
        w,dw,ddw=jets(x)
        n,_,_,_,_,_,ph,_,_=factors(w,dw,ddw)
        a=1+x*x/6; s=np.sqrt(2)*x/np.sqrt(a); ds=np.sqrt(2)/a**1.5
        hb=.21; tot_n=0.; tot_source=0.; after_T=0.; slope=0.
        for r in (-1,1):
            L=np.sqrt(.5+r*np.sin(s)/4)
            dL=r*np.cos(s)*ds/(8*L)
            wr=L*w; dwr=dL*w+L*dw
            nr=float(np.vdot(wr,wr).real)
            phr=wr/np.sqrt(nr)
            errors['conditional_amplitude']=max(errors['conditional_amplitude'],float(np.max(abs(phr-ph))))
            tot_n+=nr
            tot_source+=float(np.vdot(wr,C@C@wr).real)
            after_T+=hb*hb/2*a*float(np.vdot(dwr,dwr).real)
            slope+=dL*dL
        errors['probability_partition']=max(errors['probability_partition'],abs(tot_n-n))
        errors['source_moment']=max(errors['source_moment'],abs(tot_source-float(np.vdot(w,C@C@w).real)))
        before_T=hb*hb/2*a*float(np.vdot(dw,dw).real)
        errors['kinetic_injection']=max(errors['kinetic_injection'],abs(after_T-before_T-hb*hb/2*n*a*slope))
    assert max(errors.values())<3e-12,errors
    return dict(samples=4,maximum_absolute_errors=errors,
                original_sin_s_positive_Kraus=True,
                terminal_apparatus_autonomy_not_claimed=True,
                multi_time_equivalence_is_analytic_not_a_full_time_simulation=True)

def run():
    a=curved_equations_check();b=original_phase_check();c=original_record_check()
    deps=('research_note_643.md','research_note_704.md','research_note_706.md',
          'research_note_727.md','research_note_753.md','research_note_757.md',
          'research_note_758.md','research_note_759.md','joint_record_mass_feedback.py')
    return dict(round=760,tests_run=3,failures=0,errors=0,curved_conditional_equations=a,
                original_mass_phase_witness=b,original_record_and_resource=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='An exact conditional-amplitude representation of the existing finite-graph mixed Gauss process, retaining its full quantum slow variables. The scalar density, conditional matrix and total current alone lose a joint original-mass source current; compact finite-energy Gauss states give a witness. No spectral gap is needed for the exact representation. This does not give a reduced-complexity macroscopic approximation, a physical spatial metric dynamics, a continuum limit, or quantum gravity.')
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))

