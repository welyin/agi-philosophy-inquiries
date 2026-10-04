"""697: uniform-gap polar reduction and full-measure numerical diagnostic.

The gap and source identities are analytic. Monte Carlo errors are empirical,
not rigorous integration errors; no reflection-positivity theorem is claimed.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round697_drafts'))
import full_average_probe as probe
import importance_full_average as numerical
b=probe.b
TARGET=HERE/'joint_full_holonomy_reduction_results.json'


def polar(h,first,second):
    odd=np.diag(b.prior.rep(np.eye(3),-np.eye(2),1)).real<0
    w0=2*first*odd;w1=2*second*odd
    d=w0*w1+h;eta=d/abs(d)
    denom=np.sqrt(w0*w0+w1*w1+2+2*abs(d))
    v=np.block([[np.diag((w0+eta*w1)/denom),np.diag((-1-eta*h.conj())/denom)],
                [np.diag((h+eta)/denom),np.diag((w1+eta*w0)/denom)]])
    a=np.block([[np.diag(w0),-np.eye(16)],[np.diag(h),np.diag(w1)]])
    return a,v


def frames(v):
    # Column order: time, twofold spin, internal channel. Row order before
    # Hadamard: time, gamma4 eigenvalue, twofold spin, internal channel.
    ur=np.zeros((128,64),complex);vr=np.zeros_like(ur)
    for t in range(2):
        for s in range(2):
            for spin in range(2):
                cols=slice(32*s+16*spin,32*s+16*spin+16)
                top=slice(64*t+16*spin,64*t+16*spin+16)
                bottom=slice(64*t+32+16*spin,64*t+32+16*spin+16)
                ur[top,cols]=np.eye(16)*(t==s)/np.sqrt(2)
                vr[top,cols]=ur[top,cols]
                vr[bottom,cols]=v.conj().T[16*t:16*t+16,16*s:16*s+16]/np.sqrt(2)
                ur[bottom,cols]=-vr[bottom,cols]
    hd=np.array([[1,1],[1,-1]])/np.sqrt(2)
    rotation=np.kron(np.eye(2),np.kron(hd,np.eye(32)))
    return rotation@ur,rotation@vr


def algebra_checks():
    rng=np.random.default_rng(69781)
    angles=rng.uniform(0,2*np.pi,(3,4));hs,_=probe.torus(angles)
    fields=rng.normal(size=(3,2,10));fields/=np.linalg.norm(fields,axis=2)[:,:,None]
    rows=[]
    for first,second in ((False,False),(True,False),(False,True),(True,True)):
        for k,h in enumerate(hs):
            a,v=polar(h,first,second);u,p=frames(v)
            eig_u,eig_v,d,gap=probe.entry.kernel(first,second,np.eye(16),np.diag(h))
            gs=np.kron(np.eye(2),np.kron(b.internal.spin.G5,np.eye(16)))
            sign_h=gs@(2*d-np.eye(128))
            mat=b.fixed_matrices(fields[k],np.zeros((2,5)))
            jm,jp,m,_,_=mat
            t=np.zeros((32,32),complex)
            for site in range(2):t[16*site:16*site+16,16*site:16*site+16]=sum(x*y for x,y in zip(fields[k,site],b.internal.T))
            aux=(t+v.conj()@t@v.conj().T)/2
            det_frame=np.linalg.det(np.column_stack((u,p)))
            pf_u=b.pf(u.T@m@u)
            kl=jp.conj().T@p
            original=b.regular(eig_u,eig_v,d,mat,lam=0)['weight']
            reduced=np.linalg.det((np.eye(32)+v)/2)**2*np.linalg.det(aux)
            scalar=b.pf(u.T@m@u)*np.linalg.det(kl)/det_frame
            errors=dict(unitarity=float(np.max(abs(v.conj().T@v-np.eye(32)))),
                polar_hermitian=float(np.max(abs(v.conj().T@a-a.conj().T@v))),
                negative_frame=float(np.max(abs(sign_h@u+u))),
                positive_frame=float(np.max(abs(sign_h@p-p))),
                frame_determinant=float(abs(det_frame/np.linalg.det(v.conj().T)**2-1)),
                auxiliary_pf_absolute=float(abs(pf_u-np.linalg.det(aux))),
                physical_det_absolute=float(abs(np.linalg.det(kl)/det_frame-np.linalg.det((np.eye(32)+v)/2)**2)),
                original_relative=float(abs(original-reduced)/max(abs(original),abs(reduced),1e-28)),
                canonical_relative=float(abs(scalar-reduced)/max(abs(scalar),abs(reduced),1e-28)))
            assert max(errors.values())<1e-9,errors
            assert np.linalg.eigvalsh(v.conj().T@a).min()>0
            assert gap>=np.sqrt(2)-1-2e-13
            rows.append(dict(centres=[first,second],sample=k,gap=gap,errors=errors))
    # Deterministic endpoint checks of the exact singular-value formulas.
    spectra=[]
    for w0,w1 in ((0,0),(2,0),(0,2),(2,2)):
        z=np.exp(1j*np.linspace(0,2*np.pi,513));tr=w0*w0+w1*w1+2
        dsq=np.ones(len(z)) if w0*w1==0 else abs(w0*w1+z)**2
        low=np.sqrt((tr-np.sqrt(np.maximum(0,tr*tr-4*dsq)))/2)
        lower=1. if w0==w1 else np.sqrt(2)-1
        assert low.min()>=lower-1e-13
        spectra.append(dict(w=[w0,w1],analytic_gap='1' if w0==w1 else 'sqrt(2)-1',
            sampled_min=float(low.min()),zero_spectrum_set_empty=True))
    return dict(canonical_original_matrix_cases=rows,spectrum_crosschecks=spectra,
        universal_gap_for_this_family='sqrt(2)-1',proof_in_note_not_finite_sampling=True,
        quotient_normalization=probe.entry.weyl_normalization())


def run():
    algebra=algebra_checks()
    diagnostic=numerical.run()
    assert diagnostic==json.loads(numerical.TARGET.read_text('utf8'))
    for p,h in diagnostic['dependency_hashes'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    assert not diagnostic['fixed_rational_test']['negative_deterministic_certificate']
    deps=('research_note_673.md','research_note_674.md','research_note_675.md','research_note_696.md',
        'joint_gauss_boundary_functional.py','joint_gauss_support_marginal_results.json',
        'round697_drafts/center_holonomy_entry.py','round697_drafts/full_average_probe.py',
        'round697_drafts/importance_full_average.py','round697_drafts/importance_full_average_results.json')
    return dict(date='2026-10-02',round=697,tests_run=2,failures=0,errors=0,
        algebra=algebra,full_measure_numerical_diagnostic=diagnostic,
        scope=dict(analytic_uniform_gap_and_exact_original_source_reduction=True,
            original_full_group_and16channels_retained=True,two_S9_measures_retained=True,
            deterministic_integral_sign_proved=False,physical_RP_undecided=True,
            original_Hb_not_replaced=True,original_HF_identity_not_proved=True,
            distinct_spatial_site_propagation_not_in_this_subcase=True,
            old_space_interfaces_inherited_no_new_space_gap=True),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=697,tests=2,strict_integral_sign=False,all_checks_passed=True)))
