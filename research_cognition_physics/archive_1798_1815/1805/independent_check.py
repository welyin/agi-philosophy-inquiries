"""Conditioning of finite family-weight measurements on the frozen charged spectrum."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent


def basis():
    hs=[np.diag([1.,-1,0])/math.sqrt(2),np.diag([1.,1,-2])/math.sqrt(6)]
    for j in range(3):
        for k in range(j+1,3):
            a=np.zeros((3,3),complex);a[j,k]=1/math.sqrt(2)
            hs.extend([a+a.conj().T,1j*(a-a.conj().T)])
    return np.array(hs)


def recover(v,s):
    return (v*s)@v.conj().T


def design(hs, maps):
    def col(h):
        out=np.concatenate([fn(h).ravel() for fn in maps])
        return np.r_[out.real,out.imag]
    return np.column_stack([col(h) for h in hs])


def pinch(v):
    return lambda h:h-recover(v,np.diag(v.conj().T@h@v))


def comm(a):
    return lambda h:a@h-h@a


def conditioning(vv, yy):
    eye=np.eye(3);hs=basis()
    operators=[recover(v,y**2) for v,y in zip(vv,yy)]
    target=operators[0]/np.trace(operators[0])
    target_vec=np.real(np.einsum('aij,ji->a',hs,target))
    report={}
    for label,maps in [
            ('mass_projectors',[pinch(v) for v in vv]),
            ('dimensionless_raw_Y_Ydagger',[comm(a) for a in operators]),
            ('each_Y_Ydagger_divided_by_heaviest_eigenvalue',
             [comm(a/np.linalg.norm(a,2)) for a in operators])]:
        A=design(hs,maps)
        _,s,vt=np.linalg.svd(A,full_matrices=False)
        gain=np.linalg.norm((vt@target_vec)/s)
        direction=(vt.T@((vt@target_vec)/s**2))/gain**2
        dw=np.einsum('a,aij->ij',direction,hs)
        smallest=float(np.min(np.linalg.eigvalsh(eye+.05*dw)))
        assert smallest>0
        assert abs(np.trace(dw@target).real-1)<1e-6
        assert np.isclose(np.linalg.norm(A@direction),1/gain,rtol=1e-6)
        report[label]=dict(min_singular_value=float(s[-1]),max_trace_gain=float(gain),
                           min_eps_for_5pct=float(.05/gain),
                           min_eigenvalue_of_5pct_witness=smallest,
                           min_eigenvalue_of_minus5pct_witness=float(np.min(np.linalg.eigvalsh(eye-.05*dw))))
    principal=np.linalg.svd(abs(vv[0].conj().T@vv[1])**2,compute_uv=False)
    np.testing.assert_allclose(report['mass_projectors']['min_singular_value']**2,
                               1-principal[1],rtol=1e-10)
    return dict(transition_singular_values=principal.tolist(),certificates=report)


def top_count(high_v,low_v,y):
    high_b=high_v[:,-1]
    target=np.diag(y**2/sum(y**2))
    delta=float(.05/(1.5*(high_b.conj()@target@high_b).real-.5))
    w=[1-delta/2,1+delta]
    ordinary=abs(low_v[2,2])**2
    weighted=w[1]*ordinary/(w[1]*ordinary+w[0]*(1-ordinary))
    # Exact product Bernoulli Bhattacharyya affinity, with numerical stability.
    deficit=((math.sqrt(ordinary)-math.sqrt(weighted))**2
             +(math.sqrt(1-ordinary)-math.sqrt(1-weighted))**2)/2
    lower=math.ceil(math.log(.19)/(2*math.log1p(-deficit)))
    critical=.955
    ratio=critical*(1-ordinary)/((1-critical)*ordinary)
    max_delta=2.
    min_delta=(ratio-1)/(1+ratio/2)
    # For this extra rate hypothesis also hold the c-row CKM from its independent
    # calibration and compare the b to non-b odds in the two rows.
    charm=abs(low_v[1,2])**2
    charm_new=w[1]*charm/(w[1]*charm+w[0]*(1-charm))
    return dict(high_delta_bottom_vs_average=delta,weights_light_and_bottom=w,
                low_R_SM=ordinary,low_R_modified=weighted,
                ideal_min_top_decays_for_two_sided_5pct_errors=lower,
                admissible_delta_from_only_CMS_R_under_response=[min_delta,max_delta],
                charm_row_b_fraction_if_common_selector=charm_new,
                charm_row_b_fraction_reference=charm,
                equivalent_on_shell_Wcb_tag_amplitudes=[math.sqrt(charm),math.sqrt(charm_new)],
                ATLAS_onshell_amplitude_one_width=[.036,.061],
                bottom_to_nonbottom_odds_multiplier=w[1]/w[0])


def compute():
    old=runpy.run_path(str(HERE.parent/'1799/check.py'))
    data=old['backend']()
    up,down,_,v,_=data
    stage=runpy.run_path(str(HERE.parent/'1715/check.py'))
    flow,rep=stage['charged_background'](.1)
    low=flow['low_ckm'](rep['gamma_MZ_deg'],flow['angles']('direct'))
    scale=[]
    for s12,s13,s23,delta in [(0.309,.5*(1-math.sqrt(1-.0855)),.55,deg) for deg in [-60.,-70.,-80.]]:
        # These are actual retained low three-active angles, not high-scale quark labels.
        # flow.low_ckm expects the input angle "gamma" of its CKM construction;
        # use the PDG mixing phase delta explicitly instead.
        a,b,c=map(math.sqrt,[s12,s13,s23])
        A,B,C=map(math.sqrt,[1-s12,1-s13,1-s23])
        z=np.exp(1j*math.radians(delta))
        pmns=np.array([[A*B,a*B,b/z],[-a*C-A*c*b*z,A*C-a*c*b*z,c*B],
                       [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])
        np.testing.assert_allclose(pmns@pmns.conj().T,np.eye(3),atol=5e-15)
        sigma=np.linalg.svd(abs(pmns)**2,compute_uv=False)[1]
        scale.append(dict(delta_deg=delta,projector_gap=math.sqrt(1-sigma)))
    return dict(quark=conditioning([np.eye(3),v],[up,down]),top=top_count(v,low,up),
                three_active_lepton=scale)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    raw=json.dumps(compute(),indent=2,allow_nan=False)+'\n'
    f=HERE/'independent_results.json'
    if args.check:
        assert json.loads(raw)==json.loads(f.read_bytes())
    else:
        f.write_bytes(raw.encode())
    print('1805 independent',hashlib.sha256(raw.encode()).hexdigest()[:12])
