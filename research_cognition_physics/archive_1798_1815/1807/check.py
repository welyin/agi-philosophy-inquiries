"""Physical channel incidence, QCD number support, and a held-out gauge coupling."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def load(path):
    return json.loads((HERE/path).read_bytes())


def branch_projection(p0,p1):
    # Sharp after removing the unknown common daughter/source log-weight difference.
    return 2*math.sqrt(p0*p1/(p0+p1))


def common_without_thresholds(g1,g2,mz):
    # SU(5)-normalized U(1); alpha_i^-1(mu) = alpha_i^-1(MZ) - b_i*t/(2*pi).
    b1,b2,b3=41/10,-19/6,-7
    one,two=4*math.pi/g1**2,4*math.pi/g2**2
    ulog=2*math.pi*(one-two)/(b1-b2)
    unif=one-b1*ulog/(2*math.pi)
    a3=unif+b3*ulog/(2*math.pi)
    return dict(mu_GeV=mz*math.exp(ulog),alpha3_at_MZ=1/a3,
                inverse_alpha_common=unif)


def compute():
    cfg=load('inputs.json')
    # Columns: common lepton, physical charged pion, physical proton.
    # Rows: tau->pi nu, pi->mu nu, e+e-->p pbar; Belle mode is redundant given all.
    rows=np.array([[0,1,0],[2,-1,0],[-2,0,2],[2,1,0]])
    assert round(np.linalg.det(rows[:3]))==-4
    assert np.array_equal(rows[3]-rows[0],[2,0,0])
    p0,*s0=cfg['aleph_tau_to_pi_nu']
    p1,*s1=cfg['belle_tau_to_pi_nu_ee_M_pi_ee_above_1p05_GeV']
    working=branch_projection(p0-2*math.hypot(*s0),p1-2*math.hypot(*s1))
    central=branch_projection(p0,p1)
    # Finite color singlets with arbitrary illustrative sea, including a non-valence proton.
    n=np.array([3,5,7,11,17])
    distribution=np.array([.10,.20,.15,.25,.30])
    moments=[]
    for weight in [.8,1.,1.2]:
        wp=distribution@np.power(weight,n)
        assert abs(math.log(wp))+2e-13>=3*abs(math.log(weight))
        moments.append(float(wp))
    ew=load(cfg['electroweak_and_strong'])
    gs=[ew[key] for key in ['g1_GUT','g2','g3']]
    cen=common_without_thresholds(gs[0][0],gs[1][0],ew['mu_ref_GeV'])
    g12corners=[common_without_thresholds(gs[0][0]+i*gs[0][1],gs[1][0]+j*gs[1][1],ew['mu_ref_GeV'])
                for i,j in itertools.product([-2,2],repeat=2)]
    pred=[f['alpha3_at_MZ'] for f in g12corners]
    actual=[(gs[2][0]+i*gs[2][1])**2/(4*math.pi) for i in [-2,2]]
    assert max(pred)<min(actual)
    old=load(cfg['original_common_weights'])
    omega=old['one_loop_corner_bounds']['wl_wq']
    # Edges measured to absolute log discrepancy <=eta. With the pion dictionary:
    # h,q have equal signs and |h|>=2|q|; |h|,|2l-h|<=eta implies |l-q|<=eta.
    # Optimise the same cone with n<=N cutoff: relation h=n_eff*q, n_eff in [2,N].
    ex=[]
    for nmax in [2,6,20]:
        ex.append(dict(nmax=nmax,sharp_eta_for_omega_low=math.log(omega[0])/(1-1/nmax)))
    split=load(cfg['split_thresholds'])['shared_witness']
    return dict(channel_rank=int(np.linalg.matrix_rank(rows)),
                central_pair_log_error_gain=central,two_width_lower_pair_log_error_gain=working,
                closed_state_error_needed_to_detect_5_percent_weight=central*math.log(1.05),
                example_qcd_quark_count=n.tolist(),example_probabilities=distribution.tolist(),
                proton_weight_at_wq_0p8_1_1p2=moments,
                log_edge_lower_for_1369_one_width=math.log(omega[0]),finite_pion_cutoff_bounds=ex,
                equal_weight_one_loop=cen,alpha3_pred_two_ew_marginal_box=[min(pred),max(pred)],
                alpha3_data_two_marginal_box=actual,
                gauge_matched_split_SUSY=dict(U_GeV=split['U_GeV'],alpha_inv=split['inverse_alpha_U'],
                    soft_GeV=split['soft_GeV'],equal_residual=split['all_gauge_equal_residual']))


def stable(x):
    if isinstance(x,float):
        return float(format(x,'.11g'))
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    raw=json.dumps(stable(compute()),indent=2,allow_nan=False)+'\n'
    f=HERE/'results.json'
    if args.check:
        assert json.loads(raw)==json.loads(f.read_bytes())
    else:
        f.write_bytes(raw.encode())
    print('1807',hashlib.sha256(raw.encode()).hexdigest()[:12])
