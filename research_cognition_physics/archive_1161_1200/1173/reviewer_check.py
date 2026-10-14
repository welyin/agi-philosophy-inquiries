"""1173 independent review. No author-module import; default is read only.

The reviewer participated in admission/scope discussion, so this is not blind.
The numerical comparison uses RK4, not the author's midpoint exponential code.
Small complete site circuits independently check the Fourier and relay claims.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
FROZEN = {
    'archive_1161_/research_note_1173.md': '8400637e1ccd3aeb8c97d47703792941bac1afb5bd8b81bc6539ad81af691742',
    'archive_1161_/1173/admission.md': 'e91343edaeec04e3ff45e8cb8e7cfac1c56a5e739eaf54d619e5b9631181d8a9',
    'archive_1161_/1173/proof.md': '67aaffaa8d19962c34acdd3ad08bc3f3729a3c08a51f3ddae6edc7e37dd79be6',
    'archive_1161_/1173/check.py': 'd18e37efdb6e4167ded60d42b09d55ee1adc36cca6c4548187993bed3999b26a',
    'archive_1161_/1173/results.json': 'a634cd6f2d74ddbb4e86cac1ec561cbb7d13b8095000b74d8c4d2a4c522d9204',
    'archive_1161_/1173/NEXT.md': 'c29ab9d3855820ceaa7ce2121157582db8009221b84f4384069f4f4cf84566a8',
}
I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1., -1.]).astype(complex)


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def bindings(saved):
    assert len(FROZEN) == 6, 'author freeze is not yet registered'
    out = {p: sha(RESEARCH / p) for p in FROZEN}
    assert out == FROZEN, 'author core changed after freeze'
    for p, digest in saved['source_sha256'].items():
        assert sha(RESEARCH / p) == digest, p
    return out


def eig_exp(h, t):
    vals, vecs = np.linalg.eigh(h)
    return (vecs * np.exp(-1j * t * vals)) @ vecs.conj().T


def independent_evolution(steps=16384):
    ns = np.arange(-128, 129)
    ps = ns * (2 * math.pi / 96)
    ms = np.array([3., 3.4])
    u = np.broadcast_to(I, (2, len(ns), 2, 2)).copy()
    pp = ps[None, :, None]
    dt = 8. / steps

    def rhs(t, v):
        bb = (ms * (1 + math.sin(math.pi * t / 8.)**2))[:, None, None]
        out = np.empty_like(v)
        out[:, :, 0] = -1j * (pp * v[:, :, 0] + bb * v[:, :, 1])
        out[:, :, 1] = -1j * (bb * v[:, :, 0] - pp * v[:, :, 1])
        return out

    for k in range(steps):
        t = k * dt
        a = rhs(t, u)
        b = rhs(t + dt / 2, u + dt * a / 2)
        c = rhs(t + dt / 2, u + dt * b / 2)
        d = rhs(t + dt, u + dt * c)
        u += dt * (a + 2*b + 2*c + d) / 6
    flat = np.array([[eig_exp(p*Z + m*X, 8.) for p in ps] for m in ms])
    ph = np.exp(-(ps - 4.)**2 / (4 * .6**2))
    ph /= np.linalg.norm(ph)
    initials = []
    for m in ms:
        # Positive-frequency preparation obtained independently by diagonalization.
        spins = []
        for p in ps:
            _, v = np.linalg.eigh(p*Z + m*X)
            pos = v[:, 1]
            pos *= np.sign(pos[0])
            spins.append(pos)
        initials.append(np.array(spins) * ph[:, None])
    initials = np.array(initials)
    states = np.array([np.einsum('dnij,dnj->dni', flat, initials),
                       np.einsum('dnij,dnj->dni', u, initials)])
    defect = float(np.max(np.linalg.norm(u.conj().swapaxes(-1,-2) @ u-I, axis=(-2,-1))))
    return ns, ps, states, defect


def finite_port(ns, length, spacing, threshold):
    sites = round(length / spacing)
    first = math.floor((threshold + length/2) / spacing) + 1
    count = sites - first
    delta = 2 * math.pi * (ns[None, :] - ns[:, None]) / length
    out = np.empty(delta.shape, complex)
    mask = delta != 0
    phase = np.exp(1j * delta[mask] * (-length/2 + first*spacing))
    out[mask] = phase * np.expm1(1j * delta[mask] * spacing * count) / np.expm1(1j * delta[mask] * spacing) / sites
    out[~mask] = count / sites
    return out


def finite_sum_and_sampling():
    ns = np.arange(-7, 8)
    length, sites = 8., 64
    xs = -length/2 + np.arange(sites) * length/sites
    ff = np.exp(1j * xs[:, None] * 2 * math.pi * ns[None, :] / length) / math.sqrt(sites)
    exact = ff[xs > .73].conj().T @ ff[xs > .73]
    closed = finite_port(ns, length, length/sites, .73)
    err = float(np.linalg.norm(exact - closed, 2))
    isometry = float(np.linalg.norm(ff.conj().T @ ff - np.eye(len(ns)), 2))
    assert err < 2e-14 and isometry < 2e-14
    pmax, a = 8*math.pi/3, 1/65536
    embedding_exact = math.sqrt(2 * (1 - np.sinc(pmax*a/(2*math.pi))))
    assert embedding_exact < a*pmax
    return dict(direct_site_sum_error=err, finite_sampling_isometry_error=isometry,
                actual_band_embedding_operator_norm=embedding_exact,
                used_embedding_upper=a*pmax)


def projected_records(states, ka):
    # Project the coherent amplitude onto rank-one source and clock eigenvectors.
    eight = []
    for sx in (-1, 1):
        vx = np.array([1, sx], complex) / math.sqrt(2)
        for sy in (-1, 1):
            vals, vecs = np.linalg.eigh(-Y)
            vy = vecs[:, np.argmin(abs(vals-sy))]
            amp = np.einsum('q,d,qdnc->nc', vx.conj(), vy.conj(), states) / 2
            yes = float(np.vdot(amp, ka @ amp).real)
            total = float(np.vdot(amp, amp).real)
            eight.extend([dict(source_X=sx, clock_minusY=sy, port=0, probability=total-yes),
                          dict(source_X=sx, clock_minusY=sy, port=1, probability=yes)])
    port = [[float(np.vdot(states[q,d], ka@states[q,d]).real) for d in (0,1)] for q in (0,1)]
    _, vec = np.linalg.eigh(-Y)
    minus_y = vec[:, 1]
    joint = []
    for q in (0,1):
        amp = np.einsum('d,dnc->nc', minus_y.conj(), states[q]) / math.sqrt(2)
        joint.append(float(np.vdot(amp, ka@amp).real))
    gamma = np.array([np.vdot(states[0,d], states[1,d]) for d in (0,1)])
    return eight, port, joint, gamma


def small_complete_qw():
    # Direct site coin/shift on a whole prepared subspace, with no phase alignment.
    nsite, length, steps = 32, 4., 8
    a, duration = length/nsite, steps*length/nsite
    ns = np.arange(-3,4)
    xs = -length/2 + np.arange(nsite)*a
    ps = 2*math.pi*ns/length
    ff = np.exp(1j*xs[:,None]*ps[None,:])/math.sqrt(nsite)
    maxerr, maxbound = 0., 0.
    for q in (0,1):
        for mass in (.7,1.2):
            init = np.zeros((nsite,2,len(ns)*2), complex)
            for i in range(len(ns)):
                init[:,0,2*i] = ff[:,i]
                init[:,1,2*i+1] = ff[:,i]
            direct = init.copy()
            freq = np.broadcast_to(I,(len(ns),2,2)).copy()
            for k in range(steps):
                b = mass*(1+q*math.sin(math.pi*(k+.5)/steps)**2)
                coin = eig_exp(b*X,a)
                direct = np.einsum('cd,jde->jce',coin,direct)
                direct[:,0] = np.roll(direct[:,0],1,axis=0)
                direct[:,1] = np.roll(direct[:,1],-1,axis=0)
                for i,p in enumerate(ps):
                    freq[i] = eig_exp(p*Z,a) @ coin @ freq[i]
            expected = np.empty_like(direct)
            for i in range(len(ns)):
                expected[:,:,2*i:2*i+2] = ff[:,i,None,None]*freq[i][None,:,:]
            maxerr = max(maxerr,float(np.linalg.norm((direct-expected).reshape(nsite*2,-1),2)))
            # A bound check with a much finer numerical comparator on this small set.
            fine = np.broadcast_to(I,(len(ns),2,2)).copy()
            h = duration/4096
            for k in range(4096):
                b=mass*(1+q*math.sin(math.pi*(k+.5)/4096)**2)
                for i,p in enumerate(ps):
                    fine[i]=eig_exp(p*Z+b*X,h)@fine[i]
            actual=max(float(np.linalg.norm(freq[i]-fine[i],2)) for i in range(len(ns)))
            bound=duration*a*max(abs(ps))*mass*(1+q)+q*mass*math.pi*a/4
            assert actual < bound + q*mass*math.pi*h/4
            maxbound=max(maxbound,actual)
    assert maxerr<3e-14
    return dict(full_prepared_site_vs_Fourier_error=maxerr,
                coarse_actual_vs_fine_comparator_error=maxbound,
                full_chirality_and_unknown_internal_coherence_checked=True)


def local_relay():
    # Four actual local label qubits, source, carrier position/chirality/clock and R.
    labels, nsite = 4, 4
    # axes: source, G-bits encoded as integer, internal clock, position, chirality, R
    rng=np.random.default_rng(1173)
    source_dr=rng.normal(size=(2,2,2))+1j*rng.normal(size=(2,2,2))
    source_dr/=np.linalg.norm(source_dr)
    init=np.zeros((2,1<<labels,2,nsite,2,2),complex)
    packet=rng.normal(size=(nsite,2))+1j*rng.normal(size=(nsite,2))
    packet/=np.linalg.norm(packet)
    init[:,0]=np.einsum('qdr,jc->qdjcr',source_dr,packet)

    def cnot(v, control, target):
        out=np.empty_like(v)
        for q in (0,1):
            for g in range(1<<labels):
                bit=q if control<0 else ((g>>control)&1)
                out[q,g^(bit<<target)]=v[q,g]
        return out

    edges=[(-1,0),(0,1),(1,2),(2,3)]
    v=init.copy()
    for control,target in edges: v=cnot(v,control,target)
    prep_err=float(np.sum(abs(v[:,1:-1])**2))
    expected=init[:,0].copy()
    for k in range(3):
        factor=math.sin(math.pi*(k+.5)/3)**2
        for q in (0,1):
            for g in range(1<<labels):
                for d,mass in enumerate((.7,1.2)):
                    for j in range(nsite):
                        coin=eig_exp(mass*(1+((g>>j)&1)*factor)*X,.1)
                        v[q,g,d,j]=coin@v[q,g,d,j]
            for d,mass in enumerate((.7,1.2)):
                coin=eig_exp(mass*(1+q*factor)*X,.1)
                expected[q,d]=np.einsum('cd,jdr->jcr',coin,expected[q,d])
        v[:,:,:,:,0,:]=np.roll(v[:,:,:,:,0,:],1,axis=3)
        v[:,:,:,:,1,:]=np.roll(v[:,:,:,:,1,:],-1,axis=3)
        expected[:,:,:,0,:]=np.roll(expected[:,:,:,0,:],1,axis=2)
        expected[:,:,:,1,:]=np.roll(expected[:,:,:,1,:],-1,axis=2)
    for control,target in reversed(edges): v=cnot(v,control,target)
    cleared=float(np.sum(abs(v[:,1:])**2))
    err=float(np.linalg.norm(v[:,0]-expected))
    assert prep_err < 1e-28 and cleared<1e-28 and err<2e-14
    return dict(reset_G_leakage=cleared, full_source_clock_carrier_R_error=err,
                causal_neighbour_edges=edges, source_distribution_preserved=True,
                arbitrary_state_cloning_claimed=False,
                feedback_means_coherent_source_state_only=True)


def compare_outputs(saved, states, defect, ns):
    ka=finite_port(ns,96.,1/65536,5.8)
    eig=np.linalg.eigvalsh(ka)
    assert min(eig)>-2e-12 and max(eig)<1+2e-12
    eight,port,joint,gamma=projected_records(states,ka)
    old=saved['group2_same_moving_clock_position_and_source']
    target_gamma=np.array(old['source_gamma']['real'])+1j*np.array(old['source_gamma']['imag'])
    gammaerr=float(max(abs(gamma-target_gamma)))
    errs=[gammaerr]
    for a,b in zip(eight,old['all_eight_joint_records_plus_plus']):
        assert {k:a[k] for k in ('source_X','clock_minusY','port')} == {k:b[k] for k in ('source_X','clock_minusY','port')}
        errs.append(abs(a['probability']-b['probability']))
    for q in (0,1):
        errs.append(abs(joint[q]-old['moving_clock_joint_effects'][f'q{q}']['plus_internal_joint_probability']))
        for d in (0,1):
            errs.append(abs(port[q][d]-old['branches'][f'q{q}_d{d}']['comparator_site_port']))
    maxerr=max(errs)
    assert maxerr<3e-7 and defect<5e-9
    # Recompute all analytic constants without using author helpers.
    b=saved['group3_uniform_certificates_and_resources']
    a=1/65536; pmax=8*math.pi/3; m=3.4
    nsites=round(96/a)
    # The source is a separate factor: include its injection edge into G_0.
    assert b['relay_CNOT_count_upper'] == 2*nsites
    assert b['relay_single_trip_serial_time_upper'] == nsites*a
    ee=[8*a*pmax*m*(1+q)+q*m*math.pi*a/4+q*m*math.pi*(8/8192)/4+2e-6 for q in (0,1)]
    assert max(abs(np.array(ee)-b['actual_to_comparator_per_branch']))<1e-15
    assert abs(joint[1]-joint[0])-sum(ee)>.54
    assert abs(gamma[0].imag-gamma[1].imag)-2*sum(ee)>1.19
    assert abs(sum(p['probability'] for p in eight)-1)<5e-9
    assert min(p['probability'] for p in eight)>-1e-12
    return dict(method='16384-step RK4 matrix evolution; independent eigenvector preparation; exponential finite sum; rank-one output contractions',
                author_module_imported=False, maximum_saved_observable_difference=maxerr,
                source_gamma_difference=gammaerr, RK4_unitarity_defect=defect,
                eight_records=eight, branch_port_probabilities=port,
                moving_clock_joint_probabilities=joint,
                source_gamma_real=gamma.real.tolist(),source_gamma_imag=gamma.imag.tolist(),
                absolute_branch_error_certificates=ee)


def calculate():
    saved=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    before=bindings(saved)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    run=subprocess.run([sys.executable,'-B','-X','utf8',str(HERE/'check.py')],capture_output=True,text=True,encoding='utf-8',env=env)
    assert run.returncode==0,(run.returncode,run.stdout,run.stderr)
    author=json.loads(run.stdout)
    assert author['status']=='READ_ONLY_MATCH' and author['scientific_groups']==3
    ns,ps,states,defect=independent_evolution()
    observables=compare_outputs(saved,states,defect,ns)
    sampling=finite_sum_and_sampling()
    direct=small_complete_qw()
    relay=local_relay()
    assert bindings(saved)==before
    return dict(schema='round1173_independent_review_v1',status='PASS_WITH_SCOPE',blockers=0,
                reviewer_script_sha256=sha(Path(__file__)),author_core_sha256=before,
                author_default_run=dict(exit_code=run.returncode,stdout=author,stderr=run.stderr),
                independent_observables=observables,finite_port_and_embedding=sampling,
                direct_site_QW=direct,causal_label_relay=relay,
                independence=dict(blind_review=False,pre_freeze_contribution='admission choice and scope: temporal-profile identification, absolute phases, full moving-clock readout, positive-frequency and complete-resource limits',author_core_edited=False,author_module_imported=False),
                scope=dict(whole_finite_QW_executed=False,actual_large_grid_connection='analytic uniform prepared-subspace bound',all_unknown_Q_D_band_R=True,negative_frequency_discarded=False,one_dimensional_corridor=True,externally_specified_finite_program_permission=True,physical_device_and_energy_momentum_ledger_certified=False,source_backreaction='coherent source X/Y only, not population or geometry field equation',local_clock_and_site_records_before_causal_G_return=True,new_adopted_axioms=0,GR_hard_deliverables_completed=0,constructive_goal_completed=False))


def compare(actual,saved):
    if isinstance(saved,dict):
        assert actual.keys()==saved.keys()
        for k in saved: compare(actual[k],saved[k])
    elif isinstance(saved,list):
        assert len(actual)==len(saved)
        for a,b in zip(actual,saved): compare(a,b)
    elif isinstance(saved,float):
        assert math.isclose(actual,saved,rel_tol=2e-9,abs_tol=3e-10),(actual,saved)
    else: assert actual==saved,(actual,saved)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--write-once',action='store_true')
    args=ap.parse_args()
    result=calculate()
    path=HERE/'review_checks.json'
    if args.write_once:
        with path.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False); f.write('\n')
    else: compare(result,json.loads(path.read_text(encoding='utf-8')))
    print(json.dumps(dict(status=result['status'],mode='WRITE_ONCE' if args.write_once else 'READ_ONLY_MATCH',blockers=0,independent_maximum_difference=result['independent_observables']['maximum_saved_observable_difference']),ensure_ascii=False))
