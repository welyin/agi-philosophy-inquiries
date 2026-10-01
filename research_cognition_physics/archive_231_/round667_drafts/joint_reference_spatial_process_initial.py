"""667: actual original reference dictionary with spatial hopping and masses.

All32 modes per node enter each conditional cycle. A six-mode invariant
neutral sector independently checks the BdG thermal formula against Fock.
No whole bosonic/Gauss integral or infinite-volume dynamics is computed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_spinor_subgroup_mass as dictionary
import joint_fermion_gauss_completion as matter
import joint_gauss_fermion_influence as car
import joint_connection_matter_matching as contact

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_reference_spatial_process_results.json'
LAM=-contact.current_matrices()[0]


def norm(x):return float(np.linalg.norm(x,'fro'))


def cycle(m,theta=0.,neutral=False,stagger_phi=False):
    n=32*m;h=np.zeros((n,n),complex);d=np.zeros_like(h)
    pred_h=np.zeros_like(h);pred_d=np.zeros_like(h);naive_h=np.zeros_like(h)
    ph=dictionary.ph_matrix();u=np.kron(np.eye(m),ph[:32,:32]);v=np.kron(np.eye(m),ph[:32,32:])
    w=np.block([[u,v],[v.conj(),u.conj()]])
    linkleft=[];massnorm=0.
    for x in range(m):
        phi=car.PHI[x].copy()
        if neutral:phi[:4]=0
        if stagger_phi:phi[:4]*=(-1)**x
        hx,dx=matter.mass_matrices(phi);ix=slice(32*x,32*x+32)
        h[ix,ix]=hx;d[ix,ix]=dx
        b=ph@car.bdg(hx,dx)@ph.conj().T
        pred_h[ix,ix]=b[:32,:32];pred_d[ix,ix]=b[:32,32:]
        naive_h[ix,ix]=b[:32,:32]
        massnorm+=norm(b[:32,32:])
    hopping=np.zeros_like(h);correct_hop=np.zeros_like(h);naive_hop=np.zeros_like(h)
    for x in range(m):
        y=(x+1)%m;scale=x+1
        c=matter.gauge.group_exp(scale*np.array([.07,-.04,.03,.02,.01,-.03,.06,.05]),3)
        weak=matter.gauge.group_exp(scale*np.array([.11,-.08,.04]),2);z=np.exp(.09j*scale)
        old=matter.representation(c,weak,z)
        left=np.kron(dictionary.left_rep(c,weak,z),np.eye(2));linkleft.append(left)
        t=.23*np.exp(-2*theta);ix=slice(32*x,32*x+32);iy=slice(32*y,32*y+32)
        for out,mat in ((hopping,old),(correct_hop,LAM@left),(naive_hop,left)):
            out[ix,iy]=t*mat;out[iy,ix]=t*mat.conj().T
    h+=hopping;pred_h+=correct_hop;naive_h+=naive_hop
    return dict(h=h,d=d,w=w,correct_h=pred_h,correct_d=pred_d,naive_h=naive_h,
                source_old=car.bdg(-2*hopping,np.zeros_like(h)),
                source_correct=car.bdg(-2*correct_hop,np.zeros_like(h)),
                source_naive=car.bdg(-2*naive_hop,np.zeros_like(h)),
                links=linkleft,massnorm=massnorm)


def thermal(h,d,source,beta=.73):
    b=car.bdg(h,d);ev,v=np.linalg.eigh(b)
    logz=.5*np.logaddexp(beta*ev/2,-beta*ev/2).sum()-.5*beta*np.trace(h).real
    derivative=beta/4*np.sum(np.tanh(beta*ev/2)*np.diag(v.conj().T@source@v).real)
    derivative-=beta/2*np.trace(source[:len(h),:len(h)]).real
    return float(logz),float(derivative),ev


def spatial_transport_check():
    data=cycle(3);b=car.bdg(data['h'],data['d']);w=data['w']
    correct=car.bdg(data['correct_h'],data['correct_d'])
    transform_error=norm(w@b@w.conj().T-correct)
    source_error=norm(w@data['source_old']@w.conj().T-data['source_correct'])
    old=thermal(data['h'],data['d'],data['source_old'])
    new=thermal(data['correct_h'],data['correct_d'],data['source_correct'])
    wrong=thermal(data['naive_h'],data['correct_d'],data['source_naive'])
    thermal_error=max(abs(old[0]-new[0]),abs(old[1]-new[1]))
    spectrum_error=float(np.max(abs(old[2]-new[2])))
    assert max(transform_error,source_error,thermal_error,spectrum_error)<1e-11
    wrong_gaps=[abs(old[i]-wrong[i]) for i in (0,1)]
    assert min(wrong_gaps)>1e-5
    step=2e-5;logs=[]
    for theta in (-step,step):
        q=cycle(3,theta)
        logs.append(thermal(q['h'],q['d'],q['source_old'])[0])
    source_difference_error=abs((logs[1]-logs[0])/(2*step)-old[1])
    assert source_difference_error<1e-7
    loop=np.eye(2,dtype=complex);naive_loop=loop.copy()
    for link in data['links']:
        loop=loop@(-link[30:32,30:32]);naive_loop=naive_loop@link[30:32,30:32]
    assert norm(loop+np.eye(2))==0 and norm(naive_loop-np.eye(2))==0
    # Independent Fock check on the exact neutral invariant sector X=0.
    q=cycle(3,neutral=True);ids=[32*x+i for x in range(3) for i in (30,31)]
    ix=np.ix_(ids,ids);h=q['h'][ix];d=q['d'][ix]
    hb=car.fock(h,d);ev=np.linalg.eigvalsh(hb);beta=.73
    logf=float(np.logaddexp.reduce(-beta*ev))
    thermal_small=thermal(h,d,car.bdg(np.zeros_like(h),np.zeros_like(h)))
    fock_error=abs(logf-thermal_small[0]);assert fock_error<1e-12
    return dict(nodes=3,original_modes=96,whole_mass_and_nontrivial_SM_links=True,
                dictionary_error=transform_error,source_dictionary_error=source_error,
                BdG_spectrum_error=spectrum_error,thermal_transport_error=thermal_error,
                original_logZ=old[0],transported_logZ=new[0],naive_logZ=wrong[0],
                original_hopping_geometry_source=old[1],naive_source=wrong[1],
                source_finite_difference_error=source_difference_error,
                right_neutrino_loop_sign=-1,naive_neutrino_loop_sign=1,
                independent_neutral_Fock_modes=6,neutral_Fock_logZ_error=fock_error,
                conditional_quadratic_branch_only=True)


def joint_staggering_check():
    data=cycle(4);n=128;sign=np.ones(n)
    for x in range(4):
        sign[32*x:32*x+32][np.diag(LAM).real<0]=(-1)**x
    s=np.diag(sign);sn=np.block([[s,np.zeros_like(s)],[np.zeros_like(s),s]])
    corrected=sn@car.bdg(data['correct_h'],data['correct_d'])@sn
    staggered=cycle(4,stagger_phi=True)
    joint=car.bdg(staggered['naive_h'],staggered['correct_d'])
    unchanged_mass=car.bdg(data['naive_h'],data['correct_d'])
    joint_error=norm(corrected-joint);mass_mismatch=norm(corrected-unchanged_mass)
    assert joint_error<1e-12 and mass_mismatch>.1
    before=thermal(data['correct_h'],data['correct_d'],data['source_correct'])
    after=thermal(staggered['naive_h'],staggered['correct_d'],staggered['source_naive'])
    wrong=thermal(data['naive_h'],data['correct_d'],data['source_naive'])
    thermal_error=max(abs(before[i]-after[i]) for i in (0,1))
    assert thermal_error<1e-11 and abs(before[0]-wrong[0])>1e-5
    # Reinterpretation as a scalar-coordinate change must transform its edge action.
    x=car.PHI[0,:2]+1j*car.PHI[0,2:4];y=car.PHI[1,:2]+1j*car.PHI[1,2:4]
    edge_before=float(np.vdot(x-y,x-y).real)
    edge_reset=float(np.vdot(x+y,x+y).real)
    assert abs(edge_before-edge_reset)>.01
    return dict(nodes=4,original_modes=128,bipartite_signs=[1,-1,1,-1],
                hopping_and_mass_joint_error=joint_error,
                omitted_mass_staggering_BdG_error=mass_mismatch,
                joint_thermal_and_source_error=thermal_error,
                original_logZ=before[0],jointly_transformed_logZ=after[0],
                hopping_only_reset_logZ=wrong[0],
                original_identity_link_Higgs_edge_norm=edge_before,
                Higgs_sign_reset_edge_norm=edge_reset,
                scalar_edge_is_unnormalized_common_positive_metric_factor=False,
                Higgs_edge_witness_is_flat_numerator_not_full_curved_graph_energy=True,
                no_claim_of_bipartite_necessity_for_valid_quantum_model=True)


def run():
    deps=('joint_spinor_subgroup_mass.py','joint_fermion_gauss_completion.py',
          'joint_gauss_fermion_influence.py','joint_connection_matter_matching.py',
          'research_note_598.md','research_note_614.md','research_note_662.md','research_note_666.md',
          'round667_drafts/reference_scale_probe.py','round667_drafts/reference_scale_probe_results.json')
    return dict(date='2026-10-02',round=667,tests_run=2,failures=0,errors=0,
                spatial_reference_transport=spatial_transport_check(),
                mass_and_spatial_staggering=joint_staggering_check(),
                dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
                scope='Actual reference automorphism preserves local region/gauge structure and conditionally transports any established state/process limit without requiring an implementer in empty infinite Fock. Original hopping gains right-handed signs; odd loops and simultaneous Yukawa staggering constrain naive rewrites. Conditional96/128-mode quadratic numerics, not full bosonic/Gauss thermal integral, complete infinite dynamics, continuum chiral measure or quantum GR.',
                all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
