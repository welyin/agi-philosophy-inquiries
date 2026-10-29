"""1331: exact invariant support, finite orientation source and actual local records.

Uses the frozen 1330 complete parent builder, with its SHA256 recorded in results.
No parent reference direction is removed. Dynamics restriction follows an explicitly
checked exact conserved Z; the measurement rotation uses the complete parent.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
DEPENDENCY=HERE.parent/'1330'/'check.py'
spec=importlib.util.spec_from_file_location('frozen1330',DEPENDENCY)
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)


def norm(a):
    return float(np.linalg.norm(a,2)) if a.size else 0.


def make_operators(q,Ns):
    a=old.physical(q,Ns);d=len(a['basis']);look=a['lookup']
    jp=np.zeros((d,d),complex);jpb=np.zeros_like(jp)
    edges={k:np.zeros_like(jp) for k in ['D','U','L0','L1']}
    z=np.zeros(d,int)
    for col,(s,h0,h1) in enumerate(a['basis']):
        hs=[h0,h1]
        z[col]=h0+h1+sum(((s>>(4*x+3))&1)+((s>>(4*x+1))&1) for x in range(2))
        for x in range(2):
            b=int(a['bs'][col,x]);hh=hs[x]
            if hh:
                out=hs.copy();out[x]-=1
                jpb[look[(s,*out)],col]+=math.sqrt((b-hh+1)*hh)
            moved=old.hop(s,4*x+2,4*x+3)
            if moved:
                ss,sign=moved;jp[look[(ss,h0,h1)],col]+=sign
        for sec,mode,w in [('D',0,0),('U',1,2),('L0',2,1),('L1',3,1)]:
            moved=old.hop(s,4+mode,mode)
            if moved and a['bs'][col,1]-h1>=w:
                ss,sign=moved;edges[sec][look[(ss,h0,h1)],col]+=sign
    jp+=jpb
    a.update(Jx=(jp+jp.conj().T)/2,Jy=(jp-jp.conj().T)/(2j),
             Jxb=(jpb+jpb.conj().T)/2,Z=z,Jz=sum(Ns)/2-z,split_edges=edges)
    return a


def assemble(a,md,mu,coeff):
    d=len(a['basis']);v=np.zeros((d,d),complex);he=np.zeros_like(v)
    for x in range(2):
        t=md[x]*a['D'][x]+mu[x]*a['U'][x];v+=t+t.conj().T
    for s in coeff:
        t=coeff[s]*a['split_edges'][s];he+=t+t.conj().T
    # The common working H0 is zero after the registered independent energy dictionary.
    return v+he,v,he


def source_columns(h,j,idx):
    return 1j*(h@j[:,idx]-j@h[:,idx])


def exp_columns(eig,columns,t):
    w,v=eig
    return (v*np.exp(-1j*t*w))@(v.conj().T@columns)


def source_bound(q,Ns,c):
    m=min(Ns)-3*q
    return 2*math.sqrt(q+1)/math.sqrt(m-2)*(abs(c['L0'])+abs(c['L1'])+2*abs(c['U']))


def records(a,states,initial,h,h0,rotation,t,epsilon=.025,alpha=.01):
    rotated=rotation@states
    # Complete local reference and matter occupations, including D/U locations, label every row.
    probabilities=np.abs(rotated[:,0])**2
    z=a['Jz'];mean=float(z@probabilities)
    variance=float((z*z)@probabilities-mean*mean)
    gram=norm(rotated.conj().T@rotated-np.eye(rotated.shape[1]))
    work=-states.conj().T@(h-h0)@states
    final_energy=rotated.conj().T@h0@rotated
    energy_res=norm(final_energy-initial.conj().T@h@initial-work)
    jxmean=float(np.real(states[:,0].conj()@a['Jx']@states[:,0]))
    interval=float(np.max(z)-np.min(z))
    shots=math.ceil(interval*interval/(2*epsilon*epsilon)*math.log(2/alpha))
    pre_z=float(np.real(states[:,0].conj()@(z*states[:,0])))
    post_x=float(np.real(rotated[:,0].conj()@a['Jx']@rotated[:,0]))
    post_y=float(np.real(rotated[:,0].conj()@a['Jy']@rotated[:,0]))
    local_outcomes=[len({((s>>(4*x))&15,int(a['bs'][i,x]-hs[x]),int(hs[x]))
                        for i,(s,*hs) in enumerate(a['basis'])}) for x in range(2)]
    return dict(time=t,mean=mean,variance=variance,readout_identity_error=abs(mean-jxmean),
        complete_record_probabilities=probabilities.tolist(),probability_sum=float(probabilities.sum()),
        record_gram_error=gram,complete_switch_energy_error=energy_res,
        joint_nonzero_record_result_count=len(z),local_blank_pointer_dimensions=local_outcomes,
        local_record_capacity_bits=sum(math.ceil(math.log2(k)) for k in local_outcomes),
        free_record_compression_assumed=False,readout_range=interval,mean_precision=epsilon,
        failure_probability=alpha,hoeffding_sufficient_repetitions=shots,
        actual_reference_quanta_per_preparation=float(sum(a['bs'][a['lookup'][(1<<2,0,0)]])),
        control_rotation_angle=math.pi/2,control_orientation_transfer_Jz=mean-pre_z,
        control_orientation_transfer_Jx=post_x-jxmean,
        record_orientation_transfer_Jx=-post_x,record_orientation_transfer_Jy=-post_y,
        record_orientation_transfer_Jz=0.,record_pointer_and_control_are_registered_inputs=True,
        registered_angle_accuracy=.001,angle_systematic_bias_bound=interval*.001/2,
        autonomous_controller_derived=False,reference_or_record_reset_hidden=False,
        passed=bool(max(abs(mean-jxmean),gram,energy_res,abs(probabilities.sum()-1))<2e-9))


def run_case(q,Ns,do_records):
    a=make_operators(q,Ns);idx=np.flatnonzero(a['Z']<=q);off=np.flatnonzero(a['Z']>q)
    md=[.012,.017];mu=[.009,.013]
    c={'D':.31,'U':.43,'L0':.9,'L1':.2}
    csym=dict(c,L0=.55,L1=.55)
    h,v,he=assemble(a,md,mu,c);hs,vs,hes=assemble(a,md,mu,csym)
    hp=h[np.ix_(idx,idx)];vp=v[np.ix_(idx,idx)];ep=he[np.ix_(idx,idx)]
    jsp=a['Jx'][np.ix_(idx,idx)]
    cpcols=source_columns(h,a['Jx'],idx);cp=cpcols[idx]
    cref=source_columns(he,a['Jxb'],idx)
    symcols=source_columns(hs,a['Jx'],idx)
    symref=source_columns(hes,a['Jxb'],idx)
    z_comm=float(np.linalg.norm(h*(a['Z'][None,:]-a['Z'][:,None])))
    inv=norm(h[np.ix_(off,idx)])
    v_comm=norm(source_columns(v,a['Jx'],idx))
    refnorm=norm(cref);symnorm=norm(symcols);b=source_bound(q,Ns,c);bs=source_bound(q,Ns,csym)
    fermion_cancel=norm(symcols-symref)
    minb0=min(int(a['bs'][i,x]-a['basis'][i][1+x]) for i in idx for x in range(2))
    maxh=max(a['basis'][i][1]+a['basis'][i][2] for i in idx)
    identity=float(np.max(np.abs(a['Z']-(sum(Ns)/2-a['Jz']))))
    # All actual unknown initial left columns, including arbitrary partners.
    distribution=[1,0] if q==1 else [1,1]
    cols=[a['lookup'][(s,0,0)] for s in old.masks(q)
          if all(old.counts(s)[x,0]==old.counts(s)[x,1]==0 and
                 old.counts(s)[x,2]==distribution[x] for x in range(2))]
    initial=np.eye(len(h),dtype=complex)[:,cols]
    # First column is the physically prepared equal-direction probe; others complete an isometry.
    if q==1:
        initial=initial@np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
    else:
        transform=np.array([[1,1],[1,-1]],complex)/math.sqrt(2)
        initial=initial@np.kron(transform,transform)
    initp=initial[idx];eig=np.linalg.eigh(hp);eeig=np.linalg.eigh(ep)
    seig=np.linalg.eigh(hs[np.ix_(idx,idx)])
    hn=norm(hp);cn=norm(cp);vn=norm(vp);ccn=norm(cpcols)
    d1=initp.conj().T@cp@initp
    d2=initp.conj().T@(1j*(hp@cp-cp@hp))@initp
    curvature_error=None
    if q==1:curvature_error=abs(float(np.real(d2[0,0]))+.5*(c['L0']-c['L1'])**2)
    # All-q compressed generator checks are independent of the particular readout probe.
    support_residual=norm(h[:,idx]-np.eye(len(h),dtype=complex)[:,idx]@hp)
    times=[.65,1.1,1.65]
    holdouts=[];rec=[]
    rotation=None;rotation_identity=None
    if do_records:
        w,rr=np.linalg.eigh(a['Jy']);rotation=(rr*np.exp(1j*math.pi*w/2))@rr.conj().T
        rotation_identity=norm(rotation.conj().T@(a['Jz'][:,None]*rotation)-a['Jx'])
    for t in times:
        statep=exp_columns(eig,initp,t);edp=exp_columns(eeig,initp,t)
        symp=exp_columns(seig,initp,t)
        f=float(np.real(statep[:,0].conj()@jsp@statep[:,0]))
        fe=float(np.real(edp[:,0].conj()@jsp@edp[:,0]))
        fs=float(np.real(symp[:,0].conj()@jsp@symp[:,0]))
        f0=float(np.real(initp[:,0].conj()@jsp@initp[:,0]))
        point_bound=t*t*vn*ccn
        nodes,weights=np.polynomial.legendre.leggauss(16)
        integral=0.
        for node,weight in zip(nodes,weights):
            y=exp_columns(eig,initp[:,[0]],t*(node+1)/2)[:,0]
            integral+=weight*float(np.real(y.conj()@cp@y))*t/2
        remainder=None;remainder_bound=2/3*hn*hn*cn*t**3
        pure_formula=None
        if q==1:
            pure_formula=.5*math.cos((c['L0']-c['L1'])*t)
            remainder=abs(f-.5+t*t*(c['L0']-c['L1'])**2/4)
        holdouts.append(dict(time=t,total_Jx_mean=f,equal_left_Jx_mean=fs,pure_edge_mean=fe,
            pure_edge_formula=pure_formula,source_integral=integral,source_integral_error=abs(integral-f+f0),
            point_vertex_effect=abs(f-fe),point_vertex_effect_bound=point_bound,
            finite_curvature_remainder=remainder,finite_curvature_remainder_bound=remainder_bound,
            orientation_drop=f0-f,complete_input_gram_error=norm(statep.conj().T@statep-np.eye(len(cols))),
            passed=bool(abs(f-fe)<=point_bound+1e-9 and abs(integral-f+f0)<2e-9 and
                        (pure_formula is None or abs(fe-pure_formula)<2e-9) and
                        (remainder is None or remainder<=remainder_bound+1e-9))))
        if do_records:
            full=np.zeros_like(initial);full[idx]=statep
            record=records(a,full,initial,h,np.zeros_like(h),rotation,t)
            record['certified_drop_lower_bound']=.5-pure_formula-point_bound
            record['two_mean_total_tolerance']=2*(record['mean_precision']+record['angle_systematic_bias_bound'])
            record['two_mean_failure_probability_bound']=2*record['failure_probability']
            record['two_mean_sufficient_total_preparations']=2*record['hoeffding_sufficient_repetitions']
            record['separation_margin']=record['certified_drop_lower_bound']-record['two_mean_total_tolerance']
            rec.append(record)
    passed=max(z_comm,inv,v_comm,fermion_cancel,identity,support_residual)<2e-9
    passed=passed and refnorm<=b+1e-9 and symnorm<=bs+1e-9 and minb0>=min(Ns)-3*q and maxh<=q
    passed=passed and (curvature_error is None or curvature_error<2e-9)
    passed=passed and (rotation_identity is None or rotation_identity<2e-9)
    passed=passed and all(x['passed'] for x in holdouts+rec)
    if do_records:passed=passed and any(x['separation_margin']>0 for x in rec)
    return dict(q=q,Ns=Ns,parent_dimension=len(h),exact_invariant_support_dimension=len(idx),
        actual_input_dimension=len(cols),Z_commutator_error=z_comm,invariant_support_error=inv,
        point_SU2_error=v_comm,minimum_b0_on_support=minb0,max_b1_on_support=maxh,
        source_ref_norm=refnorm,source_ref_bound=b,equal_left_full_source_norm=symnorm,
        equal_left_source_bound=bs,equal_left_fermion_cancellation_error=fermion_cancel,
        first_source_probe=float(np.real(d1[0,0])),second_derivative_probe=float(np.real(d2[0,0])),
        exact_actual_curvature_error=curvature_error,support_H_norm=hn,support_source_norm=cn,
        point_perturbation_norm=vn,complete_source_columns_norm=ccn,
        rotation_instrument_identity_error=rotation_identity,holdouts=holdouts,records=rec,
        record_basis=[list(x) for x in a['basis']] if do_records else None,
        passed=bool(passed))


def matrix_classification_check():
    # Independent finite-dimensional source/receiver flavour calculation, no physical encoding permission inferred.
    j0=np.array([[.7,.1+.2j],[-.15j,.35]],complex)
    j1=np.array([[.25,.2j],[.08,.6]],complex)
    n=2;delta=j0-j1
    h0=np.block([[np.zeros((n,n)),j0.conj().T],[j0,np.zeros((n,n))]])
    h1=np.block([[np.zeros((n,n)),j1.conj().T],[j1,np.zeros((n,n))]])
    h=np.block([[h0,np.zeros_like(h0)],[np.zeros_like(h1),h1]])
    jx=np.block([[np.zeros_like(h0),np.eye(2*n)/2],[np.eye(2*n)/2,np.zeros_like(h0)]])
    f=np.zeros((4*n,n),complex);f[:n]=np.eye(n)/math.sqrt(2);f[2*n:3*n]=np.eye(n)/math.sqrt(2)
    c=1j*(h@jx-jx@h);d2=f.conj().T@(1j*(h@c-c@h))@f
    error=norm(d2+delta.conj().T@delta/2)
    return dict(noncommuting_flavour_example=True,curvature_matrix_error=error,
                compressed_source_norm=norm(c),expected_source_norm=norm(delta)/2,
                passed=bool(error<1e-12 and abs(norm(c)-norm(delta)/2)<1e-12))


def compute():
    cases=[run_case(1,[6,6],True),run_case(1,[9,9],True),run_case(2,[9,9],False)]
    matrix=matrix_classification_check()
    paths=[HERE/'proof.md',HERE/'check.py',DEPENDENCY,HERE.parent/'1330'/'proof.md',
           HERE.parent/'_admission'/'round1331_admission.md']
    provenance={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    return dict(round=1331,status='PASS' if all(x['passed'] for x in cases) and matrix['passed'] else 'FAIL',
                provenance=provenance,parameters=dict(MD=[.012,.017],MU=[.009,.013],JD=.31,JU=.43,
                JL0=.9,JL1=.2,equal_left=.55,working_p=[0,0]),refits=0,cases=cases,
                arbitrary_matrix_check=matrix,new_axioms=0,empirical_experiments=0,
                complete_parent_directions_kept=True,local_Gauss_claim=False,stage_complete=False,long_goal_complete=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),(path,'keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-8,abs_tol=2e-9),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group()
    g.add_argument('--check',action='store_true');g.add_argument('--write-new',action='store_true')
    args=p.parse_args();r=compute();assert r['status']=='PASS',r
    dst=HERE/'results.json'
    if args.write_new:
        payload=json.dumps(r,ensure_ascii=False,indent=2)+'\n'
        with dst.open('x',encoding='utf-8') as f:f.write(payload)
    else:compare(r,json.loads(dst.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',round=1331,complete_parents=3,time_holdouts=9,
                         complete_records=6,results_written=args.write_new)))


if __name__=='__main__':main()
