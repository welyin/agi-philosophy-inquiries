"""661: original physical Weyl/auxiliary state and a common source.

The tensor/Schur positivity fact is inherited mathematics. New checks map the
actual retained Weyl algebra, original auxiliary weight and source insertions.
No original scalar record instrument or common continuum H is constructed.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_local_mirror_process as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_physical_auxiliary_state_results.json'


def weyl(q):
    sites=q['sites'];count=len(sites);nt=max(t for t,x in sites)+1
    jm=np.kron(np.kron(np.eye(count),old.old.old.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(count),old.old.old.VP),np.eye(16))
    w=jm.conj().T@q['v'];s=w@np.linalg.inv(q['Kl']);r=len(s)
    g=np.block([[np.zeros_like(s),-s],[s.T,np.zeros_like(s)]])
    perm=[sites.index((nt-1-t,x)) for t,x in sites]
    reflect=np.eye(count)[perm]
    unit=old.old.old.VP.conj().T@old.old.old.spin.GAMMA[3]@old.old.old.VM
    a=np.kron(np.kron(reflect,unit.T),np.eye(16))
    b=np.kron(np.kron(reflect,unit),np.eye(16))
    zero=np.zeros_like(a);theta=np.block([[zero,a],[b,zero]])
    pos=np.repeat([t>=nt//2 for t,x in sites],32);pos=np.concatenate((pos,pos))
    idx=np.flatnonzero(pos)
    gram=(theta@g)[np.ix_(idx,idx)]
    assert old.err(gram-gram.conj().T)<2e-12
    assert np.min(np.linalg.eigvalsh(gram))>-3e-12
    select=np.zeros((2*r,4*r),complex)
    select[:r,r:2*r]=w;select[r:,3*r:]=np.eye(r)
    source_map=select@q['T']@q['S'].conj().T
    return dict(G=g,Theta=theta,positive=idx,Gram=gram,map=source_map)


def rising(a,n):
    value=Fraction(1)
    for i in range(n):value*=a+i
    return value


def jacobi_s9(order=10):
    # Normalized weight (1-t^2)^(7/2); Golub-Welsch, not sampled directions.
    k=np.arange(1,order,dtype=float)
    off=np.sqrt(k*(k+7)/((2*k+7)**2-1))
    matrix=np.diag(off,1)+np.diag(off,-1)
    nodes,vec=np.linalg.eigh(matrix)
    return nodes,vec[0]**2


def original_joint_gram_check():
    q=old.data(nx=1,nt=2);w=weyl(q)
    eig=np.linalg.eigvalsh(w['Gram'])
    # Whole S9xS9 integral, using original 2-time Pfaffian, not a plane measure.
    nodes,weights=jacobi_s9();z=0j;first=0j;maxweight=0.
    f=old.old.frame(1,1,2)
    for t,weight in zip(nodes,weights):
        e=np.zeros((2,10));e[0,0]=1;e[1,0]=t;e[1,1]=np.sqrt(1-t*t)
        actual=old.old.ratio(f,e);predicted=((1+t)/2)**16
        maxweight=max(maxweight,float(abs(actual-predicted)))
        z+=weight*actual;first+=weight*t*actual
    exact=rising(Fraction(9,2),16)/rising(Fraction(9),16)
    mean=first/z
    assert abs(z-float(exact))<2e-14 and abs(mean-Fraction(16,25))<2e-12
    # Features (1,sqrt(10) E_a): exact auxiliary reflected Gram.
    aux=np.diag([1.]+[mean.real]*10)
    joint=np.kron(aux,w['Gram'])
    mineig=float(np.min(np.linalg.eigvalsh(joint)))
    assert mineig>-3e-12
    # A genuinely spatial sample kernel and a physical source insertion.
    q=old.data(nx=2,nt=4);w=weyl(q);idx=w['positive']
    rng=np.random.default_rng(66101);hs=4
    fields=np.eye(10)[0]+.13*rng.normal(size=(4,hs,10))
    fields/=np.linalg.norm(fields,axis=2)[:,:,None]
    kernel=np.zeros((4,4),complex)
    frame=old.old.frame(1,2,4)
    for i in range(4):
        for j in range(4):
            negative=fields[i].reshape(2,2,10)[::-1].reshape(4,10)
            e=np.concatenate((negative,fields[j]))
            kernel[i,j]=old.old.ratio(frame,e)
    assert old.err(kernel-kernel.conj().T)<2e-12
    assert min(np.linalg.eigvalsh(kernel))>-2e-12
    # Deliberately nonorthogonal wave packets, all on the original positive half.
    p=rng.normal(size=(len(idx),4))+1j*rng.normal(size=(len(idx),4))
    p/=np.linalg.norm(p,axis=0)
    wg=p.conj().T@w['Gram']@p
    joint_sample=kernel*wg
    assert min(np.linalg.eigvalsh(joint_sample))>-2e-12
    # Direct Pfaffian insertion checks an off-diagonal kernel element, with
    # source rows transported exactly, retaining its complex phase.
    i,j=0,1;negative=fields[i].reshape(2,2,10)[::-1].reshape(4,10)
    qq=old.data(nx=2,nt=4,e=np.concatenate((negative,fields[j])))
    ww=weyl(qq)
    left=p[:,i].conj()@ww['Theta'][idx,:]@ww['map']
    right=p[:,j]@ww['map'][idx,:]
    source=np.outer(left,right)-np.outer(right,left)
    step=.001
    base=old.old.old.pfaffian(qq['N'])
    inserted=(old.old.old.pfaffian(qq['N']+step*source)-base)/(step*base)
    insertion_error=float(abs(inserted-wg[i,j]))
    assert insertion_error<3e-8
    return dict(single_spatial_site_full_S9_partition=str(exact),
        partition_error=float(abs(z-float(exact))),original_weight_error=maxweight,
        reflected_vector_feature_eigenvalue=old.cpair(mean),expected_feature_eigenvalue='16/25',
        Weyl_Gram_eigenvalue_range=[float(min(eig)),float(max(eig))],
        fully_integrated_joint_Gram_dimension=len(joint),joint_minimum_eigenvalue=mineig,
        spatial_sample_auxiliary_minimum=float(min(np.linalg.eigvalsh(kernel))),
        spatial_sample_joint_minimum=float(min(np.linalg.eigvalsh(joint_sample))),
        actual_source_Pfaffian_insertion=old.cpair(inserted),
        target_Weyl_insertion=old.cpair(wg[i,j]),insertion_error=insertion_error,
        spatial_sample_not_claimed_as_complete_spatial_sphere_integration=True,
        all_volume_positivity_depends_on_analytic_original_object_mapping=True)


def common_source_check():
    # Original615 charged flat holonomy; source algebra only, not gauge RP.
    nodes,weights=np.polynomial.legendre.leggauss(7)
    radial=(nodes+1)/2;weights=weights/2*12*radial**2*(1-radial)
    theta=.23;h=1e-4;angles=(theta-h,theta,theta+h)
    charge_index=int(np.flatnonzero(old.old.old.Q==-4)[0]);charge=-4
    totals=[];means=[];max_cov=0.
    for angle in angles:
        unnorm=0j;z=0j
        for rho,weight in zip(radial,weights):
            e=np.zeros(10);e[0]=np.sqrt(rho);e[6]=np.sqrt(1-rho)
            q=old.data('holonomy',angle,e=e);r=q['r']
            jm=np.kron(np.eye(16),old.old.old.VM)
            # Fixed local spin component, expressed through the actual chiral v.
            sel=np.zeros((2*r,4*r),complex)
            sel[:r,r:2*r]=jm.conj().T@q['v'];sel[r:,3*r:]=np.eye(r)
            transform=sel@q['T']@q['S'].conj().T
            g=transform@(-np.linalg.inv(q['N']))@transform.T
            sl=slice(2*charge_index,2*charge_index+2)
            unit=old.old.old.VP.conj().T@old.old.old.spin.GAMMA[3]@old.old.old.VM
            observed=.5j*np.trace(g[:r,r:][sl,sl]@unit)
            analytic=np.tan(charge*angle/2)
            max_cov=max(max_cov,float(abs(observed-analytic)))
            weight_pf=old.old.old.pfaffian(q['N'])
            z+=weight*weight_pf;unnorm+=weight*weight_pf*rho*observed
        totals.append(z);means.append(unnorm/z)
    a=np.cos(theta)**2;b=np.cos(1.5*theta)**2
    da=-np.sin(2*theta);db=-1.5*np.sin(3*theta)
    dda=-2*np.cos(2*theta);ddb=-4.5*np.cos(3*theta)
    x=a*radial+b*(1-radial);dx=da*radial+db*(1-radial)
    ddx=dda*radial+ddb*(1-radial)
    w0=x**8;w1=8*x**7*dx;w2=56*x**6*dx**2+8*x**7*ddx
    z0,z1,z2=[weights@w for w in (w0,w1,w2)]
    m0,m1,m2=[weights@(radial*w) for w in (w0,w1,w2)]
    mu=m0/z0;dm=m1/z0-mu*z1/z0
    ddm=m2/z0-mu*z2/z0-2*dm*z1/z0
    o=np.tan(charge*theta/2);do=charge/2*(1+o*o)
    ddo=charge**2/2*(1+o*o)*o
    expected0=o*mu;expected1=do*mu+o*dm
    expected2=ddo*mu+2*do*dm+o*ddm
    fd1=(means[2]-means[0])/(2*h);fd2=(means[2]-2*means[1]+means[0])/h**2
    assert max_cov<2e-12 and abs(means[1]-expected0)<3e-13
    assert abs(fd1-expected1)<2e-7 and abs(fd2-expected2)<2e-6
    missing=expected1-do*mu;assert abs(missing)>.005
    return dict(theta=theta,charge=charge,auxiliary_observable='sum_{a=0}^5 E_a^2',
        scalar_Weyl_bilinear='i/2 trace(G_w_bar unit)',
        conditional_Weyl_covariance_error=max_cov,full_joint_mean=old.cpair(means[1]),
        factorized_joint_mean=float(expected0),first_derivative=old.cpair(fd1),
        analytic_first=float(expected1),second_derivative=old.cpair(fd2),
        analytic_second=float(expected2),first_error=float(abs(fd1-expected1)),
        second_error=float(abs(fd2-expected2)),
        omitted_common_auxiliary_source_response=float(missing),
        full_original_measure_and_transported_observation_used=True,
        gauge_RP_or_original_scalar_record_not_inferred=True)


def run():
    deps=('joint_local_mirror_process.py','joint_auxiliary_reflection_gluing.py',
          'joint_subgroup_measure_source.py','research_note_612.md','research_note_624.md',
          'research_note_647.md','research_note_658.md','research_note_660.md',
          'round661_drafts/observable_interface_probe_results.json')
    return dict(date='2026-10-02',round=661,tests_run=2,failures=0,errors=0,
        original_joint_state=original_joint_gram_check(),common_background_source=common_source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original local physical Weyl algebra times original auxiliary E algebra share a normalized finite free reflection functional, actual signed source insertions and common background responses. Not all mirror observations, original scalar instruments, general gauge RP, original Gibbs H or quantum GR.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=661,tests=result['tests_run'],passed=result['all_checks_passed'])))
