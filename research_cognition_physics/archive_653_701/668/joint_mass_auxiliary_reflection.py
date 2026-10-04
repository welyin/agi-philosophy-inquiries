"""668: original nonzero masses, finite Weyl reflection and common sources.

Unit gauge links, constant original scalar backgrounds. No identity with the
complete original interacting Gibbs process is assumed or inferred.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_physical_auxiliary_state as prior
import joint_spinor_subgroup_mass as dictionary
import joint_gauss_fermion_influence as car
old=prior.old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_mass_auxiliary_reflection_results.json'


def mass_pairing(count,phi):
    w=dictionary.ph_matrix();b=w@dictionary.bdg(phi)@w.conj().T
    j=np.kron(dictionary.dictionary(),np.eye(2))
    delta=j@b[:32,32:]@j.T
    delta=delta.reshape(16,2,16,2).transpose(1,0,3,2).reshape(32,32)
    m=delta[:16,16:]
    assert old.err(delta-np.kron(dictionary.EPS,m))<1e-13
    assert old.err(m-m.T)<1e-13
    p=np.kron(np.eye(count),delta.conj())
    return old.diag(p,-p.conj()),np.linalg.svd(m,compute_uv=False)


def physical(nx,nt,phi,kappa=1.):
    frame=old.old.frame(kappa,nx,nt);count=nx*nt
    v=np.kron(frame['vec'][:,frame['ev']>0],np.eye(16))
    jm=np.kron(np.kron(np.eye(count),old.old.old.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(count),old.old.old.VP),np.eye(16))
    sign=(frame['vec']*np.sign(frame['ev']))@frame['vec'].conj().T
    g5=np.kron(np.eye(count),old.old.old.spin.G5)
    d=np.kron((np.eye(4*count)+g5@sign)/2,np.eye(16))
    kl=jp.conj().T@d@v;w=jm.conj().T@v;r=len(kl)
    z=np.zeros_like(kl);n0=np.block([[z,-kl.T],[kl,z]])
    select=old.diag(w,np.eye(r))
    pair,masses=mass_pairing(count,phi)
    return n0,select.T@pair@select,select,masses


def momentum_r2(nx,nt,kappa=1.):
    out=[]
    for t in range(nt):
        p4=(2*t+1)*np.pi/nt
        for x in range(nx):
            p1=2*x*np.pi/nx
            b=-1+(1-np.cos(p4))+kappa*(1-np.cos(p1))
            s2=np.sin(p4)**2+(kappa*np.sin(p1))**2
            omega=np.sqrt(b*b+s2)
            out.append(s2/(omega+b)**2)
    return np.array(out)


def massive_partition_check():
    rows=[];pf=old.old.old.pfaffian
    for nx,nt,phi,lam in ((1,2,car.PHI[0],.37),
                           (2,4,car.PHI[0],.37),
                           (3,4,car.PHI[1],1.1)):
        n0,pm,select,masses=physical(nx,nt,phi)
        actual=pf(n0+lam*pm)/pf(n0)
        r2=momentum_r2(nx,nt)
        logpred=float(np.log1p(lam**2*r2[:,None]*masses[None,:]**2).sum())
        pred=np.exp(logpred)
        relative=float(abs(actual/pred-1));assert relative<3e-11
        deriv=float((2*lam*r2[:,None]*masses[None,:]**2/
                     (1+lam**2*r2[:,None]*masses[None,:]**2)).sum())
        trace=.5*np.trace(np.linalg.solve(n0+lam*pm,pm))
        source_error=float(abs(trace-deriv));assert source_error<3e-11
        g=select@(-np.linalg.inv(n0+lam*pm))@select.T
        sites=[(t,x) for t in range(nt) for x in range(nx)]
        rr=np.eye(nx*nt)[[sites.index((nt-1-t,x)) for t,x in sites]]
        a=np.kron(rr,np.eye(32));z=np.zeros_like(a)
        theta=np.block([[z,a],[a,z]])
        idx=np.flatnonzero(np.tile(np.repeat([t>=nt//2 for t,x in sites],32),2))
        gram=(theta@g)[np.ix_(idx,idx)]
        herm=old.err(gram-gram.conj().T)
        mineig=float(min(np.linalg.eigvalsh((gram+gram.conj().T)/2)))
        assert herm<3e-12 and mineig>-3e-12
        rows.append(dict(nx=nx,nt=nt,mass_scale=lam,physical_Nambu_dimension=len(n0),
             actual_weight_ratio=old.cpair(actual),analytic_log_ratio=logpred,
             Fourier_product_relative_error=relative,mass_source_error=source_error,
             mass_source=deriv,reflection_Gram_minimum=mineig,reflection_hermiticity_error=herm))
    return dict(rows=rows,all_constant_mass_real_scales_positive_by_analytic_product=True,
                finite_Gram_diagnostics_not_used_as_full_RP_proof=True)


def common_background_check():
    nx=2;nt=2;lam=.37
    e=np.random.default_rng(66821).normal(size=(nx*nt,10))*.16;e[:,0]+=1
    e/=np.linalg.norm(e,axis=1)[:,None]
    direction=np.array([.1,-.07,.05,.03,-.12])
    def data(theta):
        kappa=1+.2*theta;phi=car.PHI[0]+theta*direction
        q=old.data(value=kappa,nx=nx,nt=nt,e=e);w=prior.weyl(q)
        p,_=mass_pairing(nx*nt,phi);f=w['map']
        n=q['N']+lam*f.T@p@f
        return q['N'],n,f,p
    step=2e-5;low=data(-step);mid=data(0);high=data(step)
    dn=(high[0]-low[0])/(2*step);df=(high[2]-low[2])/(2*step)
    dp=(high[3]-low[3])/(2*step);f,p=mid[2:]
    actual_derivative=(high[1]-low[1])/(2*step)
    complete=dn+lam*(df.T@p@f+f.T@dp@f+f.T@p@df)
    frozen_map=dn+lam*f.T@dp@f
    matrix_error=old.err(actual_derivative-complete);assert matrix_error<2e-10
    full=.5*np.trace(np.linalg.solve(mid[1],complete))
    missing=.5*np.trace(np.linalg.solve(mid[1],complete-frozen_map))
    assert abs(missing)>.001
    pf=old.old.old.pfaffian
    fd=(pf(high[1])-pf(low[1]))/(2*step*pf(mid[1]))
    error=float(abs(fd-full));assert error<3e-7
    # Independent original Weyl-mass product + original auxiliary/free weight.
    base=.5*np.trace(np.linalg.solve(mid[0],dn))
    logs=[]
    for t in (-step,step):
        _,m=mass_pairing(nx*nt,car.PHI[0]+t*direction)
        r2=momentum_r2(nx,nt,1+.2*t)
        logs.append(float(np.log1p(lam**2*r2[:,None]*m[None,:]**2).sum()))
    increment=(logs[1]-logs[0])/(2*step)
    factor_error=float(abs(full-base-increment));assert factor_error<3e-8
    return dict(source_family='kappa=1+.2 theta, phi=PHI[0]+theta*(.1,-.07,.05,.03,-.12)',
        actual_original_auxiliary_configuration=e.tolist(),full_conditional_log_source=old.cpair(full),
        original_auxiliary_and_free_source=old.cpair(base),mass_increment_source=increment,
        complete_derivative_matrix_error=matrix_error,full_source_finite_difference_error=error,
        factorized_mass_and_auxiliary_source_error=factor_error,
        omitted_observation_map_source=old.cpair(missing),
        conditional_source_not_full_bosonic_or_geometry_stress=True)


def run():
    deps=('joint_physical_auxiliary_state.py','joint_local_mirror_process.py',
          'joint_spinor_subgroup_mass.py','joint_gauss_fermion_influence.py',
          'research_note_604.md','research_note_614.md','research_note_658.md',
          'research_note_661.md','research_note_667.md',
          'round668_drafts/mass_observation_probe.py','round668_drafts/mass_observation_probe_results.json')
    return dict(date='2026-10-02',round=668,tests_run=2,failures=0,errors=0,
        constant_mass_normalization=massive_partition_check(),common_source=common_background_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Declared original constant mass pairing in retained physical Weyl algebra shares normalized finite free reflection functional and original auxiliary measure; all real strengths strictly positive by Fourier product. Exact local-candidate pullback and common background source, not bare local mass or original complete Gibbs/quantum GR equivalence.',
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
