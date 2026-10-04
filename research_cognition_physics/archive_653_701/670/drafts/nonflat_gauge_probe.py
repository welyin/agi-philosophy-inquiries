"""Actual670 entry: nonflat original hypercharge links and the same mass map.

Finite index-zero gapped patch. Algebra/covariance, not interacting gauge RP.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_mass_auxiliary_reflection as mass
old=mass.old;internal=old.old.old
SITES=[(t,x) for t in range(2) for x in range(2)]


def representation(angles,spin):
    out=np.zeros((len(angles)*16*spin,)*2,complex)
    for i,a in enumerate(angles):
        sl=slice(16*spin*i,16*spin*(i+1))
        out[sl,sl]=np.kron(np.eye(spin),np.diag(np.exp(1j*internal.Q*a)))
    return out


def vector_rotation(angle):
    r=np.diag(np.exp(1j*internal.Q*angle))
    columns=[]
    for t in internal.T:
        transformed=r.conj()@t@r.conj().T
        coeff=np.array([np.trace(a.conj().T@transformed)/16 for a in internal.T])
        assert max(abs(coeff.imag))<2e-13
        assert old.err(transformed-sum(c.real*a for c,a in zip(coeff,internal.T)))<2e-13
        columns.append(coeff.real)
    out=np.array(columns).T
    assert old.err(out.T@out-np.eye(10))<2e-13
    return out


def kernel(link_angles):
    n=256;dw=np.zeros((n,n),complex)
    for mu in range(2):
        shift=np.zeros((n,n),complex)
        for i,(t,x) in enumerate(SITES):
            dest=((t+1)%2,x) if mu==1 else (t,(x+1)%2)
            j=SITES.index(dest);sign=-1 if mu==1 and t==1 else 1
            block=sign*np.kron(np.eye(4),np.diag(np.exp(1j*internal.Q*link_angles[mu,i])))
            shift[64*i:64*i+64,64*j:64*j+64]=block
        gamma=np.kron(np.kron(np.eye(4),internal.spin.GAMMA[3 if mu==1 else 0]),np.eye(16))
        dw+=np.eye(n)-(shift+shift.conj().T)/2+gamma@(shift-shift.conj().T)/2
    g5=np.kron(np.kron(np.eye(4),internal.spin.G5),np.eye(16))
    h=g5@(dw-np.eye(n));ev,vec=np.linalg.eigh(h)
    assert old.err(h-h.conj().T)<1e-13 and min(abs(ev))>.5
    u=vec[:,ev<0];v=vec[:,ev>0]
    assert u.shape==v.shape==(256,128)
    d=(np.eye(n)+g5@(vec*np.sign(ev))@vec.conj().T)/2
    return d,u,v,float(min(abs(ev)))


def data(link_angles,e,phis):
    d,u,v,gap=kernel(link_angles)
    jm=np.kron(np.kron(np.eye(4),internal.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(4),internal.VP),np.eye(16))
    qp=jp@jp.conj().T;qm=jm@jm.conj().T
    m=np.zeros_like(d);mb=np.zeros_like(d)
    for i,ei in enumerate(e):
        t=sum(a*b for a,b in zip(ei,internal.T));sl=slice(64*i,64*i+64)
        m[sl,sl]=np.kron(internal.B,t);mb[sl,sl]=np.kron(internal.B,t.conj().T)
    n=old.diag(np.zeros_like(d),np.zeros_like(d))
    n[:256,:256]=m@qp;n[:256,256:]=-d.T;n[256:,:256]=d;n[256:,256:]=mb@qm
    s=old.diag(np.column_stack((u,v)),np.column_stack((jm.conj(),jp.conj())))
    a=u.T@m@u;bb=jm.conj().T@mb@jm.conj();k=jm.conj().T@d@u;kl=jp.conj().T@d@v
    c=u.T@m@qp@v;f=v.T@m@qp@v;r=128
    transform=np.eye(4*r,dtype=complex)
    transform[2*r:3*r,:r]=np.linalg.solve(bb,k)
    transform[3*r:,:r]=np.linalg.solve(kl.T,c.T)
    transform[3*r:,r:2*r]=-.5*np.linalg.solve(kl.T,f)
    select=np.zeros((2*r,4*r),complex)
    select[:r,r:2*r]=jm.conj().T@v;select[r:,3*r:]=np.eye(r)
    phi=select@transform@s.conj().T
    pair=np.zeros((2*r,2*r),complex)
    for i,ph in enumerate(phis):
        p,_=mass.mass_pairing(1,ph);sl=slice(32*i,32*i+32);sr=slice(r+32*i,r+32*i+32)
        pair[sl,sl]=p[:32,:32];pair[sr,sr]=p[32:,32:]
    n_mass=n+.37*phi.T@pair@phi
    expected=np.zeros_like(n);expected[:r,:r]=a;expected[2*r:3*r,2*r:3*r]=bb
    expected[r:2*r,3*r:]=-kl.T;expected[3*r:,r:2*r]=kl
    ti=np.linalg.inv(transform)
    error=old.err(ti.T@s.T@n@s@ti-expected)
    pf=internal.pfaffian
    fact=np.linalg.det(kl)*pf(a)/np.linalg.det(s)
    weight_error=float(abs(pf(n)/fact-1))
    return dict(D=d,N=n,Nmass=n_mass,Phi=phi,Pair=pair,gap=gap,
        triangular_error=error,original_weight_error=weight_error,
        signed_weight=pf(n),signed_mass_weight=pf(n_mass))


def run():
    links=np.zeros((2,4));links[0,0]=.11
    angles=np.array([-.03,.06,-.015,.03])
    rng=np.random.default_rng(67001);e=rng.normal(size=(4,10))*.1;e[:,0]+=1
    e/=np.linalg.norm(e,axis=1)[:,None]
    phis=np.array(mass.car.PHI)
    changed=links.copy()
    for mu in range(2):
        for i,(t,x) in enumerate(SITES):
            dest=((t+1)%2,x) if mu==1 else (t,(x+1)%2)
            changed[mu,i]+=angles[i]-angles[SITES.index(dest)]
    et=np.array([vector_rotation(a)@ei for a,ei in zip(angles,e)])
    pt=[]
    for a,ph in zip(angles,phis):
        x=np.exp(3j*a)*(ph[:2]+1j*ph[2:4]);pt.append(np.r_[x.real,x.imag,ph[4]])
    q=data(links,e,phis);qt=data(changed,et,np.array(pt))
    r4=representation(angles,4);r2=representation(angles,2)
    gx=old.diag(r4,r4.conj());gy=old.diag(r2,r2.conj())
    errs=dict(kernel=old.err(qt['D']-r4@q['D']@r4.conj().T),
        original_pair=old.err(gx.T@qt['N']@gx-q['N']),
        physical_observation=old.err(qt['Phi']@gx-gy@q['Phi']),
        original_mass=old.err(gy.T@qt['Pair']@gy-q['Pair']),
        full_mass_kernel=old.err(gx.T@qt['Nmass']@gx-q['Nmass']))
    assert max(errs.values())<3e-11
    assert max(q['triangular_error'],qt['triangular_error'])<3e-11
    assert max(q['original_weight_error'],qt['original_weight_error'])<3e-10
    weight=float(abs(qt['signed_weight']/q['signed_weight']-1))
    mweight=float(abs(qt['signed_mass_weight']/q['signed_mass_weight']-1))
    assert max(weight,mweight)<3e-10
    plaquette=links[0,0]+links[1,1]-links[0,2]-links[1,0]
    assert abs(plaquette-.11)<1e-14
    deps=('joint_mass_auxiliary_reflection.py','joint_local_mirror_process.py',
          'joint_subgroup_measure_source.py','research_note_660.md','research_note_668.md','research_note_669.md')
    return dict(date='2026-10-02',entry_for_round=670,formal_round_complete=False,
        base_unit_charge_plaquette_angle=float(plaquette),nonzero_curvature=True,
        original_SM_hypercharges=internal.Q.tolist(),original_kernel_gap=q['gap'],
        transformed_kernel_gap=qt['gap'],covariance_errors=errs,
        triangular_error=max(q['triangular_error'],qt['triangular_error']),
        complete_original_weight_error=max(q['original_weight_error'],qt['original_weight_error']),
        full_signed_weight_covariance_error=weight,full_mass_weight_covariance_error=mweight,
        original_mass_increment=old.cpair(q['signed_mass_weight']/q['signed_weight']),
        general_gauge_RP_or_dynamic_Gauss_process_not_claimed=True,
        dependency_hashes={p:hashlib.sha256((BASE/p).read_bytes()).hexdigest() for p in deps},
        all_checks_passed=True)


if __name__=='__main__':
    result=run();target=HERE/'nonflat_gauge_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(target.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('entry_for_round','formal_round_complete','all_checks_passed')}))
