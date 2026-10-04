"""670: joint nonflat gauge/mass measure and the massless-coordinate singularity.

Finite index-zero Wilson-gap patches only. No interacting reflection positivity,
global measure, Hamiltonian equivalence or continuum limit is inferred.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_mass_auxiliary_reflection as mass
old=mass.old
internal=old.old.old
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_nonflat_mass_measure_results.json'
SITES=[(t,x) for t in range(2) for x in range(2)]


def rep(c,w,z):
    j=mass.dictionary.dictionary()
    return j@mass.dictionary.left_rep(c,w,z)@j.conj().T


def group(seed,scale):
    rng=np.random.default_rng(seed)
    exp=mass.dictionary.old.gauge.group_exp
    return exp(scale*rng.normal(size=8),3),exp(scale*rng.normal(size=3),2),np.exp(1j*scale*rng.normal())


def vector_rotation(r):
    columns=[]
    for t in internal.T:
        target=r.conj()@t@r.conj().T
        coeff=np.array([np.trace(a.conj().T@target)/16 for a in internal.T])
        assert old.err(coeff.imag)<3e-13
        assert old.err(target-sum(x.real*a for x,a in zip(coeff,internal.T)))<3e-13
        columns.append(coeff.real)
    out=np.array(columns).T
    assert old.err(out.T@out-np.eye(10))<3e-13
    return out


def kernel(links):
    n=256;dw=np.zeros((n,n),complex)
    for mu in range(2):
        shift=np.zeros((n,n),complex)
        for i,(t,x) in enumerate(SITES):
            dest=((t+1)%2,x) if mu==1 else (t,(x+1)%2)
            j=SITES.index(dest);sign=-1 if mu==1 and t==1 else 1
            shift[64*i:64*i+64,64*j:64*j+64]=sign*np.kron(np.eye(4),links[mu,i])
        gamma=np.kron(np.kron(np.eye(4),internal.spin.GAMMA[3 if mu==1 else 0]),np.eye(16))
        dw+=np.eye(n)-(shift+shift.conj().T)/2+gamma@(shift-shift.conj().T)/2
    g5=np.kron(np.kron(np.eye(4),internal.spin.G5),np.eye(16))
    h=g5@(dw-np.eye(n));ev,vec=np.linalg.eigh(h)
    assert old.err(h-h.conj().T)<2e-13 and min(abs(ev))>.5
    u=vec[:,ev<0];v=vec[:,ev>0]
    assert u.shape==v.shape==(256,128)
    d=(np.eye(n)+g5@(vec*np.sign(ev))@vec.conj().T)/2
    return u,v,d,h


def joint(u,v,d,jm,jp,m,mb,pair,lam=.37,local_check=False):
    """Projector formula has no Kl inverse; optional legacy formula is a check."""
    r=u.shape[1];n=len(d)
    kl=jp.conj().T@d@v;w=jm.conj().T@v
    a=u.T@m@u;s=old.diag(np.column_stack((u,v)),np.column_stack((jm.conj(),jp.conj())))
    l=old.diag(w,np.eye(r));zero=np.zeros_like(kl)
    nw0=np.block([[zero,-kl.T],[kl,zero]])
    nw=nw0+lam*l.T@pair@l
    pf=internal.pfaffian
    aux=pf(a)/np.linalg.det(s);weight=aux*pf(nw)
    pu=u@u.conj().T;pv=v@v.conj().T;qp=jp@jp.conj().T
    assert old.err(d@v-qp@v)<3e-13
    assert old.err(m@qp-qp.T@m)<3e-13
    phi=np.block([[jm.conj().T@pv,np.zeros((r,n))],
                  [-jp.T@m@(pu+.5*pv),jp.T]])
    out=dict(kl=kl,L=l,A=a,Nw=nw,S=s,aux=aux,weight=weight,Phi=phi,
             covariance=-l@np.linalg.solve(nw,l.T),mass_source=.5*np.trace(np.linalg.solve(nw,l.T@pair@l)))
    if local_check:
        nc=np.block([[m@jp@jp.conj().T,-d.T],[d,mb@jm@jm.conj().T]])
        bb=jm.conj().T@mb@jm.conj();k=jm.conj().T@d@u
        c=u.T@m@jp@jp.conj().T@v;f=v.T@m@jp@jp.conj().T@v
        transform=np.eye(4*r,dtype=complex)
        transform[2*r:3*r,:r]=np.linalg.solve(bb,k)
        transform[3*r:,:r]=-jp.T@m@u
        transform[3*r:,r:2*r]=-.5*jp.T@m@v
        select=np.zeros((2*r,4*r),complex)
        select[:r,r:2*r]=w;select[r:,3*r:]=np.eye(r)
        assert old.err(phi-select@transform@s.conj().T)<4e-13
        expected=np.zeros_like(nc)
        expected[:r,:r]=a;expected[2*r:3*r,2*r:3*r]=bb
        expected[r:2*r,3*r:]=-kl.T;expected[3*r:,r:2*r]=kl
        ti=np.linalg.inv(transform)
        triangular_error=old.err(ti.T@s.T@nc@s@ti-expected)
        assert triangular_error<4e-13
        # Only compare the old inverse formula where numerically nonsingular.
        legacy_error=None
        if min(np.linalg.svd(kl,compute_uv=False))>1e-12:
            legacy=transform.copy()
            legacy[3*r:,:r]=np.linalg.solve(kl.T,c.T)
            legacy[3*r:,r:2*r]=-.5*np.linalg.solve(kl.T,f)
            legacy_error=old.err(phi-select@legacy@s.conj().T)
            assert legacy_error<3e-11
        massive=nc+lam*phi.T@pair@phi
        out.update(local_weight=pf(massive),local_weight_error=float(abs(pf(massive)/weight-1)),
                   triangular_error=triangular_error,legacy_formula_error=legacy_error,
                   local_covariance_error=old.err(-phi@np.linalg.solve(massive,phi.T)-out['covariance']))
    return out


def spatial_data(links,e,phis,rephase=False,local_check=False):
    u,v,d,h=kernel(links)
    if rephase:
        u=u*np.exp(1j*np.linspace(-.3,.4,128));v=v*np.exp(1j*np.linspace(.2,-.7,128))
    jm=np.kron(np.kron(np.eye(4),internal.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(4),internal.VP),np.eye(16))
    m=np.zeros_like(d);mb=np.zeros_like(d);pair=np.zeros((256,256),complex)
    for i,(ei,ph) in enumerate(zip(e,phis)):
        t=sum(x*a for x,a in zip(ei,internal.T));sl=slice(64*i,64*i+64)
        m[sl,sl]=np.kron(internal.B,t);mb[sl,sl]=np.kron(internal.B,t.conj().T)
        p,_=mass.mass_pairing(1,ph);a=slice(32*i,32*i+32);b=slice(128+32*i,128+32*i+32)
        pair[a,a]=p[:32,:32];pair[b,b]=p[32:,32:]
    out=joint(u,v,d,jm,jp,m,mb,pair,local_check=local_check)
    out.update(D=d,gap=float(min(abs(np.linalg.eigvalsh(h)))))
    return out


def nonflat_dictionary_check():
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    links[0,0]=rep(*group(67031,.035));links[1,1]=rep(*group(67032,.04))
    plaquette=links[0,0]@links[1,1]@links[0,2].conj().T@links[1,0].conj().T
    curvature=old.norm(plaquette-np.eye(16));commutator=old.norm(links[0,0]@links[1,1]-links[1,1]@links[0,0])
    assert curvature>.1 and commutator>.001
    gs=[group(67040+i,.045) for i in range(4)];rs=[rep(*g) for g in gs]
    changed=links.copy()
    for mu in range(2):
        for i,(t,x) in enumerate(SITES):
            dest=((t+1)%2,x) if mu==1 else (t,(x+1)%2);j=SITES.index(dest)
            changed[mu,i]=rs[i]@links[mu,i]@rs[j].conj().T
    e=np.random.default_rng(67043).normal(size=(4,10))*.1;e[:,0]+=1;e/=np.linalg.norm(e,axis=1)[:,None]
    phis=np.array(mass.car.PHI);et=np.array([vector_rotation(r)@ei for r,ei in zip(rs,e)])
    pt=[]
    for (_,weak,z),ph in zip(gs,phis):
        h=z**3*weak@(ph[:2]+1j*ph[2:4]);pt.append(np.r_[h.real,h.imag,ph[4]])
    q=spatial_data(links,e,phis,local_check=True);qt=spatial_data(changed,et,pt)
    phase=spatial_data(links,e,phis,rephase=True)
    r4=np.zeros((256,256),complex);r2=np.zeros((128,128),complex)
    for i,r in enumerate(rs):
        r4[64*i:64*i+64,64*i:64*i+64]=np.kron(np.eye(4),r)
        r2[32*i:32*i+32,32*i:32*i+32]=np.kron(np.eye(2),r)
    gy=old.diag(r2,r2.conj())
    errors=dict(kernel=old.err(qt['D']-r4@q['D']@r4.conj().T),
        signed_weight=float(abs(qt['weight']/q['weight']-1)),
        physical_two_point=old.err(qt['covariance']-gy@q['covariance']@gy.T),
        mass_source=float(abs(qt['mass_source']-q['mass_source'])),
        frame_weight=float(abs(phase['weight']/q['weight']-1)),
        frame_two_point=old.err(phase['covariance']-q['covariance']),
        local_weight=q['local_weight_error'],local_two_point=q['local_covariance_error'])
    assert max(errors.values())<3e-10
    return dict(sites=4,internal_channels=16,nonabelian_plaquette_defect=curvature,
        link_commutator_norm=commutator,Wilson_gap=q['gap'],mass_scale=.37,
        complete_signed_weight=old.cpair(q['weight']),errors=errors)


def holonomy_data(theta,lam=.37,local_check=False):
    u,v,d,h=internal.frames(theta)
    jm=np.kron(np.eye(16),internal.VM);jp=np.kron(np.eye(16),internal.VP)
    pair,_=mass.mass_pairing(1,mass.car.PHI[0])
    perm=np.array([s*16+i for i in range(16) for s in range(2)]);idx=np.r_[perm,32+perm]
    pair=pair[np.ix_(idx,idx)]
    m=np.kron(internal.T[0],internal.B);mb=np.kron(internal.T[0].conj().T,internal.B)
    out=joint(u,v,d,jm,jp,m,mb,pair,lam=lam,local_check=local_check)
    out.update(gap=float(min(abs(np.linalg.eigvalsh(h)))),mass_pair=pair)
    return out


def zero_mode_extension_check():
    center=np.pi/6;q=holonomy_data(center,local_check=True);pf=internal.pfaffian
    s0=float(min(np.linalg.svd(q['kl'],compute_uv=False)))
    s1=float(min(np.linalg.svd(q['Nw'],compute_uv=False)))
    assert s0<1e-14 and s1>.001 and abs(q['weight'])>1e-10
    rows=[]
    for eps in (1e-2,1e-3,1e-4):
        a=holonomy_data(center+eps,local_check=True)
        assert a['local_weight_error']<2e-8 and a['local_covariance_error']<2e-8
        rows.append(dict(offset=eps,Kl_minimum=float(min(np.linalg.svd(a['kl'],compute_uv=False))),
            inverse_coordinate_norm=old.norm(a['Phi']),signed_weight=old.cpair(a['weight']),
            local_weight_error=a['local_weight_error'],local_two_point_error=a['local_covariance_error']))
    assert q['local_weight_error']<2e-8 and q['local_covariance_error']<2e-8
    assert max(row['inverse_coordinate_norm'] for row in rows)<8
    # Pfaffian cofactors define unnormalized source coefficients without inverse.
    nw=q['Nw'];cofactor=np.zeros_like(nw)
    for i,j in ((0,1),(20,21),(20,52),(32,33)):
        keep=np.array([k for k in range(64) if k not in (i,j)])
        cofactor[i,j]=(-1)**(i+j+1)*pf(nw[np.ix_(keep,keep)])
    exact=pf(nw)*(-np.linalg.inv(nw))
    cofactor_error=max(abs(cofactor[i,j]-exact[i,j]) for i,j in ((0,1),(20,21),(20,52),(32,33)))
    assert cofactor_error<2e-16
    # Differentiation across the massless zero uses the nonzero massive weight.
    step=2e-7;a=holonomy_data(center-step);b=holonomy_data(center+step)
    dn=(b['Nw']-a['Nw'])/(2*step)
    source=(b['aux']-a['aux'])/(2*step*q['aux'])+.5*np.trace(np.linalg.solve(nw,dn))
    finite=(b['weight']-a['weight'])/(2*step*q['weight'])
    error=float(abs(source-finite));assert error<.003
    mstep=1e-6
    mf=(holonomy_data(center,.37+mstep)['weight']-holonomy_data(center,.37-mstep)['weight'])/(2*mstep*q['weight'])
    mass_error=float(abs(mf-q['mass_source']));assert mass_error<2e-6
    return dict(theta=center,mass_scale=.37,Wilson_gap=q['gap'],massless_Kl_minimum=s0,
        massive_Nambu_minimum=s1,exact_zero_signed_weight=old.cpair(q['weight']),
        physical_covariance_norm=old.norm(q['covariance']),nearby_rows=rows,
        zero_mode_projector_map_norm=old.norm(q['Phi']),
        zero_mode_triangular_error=q['triangular_error'],
        zero_mode_local_weight_error=q['local_weight_error'],
        zero_mode_local_two_point_error=q['local_covariance_error'],
        cofactor_error=float(cofactor_error),holonomy_log_source=old.cpair(source),
        source_finite_difference_error=error,mass_source_error=mass_error,
        inverse_Kl_not_required_for_actual_observation_map_or_generating_polynomial=True,
        global_normalization_or_index_change_not_proved=True)


def run():
    deps=('joint_mass_auxiliary_reflection.py','joint_local_mirror_process.py',
          'joint_subgroup_measure_source.py','research_note_660.md','research_note_668.md','research_note_669.md',
          'round670_drafts/nonflat_gauge_probe.py','round670_drafts/nonflat_gauge_probe_results.json')
    return dict(date='2026-10-02',round=670,tests_run=2,failures=0,errors=0,
        nonflat_joint_dictionary=nonflat_dictionary_check(),
        physical_massless_zero_extension=zero_mode_extension_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Same finite index-zero Wilson-gap patch: full non-Abelian original SM links, original masses and physical-source polynomial share covariance and frame cancellation. The apparent inverse massless Weyl kinetic block in the actual observation map cancels exactly using the original chiral pairing; the projector expression extends across its zero. An original massive holonomy witness has nonzero weight there. No all-background nonzero weight, interacting RP, complete Gauss process identity, continuum or quantum GR claim.',
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
