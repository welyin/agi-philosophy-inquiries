"""673: rectangular physical dictionary and full two-boundary Gauss functional.

The note proves existence/integrability of the finite unnormalized candidate.
Numerics check actual sewing identities, not the full Haar integral or RP.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_nonflat_mass_measure as prior
old=prior.old;internal=prior.internal;mass=prior.mass
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauss_boundary_functional_results.json'


def pf(a):
    return 0j if len(a)%2 else internal.pfaffian(a)


def fixed_matrices(e,phis):
    count=len(e);n=64*count;r=n//2
    jm=np.kron(np.kron(np.eye(count),internal.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(count),internal.VP),np.eye(16))
    m=np.zeros((n,n),complex);mb=np.zeros_like(m);pair=np.zeros_like(m)
    for i,(ei,phi) in enumerate(zip(e,phis)):
        t=sum(x*a for x,a in zip(ei,internal.T));sl=slice(64*i,64*i+64)
        m[sl,sl]=np.kron(internal.B,t)
        mb[sl,sl]=np.kron(internal.B,t.conj().T)
        p,_=mass.mass_pairing(1,phi)
        a=slice(32*i,32*i+32);b=slice(r+32*i,r+32*i+32)
        pair[a,a]=p[:32,:32];pair[b,b]=p[32:,32:]
    return jm,jp,m,mb,pair


def regular(u,v,d,matrices,lam=.37,audit=False):
    jm,jp,m,mb,pair=matrices
    n=len(d);r=n//2;ru=u.shape[1];rv=v.shape[1]
    pu=u@u.conj().T;pv=v@v.conj().T;qp=jp@jp.conj().T;qm=jm@jm.conj().T
    phi=np.block([[jm.conj().T@pv,np.zeros((r,n))],
                  [-jp.T@m@(pu+.5*pv),jp.T]])
    nc=np.block([[m@qp,-d.T],[d,mb@qm]])
    full=nc+lam*phi.T@pair@phi
    out=dict(weight=pf(full),Phi=phi,N=full)
    if not audit:return out
    kl=jp.conj().T@d@v;w=jm.conj().T@v
    a=u.T@m@u;bb=jm.conj().T@mb@jm.conj();k=jm.conj().T@d@u
    s=old.diag(np.column_stack((u,v)),np.column_stack((jm.conj(),jp.conj())))
    t=np.eye(2*n,dtype=complex)
    t[n:n+r,:ru]=np.linalg.solve(bb,k)
    t[n+r:,:ru]=-jp.T@m@u
    t[n+r:,ru:n]=-.5*jp.T@m@v
    l=np.zeros((n,rv+r),complex);l[:r,:rv]=w;l[r:,rv:]=np.eye(r)
    select=np.zeros((n,2*n),complex)
    select[:r,ru:n]=w;select[r:,n+r:]=np.eye(r)
    nw0=np.block([[np.zeros((rv,rv)),-kl.T],[kl,np.zeros((r,r))]])
    nw=nw0+lam*l.T@pair@l
    expected=np.zeros_like(nc)
    expected[:ru,:ru]=a;expected[n:n+r,n:n+r]=bb
    expected[ru:n,n+r:]=-kl.T;expected[n+r:,ru:n]=kl
    ti=np.linalg.inv(t)
    errors=dict(triangular=old.err(ti.T@s.T@nc@s@ti-expected),
        rectangular_observation=old.err(phi-select@t@s.conj().T),
        projector_identity=old.err(d@v-qp@v),bar_pair_pf=float(abs(pf(bb)-1)))
    predicted=pf(a)*pf(nw)/np.linalg.det(s)
    denom=max(abs(out['weight']),abs(predicted),1e-24)
    errors['weight_absolute']=float(abs(out['weight']-predicted))
    if ru%2==0:errors['weight_relative']=float(abs(out['weight']-predicted)/denom)
    rng=np.random.default_rng(67311)
    sources=rng.normal(size=(n,4))+1j*rng.normal(size=(n,4))
    sources/=np.linalg.norm(sources,axis=0)
    source_errors=[]
    for count in (2,4):
        z=sources[:,:count]
        left=np.block([[full,phi.T@z],[-z.T@phi,np.zeros((count,count))]])
        right=np.block([[nw,l.T@z],[-z.T@l,np.zeros((count,count))]])
        x=pf(left);y=pf(a)*pf(right)/np.linalg.det(s)
        source_errors.append(float(abs(x-y)/max(abs(x),abs(y),1e-24)) if ru%2==0 else float(abs(x-y)))
    assert max(errors.values())<2e-9 and max(source_errors)<2e-9
    out.update(errors=errors,source_errors=source_errors,ru=ru,rv=rv,
        physical_integration_dimension=rv+r,observation_rows=n,
        phi_operator_norm=float(np.linalg.norm(phi,2)),predicted_weight=predicted)
    return out


def rectangular_check():
    rng=np.random.default_rng(67312)
    z=rng.normal(size=(64,64))+1j*rng.normal(size=(64,64));v0,_=np.linalg.qr(z)
    g5=np.kron(internal.spin.G5,np.eye(16))
    matrices=fixed_matrices(np.eye(10)[0][None,:],mass.car.PHI[0][None,:])
    rows=[]
    for ru in (30,31,34):
        u=v0[:,:ru];v=v0[:,ru:]
        d=(np.eye(64)+g5@(v@v.conj().T-u@u.conj().T))/2
        q=regular(u,v,d,matrices,audit=True)
        if ru%2:assert abs(q['weight'])<1e-20
        assert q['phi_operator_norm']<(1+np.sqrt(5))/2+1e-12
        rows.append({k:q[k] for k in ('ru','rv','physical_integration_dimension','observation_rows',
                    'errors','source_errors','phi_operator_norm')}|
                    dict(weight=old.cpair(q['weight']),predicted_weight=old.cpair(q['predicted_weight'])))
    return dict(algebra_fixture_not_a_nonzero_topology_Wilson_configuration=True,
        original_16_channel_pair_and_mass=True,rows=rows)


def kernel(links):
    sites=prior.SITES;n=256;dw=np.zeros((n,n),complex)
    for mu in range(2):
        shift=np.zeros_like(dw)
        for i,(t,x) in enumerate(sites):
            dest=((t+1)%2,x) if mu==1 else (t,1-x);j=sites.index(dest)
            shift[64*i:64*i+64,64*j:64*j+64]=(-1 if mu==1 and t==1 else 1)*np.kron(np.eye(4),links[mu,i])
        gamma=np.kron(np.kron(np.eye(4),internal.spin.GAMMA[3 if mu==1 else 0]),np.eye(16))
        dw+=np.eye(n)-(shift+shift.conj().T)/2+gamma@(shift-shift.conj().T)/2
    g5=np.kron(np.kron(np.eye(4),internal.spin.G5),np.eye(16))
    h=g5@(dw-np.eye(n));ev,vec=np.linalg.eigh(h)
    gap=float(min(abs(ev)));assert gap>1e-10
    d=(np.eye(n)+g5@(vec*np.sign(ev))@vec.conj().T)/2
    return vec[:,ev<0],vec[:,ev>0],d,h,gap


def configuration(links,e,phis):
    u,v,d,h,gap=kernel(links)
    q=regular(u,v,d,fixed_matrices(e,phis))
    q.update(gap=gap,H=h)
    return q


def change_phi(phi,g):
    _,w,z=g;h=z**3*w@(phi[:2]+1j*phi[2:4])
    return np.r_[h.real,h.imag,phi[4]]


def transform(links,e,phis,groups):
    rs=[prior.rep(*g) for g in groups];out=links.copy()
    for mu in range(2):
        for i,(t,x) in enumerate(prior.SITES):
            dest=((t+1)%2,x) if mu==1 else (t,1-x);j=prior.SITES.index(dest)
            out[mu,i]=rs[i]@links[mu,i]@rs[j].conj().T
    es=np.array([prior.vector_rotation(r)@v for r,v in zip(rs,e)])
    ps=np.array([change_phi(p,g) for p,g in zip(phis,groups)])
    return out,es,ps,rs


def cosh_distance(a,b):
    fa=2-np.dot(a,a)/6;fb=2-np.dot(b,b)/6
    return (2-np.dot(a,b)/6)/np.sqrt(fa*fb)


def boundary_check():
    e=np.random.default_rng(67321).normal(size=(4,10));e/=np.linalg.norm(e,axis=1)[:,None]
    phis=mass.car.PHI.copy()
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i in range(4):links[0,i]=prior.rep(*prior.group(67322+i,.16))
    a=[prior.group(67326+i,.32) for i in range(2)]
    b=[prior.group(67328+i,.38) for i in range(2)]
    for x in range(2):
        links[1,x]=prior.rep(*a[x]);links[1,2+x]=prior.rep(*b[x])
    q=configuration(links,e,phis)
    identity=(np.eye(3),np.eye(2),1.+0j)
    reduced,er,pr,rs=transform(links,e,phis,[identity,identity,*a])
    qr=configuration(reduced,er,pr)
    weight_error=float(abs(qr['weight']/q['weight']-1))
    assert weight_error<3e-10
    r4=np.zeros((256,256),complex);r2=np.zeros((128,128),complex)
    for i,r in enumerate(rs):
        r4[64*i:64*i+64,64*i:64*i+64]=np.kron(np.eye(4),r)
        r2[32*i:32*i+32,32*i:32*i+32]=np.kron(np.eye(2),r)
    gx=old.diag(r4,r4.conj());gy=old.diag(r2,r2.conj())
    source_map_error=old.err(qr['Phi']@gx-gy@q['Phi'])
    mass_congruence_error=old.err(gx.T@qr['N']@gx-q['N'])
    rng=np.random.default_rng(67333)
    j=rng.normal(size=(256,2))+1j*rng.normal(size=(256,2));j/=np.linalg.norm(j,axis=0)
    def coefficient(data,source):
        c=data['Phi'].T@source
        return pf(np.block([[data['N'],c],[-c.T,np.zeros((2,2))]]))
    left=coefficient(q,j);right=coefficient(qr,gy.conj()@j)
    source_relative_error=float(abs(left-right)/max(abs(left),abs(right)))
    assert max(source_map_error,mass_congruence_error,source_relative_error)<3e-10
    first=max(old.err(reduced[1,x]-np.eye(16)) for x in range(2))
    closure=max(old.err(reduced[1,2+x]-prior.rep(*a[x])@prior.rep(*b[x])) for x in range(2))
    wrong=configuration(reduced,e,phis)
    wrong_relative=float(abs(wrong['weight']/qr['weight']-1))
    assert wrong_relative>.001
    # Reflection swaps slices and inverts each time link at its fixed seam.
    reflected=links[:,[2,3,0,1]].copy()
    reflected[1]=np.array([z.conj().T for z in links[1]])
    qt=configuration(reflected,e[[2,3,0,1]],phis[[2,3,0,1]])
    reflection_error=float(abs(qt['weight']/q['weight'].conjugate()-1))
    assert reflection_error<3e-10
    distances=[]
    for x in range(2):
        before=cosh_distance(phis[2+x],change_phi(phis[x],b[x]))
        after=cosh_distance(pr[2+x],change_phi(change_phi(phis[x],b[x]),a[x]))
        distances.append(abs(before-after))
    assert max(distances)<3e-13 and max(first,closure)<3e-13
    # Actual large original gauge samples: diagnostic .5 gap gate is not retained.
    samples=[]
    for seed in (67350,67351,67352):
        arbitrary=np.array([[prior.rep(*prior.group(seed+8*mu+i,.9)) for i in range(4)] for mu in range(2)])
        qs=configuration(arbitrary,e,phis)
        identity_time=arbitrary.copy();identity_time[1]=np.tile(np.eye(16),(4,1,1))
        _,_,_,_,identity_gap=kernel(identity_time)
        singular=np.linalg.svd(qs['N'],compute_uv=False)
        reflected=arbitrary[:,[2,3,0,1]].copy()
        reflected[1]=np.array([z.conj().T for z in arbitrary[1]])
        qsref=configuration(reflected,e[[2,3,0,1]],phis[[2,3,0,1]])
        resolved=bool(singular[-1]>1e-10)
        relative=float(abs(qsref['weight']/qs['weight'].conjugate()-1)) if resolved else None
        if resolved:assert relative<2e-9
        samples.append(dict(seed=seed,gap=qs['gap'],identity_seam_gap=identity_gap,
            weight=old.cpair(qs['weight']),smallest_full_singular_value=float(singular[-1]),
            numerical_near_null_modes=int(sum(singular<1e-10)),
            weight_phase_or_relative_error_resolved=resolved,
            reflected_relative_error_only_when_resolved=relative,
            reflected_absolute_difference=float(abs(qsref['weight']-qs['weight'].conjugate()))))
    return dict(original_full_group_and_variable_original_masses=True,
        unreduced_weight=old.cpair(q['weight']),reduced_weight=old.cpair(qr['weight']),
        weight_identity_error=weight_error,first_seam_identity_error=first,
        source_map_error=source_map_error,mass_congruence_error=mass_congruence_error,
        transported_two_source_relative_error=source_relative_error,
        closure_product_error=closure,omitting_field_transport_relative_error=wrong_relative,
        reflected_weight_error=reflection_error,original_H5_distance_errors=distances,
        minimum_sewing_gap=min(q['gap'],qr['gap'],wrong['gap'],qt['gap']),
        actual_large_gauge_samples=samples,
        full_Haar_integral_not_numerically_evaluated=True,
        no_reflection_positivity_or_normalization_claim=True)


def run():
    deps=('joint_nonflat_mass_measure.py','joint_gauge_history_kernel.py','research_note_643.md',
          'research_note_669.md','research_note_670.md','research_note_672.md',
          'round673_drafts/boundary_holonomy_probe_results.json')
    return dict(date='2026-10-02',round=673,tests_run=2,failures=0,errors=0,
        rectangular_source_dictionary=rectangular_check(),actual_double_boundary_sewing=boundary_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Original unnormalized physical-source dictionary extended algebraically to unequal chiral ranks; fixed even r. Actual two-slice gauge sewing transports scalar and auxiliary data together. Analytic finite double-Haar functional uses original bosonic heat kernel and configuration moments; Wilson zeros are Haar-null at two slices, with no hard admissibility cutoff imposed for integrability. No complete numerical Haar integral, positivity/nonzero partition function, identity with original ordered CAR influence, continuum limit or quantum GR claim.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=673,tests_run=2,all_checks_passed=True,
        boundary=result['actual_double_boundary_sewing'])))
