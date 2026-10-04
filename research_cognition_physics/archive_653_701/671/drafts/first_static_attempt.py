"""671: matrix spectral reflection bridge for the original static gauge field.

Analytic proof is in the numbered note. Finite diagnostics below retain full
noncommuting SM representation, signed auxiliary weights and physical fields.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_nonflat_mass_measure as prior
old=prior.old;internal=prior.internal;mass=prior.mass
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_static_gauge_reflection_results.json'


def spatial(m0=.73):
    sites=[(x,y) for x in range(2) for y in range(2)];n=64*len(sites)
    links=[prior.rep(*prior.group(67101,.13)),prior.rep(*prior.group(67102,.12))]
    ws=np.zeros((n,n),complex);c=np.zeros_like(ws)
    for mu in range(2):
        shift=np.zeros_like(ws)
        for i,(x,y) in enumerate(sites):
            dest=((x+1)%2,y) if mu==0 else (x,(y+1)%2);j=sites.index(dest)
            shift[64*i:64*i+64,64*j:64*j+64]=np.kron(np.eye(4),links[mu])
        gamma=np.kron(np.kron(np.eye(4),internal.spin.GAMMA[mu]),np.eye(16))
        ws+=np.eye(n)-(shift+shift.conj().T)/2
        c+=gamma@(shift-shift.conj().T)/2
    g0=np.kron(np.kron(np.eye(4),internal.spin.GAMMA[3]),np.eye(16))
    g5=np.kron(np.kron(np.eye(4),internal.spin.G5),np.eye(16))
    b=(1-m0)*np.eye(n)+ws;hs=g0@(b+c)
    assert max(old.err(b-b.conj().T),old.err(c+c.conj().T),old.err(hs-hs.conj().T))<3e-13
    assert old.err(g0@c+c@g0)<3e-13
    assert old.norm(b@c-c@b)>.001
    return dict(B=b,C=c,Hs=hs,g0=g0,g5=g5,n=n,links=links,
        curvature=old.norm(links[0]@links[1]@links[0].conj().T@links[1].conj().T-np.eye(16)))


def momentum(s,p):
    x=s['B']-np.cos(p)*np.eye(s['n'])+s['C']+1j*s['g0']*np.sin(p)
    h=s['g5']@x;ev,vec=np.linalg.eigh(h)
    assert old.err(h-h.conj().T)<3e-13 and min(abs(ev))>.1
    d=(np.eye(s['n'])+s['g5']@(vec*np.sign(ev))@vec.conj().T)/2
    return x,d,ev,vec


def groups(values):
    """Resolve numerical degeneracies without dropping any eigenspace."""
    out=[];start=0
    for i in range(1,len(values)):
        if abs(values[i]-values[start])>3e-10:
            out.append(np.arange(start,i));start=i
    out.append(np.arange(start,len(values)))
    assert sum(map(len,out))==len(values)
    return out


def spectral_density_check():
    s=spatial();b,c,hs,g0=(s[k] for k in ('B','C','Hs','g0'))
    ev,vec=np.linalg.eigh(b);bi=(vec/np.sqrt(ev))@vec.conj().T
    assert min(ev)>.26
    rows=[]
    for t in (.2,1.,7.):
        a=np.eye(s['n'])+hs@hs+t*np.eye(s['n'])
        lam,z=np.linalg.eigh(bi@a@bi);assert min(lam)>2
        kernels=[np.zeros_like(b) for _ in range(3)]
        maxherm=0.;minimum=0.;comm=0.;nullerror=0.;discmin=1e9
        for idx in groups(lam):
            co=float(np.mean(lam[idx])/2);en=float(np.arccosh(co));si=np.sinh(en)
            w=bi@z[:,idx];weight=w@w.conj().T;l=hs-g0*co
            basis=np.linalg.qr(w)[0]
            minus=(si*np.eye(s['n'])-l)@weight/(2*si)
            plus=(si*np.eye(s['n'])+l)@weight/(2*si)
            maxherm=max(maxherm,old.err(minus-minus.conj().T),old.err(plus-plus.conj().T))
            for rho in (minus,plus):
                small=basis.conj().T@rho@basis
                minimum=min(minimum,float(min(np.linalg.eigvalsh((small+small.conj().T)/2))))
            comm=max(comm,old.err(l@weight-weight@l))
            nullerror=max(nullerror,old.err((l@l-(si*si-t)*np.eye(s['n']))@basis))
            discmin=min(discmin,float(si*si-t))
            for k in range(3):kernels[k]+=np.exp(-(k+1)*en)*minus
        # Independent Fourier sum of the matrix resolvent, not of the formula.
        direct=[np.zeros_like(b) for _ in range(3)]
        # Use generalized spectral inverse here; original Q(p) identity checked separately.
        maxidentity=0.
        for p in 2*np.pi*np.arange(128)/128:
            x=b-np.cos(p)*np.eye(s['n'])+c+1j*g0*np.sin(p)
            q=x.conj().T@x+t*np.eye(s['n'])
            invq=bi@(z/(lam-2*np.cos(p)))@z.conj().T@bi
            maxidentity=max(maxidentity,old.err(q-(a-2*b*np.cos(p))))
            for k in range(3):direct[k]+=-g0@x@invq*np.exp(1j*p*(k+1))/128
        coefficient_error=max(old.err(x-y) for x,y in zip(kernels,direct))
        assert max(maxherm,comm,nullerror,maxidentity)<2e-9
        assert minimum>-2e-10 and coefficient_error<2e-10
        rows.append(dict(resolvent_parameter=t,degenerate_groups=len(groups(lam)),
            smallest_energy=float(np.arccosh(min(lam)/2)),minimum_discriminant=discmin,
            density_minimum=minimum,density_hermiticity_error=maxherm,
            residue_commutator_error=comm,residue_nullspace_error=nullerror,
            pencil_identity_error=maxidentity,Fourier_coefficient_error=coefficient_error))
    return dict(m0=.73,spatial_spin_internal_dimension=s['n'],
        original_plaquette_defect=s['curvature'],B_C_commutator=old.norm(b@c-c@b),rows=rows,
        integral_positivity_is_analytic_not_inferred_from_sampled_t=True)


def finite_data(nt=4,m0=1.):
    s=spatial(m0);n=s['n'];r=n//2;dp=[];up=[];vp=[];gaps=[]
    momenta=(2*np.arange(nt)+1)*np.pi/nt
    for p in momenta:
        _,d,ev,vec=momentum(s,p);dp.append(d);up.append(vec[:,ev<0]);vp.append(vec[:,ev>0]);gaps.append(float(min(abs(ev))))
    ft=np.exp(1j*np.outer(np.arange(nt),momenta))/np.sqrt(nt)
    d=np.zeros((nt*n,nt*n),complex);u=np.zeros((nt*n,nt*r),complex);v=np.zeros_like(u)
    for k in range(nt):
        for t in range(nt):
            u[t*n:(t+1)*n,k*r:(k+1)*r]=ft[t,k]*up[k]
            v[t*n:(t+1)*n,k*r:(k+1)*r]=ft[t,k]*vp[k]
            for j in range(nt):d[t*n:(t+1)*n,j*n:(j+1)*n]+=ft[t,k]*ft[j,k].conjugate()*dp[k]
    return dict(D=d,u=u,v=v,n=n,nt=nt,gap=min(gaps),spatial=s)


def joint_static_check():
    q=finite_data();d=q['D'];nt=q['nt'];ns=q['n'];count=4*nt
    # Direct action cross-kernel, including both AP images, on a two-slice half.
    positive=np.arange((nt//2)*ns,nt*ns)
    reflected=np.concatenate([np.arange((nt-1-t)*ns,(nt-t)*ns) for t in range(nt//2,nt)])
    g0half=np.kron(np.eye(nt//2),q['spatial']['g0'])
    cross=-g0half@d[np.ix_(reflected,positive)]
    herm=old.err(cross-cross.conj().T);minimum=float(min(np.linalg.eigvalsh((cross+cross.conj().T)/2)))
    assert herm<3e-12 and minimum>-3e-12
    jm=np.kron(np.kron(np.eye(count),internal.VM),np.eye(16));jp=np.kron(np.kron(np.eye(count),internal.VP),np.eye(16))
    kl=jp.conj().T@d@q['v'];w=jm.conj().T@q['v'];r=len(kl)
    l=old.diag(w,np.eye(r));z=np.zeros_like(kl);nw0=np.block([[z,-kl.T],[kl,z]])
    pair,_=mass.mass_pairing(count,mass.car.PHI[0]);nw=nw0+.37*l.T@pair@l
    covariance=-l@np.linalg.solve(nw,l.T)
    reflect=np.eye(count)[[4*(nt-1-t)+x for t in range(nt) for x in range(4)]]
    tr=np.kron(reflect,np.eye(32));theta=np.block([[np.zeros_like(tr),tr],[tr,np.zeros_like(tr)]])
    idx=np.flatnonzero(np.tile(np.repeat([t>=nt//2 for t in range(nt) for _ in range(4)],32),2))
    gram=(theta@covariance)[np.ix_(idx,idx)]
    gh=old.err(gram-gram.conj().T);gm=float(min(np.linalg.eigvalsh((gram+gram.conj().T)/2)))
    assert gh<4e-12 and gm>-4e-12
    # Full actual auxiliary Pfaffian kernel on three half-field configurations.
    fields=np.random.default_rng(67119).normal(size=(3,count//2,10))*.15;fields[:,:,0]+=1
    fields/=np.linalg.norm(fields,axis=2)[:,:,None]
    def auxiliary(e):
        mat=np.zeros_like(d)
        for i,ei in enumerate(e):
            t=sum(x*a for x,a in zip(ei,internal.T))
            mat[64*i:64*i+64,64*i:64*i+64]=np.kron(internal.B,t)
        return internal.pfaffian(q['u'].T@mat@q['u'])
    reference=auxiliary(np.tile(np.eye(10)[0],(count,1)))
    assert abs(reference)>1e-30
    ek=np.zeros((3,3),complex)
    for i in range(3):
        for j in range(3):
            negative=fields[i].reshape(nt//2,4,10)[::-1].reshape(count//2,10)
            ek[i,j]=auxiliary(np.concatenate((negative,fields[j])))/reference
    eh=old.err(ek-ek.conj().T);em=float(min(np.linalg.eigvalsh((ek+ek.conj().T)/2)))
    assert eh<4e-11 and em>-4e-11
    packets=np.random.default_rng(67120).normal(size=(len(idx),3))+1j*np.random.default_rng(67121).normal(size=(len(idx),3))
    packets/=np.linalg.norm(packets,axis=0)
    combined=ek*(packets.conj().T@gram@packets)
    cm=float(min(np.linalg.eigvalsh((combined+combined.conj().T)/2)))
    assert cm>-4e-11
    return dict(nx=2,ny=2,nt=nt,m0=1.,mass_scale=.37,Wilson_gap=q['gap'],
        action_cross_minimum=minimum,action_cross_hermiticity_error=herm,
        physical_mass_Gram_minimum=gm,physical_mass_Gram_hermiticity_error=gh,
        auxiliary_sample_minimum=em,auxiliary_sample_hermiticity_error=eh,
        combined_sample_minimum=cm,auxiliary_reference=old.cpair(reference),
        sample_not_complete_S9_integral_or_dynamic_gauge_integral=True)


def run():
    deps=('joint_nonflat_mass_measure.py','joint_mass_auxiliary_reflection.py',
          'research_note_605.md','research_note_660.md','research_note_661.md','research_note_670.md',
          'round671_drafts/magnetic_reflection_probe.py','round671_drafts/magnetic_reflection_probe_results.json')
    return dict(date='2026-10-02',round=671,tests_run=2,failures=0,errors=0,
        static_matrix_spectral_density=spectral_density_check(),
        actual_joint_static_reflection=joint_static_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Finite static spatial gauge backgrounds in temporal gauge, original Wilson/overlap at 0<m0<=1 with finite AP kernel nonsingular. Matrix spectral density extends the reflection argument beyond commuting spatial momenta. Original auxiliary/physical mass functional is conditionally unnormalized RP; strict normalization still requires nonzero partition function. No dynamic gauge path integral, original Gibbs/instrument identity, continuum or quantum GR completion.',all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
