"""940: positive physical modes + constraints for shared Dirac sources.

Canonical weak-field physics is inherited, not a new gravity derivation.
The finite Hamiltonian is an explicit mode diagnostic, not full E_rec.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, itertools, json, math
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'native_exchange_bridge_results.json'
OLD=HERE.parent/'931/charged_matter_response_results.json'
I2=np.eye(2,dtype=complex)
PAULI=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1.,-1.]).astype(complex)]
G0=np.diag([1.,1.,-1.,-1.]).astype(complex)
GAMMA=[G0]+[np.block([[np.zeros((2,2)),s],[-s,np.zeros((2,2))]]) for s in PAULI]
ETA=np.diag([1.,-1.,-1.,-1.])

def norm(a):return float(np.linalg.norm(a,2))
def real(v):
    z=complex(v);assert abs(z.imag)<3e-11
    return float(z.real)
def source(p,p1,m):
    def spinors(v):
        e=np.sqrt(m*m+v@v)
        return e,np.sqrt((e+m)/(2*e))*np.vstack([I2,sum(v[i]*PAULI[i] for i in range(3))/(e+m)])
    e,u=spinors(p);e1,u1=spinors(p1)
    pp=np.r_[e+e1,p+p1]
    j=np.array([u1.conj().T@G0@g@u for g in GAMMA])
    t=np.array([[.25*(j[a]*pp[b]+j[b]*pp[a]) for b in range(4)] for a in range(4)])
    q=np.r_[e1-e,p1-p]
    return q,j,t,e,e1

def contractions(q,t,s,j,l):
    k=q[1:];k2=float(k@k);w=q[0];den=k2-w*w
    assert k2>0 and abs(den)>1e-8
    nh=k/np.sqrt(k2);pt=np.eye(3)-np.outer(nh,nh)
    def tt(a):return pt@a@pt-.5*pt*np.trace(pt@a)
    def cn(a,b):return np.vdot(a,b)
    ttrace=np.trace(ETA@t);strace=np.trace(ETA@s)
    covg=(np.vdot(t,ETA@s@ETA)-.5*np.conj(ttrace)*strace)/den
    radg=cn(tt(t[1:,1:]),tt(s[1:,1:]))/den
    rho,rhop=t[0,0],s[0,0];sigma,sigmap=np.trace(t[1:,1:]),np.trace(s[1:,1:])
    jt,jtp=pt@t[0,1:],pt@s[0,1:]
    jl,jlp=nh@t[0,1:],nh@s[0,1:]
    instantg=(.5*np.conj(rho)*rhop+.5*(np.conj(rho)*sigmap+np.conj(sigma)*rhop)-2*cn(jt,jtp)-1.5*np.conj(jl)*jlp)/k2
    cove=-cn(j,ETA@l)/den
    rade=cn(pt@j[1:],pt@l[1:])/den
    instante=-np.conj(j[0])*l[0]/k2
    return covg,radg,instantg,cove,rade,instante

def exact_rational_check():
    # Coordinate proof at k along z with two independent conserved tensors.
    largest=F(0)
    for w,k in [(F(1,3),F(2)),(F(0),F(3,2)),(F(4,5),F(7,5))]:
        values=[]
        for shift in (0,2):
            rho,jx,jy,a,b,c=[F(i+shift,7) for i in (1,2,-3,4,5,-2)]
            t=[[F(0) for _ in range(4)] for _ in range(4)]
            t[0][0]=rho;t[0][1]=t[1][0]=jx;t[0][2]=t[2][0]=jy
            t[0][3]=t[3][0]=w*rho/k
            t[1][1]=a;t[2][2]=b;t[1][2]=t[2][1]=c
            t[3][1]=t[1][3]=w*jx/k;t[3][2]=t[2][3]=w*jy/k;t[3][3]=w*w*rho/(k*k)
            values.append(t)
        t,s=values;sgn=[1,-1,-1,-1]
        tr=lambda a:sum(sgn[i]*a[i][i] for i in range(4))
        den=k*k-w*w
        cov=(sum(sgn[i]*sgn[j]*t[i][j]*s[i][j] for i in range(4) for j in range(4))-tr(t)*tr(s)/2)/den
        rad=((t[1][1]-t[2][2])*(s[1][1]-s[2][2])/2+2*t[1][2]*s[1][2])/den
        sig=lambda a:sum(a[i][i] for i in (1,2,3))
        contact=(t[0][0]*s[0][0]/2+(t[0][0]*sig(s)+sig(t)*s[0][0])/2
                 -2*(t[0][1]*s[0][1]+t[0][2]*s[0][2])-F(3,2)*t[0][3]*s[0][3])/(k*k)
        largest=max(largest,abs(cov-rad-contact))
    assert largest==0
    return str(largest)

def fermion_sector():
    masks=[n for n in range(16) if n.bit_count()==2]
    cs=[]
    for i in range(4):
        c=np.zeros((16,16),complex)
        for m in range(16):
            if m>>i&1:c[m^(1<<i),m]=(-1)**((m&((1<<i)-1)).bit_count())
        cs.append(c)
    bs=np.array([[((cs[i].conj().T@cs[j])[np.ix_(masks,masks)]) for j in range(4)] for i in range(4)])
    dg=lambda a:np.einsum('ij,ijab->ab',a,bs)
    def normal(a,b):
        return (dg(a)@dg(b)+dg(b)@dg(a)-dg(a@b+b@a))/2
    return masks,dg,normal

def finite_hamiltonian(mass,charge):
    k=.17;px=.04;py=.02;boost=.2
    eb=np.sqrt(mass*mass+px*px+py*py+k*k/4)
    gam=1/np.sqrt(1-boost*boost)
    p=np.array([px,py,gam*(-k/2+boost*eb)])
    p1=np.array([px,py,gam*(k/2+boost*eb)])
    q,j,t,e0,e1=source(p,p1,mass);k=float(np.linalg.norm(q[1:]))
    masks,dg,normal=fermion_sector();nf=len(masks)
    hm=dg(np.diag([e0,e0,e1,e1]))
    def lift(a):
        b=np.zeros((4,4),complex);b[2:,0:2]=a
        return [(b+b.conj().T)/np.sqrt(2),(b-b.conj().T)/(1j*np.sqrt(2))]
    rho=lift(t[0,0]);sigma=lift(t[1,1]+t[2,2]+t[3,3])
    mom=[lift(t[0,i]) for i in (1,2,3)]
    re=lift(charge*j[0]);je=[lift(charge*j[i]) for i in (1,2)]
    tg=[lift((t[1,1]-t[2,2])/np.sqrt(2)),lift(np.sqrt(2)*t[1,2])]
    coulomb=sum(normal(a,a) for a in re)/(2*k*k)
    cg=np.zeros((nf,nf),complex)
    for r in range(2):
        cg+=(.5*normal(rho[r],rho[r])+normal(rho[r],sigma[r])
             -2*normal(mom[0][r],mom[0][r])-2*normal(mom[1][r],mom[1][r])
             -1.5*normal(mom[2][r],mom[2][r]))
    gravity=-cg/(8*k*k)
    currents=[dg(a) for pol in je for a in pol]
    stresses=[dg(a) for pol in tg for a in pol]
    all_sources=currents+stresses
    nb=8;shape=(2,)*nb+(nf,);dim=2**nb*nf
    aa=np.array([[0.,1.],[0.,0.]],complex)
    Q=(aa+aa.conj().T)/np.sqrt(2*k)
    P=-1j*np.sqrt(k/2)*(aa-aa.conj().T)
    top=np.diag([0.,1.])
    occupations=np.array([i.bit_count() for i in range(2**nb)],float)
    free_b=k*occupations
    def on_axis(mat,v,axis):
        return np.moveaxis(np.tensordot(mat,v.reshape(shape),axes=(1,axis)),0,axis).reshape(-1)
    def matter(a,v):return (v.reshape(-1,nf)@a.T).reshape(-1)
    def operators(kappa,em):
        hmat=hm+em*em*coulomb+kappa*kappa*gravity
        couplings=[em]*4+[kappa/2]*4
        def apply(v):
            out=matter(hmat,v)+(free_b[:,None]*v.reshape(-1,nf)).reshape(-1)
            for axis,(g,s) in enumerate(zip(couplings,all_sources)):
                out-=g*matter(s,on_axis(Q,v,axis))
            return out
        bnd=norm(hmat)+nb*k+sum(abs(g)*norm(Q)*norm(s) for g,s in zip(couplings,all_sources))
        def dk(v):
            out=matter(2*kappa*gravity,v)
            for axis,s in enumerate(stresses,4):out-=.5*matter(s,on_axis(Q,v,axis))
            return out
        return apply,bnd,dk,couplings
    kappa=.08;em=.05;T=6.
    H,bound,Dk,gs=operators(kappa,em)
    v0=np.zeros(dim,complex)
    v0[masks.index(3)]=1/np.sqrt(2);v0[masks.index(5)]=1/np.sqrt(2)
    def evolve(apply,bnd):
        steps=max(1,math.ceil(T*bnd/.3));dt=T/steps;degree=14
        error_one=math.exp(dt*bnd)*(dt*bnd)**(degree+1)/math.factorial(degree+1)
        v=v0.copy()
        for _ in range(steps):
            term=v.copy();total=v.copy()
            for order in range(1,degree+1):
                term=(-1j*dt/order)*apply(term);total+=term
            v=total
        return v,steps,error_one*steps*(1+error_one)**(steps-1)
    v,steps,remainder=evolve(H,bound)
    H0,b0,_,_=operators(0.,em);v_without,_,_=evolve(H0,b0)
    record=np.diag([float(bool(mask&12)) for mask in masks])
    probability=lambda state:real(np.vdot(state,matter(record,state)))
    force_error=0.;boundary_effect=0.
    for axis,(g,s) in enumerate(zip(gs,all_sources)):
        lhs=1j*(H(on_axis(P,v,axis))-on_axis(P,H(v),axis))
        ideal=-k*k*on_axis(Q,v,axis)+g*matter(s,v)
        correction=-2*g*matter(s,on_axis(top,v,axis))
        force_error=max(force_error,float(np.linalg.norm(lhs-ideal-correction)))
        boundary_effect=max(boundary_effect,abs(real(np.vdot(v,correction))))
    dj=1e-5;Hp,_,_,_=operators(kappa+dj,em);Hm,_,_,_=operators(kappa-dj,em)
    derivative_error=float(np.linalg.norm((Hp(v)-Hm(v))/(2*dj)-Dk(v)))
    energy_before=real(np.vdot(v0,H(v0)));energy_after=real(np.vdot(v,H(v)))
    transition=probability(v)-probability(v_without)
    assert abs(np.vdot(v,v)-1)<2e-12 and abs(energy_after-energy_before)<2e-12
    assert force_error<2e-12 and derivative_error<3e-10 and abs(transition)>1e-9
    return dict(dimension=dim,matter_sector='two particles in four electron momentum-spin modes',
        physical_bosonic_modes=8,occupation_cutoff_per_mode=2,k=k,kappa=kappa,em=em,time=T,
        generator_norm_upper=bound,time_steps=steps,Taylor_action_remainder_bound=remainder,
        norm_error=float(abs(np.vdot(v,v)-1)),energy_before=energy_before,energy_after=energy_after,
        energy_conservation_error=abs(energy_after-energy_before),
        probability_with_gravity=probability(v),probability_without_gravity=probability(v_without),
        gravity_record_difference=transition,actual_source_derivative_error=derivative_error,
        exact_projected_force_identity_error=force_error,
        largest_expected_occupation_boundary_force=boundary_effect,
        independent_physical_continuum_error_bound_claimed=False)

def run():
    old=json.loads(OLD.read_text('utf-8-sig'))
    assert old['round']==931 and old['all_scientific_checks_passed']
    errors=[];rows=[];bad=[]
    k=.23
    for spec in old['species']:
        m=spec['mass'];charge=spec['charge']
        for boost in (0.,.2,.55):
            gam=1/np.sqrt(1-boost*boost)
            eb=np.sqrt(m*m+.07**2+.03**2+k*k/4)
            p=np.array([.07,.03,gam*(-k/2+boost*eb)])
            p1=np.array([.07,.03,gam*(k/2+boost*eb)])
            q,j,t,_,_=source(p,p1,m);ql=ETA@q
            ward_j=float(np.max(np.abs(np.einsum('a,aij->ij',ql,j))))
            ward_t=float(np.max(np.abs(np.einsum('a,abij->bij',ql,t))))
            for a,b in itertools.product(range(2),repeat=2):
                jj=charge*j[:,a,b];tt=t[:,:,a,b]
                c,r,z,ce,re,ze=contractions(q,tt,tt,jj,jj)
                errors.append(max(abs(c-r-z),abs(ce-re-ze),ward_j,ward_t))
                if a==b==0:
                    rows.append(dict(species=spec['label'],boost=boost,gravitational_covariant=real(c),
                         gravitational_radiative=real(r),gravitational_constraint=real(z),
                         electromagnetic_covariant=real(ce),electromagnetic_radiative=real(re),
                         electromagnetic_constraint=real(ze),source_Ward_error=max(ward_j,ward_t)))
                    wrong=tt.copy();wrong[0,0]+=.1
                    bc,br,bz,*_=contractions(q,wrong,wrong,jj,jj)
                    bad.append(abs(bc-br-bz))
    assert max(errors)<3e-12 and max(bad)>.01
    # Static rest density is carried entirely by constraints, not TT quanta.
    ts=np.zeros((4,4),complex);ts[0,0]=1
    js=np.array([1,0,0,0],complex)
    cg,rg,zg,ce,re,ze=contractions(np.array([0.,0.,0.,1.]),ts,ts,js,js)
    assert (real(cg),real(rg),real(zg),real(ce),real(re),real(ze))==(.5,0.,.5,-1.,0.,-1.)
    electron=next(s for s in old['species'] if s['label']=='e')
    actual=finite_hamiltonian(electron['mass'],electron['charge'])
    paths=[OLD,HERE.parent/'931/charged_matter_response.py',HERE.parent/'research_note_936.md',
           HERE.parent/'939/drafts/common_model_recovery_map.md',Path(__file__)]
    return dict(round=940,date='2026-10-07',all_scientific_checks_passed=True,
        inherited_mass_charge_table=old['species'],exact_rational_decomposition_error=exact_rational_check(),
        maximum_Dirac_kernel_or_Ward_error=float(max(errors)),
        unconserved_source_negative_control_gap=float(max(bad)),exchange_rows=rows,
        static_unit_source=dict(gravity_full=.5,gravity_TT_only=0.,electromagnetic_full=-1.,electromagnetic_transverse_only=0.),
        native_joint_diagnostic=actual,
        canonical_weak_field_methods_are_inherited_physics=True,
        finite_positive_process_and_own_sources_constructed=True,
        full_original_curved_branch_or_SM_recovery_proved=False,
        cognition_generates_Einstein_or_dimension_proved=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        old=json.loads(TARGET.read_text('utf-8'))
        def check(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a:check(a[k],b[k])
            elif isinstance(a,list):
                assert len(a)==len(b)
                for x,y in zip(a,b):check(x,y)
            elif isinstance(a,float):assert abs(a-b)<1e-11
            else:assert a==b,(a,b)
        check(result,old)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_hashes','exchange_rows','inherited_mass_charge_table')},ensure_ascii=False,indent=2))
