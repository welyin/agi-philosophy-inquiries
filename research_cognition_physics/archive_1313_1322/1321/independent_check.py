"""Independent finite-mode/CAR review for 1321; imports no author code.

Full exponentially large Fock operators are defined by their actual CAR
monomials.  Low-occupation evaluations check the claimed invariant codes;
triangle estimates certify the entire finite Fock space.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
from fractions import Fraction
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE/'independent_results.json'


class Audit:
    def __init__(self): self.count=0; self.errors={}
    def close(self, name, got, want, tol=3e-10):
        g,w=np.asarray(got),np.asarray(want)
        e=float(np.linalg.norm((g-w).ravel())/max(1.,float(np.linalg.norm(w.ravel()))))
        self.count+=1; self.errors[name]=max(e,self.errors.get(name,0.))
        if not math.isfinite(e) or e>tol: raise AssertionError((name,e))
    def yes(self,name,test):
        self.count+=1; self.errors.setdefault(name,0.)
        if not bool(test): raise AssertionError(name)


def adj(x): return x.conj().T
def norm(x): return float(np.linalg.norm(x,2))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def unitary(h,t):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@adj(u)


def canonical(ops):
    """Products encountered here have distinct fermion modes, so no delta term."""
    if len({j for _,j in ops})!=len(ops):
        raise AssertionError('unexpected equal-mode CAR contraction')
    key=lambda x: (-x[0],x[1])
    inversions=sum(key(ops[i])>key(ops[j]) for i in range(len(ops)) for j in range(i+1,len(ops)))
    return (-1)**inversions,tuple(sorted(ops,key=key))


def dagger(ops):
    sign,out=canonical(tuple((1-cr,j) for cr,j in reversed(ops)))
    return sign,out


def apply_fermions(mask,ops):
    sign=1
    for cr,j in reversed(ops):
        occupied=(mask>>j)&1
        if occupied==cr: return None
        sign*=(-1)**((mask&((1<<j)-1)).bit_count())
        mask^=1<<j
    return mask,sign


def add(state,key,value):
    state[key]=state.get(key,0j)+value
    if abs(state[key])<1e-16: state.pop(key,None)


def distance(x,y): return math.sqrt(sum(abs(x.get(k,0)-y.get(k,0))**2 for k in x.keys()|y.keys()))
def state_norm(x): return math.sqrt(sum(abs(v)**2 for v in x.values()))
def scale(x,z): return {k:z*v for k,v in x.items()}


def field_modes(field,spin):
    p=-1 if spin==0 else 1
    return [(p,(0,4*field+spin)),(-p,(1,4*field+2+spin))]


def spin_integral(f,g):
    result={}
    for spin in (0,1):
        for p,op in field_modes(f,spin):
            for q,oq in field_modes(g,1-spin):
                if p+q: continue
                sg,ops=canonical((op,oq))
                result[ops]=result.get(ops,0)+(-1)**spin*sg
    return {k:v for k,v in result.items() if v}


def eps(i,j,k):
    if len({i,j,k})<3: return 0
    return (-1)**sum(x>y for x,y in ((i,j),(i,k),(j,k)))


def currents(family,color,kind):
    """Physical Eq(5) current before dividing its vacuum norm by 2."""
    result={}; offset=16*family
    if kind=='D':
        outer=[(i,j,a,eps(color,i,j)*(-1)**a/2)
               for i in range(3) for j in range(3) for a in (0,1) if eps(color,i,j)]
        for i,j,a,c in outer:
            for ops,z in spin_integral(offset+2*i+a,offset+2*j+1-a).items():
                result[ops]=result.get(ops,0)+c*z
    else:
        for a in (0,1):
            for ops,z in spin_integral(offset+2*color+a,offset+12+1-a).items():
                result[ops]=result.get(ops,0)+(-1)**a*z
    return {k:v for k,v in result.items() if v}


def pair_parts(family,color,kind):
    cur=currents(family,color,kind)
    return ({ops:c/2 for ops,c in cur.items() if sum(cr for cr,_ in ops)==0},
            {ops:c/2 for ops,c in cur.items() if sum(cr for cr,_ in ops)==2})


def pair_state(terms,do_dagger=False):
    state={}
    for ops,c in terms.items():
        if do_dagger:
            sg,ops=dagger(ops); c=sg*np.conj(c)
        applied=apply_fermions(0,ops)
        if applied is not None:
            mask,sg=applied; add(state,(0,mask),c*sg)
    return state


def portal(n):
    """(coefficient, scalar creation?, scalar label, actual CAR product)."""
    terms=[]
    for family in range(n):
        weight=4.**(-family)
        for color in range(3):
            aa,ab=pair_parts(family,color,'D')
            ba,bb=pair_parts(family,color,'L')
            for cr,sc,part in ((1,color+1,aa),(1,color+1,ab),(0,color+4,aa),(0,color+4,ab),
                               (0,color+1,ba),(0,color+1,bb),(1,color+4,ba),(1,color+4,bb)):
                for ops,c in part.items():
                    terms.append((weight*c,cr,sc,ops))
                    sg,ad=dagger(ops)
                    terms.append((weight*np.conj(c)*sg,1-cr,sc,ad))
    return terms


def shift(term):
    _,cr,_,ops=term
    return 2*(2*cr-1)+sum(2*creation-1 for creation,_ in ops)


def apply_terms(terms,state):
    out={}
    for c,cr,label,ops in terms:
        for (boson,mask),z in state.items():
            if (cr and boson!=0) or (not cr and boson!=label): continue
            applied=apply_fermions(mask,ops)
            if applied is None: continue
            new,sg=applied
            add(out,(label if cr else 0,new),c*z*sg)
    return out


def bright(n,color,kind):
    out={}; sn=sum(16.**(-i) for i in range(n))
    for i in range(n):
        p=pair_parts(i,color,'D' if kind=='A' else 'L')
        v=pair_state(p[0],True) if kind=='A' else pair_state(p[1])
        for key,z in v.items(): add(out,key,4.**(-i)*z/math.sqrt(sn))
    return out


def family_center(n,c=1.):
    dft=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
    u=np.eye(n,dtype=complex);u[:3,:3]=dft
    for j in range(2,n-1):
        r=np.eye(n);t=.01
        r[j,j]=r[j+1,j+1]=math.cos(t);r[j,j+1]=-math.sin(t);r[j+1,j]=math.sin(t)
        u=u@r
    return [np.diag(np.arange(1,2*n,2)).astype(complex),
            np.diag(np.arange(2,2*n+1,2)).astype(complex),
            np.diag(np.arange(1,n+1)).astype(complex),
            c*(u@np.diag(np.arange(1,n+1))@adj(u))]


def cp_value(me,mnu):
    he=me@adj(me);hn=mnu@adj(mnu);comm=he@hn-hn@he
    return float(np.trace(comm@comm@comm).imag)


def mass_vacuum(n,matrices):
    out={}
    for sector,mat in enumerate(matrices):
        pairs=[(2*col,6+col) for col in range(3)] if sector==0 else (
              [(2*col+1,9+col) for col in range(3)] if sector==1 else
              ([(13,14)] if sector==2 else [(12,15)]))
        for left,right in pairs:
            for i in range(n):
                for j in range(n):
                    for ops,co in spin_integral(16*i+left,16*j+right).items():
                        coef=co*np.conj(mat[i,j]) # all-left coefficient = physical M*
                        for c,opp in ((coef,ops),):
                            a=apply_fermions(0,opp)
                            if a is not None: add(out,a[0],c*a[1])
                        sg,dg=dagger(ops);a=apply_fermions(0,dg)
                        if a is not None: add(out,a[0],sg*np.conj(coef)*a[1])
    return out


def actual_mode_channel(n,u,reverse=False):
    """Kraus operators from occupation-basis partial trace, no postselection.

    Domain is vacuum plus the specified three-dimensional input code.  All
    unused receiver occupation states stay in the output Hilbert space.
    """
    masks=[]
    for i in range(n):
        fields=range(6) if reverse else list(range(6))+[12,13]
        for field in fields:
            for spin in (0,1):
                masks.append(64*i+4*field+spin+(0 if reverse else 2))
    receiving=sum(1<<j for j in masks)
    entries=[]
    for logical in range(4):
        if logical==0: state={(0,0):1.}
        else:
            color=logical-1;state={}
            column=2 if reverse else 0
            for factor,component in ((u[0,column],bright(n,color,'A')),
                                     (u[1,column],{(color+1,0):1.}),
                                     (u[2,column],bright(n,color,'C'))):
                for key,z in component.items():add(state,key,factor*z)
        for (boson,mask),z in state.items():
            entries.append(((boson,mask&~receiving),mask&receiving,logical,z))
    basis=sorted({rec for _,rec,_,_ in entries});index={m:i for i,m in enumerate(basis)}
    kraus={}
    for env,rec,logical,z in entries:
        if env not in kraus:kraus[env]=np.zeros((len(basis),4),complex)
        kraus[env][index[rec],logical]+=z
    code=np.zeros((len(basis),4),complex);code[index[0],0]=1.
    for color in range(3):
        for (_,mask),z in bright(n,color,'A' if reverse else 'C').items():code[index[mask],color+1]=z
    return list(kraus.values()),code


def run():
    a=Audit();rng=np.random.default_rng(1321019)
    # Actual Fourier quadrature of the selected spinor waves, including every
    # particle/antiparticle cross term; eight points integrate all frequencies.
    z=2*np.pi*np.arange(8)/8
    for spin in (0,1):
        for p,_ in field_modes(0,spin):
            for q,_ in field_modes(1,1-spin):
                a.close('spatial momentum integral',np.exp(1j*(p+q)*z).mean(),float(p+q==0))
    for family in (0,2,4):
        for color in range(3):
            for kind in ('D','L'):
                cur=currents(family,color,kind)
                a.yes('current has all eight pair monomials',len(cur)==8)
                a.close('current physical coefficient magnitude',list(map(abs,cur.values())),np.ones(8))
                aa,bb=pair_parts(family,color,kind)
                a.close('pair code CAR normalization',state_norm(pair_state(aa,True)),1.)
                a.close('antipair code CAR normalization',state_norm(pair_state(bb)),1.)
    # Explicit nonresonant pair creation, conjugate chain, and the whole
    # K=2 invariant input code, evaluated in the actual CAR algebra.
    families=[];vac={(0,0):1.};s3=sum(Fraction(1,16)**i for i in range(3))
    for n in (3,4,5):
        v=portal(n);v0=[t for t in v if shift(t)==0];vp=[t for t in v if shift(t)==4];vm=[t for t in v if shift(t)==-4]
        a.yes('exhaustive resonant frequencies',len(v)==len(v0)+len(vp)+len(vm))
        fsum=sum(Fraction(1,4)**i for i in range(n));sn=sum(Fraction(1,16)**i for i in range(n))
        b0=sum(abs(t[0]) for t in v0);bp=sum(abs(t[0]) for t in vp)
        a.close('global monomial norm bound resonant',b0,48*float(fsum))
        a.close('global monomial norm bound raising',bp,24*float(fsum))
        a.yes('five generation uniform full Fock bounds',b0<64 and bp<32 and len(v)>0)
        a.close('resonant vacuum dark',state_norm(apply_terms(v0,vac)),0.)
        raised=apply_terms(vp,vac)
        a.yes('counter rotating vacuum process retained',state_norm(raised)>0)
        a.close('scalar cutoff makes raising square zero',state_norm(apply_terms(vp,raised)),0.)
        for color in range(3):
            va=bright(n,color,'A');vc=bright(n,color,'C');vs={(color+1,0):1.}
            a.close('full flavor bright A normalized',state_norm(va),1.)
            a.close('full flavor bright C normalized',state_norm(vc),1.)
            a.close('actual CAR A to scalar',distance(apply_terms(v0,va),scale(vs,math.sqrt(float(sn)))),0.)
            a.close('actual CAR C to scalar',distance(apply_terms(v0,vc),scale(vs,math.sqrt(float(sn)))),0.)
            target=scale(va,math.sqrt(float(sn)))
            for key,val in scale(vc,math.sqrt(float(sn))).items():add(target,key,val)
            a.close('actual CAR scalar to both whole codes',distance(apply_terms(v0,vs),target),0.)
            # Conjugate t chain is retained, orthogonal to the intended input.
            vt={(color+4,0):1.};out=apply_terms(v0,vt)
            a.close('conjugate scalar chain norm',state_norm(out),math.sqrt(2*float(sn)))
            a.yes('conjugate chain disjoint',not(out.keys() & (va.keys()|vc.keys())))
        theta=math.pi*math.sqrt(float(sn/s3));success=math.sin(theta/2)**4
        h=np.array([[0,1,0],[1,0,1],[0,1,0]],complex)*math.sqrt(float(sn))
        time=math.pi/math.sqrt(2*float(s3));u=unitary(h,time)
        a.close('actual chain success at common time',abs(u[2,0])**2,success)
        a.close('nonpeak scalar leakage kept',abs(u[1,0])**2,math.sin(theta)**2/2)
        # Directly trace occupation states: the nonzero scalar branch for n=4,5
        # must be discarded with its actual environment, not called empty.
        ku,wc=actual_mode_channel(n,u)
        kr,wa=actual_mode_channel(n,u,True)
        a.close('occupation trace first channel is TP',sum(adj(k)@k for k in ku),np.eye(4))
        a.close('occupation trace second channel is TP',sum(adj(k)@k for k in kr),np.eye(4))
        inp=np.zeros(12,complex)
        for k in range(3):inp[3*(k+1)+k]=1/np.sqrt(3)
        rhoin=np.outer(inp,inp.conj())
        first=sum(np.kron(k,np.eye(3))@rhoin@adj(np.kron(k,np.eye(3))) for k in ku)
        decoder=np.kron(adj(wc),np.eye(3));reduced=decoder@first@adj(decoder)
        a.close('all retained receiver state in declared code and vacuum',np.trace(reduced),1.)
        second=sum(np.kron(k,np.eye(3))@reduced@adj(np.kron(k,np.eye(3))) for k in kr)
        logical=np.kron(adj(wa),np.eye(3))@second@np.kron(wa,np.eye(3))
        expected=success**2*rhoin
        expected[:3,:3]+=(1-success**2)*np.eye(3)/3
        a.close('actual two occupation traces preserve partner and all failure',logical,expected)
        # Full erased-code Choi state, preserving an untouched dimension-three partner.
        psi=np.eye(3).reshape(-1)/np.sqrt(3);rho=np.outer(psi,psi)
        output=np.zeros((12,12),complex);output[:9,:9]=success**2*rho
        output[9:,9:]=(1-success**2)*np.eye(3)/3
        ideal=np.zeros_like(output);ideal[:9,:9]=rho
        a.close('complete two stage partner erasure',np.linalg.svd(output-ideal,compute_uv=False).sum()/2,1-success**2)
        counter=Fraction(18496,2**24);mass=Fraction(196608,2**40)
        total=1-success**2+float(counter+mass)
        a.yes('one uniform actual roundtrip budget',total<.01)
        centers=[]
        for c in (1.,2.):
            ms=family_center(n,c);center_vac=mass_vacuum(n,ms)
            nmodes=[64*i+60+j for i in range(n) for j in range(4)]
            nmask=sum(1<<j for j in nmodes)
            actual=sum(abs(z0)**2 for mask,z0 in center_vac.items() if mask&nmask)
            expected=4*float(np.linalg.norm(ms[3],'fro')**2)
            a.close('actual right neutral vacuum pair norm',actual,expected)
            totalvac=4*sum(weight*np.linalg.norm(m,'fro')**2 for weight,m in zip((3,3,1,1),ms))
            a.close('all four sectors retained vacuum norm',sum(abs(x)**2 for x in center_vac.values()),totalvac)
            cp=cp_value(ms[2],ms[3]);a.yes('center CP invariant nonzero',abs(cp)>1)
            for m in ms:
                sg=np.linalg.svd(m,compute_uv=False)
                a.yes('entire radius one eighth ball full rank simple',sg.min()-.125>=.875-1e-12 and min(np.diff(sorted(sg)))>.25)
                a.yes('entire radius one eighth ball uniform norm',sg.max()+.125<11)
            # Samples verify formulae, not the analytic entire-ball claim.
            random=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));random*=.12/norm(random)
            for m in ms:
                sg=np.linalg.svd(m+random,compute_uv=False)
                a.yes('noncentral full matrices retained',sg.min()>.875 and sg.max()<11)
            centers.append(dict(multiplier=c,scaled_neutral_pair_norm_squared=actual,dimensionless_cp_invariant=cp))
        # p=1: the tiny probability is stable in float64, but no practical detector claim.
        mu=2.**-64;tau=mu/2**24;W=323.
        bounds=[]
        for item in centers:
            ay=mu*math.sqrt(item['scaled_neutral_pair_norm_squared']);bt=tau*tau*W*W/2
            bounds.append([(max(0,tau*ay-bt))**2,(tau*ay+bt)**2])
        a.yes('same short read time separates mass strengths',bounds[0][1]<bounds[1][0])
        families.append(dict(n=n,fermionic_modes=64*n,scalar_dimension=7,full_dimension_decimal=str(7*(1<<(64*n))),
                             flavor_sum=float(fsum),squared_flavor_sum=float(sn),resonant_error=1-success**2,
                             uniform_actual_roundtrip_upper=total,scaled_mass_centers=centers,neutral_probability_intervals=bounds))
    # Independent integration-by-parts checks in generic integer-gap systems.
    average_ratios=[]
    for size in (4,7,11):
        levels=rng.integers(0,5,size=size);h0=np.diag(levels.astype(float))
        x=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size));vv=(x+adj(x))/2;vv/=norm(vv)
        delta=levels[:,None]-levels[None,:];vr=np.where(delta==0,vv,0)
        primitive=np.zeros_like(vv);np.divide(vv,delta,out=primitive,where=delta!=0)
        a.close('primitive commutator exact',h0@primitive-primitive@h0,vv-vr)
        a.close('primitive anti Hermitian',primitive+adj(primitive),np.zeros_like(vv))
        for g in (.1,.03,.01):
            time=2./g;err=norm(unitary(h0+g*vv,time)-unitary(h0+g*vr,time))
            bound=g*norm(primitive)*(2+g*time*(norm(vv)+norm(vr)))
            a.yes('long time averaging uniform bound',err<=bound+1e-12);average_ratios.append(err/bound)
    # Wrongly deleting a degenerate resonant spectator never gains a g/gap error.
    h0=np.diag([0.,1.,1.]);v=np.zeros((3,3));v[1,2]=v[2,1]=1
    a.close('gap does not remove omitted resonant term',norm(unitary(h0+.001*v,math.pi/.002)-unitary(h0,math.pi/.002)),math.sqrt(2))
    # Physical normalization and independent source parameters before resonance.
    for p in (.7,1.,2.3):
        mass=2*p;volume=(2*math.pi/p)**3;kappa=2/math.sqrt(2*mass*volume)
        a.close('fixed spatial overlap dictionary',kappa,p/(2*math.pi)**1.5)
        lam=2.**-24*(2*math.pi)**1.5
        a.close('same local portal gives fixed g over p',kappa*lam/p,2.**-24)
        step=1e-5
        df=(2/math.sqrt(2*(mass+step)*volume)-2/math.sqrt(2*(mass-step)*volume))/(2*step)
        a.close('scalar mass overlap source',df,-kappa/(2*mass))
        dv=(2/math.sqrt(2*mass*volume*(1+step))-2/math.sqrt(2*mass*volume*(1-step)))/(2*step)
        a.close('volume overlap source',dv,-kappa/2)
    # The actual scalar total-occupation cutoff preserves its U(3) representation.
    x=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3));u,_=np.linalg.qr(x)
    rep=np.zeros((7,7),complex);rep[0,0]=1;rep[1:4,1:4]=u;rep[4:,4:]=u.conj()
    nb=np.diag([0,1,1,1,1,1,1]);a.close('total boson cutoff invariant',rep@nb@adj(rep),nb)
    return dict(status='PASS_WITH_SCOPE',round=1321,assertions=a.count,categories=len(a.errors),
                max_normalized_error=max(a.errors.values()),errors=a.errors,families=families,
                averaging_max_observed_ratio=max(average_ratios),
                exact_uniform_bounds=dict(vres=64,vraise=32,vfull=128,mass=17600,loose_mass=32768,
                                          hfull_over_p=323,counter_two_segment=float(Fraction(18496,2**24)),
                                          mass_two_segment=float(Fraction(196608,2**40))),
                scope=['complete finite cutoff Fock, not continuum or all gauge backgrounds',
                       'one bright flavor preparation, full internal three dimensional quantum code',
                       'mass ball full rank and simple; CP nonzero only established at centers and some neighborhoods',
                       'neutral occupation reads total mass strength, not individual masses or mixing phase',
                       'bare work vacuum is not interacting or massive vacuum',
                       'no certification of old P long time instrument or omitted spatial modes'])


def compare(x,y,path='result'):
    if isinstance(x,dict):
        if not isinstance(y,dict) or x.keys()!=y.keys():raise AssertionError(path+' fields')
        for k in x:compare(x[k],y[k],path+'.'+k)
    elif isinstance(x,list):
        if not isinstance(y,list) or len(x)!=len(y):raise AssertionError(path+' sequence')
        for i,(a,b) in enumerate(zip(x,y)):compare(a,b,f'{path}[{i}]')
    elif isinstance(x,bool):
        if type(y)is not bool or x!=y:raise AssertionError(path)
    elif isinstance(x,int):
        if type(y)is not int or x!=y:raise AssertionError(path)
    elif isinstance(x,float):
        if not isinstance(y,(int,float)) or isinstance(y,bool) or not np.isclose(x,y,rtol=2e-9,atol=0.,equal_nan=False):
            if not (abs(x)<1e-9 and abs(y)<1e-9 and 'errors.' in path):raise AssertionError((path,x,y))
    elif type(x)is not type(y) or x!=y:raise AssertionError((path,x,y))


def main():
    parser=argparse.ArgumentParser();g=parser.add_mutually_exclusive_group()
    g.add_argument('--check',action='store_true');g.add_argument('--write-new',action='store_true');args=parser.parse_args()
    result=run();result['code_sha256']=sha(__file__)
    result['admission_sha256']=sha(HERE.parent/'_admission'/'round1321_admission.md')
    if args.write_new:
        with OUT.open('x',encoding='utf-8')as stream:json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')
    else:compare(result,json.loads(OUT.read_text(encoding='utf-8')))
    print(json.dumps({k:result[k]for k in ('status','assertions','categories','max_normalized_error')},ensure_ascii=False))


if __name__=='__main__':main()
