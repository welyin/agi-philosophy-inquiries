"""901 working: integer-enclosed Fourier algebra for the original initial PDE.
No FFT is used to certify coefficients. FFT only chooses an exact dyadic witness.
"""
from pathlib import Path
from fractions import Fraction as Q
import math,json,sys,argparse,time
S=1<<80
KEEP=(12,12,6)  # third Fourier index represents physical wave number 2*k_z
DROP=1<<15
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];RESEARCH=HERE.parents[1]
WITNESS=HERE/'initial_psi_witness.json';TARGET=HERE/'wiener_certificate_results.json'

def ceilq(q):return -((-q.numerator)//q.denominator)
def bound(q):return {'exact':str(q),'value':float(q)}
class I:
    def __init__(self,a,b=None):self.a=Q(a);self.b=Q(a if b is None else b);assert self.a<=self.b
    def __add__(a,b):
        b=iv(b);return I(a.a+b.a,a.b+b.b)
    __radd__=__add__
    def __neg__(a):return I(-a.b,-a.a)
    def __sub__(a,b):return a+-iv(b)
    def __rsub__(a,b):return iv(b)+-a
    def __mul__(a,b):
        b=iv(b);v=[x*y for x in (a.a,a.b) for y in (b.a,b.b)];return I(min(v),max(v))
    __rmul__=__mul__
    def inv(a):assert a.a*a.b>0;return I(1/a.b,1/a.a)
    def __truediv__(a,b):return a*iv(b).inv()
    def __rtruediv__(a,b):return iv(b)*a.inv()
    def __pow__(a,n):
        if n<0:return a.inv()**-n
        r=I(1)
        for _ in range(n):r=r*a
        return r
    def sqrt(a):
        assert a.a>=0
        x=math.isqrt(a.a.numerator*S*S//a.a.denominator)
        y=math.isqrt(a.b.numerator*S*S//a.b.denominator)+1
        return I(Q(x,S),Q(y,S))
def iv(a):return a if isinstance(a,I) else I(a)
def pi_interval():
    def atan(n):
        v=sum((Q((-1)**j,(2*j+1)*n**(2*j+1)) for j in range(40)),Q(0))
        return I(v,v+Q(1,81*n**81))
    return 16*atan(5)-4*atan(239)

# Exact Kronecker convolution with balanced base digits, not floating FFT.
def integer_convolve(a,b):
    if not a or not b:return {}
    am=[min(k[i] for k in a) for i in range(3)];bm=[min(k[i] for k in b) for i in range(3)]
    ax=[max(k[i] for k in a) for i in range(3)];bx=[max(k[i] for k in b) for i in range(3)]
    dims=[ax[i]-am[i]+bx[i]-bm[i]+1 for i in range(3)]
    stride=(dims[1]*dims[2],dims[2],1)
    norm=lambda d:sum(abs(x)+abs(y) for x,y in d.values())
    bits=8*((max(1,(2*norm(a)*norm(b)).bit_length())+7)//8);width=bits//8;B=1<<bits
    def pack(d,offset,part):
        slots={sum((k[j]-offset[j])*stride[j] for j in range(3)): (v[0] if part==0 else v[1] if part==1 else v[0]+v[1]) for k,v in d.items()}
        length=max(slots)+1;data=bytearray((length+1)*width);carry=0
        for i in range(length):
            carry,digit=divmod(slots.get(i,0)+carry,B)
            data[i*width:(i+1)*width]=digit.to_bytes(width,'little')
        value=int.from_bytes(data[:length*width],'little')
        return value+carry*(1<<(bits*length))
    ar,ai,asum=[pack(a,am,j) for j in range(3)]
    br,bi,bsum=[pack(b,bm,j) for j in range(3)]
    real=ar*br-ai*bi;imag=asum*bsum-ar*br-ai*bi
    length=math.prod(dims);dataR=real.to_bytes((length+1)*width,'little',signed=True);dataI=imag.to_bytes((length+1)*width,'little',signed=True)
    cr=ci=0;out={}
    for index in range(length):
        r=int.from_bytes(dataR[index*width:(index+1)*width],'little')+cr
        im=int.from_bytes(dataI[index*width:(index+1)*width],'little')+ci
        cr=int(r>=B//2);ci=int(im>=B//2);r-=cr*B;im-=ci*B
        if r or im:
            x,rest=divmod(index,stride[0]);y,z=divmod(rest,stride[1]);k=(x+am[0]+bm[0],y+am[1]+bm[1],z+am[2]+bm[2])
            out[k]=(r,im)
    return out

class P:
    def __init__(self,c=None,e=0):self.c=c or {};self.e=int(e)
    @staticmethod
    def const(value):
        v=iv(value);center=((v.a+v.b)*S/2).__floor__();e=ceilq(max(abs(Q(center)-v.a*S),abs(v.b*S-Q(center))))
        return P({(0,0,0):(center,0)},e)
    @property
    def n(self):return sum(abs(a)+abs(b) for a,b in self.c.values())
    def clean(self):
        c={};e=self.e
        for k,(a,b) in self.c.items():
            if any(abs(k[j])>KEEP[j] for j in range(3)) or abs(a)+abs(b)<DROP:e+=abs(a)+abs(b)
            else:c[k]=(a,b)
        return P(c,e)
    def __add__(a,b):
        b=poly(b);c=dict(a.c)
        for k,(x,y) in b.c.items():
            u,v=c.get(k,(0,0));c[k]=(u+x,v+y)
        return P(c,a.e+b.e).clean()
    __radd__=__add__
    def __neg__(a):return P({k:(-x,-y) for k,(x,y) in a.c.items()},a.e)
    def __sub__(a,b):return a+-poly(b)
    def __rsub__(a,b):return poly(b)+-a
    def __mul__(a,b):
        b=poly(b);c=integer_convolve(a.c,b.c);e=ceilq(Q(a.n*b.e+b.n*a.e+a.e*b.e,S));out={}
        for k,(x,y) in c.items():
            qx,rx=divmod(x,S);qy,ry=divmod(y,S);out[k]=(qx,qy);e+=int(rx!=0)+int(ry!=0)
        return P(out,e).clean()
    __rmul__=__mul__
    def __truediv__(a,b):return a*P.const(iv(b).inv())
    def __pow__(a,n):
        assert n>=0
        ans=poly(1)
        for _ in range(n):ans=ans*a
        return ans
    def derivative(a,i):
        assert a.e==0,'Cannot differentiate an unweighted unknown error'
        out={}
        for k,(x,y) in a.c.items():
            f=k[i]*(2 if i==2 else 1);out[k]=(-f*y,f*x)
        return P(out)
    def inverse_power(a,p,terms=22):
        # center c is an exact rational. Includes uncertainty in input a.
        c=Q(a.c[(0,0,0)][0],S);assert c>0
        w=1-a/c;rho=Q(w.n+w.e,S);assert rho<1
        power=poly(1);out=poly(0)
        for n in range(terms+1):
            out=out+math.comb(p+n-1,n)*power
            if n<terms:power=power*w
        ratio=rho*Q(p+terms+1,terms+2)
        assert ratio<1
        tail=math.comb(p+terms,terms+1)*rho**(terms+1)/(1-ratio)
        out.e+=ceilq(tail*S)
        return out*P.const(c**-p)
    def stats(a):return {'terms':len(a.c),'center_W1':float(Q(a.n,S)),'radius_W1':float(Q(a.e,S))}
def poly(x):return x if isinstance(x,P) else P.const(x)
def trig(k,sine=False):
    assert len(k)==3
    return P({tuple(k):(0,-S//2) if sine else (S//2,0),tuple(-x for x in k):(0,S//2) if sine else (S//2,0)})
def evaluate(p,points):
    import numpy as np
    ans=np.zeros(points.shape[:-1],complex)
    for k,(r,i) in p.c.items():ans+=(float(Q(r,S))+1j*float(Q(i,S)))*np.exp(1j*(k[0]*points[...,0]+k[1]*points[...,1]+2*k[2]*points[...,2]))
    return ans.real

def parameters():
    raw=json.loads((RESEARCH/'archive_531_553/544/joint_singlet_common_mass_rg_results.json').read_text('utf-8'),parse_float=str)['examples'][2]
    z=raw['state'];lh,p,ls=map(I,(z['lambda_H'],z['p'],z['lambda_s']));cx=I(z['x'])/4;cy=I(z['y'])/4
    det=lh*ls-p*p;uh=(ls*cx-p*cy)/det;us=(lh*cy-p*cx)/det;hs=uh.sqrt();ss=us.sqrt()
    match=json.loads((RESEARCH/'archive_531_553/531/joint_gauge_matter_constraints_results.json').read_text('utf-8'),parse_float=str)['historical_running']['matched_inverse_couplings']
    pi=pi_interval();g=[(I(x)+Q(beta)*(-I(raw['u']))/(8*pi*pi)).inv() for x,beta in zip(match,('-41/6','19/6','7'))]
    gy,gw,gc=g;bc=gc/2;bw=gw/2;b0=gy/72;kc=gc.inv();kw=gw.inv();k0=36/gy
    lmax=(lh+ls+((lh-ls)**2+4*p*p).sqrt())/2
    fmin=2-((I('1.05')*hs)**2+(I('1.06')*ss)**2)/6
    db1=uh*Q('0.1025');db2=us*Q('0.1236');up=lmax*(db1*db1+db2*db2)/(4*fmin*fmin)
    return locals()

def shear_from_momentum(M):
    # Exact rational multiplier from 572. Zero mode is analytically zero.
    out=[];keys=set().union(*(p.c for p in M))
    for i in range(3):
        for j in range(3):
            c={};error=3*sum(p.e for p in M)
            for k in keys:
                w=(k[0],k[1],2*k[2]);k2=sum(x*x for x in w)
                if k2==0:continue
                vals=[p.c.get(k,(0,0)) for p in M]
                W=[]
                for a in range(3):
                    W.append(tuple(Q(vals[a][part],k2)-Q(w[a]*sum(w[b]*vals[b][part] for b in range(3)),4*k2*k2) for part in (0,1)))
                div=tuple(sum(w[b]*W[b][part] for b in range(3)) for part in (0,1))
                v=[w[i]*W[j][part]+w[j]*W[i][part]-(Q(2,3)*div[part] if i==j else 0) for part in (0,1)]
                r,im=-v[1],v[0];ir=r.__floor__();ii=im.__floor__();error+=int(r!=ir)+int(im!=ii);c[k]=(ir,ii)
            out.append(P(c,error).clean())
    return out

def source():
    z=parameters();lh,p,ls,uh,us,hs,ss,bc,bw,b0,kc,kw,k0,up=[z[k] for k in ('lh','p','ls','uh','us','hs','ss','bc','bw','b0','kc','kw','k0','up')]
    sx,cx=trig((1,0,0),True),trig((1,0,0));sy,cy=trig((0,1,0),True),trig((0,1,0));st,ct=trig((1,1,0),True),trig((1,1,0))
    # cos^2(z) is represented directly using the even-z index.
    cz2=(1+trig((0,0,1)))/2
    H=1+Q('.05')*st;V=1+Q('.06')*cy;h=poly(hs)*H;s=poly(ss)*V;h2=poly(uh)*H**2;s2=poly(us)*V**2
    dh=poly(hs)*Q('.05')*ct;ds=-poly(ss)*Q('.06')*sy
    f=1+Q('.38')*cx;ax=Q('.41')*f;ay=Q('.29')*(1+Q('.2')*sy)
    ex=Q('.17')*(f**2-(1+Q('1.5')*Q('.38')**2));r2=-Q('2')*Q('.17')*Q('.38')*f*sx;r3=-ax*ex
    e0=6*Q('.41')*Q('.17')*((2*Q('.38')-Q('.75')*Q('.38')**3)*sx+Q('.75')*Q('.38')**2*trig((2,0,0),True)+Q('.38')**3/12*trig((3,0,0),True))+Q('.13')
    f0=Q('.08')*cx+Q('.03')*sy
    ph=-poly(Q('.112')/hs)*ct;ps=Q('.025')*cx-poly(Q('.0056')/(Q('.06')*ss))*sy
    M=[-Q('.0056')*ct**2+Q('.18')*cy*r3+Q('.07')*cx*f0,
       -Q('.0056')*ct**2+ps*ds-ay*r2+Q('.48')*sx*r3-e0*f0,poly(0)]
    for m in M:assert abs(m.c.get((0,0,0),(0,0))[0])<=m.e
    shear=shear_from_momentum(M)
    FF=2-(h2+s2)/6;Fi=FF.inverse_power(1,18);Fii=FF.inverse_power(2,18)
    h2i=h2.inverse_power(1,22)
    pkp=FF*(ph**2+ps**2+4*(r2**2+r3**2)*h2i-(h*ph+s*ps)**2/12)
    A=sum((q*q for q in shear),poly(0))+pkp
    D2=2*dh**2+ds**2+h2*((ax**2+ay**2)/4+(Q('.09')*cy)**2+(Q('.24')*sx)**2+Q('.0036')*cz2)
    rad2=(h*dh)**2+(h*dh+s*ds)**2
    B=D2*Fi+rad2*Fii/6+Q('.0004')*cx**2
    deltaH=h2-poly(uh);deltaS=s2-poly(us)
    U=(poly(lh)*deltaH**2+2*poly(p)*deltaH*deltaS+poly(ls)*deltaS**2)*Fii/4
    C=2+2*poly(up)-2*U-Q('.0004')*sx**2
    Y=poly(bw)*(ex**2+Q('.0121')*cx**2+Q('.0081')*sy**2)+poly(b0)*(e0**2+Q('.0049')*cx**2)
    Y+=poly(kw)/2*((ax*ay)**2+Q('.0144')*cz2*(ax**2+ay**2))+poly(k0)/2*f0**2
    Scol=Q('.008868915');Y+=poly(bc)*(Q('.23')**2+Q('.19')**2+Q('.17')**2)+poly(kc)/2*Scol
    return dict(A=A,B=B,C=C,Y=Y,M=M,shear=shear,params=z)

def build_witness():
    import numpy as np
    sys.path.insert(0,str(ROOT/'scripts'));from research_layout import Layout,ResearchRuntime
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        q,p0,_,_,_=old.completed(24,1.)
        x=q['grid'][...,0];new=dict(q);new['B']=q['B']+Q('.0004').__float__()*np.cos(x)**2;new['C']=q['C']-(.02*np.sin(x))**2;new['U']=q['U']+.5*(.02*np.sin(x))**2
        psi,_,stats=old.geo.solve_hamiltonian(new,initial=p0)
        coeff=np.fft.fftn(psi)/24**3;waves=np.rint(old.geo.waves(24)).astype(int);data=[]
        for index in np.ndindex(coeff.shape):
            k=tuple(waves[index]);v=coeff[index]
            if max(abs(t) for t in k)>=12 or k[2]%2:continue
            if abs(v)<1e-14:continue
            data.append([int(k[0]),int(k[1]),int(k[2]//2),round(float(v.real)*S),round(float(v.imag)*S)])
        c={tuple(r[:3]):tuple(r[3:]) for r in data};sym={}
        for k in c:
            neg=tuple(-a for a in k);v=c[k];w=c.get(neg,(0,0));sym[k]=((v[0]+w[0])//2,(v[1]-w[1])//2)
        # Set opposite coefficients exactly, so the witness is rigorously real.
        for k in sorted(sym):
            if k>(0,0,0):sym[tuple(-a for a in k)]=(sym[k][0],-sym[k][1])
        sym[(0,0,0)]=(sym[(0,0,0)][0],0)
        data=[[*k,*v] for k,v in sorted(sym.items())]
        return dict(round=901,denominator=str(S),physical_z_index_multiplier=2,coefficients=data,choice='N24 old solver, small coefficient and Nyquist removal; exact dyadic polynomial thereafter',original_collocation_residual=stats['original_Hamiltonian_residual'])

def run():
    start=time.time();w=json.loads(WITNESS.read_text('utf-8'));psi=P({tuple(r[:3]):tuple(r[3:]) for r in w['coefficients']})
    assert all(psi.c.get(tuple(-a for a in k))==(v[0],-v[1]) for k,v in psi.c.items())
    center=Q(psi.c[(0,0,0)][0],S);variation=Q(psi.n-abs(psi.c[(0,0,0)][0]),S)
    assert center-variation>Q('.72') and center+variation<Q('1.125')
    print('witness',psi.stats(),flush=True)
    src=source();print('source constructed', {k:src[k].stats() for k in ('A','B','C','Y')},flush=True)
    lap=P({k:(-sum(x*x for x in (k[0],k[1],2*k[2]))*v[0],-sum(x*x for x in (k[0],k[1],2*k[2]))*v[1]) for k,v in psi.c.items()})
    psi5=psi**5;print('psi5',psi5.stats(),flush=True)
    inv7=psi.inverse_power(7,23);print('inverse7',inv7.stats(),flush=True)
    inv3=psi.inverse_power(3,21);print('inverse3',inv3.stats(),flush=True)
    residual=-8*lap+src['C']*psi5-src['B']*psi-src['A']*inv7-2*src['Y']*inv3
    r=Q(residual.n+residual.e,S)
    base=json.loads((HERE.parent/'900/continuous_reference_bounds_results.json').read_text('utf-8'))
    c0=Q(base['C0_error_per_continuum_residual']['exact']);c1=Q(base['each_C1_error_per_continuum_residual']['exact'])
    grad=[Q(sum(abs(k[j]*(2 if j==2 else 1))*(abs(v[0])+abs(v[1])) for k,v in psi.c.items()),S) for j in range(3)]
    merrors=[(Q(base['magnetic_gradient_error_residual_constant']['exact'])+Q(base['magnetic_gradient_error_residual_times_candidate_gradient']['exact'])*g)*r for g in grad]
    return dict(round=901,status='working',date='2026-10-06',formal_previous=900,cumulative_previous=3685,
        witness_terms=len(psi.c),witness_global_lower=bound(center-variation),witness_global_upper=bound(center+variation),
        source_balls={k:src[k].stats() for k in ('A','B','C','Y')},residual_ball=residual.stats(),
        continuum_residual_upper=bound(r),C0_error_upper=bound(c0*r),each_C1_error_upper=bound(c1*r),
        each_candidate_gradient_W1=[bound(g) for g in grad],each_initial_line_magnetic_gradient_error=[bound(m) for m in merrors],
        exact_integer_arithmetic_for_residual=True,FFT_only_selects_witness=True,full_time_evolution_computed=False,
        full_goal_completed=False,elapsed_seconds=time.time()-start)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--build-witness',action='store_true');a.add_argument('--write',action='store_true');arg=a.parse_args()
    if arg.build_witness:
        assert not WITNESS.exists();WITNESS.write_text(json.dumps(build_witness(),ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    result=run()
    if arg.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
