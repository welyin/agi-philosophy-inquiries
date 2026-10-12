"""Finite Fourier enclosure of the ORIGINAL 1129/1137 receiver tail.
Default read-only. Analytic midpoint, Fourier projection and cell errors are
separate from a deliberately generous IEEE binary64 arithmetic budget.
No continuum dynamics is simulated by the finite CAR diagnostic.
"""
from pathlib import Path
from decimal import Decimal as D, localcontext
import argparse, json, math, time
import numpy as np
PI=D('3.141592653589793238462643383279502884197169399375105820974944592307816406')
CUT=31
GRID=128
QUAD=16384
L=8.

def dec_cs(x):
    """Decimal Taylor; input reduced to [-pi,pi], tail < 1e-100."""
    x=x-(x/(2*PI)).to_integral_value()*(2*PI)
    c=t=D(1); s=v=x
    for j in range(1,85):
        t *= -x*x/D((2*j-1)*(2*j)); c += t
        v *= -x*x/D((2*j)*(2*j+1)); s += v
    return c,s

def coefficients():
    with localcontext() as ctx:
        ctx.prec=60
        total=[D(0) for _ in range(CUT+1)]; z=D(0)
        # even midpoint quadrature on [-1,1]; double positive half.
        for j in range(QUAD//2):
            x=D(2*j+1)/D(QUAD)
            q=(-D(1)/(1-x*x)).exp()
            z+=q*q
            c,_=dec_cs(PI*x/4)
            prev=D(1); cur=c
            total[0]+=q
            total[1]+=q*cur
            for n in range(2,CUT+1):
                prev,cur=cur,2*c*cur-prev
                total[n]+=q*cur
        scale=D(4)/QUAD; z*=scale
        raw=np.array([float(v*scale) for v in total])
        z0=float(z)
    # |q|<=1, |q'|<2, |q''|<100, hence |(q^2)''|<208.
    # Composite midpoint on length 2: bound B/(3 QUAD^2).
    dz=208/(3*QUAD**2)+1e-40
    ns=np.arange(-CUT,CUT+1); kk=np.pi*ns/4
    raw=raw[abs(ns)]
    di=(100+4*abs(kk)+kk*kk)/(3*QUAD**2)+1e-35
    center=raw/math.sqrt(L*z0)
    # interval division by sqrt(L*z), allowing either sign of raw.
    lo=np.minimum((raw-di)/math.sqrt(L*(z0-dz)),(raw-di)/math.sqrt(L*(z0+dz)))
    hi=np.maximum((raw+di)/math.sqrt(L*(z0-dz)),(raw+di)/math.sqrt(L*(z0+dz)))
    err=np.maximum(center-lo,hi-center)+1e-14
    return ns,kk,center,err,z0,dz

def phase_tables():
    with localcontext() as ctx:
        ctx.prec=60
        cos=[]; sin_over=[]
        for n in range(3*CUT**2+1):
            om=(1+PI*PI*D(n)/16).sqrt()
            c,s=dec_cs(5*om/4)
            cos.append(float(c)); sin_over.append(float(s/om))
        roots={}
        m=2
        while m<=GRID:
            vals=[]
            for j in range(m//2):
                c,s=dec_cs(2*PI*j/m)
                vals.append(complex(float(c),float(s)))
            roots[m]=np.array(vals)
            m*=2
    return np.array(cos),np.array(sin_over),roots

def inverse_fft(a,roots):
    """Explicit radix-2 butterflies, normalized inverse DFT, all axes.
    Unlike a vendor FFT, its operation count and roots are specified here.
    """
    bits=GRID.bit_length()-1
    rev=np.array([int(f'{n:0{bits}b}'[::-1],2) for n in range(GRID)])
    out=np.asarray(a,dtype=np.complex128)
    for axis in range(3):
        b=np.take(np.swapaxes(out,axis,-1),rev,axis=-1).copy()
        shape=b.shape
        m=2
        while m<=GRID:
            v=b.reshape(shape[:-1]+(-1,m))
            left=v[...,:m//2].copy()
            right=v[...,m//2:].copy()*roots[m]
            v[...,:m//2]=left+right
            v[...,m//2:]=left-right
            b=v.reshape(shape)
            m*=2
        out=np.swapaxes(b/GRID,axis,-1)
    return out

def run():
    ns,k,c,e,z,dz=coefficients()
    S=float(c@c); Slo=float(np.maximum(abs(c)-e,0)@np.maximum(abs(c)-e,0))
    assert math.sqrt(S)<1.001
    dcoeff=3.01*float(np.linalg.norm(e))
    dproj=math.sqrt(max(0,1-Slo**3))
    grad=math.sqrt(3*float((k*k*c)@c)*S*S)
    h=L/GRID
    dcell=h/math.sqrt(12)*grad
    # Explicit 21 butterfly stages, correctly rounded Decimal-derived roots;
    # IEEE binary64 arithmetic. 1e-6 exceeds the norm/positive-sum budget.
    darith=1e-6
    delta=dproj+dcoeff+dcell+darith
    phase_c,phase_d,roots=phase_tables()
    n2=ns[:,None,None]**2+ns[None,:,None]**2+ns[None,None,:]**2
    base=c[:,None,None]*c[None,:,None]*c[None,None,:]
    ds=base*phase_d[n2]
    small=[base*phase_c[n2],ds,1j*k[:,None,None]*ds,
           1j*k[None,:,None]*ds,1j*k[None,None,:]*ds]
    idx=ns%GRID; ix=np.ix_(idx,idx,idx)
    density=np.zeros((GRID,)*3)
    fft_res=0.
    for a in small:
        full=np.zeros((GRID,)*3,dtype=np.complex128); full[ix]=a
        val=inverse_fft(full,roots)*GRID**3/L**1.5
        # An independent implementation diagnostic, not the enclosure proof.
        ref=np.fft.ifftn(full)*GRID**3/L**1.5
        fft_res=max(fft_res,float(np.max(abs(val-ref))))
        density+=abs(val)**2
    x=np.fft.fftfreq(GRID)*L
    r=np.sqrt(x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2)
    def ramp(rad):
        v=np.clip((rad-1.9)/.1,0,1)
        return 35*v**4-84*v**5+70*v**6-20*v**7
    # Minimum radial mask on each closed cell, including periodic boundary cells.
    low=ramp(np.maximum(r-math.sqrt(3)*h/2,0))
    a_cell=float(np.sum(low*density)*h**3)
    b_cell=float(np.sum(low*low*density)*h**3)
    a_lower=max(0,math.sqrt(max(0,a_cell-1e-6))-delta)**2
    nu_lower=max(0,math.sqrt(max(0,b_cell-1e-6))-delta)**2
    mid=ramp(r)
    a_mid=float(np.sum(mid*density)*h**3)
    b_mid=float(np.sum(mid*mid*density)*h**3)
    # q=a^2/nu >= a since 0<=kappa<=1 and nu<=a.
    q_lower=a_lower
    assert q_lower>.035 and nu_lower>.03
    assert grad<3.1 and fft_res<1e-10
    # H^2 mode constants. Polynomial derivative bounds: 22 and 1e4.
    k1=(20+22)/math.sqrt(.03)
    lap=10000+2*22/1.9
    k2=(1800+2*22*20+lap)/math.sqrt(.03)
    assert k1<1000 and k2<1e6
    return dict(status='PASS',scope='ORIGINAL_MODE_FINITE_FOURIER_TAIL_ENCLOSURE',
        n_fourier=2*CUT+1,n_grid=GRID,n_midpoint=QUAD,
        z_midpoint=z,z_error=dz,one_dim_projected_norm_lower=Slo,
        projection_error=dproj,coefficient_error=dcoeff,cell_error=dcell,
        arithmetic_norm_budget=darith,complete_L2_error=delta,
        cell_weighted_a=a_cell,cell_weighted_nu=b_cell,
        q_lower=q_lower,nu_lower=nu_lower,q_diagnostic=a_mid*a_mid/b_mid,
        norm_diagnostic=float(np.sum(density)*h**3),gradient_norm=grad,
        independent_fft_max_error=fft_res,k1_upper=k1,k2_upper=k2,
        midpoint_enclosure_assumptions='Decimal 60 digits, exp/sqrt correctly rounded; Taylor remainder and rounding <1e-35; IEEE754 binary64 nearest arithmetic, no overflow',
        scientific_calibration_increment=0)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true');a=ap.parse_args()
    result=run(); target=Path(__file__).with_suffix('.json')
    if a.write_results:
        if target.exists(): raise FileExistsError(target)
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif target.exists():
        saved=json.loads(target.read_text(encoding='utf-8'))
        for k,v in result.items():
            if isinstance(v,float): assert math.isclose(v,saved[k],rel_tol=1e-10,abs_tol=1e-12),(k,v,saved[k])
            else: assert v==saved[k],k
    print(json.dumps(result,ensure_ascii=False,indent=2))
