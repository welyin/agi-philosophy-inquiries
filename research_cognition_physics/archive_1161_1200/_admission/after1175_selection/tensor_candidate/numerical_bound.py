"""Whole-step residual certificate under a declared binary64 gamma_n model.

Fixed physics. --bench only times a few steps on the full grid. Default computes
the whole certificate and never writes a file. See numerical_bound.md.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Fr
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time
import numpy as np

HERE = Path(__file__).resolve().parent
UR = 2.0**-53
NODES = np.array([0., .125, .25, .5, .75, .875, 1.])


def up(x):
    return np.nextafter(x, np.inf)


def gamma(n):
    return float(up(n * UR / (1 - n * UR)))


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def interpolation_matrix():
    """Exact rational inverse of degree-six Bernstein sampling at dyadic nodes."""
    mat = [[Fr(math.comb(6,j))*Fr(float(x))**j*(1-Fr(float(x)))**(6-j)
            for j in range(7)] + [Fr(int(i==j)) for j in range(7)]
           for i,x in enumerate(NODES)]
    for i in range(7):
        pivot = next(j for j in range(i,7) if mat[j][i])
        mat[i],mat[pivot]=mat[pivot],mat[i]
        d=mat[i][i]
        mat[i]=[v/d for v in mat[i]]
        for j in range(7):
            if j!=i:
                d=mat[j][i]
                mat[j]=[a-d*b for a,b in zip(mat[j],mat[i])]
    exact=[r[7:] for r in mat]
    c=np.array([[float(v) for v in r] for r in exact])
    row_sum=max(sum(abs(v) for v in r) for r in exact)
    coefficient_error=max(sum(abs(Fr(float(v))-e) for v,e in zip(r,rr))
                          for r,rr in zip(c,exact))
    return c, float(up(float(row_sum))), float(up(float(coefficient_error)))


def hermite(y0,y1,f0,f1,dt):
    diff=y1-y0
    return np.stack([y0,dt*f0,3*diff-dt*(2*f0+f1),
                     -2*diff+dt*(f0+f1)])


def evaluate(c,x,dt):
    y=c[0]+x*(c[1]+x*(c[2]+x*c[3]))
    dy=(c[1]+x*(2*c[2]+3*x*c[3]))/dt
    return y,dy


def hermite_envelope(y0,y1,f0,f1,c,dt):
    """Uniform component errors relative to exact stored-data Hermite curve.

    Works separately on real/imaginary parts; complex conversion paid below.
    Conservative Y/F envelopes make cancellation harmless.
    """
    y=max(float(np.max(np.abs(y0))),float(np.max(np.abs(y1))))
    f=max(float(np.max(np.abs(f0))),float(np.max(np.abs(f1))))
    ce=float(up(gamma(20)*(6*y+4*dt*f)))
    cs=float(up(np.max(np.sum(np.abs(c),axis=0))*(1+gamma(4))+4*ce))
    ds=float(up(np.max(np.abs(c[1])+2*np.abs(c[2])+3*np.abs(c[3]))
                /dt*(1+gamma(8))+6*ce/dt))
    ye=float(up(4*ce+gamma(20)*cs))
    de=float(up(6*ce/dt+gamma(20)*ds))
    return cs,ds,ye,de


def rhs_background(q,g,k,m):
    n,v=m.coeffs(g,m.COUPLING)
    d=np.empty((len(k),4))
    d[:,0]=n*m.MS
    d[:,1:]=k*v
    dq=np.zeros_like(q)
    for j,(a,b) in enumerate(m.PAIRS):
        dq[:,a]-=2*q[:,4+j]*d[:,b]
        dq[:,b]+=2*q[:,4+j]*d[:,a]
        dq[:,4+j]=2*(d[:,b]*q[:,a]-d[:,a]*q[:,b])
    em,et=m.source_observables(q[:,:4])
    b=np.array([m.EPS*em+m.NU*sum(et),m.COUPLING*(et[0]-et[1])])
    dg=np.empty((2,3))
    dg[:,0]=-2*b*g[:,1]
    dg[:,1]=2*(b*g[:,0]+m.OMEGA*g[:,2])
    dg[:,2]=-2*m.OMEGA*g[:,1]
    return dq,dg


def input_enclosure(m,iv):
    # Mathematical pi and exact displayed decimal parameters, not libm trust.
    pi=(np.array(3.141592653589793),np.array(3.1415926535897936))
    angle=iv.divpos(iv.scale(iv.mul(pi,iv.pt(np.arange(33))),2),iv.pt(m.GRID))
    sins=iv.sin_interval(angle)
    ell=iv.divpos(iv.pt(m.LENGTH),iv.pt(m.GRID))
    kk=iv.divpos(sins,ell)
    kval=np.sin(2*np.pi*np.arange(33)/m.GRID)/m.ELL
    kerr=float(up(max(np.max(np.abs(kval-kk[0])),np.max(np.abs(kval-kk[1])))))
    kmax=float(up(max(abs(kk[0]).max(),abs(kk[1]).max())))
    # Relative errors in exact positive rational multiplicity/V weights.
    weight_relative=gamma(2)
    eps=iv.input_interval(m.EPS);nu=iv.input_interval(m.NU);gt=iv.input_interval(m.COUPLING)
    d_eps=max(abs(float(eps[0])-m.EPS),abs(float(eps[1])-m.EPS))
    d_nu=max(abs(float(nu[0])-m.NU),abs(float(nu[1])-m.NU))
    d_gt=max(abs(float(gt[0])-m.COUPLING),abs(float(gt[1])-m.COUPLING))
    ep,nv,gv=float(eps[1]),float(nu[1]),float(gt[1])
    # d=d0+zs fs+zt ft, bounds in coefficient l1 norm.
    f_s=ep*m.MS+3*nv*kmax
    f_t=2*gv*kmax
    de0=3*kerr
    des=m.MS*d_eps+3*(nv*kerr+kmax*d_nu)
    det=2*(gv*kerr+kmax*d_gt)
    fs2=iv.add(iv.mul(iv.scale(eps,m.MS),iv.scale(eps,m.MS)),
               iv.divpos(iv.scale(iv.mul(nu,nu),3),iv.scale(iv.mul(ell,ell),2)))
    ft2=iv.divpos(iv.mul(gt,gt),iv.mul(ell,ell))
    F=float(iv.sqrt_upper(max(float(fs2[1]),float(ft2[1]))))
    return dict(kerr=kerr,kmax=kmax,weight_relative=weight_relative,fs=f_s,ft=f_t,
                de0=de0,des=des,det=det,F=F,ell_lo=float(ell[0]),
                probe_k=float(kval[1]),probe_kerr=float(up(max(
                    abs(kval[1]-kk[0][1]),abs(kval[1]-kk[1][1])))),
                eps=ep,nu=nv,g=gv)


def bernstein_samples(samples,cmat):
    return (cmat @ samples.reshape(7,-1)).reshape(samples.shape)


def background_residual(qc,gc,qenv,genv,dt,m,inp,cmat,cnorm,cerr):
    qr=np.empty((7,)+qc.shape[1:])
    gr=np.empty((7,2,3))
    qmax,qdot,qerr,qderr=qenv
    gmax,gdot,gerr,gderr=genv
    for j,x in enumerate(NODES):
        q,dq=evaluate(qc,x,dt)
        g,dg=evaluate(gc,x,dt)
        fq,fg=rhs_background(q,g,m.KABS,m)
        qr[j]=dq-fq; gr[j]=dg-fg
    qb=bernstein_samples(qr,cmat)
    gb=bernstein_samples(gr,cmat)
    # Uniform ordinary arithmetic and input enclosure, component by component.
    fs,ft=inp['fs'],inp['ft']
    dsum=m.MS+3*inp['kmax']+(fs+ft)*gmax
    derror=inp['de0']+(inp['des']+inp['det'])*gmax+(fs+ft)*gerr
    q_rhs_error=2*dsum*qerr+2*derror*(qmax+qerr)+gamma(64)*2*dsum*qmax
    q_sample_error=qderr+q_rhs_error+gamma(2)*(qdot+2*dsum*qmax)
    reduce_gamma=gamma(len(m.WEIGHTS)+32)+inp['weight_relative']
    bmax=2*max(fs,ft)*qmax
    berror=2*max(fs,ft)*qerr+2*max(inp['des'],inp['det'])*(qmax+qerr)
    berror+=reduce_gamma*2*max(fs,ft)*qmax
    g_rhs_error=2*(bmax+m.OMEGA)*gerr+2*berror*(gmax+gerr)
    g_rhs_error+=gamma(64)*2*(bmax+m.OMEGA)*gmax
    g_sample_error=gderr+g_rhs_error+gamma(2)*(gdot+2*(bmax+m.OMEGA)*gmax)
    qpad=cnorm*q_sample_error+(cerr+gamma(8)*cnorm)*float(np.max(np.abs(qr)))
    gpad=cnorm*g_sample_error+(cerr+gamma(8)*cnorm)*float(np.max(np.abs(gr)))
    # Triangle bound for uniform errors in 10 q and 6 G components.
    extra=math.sqrt(20*qpad*qpad+6*gpad*gpad)
    sq=2*np.einsum('n,jnc,jnc->j',m.WEIGHTS,qb,qb)+np.sum(gb*gb,axis=(1,2))
    arithmetic_factor=1+gamma(12*len(m.WEIGHTS)+128)+inp['weight_relative']
    raw=float(np.sqrt(np.max(sq)*arithmetic_factor))
    residual=float(up((raw+extra)*(1+gamma(100))))
    return residual,raw,float(up(extra)),gc


def h_probe(g,m):
    n,v=m.coeffs(g,m.COUPLING)
    vv=np.array([v,v[[1,0,2]]])
    kin=np.einsum('ab,pb,bij->apij',vv,m.LMOM,m.ALPHA)
    return kin[:,:,None]+n*m.MP[None,None,:,None,None]*m.BETA


def probe_residual(uc,gc,uenv,genv,dt,m,inp,cmat,cnorm,cerr):
    values=np.empty((7,)+uc.shape[1:],complex)
    for j,x in enumerate(NODES):
        u,du=evaluate(uc,x,dt)
        g,_=evaluate(gc,x,dt)
        values[j]=du+1j*h_probe(g,m)@u
    cb=bernstein_samples(values,cmat)
    raw=float(np.max(np.sqrt(np.sum(abs(cb)**2,axis=(-2,-1)))))
    umax,udot,uerr,uderr=uenv
    gmax,_,gerr,_=genv
    p=abs(inp['probe_k'])+inp['probe_kerr']; mm=.5+2e-16
    hs=inp['eps']*mm+inp['nu']*p;ht=inp['g']*p
    hmax=mm+p+(hs+ht)*gmax
    # Exact masses and decimals vs stored binary constants, liberal input term.
    herr=(hs+ht)*gerr+inp['probe_kerr']*(1+(inp['nu']+inp['g'])*gmax)
    herr+=64*UR*(1+mm+p)*(1+gmax)
    # Complex multiplication: factor2 safely covers real/imag envelopes.
    node_error=2*uderr+8*(hmax*2*uerr+herr*(umax+2*uerr))
    node_error+=gamma(128)*8*hmax*umax+gamma(4)*(udot+4*hmax*umax)
    coeff_error=cnorm*node_error+(cerr+gamma(16)*cnorm)*float(np.max(abs(values)))
    # Each 4x4 block Frobenius <=4 times a complex-entry upper bound.
    return float(up((raw+4*coeff_error)*(1+gamma(100))))


def curve_integrals(gc,dt,m):
    """Exact polynomial integral of zs; convex-hull absolute integral bound of zt."""
    # Power->cubic Bernstein, exact rational fractions paid by tiny scalar pad.
    zt=gc[:,1,2]
    b=np.array([sum(math.comb(j,k)/math.comb(3,k)*zt[k] for k in range(j+1))
                for j in range(4)])
    zs_int=dt*sum(gc[j,0,2]/(j+1) for j in range(4))
    # Triangle in Bernstein form, not an assumption that zt has fixed sign.
    abszt_int=dt*sum(abs(b))/4
    # scalar arithmetic plus Hermite coefficient error paid in run().
    return float(zs_int),float(abszt_int)


def run(steps=2048,bench=0):
    assert steps in (1024,2048,4096), "Fixed dyadic integration choices; no physics scan."
    m=module('tensor_fixed_preview',HERE/'preview.py')
    old=HERE.parent.parent/'after1174_direction/mean_field_numerical_bound.py'
    iv=module('old_basic_interval_tools',old)
    cmat,cnorm,cerr=interpolation_matrix()
    inp=input_enclosure(m,iv)
    L=float(up(up(2*up(math.sqrt(2)))*inp['F']))
    dt=m.TIME/steps
    ex,ph1,ph2=iv.positive_factors(float(up(L*dt)))
    r=np.zeros((len(m.KABS),4));r[:,1]=1
    zero=m.KABS[:,0]==0;r[zero,1]=0;r[zero,0]=1
    s=np.zeros((len(m.KABS),6))
    g=np.zeros((2,3));g[:,0]=1
    u=np.broadcast_to(np.eye(4,dtype=complex),(2,4,3,4,4)).copy()
    state=(r,s,g,np.zeros(4),u)
    err=err_int=probe_error=0.
    maxr=maxraw=maxpad=0.
    zsint=ztabs=scalar_padding=0.
    begin=time.perf_counter()
    count=bench or steps
    for step in range(count):
        f0=m.derivative(state,m.COUPLING)
        f2=m.derivative(m.plus(state,f0,dt/2),m.COUPLING)
        f3=m.derivative(m.plus(state,f2,dt/2),m.COUPLING)
        f4=m.derivative(m.plus(state,f3,dt),m.COUPLING)
        new=tuple(a+dt*(b+2*c+2*d+e)/6 for a,b,c,d,e in zip(state,f0,f2,f3,f4))
        f1=m.derivative(new,m.COUPLING)
        q0=np.concatenate(state[:2],axis=-1);q1=np.concatenate(new[:2],axis=-1)
        dq0=np.concatenate(f0[:2],axis=-1);dq1=np.concatenate(f1[:2],axis=-1)
        qc=hermite(q0,q1,dq0,dq1,dt)
        gc=hermite(state[2],new[2],f0[2],f1[2],dt)
        uc=hermite(state[4],new[4],f0[4],f1[4],dt)
        qe=hermite_envelope(q0,q1,dq0,dq1,qc,dt)
        ge=hermite_envelope(state[2],new[2],f0[2],f1[2],gc,dt)
        ue=hermite_envelope(state[4],new[4],f0[4],f1[4],uc,dt)
        rr,raw,pad,_=background_residual(qc,gc,qe,ge,dt,m,inp,cmat,cnorm,cerr)
        pr=probe_residual(uc,gc,ue,ge,dt,m,inp,cmat,cnorm,cerr)
        err_int=float(up(err_int+up(err*up(dt*ph1))+up(rr*up(dt*dt*ph2))))
        err=float(up(up(err*ex)+up(rr*up(dt*ph1))))
        probe_error=float(up(probe_error+up(dt*pr)))
        zi,ti=curve_integrals(gc,dt,m)
        zsint+=zi;ztabs+=ti
        scalar_padding=float(up(scalar_padding+dt*(4*ge[2]+gamma(64)*(1+ge[0]))))
        maxr=max(maxr,rr);maxraw=max(maxraw,raw);maxpad=max(maxpad,pad)
        state=new
        if not bench and (step+1)%256==0:
            print(json.dumps({'progress_steps':step+1,'total_steps':steps,
                              'elapsed_s':round(time.perf_counter()-begin,2),
                              'max_residual':maxr,'integrated_error':err_int}),flush=True)
    elapsed=time.perf_counter()-begin
    if bench:
        return {'status':'BENCH_ONLY','steps_measured':bench,'full_steps':steps,
                'elapsed_s':elapsed,'projected_full_s':elapsed*steps/bench,
                'maximum_residual_first_steps':maxr,'maximum_roundoff_padding':maxpad,
                'interpolation_row_norm':cnorm,'input_k_error':inp['kerr']}
    pmax=float(up(2*np.pi/m.LENGTH))
    # Trace majorizes the 2x2 probe Gram eigenvalue; masses 0,.3,.5 all covered.
    fp=float(up(math.sqrt((inp['eps']*(.5+2e-16))**2
                         +(inp['nu']**2+inp['g']**2)*pmax**2)))
    bgerr=float(up(fp*err_int))
    rawu=state[4]
    gram=max(float(np.linalg.norm(a.conj().T@a-np.eye(4),'fro'))
             for a in rawu.reshape(-1,4,4))
    gram_upper=float(up(gram+gamma(256)*16*float(np.max(abs(rawu)))**2))
    polar_shift=float(up(gram_upper/(1+np.sqrt(1-gram_upper))))
    # Exact right polar factors are the mathematical normalized comparison.
    normu=np.empty_like(rawu)
    for ii in np.ndindex(rawu.shape[:-2]):
        a=rawu[ii]; ev,vec=np.linalg.eigh(a.conj().T@a)
        normu[ii]=a@(vec@np.diag(1/np.sqrt(ev))@vec.conj().T)
    arithmetic_probability_padding=1e-10
    eta_each=float(up(bgerr+probe_error+polar_shift+arithmetic_probability_padding))
    # Source symmetry means identical error budget for both branches, but not
    # cancellation of their numerical errors.
    pair=float(up(2*eta_each))
    positive_accum_padding=gamma(steps+64)*(m.TIME+abs(zsint)+ztabs)
    total_int_padding=float(up(scalar_padding+positive_accum_padding))
    Jup=float(up(inp['g']*(ztabs+total_int_padding+err_int)))
    Vup=float(up(m.TIME+inp['nu']*zsint+inp['g']*ztabs
                 +(inp['nu']+inp['g'])*total_int_padding
                 +math.sqrt(inp['nu']**2+inp['g']**2)*err_int))
    tauerr=float(up(inp['eps']*(err_int+total_int_padding)))
    result={
        'status':'RESIDUAL_CERTIFICATE_UNDER_DECLARED_GAMMA_N_MODEL',
        'steps':steps,'source_modes_represented':m.VOL,'compressed_modes':len(m.KABS),
        'T':m.TIME,'physics_unchanged':True,'Fop_upper':inp['F'],'L_upper':L,
        'initial_error':0.,'maximum_residual':maxr,'maximum_central_residual':maxraw,
        'maximum_roundoff_residual_padding':maxpad,'terminal_background_error':err,
        'integrated_background_error':err_int,'probe_background_error':bgerr,
        'probe_own_residual_error':probe_error,'probe_gram_upper':gram_upper,
        'right_polar_operator_correction':polar_shift,
        'terminal_probability_arithmetic_padding':arithmetic_probability_padding,
        'numerical_error_each_branch':eta_each,'numerical_error_pair':pair,
        'J_integral_abs_gzt_upper':Jup,'V_integral_vbar_plus_abs_gzt_upper':Vup,
        'clock_tau_error_from_Hermite_integral':tauerr,
        'zs_Hermite_integral_center':zsint,'abs_zt_Hermite_integral_upper_center':ztabs,
        'scalar_integral_arithmetic_padding':total_int_padding,
        'input_sine_absolute_error':inp['kerr'],'interpolation_row_norm':cnorm,
        'normalized_original_menu':m.readouts(normu),
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in [Path(__file__),HERE/'preview.py',old]},
        'scope':'Declared gamma_n basic-operation model, exact polynomial residual plus explicit input/roundoff envelope; not formal IEEE verification, full-K simulation, or certified hardware.'}
    print(json.dumps({'elapsed_s':elapsed,'numerical_pair':pair}),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--steps',type=int,default=2048)
    ap.add_argument('--bench',type=int,default=0);ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    result=run(args.steps,args.bench)
    if args.output:
        assert not args.bench
        with args.output.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
