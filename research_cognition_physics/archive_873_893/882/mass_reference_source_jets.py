"""882: original64 mass family, global covariance obstruction and finite CAR jets.
No future curved PDE, full interacting Gibbs state, or field-theory sum is
numerically solved. Infinite-volume claims are NOT made: the torus is compact.
"""
from pathlib import Path
from math import comb
import argparse,functools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'mass_reference_source_jets_results.json'
BETA=2.;PARAM=.2;Q=.3;TIMES=(.12,.31)
def const(x):return np.array([np.asarray(x),np.zeros_like(x),np.zeros_like(x)])
def mul(x,y):
    return np.array([x[0]*y[0],x[1]*y[0]+x[0]*y[1],
                     x[2]*y[0]+2*x[1]*y[1]+x[0]*y[2]])
def fun(x,f,d,dd):
    return np.array([f(x[0]),d(x[0])*x[1],d(x[0])*x[2]+dd(x[0])*x[1]**2])
def sqrt(x):return fun(x,np.sqrt,lambda v:.5/np.sqrt(v),lambda v:-.25/v**1.5)
def inv(x):return fun(x,lambda v:1/v,lambda v:-1/v**2,lambda v:2/v**3)
def tanh(x):
    return fun(x,np.tanh,lambda v:1-np.tanh(v)**2,
               lambda v:-2*np.tanh(v)*(1-np.tanh(v)**2))
def sin(x):return fun(x,np.sin,np.cos,lambda v:-np.sin(v))
def cos(x):return fun(x,np.cos,lambda v:-np.sin(v),lambda v:-np.cos(v))
def mass_load():
    with ResearchRuntime(Layout()).installed():
        sys.path.insert(0,str(HERE.parent/'805'))
        import original_past_covariance_action as old
    return old
@functools.lru_cache(None)
def weights(L):
    # Radius multiplicities and smooth packet weights, exact integer geometry.
    axis=np.arange(L+1);sgn=np.where(axis==0,1.,2.)
    xy=(axis[:,None]**2+axis[None,:]**2).ravel()
    ways=(sgn[:,None]*sgn[None,:]).ravel()
    packet=(Q**(2*(axis[:,None]+axis[None,:]))*sgn[:,None]*sgn[None,:]).ravel()
    count=np.zeros(3*L*L+1);w=np.zeros_like(count)
    for k in axis:
        rad=xy+k*k
        count+=np.bincount(rad,weights=ways*sgn[k],minlength=len(w))
        w+=np.bincount(rad,weights=packet*sgn[k]*Q**(2*k),minlength=len(w))
    assert abs(count.sum()-(2*L+1)**3)<1e-8
    d=2*Q**(2*(L+1))/(1+Q*Q)
    tail=3*d-3*d*d+d**3
    w/=w.sum()
    return count,w,float(tail)
def scalar_jets(ev,r2,lam=PARAM):
    a=np.exp(lam)
    q0=r2[:,None]+a*a*ev[None,:]
    q1=np.broadcast_to(2*a*a*ev[None,:],q0.shape)
    q2=2*q1
    e=sqrt(np.array([q0,q1,q2]))
    f=mul(tanh(BETA/2*e),inv(e))
    af=a*np.array([f[0],f[0]+f[1],f[0]+2*f[1]+f[2]])
    return e,f,af
def covariance(old,k,lam=PARAM):
    r2=np.dot(k,k);_,f,af=scalar_jets(old.EV,np.array([r2]),lam)
    v=old.VEC;g=sum((old.GAMMA[j]*k[j] for j in range(3)),np.zeros_like(old.MASS))
    out=[]
    for j in range(3):
        ff=(v*f[j,0])@v.conj().T;mm=(v*af[j,0])@v.conj().T
        out.append(.5*(g@ff+old.MASS@mm+(old.I if j==0 else 0)))
    return np.array(out)
def hilbert_schmidt(old):
    rows=[]
    for L in (4,8,16,32,64):
        count,_,_=weights(L);r2=np.arange(len(count))
        _,f,af=scalar_jets(old.EV,r2,PARAM)
        _,f0,af0=scalar_jets(old.EV,r2,0.)
        vals=.25*np.sum(r2[:,None]*(f[0]-f0[0])**2+old.EV[None,:]*(af[0]-af0[0])**2,axis=1)
        d1=.25*np.sum(r2[:,None]*f[1]**2+old.EV[None,:]*af[1]**2,axis=1)
        lead=(np.exp(PARAM)-1)**2*float(old.EV.sum())/4
        asym=lead*np.dot(count[1:],1/r2[1:])
        rows.append(dict(cube_halfwidth=L,HS_squared_difference=float(count@vals),
                         HS_squared_first_derivative=float(count@d1),
                         leading_asymptotic_sum=float(asym),
                         ratio_to_leading_sum=float(count@vals/asym)))
    high=[]
    leading=(np.exp(PARAM)-1)**2*float(old.EV.sum())/4
    for r in (10,100,1000):
        p=covariance(old,np.array([r,0.,0.]),PARAM)[0]
        p0=covariance(old,np.array([r,0.,0.]),0.)[0]
        scaled=r*r*float(np.linalg.norm(p-p0,'fro')**2)
        high.append(dict(radius=r,scaled_fiber_HS_squared=scaled,leading_coefficient=leading))
    return rows,high
# Full CAR matrices, used only to reconstruct the exact two-mode reduced state.
small=np.array([[0,1],[0,0]],complex);Z=np.diag([1.,-1.])
ca=np.kron(small,np.eye(2));cb=np.kron(Z,small)
MAJ=[ca+ca.conj().T,-1j*(ca-ca.conj().T),cb+cb.conj().T,-1j*(cb-cb.conj().T)]
T=np.array([[1,0,1,0],[-1j,0,1j,0],[0,1,0,1],[0,-1j,0,1j]])
IND=[30,31,62,63]
def contraction(old,L,dt,lam=PARAM):
    _,weight,_=weights(L);r2=np.arange(len(weight))
    e,_,_=scalar_jets(old.EV,r2,lam)
    th=tanh(BETA/2*e);co=cos(dt*e);si=sin(dt*e)
    ff=co-1j*mul(th,si)
    gg=mul(mul(th,co)-1j*si,inv(e))
    a=np.exp(lam)
    ag=a*np.array([gg[0],gg[0]+gg[1],gg[0]+2*gg[1]+gg[2]])
    out=[]
    for j in range(3):
        fm=weight@ff[j];gm=weight@ag[j]
        # The odd Gamma.k terms cancel in the symmetric smearing sum;
        # full mass-square mixing is retained before selecting the modes.
        nambu=.5*((old.VEC*fm)@old.VEC.conj().T
                  +old.MASS@(old.VEC*gm)@old.VEC.conj().T)
        block=nambu[np.ix_(IND,IND)]
        out.append(T@block@T.conj().T)
    return np.array(out)
def wick_jet(word,pairs):
    if not word:return np.array([1.,0.,0.],complex)
    if len(word)%2:return np.zeros(3,complex)
    total=np.zeros(3,complex)
    for j in range(1,len(word)):
        a,b=word[0],word[j]
        contraction=pairs[round(a[0]-b[0],10)][:,a[1],b[1]]
        total+=(-1)**(j-1)*mul(contraction,wick_jet(word[1:j]+word[j+1:],pairs))
    return total
def reduced_state(pairs):
    out=np.zeros((3,4,4),complex)
    for mask in range(16):
        idx=[j for j in range(4) if mask>>j&1]
        m=np.eye(4,dtype=complex)
        for j in idx:m=m@MAJ[j]
        w=wick_jet(tuple((0.,j) for j in idx),pairs)
        out+=w.conj()[:,None,None]*m[None,:,:]/4
    return out
def polynomial_product(a,b):
    return [(x*y,w+v) for x,w in a for y,v in b]
def hist_source(old,L,lam=PARAM):
    deltas=sorted({round(a-b,10) for a in (0.,)+TIMES for b in (0.,)+TIMES})
    pairs={dt:contraction(old,L,dt,lam) for dt in deltas}
    out=[]
    for signA in (1,-1):
        for signB in (1,-1):
            aa=[(.5,()),(.5j*signA,((TIMES[0],0),(TIMES[0],1)))]
            bb=[(.5,()),(.5j*signB,((TIMES[1],2),(TIMES[1],3)))]
            p=polynomial_product(polynomial_product(aa,bb),aa)
            out.append(sum((coef*wick_jet(word,pairs) for coef,word in p),np.zeros(3,complex)))
    rho=reduced_state(pairs)
    return np.array(out),rho
def derivative_bounds(old):
    # Uniform bound on lambda in [-.2,.2], k arbitrary. Kept deliberately loose.
    amin=np.exp(-.2);amax=np.exp(.2);mm=np.sqrt(old.EV.max());gap=amin*np.sqrt(old.EV.min())
    b=BETA/2;qp=2*amax**2*mm**2;qpp=2*qp
    fprime=b/(2*gap**2)+1/(2*gap**3)
    efprime=b/(2*gap)+1/(2*gap**2)
    efsecond=b*b/(2*gap**2)+3*b/(4*gap**3)+3/(4*gap**4)
    p1=.5*(amax*mm/gap+qp*efprime)
    p2=.5*(amax*mm/gap+2*amax*mm*qp*fprime+qp*qp*efsecond+qpp*efprime)
    duration=max(TIMES)
    u1=duration*amax*mm;u2=u1+u1*u1
    return np.array([1.,u1+p1,u2+2*u1*p1+p2])
def history_tail_bound(tail,bounds):
    c=np.array([1.,2*bounds[1],2*bounds[2]])
    delta=4*tail*bounds
    estimates=[np.zeros(3)]
    power=np.array([1.,0.,0.]);pairings=1
    for m in range(1,4):
        pairings*=2*m-1
        estimates.append(pairings*m*mul(delta,power))
        power=mul(power,c)
    return np.max(estimates,axis=0)
def run():
    old=mass_load()
    hs,high=hilbert_schmidt(old)
    assert all(b['HS_squared_difference']>a['HS_squared_difference'] for a,b in zip(hs,hs[1:]))
    assert abs(high[-1]['scaled_fiber_HS_squared']/high[-1]['leading_coefficient']-1)<2e-6
    matrix_fd=[]
    for h in (2e-4,1e-4):
        err=[0.,0.]
        for k in ([0,0,0],[1,-2,3]):
            k=np.array(k,float);a=covariance(old,k,PARAM)
            # Independently diagonalize the original64 H for the differences.
            def exact(lam):
                H=sum((old.GAMMA[j]*k[j] for j in range(3)),np.exp(lam)*old.MASS.copy())
                e,v=np.linalg.eigh(H);return (v*(.5+.5*np.tanh(BETA*e/2)))@v.conj().T
            pp=exact(PARAM+h);pm=exact(PARAM-h);p0=exact(PARAM)
            err[0]=max(err[0],float(np.max(np.abs((pp-pm)/(2*h)-a[1]))))
            err[1]=max(err[1],float(np.max(np.abs((pp-2*p0+pm)/(h*h)-a[2]))))
        matrix_fd.append(dict(step=h,first_error=err[0],second_error=err[1]))
    assert matrix_fd[-1]['first_error']<1e-7 and matrix_fd[-1]['second_error']<1e-5
    reference,rhoref=hist_source(old,24)
    bounds=derivative_bounds(old);ref_tail=weights(24)[2]
    rows=[]
    for L in (0,1,2,4,8):
        actual,rho=hist_source(old,L)
        assert np.max(np.abs(actual.sum(axis=0)-np.array([1,0,0])))<1e-11
        assert actual[:,0].real.min()>0 and np.max(np.abs(actual.imag))<1e-11
        assert np.linalg.eigvalsh(rho[0]).min()>0
        assert np.max(np.abs(np.trace(rho,axis1=1,axis2=2)-np.array([1,0,0])))<1e-12
        tail=weights(L)[2]
        errors=np.max(np.abs(actual-reference),axis=0)
        certificate=history_tail_bound(tail,bounds)+history_tail_bound(ref_tail,bounds)
        assert np.all(errors<certificate+5e-13)
        rows.append(dict(cube_halfwidth=L,packet_tail=tail,
                         joint_history_source_errors=errors.tolist(),
                         two_mode_state_source_trace_errors=[float(np.linalg.svd(rho[j]-rhoref[j],compute_uv=False).sum()) for j in range(3)],
                         conservative_analytic_history_tail_bounds=certificate.tolist()))
    finite_differences=[]
    for h in (.0008,.0004):
        pp=hist_source(old,8,PARAM+h)[0][:,0]
        pm=hist_source(old,8,PARAM-h)[0][:,0]
        base=hist_source(old,8,PARAM)[0]
        finite_differences.append(dict(step=h,first_error=float(np.max(np.abs((pp-pm)/(2*h)-base[:,1]))),
                       second_error=float(np.max(np.abs((pp-2*base[:,0]+pm)/(h*h)-base[:,2])))))
    assert finite_differences[-1]['second_error']<2e-7
    assert np.max(np.abs(reference[:,1]))>1e-5 and np.max(np.abs(reference[:,2]))>1e-5
    return dict(round=882,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3667,
        argument_scope='The original uniformly mass-scaled auxiliary-past thermal references are not normal in one fixed vacuum folium under the fixed equal-time CAR identification. Nevertheless finite smooth CAR preparations and actual free ordered histories have common C2 source jets. Neither statement identifies the full interacting Q/E models.',
        original_Nambu_components=64,original_CAR_components=32,
        mass_family='All original Yukawa and singlet-Majorana matrices scaled by exp(lambda), with fixed original scalar background.',
        beta=BETA,lambda_value=PARAM,original_mass_square_trace=float(old.EV.sum()),
        high_frequency_fiber_checks=high,global_covariance_cutoff_checks=hs,
        covariance_jet_independent_spectral_checks=matrix_fd,
        finite_packet='Normalized Fourier amplitudes proportional to 0.3^(|kx|+|ky|+|kz|), original indices30,31 plus conjugates.',
        actual_original_free_times=list(TIMES),
        full_mode_propagation_preserved=True,closed_two_mode_Hamiltonian_used=False,
        source_derivative_order=[0,1,2],
        joint_history_values_and_derivatives=reference.real.tolist(),
        finite_menu_cutoff_checks=rows,history_independent_difference_checks=finite_differences,
        reference_L24_analytic_history_tail_bounds=history_tail_bound(ref_tail,bounds).tolist(),
        fixed_global_normal_state_family_proved=False,
        finite_CAR_preparation_and_history_C2_proved=True,
        smooth_background_family_843_invalidated=False,
        all_reference_transport_bridges_excluded=False,
        full_interacting_Q_E_identified=False,curved_future_PDE_numerically_solved=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
