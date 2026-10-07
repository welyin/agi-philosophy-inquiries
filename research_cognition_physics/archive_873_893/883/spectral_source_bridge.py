"""883: spectral hopping with original64 matter, joint source/history test.
A new long-range finite-graph candidate is declared. No full interacting
Gauss thermal continuum, curved target PDE, or autonomous detector is solved.
"""
from pathlib import Path
import argparse,json,sys,functools
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'spectral_source_bridge_results.json'
ETA=.1;BETA=2.;XI=.2;TIMES=(.17,.43)

def load():
    with ResearchRuntime(Layout()).installed():
        sys.path.insert(0,str(HERE.parent/'805'))
        import original_past_covariance_action as old
    return old

@functools.lru_cache(None)
def coefficients(N):
    assert N%2==1
    r=np.arange(-(N//2),N//2+1)
    a=2*np.pi/N
    coef=np.exp(-1j*a*np.outer(r,r))@r/N
    return r,coef

def polynomial(N,p,derivative=0):
    r,c=coefficients(N);a=2*np.pi/N
    return np.real(np.exp(1j*a*np.outer(np.atleast_1d(p),r))@(c*(1j*a*r)**derivative))

def smooth_step(x):
    x=np.asarray(x,float)
    out=np.zeros_like(x);out[x>=1]=1
    m=(x>0)&(x<1)
    s=x[m]
    # Stable flat C-infinity step.
    z=1/s-1/(1-s)
    out[m]=np.exp(-np.logaddexp(0,z))
    return out

def buffered(N,theta):
    a=2*np.pi/N
    theta=np.angle(np.exp(1j*np.asarray(theta)))
    return theta/a*smooth_step((np.pi-np.abs(theta))/(a*ETA))

def charges(old):
    matter=old.vertex.original.old.matter
    charge=np.zeros(32)
    for name,sl in matter.SLICES.items():
        charge[sl]=dict(Q=1,u=4,d=-2,L=-3,e=-6,nu=0)[name]
    return np.r_[charge,-charge],matter

def H(old,q,N,k,alpha,xi=XI,kappa=0.,kind='buffer'):
    p=np.array(k,float)
    x=p[0]+alpha*q
    if kind=='polynomial':dx=polynomial(N,x)
    elif kind=='buffer':dx=buffered(N,2*np.pi/N*x)
    elif kind=='continuum':dx=x
    else:raise ValueError(kind)
    kinetic=old.GAMMA[0]@np.diag(dx)+old.GAMMA[1]*p[1]+old.GAMMA[2]*p[2]
    return np.exp(-kappa)*kinetic+np.exp(xi)*old.MASS

def thermal(h):
    e,v=np.linalg.eigh(h)
    return (v*(.5+.5*np.tanh(BETA*e/2)))@v.conj().T

def unitary(h,t):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T

def wick(word,c):
    if not word:return 1.+0j
    ans=0j
    for j in range(1,len(word)):
        ans+=(-1)**(j-1)*c[word[0],word[j]]*wick(word[1:j]+word[j+1:],c)
    return ans

def history(old,q,N,alpha,kind):
    # Initial preparation xi fixed, no rethermalization during alpha source.
    p=thermal(np.exp(XI)*old.MASS)
    h=H(old,q,N,[0,0,0],alpha,kind=kind)
    fa=np.zeros(32,complex);fa[28]=1
    fb=np.zeros(32,complex);fb[28]=1/np.sqrt(2);fb[29]=1j/np.sqrt(2)
    rows=[]
    for f,t in ((fa,TIMES[0]),(fb,TIMES[1])):
        u=unitary(h,t)
        rows.extend([np.r_[f.conj(),f]@u,np.r_[-1j*f.conj(),1j*f]@u])
    rows=np.array(rows);c=rows@p@rows.conj().T
    out=[]
    for sa in (1,-1):
        for sb in (1,-1):
            aa=[(.5,()),(.5j*sa,(0,1))]
            bb=[(.5,()),(.5j*sb,(2,3))]
            value=sum(x*y*z*wick(wx+wy+wz,c) for x,wx in aa for y,wy in bb for z,wz in aa)
            out.append(value)
    out=np.array(out)
    assert np.max(abs(out.imag))<1e-11
    assert abs(out.sum()-1)<1e-11 and out.real.min()>0
    return out.real

def hermitian_exp(x):
    e,v=np.linalg.eigh(x);return (v*np.exp(1j*e))@v.conj().T

def matrix_function(N,S):
    ev,v=np.linalg.eig(S)
    f=v@np.diag(buffered(N,np.angle(ev)))@np.linalg.inv(v)
    assert np.max(abs(f-f.conj().T))<3e-12
    return (f+f.conj().T)/2

def gauge_check(old,q,matter):
    N=5;a=2*np.pi/N
    link=np.diag(np.exp(1j*a*.02*q[:32]))
    shift=np.zeros((N*32,N*32),complex)
    G=np.zeros_like(shift);rng=np.random.default_rng(883)
    for n in range(N):
        shift[n*32:(n+1)*32,((n+1)%N)*32:((n+1)%N+1)*32]=link
        group=[]
        for d in (3,2):
            x=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d))
            x=(x+x.conj().T)/2;x-=np.trace(x)/d*np.eye(d)
            group.append(hermitian_exp(.3*x))
        R=matter.representation(group[0],group[1],np.exp(.2j*rng.normal()))
        G[n*32:(n+1)*32,n*32:(n+1)*32]=R
    f=matrix_function(N,shift)
    fg=matrix_function(N,G@shift@G.conj().T)
    residual=float(np.max(abs(fg-G@f@G.conj().T)))
    assert residual<3e-11
    distant=float(np.linalg.norm(f[:32,2*32:3*32]))
    assert distant>1e-3
    return dict(nodes=N,original_particle_components=32,full_original_group=True,
                covariance_error=residual,non_neighbor_block_norm=distant,
                link_configuration='Flat hypercharge holonomy and its non-Abelian local gauge transform.')

def run():
    old=load();q,matter=charges(old)
    raw=[]
    for N in (5,7,9,11,17,19,33,35,65,67,129,131,257,259):
        k=np.arange(-N//2+1,N//2+1)
        interpolation=float(np.max(abs(polynomial(N,k)-k)))
        slope=float(polynomial(N,[0],1)[0])
        limiting=1-(-1)**(N//2)*np.pi/2
        raw.append(dict(N=N,free_symbol_error=interpolation,zero_mode_slope=slope,
            parity_subsequence_limit=limiting,limit_error=abs(slope-limiting),
            zero_mode_full64_gauge_source_error=6*abs(slope-1)))
        assert interpolation<3e-10
    matrix_rows=[]
    for N in (5,9,17):
        err=0.;reality=0.
        for k in ([0,0,0],[1,-1,0],[N//2,-N//2,N//2]):
            for alpha in (-.04,0.,.04):
                b=H(old,q,N,k,alpha,kappa=.13)
                c=H(old,q,N,k,alpha,kappa=.13,kind='continuum')
                err=max(err,float(np.max(abs(b-c))))
                hm=H(old,q,N,-np.array(k),alpha,kappa=.13)
                reality=max(reality,float(np.max(abs(old.CHARGE@b.conj()@old.CHARGE+hm))))
        assert max(err,reality)<4e-12
        matrix_rows.append(dict(N=N,full64_matrix_error=err,Nambu_reality_error=reality))
    h=2e-4
    ref0=history(old,q,5,0.,'continuum')
    refp=history(old,q,5,h,'continuum');refm=history(old,q,5,-h,'continuum')
    d1=(refp-refm)/(2*h);d2=(refp-2*ref0+refm)/h**2
    assert np.max(abs(d1))>1e-4
    histories=[]
    for N in (5,7,17,19,65,67):
        p0=history(old,q,N,0.,'polynomial')
        pp=history(old,q,N,h,'polynomial');pm=history(old,q,N,-h,'polynomial')
        slope=float(polynomial(N,[0],1)[0])
        p1=(pp-pm)/(2*h);p2=(pp-2*p0+pm)/h**2
        pointwise=max(float(np.max(abs(history(old,q,N,a,'buffer')-history(old,q,N,a,'continuum'))))
                      for a in (-.04,0.,.04))
        first_chain=float(np.max(abs(p1-slope*d1)))
        second_chain=float(np.max(abs(p2-slope*slope*d2)))
        assert pointwise<3e-13 and first_chain<2e-6 and second_chain<2e-5
        histories.append(dict(N=N,zero_source_probability_error=float(np.max(abs(p0-ref0))),
            first_source_error=float(np.max(abs(p1-d1))),
            second_source_error=float(np.max(abs(p2-d2))),
            analytic_first_chain_residual=first_chain,analytic_second_chain_residual=second_chain,
            buffered_full_source_path_history_error=pointwise))
    tails=[]
    for L in (1,2,4,8):
        qq=.3;d=2*qq**(2*(L+1))/(1+qq**2)
        tails.append(dict(cube_halfwidth=L,original_smooth_packet_squared_tail=3*d-3*d*d+d**3))
    boundary=[]
    for N in (9,17,33):
        k=N//2
        val=float(buffered(N,[2*np.pi/N*(k+.5)])[0])
        assert abs(val)<1e-12
        boundary.append(dict(N=N,shifted_momentum=k+.5,buffered_symbol_at_cut=val,
                             continuum_symbol=k+.5))
    return dict(round=883,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3668,
        argument_scope='A newly declared long-range covariant finite-graph hopping matches the original full homogeneous fermion preparation, histories and source jets in a flat holonomy patch. Minimal-path spectral gauging fails these source jets despite exact zero-source spectra. This is not the original full Q/E bridge or an interacting chiral continuum theorem.',
        original_Nambu_components=64,original_CAR_components=32,
        original_integer_hypercharges=sorted(set(q[:32].tolist())),
        smooth_buffer_eta=ETA,flat_holonomy_source_interval=[-.04,.04],
        polynomial_source_failure=raw,buffered_original_matrix_checks=matrix_rows,
        actual_ordered_history=dict(momentum=[0,0,0],times=list(TIMES),beta=BETA,preparation_log_mass=XI,
            record_modes='original right electron spin-z and spin-y wavefunctions, original full mass propagation',
            continuum_probabilities=ref0.tolist(),continuum_first_gauge_source=d1.tolist(),
            continuum_second_gauge_source=d2.tolist(),checks=histories),
        original_group_covariance=gauge_check(old,q,matter),
        smooth_menu_tail_examples=tails,branch_boundary_failure=boundary,
        finite_graph_hopping_is_new_declared_choice=True,
        exact_finite_graph_locality=False,uniform_dynamical_gauge_continuum_proved=False,
        finite_homogeneous_CAR_menu_source_bridge=True,
        full_interacting_Gauss_Gibbs_state_matched=False,original_Q_continuous_E_bridge_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
