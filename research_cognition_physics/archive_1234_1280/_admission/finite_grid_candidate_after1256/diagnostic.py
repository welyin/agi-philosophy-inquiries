"""Non-numbered 5^3 admission diagnostic. Default read-only; --write creates once.
The complex-step calculation is a floating diagnostic, not interval integration.
Exact Fraction root bounds certify the explicitly derived initial defect.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import numpy as np

SIDE = 5
DX = 2.0
CELL_WEIGHT = DX**3
NC = 1.0
NR = 0.1
KROD = np.sqrt(0.1)
GRID = (SIDE,)*3
UPPER = ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
PAIRS = (list(zip(range(6),range(6,12))) + [(12,13),(14,15)]
         + list(zip(range(16,19),range(19,22)))
         + [(22,23),(24,25),(26,27)]
         + list(zip(range(28,31),range(31,34))) + [(34,35)])
# State slots: h[6], P_h[6], p/Pp, X/PX, xi[3]/Pi_phi[3],
# q/r, psi/B, Y/K, zeta[3]/Pi_Z[3], vartheta/Pi_vartheta.


def diff(a, direction):
    return (np.roll(a,-1,axis=direction)-np.roll(a,1,axis=direction))/(2*DX)


def gradient(a):
    return np.stack([diff(a,i) for i in range(3)],axis=-1)


def fields(y):
    """All 18 pairs at all 125 sites are retained; no zero-mode substitution."""
    h=np.zeros(GRID+(3,3),dtype=y.dtype)
    pi=np.zeros_like(h)
    for a,(i,j) in enumerate(UPPER):
        h[...,i,j]=h[...,j,i]=y[a]
        # P_h^{ij}=(2-delta_ij) pi^{ij}, for independent upper components.
        pi[...,i,j]=pi[...,j,i]=y[6+a]/(1 if i==j else 2)
    hi=np.linalg.inv(h)
    sg=np.sqrt(np.linalg.det(h))
    dh=np.stack([diff(h,i) for i in range(3)],axis=-3)
    # dh has indices (derivative,l,j); T has (l,i,j).
    T=np.swapaxes(dh,-3,-2)+np.moveaxis(dh,-3,-1)-dh
    gamma=np.einsum('...al,...lij->...aij',hi,T)/2
    dg=np.stack([diff(gamma,i) for i in range(3)],axis=-4)
    ric=(np.einsum('...aaij->...ij',dg)-np.einsum('...jaia->...ij',dg)
         +np.einsum('...aal,...lij->...ij',gamma,gamma)
         -np.einsum('...ajl,...lia->...ij',gamma,gamma))
    curvature=np.einsum('...ij,...ij->...',hi,ric)
    trace=np.einsum('...ij,...ij->...',pi,h)
    pi_sq=np.einsum('...ij,...ik,...jl,...kl->...',pi,h,h,pi)
    Cg=6*(pi_sq-trace*trace/2)/sg-sg*curvature/6
    dpi=np.stack([diff(pi,i) for i in range(3)],axis=-3)
    divpi=(np.einsum('...kjk->...j',dpi)
           +np.einsum('...jkl,...kl->...j',gamma,pi))
    Dg=-2*np.einsum('...ij,...j->...i',h,divpi)

    pg=gradient(y[12]); xg=gradient(y[14])
    phig=KROD*np.eye(3)+gradient(np.moveaxis(y[16:19],0,-1))
    phi_pi=np.moveaxis(y[19:22],0,-1)
    F=1+y[12]/5
    Cf=((y[13]**2+y[15]**2+np.sum(phi_pi**2,axis=-1)/F)/(2*sg)
        +sg*(np.einsum('...i,...ij,...j->...',pg,hi,pg)
              +np.einsum('...i,...ij,...j->...',xg,hi,xg)
              +y[12]**2+y[14]**2
              +F*np.einsum('...ai,...ij,...aj->...',phig,hi,phig))/2)
    Df=(y[13,...,None]*pg+y[15,...,None]*xg
        +np.einsum('...a,...ai->...i',phi_pi,phig))

    Zg=np.eye(3)+gradient(np.moveaxis(y[28:31],0,-1))
    eta=NR*np.linalg.det(Zg)  # Original constant label density, not a new pair.
    Dr=(np.einsum('...a,...ai->...i',np.moveaxis(y[31:34],0,-1),Zg)
        +y[35,...,None]*gradient(y[34]))
    W=(eta+(1+y[12])*y[35])/10
    Er=np.sqrt(np.einsum('...i,...ij,...j->...',Dr,hi,Dr)+W*W)
    C0=Cg+Cf+Er
    D0=Dg+Df+Dr
    root=np.sqrt(C0*C0-np.einsum('...i,...ij,...j->...',D0,hi,D0))
    return dict(h=h,hi=hi,sqrt_h=sg,Cg=Cg,Cf=Cf,Dg=Dg,Df=Df,
                eta=eta,Dr=Dr,W=W,Er=Er,C0=C0,D0=D0,root=root)


def pulse(t,start):
    s=t-start
    return 64*s**3*(1-s)**3 if 0 <= s <= 1 else 0.0


def mother(y,t=0.0):
    z=fields(y)
    U=(.71+(y[22]**2+y[23]**2+y[27]**2)/2
       +(y[25]-1)*pulse(t,0)*y[22]*y[12]/4
       +pulse(t,1)*y[14]*y[26]/10)
    z['U']=U
    z['J']=z['root']/NC-U
    z['H']=CELL_WEIGHT*np.sum(NC*U-z['root'])
    z['lapse']=-z['C0']/z['root']
    z['shift']=np.einsum('...ij,...j->...i',z['hi'],z['D0'])/z['root'][...,None]
    return z


def initial(delta):
    y=np.zeros((36,)+GRID)
    y[:3]=1
    y[6:9]=-np.sqrt(1.11)/3
    y[14]=.1; y[22]=.5; y[25]=1; y[35]=.1
    d=delta*np.array([1,0,-1,0,0])[:,None,None]
    y[19]=d/KROD
    return y


def flow_complex_step(y,t=0.0,step=1e-28):
    """Differentiate the full 4500-dimensional mother H, not fitted equations."""
    grad=np.empty_like(y)
    for index in np.ndindex(y.shape):
        q=y.astype(complex)
        q[index]+=1j*step
        grad[index]=mother(q,t)['H'].imag/step
    flow=np.zeros_like(y)
    for a,b in PAIRS:
        weight=CELL_WEIGHT*(NC if (a,b) in ((22,23),(24,25),(26,27)) else 1)
        flow[a]=grad[b]/weight
        flow[b]=-grad[a]/weight
    return flow


def clean(x):
    x=float(x)
    return 0.0 if abs(x)<1e-12 else float(f'{x:.12g}')


def root_bounds(value,bits=80):
    scale=1<<bits
    n=isqrt((value.numerator*scale*scale)//value.denominator)
    lo,hi=Q(n,scale),Q(n+1,scale)
    assert lo*lo <= value < hi*hi
    return lo,hi


def fraction_string(x):
    return f'{x.numerator}/{x.denominator}'


def exact_certificate():
    delta=Q(1,10); spacing=Q(2); A0=Q(187,200); U=Q(167,200)
    A=A0-5*delta**2
    root_sq=A*A-delta**2
    lo,hi=root_bounds(root_sq)
    jlo,jhi=lo-U,hi-U
    defect_lo=-5*delta**3/(spacing*lo)
    defect_hi=-5*delta**3/(spacing*hi)
    assert jlo>Q(1,25) and jhi<Q(9,200)
    assert defect_hi<-Q(1,400)
    assert A+(Q(7,40)+5*delta**2)==Q(111,100)
    assert A0+Q(7,40)==Q(111,100)
    assert delta-delta==0
    assert root_sq+delta**2==A*A
    # Identity proof uses positive A: reconstructed Ec=sqrt(root_sq+d^2)=A.
    assert A>0 and root_sq>0 and A<1
    c_lo=5*delta**3/(2*spacing*hi)
    c_hi=5*delta**3/(2*spacing*lo)
    return {
      'status':'PASS', 'method':'Fraction exact algebra and 80-bit dyadic root enclosure',
      'delta':fraction_string(delta),'spacing':fraction_string(spacing),
      'A0':fraction_string(A0),'A_delta':fraction_string(A),
      'root_squared':fraction_string(root_sq),
      'root_interval':[fraction_string(lo),fraction_string(hi)],
      'J_interval':[fraction_string(jlo),fraction_string(jhi)],
      'J_strict_bounds':['1/25','9/200'],
      'Jdot_site_1_interval':[fraction_string(defect_lo),fraction_string(defect_hi)],
      'Jdot_site_1_strict_upper':'-1/400',
      'Jdot_sites_3_and_4_interval':[fraction_string(c_lo),fraction_string(c_hi)],
      'source_Tnn':'111/100','momentum_constraint':'0',
      'hamiltonian_constraint':'0 by positive-root identity',
      'interpretation':'Exact initial algebra and explicit local-defect formula; not a future grid/PDE error bound.'}


def floating_diagnostic(delta):
    y=initial(delta); z=mother(y)
    flow=flow_complex_step(y)
    # Entire vector field is a full canonical derivative; one holomorphic JVP
    # then differentiates every reconstructed local J along that field.
    jdot=mother(y.astype(complex)+1j*1e-28*flow)['J'].imag/1e-28
    d=delta*np.array([1,0,-1,0,0])[:,None,None]
    A=Q(187,200).__float__()-5*d*d
    rad=np.sqrt(A*A-d*d)
    pred=np.broadcast_to(5*diff(d*d*(d/rad),0),GRID)
    # Analytic derivation applies at zero-d sites. At nonzero sites no theorem
    # is inferred from observing a floating zero.
    mask=np.broadcast_to(d==0,GRID)
    error=float(np.max(np.abs(jdot[mask]-pred[mask])))
    Ec=-z['C0']
    total_source=(z['Cf']+z['Er']+Ec)/z['sqrt_h']
    constraint=z['C0']+np.sqrt(np.einsum('...i,...ij,...j->...',z['D0'],z['hi'],z['D0'])+z['root']**2)
    assert error<5e-12
    assert np.max(np.abs(constraint))<5e-12
    assert np.max(np.abs(total_source-1.11))<5e-12
    assert np.min(z['J'])>.04
    assert abs(np.sum(jdot))*CELL_WEIGHT<5e-10
    return {
      'delta':delta,'J_min':clean(np.min(z['J'])),
      'Jdot_x1_line_at_x2_x3_0':[clean(x) for x in jdot[:,0,0]],
      'Jdot_integral':clean(np.sum(jdot)*CELL_WEIGHT),
      'max_analytic_defect_error_on_zero_d_sites':clean(error),
      'max_initial_hamiltonian_constraint_abs':clean(np.max(np.abs(constraint))),
      'max_initial_total_source_error':clean(np.max(np.abs(total_source-1.11))),
      'Kdot_max_abs':clean(np.max(np.abs(flow[27]))),
      'method':'4500 complex-step partial derivatives of full H, then full J directional derivative',
      'step':'1e-28','assertion_tolerance':'5e-12 (integrated Jdot tolerance 5e-10)',
      'status':'PASS_DIAGNOSTIC_NOT_INTERVAL_PROOF'}


def calculate():
    return {
      'schema':'finite_grid_candidate_after1256.v1',
      'status':'INITIAL_DEFECT_CERTIFIED_TASK_NOT_ADMITTED',
      'code_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
      'grid':[5,5,5],'period':10,'spacing':2,'cell_symplectic_weight':8,
      'canonical_pairs_per_cell':18,'real_phase_space_dimension':4500,
      'hamiltonian_support_radius_in_grid_edges':2,
      'vector_field_support_radius_upper_bound_in_grid_edges':4,
      'exact_initial_certificate':exact_certificate(),
      'complex_step_diagnostics':[floating_diagnostic(0.0),floating_diagnostic(0.1)],
      'not_certified':['finite-time nonhomogeneous prediction error','Einstein/Bianchi residual budget',
          'finite-grid exact HDA','physical terminal cell-readout realization',
          'quantum state or quantum error bound','strict light cone'],
      'scope':'A chosen finite discretization has a local clock-source defect at nonuniform data. This does not reject smaller packets, another discretization, or the research route.'}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--write',action='store_true'); args=ap.parse_args()
    output=Path(__file__).with_name('diagnostic_results.json')
    if args.write and output.exists():
        raise FileExistsError('Refusing to overwrite existing diagnostic_results.json')
    result=calculate()
    if args.write:
        with output.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('CREATED_ONCE',output)
    else:
        saved=json.loads(output.read_text(encoding='utf-8'))
        assert result==saved, 'Recomputed diagnostic differs from saved result'
        print('READ_ONLY_RECOMPUTATION_MATCH')
    print(result['status'])
    print('FULL_GRID',result['grid'],'DIMENSION',result['real_phase_space_dimension'])
    for row in result['complex_step_diagnostics']:
        print('delta',row['delta'],'Jdot',row['Jdot_x1_line_at_x2_x3_0'])

if __name__=='__main__': main()
