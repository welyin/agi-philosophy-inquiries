"""754: all-color source completion on the original 753 constrained background.

Finite-strength numerical sources are explicitly prescribed diagnostics.
The conditional quantum-source application is analytic in note754.
"""
import argparse,hashlib,json
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import numpy as np
import joint_reference_constraint_strata as bg
import joint_source_constraint_response as old
geo=bg.geo
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_irreducible_source_completion_results.json'

def maxabs(a):return float(np.max(abs(a)))
def components(a):return np.array([bg.inner(t,a) for t in bg.T])
def connection_matrices():
    c=bg.color(1.)
    C=np.array([np.column_stack([components(bg.cross(ai,t)) for t in bg.T]) for ai in c['A']])
    return c,C

def exact_gap_certificate():
    # 2*T8 is unnormalized here, so its squared basis norm is three.
    G=2*bg.T.copy();G[7]=np.diag([1.,1.,-2.])
    assert np.all(G.real==np.round(G.real)) and np.all(G.imag==np.round(G.imag))
    V=np.array([[1j*(G[k]@g-g@G[k]) for g in G] for k in (0,1,3)])
    ns=(31,27,21);mass=[[Q(0) for j in range(8)] for i in range(8)]
    for i in range(8):
        for j in range(8):
            for a,n in enumerate(ns):
                tr=2*np.trace(V[a,i].conj().T@V[a,j])
                assert tr.imag==0 and tr.real==round(tr.real)
                mass[i][j]+=Q(n*n*int(tr.real),160000)
    # Congruence of L(0)-.01 I in this exact rational basis.
    for i in range(8):mass[i][i]-=Q(3 if i==7 else 1,100)
    lower=[[Q(int(i==j)) for j in range(8)] for i in range(8)];pivots=[]
    for j in range(8):
        pivot=mass[j][j]-sum(lower[j][k]**2*pivots[k] for k in range(j))
        assert pivot>0;pivots.append(pivot)
        for i in range(j+1,8):
            lower[i][j]=(mass[i][j]-sum(lower[i][k]*lower[j][k]*pivots[k] for k in range(j)))/pivot
    assert Q(31**2+27**2+21**2,10000)<Q(1,4)
    return dict(zero_mode_gap_strict_lower='1/100',exact_LDL_pivots=[str(p) for p in pivots],
                connection_squared_norm_bound='2131/10000',nonzero_integer_mode_gap_lower='1/4',
                all_fourier_modes_covered_analytically=True)

@lru_cache(maxsize=3)
def color_operator(N):
    c,C=connection_matrices();k=geo.waves(N);I=np.eye(8)
    symbol=np.zeros(k.shape[:-1]+(8,8),complex)
    for i in range(3):
        D=1j*k[...,i,None,None]*I+C[i]
        symbol-=D@D
    assert maxabs(symbol-symbol.swapaxes(-1,-2).conj())<1e-12
    return c,C,symbol

def D(v,i,C):
    return geo.derivative(v,i)+np.einsum('ab,...b->...a',C[i],v)

def color_inverse(sigma):
    c,C,symbol=color_operator(sigma.shape[0])
    xi=geo.ifft(np.linalg.solve(symbol,geo.fft(sigma)[...,None])[...,0])
    dE=np.stack([D(xi,i,C) for i in range(3)],axis=-2)
    residual=sum(D(dE[...,i,:],i,C) for i in range(3))+sigma
    assert maxabs(residual)<2e-13
    return dE,dict(residual=maxabs(residual),xi_max=maxabs(xi),
                   sampled_symbol_min_eigenvalue=float(np.linalg.eigvalsh(symbol).min()))

@lru_cache(maxsize=1)
def setup():
    q,psi,tensor,info,c=bg.completed(15,1.)
    _,C,_=color_operator(15);k=old.tangent_matrix(q)
    x,y,z=np.moveaxis(q['grid'],-1,0)
    sigma=np.stack((.012*np.cos(x+y),.008*np.sin(y+z),.009*np.cos(z)-.005,
                    8/geo.VOL+.004*np.sin(x)),axis=-1)
    sigma_c=np.zeros(x.shape+(8,))
    sigma_c[...,0]=.003*np.sin(x);sigma_c[...,4]=.002*np.cos(y+z)
    sigma_c[...,2]=1/geo.VOL;sigma_c[...,7]=1/(np.sqrt(3)*geo.VOL)
    J=np.stack((.003+.004*np.sin(x),-.002+.003*np.cos(y),.001+.002*np.sin(z)),axis=-1)
    rho=.01+.002*np.cos(x-y)+.003*np.sin(z)
    dp,dE,dE0,ew_info=old.inverse_gauss(q,k,sigma)
    dEc,col_info=color_inverse(sigma_c)
    Ec=np.array([components(e) for e in c['E']]);Fc=np.array([[components(f) for f in row] for row in c['F']])
    color_mom=np.einsum('...ja,ija->...i',dEc,Fc)
    particular=old.momentum(q,dp,dE,dE0)+color_mom+J
    dh=np.stack([geo.derivative(q['f']['h'],i) for i in range(3)],axis=-1)
    ds=np.stack([geo.derivative(q['f']['s'],i) for i in range(3)],axis=-1)
    menus=[]
    for i in (0,1):
        p=np.zeros_like(dp);p[...,1]=dh[...,i];p[...,4]=ds[...,i]
        menus.append((p,np.zeros_like(dE),np.zeros_like(dE0)))
    fw,_=old.curvature(q);weighted=np.sin(2*z)[...,None]*fw[...,0,2,:]
    electric=np.zeros_like(dE)
    electric[...,0,:]=old.covariant(weighted,2,q)
    electric[...,2,:]=-old.covariant(weighted,0,q)
    menus.append((np.zeros_like(dp),electric,np.zeros_like(dE0)))
    matrix=np.column_stack([q['dx']**3*old.momentum(q,*m).sum(axis=(0,1,2)) for m in menus])
    before=q['dx']**3*particular.sum(axis=(0,1,2))
    coeff=np.linalg.solve(matrix,-before)
    for value,(pm,em,e0m) in zip(coeff,menus):
        dp+=value*pm;dE+=value*em;dE0+=value*e0m
    dmom=old.momentum(q,dp,dE,dE0)+color_mom+J
    return dict(q=q,psi=psi,tensor=tensor,k=k,c=c,C=C,Ec=Ec,Fc=Fc,
                sigma=sigma,sigma_c=sigma_c,J=J,rho=rho,
                dp=dp,dE=dE,dE0=dE0,dEc=dEc,dmom=dmom,
                before=before,columns=matrix,coeff=coeff,ew_info=ew_info,color_info=col_info)

def gauss_and_quantum_charge_check():
    d=setup();q=d['q'];certificate=exact_gap_certificate()
    ew=maxabs(old.gauss(q,d['k'],d['dp'],d['dE'],d['dE0'])+d['sigma'])
    color=maxabs(sum(D(d['dEc'][...,i,:],i,d['C']) for i in range(3))+d['sigma_c'])
    total=q['dx']**3*d['sigma_c'].sum(axis=(0,1,2))
    internal=q['dx']**3*sum(np.einsum('ab,...b->...a',d['C'][i],d['dEc'][...,i,:]) for i in range(3)).sum(axis=(0,1,2))
    momentum=q['dx']**3*d['dmom'].sum(axis=(0,1,2))
    assert ew<3e-12 and color<2e-13 and maxabs(total+internal)<2e-12 and maxabs(momentum)<2e-12
    expected=np.zeros(8);expected[2]=1;expected[7]=1/np.sqrt(3)
    assert maxabs(total-expected)<2e-12
    # Cauchy CAR charge difference for two orthonormal u_R spin modes, color1.
    F=np.zeros((6,2),complex);F[0,0]=1;F[1,1]=1;P=F@F.conj().T
    pair=np.array([np.trace(P@np.kron(t,np.eye(2))).real for t in bg.T])
    assert maxabs(pair-expected)<1e-14
    return dict(exact_gap_certificate=certificate,color_solver=d['color_info'],
                electroweak_Gauss_error=ew,color_Gauss_error=color,
                prescribed_total_color=total.tolist(),internal_gauge_total_color=internal.tolist(),
                total_momentum_before=d['before'].tolist(),total_momentum_after=momentum.tolist(),
                inherited_counterflow_coefficients=d['coeff'].tolist(),
                inherited_counterflow_determinant=float(np.linalg.det(d['columns'])),
                two_mode_quantum_charge_difference=pair.tolist(),
                numeric_J_rho_are_prescribed_not_claimed_quantum=True)

def updated(eps):
    d=setup();q=d['q'];z=dict(q)
    p=q['p']+eps*d['dp'];E=q['f']['E']+eps*d['dE'];E0=q['f']['E0']+eps*d['dE0']
    Ec=d['Ec']+eps*d['dEc']
    bc,bw,b0=geo.old.PAR['b']
    z['pKp']=np.sum(p*old.kinverse(q,p),axis=-1)
    z['Y']=(q['Y']+bw*(np.sum(E*E,axis=(-1,-2))-np.sum(q['f']['E']**2,axis=(-1,-2)))
            +b0*(np.sum(E0*E0,axis=-1)-np.sum(q['f']['E0']**2,axis=-1))
            +bc*(np.sum(Ec*Ec,axis=(-1,-2))-np.sum(d['Ec']**2)))
    z['C']=q['C']-2*eps*d['rho']
    z['mom']=old.momentum(q,p,E,E0)+np.einsum('...ja,ija->...i',Ec,d['Fc'])+eps*d['J']
    assert maxabs(z['mom']-(q['mom']+eps*d['dmom']))<1e-13
    assert np.min(z['C'])>0 and np.min(z['Y'])>q['Y_lower']
    return z,p,E,E0,Ec

def nonlinear_check():
    d=setup();rows=[]
    for eps in (.1,.05,-.05):
        z,p,E,E0,Ec=updated(eps);psi,tensor,info=geo.solve_hamiltonian(z,initial=d['psi'])
        ew=maxabs(old.gauss(z,d['k'],p,E,E0)+eps*d['sigma'])
        col=maxabs(sum(D(Ec[...,i,:],i,d['C']) for i in range(3))+eps*d['sigma_c'])
        mom=maxabs(sum(geo.derivative(tensor[..., :,i],i) for i in range(3))+z['mom'])
        rho=(.5*psi**-12*z['pKp']+.5*psi**-4*z['B']+z['U']+psi**-8*z['Y']+eps*d['rho'])
        h=maxabs(-8*psi**-5*geo.laplace(psi)-psi**-12*np.sum(tensor*tensor,axis=(-1,-2))+2*z['tau2']/3-2*rho)
        assert ew<3e-12 and col<3e-13 and mom<3e-12 and h<3e-8
        rows.append(dict(epsilon=eps,electroweak_Gauss_error=ew,color_Gauss_error=col,
                         total_momentum_error=mom,full_Einstein_constraint_error=h,
                         psi_change=maxabs(psi-d['psi']),C_min=float(z['C'].min()),
                         full_electric_Y_change=maxabs(z['Y']-d['q']['Y']),
                         error_if_extra_energy_omitted=info['original_Hamiltonian_residual']))
    return dict(rows=rows,fixed_original_nonzero_color_connection=True,
                finite_strength_sources_prescribed_not_self_consistent_quantum=True)

def first_order_check():
    d=setup();q=d['q'];psi=d['psi'];tensor=d['tensor']
    dtensor,_=geo.solve_momentum(dict(q,mom=d['dmom']))
    bc,bw,b0=geo.old.PAR['b']
    A=np.sum(tensor*tensor,axis=(-1,-2))+q['pKp']
    dA=2*np.sum(tensor*dtensor,axis=(-1,-2))+2*np.sum(old.kinverse(q,q['p'])*d['dp'],axis=-1)
    dY=(2*bw*np.sum(q['f']['E']*d['dE'],axis=(-1,-2))
        +2*b0*np.sum(q['f']['E0']*d['dE0'],axis=-1)+2*bc*np.sum(d['Ec']*d['dEc'],axis=(-1,-2)))
    potential=5*q['C']*psi**4-q['B']+7*A*psi**-8+6*q['Y']*psi**-4
    op=lambda v:-8*geo.laplace(v)+potential*v
    rhs=dA*psi**-7+2*dY*psi**-3+2*d['rho']*psi**5
    k2=np.sum(geo.waves(q['N'])**2,axis=-1)
    v,it=geo.cg(op,rhs,lambda x:geo.ifft(geo.fft(x)/(8*k2+potential.mean())))
    assert maxabs(op(v)-rhs)<1e-10
    rows=[]
    for eps in (.02,.01,.005):
        zp,*_=updated(eps);zm,*_=updated(-eps)
        pp,*_=geo.solve_hamiltonian(zp,initial=psi);pm,*_=geo.solve_hamiltonian(zm,initial=psi)
        error=maxabs((pp-pm)/(2*eps)-v)
        remainder=maxabs(pp-psi-eps*v)
        rows.append(dict(epsilon=eps,centered_derivative_error=error,
                         remainder_over_epsilon_squared=remainder/eps**2))
    assert rows[-1]['centered_derivative_error']<rows[0]['centered_derivative_error']/10
    return dict(response_max=maxabs(v),equation_residual=maxabs(op(v)-rhs),rows=rows,
                absolute_quantum_feedback_and_state_preparation_not_computed=True)

def run():
    a=gauss_and_quantum_charge_check();b=nonlinear_check();c=first_order_check()
    deps=('research_note_731.md','research_note_732.md','research_note_735.md','research_note_741.md',
          'research_note_753.md','joint_source_constraint_response.py','joint_reference_constraint_strata.py',
          'joint_gauss_einstein_initial_data.py')
    return dict(round=754,tests_run=3,failures=0,errors=0,
                all_color_Gauss_and_internal_balance=a,full_prescribed_source_constraints=b,
                same_source_first_order_response=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original753 nonzero color connection admits a proven all-mode covariant Gauss inverse, including nonzero color means. Original731 counterflows and conformal solve complete all prescribed smooth initial sources, with full internal cost. Numerics test specified smooth diagnostics, not a computed continuum quantum stress. A separately proven finite-rank Hadamard state comparison gives a genuine nonzero color source; inherited Ward and first-order response tools apply conditionally. No full quantum-gravity state, exact semiclassical feedback or independent Einstein emergence.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))

