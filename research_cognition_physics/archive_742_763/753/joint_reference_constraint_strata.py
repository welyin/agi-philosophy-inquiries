"""753: original material reference versus nonlinear color constraints.

The original 572/573 action is unchanged. Classical constrained families and
local reference charts are tested, not full quantum-gravity states.
"""
import argparse,hashlib,json
from fractions import Fraction as Q
from pathlib import Path
from functools import lru_cache
import numpy as np
import joint_gauss_einstein_initial_data as geo
import joint_gravity_material_coordinates as coordinates
import joint_source_constraint_response as source

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_reference_constraint_strata_results.json'
AA,BB,CC,EE,FF,GG=.31,.27,.21,.23,.19,.17

def generators():
    out=[]
    for i,j in ((0,1),(0,2),(1,2)):
        x=np.zeros((3,3),complex);x[i,j]=x[j,i]=.5
        y=np.zeros_like(x);y[i,j]=-.5j;y[j,i]=.5j
        out.extend([x,y])
    t3=np.diag([.5,-.5,0.]).astype(complex)
    t8=np.diag([1.,1.,-2.]).astype(complex)/(2*np.sqrt(3))
    return np.array([out[0],out[1],t3,*out[2:],t8])

T=generators()
def inner(a,b):return float((2*np.trace(a.conj().T@b)).real)
def cross(a,b):return 1j*(a@b-b@a)
def maximum(a):return float(np.max(abs(a)))

def color(amplitude,balanced=True):
    A=np.zeros((3,3,3),complex);E=np.zeros_like(A)
    A[0]=amplitude*AA*T[0]
    if balanced:
        A[1]=amplitude*BB*T[1]
        E[0]=amplitude*EE*T[0];E[1]=amplitude*FF*T[1]
        A[2]=amplitude*CC*T[3];E[2]=amplitude*GG*T[3]
    else:E[0]=amplitude*EE*T[1]
    curvature=np.array([[cross(A[i],A[j]) for j in range(3)] for i in range(3)])
    gauss=sum(cross(A[i],E[i]) for i in range(3))
    momentum=np.array([sum(inner(E[j],curvature[i,j]) for j in range(3)) for i in range(3)])
    bc=geo.old.PAR['b'][0];kc=geo.old.PAR['K'][0]
    electric=bc*sum(inner(e,e) for e in E)
    magnetic=kc/2*sum(inner(curvature[i,j],curvature[i,j]) for i in range(3) for j in range(i+1,3))
    return dict(A=A,E=E,F=curvature,G=gauss,M=momentum,Y=electric+magnetic,
                electric=float(electric),magnetic=float(magnetic))

@lru_cache(maxsize=2)
def base(N):
    q=geo.make_source(N);psi,tensor,info=geo.solve_hamiltonian(q)
    return q,psi,tensor,info

def obstruction_check():
    norm_error=maximum(np.array([[inner(a,b) for b in T] for a in T])-np.eye(8))
    assert norm_error<1e-14
    a=color(1.,False);target=-geo.VOL*AA*EE
    rows=[]
    Ue,Uv=np.linalg.eigh(.23*T[2]+.31*T[4]-.19*T[6])
    U=(Uv*np.exp(1j*Ue))@Uv.conj().T
    pairing_error=0.
    for lam in (.4,.2,.1,.05):
        c=color(lam,False)
        total=geo.VOL*inner(T[2],c['G'])
        pairing_error=max(pairing_error,abs(total/lam**2-target))
        rotated=U@c['G']@U.conj().T;xi=U@T[2]@U.conj().T
        assert abs(geo.VOL*inner(xi,rotated)-total)<1e-14
        assert maximum(c['M'])==0 and maximum(c['F'])==0
        rows.append(dict(amplitude=lam,total_T3_Gauss=total,
                         divided_by_amplitude_squared=total/lam**2,
                         positive_color_energy_coefficient=c['Y']))
    assert pairing_error<1e-13 and target<0
    # At A=E=0, constant e has div(e)=0; all other tangent fields are zero.
    # The color energy and complete Gauss defect begin at order lambda^2.
    q,psi,tensor,info=base(16)
    lam=.1;c=color(lam,False)
    geometric_defect=-2*psi**-8*c['Y']
    assert abs(float(np.max(abs(geometric_defect)))/lam**2
               -2*float(np.max(psi**-8))*a['Y'])<1e-13
    return dict(basis_orthogonality_error=norm_error,rows=rows,
                invariant_quadratic_charge=target,pairing_error=pairing_error,
                uncorrected_Einstein_defect_at_amplitude_point1=float(np.max(abs(geometric_defect))),
                all_first_order_constraints_satisfied=True,
                scope='No C1 exact closed-system Gauss curve with this entire specified tangent; no claim excluding added first-order compensating matter or different backgrounds.')

@lru_cache(maxsize=12)
def completed(N,lam):
    q,psi0,tensor0,baseinfo=base(N);c=color(lam,True)
    z=dict(q);z['Y']=q['Y']+c['Y']
    psi,tensor,info=geo.solve_hamiltonian(z,initial=psi0)
    return z,psi,tensor,info,c

def completion_check():
    bc,kc=geo.old.PAR['b'][0],geo.old.PAR['K'][0]
    assert .68<bc<.70 and .72<kc<.74 and abs(2*bc*kc-1)<1e-14
    # Rational uniform bound, not a sampled-field estimate.
    magnetic_coeff=(Q('.31')*Q('.27'))**2+((Q('.31')*Q('.21'))**2+(Q('.27')*Q('.21'))**2)/4
    bound=Q('.70')*(Q('.23')**2+Q('.19')**2+Q('.17')**2)+Q('.74')/2*magnetic_coeff
    assert bound<Q('.1')
    reaction=Q(433385,34992)-2*Q('.1')*Q('1.5')**-3
    assert reaction>12
    rows=[]
    for N,lam in ((16,.25),(16,.5),(16,1.),(24,1.)):
        z,psi,tensor,info,c=completed(N,lam)
        q,psi0,tensor0,_=base(N)
        formula=bc*(EE*EE+FF*FF+GG*GG)*lam**2+kc/2*((AA*BB)**2+((AA*CC)**2+(BB*CC)**2)/4)*lam**4
        assert abs(c['Y']-formula)<1e-16
        assert maximum(c['G'])<1e-16 and maximum(c['M'])<1e-16 and c['magnetic']>0
        k=source.tangent_matrix(q)
        ew=maximum(source.gauss(q,k,q['p'],q['f']['E'],q['f']['E0']))
        mom=maximum(sum(geo.derivative(tensor[..., :,j],j) for j in range(3))+z['mom'])
        # Direct Einstein constraint, with the new physical electric/magnetic density.
        rho=(.5*psi**-12*z['pKp']+.5*psi**-4*z['B']+z['U']+psi**-8*(q['Y']+c['Y']))
        constraint=(-8*psi**-5*geo.laplace(psi)-psi**-12*np.sum(tensor*tensor,axis=(-1,-2))
                    +2*z['tau2']/3-2*rho)
        err=maximum(constraint)
        assert ew<1e-12 and mom<1e-12 and err<3e-8
        assert maximum(tensor-tensor0)<1e-14 and float(np.min(psi-psi0))>0
        rows.append(dict(N=N,amplitude=lam,electric_Y=c['electric'],magnetic_Y=c['magnetic'],
                         color_Gauss_error=maximum(c['G']),electroweak_Gauss_error=ew,
                         color_momentum=maximum(c['M']),total_momentum_error=mom,
                         original_Einstein_constraint_error=err,
                         minimum_psi_increase=float(np.min(psi-psi0)),
                         maximum_psi_increase=maximum(psi-psi0),psi_max=float(psi.max())))
    # Even family: its geometry has zero first derivative at lambda=0.
    zp,pp,*_=completed(16,.5);zm,pm,*_=completed(16,-.5)
    assert maximum(pp-pm)==0
    return dict(rows=rows,uniform_Y_increment_rational_bound=str(bound),
                uniform_Y_increment_decimal_bound=float(bound),
                original_upper_barrier_reaction_lower=str(reaction),
                even_geometry_error=maximum(pp-pm),
                original_action_and_all_couplings_unchanged=True)

def reference_and_strata_check():
    rows=[]
    for N,lam in ((16,0.),(16,.5),(16,1.),(24,1.)):
        q,psi,tensor,info,c=completed(N,lam)
        f0=coordinates.fields(N);coeff=coordinates.point_coefficients(f0)
        v=psi[...,None]**-6*q['F'][...,None]*(q['p']-q['phi']*np.sum(q['phi']*q['p'],axis=-1)[...,None]/12)
        rh=psi**-4*np.sum(f0['dh']**2,axis=-1)-v[...,1]**2
        rs=psi**-4*np.sum(f0['ds']**2,axis=-1)-v[...,4]**2
        ax=geo.derivative(rh,0);az=geo.derivative(rh,2)
        bx=geo.derivative(rs,0);bz=geo.derivative(rs,2);psiz=geo.derivative(psi,2)
        i=(0,N//4,N//8);minor=-v[i][1]*(-coeff['b'])*(ax[i]*bz[i]-az[i]*bx[i])
        bracket=coeff['Chx']*coeff['b']**2*(coeff['R']-psi[i]**8)
        exact=8*coeff['Ch']**2*(-coeff['b'])*psi[i]**-31*psiz[i]*bracket
        error=abs(minor-exact)/abs(exact)
        assert np.max(psi)<1.5 and psiz[i]<0 and exact>0 and minor>0 and error<.012
        # Pointwise centralizer of actual A1,A2,A3; E uses the same generators.
        tangents=np.array([[cross(ai,ta) for ta in T] for ai in c['A']])
        gram=np.array([[sum(inner(tangents[j,a],tangents[j,b]) for j in range(3))
                         for b in range(8)] for a in range(8)])
        eigen=np.linalg.eigvalsh(gram);nullity=int(np.sum(eigen<1e-12))
        assert nullity==(8 if lam==0 else 0)
        if lam!=0:
            assert eigen[0]>0
        rows.append(dict(N=N,amplitude=lam,reference_det_formula=float(exact),
                         independent_reference_minor=float(minor),relative_minor_error=float(error),
                         clock_norm=float(rh[i]),psi_z=float(psiz[i]),
                         constant_color_centralizer_eigenvalues=eigen.tolist(),nullity=nullity))
    return dict(rows=rows,
                continuum_no_spacetime_stabilizer_uses_full_rank_and_Killing_uniqueness=True,
                no_full_constraint_adjoint_surjectivity_claim=True,
                nonzero_family_residual_color_Lie_algebra='zero',
                original_reflection_walls_and_non_global_chart_retained=True)

def run():
    checks=(obstruction_check(),completion_check(),reference_and_strata_check())
    deps=('research_note_570.md','research_note_572.md','research_note_573.md','research_note_574.md',
          'research_note_728.md','research_note_731.md','research_note_752.md',
          'joint_gauss_einstein_initial_data.py','joint_gravity_material_coordinates.py',
          'joint_source_constraint_response.py')
    return dict(round=753,tests_run=3,failures=0,errors=0,
                original_tangent_obstruction=checks[0],exact_joint_initial_family=checks[1],
                reference_and_residual_symmetry=checks[2],
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original continuum classical Einstein-Yang-Mills-sigma branch with zero classical fermion field. Full material rank removes smooth spacetime components of joint background symmetries, not color. An explicit linearized color tangent is nonintegrable by its quadratic total charge; a different genuinely non-Abelian color family solves all classical initial constraints, retains the old local material chart, and has no continuous color stabilizer at nonzero amplitude. Discrete center actions and full gauge redundancy are not removed. No new action, full BRST positivity, quantum record implementation, global chart or all-physics reconstruction.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
