"""Round 549: reproducible jets for a local analytic Einstein--matter construction.

Existence is proved in research_note_549.md, not inferred from finite jets.
The two-derivative model is NOT asserted to approximate the full spectral action.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_reference_gravity_constraints_results.json'
ETA = np.diag([-1., 1., 1., 1.])


def bare_parameters():
    r, q, s, t, c0, ng = Q(7,3), Q(1,4), Q(3,28), Q(1), Q(1,4), 3
    n = 8*ng*(3+r)  # zero-Yukawa generations still enter Tr Z
    return dict(r=r,q=q,S=s,T=t,C0=c0,ng=ng,n=n,
                M02=n*c0/(24*t),a=c0,b=c0,
                lh=(3*q*q+r*s*s)/t,p=s,ls=t/r,
                V0=n*c0*c0/(8*t),cutoff_squared=c0/2)


def fields(phi, pars):
    """Analytic target metric, derivatives, potential and its full derivative."""
    phi = np.asarray(phi, dtype=float)
    h,s = phi
    f = float(pars['M02'])-phi@phi/6
    assert f > 0
    df, ddf = -phi/3, -np.eye(2)/3
    g = np.eye(2)/f + 1.5*np.outer(df,df)/f**2
    dg = np.empty((2,2,2))  # derivative, row, column
    for k in range(2):
        dg[k] = (-np.eye(2)*df[k]/f**2
                 +1.5*(np.outer(ddf[:,k],df)+np.outer(df,ddf[:,k]))/f**2
                 -3*np.outer(df,df)*df[k]/f**3)
    gi = np.linalg.inv(g)
    dgi = np.array([-gi@dgk@gi for dgk in dg])
    a,b,lh,p,ls,v0 = (float(pars[k]) for k in ('a','b','lh','p','ls','V0'))
    v = v0-(a*h*h+b*s*s)/2+(lh*h**4+2*p*h*h*s*s+ls*s**4)/4
    dv = np.array([-a*h+lh*h**3+p*h*s*s, -b*s+ls*s**3+p*s*h*h])
    u,du = v/f**2, dv/f**2-2*v*df/f**3
    conn = np.zeros((2,2,2))
    for i,j,k,l in itertools.product(range(2),repeat=4):
        conn[i,j,k] += .5*gi[i,l]*(dg[j,l,k]+dg[k,l,j]-dg[l,j,k])
    return dict(F=f,dF=df,G=g,Gi=gi,dG=dg,dGi=dgi,V=v,dV=dv,U=u,dU=du,connection=conn)


def origin_jets(phi, eps, pars):
    data=fields(phi,pars)
    g,gi=data['G'],data['Gi']
    p0=eps/gi[0,0]
    grad=np.zeros((2,4)); grad[:,0]=gi[:,0]*p0; grad[1,1]=eps
    hess=np.zeros((2,4,4))
    hess[:,0,1]=data['dGi'][1,:,0]*eps*p0
    hess[:,0,2]=gi[:,0]*p0
    hess[:,1,0]=hess[:,0,1]; hess[:,2,0]=hess[:,0,2]
    hess[1,1,3]=hess[1,3,1]=eps
    contraction=grad@ETA@grad.T
    # U need not be stationary even when V is stationary.
    acceleration=np.einsum('ijk,jk->i',data['connection'],contraction)-gi@data['dU']
    hess[:,0,0]=acceleration
    jj=np.vstack([grad,2*grad[0]@ETA@hess[0],2*grad[1]@ETA@hess[1]])
    invariants=np.diag(contraction)
    menu_map=np.eye(4)
    menu_map[2:,:2]=np.outer(invariants,data['dF'])
    menu_map[2:,2:]=data['F']*np.eye(2)
    jj_j=menu_map@jj
    rho2=float(grad[:,0]@g@grad[:,0]+grad[:,1]@g@grad[:,1]+2*data['U'])
    return data,grad,hess,jj,jj_j,-rho2/8


def exact_det(matrix):
    ans=Q(0)
    for perm in itertools.permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        product=Q(sign)
        for i,j in enumerate(perm): product*=matrix[i][j]
        ans+=product
    return ans


def scalar_curvature_from_metric_hessian(psi_zz):
    """At psi=1, dpsi=0, gamma_ij=psi^4 delta_ij; use Ricci definition."""
    ddmetric=np.zeros((3,3,3,3))  # derivative, derivative, metric, metric
    ddmetric[2,2]=4*psi_zz*np.eye(3)
    dgamma=np.zeros((3,3,3,3))  # derivative, upper, lower, lower
    for a,k,i,j in itertools.product(range(3),repeat=4):
        dgamma[a,k,i,j]=.5*(ddmetric[a,i,k,j]+ddmetric[a,j,k,i]-ddmetric[a,k,i,j])
    ricci=np.array([[sum(dgamma[k,k,i,j]-dgamma[j,k,i,k] for k in range(3))
                     for j in range(3)] for i in range(3)])
    return float(np.trace(ricci))


def run():
    checks=[]; pars=bare_parameters()
    assert pars['T']==3*pars['q']+pars['r']*pars['S']
    assert pars['n']==128 and pars['M02']==Q(4,3)
    h2=pars['C0']/pars['q']
    s2=pars['C0']*pars['r']/pars['T']*(1-pars['S']/pars['q'])
    fstar=pars['M02']-(h2+s2)/6
    vstar=pars['V0']-pars['C0']*(h2+s2)/2+(pars['lh']*h2*h2+2*pars['p']*h2*s2+pars['ls']*s2*s2)/4
    assert (h2,s2,fstar,vstar,vstar/fstar**2)==(Q(1),Q(1,3),Q(10,9),Q(11,12),Q(297,400))
    checks.append('exact_same_scale_coefficients_include_all_generations_and_vacuum_constant')

    phi=np.array([1.,np.sqrt(1/3)])
    fd_error=0.; inverse_error=0.; metric_floor=10.
    for state in (phi,np.array([.7,.4]),np.array([1.2,-.2])):
        data=fields(state,pars)
        metric_floor=min(metric_floor,float(np.linalg.eigvalsh(data['G'])[0]))
        for k in range(2):
            step=np.eye(2)[k]*1e-5
            plus,minus=fields(state+step,pars),fields(state-step,pars)
            fd_error=max(fd_error,float(np.max(abs((plus['G']-minus['G'])/2e-5-data['dG'][k]))),
                         abs((plus['U']-minus['U'])/2e-5-data['dU'][k]))
            inverse_error=max(inverse_error,float(np.max(abs((plus['Gi']-minus['Gi'])/2e-5-data['dGi'][k]))))
    assert metric_floor>0 and max(fd_error,inverse_error)<1e-8
    checks.append('positive_mixed_target_metric_and_independent_derivative_differences')

    eps=.25; anchor=fields(phi,pars); momentum_error=0.; constraint_error=0.
    for x,y,z,psi in ((0.,0.,0.,1.),(.1,-.2,.1,.91),(-.2,.3,-.1,1.1)):
        state=phi+[0,eps*(x+x*z)]
        data=fields(state,pars); g,gi=data['G'],data['Gi']
        p=eps*(1+y)/anchor['Gi'][0,0]
        velocity=psi**-6*gi[:,0]*p
        spacegrad=np.array([[0.,0.,0.],[eps*(1+z),0,eps*x]])
        j=-velocity@g@spacegrad
        momentum_error=max(momentum_error,float(np.max(abs(j))))
        lapseq=-(psi**-7*gi[0,0]*p*p+psi*g[1,1]*(spacegrad[1]@spacegrad[1])+2*data['U']*psi**5)/8
        ricci=-8*psi**-5*lapseq
        source=velocity@g@velocity+psi**-4*np.einsum('ij,ik,jk',g,spacegrad,spacegrad)+2*data['U']
        constraint_error=max(constraint_error,abs(ricci-source))
    naive_jx=-anchor['G'][0,1]*eps*eps
    assert momentum_error<1e-15 and constraint_error<1e-13 and abs(naive_jx)>1e-4
    checks.append('covariant_momentum_solves_constraint_naive_zero_s_velocity_does_not')

    exact=[]
    for e,c in ((Q(1,4),Q(-2,7)),(Q(1,8),Q(5,2)),(Q(1,3),Q(0))):
        jj=[[e,0,0,0],[c*e,e,0,0],[Q(3,7),Q(-2,9),-2*e*e,0],
            [Q(-8,3),Q(4,5),-2*c*c*e*e,2*e*e]]
        det=exact_det(jj); assert det==-4*e**6
        exact.append(dict(epsilon=str(e),mixed_ratio=str(c),determinant=str(det)))
    checks.append('exact_full_rank_certificate_with_arbitrary_front_columns')

    examples=[]; wave_error=0.; det_error=0.; ricci_error=0.
    for state in (phi,np.array([.7,.4]),np.array([1.2,-.2])):
        for e in (.5,.25,.125):
            data,grad,hess,jj,jj_j,psi_zz=origin_jets(state,e,pars)
            residual=np.einsum('mn,imn->i',ETA,hess)+np.einsum('ijk,jk->i',data['connection'],grad@ETA@grad.T)-data['Gi']@data['dU']
            wave_error=max(wave_error,float(np.max(abs(residual))))
            det_error=max(det_error,abs(np.linalg.det(jj)/(-4*e**6)-1))
            ricci_error=max(ricci_error,abs(scalar_curvature_from_metric_hessian(psi_zz)+8*psi_zz))
            assert data['F']>0 and grad[0]@ETA@grad[0]<0
            assert np.max(abs(jj[:2,2:]))<1e-15
            examples.append(dict(fields=state.tolist(),epsilon=e,F=data['F'],V=data['V'],U=data['U'],
                s_velocity=float(grad[1,0]),field_accelerations=hess[:,0,0].tolist(),
                determinant_E=float(np.linalg.det(jj)),determinant_J=float(np.linalg.det(jj_j)),
                psi_zz=psi_zz,initial_spatial_scalar_curvature=-8*psi_zz))
    assert wave_error<1e-14 and det_error<1e-13 and ricci_error<1e-14
    assert np.max(abs(anchor['dV']))<1e-15 and np.linalg.norm(anchor['dU'])>.1
    checks.append('complete_sigma_equation_jets_and_Ricci_definition_match_initial_constraints')

    data,grad,hess,jj,jj_j,_=origin_jets(phi,.25,pars)
    def menu(x,jordan=False):
        values=phi+grad@x+.5*np.einsum('imn,m,n->i',hess,x,x)
        deriv=grad+np.einsum('imn,n->im',hess,x)
        invariants=np.diag(deriv@ETA@deriv.T).copy()
        if jordan: invariants*=fields(values,pars)['F']
        return np.r_[values,invariants]
    # The frozen tangent metric differs from the actual metric by O(|x|^2).
    # These finite Taylor representatives validate first derivatives only.
    menu_error=0.; jordan_error=0.
    for k in range(4):
        step=np.eye(4)[k]*2e-5
        menu_error=max(menu_error,float(np.max(abs((menu(step)-menu(-step))/4e-5-jj[:,k]))))
        jordan_error=max(jordan_error,float(np.max(abs((menu(step,True)-menu(-step,True))/4e-5-jj_j[:,k]))))
    assert max(menu_error,jordan_error)<1e-8
    assert abs(np.linalg.det(jj_j)/(data['F']**2*np.linalg.det(jj))-1)<1e-13
    checks.append('independent_reference_menu_differences_and_Jordan_frame_factor')

    shifted=dict(pars); shifted['V0']+=Q(2,5)
    d2,g2,h2,j2,jj2,zz2=origin_jets(phi,.25,shifted)
    _,_,_,_,_,zz1=origin_jets(phi,.25,pars)
    curvature_shift=-8*(zz2-zz1)
    expected_shift=2*float(Q(2,5))/data['F']**2
    assert abs(curvature_shift-expected_shift)<1e-14
    assert abs(np.linalg.det(j2)-np.linalg.det(jj))<1e-15
    assert np.linalg.norm(h2[:,0,0]-hess[:,0,0])>.1
    checks.append('vacuum_constant_changes_curvature_and_acceleration_not_rank')

    _,_,_,zero_j,_,_=origin_jets(phi,0.,pars)
    assert np.linalg.matrix_rank(zero_j)<4
    failed=False
    try: fields([3.,3.],pars)
    except AssertionError: failed=True
    assert failed
    checks.append('zero_reference_amplitude_and_nonpositive_F_are_outside_positive_conclusion')

    dependencies=('research_note_533.md','research_note_536.md','research_note_543.md',
        'research_note_548.md','research_round_548_checks.json','joint_condition_compression_table.md')
    return dict(round=549,tests_run=len(checks),failures=0,errors=0,checks=checks,
        bare_parameters={k:str(v) for k,v in pars.items()},exact_stationary_values=dict(
            h_squared='1',s_squared='1/3',F='10/9',V='11/12',U='297/400'),
        exact_determinant_certificates=exact,mixed_target_metric=anchor['G'].tolist(),
        metric_derivative_error=fd_error,inverse_metric_derivative_error=inverse_error,
        momentum_residual=momentum_error,naive_momentum_jx=float(naive_jx),
        Hamiltonian_algebra_error=constraint_error,wave_equation_residual=wave_error,
        determinant_relative_error=det_error,Ricci_definition_error=ricci_error,
        menu_difference_error=menu_error,Jordan_menu_difference_error=jordan_error,
        retained_constant_curvature_shift=curvature_shift,coupled_origin_examples=examples,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
        scope=dict(local_analytic_constraint_solution_proved_in_note=True,
          local_coupled_development_uses_harmonic_reduction=True,
          numerics_validate_jets_not_full_PDE=True,same_nonminimal_Higgs_and_explicit_singlet=True,
          declared_two_derivative_classical_branch=True,full_spectral_action_solution=False,
          controlled_spectral_derivative_expansion=False,constant_potential_retained=True,
          prescribed_four_dimensions_and_action=True,no_extra_four_reference_fields=True,
          global_extension_or_finite_total_energy_proved=False,actual_quantum_instrument=False,
          unique_dimension_selected=False,completed_unification=False))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else: assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','momentum_residual',
        'wave_equation_residual','menu_difference_error','Jordan_menu_difference_error')}))
