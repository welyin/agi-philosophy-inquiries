"""863: a constraint-preserving original-color witness for a relational loop.
The proof of the full mixed physical quantum pairing is in research_note_863.
Numerics calibrate initial constraints and color holonomy, not a Hadamard kernel.
"""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'physical_relational_loop_bridge_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
import smeared_loop_current_probe as ward


def constants():
    a,b,c=map(Q,('.31','.27','.21'));v=b*b;w=c*c
    S=a*a*v+(a*a+v)*w/4
    Rprime=v*(v*v+S)/(4*(v+w/4)**2)*(2*w)
    uprime=-(v*v+S)/(4*(v+w/4)**2)*(2*w)
    assert Rprime>0 and v*uprime+Rprime==0
    exact=[]
    for t in map(Q,('-.1','0','.1')):
        ct=c*(1+t);ut=(S-v*ct*ct/4)/(v+ct*ct/4)
        assert ut>0 and ut*v+(ut+v)*ct*ct/4==S
        exact.append(dict(parameter=str(t),a_squared=str(ut),c_squared=str(ct*ct),total_magnetic_norm=str(S)))
    return dict(a=float(a),b=float(b),c=float(c),S=float(S),Rprime=float(Rprime),
        a_prime=float(uprime/(2*a)),rational_S=str(S),rational_Rprime=str(Rprime),exact_family=exact)


def color_at(old,t,c):
    z=c['c']*(1+t);u=(c['S']-c['b']**2*z*z/4)/(c['b']**2+z*z/4)
    A=np.array([np.sqrt(u)*old.T[0],c['b']*old.T[1],z*old.T[3]])
    E=old.color(1.)['E']
    F=np.array([[old.cross(A[i],A[j]) for j in range(3)] for i in range(3)])
    G=sum(old.cross(A[i],E[i]) for i in range(3))
    M=np.array([sum(old.inner(E[j],F[i,j]) for j in range(3)) for i in range(3)])
    S=sum(old.inner(F[i,j],F[i,j]) for i in range(3) for j in range(i+1,3))
    energy=old.geo.old.PAR['b'][0]*sum(old.inner(e,e) for e in E)+old.geo.old.PAR['K'][0]*S/2
    return dict(A=A,E=E,F=F,G=G,M=M,S=S,Y=float(energy))


def exp_with_tangent(A,dA):
    # A, dA are anti-Hermitian. Stable divided differences include degeneracies.
    e,V=np.linalg.eigh(1j*A)
    U=(V*np.exp(-1j*e))@V.conj().T
    divided=-1j*np.exp(-.5j*(e[:,None]+e[None,:]))*np.sinc((e[:,None]-e[None,:])/(2*np.pi))
    dU=V@(divided*(V.conj().T@(1j*dA)@V))@V.conj().T
    return U,dU


def loop_check(old,c):
    data=color_at(old,0.,c);A=-1j*data['A']
    dA=-1j*np.array([c['a_prime']*old.T[0],np.zeros((3,3)),c['c']*old.T[3]])
    V=A[0]-A[1];Z=A[2];dV=dA[0]-dA[1];dZ=dA[2]
    expected=-c['Rprime']/4
    rows=[]
    for L in (.8,.4,.2,.1):
        U=np.eye(3,dtype=complex);dU=np.zeros_like(U)
        for B,dB in ((-L*V,-L*dV),(-L*Z,-L*dZ),(L*V,L*dV),(L*Z,L*dZ)):
            C,dC=exp_with_tangent(B,dB);dU=dC@U+C@dU;U=C@U
        derivative=float(np.trace(dU).real);ratio=derivative/L**4
        assert derivative<0
        rows.append(dict(side=L,real_color_trace=float(np.trace(U).real),
            tangent_derivative=derivative,derivative_over_side_fourth=ratio,
            relative_error_to_small_loop_limit=abs(ratio/expected-1)))
    errors=[r['relative_error_to_small_loop_limit'] for r in rows]
    assert all(a>b for a,b in zip(errors,errors[1:])) and errors[-1]<.001
    # Independent central finite difference checks the Frechet derivative at L=.8.
    traces=[];dt=1e-4;L=.8
    for t in (-dt,dt):
        At=-1j*color_at(old,t,c)['A'];U=np.eye(3,dtype=complex)
        for B in (-L*(At[0]-At[1]),-L*At[2],L*(At[0]-At[1]),L*At[2]):
            C,_=exp_with_tangent(B,np.zeros_like(B));U=C@U
        traces.append(float(np.trace(U).real))
    fd=(traces[1]-traces[0])/(2*dt)
    error=abs(fd-rows[0]['tangent_derivative'])
    assert error<2e-10
    return dict(rows=rows,exact_leading_coefficient=str(-Q(c['rational_Rprime'])/4),
        finite_difference_crosscheck_error=error,
        original_color_connection_not_flat_replacement=True,
        quotient_representation='d_R=(3,1)_{-2}; numeric color factor only, fixed U(1) factor treated analytically')


def background_checks(old,c):
    geo=old.geo;families=[color_at(old,t,c) for t in (-.1,0.,.1)]
    ref=old.color(1.);rows=[]
    assert np.max(abs(families[1]['A']-ref['A']))<1e-15
    for N in (16,24):
        q,psi0,tensor0,_,_=old.completed(N,1.)
        x,y,z=np.moveaxis(q['grid'],-1,0);epsilon=.02;probe=epsilon*np.sin(x)
        new=dict(q);new['B']=q['B']+epsilon**2*np.cos(x)**2
        new['U']=q['U']+.5*probe**2;new['C']=q['C']-probe**2
        psi,tensor,stats=geo.solve_hamiltonian(new,initial=psi0)
        assert np.max(abs(tensor-tensor0))<1e-14
        hs=np.sqrt(geo.old.PAR['h2']);_,pars,_=geo.old.scalar.parameters();ss=np.sqrt(pars[1])
        ix=np.array([N-1,0,1]);iy=(N//4-ix)%N;iz=np.full(3,N//8)
        points=(ix,iy,iz)
        plane_h=q['phi'][points][...,1];plane_s=q['phi'][points][...,4]
        h_error=float(np.max(abs(plane_h-1.05*hs)))
        s_error=float(np.max(abs(plane_s-ss-.06*ss/epsilon*probe[points])))
        dh_error=float(np.max(abs(.05*hs*np.cos(x[points]+y[points]))))
        vel=psi[...,None]**-6*q['F'][...,None]*(q['p']-q['phi']*np.sum(q['phi']*q['p'],axis=-1)[...,None]/12)
        assert h_error<1e-14 and s_error<1e-14 and dh_error<1e-14
        assert float(np.min(vel[points][...,1]))>0
        point=(0,N//4,N//8);magz=-8*c['S']*psi[point]**-9*geo.derivative(psi,2)[point]
        J=float(vel[point][1]*.06*ss*epsilon*magz);assert J>0
        sub=[]
        for t,data in zip((-.1,0.,.1),families):
            Y=new['Y']+data['Y']-ref['Y']
            rho=.5*psi**-12*q['pKp']+.5*psi**-4*new['B']+new['U']+psi**-8*Y
            residual=-8*psi**-5*geo.laplace(psi)-psi**-12*np.sum(tensor*tensor,axis=(-1,-2))+2*q['tau2']/3-2*rho
            R=old.inner(data['F'][0,2]-data['F'][1,2],data['F'][0,2]-data['F'][1,2])
            magnetic_error=float(np.max(abs((data['S']-c['S'])*psi[points]**-8)))
            r=dict(parameter=t,Gauss_residual=float(np.max(abs(data['G']))),
                color_momentum_residual=float(np.max(abs(data['M']))),
                total_magnetic_norm_error=abs(data['S']-c['S']),energy_change=abs(data['Y']-ref['Y']),
                original_Hamiltonian_residual=float(np.max(abs(residual))),
                surface_Mh_change=magnetic_error,loop_plane_curvature_norm=R)
            assert r['Gauss_residual']<1e-14 and r['color_momentum_residual']<1e-14
            assert r['energy_change']<1e-14 and r['surface_Mh_change']<1e-14
            assert r['original_Hamiltonian_residual']<3e-8
            sub.append(r)
        rows.append(dict(N=N,probe_amplitude=epsilon,h_surface_error=h_error,
            s_p_surface_relation_error=s_error,h_spatial_gradient_surface_error=dh_error,
            central_background_reference_Jacobian=J,rows=sub))
    return rows


def run():
    w=ward.run();assert w==json.loads(ward.TARGET.read_text('utf-8'))
    c=constants()
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        loops=loop_check(old,c);background=background_checks(old,c)
    return dict(round=863,date='2026-10-06',formal_reports=863,fresh_numbered_groups=1,
        cumulative_numbered_groups=3648,all_checks_passed=True,constants=c,
        same_original_constraint_background_family=background,original_loop_response=loops,
        inherited_current_Ward_maximum=w['full_ward_maximum'],
        omitted_commutator_defect=w['omitted_commutator_maximum'],
        new_species_or_couplings_added=False,
        finite_resolution_full_mixed_linear_physical_current='analytic proof; nonzero physical pairing by constrained family',
        strictly_positive_existing_free_state_variance='analytic CCR inequality with explicit compact conserved companion source',
        full_Hadamard_variance_numerically_evaluated=False,
        spacetime_nonlinear_Cauchy_evolution_numerically_integrated=False,
        all_order_interacting_relational_loop_defined=False,
        finite_graph_quantum_state_matching_proved=False,full_goal_completed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
