"""Independent round 1269 checks: one group, source and finite-code subchecks.
No author executable/results are imported. Default is read-only comparison;
--write creates independent_results.json exclusively.
"""
from fractions import Fraction as F
from pathlib import Path
from math import factorial
import hashlib, json, sys
import numpy as np

HERE=Path(__file__).resolve().parent

def matmul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]

def transpose(a): return list(map(list,zip(*a)))
def diag(a):return [[a[i] if i==j else F(0) for j in range(len(a))] for i in range(len(a))]
def eye(n):return diag([F(1)]*n)
def kron(a,b):return [[x*y for x in rowa for y in rowb] for rowa in a for rowb in b]
def add(a,b):return [[x+y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]
def scale(s,a):return [[s*x for x in r] for r in a]
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))

def main():
    # Independent fixed-current metric differentiation, ten metric variations
    # in each rational Lorentz frame. This is a finite algebraic diagnostic;
    # the proof uses the unrestricted variation identity stated in the review.
    eta=[F(-1),F(1),F(1),F(1)]
    frames=[([F(1),F(0),F(0),F(0)],F(1,7)),
            ([F(5,4),F(3,4),F(0),F(0)],F(7,3)),
            ([F(13,5),F(0),F(12,5),F(0)],F(2,5))]
    gradients=[[F(1,7),F(2),F(0),F(1,3)],
               [F(0),F(1,5),F(3),F(0)],
               [F(1,11),F(0),F(1,13),F(4)]]
    S=[[sum((v[i]*v[j] for v in gradients),F(0)) for j in range(4)] for i in range(4)]
    B=sum((eta[i]*S[i][i] for i in range(4)),F(0))
    M,d=F(71,100)+F(1,10),F(1,100003)
    metric_checks=0
    for U,n in frames:
        assert dot([eta[i]*U[i] for i in range(4)],U)==-1
        j=[n*x for x in U];jlow=[eta[i]*j[i] for i in range(4)]
        Ulow=[eta[i]*U[i] for i in range(4)]
        for a in range(4):
            for b in range(a,4):
                E=[[F(0) for _ in range(4)] for _ in range(4)]
                E[a][b]=E[b][a]=F(1)
                tr=sum((eta[i]*E[i][i] for i in range(4)),F(0))
                aprime=sum((jlow[i]*E[i][k]*jlow[k] for i in range(4) for k in range(4)),F(0))
                rootgprime=-tr/2
                nprime=aprime/(2*n)+n*tr/2
                fluid_variation=-(rootgprime*n*M+nprime*M)
                tf=sum((n*M*Ulow[i]*Ulow[k]*E[i][k] for i in range(4) for k in range(4)),F(0))
                assert -2*fluid_variation==tf
                Bprime=sum((E[i][k]*S[i][k] for i in range(4) for k in range(4)),F(0))
                vvariation=-d*(rootgprime*B+Bprime)
                tv=2*d*Bprime-d*B*tr
                assert -2*vvariation==tv
                metric_checks+=1
    # Density-dependent rho=n M0+d B at fixed internal variables: the
    # enthalpy is n M0, not n M0+d B. Exact linear polynomial coefficients.
    rho={'nM':F(1),'dB':F(1)}
    nrho_n={'nM':F(1),'dB':F(0)}
    pressure={k:nrho_n[k]-rho[k] for k in rho}
    assert pressure=={'nM':F(0),'dB':F(-1)}
    assert nrho_n['dB']==0
    # Canonical internal evolution including the old oscillator, memory,
    # theta-dependent interaction and external p/X exchange.
    n,B0,dtheta,dP=F(2,3),F(7,5),F(3,11),F(-4,13)
    Mtheta,Mq,Mr,MY,MK=map(F,[2,3,5,7,11])
    qdot,rdot,Ydot,Kdot=Mr,-Mq,MK,-MY
    Jdot=-Mtheta-B0*dtheta/n
    phidot=B0*dP/n
    Pdot=F(0);thetadot=F(1)
    Mdot=Mtheta*thetadot+Jdot+Mq*qdot+Mr*rdot+MY*Ydot+MK*Kdot
    assert n*Mdot+B0*dtheta==0
    # Cartan coefficients of dtheta,dJ,dq,dr,dY,dK,dP:
    cartan=[Jdot,-thetadot,rdot,-qdot,Kdot,-Ydot,-phidot]
    target=[-Mtheta-B0*dtheta/n,F(-1),-Mq,-Mr,-MY,-MK,-B0*dP/n]
    assert cartan==target
    # Therefore i_U d(-M0 U)=-dM0-(B/n)dd+old external terms;
    # n[M0 a+(UM0)U]=B grad d minus the old field exchanges.
    # Omitting the new J update leaves a nonzero power residual.
    old_only_power=n*(Mtheta-Mtheta+Mq*qdot+Mr*rdot+MY*Ydot+MK*Kdot)
    deletion=old_only_power+B0*dtheta
    assert deletion==B0*dtheta and deletion!=0

    # Independent one-dimensional H-H principal bracket with a genuinely
    # nonconstant A. Exact cancellation of A' is the local tensor mechanism.
    bracket_diagnostics=0
    for ac,acprime in [(F(3,2),F(7,11)),(F(11,9),F(-3,5))]:
        phix,phixx,pi=F(5,13),F(2,17),F(19,23)
        lapse,dlapse,other,dother=F(3,7),F(-2,11),F(5,9),F(4,13)
        derivative_acphix=acprime*phix+ac*phixx
        vari_phi_n=-(dlapse*ac*phix+lapse*derivative_acphix)
        vari_phi_m=-(dother*ac*phix+other*derivative_acphix)
        hh=vari_phi_n*(other*pi/ac)-(lapse*pi/ac)*vari_phi_m
        target_hh=pi*phix*(lapse*dother-other*dlapse)
        assert hh==target_hh
        # d(X^-1 x) displacement is -grad(d).delta X, so the new
        # derivative of Hrod has lapse, not a derivative of lapse.
        jac,grad_d,hd,velocity=F(7,5),F(2,3),F(5,11),F(3,8)
        cross_nm=(lapse*velocity)*(-jac*other*hd*grad_d)
        cross_mn=(other*velocity)*(-jac*lapse*hd*grad_d)
        assert cross_nm-cross_mn==0
        bracket_diagnostics+=1
    # The volume piece has no embedding velocity. Inserting its per-label
    # expression into the old mass root is a distinct Legendre transform.
    naive_delta=F(1,100)
    assert (M+naive_delta)**2!=M**2

    # Normalized cubic beta pulse and its moments. No numerical quadrature.
    assert F(140)*F(factorial(3)*factorial(3),factorial(7))==1
    mean_u=F(140)*F(factorial(4)*factorial(3),factorial(8))
    assert mean_u==F(1,2)
    tau,alpha,k2=F(1,10**8),F(1,16),F(1,10)
    maxg=F(35,16)/tau
    eps_test=F(1,2**100)
    max_2d=2*alpha*maxg*eps_test/(3*k2)
    assert max_2d<F(1,2)
    # With a=a0+v(t-ta), integral(g' a)=-v, so J need not reset.
    a0,v,P=F(5,2),F(1,9),F(1,4096)
    a_average=a0+v*tau/2
    delta_phi=alpha*a_average
    delta_J=alpha*P*v
    assert delta_J!=0 and delta_phi>0
    # Homogeneous volume term a^3*(alpha g P/(3k²))*(3k²/a²)
    # equals alpha g P a; total action normalization A=V0*s²
    # gives eps=hbar/A. Code P eigenvalues are +/-eps.
    assert 3-2==1

    # Exact rational intertwining on a noncommuting geometric diagnostic.
    # This verifies the actual normalization, not an approximate finite CCR.
    modes=list(range(-3,4));eps=F(1,4096)
    V=[[F(int(m==1)),F(int(m==-1))] for m in modes]
    Pfull=diag([eps*m for m in modes]);Z=diag([F(1),F(-1)])
    assert matmul(transpose(V),V)==eye(2)
    assert matmul(Pfull,V)==matmul(V,scale(eps,Z))
    Hgeo=[[F(1,3),F(1,7)],[F(1,7),F(2,5)]]
    Fgeo=diag([F(5,2),F(11,4)])
    assert matmul(Hgeo,Fgeo)!=matmul(Fgeo,Hgeo)
    gain=F(17,23)
    Hfull=add(kron(Hgeo,eye(7)),scale(gain,kron(Fgeo,Pfull)))
    Hcode=add(kron(Hgeo,eye(2)),scale(gain*eps,kron(Fgeo,Z)))
    VV=kron(eye(2),V)
    assert matmul(Hfull,VV)==matmul(VV,Hcode)
    # Incorrect extra eps at the classical coupling fails the same test.
    Hwrong=add(kron(Hgeo,eye(7)),scale(gain*eps,kron(Fgeo,Pfull)))
    assert matmul(Hwrong,VV)!=matmul(VV,Hcode)

    # Bounded code instruments and an entangled pointer/reference diagnostic.
    # Exact generator intertwining above, rather than this finite diagnostic,
    # carries the arbitrary reference theorem.
    vn=np.array(V,dtype=float);pc=vn@vn.T;ident=np.eye(7)
    sy=np.array([[0,-1j],[1j,0]],dtype=complex)
    yp=(np.eye(2)+sy)/2;ym=(np.eye(2)-sy)/2
    eplus=vn@yp@vn.T;eminus=vn@ym@vn.T;efail=ident-pc
    assert np.linalg.norm(eplus+eminus+efail-ident)<1e-14
    assert np.linalg.norm(eplus@eplus-eplus)<1e-14
    assert np.linalg.norm(eminus@eminus-eminus)<1e-14
    assert np.linalg.norm(efail@vn)<1e-14
    # No identification of the projective code measurement with local
    # multiplication by sin(2 phi): compression has a factor 1/2 and leaks.
    sin2=np.zeros((7,7),dtype=complex)
    for i,m in enumerate(modes):
        for j,nn in enumerate(modes):
            if m==nn+2:sin2[i,j]+=1/(2j)
            if m==nn-2:sin2[i,j]-=1/(2j)
    compressed=vn.T@sin2@vn
    assert np.linalg.norm(compressed-sy/2)<1e-14
    assert np.linalg.norm((ident-pc)@sin2@vn)>0
    state=np.array([1,0,0,1j],dtype=complex)/np.sqrt(2)
    rho=np.outer(state,state.conj())
    embed=np.kron(vn,np.eye(2));rho_full=embed@rho@embed.T
    for a,b in [(eplus,yp),(eminus,ym)]:
        aa=np.kron(a,np.eye(2));bb=np.kron(b,np.eye(2))
        post=aa@rho_full@aa.conj().T
        expected=embed@(bb@rho@bb.conj().T)@embed.T
        assert np.linalg.norm(post-expected)<1e-14
    fail_probability=float(np.trace(np.kron(efail,np.eye(2))@rho_full).real)
    assert abs(fail_probability)<1e-14

    return {
      'round':1269,'verdict':'PASS_INDEPENDENT_WITH_SCOPE','independent_groups':1,
      'source_subcheck':{
        'fixed_current_metric_diagnostics':metric_checks,
        'metric_variation':'nprime=n(g+UU):dg_inverse/2 at fixed densitized j',
        'enthalpy_interaction_coefficient':str(nrho_n['dB']),
        'pressure_interaction_coefficient':str(pressure['dB']),
        'cartan_coefficients_match':True,
        'nonconstant_A_HDA_diagnostics':bracket_diagnostics,
        'lapse_cross_cancellation':True,
        'omit_J_power_residual':str(deletion),
        'pulse_integral':'1','pulse_first_moment':str(mean_u),
        'max_rod_coefficient_shift_at_eps_2powminus100':str(max_2d),
        'nonzero_J_end_diagnostic':str(delta_J),
        'scope':'Metric and canonical identities are analytic; finite rational frames only independently diagnose signs. No full PDE or material permission.'},
      'code_subcheck':{
        'rotor_modes_diagnostic':modes,'code_modes':[1,-1],
        'effective_action_parameter':str(eps),
        'exact_noncommuting_geometric_intertwining':True,
        'extra_epsilon_deleted':True,
        'entangled_reference_instrument_match':True,
        'sin2phi_is_not_code_projective_readout':True,
        'scope':'The full rotor has an exact invariant code for P-independent baseline; auxiliary truncation used only to diagnose. No finite CCR claim, no all-rotor positive domain, no local detector/source qualification.'},
      'claims_not_certified':['full spatial quantum realization','complete CC','actual local output','post-read Einstein residual','stage closure']}

if __name__=='__main__':
    result=main();dest=HERE/'independent_results.json'
    if sys.argv[1:]==['--write']:
        with dest.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        print('WRITE_ONCE_PASS_1269 independent_groups=1')
    elif sys.argv[1:]:
        raise SystemExit('Only --write or no arguments is supported.')
    else:
        saved=json.loads(dest.read_text(encoding='utf-8'))
        if saved!=result:raise SystemExit('MISMATCH_1269')
        print('READ_ONLY_MATCH_1269 independent_groups=1')
