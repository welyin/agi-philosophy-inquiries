"""1297: gauge-orbit transport of the fixed 1296 regional record contract.
Default/--check is read only. --write exclusively creates results.json.
Finite diagnostics do not recalculate the inherited continuum signal theorem.
"""
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'
ALPHA=F(31,2**27)
EPS=ALPHA/128
DELTA=EPS/2
TOL=5e-11
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
JSP=1j*SY


def enc(x):
    if isinstance(x,F):return f'{x.numerator}/{x.denominator}'
    if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v) for v in x]
    if isinstance(x,np.generic):return x.item()
    return x


def err(a,b=0):
    e=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    assert e<TOL,e
    return e


def bd(a,b):
    return np.block([[a,np.zeros((len(a),len(b)),complex)],
                     [np.zeros((len(b),len(a)),complex),b]])


def bdg(h,d):return np.block([[h,d],[d.conj().T,-h.conj()]])
def gamma(n):return np.block([[np.zeros((n,n)),np.eye(n)],[np.eye(n),np.zeros((n,n))]])

def unitary(h,t):
    w,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*w))@v.conj().T


def negative(h):
    w,v=np.linalg.eigh(h);assert min(abs(w))>1e-7
    return (v*(w<0))@v.conj().T


def exact_inheritance_and_charges():
    e=7*ALPHA/256;b=ALPHA/256
    gap=ALPHA**2-2*ALPHA*e-b*b
    hs=(ALPHA-e)**2-b*b
    assert DELTA==F(31,2**35)
    assert 3*EPS+DELTA/2<=e and DELTA/2<=b
    assert gap==F(61951,65536)*ALPHA**2>F(15,16)*ALPHA**2
    assert hs==F(3875,4096)*ALPHA**2
    assert 4*(524288+4)<1450**2
    names=['Q','uc','dc','L','ec','N','H','S']
    y0=list(map(F,[F(1,6),F(-2,3),F(1,3),F(-1,2),1,0,F(1,2),0]))
    bl=list(map(F,[F(1,3),F(-1,3),F(-1,3),-1,1,1,0,-2]))
    trial=[1,-1,-1,0,0,0,0,0];weak=[1,0,0,1,0,0,1,0]
    rows=[]
    for k in (0,1):
        y=[x+3*k*z for x,z in zip(y0,bl)]
        assert y[3]+y[6]==-3*k and -y[5]==-3*k and y[7]==-6*k
        assert y[7]+2*y[5]==0
        assert y[3]+y[6]+y[5]==0 and y[3]-y[6]+y[4]==0
        center=[F(t,3)+F(w,2)+yy for t,w,yy in zip(trial,weak,y)]
        assert all(x.denominator==1 for x in center)
        rows.append({'k':k,'hypercharges':dict(zip(names,y)),'G6_center_exponents':center,
                     'chi_hypercharge':-3*k,'S_hypercharge':-6*k,
                     'scalar_hypercharge_contact_second_derivative':F(1,2)+72*k*k})
    return {'alpha_star':ALPHA,'epsilon':EPS,'fixed_delta':DELTA,
            'time_window':[F(7,16),F(1,2)],'probability_strict_lower':gap,
            'HS_nonzero_input_strict_lower_coefficient':F(3875,4096),
            'relative_energy_strict_upper':1450,'charge_branches':rows,
            'constants_inherited_not_new_continuum_computation':True}


def wilson_and_time_transport():
    x0=np.array([.1,-.2,.15]);x=np.array([.13,-.21,.14]);d=x-x0
    c=np.array([.2,-.4,.3]);M=np.array([[.3,.2,-.1],[-.3,.4,.2],[.1,.5,-.2]])
    l=np.array([.4,.2,-.3]);Q=np.array([[.2,.1,.05],[.1,-.3,.04],[.05,.04,.1]])
    phase=lambda y:float(l@y+y@Q@y/2)
    integral=float(c@d+d@M@x0+d@M@d/2)
    shifted=float((c+l)@d+d@(M+Q)@x0+d@(M+Q)@d/2)
    maximum=err(shifted-integral,phase(x)-phase(x0))
    h=np.array([[-1.4,.2,.5,.1j],[.2,-.7,.2j,.4],[.5,-.2j,1.1,.2],[-.1j,.4,.2,1.8]],complex)
    pair=np.array([[0,.1,.3,0],[-.1,0,0,.2],[-.3,0,0,.1j],[0,-.2,-.1j,0]],complex)
    HH=bdg(h,pair);P=negative(HH);U=unitary(HH,.4);Gam=gamma(4)
    angle0=np.array([.13,.13,-.2,-.2]);angle1=np.array([.4,.4,.17,.17])
    rate=np.array([.8,.8,-.6,-.6]);rows=[]
    for k in (0,1):
        q=-3*k;W=np.exp(-1j*q*integral);Wp=np.exp(-1j*q*shifted)
        maximum=max(maximum,err(Wp*np.exp(1j*q*phase(x)),np.exp(1j*q*phase(x0))*W))
        g0=np.diag(np.exp(1j*q*angle0));g1=np.diag(np.exp(1j*q*angle1))
        G0=bd(g0,g0.conj());G1=bd(g1,g1.conj());Dq=bd(np.diag(q*rate),np.diag(-q*rate))
        HHp=G1@HH@G1.conj().T-Dq;Pp=G1@P@G1.conj().T;Up=G1@U@G0.conj().T
        maximum=max(maximum,err(Gam@G1.conj()@Gam,G1),err(Gam@HHp.conj()@Gam,-HHp),
                    err(Gam@Pp.conj()@Gam,np.eye(8)-Pp))
        # Exact derivative of G(t) U(t,0) G(0)^dagger, with the time connection.
        derivative=1j*Dq@Up+G1@(-1j*HH@U)@G0.conj().T
        maximum=max(maximum,err(1j*derivative,HHp@Up))
        A=U[2:4,:2];B=U[2:4,4:6];Ap=Up[2:4,:2];Bp=Up[2:4,4:6]
        maximum=max(maximum,err(Ap,np.exp(1j*q*(angle1[2]-angle0[0]))*A),
                    err(Bp,np.exp(1j*q*(angle1[2]+angle0[0]))*B))
        X=np.array([[.4,.2j],[-.2j,-.4]])
        L=A@X@A.conj().T-B@X.T@B.conj().T
        Lp=Ap@X@Ap.conj().T-Bp@X.T@Bp.conj().T
        maximum=max(maximum,err(Lp,L))
        instantaneous=float(np.linalg.norm(Pp-negative(HHp),2))
        missing_time=float(np.linalg.norm(1j*derivative-G1@HH@G1.conj().T@Up,2))
        if k:assert instantaneous>.05 and missing_time>1
        rows.append({'k':k,'instantaneous_negative_projection_error':instantaneous,
                     'omitted_time_connection_defect':missing_time,
                     'normal_and_anomalous_endpoint_phases_checked':True})
    return {'max_identity_residual':maximum,'rows':rows,
            'finite_matrix_orbit_diagnostic_not_new_spacetime_model':True}


def source_ward():
    n=3;D=np.diag(np.arange(1,n+1));path=np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
    yn=np.eye(n)+float(EPS)*(D+path);ye=2*np.eye(n)+float(EPS)*D
    def mass(hp,h0):
        return np.block([[np.kron(np.eye(2),yn*np.conj(h0)),np.kron(np.eye(2),ye*hp)],
                         [-np.kron(np.eye(2),yn*np.conj(hp)),np.kron(np.eye(2),ye*h0)]])
    m=mass(0,1);z=np.zeros_like(m);h=np.block([[z,m],[m.conj().T,z]])
    dr=np.zeros_like(h,dtype=complex);dr[4*n:6*n,4*n:6*n]=np.kron(JSP,np.eye(n))
    dd=float(DELTA);HH=bdg(h,dd*dr);maximum=0.;rows=[]
    def var(hp=0,h0=0,s=0):
        dm=mass(hp,h0);dh=np.block([[z,dm],[dm.conj().T,z]])
        return bdg(dh,dd*s*dr)
    for k in (0,1):
        yl=(-.5-3*k)*np.eye(4*n);yr=np.kron(np.diag([-3*k,-1-3*k]),np.eye(2*n))
        QY=bdg(bd(yl,yr),np.zeros_like(h))
        Hsource=var(h0=.5j);Ssource=var(s=-6j*k)
        defect=1j*(HH@QY-QY@HH)
        maximum=max(maximum,err(defect+Hsource+Ssource))
        omitS=float(np.linalg.norm(defect+Hsource,2))/dd
        if k:assert math.isclose(omitS,6.,rel_tol=1e-7)
        for dH,d0,dS in [(1,0,0),(1j,0,0),(0,1,0),(0,1j,0),(0,0,1),(0,0,1j)]:
            J=var(dH,d0,dS);maximum=max(maximum,err(J,J.conj().T),err(gamma(len(h))@J.conj()@gamma(len(h)),-J))
        # Finite gauge transformation of the original Yukawa and SNN mass matrix.
        theta=.19
        gp=np.diag(np.exp(1j*theta*np.diag(bd(yl,yr))))
        GG=bd(gp,gp.conj());mp=mass(0,np.exp(.5j*theta))
        hp=np.block([[z,mp],[mp.conj().T,z]])
        direct=bdg(hp,dd*np.exp(-6j*k*theta)*dr)
        maximum=max(maximum,err(direct,GG@HH@GG.conj().T))
        rows.append({'k':k,'hypercharge_Ward_residual':float(np.max(abs(defect+Hsource+Ssource))),
                     'omit_S_phase_defect_divided_by_fixed_delta':omitS,
                     'absolute_source_means_not_asserted_distinct':True})
    # All three weak directions and Higgs phase were varied before the background.
    for T in [SX/2,SY/2,SZ/2]:
        dH=1j*T@np.array([0,1.],complex)
        Q=bdg(bd(np.kron(T,np.eye(2*n)),np.zeros((4*n,4*n))),np.zeros_like(h))
        maximum=max(maximum,err(1j*(HH@Q-Q@HH)+var(dH[0],dH[1])))
    return {'fixed_delta':DELTA,'max_identity_residual':maximum,'rows':rows,
            'all_four_H_two_S_and_weak_sources_retained':True,
            'g_source_and_general_insertions_use_covariant_pullback':True,
            'not_a_full_differential_spacetime_Ward_simulation':True}


def annihilation(i,n):
    a=np.zeros((2**n,2**n),complex);bit=1<<(n-1-i)
    for j in range(2**n):
        if j&bit:a[j^bit,j]=(-1)**((j>>(n-i)).bit_count())
    return a


def trace_out(rho,dims,keep):
    a=rho.reshape(tuple(dims)*2);cur=list(dims)
    for i in reversed(range(len(cur))):
        if i not in keep:a=np.trace(a,axis1=i,axis2=i+len(cur));cur.pop(i)
    return a.reshape(math.prod(cur),math.prod(cur))


def covariance(rho,ops):
    return np.array([[np.trace(rho@ops[j].conj().T@ops[i]) for j in range(len(ops))] for i in range(len(ops))])


def quadratic(h,d,ops):
    H=np.zeros_like(ops[0]);P=np.zeros_like(H)
    for i in range(len(ops)):
        for j in range(len(ops)):
            H+=h[i,j]*ops[i].conj().T@ops[j]
            P+=d[i,j]*ops[i].conj().T@ops[j].conj().T/2
    return H+P+P.conj().T-.5*np.trace(h)*np.eye(len(H))


def complete_code_reference_and_source():
    ops=[annihilation(i,4) for i in range(4)];nops=ops+[a.conj().T for a in ops]
    h=np.array([[-2,.2,.7,.1j],[.2,-1,.3j,.6],[.7,-.3j,1,.4],[-.1j,.6,.4,2]],complex)
    d=np.zeros((4,4),complex)
    for i,j,v in [(0,1,.1),(0,2,.35),(1,3,.12),(2,3,.25j)]:d[i,j]=v;d[j,i]=-v
    H=quadratic(h,d,ops);ev,V=np.linalg.eigh(H);vac=V[:,0];ref=np.outer(vac,vac.conj())
    P=negative(bdg(h,d));maximum=err(covariance(ref,nops),P)
    env=trace_out(ref,[4,4],[1]);rf=trace_out(ref,[4,4],[0])
    nonfactor=float(np.linalg.norm(ref-np.kron(rf,env),2));assert nonfactor>.05
    eo=[annihilation(i,2) for i in range(2)]
    anomaly=abs(np.trace(env@eo[0]@eo[1]));assert anomaly>1e-3
    W=np.zeros((4,2),complex);W[2,0]=1;W[1,1]=1
    rho=np.array([[.6,.2j],[-.2j,.4]])
    state=np.kron(W@rho@W.conj().T,env)
    bell=np.array([1,0,0,1],complex)/math.sqrt(2);QR=np.outer(bell,bell.conj())
    WR=np.kron(W,np.eye(2));FR=WR@QR@WR.conj().T
    total=np.kron(FR,env).reshape(4,2,4,4,2,4).transpose(0,2,1,3,5,4).reshape(32,32)
    maximum=max(maximum,err(trace_out(total,[4,4,2],[1,2]),np.kron(env,np.eye(2)/2)))
    Number=sum(a.conj().T@a for a in ops);Nf=eo[0].conj().T@eo[0]+eo[1].conj().T@eo[1]
    charge_variance=float((np.trace(state@Number@Number)-np.trace(state@Number)**2).real)
    assert charge_variance>1e-3
    reset=np.outer(W[:,0],np.array([1.,0,0,0]));assert np.linalg.norm(Nf@reset-reset@Nf)>0.9
    double=float(np.trace(state@ops[0].conj().T@ops[0]@ops[1].conj().T@ops[1]).real)
    assert abs(double)<TOL and np.linalg.det(rho).real>.1
    js=np.array([[1,.2j,.3,0],[-.2j,-.4,.1,.2j],[.3,.1,.2,-.5],[0,-.2j,-.5,.8]],complex)
    ds=np.zeros((4,4),complex)
    for i,j,z in [(0,2,.3j),(1,3,.2),(2,3,.4j)]:ds[i,j]=z;ds[j,i]=-z
    J=quadratic(js,ds,ops);jj=bdg(js,ds)
    vf=np.eye(8)[:,:2];vc=np.eye(8)[:,4:6];fh=vf@vf.conj().T+vc@vc.conj().T
    fc=vc@vc.conj().T
    c=np.trace(ref@J)+.5*np.trace(jj@(-fh@P-P@fh+fh@P@fh+fc));jm=vf.conj().T@jj@vf
    maximum=max(maximum,err(np.trace(state@J),c+np.trace(rho@jm)))
    corr=np.trace(total@np.kron(J,SZ));rhs=c*np.trace(QR@np.kron(np.eye(2),SZ))+np.trace(QR@np.kron(jm,SZ))
    maximum=max(maximum,err(corr,rhs));assert abs(corr)>.1
    U=unitary(H,.4);rows=[]
    a0=np.array([.11,.11,-.23,-.23]);a1=np.array([.32,.32,.07,.07])
    for k in (0,1):
        q=-3*k;V0=unitary(-sum(q*a0[i]*ops[i].conj().T@ops[i] for i in range(4)),1)
        V1=unitary(-sum(q*a1[i]*ops[i].conj().T@ops[i] for i in range(4)),1)
        codephase=np.diag(np.exp(1j*q*a0[0]*np.diag(Nf)))
        maximum=max(maximum,err(codephase@(W@rho@W.conj().T)@codephase.conj().T,W@rho@W.conj().T))
        initial=np.kron(V0,np.eye(2))@total@np.kron(V0.conj().T,np.eye(2))
        envphase=unitary(-q*a0[2]*Nf,1);envp=envphase@env@envphase.conj().T
        maximum=max(maximum,err(trace_out(initial,[4,4,2],[1,2]),np.kron(envp,np.eye(2)/2)))
        Uprime=V1@U@V0.conj().T
        result=np.kron(Uprime,np.eye(2))@initial@np.kron(Uprime.conj().T,np.eye(2))
        original=np.kron(U,np.eye(2))@total@np.kron(U.conj().T,np.eye(2))
        maximum=max(maximum,err(result,np.kron(V1,np.eye(2))@original@np.kron(V1.conj().T,np.eye(2))),
                    err(trace_out(result,[4,4,2],[2]),np.eye(2)/2),
                    err(np.trace(result@np.kron(V1@J@V1.conj().T,SZ)),np.trace(original@np.kron(J,SZ))))
        for i in (2,3):
            effect=ops[i].conj().T@ops[i]
            maximum=max(maximum,err(np.trace(result@np.kron(effect,SZ)),np.trace(original@np.kron(effect,SZ))))
        rows.append({'k':k,'full_charge_variance':q*q*charge_variance,
                     'initial_remote_R_product_and_full_output_transport':True})
    return {'Fock_dimension':16,'physical_modes':4,'max_identity_residual':maximum,
            'reference_factorization_defect':nonfactor,'environment_anomaly':float(anomaly),
            'actual_double_occupation':double,'incorrect_Gaussian_double_occupation':float(np.linalg.det(rho).real),
            'source_R_correlation_real':float(corr.real),'reset_Kraus_charge_change_nonzero':True,
            'rows':rows,'finite_algebra_diagnostic_not_new_continuum_parameters':True}


def bindings():
    paths={'proof':'proof.md','check':'check.py','admission':'../_shared/round1297_admission.md',
           '1296proof':'../1296/proof.md','1296results':'../1296/results.json',
           '1295proof':'../1295/proof.md','1293proof':'../1293/proof.md','1292proof':'../1292/proof.md',
           'background_dressing':'../_admission/background_dressing_screen_after1293.md',
           'physical_subsystems':'../../archive_1046_1062/_admission/physical_subsystems_after1055/selection.md',
           '1060':'../../archive_1046_1062/1060/proof.md','730':'../../archive_702_741/research_note_730.md',
           '805':'../../archive_797_818/research_note_805.md','813':'../../archive_797_818/research_note_813.md',
           '818':'../../archive_797_818/research_note_818.md'}
    return {k:{'path':p,'sha256':hashlib.sha256((HERE/p).read_bytes()).hexdigest()} for k,p in paths.items()}


def compare(a,b):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isfinite(a) and math.isfinite(b) and math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-14),(a,b)
    else:assert a==b,(a,b)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.write and args.check:parser.error('choose --write or --check')
    payload=enc({'round':1297,'author_check_groups':1,'verdict':'PASS_AUTHOR_SCOPED',
                 'fixed_contract':exact_inheritance_and_charges(),'Wilson_and_Nambu_transport':wilson_and_time_transport(),
                 'same_action_sources':source_ward(),'complete_CAR_and_partner':complete_code_reference_and_source(),
                 'scope':{'k_branches':[0,1],'fixed_delta_no_scan':True,'background_orbit_only':True,
                          'instantaneous_negative_reference_reselection':False,'full_g_A_H_S_source_identity':True,
                          'mode_background_derivative_is_not_physical_source':True,
                          'source_matrices_differ_not_every_state_mean':True,'code_charge_is_not_whole_field_charge':True,
                          'normal_anomalous_environment_retained':True,'independent_partner_full_state':True,
                          'dynamic_quantum_Gauss_certified':False,'autonomous_preparation_certified':False,
                          'actual_terminal_readout_device_certified':False,'full_interacting_SM_certified':False,
                          'all_k_uniform_source_noise_budget':False},'sources':bindings()})
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:json.dump(payload,f,ensure_ascii=False,sort_keys=True,indent=2);f.write('\n')
        print('PASS_AUTHOR_SCOPED EXCLUSIVE_WRITE')
    else:
        compare(payload,json.loads(OUT.read_text(encoding='utf-8')));print('PASS_AUTHOR_SCOPED READ_ONLY_MATCH')
    print('k=0,1; delta=epsilon/2; inherited whole-window gap>61951 alpha_star^2/65536; Erel<1450')
    print('gauge-orbit transport, both anomalous endpoints, complete Q/R state and all-source identity checked')


if __name__=='__main__':main()
