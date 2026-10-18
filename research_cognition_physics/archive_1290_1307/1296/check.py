"""1296: exact budgets and finite self-dual CAR/source cross-checks.
Default and --check read only; --write exclusively creates results.json.
Finite matrices do not replace the analytical continuum theorem.
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
TOL=5e-11
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
JSP=1j*SY


def enc(v):
    if isinstance(v,F):return f'{v.numerator}/{v.denominator}'
    if isinstance(v,dict):return {str(k):enc(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [enc(x) for x in v]
    if isinstance(v,np.generic):return v.item()
    return v


def err(a,b=0):
    x=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    assert x<TOL,x
    return x


def unitary(h,t):
    w,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*w))@v.conj().T


def negative_and_absolute(h):
    w,v=np.linalg.eigh(h)
    assert min(abs(w))>1e-6
    return (v*(w<0))@v.conj().T,(v*abs(w))@v.conj().T


def blockdiag(a,b):
    return np.block([[a,np.zeros((a.shape[0],b.shape[1]),complex)],
                     [np.zeros((b.shape[0],a.shape[1]),complex),b]])


def bdg(h,d):
    assert err(d.T,-d)<TOL
    return np.block([[h,d],[d.conj().T,-h.conj()]])


def gamma(n):
    z=np.zeros((n,n));return np.block([[z,np.eye(n)],[np.eye(n),z]])


def exact_budget():
    e=7*ALPHA/256;b=ALPHA/256
    prob=ALPHA**2-2*ALPHA*e-b*b
    hs=(ALPHA-e)**2-b*b
    assert prob==F(61951,65536)*ALPHA**2>F(15,16)*ALPHA**2
    assert hs==F(3875,4096)*ALPHA**2
    assert 3*EPS+EPS/2==e
    assert 1+5*EPS+EPS<2
    assert 4*(524288+4)<1450**2
    return {'alpha_star':ALPHA,'epsilon':EPS,'delta_interval':[0,EPS],
            'time_window':[F(7,16),F(1,2)],'A_error_upper':e,'B_error_upper':b,
            'probability_strict_lower':prob,'probability_coefficient':F(61951,65536),
            'HS_nonzero_input_strict_lower_coefficient':F(3875,4096),
            'zero_input_only_nonstrict':True,'energy_square_upper':4*524292,
            'relative_energy_strict_upper':1450,'physical_modes':2,'Nambu_code_rank':4}


def yukawa(n=3):
    d=np.diag(np.arange(1,n+1,dtype=float))
    a=d+np.diag(np.ones(n-1),1)+np.diag(np.ones(n-1),-1)
    return np.eye(n)+float(EPS)*a,2*np.eye(n)+float(EPS)*d


def dirac(m,p):
    n=len(m);sp=p[0]*SX+p[1]*SY+p[2]*SZ
    k=np.kron(sp,np.eye(n));mm=np.kron(np.eye(2),m)
    return np.block([[-k,mm],[mm.conj().T,k]])


def continuous_symbol_diagnostics():
    n=3;m,_=yukawa(n);p=np.array([.3,-.2,.4])
    hp=dirac(m,p);hm=dirac(m,-p);h=blockdiag(hp,hm);dim=len(h)
    dr=np.zeros_like(hp);dr[2*n:,2*n:]=np.kron(JSP,np.eye(n))
    pair=np.kron(SX,dr)  # Local pairing connects opposite Fourier momenta.
    g=gamma(dim)
    v=np.zeros((dim,2),complex)
    for i in range(2):v[2*n+i,i]=1/math.sqrt(2);v[len(hp)+2*n+i,i]=1/math.sqrt(2)
    vf=np.vstack([v,np.zeros_like(v)]);vc=np.vstack([np.zeros_like(v),v.conj()])
    ff=vf@vf.conj().T;fc=vc@vc.conj().T;fh=ff+fc
    h0=bdg(h,np.zeros_like(h));maxerr=0.;rows=[]
    charge=blockdiag(-np.eye(dim),np.eye(dim))
    # The last value is a finite-matrix stress test, not an enlarged theorem domain.
    for delta in (0.,float(EPS),.125):
        hh=bdg(h,delta*pair);pref,absolute=negative_and_absolute(hh)
        maxerr=max(maxerr,err(hh,hh.conj().T),err(g@hh.conj()@g,-hh),
                   err(g@pref.conj()@g,np.eye(2*dim)-pref),err(g@fh.conj()@g,fh),
                   err(fh@hh@fh))
        comm=float(np.linalg.norm(hh@charge-charge@hh,2))
        maxerr=max(maxerr,err(comm,2*delta))
        expected=float((np.trace(fh@absolute)/2).real);energies=[]
        for rho in (np.diag([1.,0.]),np.array([[.6,.2j],[-.2j,.4]])):
            code=vf@rho@vf.conj().T+vc@(np.eye(2)-rho.T)@vc.conj().T
            q=np.eye(2*dim)-fh;state=q@pref@q+code;dc=state-pref
            energy=float((np.trace(hh@dc)/2).real)
            maxerr=max(maxerr,err(energy,expected),err(g@state.conj()@g,np.eye(2*dim)-state))
            energies.append(energy)
        assert expected>0
        u=unitary(hh,.5);u0=unitary(h0,.5)
        diff=float(np.linalg.norm(u-u0,2));anti=float(np.linalg.norm(u[:dim,dim:],2))
        assert diff<=.5*delta+2e-14 and anti<=.5*delta+2e-14
        maxerr=max(maxerr,err(g@u.conj()@g,u))
        rows.append({'delta':delta,'in_admitted_interval':delta<=float(EPS),
                     'B_minus_L_commutator_norm':comm,'Duhamel_norm':diff,
                     'Duhamel_upper':.5*delta,'anomalous_evolution_norm':anti,
                     'relative_energies':energies,'half_trace_Fhat_absolute_H':expected,
                     'omitting_Nambu_half_double_counts':2*expected})
    return {'Nambu_dimension':2*dim,'max_identity_residual':maxerr,'rows':rows,
            'opposite_momenta_retained':True,'not_continuum_tail_computation':True}


def background_ward():
    n=3;yn,ye=yukawa(n);z=np.zeros((2*n,2*n),complex)
    def mass(hp,h0):
        return np.block([[np.kron(np.eye(2),yn*np.conj(h0)),np.kron(np.eye(2),ye*hp)],
                         [-np.kron(np.eye(2),yn*np.conj(hp)),np.kron(np.eye(2),ye*h0)]])
    m=mass(0,1);zero=np.zeros_like(m)
    h=np.block([[zero,m],[m.conj().T,zero]]);dim=len(h)
    dunit=np.zeros_like(h,dtype=complex);dunit[4*n:6*n,4*n:6*n]=np.kron(JSP,np.eye(n))
    delta=.125;hh=bdg(h,delta*dunit)
    def variation(dhp=0,dh0=0,ds=0):
        dm=mass(dhp,dh0);dh=np.block([[zero,dm],[dm.conj().T,zero]])
        return bdg(dh,delta*ds*dunit)
    rows=[];maxerr=0.
    for name,th in [('W1',SX/2),('W2',SY/2),('W3',SZ/2),('Y',np.eye(2)/2)]:
        dh=1j*th@np.array([0,1.],complex)
        if name=='Y':
            tl=-.5*np.eye(4*n);tr=np.kron(np.diag([0.,-1.]),np.eye(2*n))
        else:
            tl=np.kron(th,np.eye(2*n));tr=np.zeros((4*n,4*n),complex)
        charge=bdg(blockdiag(tl,tr),np.zeros_like(h))
        defect=1j*(hh@charge-charge@hh)
        source=variation(dhp=dh[0],dh0=dh[1])
        er=err(defect+source);maxerr=max(maxerr,er)
        rows.append({'generator':name,'Ward_residual':er,
                     'drop_Higgs_source_defect':float(np.linalg.norm(defect,2))})
    bl=blockdiag(-np.eye(dim),np.eye(dim))
    defect=1j*(hh@bl-bl@hh);ssource=variation(ds=-2j)
    er=err(defect+ssource);maxerr=max(maxerr,er)
    assert np.linalg.norm(defect,2)>.2
    for hp,h0,ss in [(1,0,0),(1j,0,0),(0,1,0),(0,1j,0),(0,0,1),(0,0,1j)]:
        source=variation(hp,h0,ss);maxerr=max(maxerr,err(source,source.conj().T))
    return {'diagnostic_delta':delta,'all_four_H_and_two_S_real_sources':True,
            'electroweak_Wards':rows,'total_B_minus_L_Ward_residual':er,
            'drop_S_phase_source_defect':float(np.linalg.norm(defect,2)),
            'max_identity_residual':maxerr,'not_spacetime_differential_Ward_test':True,
            'fixed_couplings_no_arbitrary_Yukawa_spurion':True}


def annihilation(i,n):
    a=np.zeros((2**n,2**n),complex);bit=1<<(n-1-i)
    for c in range(2**n):
        if c&bit:a[c^bit,c]=(-1)**((c>>(n-i)).bit_count())
    return a


def trace_out(rho,dims,keep):
    a=rho.reshape(tuple(dims)*2);cur=list(dims)
    for i in reversed(range(len(cur))):
        if i not in keep:a=np.trace(a,axis1=i,axis2=i+len(cur));cur.pop(i)
    d=math.prod(cur);return a.reshape(d,d)


def covariance(rho,ops):
    return np.array([[np.trace(rho@ops[j].conj().T@ops[i]) for j in range(len(ops))]
                     for i in range(len(ops))])


def quadratic(h,d,ops):
    out=np.zeros_like(ops[0]);pair=np.zeros_like(out)
    for i in range(len(ops)):
        for j in range(len(ops)):
            out+=h[i,j]*ops[i].conj().T@ops[j]
            pair+=.5*d[i,j]*ops[i].conj().T@ops[j].conj().T
    return out+pair+pair.conj().T-.5*np.trace(h)*np.eye(len(out))


def full_CAR_and_partner():
    ops=[annihilation(i,4) for i in range(4)];nops=ops+[x.conj().T for x in ops]
    h=np.array([[-2,.2,.7,.1j],[.2,-1,.3j,.6],[.7,-.3j,1,.4],[-.1j,.6,.4,2]],complex)
    d=np.zeros((4,4),complex)
    for i,j,v in [(0,1,.1),(0,2,.35),(1,3,.12),(2,3,.25j)]:d[i,j]=v;d[j,i]=-v
    hh=bdg(h,d);H=quadratic(h,d,ops);vals,vec=np.linalg.eigh(H)
    psi=vec[:,0];ref=np.outer(psi,psi.conj());pref,_=negative_and_absolute(hh)
    maxerr=err(covariance(ref,nops),pref)
    env=trace_out(ref,[4,4],[1]);code_ref=trace_out(ref,[4,4],[0])
    factorization_defect=float(np.linalg.norm(ref-np.kron(code_ref,env),2))
    assert factorization_defect>.05
    envops=[annihilation(i,2) for i in range(2)]
    anomaly=complex(np.trace(env@envops[0]@envops[1]));assert abs(anomaly)>1e-3
    w=np.zeros((4,2),complex);w[2,0]=1;w[1,1]=1
    rho=np.array([[.6,.2j],[-.2j,.4]],complex)
    state=np.kron(w@rho@w.conj().T,env)
    vf=np.eye(8)[:,:2];vc=np.eye(8)[:,4:6]
    ff=vf@vf.conj().T;fc=vc@vc.conj().T;fh=ff+fc;q=np.eye(8)-fh
    code=vf@rho@vf.conj().T+vc@(np.eye(2)-rho.T)@vc.conj().T
    maxerr=max(maxerr,err(covariance(state,nops),q@pref@q+code))
    n01=ops[0].conj().T@ops[0]@ops[1].conj().T@ops[1]
    actual_double=float(np.trace(state@n01).real);assert abs(actual_double)<TOL
    wrong_wick=float(np.linalg.det(rho).real);assert wrong_wick>.1
    # The complete state, not only its two-point covariance, carries the reference partner.
    bell=np.array([1,0,0,1],complex)/math.sqrt(2);qr=np.outer(bell,bell.conj())
    wr=np.kron(w,np.eye(2));fr=wr@qr@wr.conj().T
    full=np.kron(fr,env).reshape(4,2,4,4,2,4).transpose(0,2,1,3,5,4).reshape(32,32)
    maxerr=max(maxerr,err(trace_out(full,[4,4,2],[1,2]),np.kron(env,np.eye(2)/2)),
               err(trace_out(full,[4,4,2],[0,2]),fr))
    t=.4;u=unitary(H,t);ub=unitary(hh,t)
    maxerr=max(maxerr,err(covariance(u@state@u.conj().T,nops),ub@covariance(state,nops)@ub.conj().T))
    ar=ub[2:4,:2];br=ub[2:4,4:6]
    X=np.diag([1.,-1.]);theory=ar@X@ar.conj().T-br@X.T@br.conj().T
    s0=np.kron(np.outer(w[:,0],w[:,0].conj()),env)
    s1=np.kron(np.outer(w[:,1],w[:,1].conj()),env)
    diff=covariance(u@(s0-s1)@u.conj().T,nops)[2:4,2:4]
    maxerr=max(maxerr,err(diff,theory))
    omit_B=float(np.linalg.norm(diff-ar@X@ar.conj().T,2));assert omit_B>1e-4
    full_u=np.kron(u,np.eye(2));out=full_u@full@full_u.conj().T
    maxerr=max(maxerr,err(trace_out(out,[4,4,2],[2]),np.eye(2)/2))
    # General normal + anomalous source, with reference baseline and Nambu half.
    j=np.array([[1,.2j,.3,0],[-.2j,-.4,.1,.2j],[.3,.1,.2,-.5],[0,-.2j,-.5,.8]],complex)
    dj=np.zeros((4,4),complex)
    for i,k,v in [(0,2,.3j),(1,3,.2),(2,3,.4j)]:dj[i,k]=v;dj[k,i]=-v
    js=bdg(j,dj);J=quadratic(j,dj,ops)
    common=-fh@pref-pref@fh+fh@pref@fh+fc
    c=np.trace(ref@J)+.5*np.trace(js@common);bm=vf.conj().T@js@vf
    maxerr=max(maxerr,err(np.trace(state@J),c+np.trace(rho@bm)),
               err(np.trace(state@J),.5*np.trace(js@covariance(state,nops))))
    corr=np.trace(full@np.kron(J,SZ))
    rhs=c*np.trace(qr@np.kron(np.eye(2),SZ))+np.trace(qr@np.kron(bm,SZ))
    maxerr=max(maxerr,err(corr,rhs));assert abs(corr)>.1
    no_anomaly=np.array(covariance(state,nops),copy=True);no_anomaly[:4,4:]=0;no_anomaly[4:,:4]=0
    omitted_env=float(abs(np.trace(state@J)-.5*np.trace(js@no_anomaly)));assert omitted_env>1e-4
    return {'physical_modes':4,'Fock_dimension':16,'max_identity_residual':maxerr,
            'reference_factorization_defect':factorization_defect,
            'environment_anomalous_covariance_magnitude':float(abs(anomaly)),
            'actual_double_occupation':actual_double,'wrong_Gaussian_double_occupation':wrong_wick,
            'remote_R_product_at_initial_slice':True,'R_marginal_after_full_evolution':True,
            'omit_anomalous_transfer_leg_defect':omit_B,
            'omit_environment_anomalous_source_defect':omitted_env,
            'full_source_R_correlation_real':float(corr.real),
            'difference_map_not_claimed_CP':True,'finite_diagnostic_not_new_physical_parameters':True}


def bindings():
    paths={'proof':'proof.md','check':'check.py','admission':'../_shared/round1296_admission.md',
           '1295proof':'../1295/proof.md','1295results':'../1295/results.json',
           '1293proof':'../1293/proof.md','1293results':'../1293/results.json',
           '1292proof':'../1292/proof.md','730':'../../archive_702_741/research_note_730.md',
           '805':'../../archive_797_818/research_note_805.md','813':'../../archive_797_818/research_note_813.md',
           '818':'../../archive_797_818/research_note_818.md','1060':'../../archive_1046_1062/1060/proof.md',
           '1034':'../../archive_1009_1043/research_note_1034.md'}
    return {k:{'path':p,'sha256':hashlib.sha256((HERE/p).read_bytes()).hexdigest()} for k,p in paths.items()}


def compare(a,b):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isfinite(a) and math.isfinite(b)
        assert math.isclose(a,b,rel_tol=2e-8,abs_tol=2e-14),(a,b)
    else:assert a==b,(a,b)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.check and args.write:parser.error('choose --check or --write')
    payload=enc({'round':1296,'author_check_groups':1,'verdict':'PASS_AUTHOR_SCOPED',
                 'exact_budget':exact_budget(),'BdG_energy_and_evolution':continuous_symbol_diagnostics(),
                 'same_action_background_Wards':background_ward(),'complete_CAR_and_partner':full_CAR_and_partner(),
                 'scope':{'own_selfdual_Hadamard_reference_per_model':True,'code_actual_not_Gaussian':True,
                          'full_source_includes_g_A_H_S':True,'Nambu_half_retained':True,
                          'finite_effect_bound_controls_unbounded_source':False,
                          'difference_map_itself_is_CP':False,'fermionic_B_minus_L_vs_total_charge_distinguished':True,
                          'fixed_nonzero_S_background_r_broken_in_both_branches':True,
                          'full_interacting_SM_certified':False,'actual_readout_apparatus_certified':False,
                          'all_cognitive_roles_independence_proved':False},'sources':bindings()})
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:json.dump(payload,f,ensure_ascii=False,indent=2,sort_keys=True);f.write('\n')
        print('PASS_AUTHOR_SCOPED EXCLUSIVE_WRITE')
    else:
        compare(payload,json.loads(OUT.read_text(encoding='utf-8')));print('PASS_AUTHOR_SCOPED READ_ONLY_MATCH')
    print('n=3; 0<=delta<=alpha*/128; gap>61951 alpha*^2/65536; HS>3875 alpha*^2/4096 for X!=0')
    print('Gamma-invariant rank-4 code; Nambu half-energy<1450; all fixed-background sources retained')


if __name__=='__main__':main()
