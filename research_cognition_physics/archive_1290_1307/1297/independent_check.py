"""1297 independent finite diagnostics and exact inherited certificates.
No import of author or historical calculation code.  Continuous field claims
are analytical: finite matrices below test the dictionary, not a lattice QFT.
Default/--check is read-only; --write exclusively creates independent_results.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if (p/'research_cognition_physics').is_dir())
OUT=HERE/'independent_results.json'
ALPHA=F(31,2**27); EPS=ALPHA/128

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def residual(a,b): return float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
def norm(a): return float(np.linalg.norm(a,2))
def sp(h,fn):
    x,v=np.linalg.eigh(h)
    return (v*fn(x))@v.conj().T
def unitary(h,t): return sp(h,lambda x:np.exp(-1j*t*x))
def bdg(h,d):
    assert residual(h,h.conj().T)<1e-14 and residual(d,-d.T)<1e-14
    return np.block([[h,d],[-d.conj(),-h.conj()]])
def diag_blocks(a,b):
    z=np.zeros((len(a),len(b)),complex)
    return np.block([[a,z],[z.T,b]])
def conjugator(m):
    z=np.zeros((m,m));i=np.eye(m)
    return np.block([[z,i],[i,z]])
def lift(g): return diag_blocks(g,g.conj())
def clean(x):
    if isinstance(x,F):return {'numerator':x.numerator,'denominator':x.denominator}
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(np.floating,)):return float(x)
    return x

def exact_data():
    records=[]
    for k in (0,1):
        charges={'Q':F(1,6)+k,'uc':-F(2,3)-k,'dc':F(1,3)-k,
                 'L':-F(1,2)-3*k,'ec':1+3*k,'N':3*k,'H':F(1,2),'S':-6*k}
        assert charges['L']+charges['H']+charges['N']==0
        assert charges['L']-charges['H']+charges['ec']==0
        assert charges['S']+2*charges['N']==0
        assert charges['Q']+charges['H']+charges['uc']==0
        assert charges['Q']-charges['H']+charges['dc']==0
        triality={'Q':1,'uc':-1,'dc':-1,'L':0,'ec':0,'N':0,'H':0,'S':0}
        weak={'Q':1,'uc':0,'dc':0,'L':1,'ec':0,'N':0,'H':1,'S':0}
        center={n:F(triality[n],3)+F(weak[n],2)+v for n,v in charges.items()}
        assert all(x.denominator==1 for x in center.values())
        q=charges['L']+charges['H']
        assert q==-charges['N']==-3*k and charges['S']==2*q
        records.append(dict(k=k,hypercharges=charges,G6_center_exponents=center,
                            physical_two_track_charge=q,integer_cover_line_charge=6*q))
    e=7*ALPHA/256;b=ALPHA/256
    assert ALPHA**2-2*ALPHA*e-b*b==F(61951,65536)*ALPHA**2
    assert (ALPHA-e)**2-b*b==F(3875,4096)*ALPHA**2
    assert 4*(524288+4)<1450**2
    return dict(alpha=ALPHA,epsilon=EPS,delta=EPS/2,
       whole_time_window=[F(7,16),F(1,2)],probability_lower=F(61951,65536)*ALPHA**2,
       HS_lower=F(3875,4096)*ALPHA**2,relative_energy_upper=1450,
       original_continuous_proof_not_recomputed=True,charges=records)

def wilson_tests():
    # Polynomial connection and gauge parameter; line integral is exact rational.
    def alpha(x):return F(1,7)+2*x[0]-3*x[1]+x[0]*x[1]+2*x[2]**2
    def integral(x0,x):
        d=[x[j]-x0[j] for j in range(3)]
        # a=(y,2z,-x) on the straight radial path.
        return d[0]*(x0[1]+d[1]/2)+2*d[1]*(x0[2]+d[2]/2)-d[2]*(x0[0]+d[0]/2)
    centers=[(F(0),F(0),F(0)),(F(1,4),F(0),F(0))]
    offsets=[(F(1,128),F(1,160),F(0)),(-F(1,128),F(0),F(1,192)),
             (F(0),-F(1,160),-F(1,192))]
    radius=F(1,32);checks=[];mx=0.
    for k in (0,1):
        q=-3*k
        for base in centers:
            points=[tuple(base[j]+o[j] for j in range(3)) for o in offsets]
            assert all(sum(o[j]**2 for j in range(3))<radius**2 for o in offsets)
            for x in points:
                a=integral(base,x); da=alpha(x)-alpha(base)
                lhs=np.exp(-1j*q*float(a+da))*np.exp(1j*q*float(alpha(x)))
                rhs=np.exp(1j*q*float(alpha(base)))*np.exp(-1j*q*float(a))
                mx=max(mx,float(abs(lhs-rhs)))
            # A spatially varying mode projection must be transported, not fixed.
            weight=np.sqrt(np.array([.2,.3,.5]))
            v=np.kron(weight[:,None],np.eye(2))
            g=np.diag(np.repeat([np.exp(1j*q*float(alpha(x))) for x in points],2))
            ff=v@v.conj().T;vg=g@v; fg=vg@vg.conj().T
            mx=max(mx,residual(fg,g@ff@g.conj().T),residual(fg@fg,fg))
            checks.append(dict(k=k,base=[str(y) for y in base],
               original_projection_change=norm(fg-ff)))
    # Higgs dressing: invariant SU(2) epsilon tensor and hypercharge phase.
    eps2=np.array([[0,1],[-1,0]],complex)
    a=.4;axis=np.array([.2,-.3,.6]);axis/=np.linalg.norm(axis)
    sig=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    su=math.cos(a)*np.eye(2)+1j*math.sin(a)*sum(axis[j]*sig[j] for j in range(3))
    H=np.array([.2+.1j,1.1-.3j]); L=np.array([.7-.2j,.1+.6j]);angle=.37
    dress=lambda hh,ll: (eps2@hh.conj()).conj()@ll/np.linalg.norm(hh)
    for k in (0,1):
        hp=np.exp(1j*angle/2)*su@H
        lp=np.exp(1j*(-.5-3*k)*angle)*su@L
        mx=max(mx,float(abs(dress(hp,lp)-np.exp(-3j*k*angle)*dress(H,L))))
    assert mx<2e-13
    return dict(max_residual=mx,radial_trees_stay_in_declared_balls=True,
       task_support_not_point_local=True,projection_diagnostics=checks,
       no_quantized_Wilson_carrier=True)

def car(m):
    aa=[]
    for j in range(m):
        mat=np.zeros((2**m,2**m),complex)
        bit=1<<(m-1-j)
        for c in range(2**m):
            if c&bit:
                sign=(-1)**sum((c>>(m-1-l))&1 for l in range(j))
                mat[c^bit,c]=sign
        aa.append(mat)
    return aa

def quantize(j,aa):
    zz=aa+[a.conj().T for a in aa]
    out=np.zeros_like(aa[0])
    for r in range(len(zz)):
        for s in range(len(zz)):out+=j[r,s]*(zz[r].conj().T@zz[s])/2
    return out

def cov(state,aa):
    zz=aa+[a.conj().T for a in aa]
    return np.array([[np.trace(state@zz[s].conj().T@zz[r]) for s in range(len(zz))]
                     for r in range(len(zz))])

def finite_dictionary():
    # Deliberately different finite coefficients from the continuum / author.
    h=np.array([[0,0,.7,.2j],[0,0,-.3,.9],[.7,-.3,.4,.13j],[-.2j,.9,-.13j,-.2]],complex)
    dd=np.zeros((4,4),complex)
    for i,j,v in [(0,2,.11),(1,3,-.07j),(2,3,.23)]:dd[i,j]=v;dd[j,i]=-v
    hh=bdg(h,dd);cc=conjugator(4)
    pp=sp(hh,lambda x:(x<0).astype(float));aa=car(4)
    ham=quantize(hh,aa);e,vec=np.linalg.eigh(ham);vac=vec[:,0]
    reference=np.outer(vac,vac.conj())
    env=np.einsum('aiaj->ij',reference.reshape(4,4,4,4))
    env_anom=norm(cov(env,car(2))[:2,2:])
    assert env_anom>.01
    ff=np.diag([1,1,0,0]);fh=diag_blocks(ff,ff);qc=np.eye(8)-fh
    fc=diag_blocks(np.zeros((4,4)),ff)
    w=np.zeros((4,2),complex);w[2,0]=1;w[1,1]=1
    psiQR=np.array([math.sqrt(.41),0,0,1j*math.sqrt(.59)])
    rhoQR=np.outer(psiQR,psiQR.conj())
    qr4=rhoQR.reshape(2,2,2,2)
    rQ=np.einsum('aibi->ab',qr4);rR=np.einsum('aiaj->ij',qr4)
    codeunits=[[np.kron(np.outer(w[:,i],w[:,j].conj()),env) for j in range(2)] for i in range(2)]
    joint=sum((np.kron(codeunits[i][j],qr4[i,:,j,:]) for i in range(2) for j in range(2)),
              np.zeros((32,32),complex))
    state=np.einsum('aibi->ab',joint.reshape(16,2,16,2))
    remoteR=np.einsum('aijakl->ijkl',joint.reshape(4,4,2,4,4,2)).reshape(8,8)
    mx=max(residual(cov(reference,aa),pp),residual(cc@pp.conj()@cc,np.eye(8)-pp),
           residual(remoteR,np.kron(env,rR)),residual(fh@hh@fh,0))
    actualDouble=float(np.trace(state@aa[0].conj().T@aa[0]@aa[1].conj().T@aa[1]).real)
    assert abs(actualDouble)<1e-14 and np.linalg.det(rQ).real>.2
    t=.39;u=unitary(hh,t);uf=unitary(ham,t)
    A=u[2:4,:2];B=u[2:4,4:6]
    assert norm(B)>1e-3
    baseline=qc@pp@qc+fc
    for i in range(2):
        for j in range(2):
            x=np.zeros((2,2),complex);x[i,j]=1
            actual=cov(uf@codeunits[i][j]@uf.conj().T,aa)[2:4,2:4]
            expected=np.trace(x)*(u@baseline@u.conj().T)[2:4,2:4]+A@x@A.conj().T-B@x.T@B.conj().T
            mx=max(mx,residual(actual,expected))
    # Full source pullbacks, including arbitrary partner matrix effects.
    nn=np.array([[.4,.09j,.1,0],[-.09j,-.2,0,.3],[.1,0,.7,.05],[0,.3,.05,-.4]])
    jlist=[('normal',bdg(nn,np.zeros((4,4)))),('S_real',bdg(np.zeros((4,4)),dd)),
           ('S_imag',bdg(np.zeros((4,4)),1j*dd)),('energy',hh)]
    partner=np.array([[.2,.1+.07j],[.1-.07j,-.6]])
    sources=[]
    for name,js in jlist:
        jt=u.conj().T@js@u
        obs=uf.conj().T@quantize(js,aa)@uf+F(2,7).__float__()*np.eye(16)
        pulled=np.array([[np.trace(codeunits[j][i]@obs) for j in range(2)] for i in range(2)])
        constant=.5*np.trace(jt@baseline)+2/7
        pred=constant*np.eye(2)+jt[:2,:2]
        er=max(residual(pulled,pred),float(abs(np.trace(joint@np.kron(obs,partner))-
                                             np.trace(rhoQR@np.kron(pred,partner)))))
        mx=max(mx,er)
        sources.append(dict(name=name,max_residual=er,
            logical_traceless_norm=norm(pulled-np.trace(pulled)*np.eye(2)/2)))
    # Energy identity is relative to this model's own reference, with Nambu 1/2.
    rel=float(np.trace((state-reference)@ham).real)
    half=float((np.trace(fh@sp(hh,np.abs))/2).real)
    mx=max(mx,abs(rel-half));assert rel>0
    gauges=[]
    for k in (0,1):
        q=-3*k
        qa=np.diag([q]*4+[-q]*4)
        ward=residual(1j*(hh@qa-qa@hh)+2*q*bdg(np.zeros((4,4)),1j*dd),0)
        mx=max(mx,ward)
        def phases(s):return np.array([.17+1.7*s+.13*s*s]*2+[-.23+1.3*s-.09*s*s]*2)
        def rates(s):return np.array([1.7+.26*s]*2+[1.3-.18*s]*2)
        def G(s):return lift(np.diag(np.exp(1j*q*phases(s))))
        g0=G(0);gt=G(t); rate=np.diag(np.r_[q*rates(t),-q*rates(t)])
        hg=gt@hh@gt.conj().T-rate
        ug=gt@u@g0.conj().T;pg=gt@pp@gt.conj().T
        dotg=1j*rate@gt
        dotu=dotg@u@g0.conj().T-1j*gt@hh@u@g0.conj().T
        er=max(residual(cc@gt.conj()@cc,gt),residual(1j*dotu,hg@ug),
               residual(cc@pg.conj()@cc,np.eye(8)-pg))
        ag=ug[2:4,:2];bg=ug[2:4,4:6]
        pa=np.exp(1j*q*phases(0)[0]);pb=np.exp(1j*q*phases(t)[2])
        er=max(er,residual(ag,pb*A*pa.conjugate()),residual(bg,pb*B*pa))
        for x in (np.array([[1,0],[0,-1]]),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]])):
            er=max(er,residual(ag@x@ag.conj().T-bg@x.T@bg.conj().T,
                               A@x@A.conj().T-B@x.T@B.conj().T))
        # Full prepared state: local gauge unitaries also transport its anomalous environment.
        chargeF=sum((q*phases(t)[j]*aa[j].conj().T@aa[j] for j in range(4)),np.zeros((16,16),complex))
        vf=unitary(chargeF,-1)
        gauge_state=vf@state@vf.conj().T
        er=max(er,residual(cov(gauge_state,aa),gt@cov(state,aa)@gt.conj().T))
        for _,js in jlist:
            er=max(er,float(abs(np.trace(gauge_state@quantize(gt@js@gt.conj().T,aa))-
                                  np.trace(state@quantize(js,aa)))))
        inst=sp(hg,lambda x:(x<0).astype(float))
        wrongP=norm(inst-pg)
        omitted=norm(1j*dotu-(gt@hh@gt.conj().T)@ug)
        fixedS=norm(gt@hh@gt.conj().T-bdg((gt[:4,:4]@h@gt[:4,:4].conj().T),dd))
        if k:assert wrongP>.5 and omitted>1 and fixedS>.1
        mx=max(mx,er)
        gauges.append(dict(k=k,max_residual=er,instantaneous_negative_projection_mismatch=wrongP,
            omitted_temporal_connection_residual=omitted,untransported_S_pairing_defect=fixedS,
            hypercharge_S_Ward_residual=ward))
    # Chain rule: a background-dependent representation adds Euler/contact terms.
    generator=np.diag([.2,-.1,.3,-.4]);D=lift(np.exp(1j*.23*generator.diagonal())*np.eye(4))
    K=diag_blocks(generator,-generator)
    physical_source=jlist[0][1]
    transformedH=D.conj().T@hh@D
    transformedSource=D.conj().T@physical_source@D
    chain=transformedSource+1j*(transformedH@K-K@transformedH)
    eta=1e-5
    def changed(z):
        ddg=lift(np.diag(np.exp(1j*(.23+z)*generator.diagonal())))
        return ddg.conj().T@(hh+z*physical_source)@ddg
    derivative=(changed(eta)-changed(-eta))/(2*eta)
    chainerr=residual(chain,derivative)
    assert chainerr<2e-9 and norm(chain-transformedSource)>.05
    assert mx<4e-12
    return dict(finite_CAR_modes=4,Nambu_dimension=8,Fock_dimension=16,partner_dimension=2,
       diagnostic_coefficients_not_continuum_model=True,environment_anomalous_norm=env_anom,
       initial_remote_R_product_residual=residual(remoteR,np.kron(env,rR)),
       exact_double_occupation=actualDouble,wrong_quasifree_double=float(np.linalg.det(rQ).real),
       normal_leg_norm=norm(A),anomalous_leg_norm=norm(B),source_pullbacks=sources,
       relative_energy=rel,half_weight_energy=half,gauge_transports=gauges,
       chain_rule_finite_difference_residual=chainerr,
       omitted_dictionary_derivative_norm=norm(chain-transformedSource),max_residual=mx)

def original_field_Ward():
    # Undressed L/N-dagger components: different original hypercharges.
    # Two spin slots per chirality; pairing only in the right block.
    spin_eps=np.array([[0,1],[-1,0]],complex)
    m=np.eye(2)*1.17
    h=np.block([[np.zeros((2,2)),m],[m,np.zeros((2,2))]])
    pair=np.zeros((4,4),complex);pair[2:,2:]=.19*spin_eps
    hh=bdg(h,pair);data=[]
    for k in (0,1):
        qL=-.5-3*k;qR=-3*k
        q=np.diag([qL,qL,qR,qR]);Q=diag_blocks(q,-q)
        # H transforms with +1/2; this Hamiltonian's LR entry uses its conjugate.
        dhH=np.block([[np.zeros((2,2)),-.5j*m],[.5j*m,np.zeros((2,2))]])
        dirac=bdg(dhH,np.zeros((4,4)))
        majorana=bdg(np.zeros((4,4)),2j*qR*pair)
        total=1j*(hh@Q-Q@hh)+dirac+majorana
        err=residual(total,0)
        missingH=norm(total-dirac);missingS=norm(total-majorana)
        assert err<1e-13 and missingH>.5
        if k:assert missingS>1
        theta=.14
        G=lift(np.diag(np.exp(1j*theta*np.diag(q))))
        # Finite background phase includes both H and S.
        hp=np.block([[np.zeros((2,2)),np.exp(-.5j*theta)*m],
                     [np.exp(.5j*theta)*m,np.zeros((2,2))]])
        hp=bdg(hp,np.exp(2j*qR*theta)*pair)
        finiteerr=residual(hp,G@hh@G.conj().T)
        assert finiteerr<2e-14
        data.append(dict(k=k,left_hypercharge=qL,right_hypercharge=qR,
                         combined_Ward_residual=err,finite_covariance_residual=finiteerr,
                         omit_H_phase_defect=missingH,omit_S_phase_defect=missingS))
    return dict(original_undressed_chiral_source_normalization=True,diagnostics=data,
                generic_finite_coefficients_not_continuum_operator_bounds=True)

def charge_budget_control():
    # Dirac control only: vacuum fills two negative modes of bipartite h.
    # Resetting F to one occupation does NOT fix the whole charge sector.
    aa=car(4);h=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]])
    ham=quantize(bdg(h,np.zeros((4,4))),aa)
    _,v=np.linalg.eigh(ham);vac=np.outer(v[:,0],v[:,0].conj())
    env=np.einsum('aiaj->ij',vac.reshape(4,4,4,4))
    code=np.diag([0,.5,.5,0]);state=np.kron(code,env)
    n=sum((a.conj().T@a for a in aa),np.zeros((16,16),complex))
    diag=np.diag(n).real.astype(int);probs=[]
    for j in (1,2,3):probs.append(float(sum(state[r,r].real for r in range(16) if diag[r]==j)))
    assert residual(probs,[.25,.5,.25])<1e-13
    variance=float(np.trace(state@(n-2*np.eye(16))@(n-2*np.eye(16))).real)
    assert abs(variance-.5)<1e-13
    return dict(Dirac_only_negative_control=True,relative_charge_values_in_units_q=[-1,0,1],
      probabilities=probs,variance_in_units_q_squared=variance,
      no_claim_same_distribution_for_nonzero_pairing=True,
      logical_fixed_charge_not_full_reset_sharp_charge=True)

def bindings():
    paths={'admission':HERE.parent/'_shared/round1297_admission.md',
      'author_proof':HERE/'proof.md','author_check':HERE/'check.py','author_results':HERE/'results.json',
      'independent_check':Path(__file__), 'inherited_1296_proof':HERE.parent/'1296/proof.md',
      'background_dressing':HERE.parent/'_admission/background_dressing_screen_after1293.md'}
    return {key:{'path':str(p.relative_to(ROOT)).replace('\\','/'),'sha256':sha(p)} for key,p in paths.items()}
def compare(x,y):
    if isinstance(x,dict):
        assert isinstance(y,dict) and x.keys()==y.keys()
        for k in x:compare(x[k],y[k])
    elif isinstance(x,list):
        assert isinstance(y,list) and len(x)==len(y)
        for a,b in zip(x,y):compare(a,b)
    elif isinstance(x,float):assert math.isfinite(x) and math.isclose(x,y,rel_tol=2e-8,abs_tol=2e-12),(x,y)
    else:assert x==y,(x,y)
def main():
    parser=argparse.ArgumentParser();mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',action='store_true');mode.add_argument('--check',action='store_true')
    mode.add_argument('--diagnose',action='store_true');args=parser.parse_args()
    data=clean(dict(round=1297,independent_groups=1,status='PASS_INDEPENDENT_SCOPED',
       exact=exact_data(),Wilson=wilson_tests(),finite_dictionary=finite_dictionary(),
       original_field_Ward=original_field_Ward(),charge_control=charge_budget_control(),scope=dict(
       declared_classical_background_gauge_orbits_only=True,
       continuum_constants_inherited_analytically=True,
       sources_varied_before_fixing_background=True,
       own_reference_and_anomalous_environment_preserved=True,
       no_full_field_sharp_total_charge_claim=True,
       no_dynamic_Gauss_or_apparatus_certificate=True,
       no_uniform_all_k_source_variance_or_resource_bound=True)))
    if args.diagnose:
        print('DIAGNOSTIC_ONLY_NO_WRITE');print(json.dumps(data,ensure_ascii=False,indent=2));return
    data['source_bindings']=bindings()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:json.dump(data,f,ensure_ascii=False,indent=2,sort_keys=True);f.write('\n')
        print('PASS_INDEPENDENT_SCOPED EXCLUSIVE_WRITE')
    else:
        compare(data,json.loads(OUT.read_text(encoding='utf-8')));print('PASS_INDEPENDENT_SCOPED READ_ONLY_MATCH')
    print('Max finite algebra residual:',max(data['Wilson']['max_residual'],data['finite_dictionary']['max_residual']))
    print('Whole continuum window inherited; finite diagnostics are not a spatial discretization.')
if __name__=='__main__':main()
