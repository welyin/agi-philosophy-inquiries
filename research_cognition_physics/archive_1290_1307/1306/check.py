# 1306: fixed-background Majorana routing into finite left-hand spacetime ports.
# Default/--check read only; --write creates results.json exclusively.
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'
MU=F(1,2**30)
A=F(1,32)
L=F(1,4)
M=np.array([[0,1,0,0],[1,1,1,0],[0,1,2,1],[0,0,1,3]],dtype=np.int64)
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.]).astype(complex)
J=np.array([[0,1],[-1,0]],complex)

def det(a):
    a=[[F(int(v)) for v in row] for row in a]
    d=F(1)
    for i in range(len(a)):
        j=next((j for j in range(i,len(a)) if a[j][i]),None)
        if j is None: return F(0)
        if j!=i: a[i],a[j]=a[j],a[i]; d=-d
        p=a[i][i]; d*=p
        for j in range(i+1,len(a)):
            q=a[j][i]/p
            for k in range(i+1,len(a)): a[j][k]-=q*a[i][k]
    return d

def exact_structure():
    e=np.eye(4,dtype=np.int64)[:,0]
    C=np.column_stack([np.linalg.matrix_power(M,2*j)@e for j in range(4)])
    ps=[[F(1),F(0),F(0),F(0)],
        [F(-85,9),F(18),F(-13,3),F(2,9)],
        [F(76,9),F(-17),F(13,3),F(-2,9)],
        [F(-32,9),F(22,3),F(-2),F(1,9)]]
    for j,p in enumerate(ps):
        assert [sum(F(int(C[r,k]))*p[k] for k in range(4)) for r in range(4)]==[F(int(r==j)) for r in range(4)]
    assert det(M)==-5 and det(C)==27
    d=np.diag([1,0,0])
    ra=np.diag([1,2,3]); rb=np.array([[1,1,0],[1,2,1],[0,1,3]])
    ma=np.block([[np.zeros((3,3),dtype=int),d],[d,ra]])
    mb=np.block([[np.zeros((3,3),dtype=int),d],[d,rb]])
    pdark=np.diag([0,0,0,0,1,1]); pleft=np.diag([1,1,1,0,0,0])
    assert np.array_equal(ma@pdark,pdark@ma) and not np.any(pleft@pdark)
    idx=[0,3,4,5]
    assert np.array_equal(mb[np.ix_(idx,idx)],M)
    assert np.linalg.matrix_rank(mb)==4
    assert max(np.sum(np.abs(M),axis=1))<=5
    q=np.diag([-1,1,1,1])
    dr=M.copy(); dr[0,:]=0; dr[:,0]=0
    assert not np.any(q@M+M@q-2*dr)
    return {'Dirac_rank_both':1,'M0':M.tolist(),'det_M0':str(det(M)),
        'mass_squared_Krylov':C.tolist(),'det_Krylov':str(det(C)),
        'polynomials_ascending_z':[[str(v) for v in row] for row in ps],
        'all_polynomial_identities_exact':True,'A_dark_projector':pdark.tolist(),
        'A_dark_invariant_and_left_orthogonal':True,'B_full_six_Weyl_rank':4,
        'full_spectrum_has_two_massless_left_fields':True,
        'B_minus_L_joint_S_source_identity_exact':True,
        'sigma_M0_lower':'1/25','real_mass_background_only':True}

def exact_bounds():
    tmin,tmax=F(7,16),F(1,2)
    rmin,rmax=L-2*A,L+2*A
    smin,smax=tmin*tmin-rmax*rmax,tmax*tmax-rmin*rmin
    assert (smin,smax)==(F(3,32),F(55,256))
    # exp(m^2 S/4) < 1/(1-1/16)=16/15 for m<=1.
    eupper=F(16,15)
    off_coefficient=smax*eupper/F(64)
    normal_coefficient=(tmax+rmax)*eupper/F(32)
    assert off_coefficient<F(1,128)
    assert normal_coefficient<F(1,32)
    # mmax=5mu<=1; mmax^5<=mmax^4.
    norm_error_coef=F(5**4)*(F(1,128)+F(1,32))
    assert norm_error_coef<25<100
    ratio=F(12500000)*MU
    assert ratio<F(1,2)
    beta=A**3*MU**3/F(3000000)
    t0=F(15,32); eta=F(1,128); radius=A+eta; deadline=F(3,4)
    assert t0-eta>L+A+radius
    assert t0+eta+radius<deadline
    energy_square=9*(F(524288)+25*MU**2)
    assert energy_square<2173**2
    return {'mu':str(MU),'a':str(A),'L':str(L),'time_window':[str(tmin),str(tmax)],
        'timelike_s_squared_range':[str(smin),str(smax)],
        'exp_upper':'16/15','off_tail_coefficient_without_mass_W_pi':str(off_coefficient),
        'normal_tail_coefficient_without_mass_W_pi':str(normal_coefficient),
        'sufficient_error_coef':100,'sharper_error_coef':str(norm_error_coef),
        'error_to_leading_gap':str(ratio),'error_to_leading_gap_decimal':float(ratio),
        'beta0_exact':str(beta),'beta0_decimal':float(beta),
        'left_port_support':{'t0':str(t0),'eta':str(eta),'radius':str(radius),
            'deadline':str(deadline),'deadline_causal_slack':str(deadline-t0-eta-radius)},
        'field_relative_energy_strict_upper':2173,
        'port_access_and_cost_are_additional_inputs':True}

def nambu_symbol():
    p=np.array([.13,-.2,.31])
    dm=np.kron(M,J)
    def h(p):
        h2=-(p[0]*SX+p[1]*SY+p[2]*SZ)
        return np.block([[np.kron(np.eye(4),h2),dm],
                         [-dm,np.kron(np.eye(4),h2.conjugate())]])
    hp=h(p)
    expected=np.kron(np.eye(2),np.kron(p@p*np.eye(4)+M@M,np.eye(2)))
    swap=np.block([[np.zeros((8,8)),np.eye(8)],[np.eye(8),np.zeros((8,8))]])
    square=float(np.linalg.norm(hp@hp-expected))
    reality=float(np.linalg.norm(swap@h(-p).conjugate()@swap+hp))
    herm=float(np.linalg.norm(hp-hp.conjugate().T))
    src=[2,4,6,10,12,14]
    mass_compression=float(np.linalg.norm(h(np.zeros(3))[np.ix_(src,src)]))
    assert max(square,reality,herm,mass_compression)<1e-12
    return {'square_residual':square,'Gamma_momentum_reversal_residual':reality,
        'Hermiticity_residual':herm,'source_Nambu_mass_compression':mass_compression,
        'momentum_average_zero_uses_real_radial_packet':True}

def radial(n):
    z,w=np.polynomial.legendre.leggauss(n)
    r=(z+1)*float(A)/2; w=w*float(A)/2
    b=np.exp(-1/(1-(r/float(A))**2))
    norm=np.sqrt(4*np.pi*np.sum(w*r*r*b*b))
    b0=4*np.pi*np.sum(w*r*r*b)/norm
    u2=np.sum(w*r**4*b)/np.sum(w*r*r*b)
    u4=np.sum(w*r**6*b)/np.sum(w*r*r*b)
    return b0*b0,u2,u4

def flavor_basis():
    rows=[]
    for a in range(3):
        for b in range(a+1,3):
            x=np.zeros((3,3),complex); x[a,b]=x[b,a]=1/np.sqrt(2); rows.append(x)
            x=np.zeros((3,3),complex); x[a,b]=-1j/np.sqrt(2); x[b,a]=1j/np.sqrt(2); rows.append(x)
    rows.append(np.diag([1,-1,0])/np.sqrt(2))
    rows.append(np.diag([1,1,-2])/np.sqrt(6))
    return rows

def continuous_tail_diagnostic():
    w80,u280,u480=radial(80); W,u2,u4=radial(120)
    radial_error=max(abs(w80-W)/W,abs(u280-u2)/u2,abs(u480-u4)/u4)
    assert radial_error<1e-10
    assert W>np.pi*float(A)**3/12
    mu=float(MU); t=float(F(15,32)); ell=float(L)
    I1=W*(t*t-ell*ell-2*u2)
    I2=W*((t*t-ell*ell)**2-4*u2*(t*t-ell*ell)+(8/3)*ell*ell*u2+2*u4+(10/3)*u2*u2)
    dt1=2*t*W; dl1=-2*ell*W
    dt2=W*(4*t*(t*t-ell*ell)-8*t*u2)
    dl2=W*(-4*ell*(t*t-ell*ell)+(40/3)*ell*u2)
    m2=(mu**2)*(M@M); m4=m2@m2; m6=m4@m2
    k=-W*m2/(8*np.pi)+I1*m4/(64*np.pi)-I2*m6/(1536*np.pi)
    dt=dt1*m4/(64*np.pi)-dt2*m6/(1536*np.pi)
    dl=dl1*m4/(64*np.pi)-dl2*m6/(1536*np.pi)
    normal=np.kron(dt,np.eye(2))+np.kron(dl,SZ)
    off=-1j*np.kron((mu*M)@k,J)
    T=np.block([[normal,off],[-off,normal]])
    scale=W*mu**3/(8*np.pi)
    Tn=T/scale
    swap=np.block([[np.zeros((8,8)),np.eye(8)],[np.eye(8),np.zeros((8,8))]])
    assert np.linalg.norm(swap@Tn.conjugate()@swap-Tn)<1e-12
    sn=np.linalg.svd(Tn,compute_uv=False)
    assert min(sn)>1/(2*15625)
    indices=[2,4,6,10,12,14]
    Tc=Tn[:,indices]
    conv=np.block([[np.eye(8),np.eye(8)],[-1j*np.eye(8),1j*np.eye(8)]])
    columns=[]; covariance_errors=[]
    for x in flavor_basis():
        inp=np.block([[x,np.zeros((3,3))],[np.zeros((3,3)),-x.T]])
        d=Tc@inp@Tc.conjugate().T
        G=1j*conv@d@conv.conjugate().T
        covariance_errors.extend([np.max(np.abs(G.imag)),np.max(np.abs(G+G.T))])
        columns.append(G.real[np.triu_indices(16,1)]/2)
    tomography=np.column_stack(columns)
    tv=np.linalg.svd(tomography,compute_uv=False)
    assert np.linalg.matrix_rank(tomography)==8 and min(tv)>0
    assert max(covariance_errors)<1e-10
    # These are scaled field effects; C_port and actual devices are not computed here.
    return {'radial_orders':[80,120],'radial_relative_difference':radial_error,
        'W':W,'radial_u2':u2,'radial_u4':u4,'I1':I1,'I2':I2,
        'time':t,'series_terms_kept':3,'remaining_tail_certified_analytically':True,
        'scaled_Nambu_transfer_singular_values':sn.tolist(),
        'actual_transfer_scale':scale,'tail_truncation_diagnostic_min_singular':float(min(sn)*scale),
        'qutrit_traceless_tomography_rank':8,'effects':120,
        'scaled_tomography_singular_values':tv.tolist(),
        'covariance_reality_and_antisymmetry_residual':float(max(covariance_errors)),
        'port_C_and_device_statistics_not_numerically_claimed':True}

def annihilators(n):
    low=np.array([[0,1],[0,0]],complex); z=np.diag([1,-1]).astype(complex)
    ops=[]
    for j in range(n):
        op=np.array([[1]],complex)
        for k in range(n):
            op=np.kron(op,z if k<j else low if k==j else np.eye(2))
        ops.append(op)
    return ops

def covariance(rho,ops):
    z=ops+[a.conjugate().T for a in ops]
    return np.array([[np.trace(rho@b.conjugate().T@a) for b in z] for a in z])

def finite_car_partner_source():
    # A diagnostic finite CAR reference with cross and anomalous environment correlations.
    # It is not a replacement for either continuum Hadamard reference.
    ops=annihilators(5); dim=32
    number=[a.conjugate().T@a for a in ops]
    hf=sum((j+1)/7*number[j] for j in range(5))
    hop=ops[0].conjugate().T@ops[3]+.3j*ops[2].conjugate().T@ops[4]
    pair=.31*ops[3].conjugate().T@ops[4].conjugate().T+.19j*ops[1].conjugate().T@ops[3].conjugate().T
    hf=hf+hop+hop.conjugate().T+pair+pair.conjugate().T
    vals,vec=np.linalg.eigh(hf)
    weights=np.exp(-vals); weights/=sum(weights)
    ref=(vec*weights)@vec.conjugate().T
    P=covariance(ref,ops)
    env=np.einsum('aeaf->ef',ref.reshape(8,4,8,4))
    code=np.zeros((8,3),complex)
    code[4,0]=code[2,1]=code[1,2]=1
    psi=np.zeros(6,complex); psi[0]=.5; psi[3]=.5j; psi[4]=1/np.sqrt(2)
    qr=np.outer(psi,psi.conjugate())
    qr4=qr.reshape(3,2,3,2)
    rhoq=np.einsum('arbr->ab',qr4)
    w=np.kron(code,np.eye(2))
    source_r=w@qr@w.conjugate().T
    joint=np.einsum('arbs,ef->aer bfs'.replace(' ',''),source_r.reshape(8,2,8,2),env).reshape(64,64)
    prep=np.einsum('a r b r -> a b'.replace(' ',''),joint.reshape(32,2,32,2))
    Fhat=np.diag([1,1,1,0,0,1,1,1,0,0])
    Q=np.eye(10)-Fhat
    V=np.eye(10)[:,:3]; Vc=np.eye(10)[:,5:8]
    pred=Q@P@Q+V@rhoq@V.conjugate().T+Vc@(np.eye(3)-rhoq.T)@Vc.conjugate().T
    cov_error=float(np.linalg.norm(covariance(prep,ops)-pred))
    eigmin=float(min(np.linalg.eigvalsh(joint)))
    # General Hermitian, self-dual-odd source kernel with pairing and cross legs.
    h=np.diag(np.arange(1,6)/9).astype(complex)
    h[0,3]=.2j; h[3,0]=-.2j; h[1,2]=h[2,1]=.17
    delta=np.zeros((5,5),complex)
    for i,j,v in [(0,4,.11j),(2,3,.23),(3,4,.29j)]:
        delta[i,j]=v; delta[j,i]=-v
    kern=np.block([[h,delta],[-delta.conjugate(),-h.conjugate()]])
    z=ops+[a.conjugate().T for a in ops]
    obs=sum(.5*kern[i,j]*(z[i].conjugate().T@z[j]) for i in range(10) for j in range(10))
    const=.5*np.trace(kern@(Q@P@Q+Vc@Vc.conjugate().T))
    logical=V.conjugate().T@kern@V
    partner=np.array([[.4,.2-.13j],[.2+.13j,-.3]],complex)
    rhor=np.einsum('aras->rs',qr4)
    direct=np.trace(joint@np.kron(obs,partner))
    pull=const*np.trace(rhor@partner)+np.trace(qr@np.kron(logical,partner))
    source_error=float(abs(direct-pull))
    env_pair=float(abs(np.trace(prep@ops[3]@ops[4])))
    env_pair_ref=float(abs(np.trace(ref@ops[3]@ops[4])))
    double=float(np.trace(prep@number[0]@number[1]).real)
    wrong_wick=float(np.linalg.det(rhoq[:2,:2]).real)
    assert cov_error<1e-12 and source_error<1e-12
    assert eigmin>-1e-12 and abs(np.trace(joint)-1)<1e-12
    assert abs(env_pair-env_pair_ref)<1e-12 and env_pair>1e-3
    assert abs(double)<1e-12 and wrong_wick>1e-3
    return {'finite_modes':5,'logical_dimension':3,'partner_dimension':2,
        'covariance_CP_formula_residual':cov_error,'partner_source_pullback_residual':source_error,
        'joint_state_min_eigenvalue':eigmin,'environment_anomalous_pair_abs':env_pair,
        'actual_source_double_occupation':double,'incorrect_Gaussian_double_prediction':wrong_wick,
        'full_environment_restriction_preserved':True,
        'finite_reference_is_algebra_diagnostic_only':True}

def payload():
    sources=['proof.md','NEXT.md','../research_note_1306.md',
        '../_shared/round1306_admission.md','../_shared/comparison_scope.md',
        '../_admission/positive_cognitive_selection_screen.md','../1296/proof.md','../1303/proof.md',
        '../../archive_797_818/research_note_799.md',
        '../../archive_1009_1043/research_note_1011.md',
        '../../archive_1009_1043/research_note_1036.md',
        '../../archive_467_530/research_note_512.md',
        '../../../猜想/层级认知接口与尺度依赖空间维数猜想.md',
        '../../../猜想/从认知问题到社会机制再到物理检验的研究方法论.md']
    data={'status':'PASS_AUTHOR_LEFT_FUTURE_PORT_CONDITIONAL',
        'scope':'fixed-background quadratic S/Q plus O conditional on declared finite spacetime ports',
        'exact_structure':exact_structure(),'exact_full_window_bounds':exact_bounds(),
        'Nambu_symbol_diagnostic':nambu_symbol(),
        'continuous_tail_diagnostic':continuous_tail_diagnostic(),
        'finite_CAR_partner_sources':finite_car_partner_source(),
        'new_axioms':0,'autonomous_receiver_proven':False,
        'old_single_slice_port_automatically_grants_new_menu':False,
        'full_CC_or_M4_claimed':False,'thermal_geometry_window_inherited':False,
        'Dirac_rank_one_necessary_darkness_rejected_under_declared_port_condition':True,
        'source_bindings':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
        'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    assert data['source_bindings']['../_shared/round1306_admission.md']=='6349122eef0bb61b04e975e9024df4e70f8f62a998f4386f71586269c5621194'
    return json.loads(json.dumps(data,ensure_ascii=False))

def main():
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--write',action='store_true'); mode.add_argument('--check',action='store_true')
    args=parser.parse_args(); result=payload()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2); stream.write('\n')
        print('CREATED',result['status'])
    else:
        saved=json.loads(OUT.read_text(encoding='utf-8'))
        if result!=saved: raise SystemExit('MISMATCH: results or bound source changed')
        print('READ_ONLY_MATCH',result['status'])
    print('Krylov_det',result['exact_structure']['det_Krylov'])
    print('beta0',result['exact_full_window_bounds']['beta0_exact'])
    print('tomography_rank',result['continuous_tail_diagnostic']['qutrit_traceless_tomography_rank'])
    print('partner_source_residual',result['finite_CAR_partner_sources']['partner_source_pullback_residual'])

if __name__=='__main__': main()
