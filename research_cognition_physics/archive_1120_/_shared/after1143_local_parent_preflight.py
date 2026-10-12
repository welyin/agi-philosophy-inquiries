"""Unnumbered local-parent preflight, default read-only.
Matrix tests diagnose exact proofs; no full sea, autonomous apparatus or continuum simulation.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCES={
 'after1143_direction_screen.md':'43fa1cdf9799a53ef2c4e98dd1295067afedbcab838cbe2e12595f491bb369c0',
 '../1127/proof.md':'d1c8827cc12563b996c951da6039e041e08aa711de18e84a56bb240660660c27',
 '../../archive_742_763/research_note_742.md':'01382af9384d07bb0847bfaf6690d4cd051aafce9436c25b9e96e87263403329',
 '../../archive_742_763/research_note_743.md':'9d31e9ec2ad8f9d02ec74cd99b18ba27584ee0aaba8cbfdd9bf0d7a084bac83e',
}
def norm(A):return float(np.linalg.norm(A,2))
def expH(H,t=1):
    ev,V=np.linalg.eigh(H)
    return (V*np.exp(-1j*t*ev))@V.conj().T

def inverse_check():
    rng=np.random.default_rng(114304)
    residual=[];hermitian=[];bounds=[];lower=[]
    for rank in (0,1,2,4,5):
        n,k=7,5
        U=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))[0]
        W=np.linalg.qr(rng.normal(size=(k,k))+1j*rng.normal(size=(k,k)))[0]
        D=np.zeros((n,k));D[np.arange(rank),np.arange(rank)]=np.linspace(.4,1.7,rank)
        F=U@D@W.conj().T
        seed=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));seed=(seed+seed.conj().T)/2
        R=seed@F;Fp=np.linalg.pinv(F,rcond=1e-12)
        P=F@Fp;A=R@Fp;Q=np.eye(n)-P
        V=A+A.conj().T-P@A@P
        residual.append(norm(V@F-R));hermitian.append(norm(V-V.conj().T))
        bounds.append(norm(V)-norm(A)-norm(Q@A));lower.append(norm(A)-norm(seed))
        assert norm(R@(np.eye(k)-Fp@F))<1e-12
    assert max(residual+hermitian)<1e-12
    assert max(bounds+lower)<1e-12
    F=np.array([[1.,1.],[0.,0.]])
    R=np.array([[0.,0.],[0.,1.]])
    assert norm(F.T@R-R.T@F)==0
    kernel=np.array([1.,-1.])/math.sqrt(2)
    assert np.linalg.norm(F@kernel)==0 and np.linalg.norm(R@kernel)>.7
    return dict(ranks=[0,1,2,4,5],max_equation_residual=max(residual),
                max_Hermitian_residual=max(hermitian),max_upper_bound_excess=max(bounds),
                max_lower_bound_excess=max(lower),Hermitian_without_kernel_defect=float(np.linalg.norm(R@kernel)))

def boundary_check():
    ts=np.array([.1,.01,.001])
    forced=ts/(ts*ts)
    assert np.allclose(forced,[10,100,1000])
    plus=[];minus=[]
    for n in range(3):
        tp=1/(math.pi/2+2*math.pi*n);tm=1/(3*math.pi/2+2*math.pi*n)
        plus.append(math.sin(1/tp));minus.append(math.sin(1/tm))
    assert max(abs(v-1) for v in plus)<1e-12
    assert max(abs(v+1) for v in minus)<1e-12
    return dict(unbounded_forced_V=forced.tolist(),
                flat_smooth_data_plus_sequence=plus,flat_smooth_data_minus_sequence=minus,
                canonical_pseudoinverse_jump='F=t,R=t: A=1 off zero,A(0)=0, but V=1 is smooth',
                scope='analytic examples proved in Markdown; finite samples diagnose them')

def positive_check():
    I2=np.eye(2);sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
    PZ=np.diag([0.,1.]);PX=(I2+sx)/2
    records=np.zeros((12,12),complex);continuation=np.zeros((12,12),complex)
    for j,P in ((1,PZ),(2,PX)):
        records[4*j:4*j+4,4*j:4*j+4]=np.kron(P,sx)
        continuation[4*j:4*j+4,4*j:4*j+4]=np.kron(sy,I2)
    assert abs(norm(records)-1)<1e-12 and abs(norm(continuation)-1)<1e-12
    comm=norm(records@continuation-continuation@records);assert comm>.9
    Ur=expH(records,math.pi/2);Uc=expH(continuation,math.pi/4);U=Uc@Ur
    errors=[];isometries=[]
    for j,P in ((0,np.zeros((2,2))),(1,PZ),(2,PX)):
        Ein=np.zeros((12,2),complex);Ein[4*j,0]=1;Ein[4*j+2,1]=1
        W=U@Ein;expected=np.zeros_like(W)
        Qrot=I2 if j==0 else expH(sy,math.pi/4)
        expected[[4*j,4*j+2],:]=Qrot@(I2-P)
        expected[[4*j+1,4*j+3],:]=-1j*Qrot@P
        errors.append(norm(W-expected));isometries.append(norm(W.conj().T@W-I2))
    assert max(errors+isometries)<1e-12
    # Independent free-Dirac PDE diagnostic using a finite plane-wave sum.
    # Compact support of the actual Cauchy packet is proved analytically, not simulated here.
    zero=np.zeros((2,2));alphas=[np.block([[zero,s],[s,zero]]) for s in (sx,sy,sz)]
    beta=np.diag([1.,1.,-1.,-1.]);xi=np.array([1,0,1,0],complex)/math.sqrt(2)
    t=.37;x=np.array([.12,-.31,.17]);u=np.zeros(4,complex);du=np.zeros(4,complex);grad=np.zeros((3,4),complex)
    for amp,k in [(1.,np.array([.3,.2,-.4])),(.4,np.array([-.6,.7,.2])),(-.2,np.array([.2,-.5,.8]))]:
        hp=beta+sum(k[j]*alphas[j] for j in range(3))
        term=amp*np.exp(1j*k@x)*(expH(hp,t)@xi)
        u+=term;du+=-1j*hp@term
        for j in range(3):grad[j]+=1j*k[j]*term
    hu=-1j*sum(alphas[j]@grad[j] for j in range(3))+beta@u
    free_error=float(np.linalg.norm(1j*du-hu))
    K=.8*continuation;Ut=expH(continuation,.2)@Ur
    F=np.kron(u[:,None],Ut)
    dF=np.kron(du[:,None],Ut)-1j*np.kron(u[:,None],K@Ut)
    hF=np.kron(hu[:,None],Ut);R=1j*dF-hF;V=np.kron(np.eye(4),K)
    pde=norm(V@F-R);compat=norm(F.conj().T@R-R.conj().T@F)
    assert max(free_error,pde,compat)<1e-11
    return dict(internal_flavors=12,record_continue_commutator_norm=comm,
                isometry_max_error=max(isometries),instrument_continue_max_error=max(errors),
                free_Dirac_diagnostic_error=free_error,local_inverse_PDE_error=pde,
                Hermitian_compatibility_error=compat,
                task='same-carrier internal record, not independent receiver or pointer',
                scope='mode-level exact construction; full normal Fock/source/clock bridge not claimed')

def slater_check():
    n=4;dim=16;cs=[]
    for j in range(n):
        c=np.zeros((dim,dim),complex)
        for bits in range(dim):
            if bits&(1<<j):c[bits^(1<<j),bits]=(-1)**((bits&((1<<j)-1)).bit_count())
        cs.append(c)
    vac=np.zeros(dim,complex);vac[0]=1
    a0p0=cs[0].conj().T@cs[2].conj().T@vac
    a1p1=cs[1].conj().T@cs[3].conj().T@vac
    target=(a0p0-1j*a1p1)/math.sqrt(2)
    initial=(cs[0].conj().T+cs[1].conj().T)/math.sqrt(2)@cs[2].conj().T@vac
    n1=cs[1].conj().T@cs[1]
    pointer_x=cs[3].conj().T@cs[2]+cs[2].conj().T@cs[3]
    quartic=n1@pointer_x
    recording_error=float(np.linalg.norm(expH(quartic,math.pi/2)@initial-target))
    assert recording_error<1e-12
    A=np.zeros((4,4),complex);A[0,2]=1/math.sqrt(2);A[1,3]=-1j/math.sqrt(2);A=A-A.T
    plucker=A[0,1]*A[2,3]-A[0,2]*A[1,3]+A[0,3]*A[1,2]
    operator_bound=norm(A)**2
    assert abs(abs(plucker)-.5)<1e-12 and abs(operator_bound-.5)<1e-12
    saturating=float(abs(np.vdot(target,a0p0))**2)
    rng=np.random.default_rng(1143742);overlaps=[]
    for _ in range(256):
        Q=np.linalg.qr(rng.normal(size=(4,2))+1j*rng.normal(size=(4,2)))[0]
        x,y=Q[:,0],Q[:,1]
        dx=sum(x[j]*cs[j].conj().T for j in range(4))
        dy=sum(y[j]*cs[j].conj().T for j in range(4))
        psi=dx@dy@vac;overlaps.append(float(abs(np.vdot(target,psi))**2))
    assert max(overlaps)<=.5+1e-12
    covariance=np.array([[np.vdot(target,cs[j].conj().T@cs[i]@target) for j in range(4)] for i in range(4)])
    assert norm(covariance-.5*np.eye(4))<1e-12
    return dict(quartic_recording_error=recording_error,Plucker_abs=abs(plucker),
                squared_max_Slater_overlap=operator_bound,saturating_overlap=saturating,
                sampled_Slater_max=max(overlaps),target_covariance_eigenvalues=np.linalg.eigvalsh(covariance).tolist(),
                reused_result='742 Gaussian realization restriction;743 quantum scalar removes premise',
                exclusions='no assertion against nonGaussian preparation, measurement/feedback, or operator-valued sources',
                scope='specific independent two-dual-rail record; not a new general Gaussian theorem')

def run():
    for name,want in SOURCES.items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==want,name
    return dict(status='PASS_PREFLIGHT_ONLY',numbered_round=None,source_hashes=SOURCES,
                inverse=inverse_check(),boundary=boundary_check(),positive=positive_check(),
                independent_record=slater_check(),scientific_calibration_increment=0)
def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for j,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(j))
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write-results',action='store_true');args=ap.parse_args()
    result=run();target=Path(__file__).with_suffix('.json')
    if args.write_results:
        if target.exists():raise FileExistsError(target)
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    elif target.exists():compare(result,json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps(result,ensure_ascii=False,indent=2))
