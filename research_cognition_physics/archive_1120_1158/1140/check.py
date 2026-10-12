"""1140: original-mode finite collection bridge, three verification groups.
Default replay is read-only. Numerical finite CAR diagnostics do not simulate the continuous sea.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FROZEN={'1121/proof.md': '300c1a401edc2789059b7773013a12a2ad8eb225ad3d8a14dce0235a8141deb4', '1126/proof.md': '06644f8cc78f2075d62b67dee961a4064208c5226047405652201c337b5f3a76', '1128/proof.md': 'fc9963b8f77d0a60e2233acbcc68e02e583af345a5653d1d4431ea619be0fcd2', '1133/proof.md': '476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d', '1137/proof.md': 'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb', '1139/acceptance.json': '9e32524e342b65798f0edc81a4465bb82a92341531adbfbcba7509bc8f7d5c81', '_shared/after1139_collection_admission.md': '8fb7011e1bc4d08a78cdc7a891eeabefa302bedbac041d944fe4d0fb11268b86'}

def analytic():
    n=257;N=n**3;a=F(1,32);b=F(1,10**21);L=n*a
    volume_coefficient=6*(n+1)*L*L
    linfty_squared=F(1801**2,24) # pi>3
    assert volume_coefficient<101000 and linfty_squared<136000
    loss=F(136000*101000)*b
    delta=F(1,10**5)
    assert 4*loss<delta*delta
    unitary_error=F(44,7)*delta; isometry_error=6*delta
    assert unitary_error<F(1,1000)
    slot=3*a
    assert 3*a*a/4<F(1,20)**2 # center cube contained in readout ball .05
    assert slot*slot>6*a*a # declared slot > two face-adjacent cubes diameter
    rows={}
    # Actual schedule construction on a smaller three-dimensional grid, no spectral physics truncation.
    for m in (3,5,7):
        k=(m-1)//2; vertices=[(x,y,z) for x in range(-k,k+1) for y in range(-k,k+1) for z in range(-k,k+1)]
        ordered=sorted((v for v in vertices if v!=(0,0,0)),key=lambda v:(-sum(abs(c) for c in v),v))
        processed=set(); maxdist=0
        for v in ordered:
            parent=list(v)
            j=next(j for j in range(3) if v[j]!=0)
            parent[j]-=1 if v[j]>0 else -1; parent=tuple(parent)
            assert sum(abs(v[j]-parent[j]) for j in range(3))==1
            assert parent not in processed
            processed.add(v);maxdist=max(maxdist,sum(abs(c) for c in v))
        assert len(processed)==m**3-1 and maxdist==3*k
        rows[str(m)]={'vertices':m**3,'collection_events':len(ordered),'tree_depth':maxdist}
    return dict(grid_side=n,cell_modes=N,rails=2,collection_events=N-1,return_events=N-1,
        central_events=2,total_events=2*N,tree_depth=384,
        fractions={'cell_side':str(a),'corridor_width':str(b),'volume_coefficient':str(volume_coefficient),
        'linfty_squared_upper':str(linfty_squared),'lost_L2_squared_upper':str(loss),
        'mode_norm_error_upper':str(delta),'full_CP_unitary_error_upper':str(unitary_error),
        'blank_pointer_isometry_error_upper':str(isometry_error),'slot':str(slot),
        'record_time_after_new_start':str(N*slot),'complete_time_after_new_start':str(2*N*slot)},
        schedule_checks=rows,old_deadline_preserved=False,actual_clock_schedule_certified=False,
        continuous_microcausality_claimed=False)

def annihilate(n,i):
    a=np.zeros((2**n,2**n),complex)
    for s in range(2**n):
        if (s>>i)&1:a[s^(1<<i),s]=(-1)**((s&((1<<i)-1)).bit_count())
    return a

def exp_h(h,t):
    w,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*w))@v.conj().T

def qmat(a,c):
    out=np.zeros_like(c[0])
    for i,j in zip(*np.nonzero(abs(a)>1e-15)):
        out+=a[i,j]*c[i].conj().T@c[j]
    return out

def opnorm(a):return float(np.linalg.svd(a,compute_uv=False)[0])

def car():
    # Three cell modes plus an unresolved direction in each of two rails: complete 8-mode Fock.
    c=[annihilate(8,j) for j in range(8)]
    w=np.array([1.,2.,3.])/math.sqrt(14); z=np.r_[w,0.]
    C=np.eye(4);current=z.copy();G=np.eye(256,dtype=complex)
    for parent,child in ((1,2),(0,1)):
        a,b=current[parent],current[child];s=math.hypot(a,b)
        cosine,sine=a/s,b/s;theta=math.atan2(b,a)
        one=np.eye(4);one[parent,parent]=cosine;one[child,child]=cosine
        one[parent,child]=sine;one[child,parent]=-sine
        A=np.zeros((8,8),complex)
        for rail in (0,1):
            p=4*rail+parent;k=4*rail+child
            A[p,k]=1j;A[k,p]=-1j
        G=exp_h(qmat(A,c),theta)@G
        C=one@C;current=one@current
    assert np.linalg.norm(current-np.eye(4)[0])<1e-14
    X=np.array([[0,1],[1,0]],complex);I=np.eye(256);IP=np.eye(512)
    def V(f,setting):
        cs=[sum(np.conj(f[j])*c[4*r+j] for j in range(4)) for r in (0,1)]
        meas=cs[0] if setting=='Z' else (cs[0]+cs[1])/math.sqrt(2)
        N=meas.conj().T@meas
        Y=1j*(cs[1].conj().T@cs[0]-cs[0].conj().T@cs[1])
        record=IP-np.kron(N,np.eye(2))-1j*np.kron(N,X)
        return np.kron(exp_h(Y,math.pi/4),np.eye(2))@record
    GP=np.kron(G,np.eye(2)); center=np.eye(4)[0];out=[]
    # Full operator norms include every occupation sector, pointer phase and arbitrary tensor references.
    for setting in ('Z','X'):
        approx=V(z,setting);compiled=GP.conj().T@V(center,setting)@GP
        residual=opnorm(compiled-approx);assert residual<1e-12
        for tail in (1e-5,.02):
            f=math.sqrt(1-tail*tail)*z+tail*np.eye(4)[3]
            delta=np.linalg.norm(f-z);error=opnorm(compiled-V(f,setting))
            assert error<=2*math.pi*delta+1e-12
            # P starts at |0>: restrict input columns; no postselection.
            isometry_error=opnorm((compiled-V(f,setting))[:,::2])
            assert isometry_error<6*delta+1e-12
            out.append(dict(setting=setting,tail=tail,mode_error=float(delta),
                complete_CP_unitary_error=error,blank_P_isometry_error=isometry_error,compiler_residual=residual))
    survival=float(np.sum(w**4));assert abs(survival-.5)<1e-14
    return dict(full_Fock_dimension=256,pointer_dimension=2,
        exact_compilation_and_error=out,position_dephasing_g_occupation=survival,
        scope='Full finite CAR algebra diagnostic; no truncation of the original continuous sea in the proof.')

def background():
    c=[annihilate(4,j) for j in range(4)]
    h=np.zeros((4,4));h[0,1]=h[1,0]=h[2,3]=h[3,2]=1
    g=np.array([1,0,1,0])/math.sqrt(2);v=np.array([1,0,-1,0])/math.sqrt(2)
    bm=np.array([0,1,0,-1])/math.sqrt(2);P=np.diag([1,0,1,0]);Pg=np.outer(g,g)
    Qg=h@Pg+Pg@h-Pg@h@Pg;Qw=h@P+P@h-P@h@P
    H=qmat(h,c)+2*np.eye(16) # original negative-energy normal-order subtraction
    Hg=H-qmat(Qg,c);Hw=H-qmat(Qw,c)
    vac=np.eye(16)[0];cv=sum(v[j]*c[j].conj().T for j in range(4));psi=cv@vac
    Nbm=qmat(np.outer(bm,bm),c);t=math.pi/2
    old=exp_h(Hg,t)@psi;new=exp_h(Hw,t)@psi
    po=float(np.real(np.vdot(old,Nbm@old)));pn=float(np.real(np.vdot(new,Nbm@new)))
    assert po>1-1e-12 and abs(pn)<1e-12 and abs(np.vdot(old,new))<1e-12
    assert opnorm(Qw-Qg)>=1-1e-12
    # One cell-to-cell finite even gate: exact physical-field energy defect, with sea term retained.
    A=np.zeros((4,4),complex);A[0,2]=1j;A[2,0]=-1j
    theta=.7;u=exp_h(A,theta);U=exp_h(qmat(A,c),theta)
    B=u.conj().T@h@u-h
    defect=U.conj().T@H@U-H
    identity_residual=opnorm(defect-qmat(B,c));assert identity_residual<1e-12
    ew,ev=np.linalg.eigh(h);neg=ev[:,ew<0]
    omega=vac.astype(complex)
    for j in range(2):omega=sum(neg[k,j]*c[k].conj().T for k in range(4))@omega
    sea=np.real(np.trace(neg@neg.conj().T@B))
    work=float(np.real(np.vdot(omega,defect@omega)))
    assert work>.1 and abs(work-sea)<1e-12
    return dict(mode_error=0.,time=t,old_bminus_probability=po,new_bminus_probability=pn,
        background_state_overlap=float(abs(np.vdot(old,new))),
        single_particle_holding_difference_norm=opnorm(Qw-Qg),full_CAR_holding_difference_norm=opnorm(qmat(Qw-Qg,c)),
        gate_energy_defect_identity_residual=identity_residual,
        normal_sea_gate_work=work,normal_order_scalar=float(sea),
        scope='Finite normal CAR counterexample to a universal background-error inference; not a bound on the original fixed input at the new long time.',
        autonomous_controller_work_proved=False,total_four_momentum_conservation_proved=False)

def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+'/'+str(i))
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-7,abs_tol=2e-12),(path,a,b)
    else:assert a==b,(path,a,b)

def main():
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    for name,digest in FROZEN.items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    result=dict(status='PASS',scope='ORIGINAL_RECEIVER_FINITE_COHERENT_COLLECTION_BOUNDARY_BRIDGE_AND_BACKGROUND_COST',
        frozen_sources=FROZEN,numerical_check_groups=3,analytic=analytic(),CAR=car(),background=background(),
        new_axioms=0,science_increment=0,old_complete_dynamics_localization_claimed=False,
        actual_Y_clock_and_source_implementation_claimed=False,stage_completed=False)
    target=HERE/'results.json'
    if args.write_results:
        assert not target.exists(),'Results frozen after initial registration'
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:compare(result,json.loads(target.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',groups=3,frozen_sources=len(FROZEN),
        old_complete_dynamics_localized=False,results_registered=args.write_results)))
if __name__=='__main__':main()
