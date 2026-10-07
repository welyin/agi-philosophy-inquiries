"""972: native charge protection, existing exchange, and finite thermal contact.
The bath is one physical oscillator with explicit energy and preparation.
No full erasure, Markov thermalization or irreversible arrow is claimed.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json,math
from decimal import Decimal,localcontext
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/'native_thermal_record_results.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text('utf-8-sig'))
def norm(x):return float(np.linalg.norm(x,2))
def fro(x):return float(np.linalg.norm(x))
def entropy(r):
    p=np.linalg.eigvalsh((r+r.conj().T)/2)
    assert min(p)>-2e-13
    p=p[p>1e-15]
    return float(-p@np.log(p))
def hb(x):
    if x==0:return 0.
    return -x*math.log(x)-(1-x)*math.log1p(-x)
def continuity4(e):return e*math.log(3)+hb(e)
def gamma(n):
    u=np.finfo(float).eps/2
    return n*u/(1-n*u)
def construction_error(H,N,eta):
    with localcontext() as ctx:
        ctx.prec=65
        gd=Decimal(3)/1000;om=Decimal(1)/2;et=Decimal(0) if eta==0 else Decimal(1)/20
        hs={(0,1):-Decimal(1)/5,(1,0):-Decimal(1)/5,(1,1):Decimal(1),
            (2,2):Decimal(1),(0,3):-et*Decimal(3).sqrt()/4,
            (3,0):-et*Decimal(3).sqrt()/4,(3,3):-et/2}
        errors=np.zeros(4*N)
        for i in range(4*N):
            x,n=divmod(i,N)
            for j in range(4*N):
                y,k=divmod(j,N);target=Decimal(0)
                if n==k:target+=hs.get((x,y),Decimal(0))
                if x==y and n==k:
                    target+=om*n
                    if x in (1,2):target+=gd*gd/om
                if (x,y) in ((1,2),(2,1)) and abs(n-k)==1:
                    target+=gd*Decimal(max(n,k)).sqrt()
                errors[i]+=float(abs(Decimal.from_float(float(H[i,j]))-target))+1e-60
    return float(max(errors))*(1+1e-10)+1e-20
def certificate(H,ev,V,q,w,pr,N,g,t,definition_guard):
    # Sparse H V: each row has at most four nonzero terms. Bound the actual
    # dense BLAS accumulation conservatively by the count of nonzero products.
    count=int(max(np.count_nonzero(H,axis=1)))
    gram=fro(V.T@V-np.eye(len(H)))+gamma(len(H))*fro(abs(V).T@abs(V))+1e-13
    res=fro(H@V-V*ev)+gamma(count)*fro(abs(H)@abs(V))
    res+=gamma(3)*fro(abs(V)*abs(ev))+1e-13
    # Exact rational/radical entries were reconstructed using Decimal.
    res+=definition_guard*math.sqrt(1+gram)
    top=q@V[np.arange(4)*N+N-1,:]
    topnorm=np.linalg.norm(top,axis=0)
    boundary=[]
    for k in (0,1):
        starts=np.kron(w[:,k:k+1],np.eye(N))
        coeff=V.T@starts
        bn=topnorm@abs(coeff)
        boundary.append(float(math.sqrt(pr@(bn*bn)))*(1+1e-10)+1e-13)
    # The ONLY discarded link is g sqrt(N) Q |N><N-1|.
    leak=g*math.sqrt(N)*max(boundary)
    phase_guard=2e-10 # |t ev| <= 3.3e5 on the declared grid.
    state=gram+t*(res*math.sqrt(1+gram)+leak)+phase_guard
    # Thermal normalized truncation has exact trace distance r^N, r=1/2.
    trace=2*state+2.**(-N)
    return dict(nonzero_products_per_row=count,gram_bound=gram,
        eigen_residual_bound=res,defining_H_guard=definition_guard,
        boundary_amplitude_bounds=boundary,boundary_residual_bound=leak,
        phase_guard=phase_guard,state_vector_error_bound=state,
        material_trace_distance_bound=trace,
        material_record_information_error_bound=2*continuity4(trace)+1e-10)
def run():
    spec=importlib.util.spec_from_file_location('old968',STAGE/'968/internal_relay.py')
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    specp=importlib.util.spec_from_file_location('phases965',STAGE/'965/material_field_window.py')
    oldphase=importlib.util.module_from_spec(specp);specp.loader.exec_module(oldphase)
    material,h,q,f,w,intertwining=old.material()
    ps=np.diag([1.,1.,1.,0.]);N=64;om=.5;g=.003;beta=2*math.log(2)
    a=np.diag(np.sqrt(np.arange(1,N)),1);num=np.diag(np.arange(N))
    pr=2.**(-np.arange(N));pr/=sum(pr);tau=np.diag(pr)
    S0=math.log(2);SB0=entropy(tau);logZ=math.log(sum(2.**(-np.arange(N))))
    init=[np.kron(np.outer(w[:,k],w[:,k]),tau) for k in (0,1)]
    rho0=(init[0]+init[1])/2;result=[];algebra=[]
    for eta in (0.,.05):
        hm=h+eta*f;HM=np.kron(hm,np.eye(N));HB=np.kron(np.eye(4),om*num)
        HI=g*np.kron(q,a+a.T)+g*g/om*np.kron(q@q,np.eye(N))
        H=HM+HB+HI;ev,V=np.linalg.eigh(H)
        definition_guard=construction_error(H,N,eta)
        assert definition_guard<1e-13
        stack=np.vstack([np.kron(np.eye(4),hm)-np.kron(hm.T,np.eye(4)),
                         np.kron(np.eye(4),q)-np.kron(q.T,np.eye(4))])
        sv=np.linalg.svd(stack,compute_uv=False)
        nullity=int(sum(sv<1e-12));assert nullity==(2 if eta==0 else 1)
        P=np.kron(ps,np.eye(N))
        comm=norm(H@P-P@H)
        assert abs(comm-eta*math.sqrt(3)/4)<1e-13
        # In basis S,D+,T, D+ is cyclic for h_eta if eta != 0.
        block=hm[np.ix_([0,1,3],[0,1,3])];v=np.array([0.,1.,0.])
        krylov=float(np.linalg.det(np.stack([v,block@v,block@block@v],axis=1)))
        assert abs(krylov-eta*.1**2*math.sqrt(3))<1e-16
        algebra.append(dict(eta=eta,commutant_dimension=nullity,sector_commutator=comm,
            charge_commutator=norm(q@ps-ps@q),cyclic_determinant=krylov,
            singular_values=sv.tolist()))
        for t in (100.,1000.,10000.):
            U=(V*oldphase.phases(ev,t))@V.T
            full=[U@r@U.conj().T for r in init]
            # Renormalize at O(machine epsilon); the certificate uses 2*state.
            full=[r/np.trace(r).real for r in full]
            sm=[np.einsum('anbn->ab',r.reshape(4,N,4,N)) for r in full]
            rb=[np.einsum('aman->mn',r.reshape(4,N,4,N)) for r in full]
            rho=(full[0]+full[1])/2;rm=(sm[0]+sm[1])/2;rf=(rb[0]+rb[1])/2
            ent=entropy(rm);eb=float(np.trace(rf@(om*num)).real)
            heat=eb-float(pr@(om*np.arange(N)))
            I_record=ent-(entropy(sm[0])+entropy(sm[1]))/2
            I_field=entropy(rf)-(entropy(rb[0])+entropy(rb[1]))/2
            I_joint=ent+entropy(rf)-S0-SB0
            rel=-entropy(rf)+beta*eb+logZ
            sigma=ent-S0+beta*heat
            assert abs(sigma-I_joint-rel)<1e-13
            assert min(I_field,I_joint,rel,sigma)>-2e-12
            delta=[float(np.trace((rho-rho0)@A).real) for A in (HM,HB,HI)]
            assert abs(sum(delta))<1e-11
            currents=[float(np.trace(rho@(1j*(H@A-A@H))).real) for A in (HM,HB,HI)]
            assert abs(sum(currents))<1e-12
            cert=certificate(H,ev,V,q,w,pr,N,g,t,definition_guard)
            eps=cert['material_trace_distance_bound']
            # For infinite and projected evolution, energy positivity bounds
            # sqrt(<N>) <= sqrt((E-hmin)/omega)+g/omega. Use common upper 2;
            # it is verified from the declared initial energy below.
            e0=float(np.trace(rho0@H).real)
            nroot=math.sqrt((e0-min(np.linalg.eigvalsh(hm)))/om)+g/om
            assert nroot<2
            # |delta <a>| <= eps (sqrt(<N>_1)+sqrt(<N>_2+1))
            # together with bounded hm and Q^2 and conserved total energy.
            heat_error=(2*norm(hm)+2*g*(2+math.sqrt(5))+2*g*g/om)*eps
            heat_error+=om*(N+1)*2.**(-N)+1e-12
            entropy_error=continuity4(eps)+beta*heat_error+1e-10
            contrast=sum(abs(np.linalg.eigvalsh(sm[0]-sm[1])))/2
            row=dict(eta=eta,time=t,record_information=I_record,
                record_information_loss=S0-I_record,field_label_information=I_field,
                global_label_information=S0,material_trace_distance=float(contrast),
                material_entropy=ent,material_entropy_change=ent-S0,
                field_heat=heat,material_field_mutual_information=I_joint,
                field_relative_entropy_to_initial_Gibbs=rel,entropy_production=sigma,
                energy_changes=delta,energy_currents=currents,
                landauer_identity_error=abs(sigma-I_joint-rel),
                oscillator_number_root_bound=nroot,heat_error_bound=heat_error,
                entropy_production_error_bound=entropy_error,certificate=cert)
            if eta==0:
                assert abs(I_record-S0)<1e-12 and abs(contrast-1)<1e-12
            result.append(row)
    witness=result[-1]
    assert witness['record_information_loss']>witness['certificate']['material_record_information_error_bound']
    # Internal unitary alone cannot reduce label information (not merely Ps population).
    vals,B=np.linalg.eigh(h+.05*f);U=(B*np.exp(-1j*10000*vals))@B.T
    pure=[np.outer(U@w[:,k],(U@w[:,k]).conj()) for k in (0,1)]
    control=entropy((pure[0]+pure[1])/2)-sum(entropy(r) for r in pure)/2
    assert abs(control-S0)<1e-12
    files=[STAGE/'968/internal_relay.py',STAGE/'965/material_field_window.py',
        STAGE/'research_note_971.md',HERE/'drafts/thermal_record_decision.md']
    return dict(round=972,all_scientific_checks_passed=True,
        inputs=dict(U=1.,v=.1,omega=om,g=g,beta=beta,occupation_ratio=.5,
            oscillator_computational_dimension=N,thermal_tail=2.**(-N),
            old_exact_material_intertwining=intertwining),
        algebra=algebra,rows=result,isolated_material_control_information=control,
        witness=dict(eta=.05,time=10000.,information_loss_lower=witness['record_information_loss']-
            witness['certificate']['material_record_information_error_bound']),
        scope=dict(same_native_material_and_old_internal_exchange=True,
            finite_time_information_transfer=True,arbitrary_electric_bath_reset_obstruction=True,
            finite_thermal_reservoir_energy_accounted=True,
            complete_reset_proved=False,thermalization_proved=False,
            monotone_irreversible_arrow_proved=False,cosmology_common_window_proved=False,
            full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=['https://arxiv.org/abs/1306.4352'])
def compare(a,b,path=''):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'/'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'/{i}')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=2e-11),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=run()
    if args.write:
        with OUT.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(OUT))
    print(json.dumps(dict(round=972,witness=out['witness'],algebra=out['algebra'],
        last_row=out['rows'][-1],scope=out['scope']),ensure_ascii=False,indent=2))
