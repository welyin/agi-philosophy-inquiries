"""1139: three reproducible verification groups; default read-only.
Continuous domain and finite-parent scope are proved in proof.md.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math, types
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
FROZEN={'1129/proof.md': 'c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa', '1133/proof.md': '476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d', '1137/proof.md': 'b665fe6ac4a8f326b5289b93321206666c1c24d5684e7277d5b683f0ebad1beb', '1137/results.json': '1fb7f63266bc240d80ab880af707414d4e41a92c2075e14b69459e53b819ceb4', '1137/acceptance.json': '8101ad968bca9dcc40ccb19c35280a69d64c23cf84603e11788cea7b1b69adad', '1138/acceptance.json': 'e714abf633301cd7b0bbd288781c0c8b16ae50bcc2b0c92be9cba004a801d342', '_shared/local1138_dynamics_preflight.md': '56a559135a1288336d38907cfdf14e7fa5d23b4f3d1ec7268f83b97936d7f1b9', '_shared/local1138_dynamics_preflight.json': '3482a6c485289c414e4bae28243ac052902ae067c4eb987609c14aaf3f1a13db', '_shared/local1139_tail_certificate.py': 'fbc400a14c92788f00efd4579f4197e96ff579bacc2f42e1ffb92323ee478dd3', '_shared/local1139_tail_certificate.md': 'a502e98a6bcf634e6976c394a57bcf226b402447cfad1b61b1fe5d5507427a90', '_shared/local1139_tail_certificate.json': '27a7ae96911148f5073c633b0f0b2285d92b9c13acab87408bc3b48e756664cd', '_shared/after1138_refinement_admission.md': '07d182d478a1ac05ef104b9c5f315154319dad97c2a10dde6c3fcd5e0f133b7a'}

def load_tail():
    p=ROOT/'_shared/local1139_tail_certificate.py'
    m=types.ModuleType('tail_certificate');m.__file__=str(p)
    exec(compile(p.read_bytes(),str(p),'exec'),m.__dict__)
    data=m.run()
    compare(data,json.loads(p.with_suffix('.json').read_text(encoding='utf-8')))
    return data

def exact():
    T,s=F(4),F(11,5);tau=T-s
    M,A,sig=F(10**80),F(10**34),F(1,10**30)
    cy=F(24780800000000000000000,49)+16316
    v=16*F(11,7)*33536000000000000000+F(558,7)
    mstar=3*M-v
    p1=268+94*A; a1=1/sig+p1
    theta=F(11,7);u1=F(1000);u2=F(10**6)
    c1=4*theta*u1
    c2=4*theta*u2+(4*theta+16*theta*theta)*u1*u1
    assert c1<7000 and c2<6*10**7
    dk=2*c1*a1+c2
    dpf_initial=2*c1*(p1+cy*s)+c2
    added=tau*dk/(2*mstar)+(tau*dpf_initial+cy*c1*tau*tau)/(2*M)
    assert added<F(2,10**40)
    source=F(12)/(F(1,10**5)*A)*(theta*33536000000000000000*512+42*8*10**6+(theta+F(11,14))*512)
    y_compare=cy*(8*sig+8/(M*sig))
    full=F(11,8750)+F(42,10**6)+source+y_compare+F(4000024,10**65)+F(144,10**6)+F(3201,10**11)+added+F(1,10**19)
    e=F(1447,10**6); q=F(35,1000)
    assert full<F(1446416765880,10**15)<e
    eta=F(2,10**6)
    response=q-4*e
    floor=q/2-2*e-eta/2
    assert floor>F(1,100)
    y=2*sig+s*a1/mstar+tau*(a1+c1)/mstar
    tail=(y/F(1,100))**2
    assert tail<F(1,10**53)
    spatial_gap=F(19,10)-F(6,100)-tau
    assert spatial_gap==F(1,25)
    data=dict(t=T,s=s,tau=tau,C1=c1,C2=c2,extra_HL_error=added,
        original_safe_HL_error=F(3201,10**11),full_recomputed_error=full,per_branch_complete_error=e,
        q_used=q,source_support_inner_radius=F(19,10),read_lab_radius=F(6,100),
        causal_gap=spatial_gap,Y_rms_upper=y,normal_Y_tail_upper=tail,
        actual_response_lower=response,parent_response_tolerance=eta,
        parent_max_probability_error_lower=floor,target=F(1,100))
    return dict(fractions={k:str(v) for k,v in data.items()},decimal={k:float(v) for k,v in data.items()})

def annihilate(n,i):
    out=np.zeros((2**n,2**n),complex)
    for state in range(2**n):
        if (state>>i)&1:
            out[state^(1<<i),state]=(-1)**((state&((1<<i)-1)).bit_count())
    return out

def exp_h(h,t):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T

def car():
    c=[annihilate(4,j) for j in range(4)]
    ident=np.eye(16);N=c[0].conj().T@c[0]
    # Common even environment (empty and doubly occupied), purified by passive R.
    xi=np.zeros((16,2),complex)
    xi[0,0]=math.sqrt(.31);xi[12,0]=1j*math.sqrt(.19)
    xi[0,1]=math.sqrt(.23);xi[12,1]=math.sqrt(.27)
    assert abs(np.linalg.norm(xi)-1)<1e-12
    j=np.stack([c[i].conj().T@xi for i in (0,1)],axis=-1)
    assert np.linalg.norm(j.reshape(32,2).conj().T@j.reshape(32,2)-np.eye(2))<1e-12
    sx=np.array([[0.,1.],[1.,0.]])
    rec=exp_h(np.kron(N,sx),math.pi/2)
    Y=1j*(c[1].conj().T@c[0]-c[0].conj().T@c[1])
    cont=exp_h(Y,math.pi/4)
    result=[]
    for q in (0.,.035,.109,.5,1.):
        u1=math.sqrt(q)*c[0]+math.sqrt(1-q)*c[2]
        u2=math.sqrt(q)*c[1]+math.sqrt(1-q)*c[3]
        W=exp_h(u1.conj().T@u2+u2.conj().T@u1,math.pi/2)
        ps=[];full=[]
        for op in (ident,W):
            row=[];cols=[]
            for col in range(2):
                phi=op@j[:,:,col]
                withp=np.zeros((16,2,2),complex);withp[:,0,:]=phi
                recorded=(rec@withp.reshape(32,2)).reshape(16,2,2)
                final=np.einsum('ab,bpr->apr',cont,recorded)
                row.append(float(np.sum(abs(final[:,1,:])**2)))
                cols.append(final.reshape(64))
            columns=np.stack(cols,axis=1)
            assert np.linalg.norm(columns.conj().T@columns-np.eye(2))<1e-12
            ps.append(row);full.append(columns)
        contrast0=ps[0][0]-ps[0][1];contrastW=ps[1][0]-ps[1][1]
        assert abs(contrast0-1)<1e-12
        residual=abs(contrastW-(1-2*q))
        assert residual<1e-12
        response=max(abs(ps[1][i]-ps[0][i]) for i in range(2))
        assert response>=q-1e-12
        result.append(dict(q=q,no_W=ps[0],with_W=ps[1],contrast=contrastW,
            contrast_identity_residual=residual,maximum_input_response=response))
    return dict(full_CAR_dimension=16,passive_reference_dimension=2,pointer_dimension=2,
        common_background_not_empty=True,continuation_included=True,rows=result,
        scope='Finite identity diagnostic; the original normal infinite sea is handled analytically.')

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-12),(a,b)
    else:assert a==b,(a,b)

def calculate():
    sha={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert sha==FROZEN
    return dict(status='PASS',scope='CONDITIONAL_FINITE_LOCAL_READOUT_REFINEMENT_OBSTRUCTION',
        frozen_sources=sha,numerical_check_groups=3,tail=load_tail(),exact=exact(),CAR=car(),
        inherited_coarse_receiver_refuted=False,all_cognitive_axioms_refuted=False,
        added_late_control_and_small_readout_contract=True,
        local_control_hardware_implementation_claimed=False,new_axioms=0,science_increment=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    data=calculate();out=HERE/'results.json'
    if args.write_results:
        if out.exists():raise FileExistsError(out)
        out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:compare(data,json.loads(out.read_text(encoding='utf-8')))
    print(json.dumps(dict(status='PASS',groups=3,frozen_sources=len(FROZEN),
        q_lower=data['tail']['q_lower'],nu_lower=data['tail']['nu_lower'],
        matching_error_floor=data['exact']['decimal']['parent_max_probability_error_lower'],
        old_contract_refuted=False),ensure_ascii=False))
