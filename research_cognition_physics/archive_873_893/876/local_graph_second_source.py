"""876: same product cutoff, own Gibbs preparation, ordered records and CTP jets.
Reuses704's derivative jet algebra on706/875's original declared finite diagnostic.
Analytic claims concern the full fixed graph, not a spatial continuum simulation.
"""
from pathlib import Path
import argparse,itertools,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'local_graph_second_source_results.json'
KEYS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
BETA=1.2;HBAR=.7;TIMES=(.12,.19)

def local(d,bands):
    rank=4*bands;ix=np.array([a*d['nloc']+b for a in range(rank) for b in range(rank)])
    parts=[m[np.ix_(ix,ix)] for m in d['parts']]
    H=sum(parts);G=sum(k*m for k,m in zip(d['powers'],parts))
    C=sum(k*k*m for k,m in zip(d['powers'],parts))
    n=np.diag(np.tile([0.,1.,0.,1.],bands));I=np.eye(rank)
    p=np.kron(n,I);q=np.kron(I,n);eye=np.eye(len(H))
    return dict(H=H,G=G,C=C,parts=parts,ix=ix,records=((p,eye-p),(q,eye-q)))
def prep(x,jets):
    shift=float(np.linalg.eigvalsh(x['H']).min())
    return jets.geometry.statejets(jets.expjet(x['H']-shift*np.eye(len(x['H'])),x['G'],x['C'],BETA))[0]
def history_jets(x,rho,jets,plus=1.,minus=-.4,contact=True):
    C=x['C'] if contact else np.zeros_like(x['C'])
    evolution=[(jets.expjet(x['H'],plus*x['G'],plus**2*C,1j*t/HBAR),
                jets.expjet(x['H'],minus*x['G'],minus**2*C,1j*t/HBAR)) for t in TIMES]
    out=[]
    for hist in itertools.product((0,1),repeat=2):
        state={key:rho[key[0]].copy() if key[1]==0 else np.zeros_like(rho[0]) for key in KEYS}
        for step,(r,(up,um)) in enumerate(zip(hist,evolution)):
            updated={};L=x['records'][step][r]
            for a,b in KEYS:
                total=np.zeros_like(rho[0])
                for left in range(b+1):
                    for middle in range(b-left+1):
                        right=b-left-middle
                        factor=math.factorial(b)/(math.factorial(left)*math.factorial(middle)*math.factorial(right))
                        total+=factor*up[left]@state[(a,middle)]@um[right].conj().T
                updated[(a,b)]=L@total@L.conj().T
            state=updated
        out.append([complex(np.trace(state[key])) for key in KEYS])
    return np.array(out)
def expm(H,z):
    e,v=np.linalg.eigh(H);return (v*np.exp(z*e))@v.conj().T
def values(d,x,xi,epsilon,plus=1.,minus=-.4):
    def h(s):return sum(np.exp(k*s)*a for k,a in zip(d['powers'],x['parts']))
    H=h(xi);shift=float(np.linalg.eigvalsh(H).min())
    rho=expm(H-shift*np.eye(len(H)),-BETA);rho/=np.trace(rho)
    steps=[(expm(h(plus*epsilon),-1j*t/HBAR),expm(h(minus*epsilon),-1j*t/HBAR)) for t in TIMES]
    out=[]
    for hist in itertools.product((0,1),repeat=2):
        y=rho.copy()
        for i,(r,(up,um)) in enumerate(zip(hist,steps)):
            L=x['records'][i][r];y=L@up@y@um.conj().T@L.conj().T
        out.append(np.trace(y))
    return np.array(out)
def carray(a):
    return np.stack((np.asarray(a).real,np.asarray(a).imag),axis=-1).tolist()
def tnorm(a):return float(np.linalg.svd(a,compute_uv=False).sum())
def run():
    with ResearchRuntime(Layout()).installed():
        import joint_local_history_source_limit as old
        import joint_preparation_history_limit as jets
        d=old.data()
        fullx=local(d,d['r']);fullrho=prep(fullx,jets)
        full=history_jets(fullx,fullrho,jets)
        S=sum(d['parts'][:2])+np.eye(len(d['H']))
        rows=[];medium=None
        for bands in range(1,d['r']+1):
            x=fullx if bands==d['r'] else local(d,bands)
            rho=fullrho if bands==d['r'] else prep(x,jets)
            actual=full if bands==d['r'] else history_jets(x,rho,jets)
            weighted=[]
            for j in range(3):
                embed=np.zeros_like(fullrho[j]);embed[np.ix_(x['ix'],x['ix'])]=rho[j]
                weighted.append(tnorm(S@(embed-fullrho[j])@S))
            same=history_jets(x,rho,jets,plus=1.,minus=1.)
            normalization=float(np.max(np.abs(same.sum(axis=0)-np.array([1,0,0,0,0,0]))))
            assert normalization<3e-9 and np.min(same[:,0].real)>0
            rows.append(dict(bands=bands,dimension=len(x['H']),
                all_ctp_source_jet_errors=np.max(np.abs(actual-full),axis=0).tolist(),
                comparison_weight_sandwiched_preparation_derivative_errors=weighted,
                same_history_joint_derivative_normalization_error=normalization,
                own_preparation_trace_derivative_errors=[float(abs(np.trace(rho[j])-(1 if j==0 else 0))) for j in range(3)]))
            if bands==2:medium=(x,rho,actual)
        assert max(rows[-1]['all_ctp_source_jet_errors'])==0
        x,rho,actual=medium
        zero=values(d,x,0.,0.)
        fd=[]
        for h in (.0008,.0004):
            pp=values(d,x,h,h);pm=values(d,x,h,-h)
            mp=values(d,x,-h,h);mm=values(d,x,-h,-h)
            mixed=(pp-pm-mp+mm)/(4*h*h)
            dyn=(values(d,x,0.,h)-2*zero+values(d,x,0.,-h))/(h*h)
            pre=(values(d,x,h,0.)-2*zero+values(d,x,-h,0.))/(h*h)
            fd.append(dict(step=h,mixed_derivative_error=float(np.max(np.abs(mixed-actual[:,4]))),
                dynamical_second_derivative_error=float(np.max(np.abs(dyn-actual[:,5]))),
                preparation_second_derivative_error=float(np.max(np.abs(pre-actual[:,3])))))
        assert max(fd[-1][k] for k in fd[-1] if k!='step')<2e-3
        for k in ('mixed_derivative_error','dynamical_second_derivative_error','preparation_second_derivative_error'):
            assert fd[-1][k]<fd[0][k]/3
        wrong=history_jets(x,rho,jets,contact=False)
        contact_error=float(np.max(np.abs(wrong[:,5]-actual[:,5])))
        mixed_size=float(np.max(np.abs(actual[:,4])))
        assert contact_error>1e-3 and mixed_size>1e-3
        return dict(round=876,date='2026-10-06',formal_reports=876,
            cumulative_numbered_groups=3661,fresh_numbered_groups=1,all_checks_passed=True,
            original_fixed_graph_theorem_analytic=True,
            numerical_scope='Inherited706/875 frozen-link radial and neutral-CAR 256-dimensional process; all original five geometry-dependent terms retained.',
            source_jet_order=[list(k) for k in KEYS],derivatives_not_Taylor_coefficients=True,
            original_geometry_powers=list(d['powers']),product_cutoff_checks=rows,
            full_ordered_history_source_jets=carray(full),
            independent_eigen_exponential_finite_difference_checks=fd,
            omitted_dynamical_Hessian_contact_second_response_error=contact_error,
            frozen_preparation_would_lose_mixed_response=mixed_size,
            uniform_graph_norm_and_sandwiched_preparation_transport_proved=True,
            total_order_two_same_source_history_match_proved=True,
            exact_all_scale_dynamics_or_original_Q_E_identification_proved=False,
            full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    result=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
