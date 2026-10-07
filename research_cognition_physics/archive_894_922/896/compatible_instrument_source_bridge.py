"""896: physical CAR instrument approximation, original joint-source diagnostic.
Analytic Sobolev and constrained-development proof is in the report. Numeric
checks reuse original732's full two-momentum128 local symbol, not a curved PDE.
"""
from pathlib import Path
import json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
with ResearchRuntime(Layout()).installed():
    import joint_relative_source_development as old
TARGET=HERE/'compatible_instrument_source_bridge_results.json'

def projector(f,C):
    cf=C@f.conj();F=np.column_stack([f,cf]);Q=F@F.conj().T
    assert np.max(abs(F.conj().T@F-np.eye(2)))<1e-13
    return Q

def trace_norm(A):return float(np.sum(abs(np.linalg.eigvalsh((A+A.conj().T)/2))))

def source_vertices(bg):
    zero=[np.zeros_like(x) for x in bg];vertices=[]
    for i,j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2)):
        v=[x.copy() for x in zero];v[0][i,j]=1;v[0][j,i]=1
        if i!=j:v[0]/=np.sqrt(2)
        vertices.append((f'metric_{i}{j}',old.matrix_and_source(*bg,*v)[1]))
    for j in range(5):
        v=[x.copy() for x in zero];v[1][j]=1
        vertices.append((f'scalar_{j}',old.matrix_and_source(*bg,*v)[1]))
    for i in range(3):
        for j in range(3):
            v=[x.copy() for x in zero];v[2][i,j]=1
            vertices.append((f'weak_{i}{j}',old.matrix_and_source(*bg,*v)[1]))
    for j in range(3):
        v=[x.copy() for x in zero];v[3][j]=1
        vertices.append((f'hypercharge_{j}',old.matrix_and_source(*bg,*v)[1]))
    B=old.matrix_and_source(*bg,*zero)[0]
    vertices.append(('energy',B))
    return B,vertices

def gauge_pairs(bg,B):
    _,phi,a,_=bg;pairs=[]
    for kind,index,H in old.physical_generators():
        charge=old.assemble([old.old.old.bdg(H,np.zeros_like(H))]*2)
        zero=[np.zeros_like(x) for x in bg];dphi=np.zeros(5);da=np.zeros_like(a)
        X=phi[:2]+1j*phi[2:4]
        if kind=='weak':
            change=1j*(old.old.SIG[index]/2)@X
            dphi[:4]=np.r_[change.real,change.imag];da=np.cross(a,np.eye(3)[index])
        elif kind=='circle':
            change=3j*X;dphi[:4]=np.r_[change.real,change.imag]
        D=old.matrix_and_source(*bg,zero[0],dphi,da,zero[3])[1]
        assert np.max(abs(D-1j*(charge@B-B@charge)))<5e-13
        pairs.append((charge,D))
    return pairs

def run():
    P0,_,_,_,C=old.initial_modes()
    f=np.zeros(128,complex);f[30]=f[62]=1/np.sqrt(2)
    g=np.zeros(128,complex);g[31]=1/np.sqrt(2);g[63]=1j/np.sqrt(2)
    assert abs(np.vdot(f,g))<1e-15
    U=np.eye(128,dtype=complex);T=.04;steps=64
    for n in range(steps):U=old.unitary(n*T/steps,T/steps)@U
    P=U@P0@U.conj().T;bg=old.background(T);B,vertices=source_vertices(bg)
    vertices += [(f'charge_{j}',a) for j,(a,_) in enumerate(gauge_pairs(bg,B))]
    vp=np.array([a for _,a in vertices]);vnorm=np.array([np.linalg.norm(a,2) for a in vp])
    Q=projector(f,C);R=np.eye(128)-2*Q;d0=(R@P0@R-P0)/2;d0=U@d0@U.conj().T
    reference=np.einsum('aij,ji->a',vp,d0).real/2
    read=np.zeros(128,complex);read[26]=1
    # A fixed physical electron readout after the original actual propagation.
    rows=[]
    for angle in (.1,.03,.01,.003):
        fl=np.cos(angle)*f+np.sin(angle)*g;delta=float(np.linalg.norm(fl-f))
        ql=projector(fl,C);rl=np.eye(128)-2*ql
        dl=(rl@P0@rl-P0)/2;dl=U@dl@U.conj().T
        post=P+dl;ev=np.linalg.eigvalsh(post)
        assert ev.min()>-2e-12 and ev.max()<1+2e-12
        difference=dl-d0;dn=trace_norm(difference)
        # Rank-two Qs give ||Q_L-Q||_1 <=4||f_L-f||, hence <=8 delta
        # for the covariance difference of the actual nonselective channels.
        assert dn<=8*delta*(1+1e-12)
        sources=np.einsum('aij,ji->a',vp,dl).real/2
        errors=abs(sources-reference);bounds=4*vnorm*delta
        assert np.all(errors<=bounds+1e-13)
        ddot=-1j*(B@difference-difference@B)
        ward=max(abs(old.source(ddot,charge)+old.source(difference,D)) for charge,D in gauge_pairs(bg,B))
        assert ward<3e-12
        observation=abs(np.vdot(read,difference@read).real)
        assert observation<=4*delta+1e-13
        # Joint background-work jet, fixed independently of the record change.
        bp=old.old.point();tangent=(bp['dg'],bp['v'],bp['da'],bp['da0'])
        _,Dwork=old.matrix_and_source(*bg,*tangent)
        work=old.source(difference,Dwork)
        rows.append(dict(mode_angle=angle,mode_L2_error=delta,covariance_trace_error=dn,
            covariance_trace_error_bound=8*delta,source_errors={name:float(e) for (name,_),e in zip(vertices,errors)},
            largest_normalized_source_error=float(np.max(errors/np.maximum(vnorm,1e-15))),
            normalized_source_bound=4*delta,original_gauge_exchange_defect=ward,
            fixed_final_record_error=observation,
            joint_background_work_difference=work,
            post_covariance_lower=float(ev.min()),post_covariance_upper=float(ev.max())))
    assert rows[-1]['covariance_trace_error']<rows[0]['covariance_trace_error']
    assert max(rows[0]['source_errors'].values())>1e-4
    return dict(round=896,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3681,
        argument_scope='Approximate the same sterile CAR instrument while retaining the actual background and full Dirac evolution; source Sobolev error and constrained first-order response are proved analytically. Numerical calibration uses original732 full128 local symbol only.',
        numeric_modes=64,numeric_Nambu_dimension=128,physical_source_vertices=len(vertices),
        time=T,source_labels=[name for name,_ in vertices],
        actual_compatible_instrument_rows=rows,
        instrument_is_exact_CAR_channel=True,post_state_not_assumed_Gaussian=True,
        same_initial_reference_and_full_evolution=True,
        finite_mode_instrument_not_closed_finite_Hamiltonian=True,
        curved_future_PDE_numerically_solved=False,
        absolute_reference_source_approximated=False,
        nonlinear_feedback_completed=False,full_goal_completed=False)

if __name__=='__main__':
    result=run();assert not TARGET.exists()
    TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
