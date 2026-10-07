"""967: collective charge response of the unchanged 956/958 material.
A triangle and its sign-flipped control; no fitted three-body interaction.
Default read-only. --write creates the result once.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, importlib.util, itertools, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
OUT=HERE/"collective_material_results.json"
EDGES=((0,1),(1,2),(2,0))

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def tensor(items):
    out=np.ones((1,1))
    for a in items:out=np.kron(out,a)
    return out
def ground(a):
    w,v=np.linalg.eigh(a);return float(w[0]),v[:,0]
def operators(n,h,q):
    I=np.eye(len(h))
    H0=sum(tensor([h if i==j else I for j in range(n)]) for i in range(n))
    vs={(i,j):tensor([q if k in (i,j) else I for k in range(n)])
        for i in range(n) for j in range(i+1,n)}
    return H0,vs
def ldl_inertia(a,shift):
    """Exact rational unpivoted LDL; no zero pivots for stored endpoints."""
    n=len(a);L=[[F(0) for _ in range(n)] for _ in range(n)];d=[]
    for i in range(n):
        pivot=a[i][i]-shift-sum(L[i][k]**2*d[k] for k in range(i))
        assert pivot!=0
        d.append(pivot);L[i][i]=F(1)
        for j in range(i+1,n):
            L[j][i]=(a[j][i]-sum(L[j][k]*L[i][k]*d[k] for k in range(i)))/pivot
    return [sum(x<0 for x in d),sum(x>0 for x in d)]
def certified_ground(a):
    # Original singlet matrices have entries in (1/5) Z at the chosen inputs.
    ints=np.rint(5*a).astype(int)
    assert np.max(abs(a-ints/5))<1e-14
    exact=[[F(int(x),5) for x in row] for row in ints]
    eig,_=ground(a);scale=10**12;base=math.floor(eig*scale)
    lo=F(base-2,scale);hi=F(base+3,scale)
    il=ldl_inertia(exact,lo);ih=ldl_inertia(exact,hi)
    assert il==[0,len(a)] and ih==[1,len(a)-1]
    return dict(lo=str(lo),hi=str(hi),lower_inertia=il,upper_inertia=ih),lo,hi
def load():
    p=STAGE/"956/native_material_interface.py"
    sp=importlib.util.spec_from_file_location("m956",p)
    m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m)
    mat=m.material(1.,.1)
    cs=[m.annihilate(i) for i in range(4)]
    ns=[c.T@c for c in cs]
    sub=np.eye(16)[:,[k for k in range(16) if k.bit_count()==2]]
    q=sub.T@(ns[0]+ns[1]-ns[2]-ns[3])@sub/2
    S=np.stack([mat["s"],mat["d"],q@mat["d"]],axis=1)
    h=S.T@mat["H"]@S;qs=S.T@q@S
    h_exact=np.array([[0.,-.2,0.],[-.2,1.,0.],[0.,0.,1.]])
    q_exact=np.array([[0.,0.,0.],[0.,0.,1.],[0.,1.,0.]])
    assert norm(h-h_exact)<1e-14 and norm(qs-q_exact)<1e-14
    g=mat["W"]@np.array([0.,1.,-1.,0.])/math.sqrt(2)
    triplet=mat["single"]@np.array([[1,0,0],[0,1/math.sqrt(2),0],
                                  [0,1/math.sqrt(2),0],[0,0,1]])
    L=np.column_stack([g,triplet])
    return m,mat,q,S,h_exact,q_exact,g,triplet,L

def run():
    old,mat,q,S,h,qs,g,triplet,L=load()
    J=mat["J"];d=mat["occupancy"];delta=1+J
    cert1,lo1,hi1=certified_ground(h);E1,gs=ground(h)
    # Rational inequalities justify the conservative all-input bound.
    assert F(269,250)**2<F(29,25)<F(125,116)**2
    assert F(347,1000)**2>F(3,25)
    bound=2*F(9,250)*F(347,1000)/F(219,500)
    H02,V2=operators(2,h,qs);H2=H02+.2*V2[(0,1)]
    E2,phi2=ground(H2);cert2,lo2,hi2=certified_ground(H2)
    c2=E2-2*E1
    saved=json.loads((STAGE/"958/capacitive_material_write_results.json").read_text("utf-8"))
    assert abs(c2+saved["exact_interface"]["conditional_energy_chi"])<1e-14
    H03,V3=operators(3,h,qs)
    psi3=np.kron(np.kron(gs,gs),gs)
    P0=np.outer(psi3,psi3)
    R=np.linalg.inv(H03-3*E1*np.eye(27)+P0)-P0
    coefficient=1.5*d**3/delta**2
    Vsum=sum(V3.values())
    third=float(psi3@Vsum@R@Vsum@R@Vsum@psi3)
    assert abs(third-coefficient)<2e-14
    # Local orbital swap flips Q and leaves h invariant.
    parity=np.diag([1.,1.,-1.])
    assert np.array_equal(parity@h@parity,h)
    assert np.array_equal(parity@qs@parity,-qs)
    rows=[]
    fullH0,fullV=operators(3,mat["H"],q)
    W0=tensor([L,L,L])
    for signs in ((1,1,1),(-1,1,1)):
        ks={tuple(sorted(e)):.2*s for e,s in zip(EDGES,signs)}
        H=H03+sum(ks[e]*v for e,v in V3.items())
        E3,phi3=ground(H);cert3,lo3,hi3=certified_ground(H)
        zeta=E3-3*E2+3*E1
        zlo=lo3-3*hi2+3*lo1;zhi=hi3-3*lo2+3*hi1
        assert zlo*zhi>0 and float(zlo)<=zeta<=float(zhi)
        HF=fullH0+sum(ks[e]*v for e,v in fullV.items())
        # Build all 64 low-energy eigenvectors using exact spin sectors.
        columns=[];energies=[];bare_distances=[];subset={}
        for labels in itertools.product(range(4),repeat=3):
            active=tuple(i for i,x in enumerate(labels) if x==0);n=len(active)
            if active not in subset:
                if n==0:subset[active]=(0.,np.ones(1))
                else:
                    hh,vv=operators(n,h,qs)
                    hh=hh+sum(ks[tuple(sorted((active[i],active[j])))]*op
                              for (i,j),op in vv.items())
                    subset[active]=ground(hh)
            energy,phi=subset[active]
            factors=[S if x==0 else triplet[:,x-1:x] for x in labels]
            vec=tensor(factors)@phi
            bare=np.kron(np.kron(L[:,labels[0]],L[:,labels[1]]),L[:,labels[2]])
            if np.vdot(bare,vec).real<0:vec=-vec
            columns.append(vec);energies.append(energy)
            bare_distances.append(float(np.linalg.norm(vec-bare)))
        W=np.column_stack(columns);K=np.diag(energies)
        iso=norm(W.T@W-np.eye(64));inter=norm(HF@W-W@K)
        distance=norm(W-W0)
        assert max(iso,inter)<2e-13 and 2*distance<float(bound)
        # Connected energy in the SAME H's record sectors.
        energy000=energies[0]
        assert abs(energy000-E3)<1e-14
        Kadd=[]
        for labels in itertools.product(range(4),repeat=3):
            active=[i for i,x in enumerate(labels) if x==0]
            Kadd.append(len(active)*E1+len(active)*(len(active)-1)/2*c2)
        residual=K-np.diag(Kadd)
        target=np.zeros((64,64));target[0,0]=zeta
        assert norm(residual-target)<2e-14
        # Hellmann-Feynman sources: keep exact virtual correlations.
        sources=[];dk=2e-5
        for edge in sorted(V3):
            op=V3[edge]
            current=float(phi3@op@phi3)
            pairH=H02+ks[edge]*V2[(0,1)]
            _,pairphi=ground(pairH)
            pairsource=float(pairphi@V2[(0,1)]@pairphi)
            fd=(ground(H+dk*op)[0]-ground(H-dk*op)[0])/(2*dk)
            assert abs(fd-current)<2e-9
            sources.append(dict(edge=list(edge),collective_source=current,
                                isolated_pair_source=pairsource,
                                connected_source=current-pairsource,
                                derivative_check=fd))
        # Same old relation code, two singlet controls and one + receiver.
        b=np.eye(16)
        code=np.stack([(b[5]-b[6]-b[9]+b[10])/2,
          (2*b[3]+2*b[12]-b[5]-b[6]-b[9]-b[10])/math.sqrt(12)],axis=1)
        local=np.kron(mat["W"],np.eye(4))@code
        initial=np.kron(np.kron(g,g),(local[:,0]+local[:,1])/math.sqrt(2))
        time=math.pi/abs(zeta)
        eig,ev=np.linalg.eigh(HF)
        assert abs(eig[0]-E3)<2e-14
        assert abs(np.trace(HF)-216)<1e-12
        evolution=(ev*np.exp(-1j*np.remainder(eig*time,2*math.pi)))@ev.T
        out=(evolution@initial.reshape(216,4)).reshape(6,6,6,4)
        targetC=(np.exp(-1j*np.remainder((E1+2*c2)*time,2*math.pi))*local[:,0]
                  +local[:,1])/math.sqrt(2)
        amplitude=np.einsum("abcd,cd->ab",out,targetC.reshape(6,4).conj())
        probability=float(np.linalg.norm(amplitude)**2)
        # Exact ideal time and effect define the analytic statement. The
        # separately saved rational energy intervals quantify phase calibration.
        phase_guard=time*(float(zhi-zlo)+float(hi1-lo1)+2*float(hi2-lo2)
                          +4*float(hi1-lo1))
        contrast_lower=1-2*float(bound)-phase_guard
        assert probability<1-contrast_lower
        exact_operator_error=norm(evolution@W0-W0@np.diag(
            np.exp(-1j*np.remainder(np.array(energies)*time,2*math.pi))))
        assert exact_operator_error<float(bound)+1e-7
        # Threefold orbital role flips change two edges, preserving the loop.
        conjugator=tensor([parity,np.eye(3),np.eye(3)])
        ks_flip={edge:(-k if 0 in edge else k) for edge,k in ks.items()}
        Hflip=H03+sum(ks_flip[e]*v for e,v in V3.items())
        assert norm(conjugator@H@conjugator-Hflip)<1e-14
        # Finite energy exchange is from the same H, not an extra force fit.
        parts=[tensor([mat["H"] if i==j else np.eye(6) for j in range(3)])
               for i in range(3)]
        parts.extend(ks[e]*v for e,v in fullV.items())
        balance=norm(sum(1j*(HF@p-p@HF) for p in parts))
        assert balance<1e-12
        rows.append(dict(signs=list(signs),loop_product=float(np.prod(signs)),
            full_trace=float(np.trace(HF)),full_ground_energy=float(eig[0]),
            all_singlet_energy=E3,ground_interval=cert3,connected_energy=zeta,
            connected_energy_interval=dict(lo=str(zlo),hi=str(zhi)),
            third_order_prediction=coefficient*math.prod(ks.values()),
            sources=sources,dressed_embedding_distance=distance,
            all_time_uniform_input_error_bound=str(bound),intertwining_error=inter,
            isometry_error=iso,direct_uniform_error=exact_operator_error,
            pair_only_prediction_probability=1.,actual_receiver_probability=probability,
            read_time=time,energy_interval_phase_guard=phase_guard,
            certified_prediction_discrepancy_lower=contrast_lower,
            energy_current_balance=balance))
    assert rows[0]["connected_energy"]>0 and rows[1]["connected_energy"]<0
    difflo=F(rows[0]["ground_interval"]["lo"])-F(rows[1]["ground_interval"]["hi"])
    diffhi=F(rows[0]["ground_interval"]["hi"])-F(rows[1]["ground_interval"]["lo"])
    assert difflo>0
    files=[STAGE/"956/native_material_interface.py",
           STAGE/"958/capacitive_material_write.py",
           STAGE/"958/capacitive_material_write_results.json",
           STAGE/"research_note_960.md",STAGE/"research_note_966.md",
           STAGE.parent/"archive_429_466/research_note_435.md",
           HERE/"drafts/collective_material_decision.md"]
    return dict(round=967,all_scientific_checks_passed=True,
        parameters=dict(U=1.,v=.1,absolute_kappa=.2,active_dimension=216,
                        full_with_old_spectator_spins=13824),
        one_body_energy=E1,pair_energy=E2,pair_conditional_energy=c2,
        one_body_interval=cert1,pair_interval=cert2,
        third_order_loop_coefficient=coefficient,independent_resolvent_coefficient=third,
        fixed_menu_rows=rows,loop_sign_energy_difference_interval=dict(lo=str(difflo),hi=str(diffhi)),
        all_time_input_bound_float=float(bound),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files},
        references=["https://arxiv.org/abs/1201.1532"],
        scope=dict(same_material_family_and_density_rule=True,
            independent_three_body_coupling_introduced=False,
            exact_pair_spectrum_reused=True,
            all_declared_low_inputs_and_passive_references=True,
            unchanged_two_body_read_protocol_under_triangle=False,
            autonomous_geometry_or_binding_proved=False,
            actual_field_process_965_embedded=False,
            natural_gauge_group_or_spatial_dimension_derived=False,
            full_goal_completed=False,visual_checks_performed=False))
def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f"/{i}")
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-6,abs_tol=2e-8),(path,a,b)
    else:assert a==b,(path,a,b)
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true")
    args=parser.parse_args();result=run()
    if args.write:
        with OUT.open("x",encoding="utf-8") as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else:compare(result,json.loads(OUT.read_text("utf-8")))
    print(json.dumps(dict(round=967,passed=True,uniform_bound=result["all_time_input_bound_float"],
        rows=[{k:r[k] for k in ("signs","connected_energy","third_order_prediction","read_time",
           "actual_receiver_probability","certified_prediction_discrepancy_lower")}
           for r in result["fixed_menu_rows"]]),ensure_ascii=False))

