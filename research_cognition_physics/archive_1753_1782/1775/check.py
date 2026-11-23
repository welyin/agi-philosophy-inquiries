"""Count chiral weak doublets from windowed charged and neutral currents."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE / "inputs.json").read_bytes())


def weight(matrix):
    return float(np.vdot(matrix,matrix).real)


def visible(A, B):
    qu, qd = A@A.conj().T, B@B.conj().T
    w = A@B.conj().T
    gram = np.block([[qu,w],[w.conj().T,qd]])
    return qu,qd,w,gram


def strength_budget(U, D, nu, nd):
    """U,D are the full isometric embeddings of the doublet multiplicity."""
    C = U@D.conj().T
    qu,qd,w,gram = visible(U[:nu],D[:nd])
    blocks = [weight(C[:nu,:nd]),weight(C[:nu,nd:]),
              weight(C[nu:,:nd]),weight(C[nu:,nd:])]
    vu,vd = np.trace(qu).real, np.trace(qd).real
    np.testing.assert_allclose([vu,vd,sum(blocks)], [sum(blocks[:2]),blocks[0]+blocks[2],
                                U.shape[1]],atol=2e-14,rtol=0)
    gap_u=qu-w@w.conj().T
    gap_d=qd-w.conj().T@w
    np.testing.assert_allclose(gap_u,C[:nu,nd:]@C[:nu,nd:].conj().T,atol=4e-15)
    np.testing.assert_allclose(gap_d,C[nu:,:nd].conj().T@C[nu:,:nd],atol=4e-15)
    eig = np.linalg.eigvalsh(gram)
    trace2=weight(qu)+weight(qd)+2*weight(w)
    effective=(vu+vd)**2/trace2 if trace2 else 0.
    return dict(visible_visible=blocks[0],visible_up_only=blocks[1],visible_down_only=blocks[2],
                both_hidden=blocks[3],sum_touching_visible=sum(blocks[:3]),
                up_neutral_trace=vu,down_neutral_trace=vd,
                Gram_rank=int(sum(eig>1.e-10)),phase_independent_effective_rank=effective)


def complete_rows(A, n, hidden_count):
    """Append hidden rows to turn an observed contraction into an isometry."""
    A=np.pad(A,((0,0),(0,n-A.shape[1])))
    gap=np.eye(n)-A.conj().T@A
    eig,vec=np.linalg.eigh(gap)
    assert min(eig)>-1.e-12
    nonzero=eig>1.e-12
    rows=np.sqrt(eig[nonzero,None])*vec[:,nonzero].conj().T
    assert len(rows)<=hidden_count
    full=np.vstack([A,rows,np.zeros((hidden_count-len(rows),n))])
    np.testing.assert_allclose(full.conj().T@full,np.eye(n),atol=1.e-12)
    return full


def run():
    old=runpy.run_path(str(HERE/CFG["registered_threshold_program"]))
    old_data=json.loads((HERE/CFG["registered_threshold_data"]).read_bytes())["matching"]
    _,yd,V,_,_=old["threshold"]()
    c=abs(V[0,0])
    m=yd*old["VEV"]/math.sqrt(2)
    eps=np.array([old["pole"](m[0],c)["lost_doublet_probability"],0.,
                  old["pole"](m[2],math.sqrt((1+c)/2))["lost_doublet_probability"]])
    np.testing.assert_allclose(eps[[0,2]],[old_data[k]["lost_doublet_probability"]
                               for k in ("gate_d","gate_b")],rtol=1.e-9,atol=0)
    # Up mass rows carry the native V; the two omitted down mass rows recover the loss.
    D=np.vstack([np.diag(np.sqrt(1-eps)),[math.sqrt(eps[0]),0,0],[0,0,math.sqrt(eps[2])]])
    np.testing.assert_allclose(D.conj().T@D,np.eye(3),atol=2.e-15)
    windows={}
    for name in ("ordinary_visible_window","artificial_endpoint_window"):
        win=CFG[name]
        windows[name]=strength_budget(V,D,win["up"],win["down"])
    np.testing.assert_allclose(windows["ordinary_visible_window"]["visible_up_only"],
                               eps.sum(),rtol=1.e-12,atol=0)
    # Add a hidden ordinary or mirror quark sector, or one vectorlike doublet, at current level.
    plus_up=np.eye(4,dtype=complex)
    plus_up[:3,:3]=V
    plus_down=np.zeros((6,4),complex)
    plus_down[:5,:3]=D
    plus_down[-1,3]=1
    right_up=np.zeros((4,1))
    right_down=np.zeros((6,1))
    right_up[-1,0]=right_down[-1,0]=1
    ordinary=windows["ordinary_visible_window"]
    appended_L=strength_budget(plus_up,plus_down,3,3)
    appended_R=strength_budget(right_up,right_down,3,3)
    for field in ("visible_visible","visible_up_only","visible_down_only","up_neutral_trace","down_neutral_trace","Gram_rank"):
        assert abs(ordinary[field]-appended_L[field])<2.e-14
        assert appended_R[field]==0
    assert appended_L["both_hidden"]==appended_R["both_hidden"]==1
    config=CFG["two_rows_completion_example"]
    a=np.array([config["visible_up_row"]])
    b=np.array([config["visible_down_row"]])
    slots=config["hidden_mass_slots_per_charge"]
    demo={}
    for n in range(2,slots+1):
        A=complete_rows(a,n,slots)
        B=complete_rows(b,n,slots)
        demo[n]=strength_budget(A,B,1,1)
        assert abs(demo[n]["sum_touching_visible"]-1)<1.e-14 and demo[n]["Gram_rank"]==2
    # Spin 1, Y=2/3: the down-to-up link has Clebsch sqrt(2) and connects m=-1 to m=0.
    jump=np.diag([math.sqrt(2),math.sqrt(2)],k=-1)
    np.testing.assert_allclose(jump@jump.T-jump.T@jump,np.diag([-2,0,2]),atol=1.e-14)
    return dict(matching_scale_GeV=old_data["mass_GeV"],registered_missing_down_weight=eps.tolist(),
                windows_L=windows,ordinary_window_R_all_strengths=0,
                direct_sum_fixed_observations=dict(baseline_index=3,extra_sequential=4,extra_mirror=2,
                                                   extra_vectorlike=3,single_new_hidden_current_norm=1),
                scalar_counterexample=demo,triplet_middle_T3=0,triplet_ordinary_link_weight=2,
                top_only_conditional_R_strength_max=max(abs(x) for x in CFG["Wtb_single_real_parameter_95"]["range"])**2)


def tidy(x):
    if isinstance(x,dict):
        return {str(k):tidy(v) for k,v in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.13g}")
    return x


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    verify=cli.parse_args().check
    result=tidy(run())
    path=HERE/"results.json"
    if verify:
        assert json.loads(path.read_bytes())==result
        print("1775 verified: two finite windows, partial-isometry completions and hidden chiral partners")
    else:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
