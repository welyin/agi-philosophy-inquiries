"""863 working: original-color loop-current Ward identity with smooth anchors.
A compact bump Fourier multiplier calibrates smearing; the full gravitational
reference current is treated in the working proof, not simulated by this code.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'smeared_loop_current_probe_results.json'
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(HERE.parent/'862'))
from research_layout import Layout,ResearchRuntime
from material_link_path_probe import exp_anti as exp

def run():
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        color=old.color(1.);generators=-1j*old.T
    Ay=-1j*color['A'][1];Az=-1j*color['A'][2];ell=.8
    vertices=np.array([[0,0],[ell,0],[ell,ell],[0,ell],[0,0]],float)
    increments=np.diff(vertices,axis=0)
    Alist=[Ay*d[0]+Az*d[1] for d in increments]
    pieces=[exp(-a) for a in Alist]
    def product(items):
        out=np.eye(3,dtype=complex)
        for A in items:out=A@out
        return out
    nodes,weights=np.polynomial.legendre.leggauss(64);us=(nodes+1)/2;weights=weights/2
    def bump_hat(k):
        # C-infinity bump exp(-1/(1-x^2)) on [-1,1], normalized.
        v=np.exp(-1/(1-nodes*nodes))
        return sum((2*weights)*v*np.exp(1j*k*nodes))/sum((2*weights)*v)
    rows=[];max_good=0.;max_bad=0.
    for omega,ky,kz in ((.5,.7,-.4),(1.,1.3,.9),(.2,-.3,1.7)):
        multiplier=bump_hat(.15*omega)*bump_hat(.15*ky)*bump_hat(.15*kz)
        good=[];bad=[]
        for T in generators:
            total=0j;partial_only=0j
            for j,Adot in enumerate(Alist):
                before=product(pieces[:j]);after=product(pieces[j+1:])
                dx=increments[j]
                for u,w in zip(us,weights):
                    U0=exp(-Adot*u)@before;U1=after@exp(-Adot*(1-u))
                    x=vertices[j]+u*dx;phase=np.exp(1j*(ky*x[0]+kz*x[1]))
                    c=float(np.trace(U1@T@U0).real)
                    cp=float(np.trace(U1@(Adot@T-T@Adot)@U0).real)
                    partial=1j*(ky*dx[0]+kz*dx[1])*c
                    total+=w*phase*(partial+cp);partial_only+=w*phase*partial
            good.append(abs(multiplier*total));bad.append(abs(multiplier*partial_only))
        max_good=max(max_good,max(good));max_bad=max(max_bad,max(bad))
        rows.append(dict(frequency=[omega,ky,kz],compact_anchor_multiplier_real=float(multiplier.real),
            full_covariant_Ward_residual=float(max(good)),omitting_commutator_defect=float(max(bad))))
    assert max_good<1e-13 and max_bad>1e-5
    return dict(kind='round_863_working_loop_current',formal_reports=862,new_numbered_scientific_groups=0,
        all_checks_passed=True,original_753_color_background_used=True,compact_smooth_bump_anchor=True,
        rows=rows,full_ward_maximum=float(max_good),omitted_commutator_maximum=float(max_bad),
        full_gravitational_reference_current_numerically_implemented=False,
        nonzero_original_physical_quantum_variance_proved=False,
        all_order_interacting_relational_Wilson_operator_proved=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
