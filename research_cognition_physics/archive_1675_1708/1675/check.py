"""Invert the three-family charged mass conditions for the CKM phase."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())
LOW = runpy.run_path(str(HERE/CFG["matrix"]))
IV = runpy.run_path(str(HERE/CFG["interval"]))
I = IV["I"]


def read(path):
    return json.loads((HERE/path).read_bytes())


qdata = read(CFG["quarks"])["high_scale"]
qrows = np.array(qdata["down"] + qdata["angles_theta12_theta13_theta23_delta"][:3])
lrows = np.array(read(CFG["charged_leptons"])["scales"]["1e16"][3:])
qcenter, lcenter = qrows[:,0], lrows[:,0]


def ratio_shapes(x, quarks, leptons):
    edges = LOW["texture"]([*quarks, math.acos(x)])
    a,b,c,e,f,zr,zi = edges
    K = a*a+b*b+9*c*c
    L = e**2+f**2+zr**2+zi**2
    h = abs(a*b*f)
    observed = np.array(leptons)**2
    target = float(np.prod(leptons)/sum(observed)**1.5)
    target_p = sum(observed[i]*observed[j] for i,j in ((0,1),(1,2),(0,2)))/sum(observed)**2
    alpha = (target * K * np.sqrt(L)/h)**2
    roots = np.roots([1,-2,1,-alpha])
    ws = sorted(z.real for z in roots if abs(z.imag)<1e-10 and 0<z.real<1)
    assert len(ws)==2
    out=[]
    for w in ws:
        for sign in (-1,1):
            t=sign*np.sqrt(w*K/((1-w)*L))
            m=LOW["mtx"](edges,[1,1,3],t)
            vals=np.linalg.svd(m,compute_uv=False)[::-1]
            P=sum((vals[i]*vals[j])**2 for i,j in ((0,1),(1,2),(0,2)))/sum(vals**2)**2
            scale=math.sqrt(sum(observed)/sum(vals**2))
            out.append((P-target_p,t,scale,(vals*scale).tolist()))
    return out


def gamma_deg(x, quarks):
    u,v,w = quarks[3:6]
    delta = math.acos(x)
    k = math.cos(u)*math.tan(w)*math.sin(v)/math.sin(u)
    return np.degrees(delta - np.angle(1+k*np.exp(1j*delta)))


def zero(f, a, b):
    fa = f(a)
    assert fa*f(b) < 0
    for _ in range(48):
        mid=(a+b)/2
        if fa*f(mid)>0:
            a=mid
        else:
            b=mid
    return (a+b)/2


def interval_domain(quarks):
    six = [I.point(v) for v in quarks]
    return [six[0]/six[2],six[1]/six[2]]+[IV["trig"](v,True) for v in six[3:]]


def central_cover():
    base = interval_domain(qcenter)
    target=IV["spectral_target"]([[float(v),0] for v in lcenter],0)
    # Search all cos(delta), both signs and both determinant branches with outward intervals.
    todo=[IV["bounds"](-1,1)]
    survivors=[[] for _ in range(4)]
    cells=0
    while todo:
        box=todo.pop()
        cells+=1
        try:
            gaps=IV["classifications"](base+[box],(1,1,3),*target)
        except (AssertionError,ZeroDivisionError):
            # Enclosures crossing a denominator or a square-root boundary are refined.
            gaps=[0.]*4
        if min(gaps)>0:
            continue
        if box.hi-box.lo > IV["iv"].SCALE//250000:
            todo.extend(box.parts(2))
        else:
            for j,gap in enumerate(gaps):
                if gap==0:
                    survivors[j].append(box.values())
    assert not any(survivors[:3]) and survivors[3]
    intervals=sorted(survivors[3])
    for a,b in zip(intervals[:-1],intervals[1:]):
        assert a[1]>=b[0]
    return [intervals[0][0],intervals[-1][1]],cells


def cp_radius_certificate():
    """Exclude both signs of a real CKM phase with nine inputs allowed to vary."""
    radius=I.point("0.58")
    def vary(row):
        center,width=map(I.point,row)
        extent=width*radius
        return center+I(-extent.hi,extent.hi)
    dom=interval_domain([vary(row) for row in qrows])
    leptons=[vary(row) for row in lrows]
    e,m=leptons[0]/leptons[2],leptons[1]/leptons[2]
    tr=1+e.square()+m.square()
    D=e*m/(tr*tr.sqrt())
    P=(e.square()*m.square()+e.square()+m.square())/tr.square()
    counts={}
    gaps={}
    for x in (-1,1):
        pending=[(dom,[0]*5)]
        count=0
        min_gap=1.
        while pending:
            box,depths=pending.pop()
            count+=1
            residual=min(IV["classifications"](box+[I.point(x)],(1,1,3),D,P))
            if residual>0:
                min_gap=min(min_gap,residual)
                continue
            assert sum(depths)<35
            coord=min(range(5),key=lambda i:depths[i]+[0,-.4,0,-.3,-.1][i])
            next_depths=depths.copy()
            next_depths[coord]+=1
            for half in box[coord].parts(2):
                sub=box.copy()
                sub[coord]=half
                pending.append((sub,next_depths))
        counts[str(x)]=count
        gaps[str(x)]=min_gap
    return {"closed_common_pull_radius":0.58,"tree_nodes_by_cos_delta":counts,
            "smallest_disjoint_interval_gap_by_cos_delta":gaps}


def witness(x,q,l,branch):
    diff,t,k,actual = ratio_shapes(x,q,l)[branch]
    assert abs(diff)<1e-11
    np.testing.assert_allclose(actual,l,rtol=1e-9,atol=1e-16)
    gamma=gamma_deg(x,q)
    return {"cos_delta":x,"delta_abs_deg":math.degrees(math.acos(x)), "t":t,
            "k_times_yb":k,"tree_gamma_mod180_for_two_signs_deg":[gamma%180,(-gamma)%180],
            "quark_yukawas_and_angles":list(q),"three_charged_yukawas":actual,
            "largest_marginal_pull":
                max(np.abs((q-qcenter)/qrows[:,1]).max(),np.abs((l-lcenter)/lrows[:,1]).max())}


def signed_bracket(candidates, branch):
    """Outward signed residual at two endpoints; double determinant roots excluded on the path."""
    def as_interval(a,b):
        return IV["bounds"](min(float(a),float(b)),max(float(a),float(b)))
    (x1,q1,l1),(x2,q2,l2)=candidates
    def calc(x,q,l,endpoint):
        tr=sum(v.square() for v in l)
        D=l[0]*l[1]*l[2]/(tr*tr.sqrt())
        P=sum(l[i].square()*l[j].square() for i,j in ((0,1),(1,2),(0,2)))/tr.square()
        c,u0,u2,u3,u4,_=IV["coefficients"](interval_domain(q)+[x],(1,1,3))
        critical=(D/c).square()
        assert critical.lo>0 and 27*critical.hi<4*IV["iv"].SCALE
        if endpoint:
            w=list(IV["level_set"](D,c))[branch//2]
            b=1-w
            # branch 2 has negative t, branch 3 has positive t.
            return (u0*b.square()+u2*b*w+u4*w.square()
                    +(1 if branch%2==0 else -1)*u3*(b*w).sqrt()*w-P).values()
    calc(as_interval(x1,x2),[as_interval(a,b) for a,b in zip(q1,q2)],
         [as_interval(a,b) for a,b in zip(l1,l2)],False)
    ends=[calc(I.point(x),[I.point(v) for v in q],[I.point(v) for v in l],True)
          for x,q,l in candidates]
    assert ends[0][1]<0<ends[1][0] or ends[1][1]<0<ends[0][0]
    return ends


def main():
    enclosure,cells=central_cover()
    x=zero(lambda z:ratio_shapes(z,qcenter,lcenter)[3][0],*enclosure)
    center=witness(x,qcenter,lcenter,3)
    # Interpolate only CP-even inputs to construct an exact real CKM counterexample.
    signs=np.array([1,1,-1,1,-1,-1,-1,1,-1])
    def walk(v):
        return qcenter+v*signs[:6]*qrows[:,1],lcenter+v*signs[6:]*lrows[:,1]
    pull=zero(lambda v:ratio_shapes(-1,*walk(v))[2][0],0,2)
    real=witness(-1,*walk(pull),2)
    real_signs=signed_bracket([(-1,*walk(.67)),(-1,*walk(.672))],2)
    # Establish simultaneous feasibility before comparison to the direct phase.
    signs=np.array([1,-1,-1,-1,1,-1,-1,1,1])
    q=qcenter+1.99*signs[:6]*qrows[:,1]
    l=lcenter+1.99*signs[6:]*lrows[:,1]
    x=zero(lambda z:ratio_shapes(z,q,l)[3][0],.3,.5)
    positive=witness(x,q,l,3)
    positive_signs=signed_bracket([(.3888,q,l),(.3889,q,l)],3)
    central_signs=signed_bracket([(.2398,qcenter,lcenter),(.2399,qcenter,lcenter)],3)
    peak,lower,upper=CFG["LHCb_gamma_2026_deg"]
    gamma=center["tree_gamma_mod180_for_two_signs_deg"][0]
    p=positive["tree_gamma_mod180_for_two_signs_deg"][0]
    assert min(center["tree_gamma_mod180_for_two_signs_deg"])>peak+2*upper
    assert peak-2*lower<p<peak+2*upper
    assert abs(real["tree_gamma_mod180_for_two_signs_deg"][0])<1e-10
    return {"CP_conservation_excluded_box":cp_radius_certificate(),
            "central_all_four_branches_cover_cells":cells,
            "central_only_surviving_branch":"larger w at fixed determinant, positive t",
            "certified_signed_shape_residuals_at_brackets":{
                "central_x_0.2398_0.2399":central_signs,
                "real_x_minus1_pull_0.670_0.672":real_signs,
                "two_width_x_0.3888_0.3889":positive_signs},
            "central_all_cos_delta_enclosure":enclosure, "central_equality_witness":center,
            "central_direct_phase_windows_by_modulus_sign_deg":
                [[gamma_deg(enclosure[1],qcenter),gamma_deg(enclosure[0],qcenter)],
                 [180-gamma_deg(enclosure[0],qcenter),180-gamma_deg(enclosure[1],qcenter)]],
            "central_positive_margin_to_two_width_deg":gamma-(peak+2*upper),
            "two_width_CP_conserving":real,"two_width_overlapping_direct_phase":positive}


def tidy(x):
    if isinstance(x,dict):return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list):return [tidy(y) for y in x]
    if isinstance(x,(float,np.floating)):return float(f"{x:.12g}")
    return x


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=tidy(main())
    path=HERE/"results.json"
    if args.write:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    else:
        assert json.loads(path.read_bytes())==result
    print(json.dumps(result,indent=2))
