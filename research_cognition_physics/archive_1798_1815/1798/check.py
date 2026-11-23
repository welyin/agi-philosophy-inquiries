"""Joint Takagi and rephasing-free up-quark moment tests at the tagged three-family boundary."""
import argparse
import hashlib
import json
import math
import runpy
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
def pack(z):
    return [float(z.real),float(z.imag)]


def frames(y):
    u,m,vh=np.linalg.svd(y)
    return u[:,::-1],m[::-1],vh.conj().T[:,::-1]


def symmetry(d,e,cells):
    u,m,v=frames(d)
    c=u.conj().T@e@v
    edges=[(0,1),(0,2),(1,2)]
    discrepancy=[abs(abs(c[i,j])-abs(c[j,i])) for i,j in edges]
    # skew Frobenius norm squared = fixed - max Re sum_ij c_ij conj(c_ji) p_j conj(p_i).
    W=[c[i,j]*c[j,i].conj() for i,j in edges]
    const=.5*sum(abs(c[i,j])**2+abs(c[j,i])**2 for i,j in edges)
    theta=np.arange(cells)*2*math.pi/cells
    p=np.exp(1j*theta)
    objective=np.real(W[0]*p)+abs(W[1]+W[2]*p.conj())
    n=int(np.argmax(objective))
    best=objective[n]
    slack=(abs(W[0])+abs(W[2]))*math.pi/cells
    norm=np.linalg.norm(e)
    return {"ordered_cross_modulus_pairs":[[abs(c[i,j]),abs(c[j,i])] for i,j in edges],
            "modulus_Frobenius_lower_as_fraction_E":np.linalg.norm(discrepancy)/(math.sqrt(2)*norm),
            "minimum_change_to_E_at_fixed_D_as_fraction_E":[max(np.linalg.norm(discrepancy)/(math.sqrt(2)*norm),math.sqrt(max(0,const-best-slack))/norm),
                                                           math.sqrt(max(0,const-best))/norm],
            "cycle_product_relative_phase_deg":np.degrees(np.angle(np.prod([c[0,1],c[1,2],c[2,0]])/
                                                                     np.prod([c[1,0],c[2,1],c[0,2]])))}


def bloch(h):
    return np.trace(h).real/2,np.array([h[0,1].real,-h[0,1].imag,(h[0,0]-h[1,1]).real/2])


def conditional(forms,j,k,base,t):
    b,vec=np.linalg.eigh(forms[base])
    assert min(b)>0
    inv=(vec/np.sqrt(b))@vec.conj().T
    s0,s=bloch(inv@forms[j]@inv)
    d0,d=bloch(inv@forms[k]@inv)
    sn=np.linalg.norm(s)
    co=(t-s0)/sn
    if abs(co)>1:
        return {"unconditional_condition_axis_range":[s0-sn,s0+sn],"input":t,"output_range":None}
    dot=d@s/sn
    mid=d0+dot*co
    wid=np.linalg.norm(d-dot*s/sn)*math.sqrt(1-co*co)
    return {"unconditional_condition_axis_range":[s0-sn,s0+sn],"input":t,"output_range":[mid-wid,mid+wid]}


def inverse_rows(forms,actual,D,C):
    mu,vec=np.linalg.eigh(forms[2])
    inv=(vec/np.sqrt(mu))@vec.conj().T
    terms=[bloch(inv@m@inv) for m in forms[:2]]
    planes=np.array([v for _,v in terms])
    offsets=actual[:2]/actual[2]-np.array([x for x,_ in terms])
    n0=np.linalg.lstsq(planes,offsets,rcond=None)[0]
    axis=np.cross(*planes)
    assert np.linalg.norm(axis)>1e-12 and n0@n0<1
    move=axis/np.linalg.norm(axis)*math.sqrt(1-n0@n0)
    roots=[]
    for sign in (1,-1):
        x,y,z=n0+sign*move
        P=np.array([[1+z,x-1j*y],[x+1j*y,1-z]])/2
        _,vv=np.linalg.eigh(P)
        a,b=inv@vv[:,-1]*math.sqrt(actual[2])
        Y=a*D+b*C
        ul,s,_=frames(Y)
        np.testing.assert_allclose(np.sum(abs(Y)**2,axis=1),actual,rtol=3e-7)
        roots.append({"b_over_a":pack(b/a),"up_singular_values":s.tolist(),
                      "up_to_top":s[0]/s[2],"charm_to_top":s[1]/s[2],
                      "CKM_moduli":abs(ul.conj().T).tolist()})
    return roots


def up_test(d,e,up):
    l,md,r=frames(d)
    # Coefficient spinor (a,b). Each down-left row norm is a 2x2 Hermitian quadratic form.
    D=np.diag(md)
    C=l.conj().T@e@r
    forms=[]
    for i in range(3):
        arms=np.column_stack((D[i],C[i]))
        forms.append(arms.conj().T@arms)
    actual=up**2 @ abs(l)**2
    out={}
    for j,k in [(1,0),(0,1)]:
        row=conditional(forms,j,k,2,actual[j]/actual[2])
        row["actual_output"]=actual[k]/actual[2]
        out[str(j)+"_to_"+str(k)]=row
    rownorm=abs(d[1,2])/np.linalg.norm(d[2,:])
    k=(rownorm**2-(up[1]/up[2])**2)
    return {"CKM_moduli":abs(l).tolist(),"high_up_Yukawas":up.tolist(),
            "actual_charm_to_top":up[1]/up[2],
            "down_axis_moments_normalized_to_b":(actual/actual[2]).tolist(),
            "conditional_ranges":out,"tagged_row_c_to_t_norm_floor":rownorm,
            "all_solutions_to_three_actual_norms":inverse_rows(forms,actual,D,C),
            "uniform_strange_moment_and_charm_test":whole_pencil(D,C),
            "minimum_top_tilt_in_tagged_up_frame_deg":np.degrees(np.arcsin(math.sqrt(max(0,k)/(1+k))))}


def whole_pencil(d,e):
    """Certify disjoint windows on the full complex line by a finite Lipschitz cover."""
    matrices=np.stack([d,e]).reshape(2,-1)
    eigen,axes=np.linalg.eigh(matrices.conj()@matrices.T)
    whiten=(axes/np.sqrt(eigen))@axes.conj().T
    basis=(whiten.T@matrices).reshape(2,3,3)
    gram=basis.reshape(2,9).conj()@basis.reshape(2,9).T
    np.testing.assert_allclose(gram,np.eye(2),atol=2e-12)
    lipschitz=math.sqrt(np.linalg.eigvalsh(gram)[-1])+1e-12
    lower,upper=np.sqrt(CFG["robust_strange_to_bottom_moment_window"])
    charm=CFG["robust_charm_to_top_ceiling"]
    cells=[(math.pi/4,math.pi,math.pi/4,math.pi)]
    counts=[]
    least_margin=math.inf
    for level in range(22):
        t,p,a,b=np.array(cells).T
        w=np.cos(t)[:,None,None]*basis[0]+(np.sin(t)*np.exp(1j*p))[:,None,None]*basis[1]
        # The operator-norm changes of W and of a row are at most the Frobenius change.
        err=lipschitz*(a+b*np.sin(np.minimum(math.pi/2,t+a)))+1e-12
        s=np.linalg.svd(w,compute_uv=False)
        rn=np.linalg.norm(w,axis=2)
        gap_charm=s[:,1]-err-charm*(s[:,0]+err)
        gap_row=np.maximum(rn[:,1]-err-upper*(rn[:,2]+err),lower*(rn[:,2]-err)-rn[:,1]-err)
        fails_charm=gap_charm>0
        fails_row=gap_row>0
        keep=~(fails_charm|fails_row)
        if (~keep).any():
            least_margin=min(least_margin,np.maximum(gap_charm,gap_row)[~keep].min())
        counts.append([len(cells),int(sum(fails_charm)),int(sum(fails_row&~fails_charm))])
        if not keep.any():
            return {"window":CFG["robust_strange_to_bottom_moment_window"],"excluded_charm_ceiling":charm,
                    "tested_cells":sum(x[0] for x in counts),"depth":level,
                    "smallest_certified_margin_after_roundoff_allowance":least_margin,
                    "implied_window_if_Vts_033_to_05_Vtb_at_least_099":[.033**2,(.05**2+charm**2)/.99**2],
                    "cells_charm_reject_then_row_reject":counts}
        cells=[]
        for x,y,c,f in zip(t[keep],p[keep],a[keep]/2,b[keep]/2):
            cells.extend((x+i*c,y+j*f,c,f) for i in (-1,1) for j in (-1,1))
    raise AssertionError("The declared windows overlap or need a finer cover")


def calculate():
    old=runpy.run_path(str(HERE/CFG["charged_frame"]))
    _,rep=old["charged_background"](CFG["step"])
    unpack=old["unpack"]
    d=unpack(rep["down_texture_10PeV"])
    e=unpack(rep["charged_texture_10PeV"]).T
    up=np.array(rep["up_texture_diag_10PeV"])
    result=symmetry(d,e,CFG["cycle_grid_cells"])
    rng=np.random.default_rng(1798)
    def unitary():
        z=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        return np.linalg.qr(z)[0]
    l,r=unitary(),unitary()
    again=symmetry(l@d@r,l@e@r,CFG["cycle_grid_cells"])
    np.testing.assert_allclose(result["ordered_cross_modulus_pairs"],again["ordered_cross_modulus_pairs"],rtol=2e-11)
    symmetric=unitary()
    synthetic=symmetry(l@np.diag([.1,.3,1])@r,l@(symmetric+symmetric.T)@r,4096)
    assert max(abs(a-b) for a,b in synthetic["ordered_cross_modulus_pairs"])<1e-12
    pair=result["ordered_cross_modulus_pairs"][-1]
    result["sb_chirality_asymmetry_at_tagged_scale"]=(pair[0]**2-pair[1]**2)/(pair[0]**2+pair[1]**2)
    return {"symmetric_pair":result,"up_pencil":up_test(d,e,up)}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    result=clean(calculate())
    data=json.dumps(result,indent=2)+"\n"
    path=HERE/"results.json"
    if args.check:
        assert result==json.loads(path.read_bytes())
    else:
        path.write_text(data,encoding="utf-8")
    print("1798",hashlib.sha256(data.encode()).hexdigest()[:12])
