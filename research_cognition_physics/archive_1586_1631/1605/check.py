"""Two specified input rays and separately minimal rank-one Yukawa corrections."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
TAGS = "uct"


def read(rel):
    return json.loads((HERE/rel).read_text(encoding="utf-8-sig"))


def unit(a):
    return a/np.linalg.norm(a)


def plane(v):
    return np.linalg.svd(v.conj()[None,:], full_matrices=True)[2].conj().T[:,1:]


def pmatrix(t12,t13,t23,delta):
    a,b,c = np.sin([t12,t13,t23])
    A,B,C = np.cos([t12,t13,t23])
    z = np.exp(1j*delta)
    return np.array([[A*B,a*B,b/z],[-a*C-A*c*b*z,A*C-a*c*b*z,c*B],
                     [a*c-A*C*b*z,-A*c-a*C*b*z,C*B]])


def central_inputs(cfg,masses,angles,osc,atm,order):
    V = pmatrix(*angles[:,0])
    D = V*masses[:3,0]
    C = np.diag(masses[3:,0])
    params = [osc[order]["s12"][0],osc[order]["s13"][0],atm[order][0]]
    neutral = pmatrix(*np.arcsin(np.sqrt(params)),
                      math.radians(cfg["example_neutrino_phase_deg"][order]))
    n = neutral[:,0 if order=="NO" else 2]
    return V,D,C,n


def basis_from_two(z,w):
    z = unit(z)
    w = unit(w-z*np.vdot(z,w))
    third = unit(np.cross(z,w).conj())
    return np.column_stack((z,w,third))


def component_norms(vector,z):
    return np.array([abs(np.vdot(z,vector)),np.linalg.norm(vector-z*np.vdot(z,vector))])


def single_point(D,C,n,j,l,k,grid,refinements):
    eye=np.eye(3)
    m=3-j-l
    b=unit(np.linalg.solve(D,eye[:,j]))
    b2=unit(np.linalg.solve(D,eye[:,l]))
    b2*=np.exp(-1j*np.angle(np.vdot(b,b2)))
    h=float(np.vdot(b,b2).real)
    s=float(np.linalg.norm(b2-h*b))
    u=(b2-h*b)/s
    v=unit(D.conj().T[:,m])
    tag_first = s*b-h*u
    base_input=np.column_stack((b,u,v))
    assert np.linalg.norm(base_input.conj().T@base_input-eye)<1e-12

    F=plane(n)
    eigen,turn=np.linalg.eigh(F.conj().T@C@C@F)
    F=F@turn
    g=(C@F).conj().T@(C@n)
    F=F*np.exp(1j*np.angle(g))
    g=np.abs(g)
    c1,c2=np.sqrt(eigen)
    z=unit(D@v)
    r=k*np.linalg.norm(D@v)
    assert c1<r<c2
    t=math.sqrt((r*r-c1*c1)/(c2*c2-c1*c1))
    c=math.sqrt(1-t*t)
    target_second=component_norms(k*D@u,z)
    allowed_second=np.array([math.sqrt((r*r-c1*c1)*(c2*c2-r*r))/r,c1*c2/r])
    min_second=np.linalg.norm(target_second-allowed_second)/s

    K=c*t*(c2*c2-c1*c1)/r
    norm_Cu2=t*t*c1*c1+c*c*c2*c2
    target_first=component_norms(k*D@tag_first,z)

    def objective(delta,gamma):
        ed=np.exp(1j*delta)
        eg=np.exp(1j*gamma)
        zn=(c*g[0]+t*g[1]*ed.conj())/r
        nu=-t*g[0]+c*g[1]*ed
        proj=np.abs(s*zn-h*K*eg)
        norm2=s*s*np.linalg.norm(C@n)**2+h*h*norm_Cu2-2*s*h*np.real(eg*nu)
        perp=np.sqrt(np.maximum(0.,norm2-proj**2))
        return np.hypot(proj-target_first[0],perp-target_first[1])/s

    phase=np.arange(grid)*(2*math.pi/grid)
    scores=objective(phase[:,None],phase[None,:])
    ii,jj=np.unravel_index(np.argmin(scores),scores.shape)
    d0,g0=phase[ii],phase[jj]
    grid_min=float(scores[ii,jj])
    step=2*math.pi/grid
    for _ in range(refinements):
        shifts=np.linspace(-step,step,17)
        patch=objective(d0+shifts[:,None],g0+shifts[None,:])
        ii,jj=np.unravel_index(np.argmin(patch),patch.shape)
        d0,g0=d0+shifts[ii],g0+shifts[jj]
        step/=8
    # Two components of the reduced distance are 1-Lipschitz in T.
    # Varying the output unit z adds at most ||T|| times the change in z
    # to each component (projector norm is <= distance between unit vectors).
    lip_delta=math.sqrt(2)*(h*c*c2+np.linalg.norm(C,2)*t*c2/r)/s
    lip_gamma=math.sqrt(2)*h*math.sqrt(norm_Cu2)/s
    # Cover the torus with cells and split only cells that still contain
    # a possible improvement. Every retained lower bound covers the
    # entire torus, even if the memory cap stops the subdivision early.
    radius=math.pi/grid
    lipschitz=lip_delta+lip_gamma
    upper=float(objective(d0,g0))+1e-12
    certified_lower=max(0.,grid_min-radius*lipschitz-1e-12)
    active=np.nonzero(scores-radius*lipschitz<=upper)
    ds,gs=phase[active[0]],phase[active[1]]
    subdivisions=0
    for _ in range(8):
        if ds.size>65536 or ds.size==0:
            break
        radius/=2
        ds=np.concatenate((ds-radius,ds-radius,ds+radius,ds+radius))
        gs=np.concatenate((gs-radius,gs+radius,gs-radius,gs+radius))
        mids=objective(ds,gs)
        best=int(np.argmin(mids))
        if mids[best]+1e-12<upper:
            d0,g0,upper=ds[best],gs[best],float(mids[best])+1e-12
        certified_lower=max(certified_lower,float(mids.min())-radius*lipschitz-1e-12)
        keep=mids-radius*lipschitz<=upper
        ds,gs=ds[keep],gs[keep]
        subdivisions+=1

    def full_witness(delta,gamma,which):
        ed=np.exp(1j*delta)
        vv=c*F[:,0]+t*ed*F[:,1]
        uu=-t*F[:,0]+c*ed*F[:,1]
        nn=n.copy()
        zz=C@vv/r
        if which=="second":
            w=k*D@u
            common=np.exp(1j*(np.angle(np.vdot(z,w))-np.angle(np.vdot(zz,C@uu))))
            uu=common*uu
            ww=C@uu
        else:
            w=k*D@tag_first
            uu=np.exp(1j*gamma)*uu
            ww=C@(s*nn-h*uu)
            common=np.exp(1j*(np.angle(np.vdot(z,w))-np.angle(np.vdot(zz,ww))))
            nn,uu,ww=common*nn,common*uu,common*ww
        L=basis_from_two(z,w)@basis_from_two(zz,ww).conj().T
        R=np.column_stack((nn,uu,vv))@base_input.conj().T
        E=L@C@R
        deltaE=E-k*D
        a=deltaE@tag_first/s
        f=deltaE@u/s
        err=np.linalg.norm(deltaE-np.outer(a,b.conj())-np.outer(f,b2.conj()),2)
        assert err<1e-12 and abs(abs(np.vdot(n,R@b))-1)<1e-12
        assert np.max(abs(np.linalg.svd(E,compute_uv=False)[::-1]-np.diag(C)))<1e-13
        return {"first":float(np.linalg.norm(a)),"second":float(np.linalg.norm(f)),
                "error":float(err)}

    opt2=full_witness(0.,0.,"second")
    opt1=full_witness(d0,g0,"first")
    assert abs(opt2["second"]-min_second)<1e-13
    assert abs(opt1["first"]-objective(d0,g0))<1e-10
    assert certified_lower<=opt1["first"]+1e-13
    return {"k":k,"unused":TAGS[m],"sin_between_inputs":s,"r":r,
            "admissible_k":[k*c1/r,k*c2/r],"second_exact_min":min_second,
            "first_min_bracket":[certified_lower,opt1["first"]+1e-12],
            "phase_subdivisions":subdivisions,"first_witness_phases_rad":[float(d0),float(g0)],
            "second_minimizing_witness":opt2,"first_witness":opt1}


def full_box(m,a,osc,atm):
    lo,hi=m[:,0]-2*m[:,1],m[:,0]+2*m[:,1]
    amin,amax=a[:,0]-2*a[:,1],a[:,0]+2*a[:,1]
    s12,s13,s23=np.sin(amin[:3]); S12,S13,S23=np.sin(amax[:3])
    c12,c13,c23=np.cos(amax[:3])
    # Each row lists positive |V_jd| lower and |V_js| upper, then
    # |V_jb| lower and upper. All CP phases are included.
    row=np.array([[c12*c13,S12,s13,S13],
                  [s12*c23-S23*S13,1.,s23*c13,S23],
                  [s12*s23-S13,S23+S12*S13,c23*c13,1.]])
    ratio=lo[2]/np.sqrt(1+(hi[1]/lo[2]*row[:,1]/row[:,2])**2+
                          (hi[0]/lo[2]/row[:,2])**2)
    slopes=np.hypot(hi[0]/lo[1]*row[:,1]/row[:,0],hi[0]/lo[2]*row[:,3]/row[:,0])
    tilt=np.arctan(slopes)
    pb=(hi[0]/lo[2]*row[:,3]/row[:,0])**2
    s2lo=lo[2]*np.sqrt(1-pb)
    out={}
    for order in ("NO","IO"):
        n=osc[order]
        if order=="NO":
            tau_max=math.sqrt(n["s12"][2]*atm[order][2])+math.sqrt((1-n["s12"][1])*(1-atm[order][1])*n["s13"][2])
            ne=math.sqrt((1-n["s12"][2])*(1-n["s13"][2]))
        else:
            tau_max=math.sqrt((1-atm[order][1])*(1-n["s13"][1]))
            ne=math.sqrt(n["s13"][1])
        c2lo=lo[5]*math.sqrt(1-tau_max*tau_max)
        all_k=[hi[4]/min(ratio),c2lo/hi[2]]
        assert all_k[0]<1 and all_k[1]>1.02
        second={}
        for j,l in itertools.permutations(range(3),2):
            sin_max=math.sin(min(math.pi/2,tilt[j]+tilt[l]))
            free=max(0.,ne*lo[4]-hi[1]*hi[5]/s2lo[j])/sin_max
            fixed=(max(0.,ne*lo[4]-hi[1])+max(0.,c2lo-hi[2]))/sin_max
            cos_min=math.sqrt(max(0.,1-sin_max**2))
            other_e=max(0.,cos_min*ne-sin_max*math.sqrt(1-ne**2))
            other_tau_plane=max(0.,cos_min*math.sqrt(1-tau_max**2)-sin_max*tau_max)
            first=(max(0.,other_e*lo[4]-hi[1])+max(0.,other_tau_plane*lo[5]-hi[2]))/sin_max
            second[TAGS[j]+TAGS[l]]={"sin_max":sin_max,"free_k_floor":free,
                                      "k1_floor":fixed,"first_k1_floor":first}
        out[order]={"all_k_sufficient_band":all_k,"r_unscaled_lower_by_unused":ratio.tolist(),
                    "c2_lower":c2lo,"second_source":second}
    return out


def compute():
    cfg=read("inputs.json")
    m=np.array(read(cfg["spectra"])["scales"]["1e16"])
    ang=np.array(read(cfg["angles"])["high_scale"]["angles_theta12_theta13_theta23_delta"])
    osc=read(cfg["oscillations"]); atm=read(cfg["atmospheric"])["theta23_sin2"]
    out={"box":full_box(m,ang,osc,atm),"central":{}}
    for order in ("NO","IO"):
        _,D,C,n=central_inputs(cfg,m,ang,osc,atm,order)
        out["central"][order]={}
        for j,l in itertools.permutations(range(3),2):
            name=TAGS[j]+TAGS[l]
            out["central"][order][name]=[single_point(D,C,n,j,l,k,cfg["phase_grid"],
                                         cfg["phase_refinements"]) for k in cfg["central_gains"]]
    return out


def close(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(close(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    if isinstance(a,str):
        return a==b
    return math.isclose(a,b,rel_tol=3e-8,abs_tol=2e-11)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=compute()
    path=HERE/"results.json"
    if args.check:
        assert close(output,read("results.json"))
        print("1605: two orderings, six ordered source pairs, sharp distances and reconstructions passed")
    else:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(output,stream,indent=2,allow_nan=False)
            stream.write("\n")
