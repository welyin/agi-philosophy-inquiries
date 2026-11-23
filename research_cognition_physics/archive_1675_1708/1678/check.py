"""Joint charged phases and registered neutral interfaces with one frozen flavon direction."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
OLD=runpy.run_path(str(HERE/CFG["old_charged"]))
NU=runpy.run_path(str(HERE/CFG["neutrino"]))
E0=OLD["qcenter"]
L0=OLD["lcenter"]
phi=CFG["phi_over_pi"]*np.pi


def inputs_at_gamma(g,q=E0):
    u,v,w=q[3:6]
    k=np.cos(u)*np.tan(w)*np.sin(v)/np.sin(u)
    gamma=np.radians(g)
    delta=gamma+np.arcsin(k*np.sin(gamma))
    return OLD["LOW"]["texture"]([*q,delta]),delta


def parameters(g,q=E0,l=L0):
    edges,delta=inputs_at_gamma(g,q)
    a,b,c,e,f,zr,zi=edges
    z=complex(zr,zi)
    K=a*a+b*b+9*c*c
    L=e*e+f*f+abs(z)**2
    h=abs(a*b*f)
    obs=np.array(l)**2
    P=sum(obs[i]*obs[j] for i,j in ((0,1),(1,2),(0,2)))/sum(obs)**2
    D=np.prod(l)/sum(obs)**1.5
    alpha=(D*K*np.sqrt(L)/h)**2
    roots=sorted(r.real for r in np.roots([1,-2,1,-alpha]) if abs(r.imag)<1e-9 and 0<r.real<1)
    assert len(roots)==2
    out=[]
    for w in roots:
        s=np.sqrt(w*K/((1-w)*L))
        T=K+L*s*s
        base=(a*b)**2+s*s*(a*a*(abs(z)**2+f*f)+b*b*(e*e+f*f)+(3*c*f)**2)+s**4*(e*abs(z))**2
        lever=6*c*f*e*abs(z)*s**3
        cos=(base-P*T*T)/lever
        out.append((s,np.sqrt(sum(obs)/T),cos))
    return edges,delta,out


def residual(g,k,branch=1,q=E0,l=L0):
    edges,delta,states=parameters(g,q,l)
    z=complex(*edges[-2:])
    return np.cos(np.angle(z)+k*phi)-states[branch][2]


OV=OLD["IV"]
I=OV["I"]
ib=OV["bounds"]
charged_base=OLD["interval_domain"](E0)
D,P=OV["spectral_target"]([[v,0] for v in L0],0)


def phase_pair(k):
    # All twelfth-root sines/cosines from exact rationals and interval sqrt(3).
    one,half,r=I.point(1),I.point(".5"),I.point(3).sqrt()/2
    zero=I.point(0)
    cos=[one,r,half,zero,-half,-r,-one,-r,-half,zero,half,r]
    return cos[k%12],cos[(k+3)%12]


def interval_residuals(x,sg,keys,base_quarks=charged_base):
    co,u0,_,_,_,aux=OV["coefficients"](base_quarks+[x],(1,1,3))
    h,L,base,wf,er,ei,K=aux
    ei=sg*ei
    levels=list(OV["level_set"](D,co))
    out={}
    for k,b in keys:
        c,s=phase_pair(k)
        real,imag=c*er-s*ei,s*er+c*ei
        w=levels[b]
        rat=D*w.sqrt()/(h*L.sqrt())
        v=w/L
        residual=(u0*(1-w).square()+base*rat+(imag*v).square()
                  +((wf*rat).sqrt()-real*v).square()-P)
        out[(k,b)]=residual
    return out


def covers():
    # Enlargements of both direct two-width arcs, |delta-gamma|<0.001 radians.
    # More than 0.02 of slack in cos at either endpoint.
    allkeys=[(k,b) for k in range(12) for b in (0,1)]
    enclosures={}
    runs=[]
    for sg,domain in [(1,ib(".35",".54")),(-1,ib("-.54","-.35"))]:
        todo=[(domain,allkeys)]
        count=0
        while todo:
            x,keys=todo.pop()
            count+=1
            rr=interval_residuals(x,sg,keys)
            keys=[key for key,val in rr.items() if val.lo<=0<=val.hi]
            if not keys:
                continue
            if x.hi-x.lo>OV["iv"].SCALE//100000000:
                todo.extend((xx,keys) for xx in x.parts(2))
            else:
                for k,b in keys:
                    bounds=enclosures.setdefault((sg,k,b),[1.,-1.])
                    bounds[0]=min(bounds[0],x.values()[0])
                    bounds[1]=max(bounds[1],x.values()[1])
        runs.append(count)
    return enclosures,runs


def active_columns(xi, sign, k, branch):
    """Complex entries enclosed in the same row/rephasing convention as 1657."""
    d,sm,p,u,t=charged_base
    cp,cu,ct=[(1-y.square()).sqrt() for y in (p,u,t)]
    sine=sign*(1-xi.square()).sqrt()
    base=d.square()*cp.square()+sm.square()*p.square()
    mu=(sm.square()-d.square())*p*cp
    a2=base*cu.square()+u.square()
    h12=(cu*(ct*mu+t*u*(1-base)*xi),cu*t*u*(1-base)*sine)
    h13=(cu*(-t*mu+ct*u*(1-base)*xi),cu*ct*u*(1-base)*sine)
    c=(OV["abs2"](h12)/a2).sqrt()
    e=(OV["abs2"](h13)/a2).sqrt()
    v2=[(-p*ct-cp*t*u*xi,-cp*t*u*sine),(cp*ct-p*t*u*xi,-p*t*u*sine),(t*cu,I.point(0))]
    minor=sum(q*OV["abs2"](z) for q,z in zip((sm.square(),d.square(),sm.square()*d.square()),v2))
    f=(minor/a2).sqrt()
    b=-d*sm/minor.sqrt()
    co,_,_,_,_,aux=OV["coefficients"](charged_base+[xi],(1,1,3))
    h,L,_,_,er,ei,K=aux
    w=list(OV["level_set"](D,co))[branch]
    abs_t=(w*K/((1-w)*L)).sqrt()
    re,im=phase_pair(k)
    re,im=re*abs_t,im*abs_t
    z=(er/e,sign*ei/e)
    return {(0,1):(a2.sqrt(),I.point(0)),(1,0):(b,I.point(0)),(1,1):(3*c,I.point(0)),
            (1,2):OV["cmul"]((re,im),z),(2,1):(re*e,im*e),(2,2):(re*f,im*f)}


def mass_basis(g,k,b,fullx):
    edges,delta,states=parameters(g)
    mag,scale,_=states[b]
    M=OLD["LOW"]["mtx"](edges,[1,1,3],mag*np.exp(1j*k*phi))
    y=M.T
    eig,charged=np.linalg.eigh(y@y.conj().T)
    np.testing.assert_allclose(np.sqrt(eig)*scale,L0,rtol=1e-8)
    delta2=0.
    for ij,(real,imag) in active_columns(ib(str(fullx[0]-1e-12),str(fullx[1]+1e-12)),1,k,b).items():
        for num,component in ((M[ij].real,real),(M[ij].imag,imag)):
            low,high=component.values()
            assert low-1e-11 <= num <= high+1e-11
            delta2+=max(abs(low-num),abs(high-num))**2
    dm=np.sqrt(delta2)+1e-10
    dh=(2*np.linalg.norm(M,2)+dm)*dm
    gap=float(np.diff(eig).min())
    # Davis--Kahan per spectral line of Y_e Y_e^dagger.
    ep=dh/(gap-dh)+1e-8
    return charged,mag,scale,ep


def neutral(charged,ep,ph,sg):
    # Positive r on each of the two CSD(3) determinant branches.
    rbox=NU["gap_box"](2)
    g=np.linspace(*rbox,4097)
    R=(g[:-1]+g[1:])/2
    h=(g[-1]-g[0])/(2*len(R))
    phs=NU["PHASE"]
    co=phs["coefficients"](charged,ph)
    rsmall=phs["root"](R)
    windows=NU["read"](NU["INPUT"]["angle_windows"])["joint_windows"]
    reports={}
    for tag,r in [("small",rsmall),("large",4/(121*rsmall))]:
        root_phase=np.exp(2j*np.pi*sg/3)
        pe=phs["projector"](co,root_phase,R,r,0,0).real
        pm=phs["projector"](co,root_phase,R,r,1,1).real
        p1=abs(co[2][0])**2
        # Direct norm bound using r' and the minimum R in the whole interval.
        rmin, rmax=(.08,.13) if tag=="small" else (.25,.42)
        assert rmin < min(r) <= max(r) < rmax
        Rmin=rbox[0]
        coeff=6*(Rmin+1/Rmin)+16
        dr=6*rmax/np.sqrt(coeff*coeff-1936)
        slope=(4*dr/rmin**2+121*dr+6/Rmin**2+6*(1+1/Rmin**2))/(6*(Rmin-1/Rmin))
        err=ep+slope*h+2e-8
        pl,pu=pe-err,pe+err
        assert np.min(pl)>0 and np.max(pu)<.5
        lows=np.array([4*pl*(1-pl), 1-(p1+ep)/(1-pu), (pm-err)/(1-pl)])
        highs=np.array([4*pu*(1-pu), 1-(p1-ep)/(1-pl), (pm+err)/(1-pu)])
        win=np.array(windows)
        ok=((lows<=win[:,1:])&(highs>=win[:,:1])).all(axis=0)
        reports[tag]={"R_grid_cells_passing_outer":int(sum(ok)),
                      "three_angles_outer":[[lo.min(),hi.max()] for lo,hi in zip(lows,highs)]}
        if np.any(ok):
            actual_r=float(np.median(R[ok]))
            mix=NU["neutrino_mix"](charged,ph,sg,actual_r) if tag=="small" else None
            assert mix is not None
            angles=NU["OLD"]["angles"](mix)
            assert np.all((angles>win[:,0]+5*ep)&(angles<win[:,1]-5*ep))
            solar=NU["window"](NU["INPUT"]["juno_NO_dm21"],2)
            heavy=NU["window"](NU["INPUT"]["juno_NO_dm31"],2)
            ds=(max(solar[0],heavy[0]/actual_r**2)+min(solar[1],heavy[1]/actual_r**2))/2
            assert solar[0]<ds<solar[1] and heavy[0]<ds*actual_r**2<heavy[1]
            reports[tag]["common_point"]={"R":actual_r,"d21":ds,"d31":ds*actual_r**2,"three_angles":angles.tolist()}
    return reports


def main():
    allx,nodes=covers()
    assert set(allx)=={(1,1,1),(1,2,1)}
    out={}
    for (sg,k,b),xs in allx.items():
        g=OLD["zero"](lambda gg:residual(gg,k,b),58.8,67.8)
        # Endpoint interval signs prove existence independently of the point root finder.
        witness_x=np.cos(inputs_at_gamma(g)[1])
        signed=[]
        for xx in (witness_x-2e-5,witness_x+2e-5):
            signed.append(interval_residuals(I.point(xx),sg,[(k,b)])[(k,b)].values())
        assert signed[0][1]<0<signed[1][0] or signed[1][1]<0<signed[0][0]
        charged,mag,scale,ep=mass_basis(g,k,b,xs)
        phases={"third_3_2":([2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))],1),
                "second_2_minus1":([10*phi+np.angle(2-np.exp(4j*phi)),7*phi],-1),
                "second_original_11_5deg":([10*phi+np.radians(11.5),7*phi],-1)}
        out[str(k)]={"x":xs,"gamma_deg":g,"g_interval_deg":
                     [OLD["gamma_deg"](xs[1],E0),OLD["gamma_deg"](xs[0],E0)],
                     "abs_t":mag,"k_times_yb":scale,"signed_residuals":signed,
                     "charged_projector_error":ep,
                     "oscillation":{n:neutral(charged,ep,ph,s) for n,(ph,s) in phases.items()}}
    free={}
    for g in [63.2,-116.8]:
        edges,_,states=parameters(g)
        z=complex(*edges[-2:])
        c=states[1][2]
        free[str(g)]={"cos_required_small_large":[s[2] for s in states],
                      "large_phase_angles_deg":[float(np.degrees(-np.angle(z)+sign*np.arccos(c))) for sign in (-1,1)]}
    q=E0.copy()
    q[0]-=.85*OLD["qrows"][0,1]
    g=OLD["zero"](lambda gg:residual(gg,1,1,q,L0),61.,65.5)
    edges,delta,states=parameters(g,q,L0)
    x=np.cos(delta)
    shifted_intervals=OLD["interval_domain"](q)
    signs=[interval_residuals(I.point(x+step),1,[(1,1)],shifted_intervals)[(1,1)].values()
           for step in (-2e-5,2e-5)]
    assert signs[0][1]<0<signs[1][0] or signs[1][1]<0<signs[0][0]
    s,scale,_=states[1]
    y=scale*OLD["LOW"]["mtx"](edges,[1,1,3],s*np.exp(1j*phi)).T
    eig,uc=np.linalg.eigh(y@y.conj().T)
    np.testing.assert_allclose(np.sqrt(eig),L0,rtol=1e-10)
    R,ds=5.83,7.388e-5
    p3=[2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))]
    angles=NU["OLD"]["angles"](NU["neutrino_mix"](uc,p3,1,R))
    windows=np.array(NU["read"](NU["INPUT"]["angle_windows"])["joint_windows"])
    assert np.all((angles>windows[:,0]) & (angles<windows[:,1]))
    assert NU["window"](NU["INPUT"]["juno_NO_dm31"],1)[0]<R*R*ds<NU["window"](NU["INPUT"]["juno_NO_dm31"],1)[1]
    return {"cover_nodes_positive_negative_orientation":nodes,"all_discrete_solutions":out,"continuous_phase":free,
            "one_shift_full_witness":{"down_yukawa_pull":-.85,"other_eight_CP_even_pulls":0,
                                      "k":1,"gamma_deg":g,"abs_t":s,"k_times_yb":scale,"R":R,
                                      "d21":ds,"d31":R*R*ds,"three_angles":angles.tolist(),
                                      "signed_bracket":signs}}


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    mode=cli.add_mutually_exclusive_group()
    mode.add_argument("--write",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=cli.parse_args()
    result=OLD["tidy"](main())
    path=HERE/"results.json"
    if args.write:
        path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    else:
        assert json.loads(path.read_bytes())==result
    print(json.dumps(result,indent=2))
