"""Neutral oscillation holdouts with the tagged charged CP basis kept fixed."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy
import numpy as np
from numpy.polynomial import Polynomial

HERE=Path(__file__).resolve().parent
def read(f): return json.loads((HERE/f).read_bytes())
cfg=read("inputs.json")
sim=runpy.run_path(str(HERE/cfg["parent"]))
 
def bounds(row):
    x,d,u=row
    return np.array([x-d,x+u])

def ratio_window():
    f=read(cfg["clocks"])
    s,d3,d2=map(bounds,(f["juno_NO_dm21"],f["juno_NO_dm31"],f["accelerator_NO_dm32"]))
    # Vertices of the joint frequency rectangle cut by the d31-d21 strip.
    pts=list(itertools.product(s,d3))+[(x,x+z) for x,z in itertools.product(s,d2)]
    pts += [(y-z,y) for y,z in itertools.product(d3,d2)]
    ratios=[math.sqrt(y/x) for x,y in pts if s[0]-1e-16<=x<=s[1]+1e-16 and
            d3[0]-1e-16<=y<=d3[1]+1e-16 and d2[0]-1e-16<=y-x<=d2[1]+1e-16]
    center=math.sqrt(f["juno_NO_dm31"][0]/f["juno_NO_dm21"][0])
    return min(ratios),center,max(ratios)

def real_roots(p):
    c=np.array(p.coef,float)
    p=Polynomial(c/max(abs(c)))
    roots=p.roots()
    return sorted(z.real for z in roots if abs(z.imag)<1e-8)

def features(Q,O,V,x):
    M=O+x*V
    u,s,_=np.linalg.svd(M)
    mixing=Q@u[:,::-1].conj()
    p=abs(mixing)**2
    z=p[0,1]
    n=np.linalg.svd(Q.conj().T)[2][-1].conj()
    pm=np.column_stack([n,mixing])
    J=(pm[0,0]*pm[1,1]*np.conj(pm[0,1]*pm[1,0])).imag
    sol,atm=p[0,0]/(1-z),p[1,1]/(1-z)
    cos=(abs(n[1])**2-sol*(1-atm)-(1-sol)*atm*z)/(2*math.sqrt(sol*(1-sol)*atm*(1-atm)*z))
    sin=J/((1-z)*math.sqrt(sol*(1-sol)*atm*(1-atm)*z))
    return s[0]/s[1],np.array([4*z*(1-z),sol,atm]),np.arctan2(sin,cos),J

def compact_certificate(n,Q,O,V,rmin,rmax,angwin,max_depth=24,only_reactor=False,necessary_J_upper=None):
    # All real coupling strengths, including both limiting unit dyads.
    edges=np.linspace(-math.pi/2,math.pi/2,17)
    L,R=edges[:-1],edges[1:]
    zmin,zmax=(1-np.sqrt(1-angwin[0]))/2
    pe=abs(n[0])**2
    if not only_reactor:
        zmin=max(zmin,1-pe/(1-angwin[1,1]))
        zmax=min(zmax,1-pe/(1-angwin[1,0]))
    visited=0
    for depth in range(max_depth):
        mid=(L+R)/2
        err=math.sqrt(2)*(R-L)/2+2e-12
        mat=np.cos(mid)[:,None,None]*O+np.sin(mid)[:,None,None]*V
        left,s,_=np.linalg.svd(mat)
        visited+=len(L)
        viable=(s[:,0]+err>=rmin*(s[:,1]-err)) & (s[:,0]-err<=rmax*(s[:,1]+err))
        eta=(2*s[:,0]+err)*err
        gap=s[:,0]**2-s[:,1]**2
        sine=np.minimum(1,eta/np.maximum(eta,gap-eta))
        amp=abs(Q@left.conj())
        def envelope(row):
            a,b=amp[:,row,0],amp[:,row,1]
            return np.maximum(0,a*np.sqrt(np.maximum(0,1-sine*sine))-b*sine)**2,np.minimum(1,(a+b*sine)**2)
        lo,hi=envelope(0)
        ml,mh=envelope(1)
        viable &= (hi>=zmin-1e-10)&(lo<=zmax+1e-10)
        if not only_reactor:
            viable &= (mh>=angwin[2,0]*(1-np.minimum(zmax,hi))-1e-10)
            viable &= (ml<=angwin[2,1]*(1-np.maximum(zmin,lo))+1e-10)
        if necessary_J_upper is not None:
            H=mat.conj().transpose(0,2,1)@mat
            j_num=-(n[0]*n[1].conj() * (Q[1]@H@Q[0].conj())).imag
            jerr=abs(n[0]*n[1])*np.linalg.norm(Q[1])*np.linalg.norm(Q[0])*eta
            # gap' >= gap-2*eta, and the required J upper is negative.
            viable &= j_num-jerr <= necessary_J_upper*np.maximum(0,gap-2*eta)+1e-12
        if not np.any(viable):
            return {"cells":visited,"depth":depth+1}
        L,R=L[viable],R[viable]
        mid=(L+R)/2
        L,R=np.r_[L,mid],np.r_[mid,R]
    raise AssertionError("Unresolved global compact interval")

def fc_read(theta,delta):
    p=read(cfg["accelerator_FC"])["official_coordinates_rad"]
    inside=runpy.run_path(str(HERE/"../1705/independent_check.py"))["interior"]
    return [any(inside(p[key],theta,delta) for key in keys) for keys in
            [["conf68_0"],["conf90_0","conf90_1"]]]

def shape(n,Q,W,A,b,k,psi):
    z=np.exp(1j*np.pi/3*np.array([0,*k]))
    C=W.T@np.diag(1/(z*A))@W
    Y=C-np.outer(C@n,C@n)/(n@C@n)
    unit=np.linalg.norm(Y,2)
    O=Q.T@Y@Q/unit
    hb=Q.T@b
    weight=np.linalg.norm(hb)**2
    V=np.exp(1j*psi)*np.outer(hb,hb)/weight
    F=Polynomial([np.vdot(O,O).real,2*np.vdot(O,V).real,np.vdot(V,V).real])
    det0=np.linalg.det(O)
    det1=O[0,0]*V[1,1]+O[1,1]*V[0,0]-2*O[0,1]*V[0,1]
    cross=Polynomial([det0,det1])
    detabs=Polynomial((cross*Polynomial(cross.coef.conj())).coef.real)
    def roots_at(r):
        p=F*F-(r+1/r)**2*detabs
        roots=real_roots(p)
        for x in roots:
            assert abs(features(Q,O,V,x)[0]-r)<1e-7
        return roots
    lo,center,hi=ratio_window()
    endpoints=sorted(roots_at(lo)+roots_at(hi))
    intervals=[(a,b) for a,b in zip(endpoints,endpoints[1:])
               if lo<=features(Q,O,V,(a+b)/2)[0]<=hi]
    angles=np.array(read(cfg["angles"])["joint_windows"])
    def row(alpha):
        q=Q[alpha]
        H=[O.conj().T@O,O.conj().T@V+V.conj().T@O,V.conj().T@V]
        return Polynomial([2*np.vdot(q.conj(),v@q.conj()).real for v in H])-np.linalg.norm(q)**2*F
    e,m=row(0),row(1)
    pe,pm=np.sum(abs(Q[:2])**2,axis=1)
    G=F*F-4*detabs
    # Each threshold crossing is among the real roots of a quartic.
    def cut(poly,c):
        return real_roots(poly*poly-c*c*G)
    zr=(1-np.sqrt(1-angles[0]))/2
    zsolar=1-abs(n[0])**2/(1-angles[1])
    crossings=[x for z in [*zr,*zsolar] for x in cut(e,2*z-pe)]
    crossings += [x for w in angles[2] for x in cut(m+w*e,2*w-pm-w*pe)]
    good=[]
    for a,b in intervals:
        sep=[a]+sorted(x for x in crossings if a<x<b)+[b]
        for l,r in zip(sep,sep[1:]):
            test=(l+r)/2
            ratio,ang,delta,J=features(Q,O,V,test)
            if np.all((angles[:,0]<=ang)&(ang<=angles[:,1])):
                good.append({"interval":[l,r],"mid":test,"angles":ang,"ratio":ratio,"delta_deg":np.degrees(delta),
                             "FC68_90_at_mid":fc_read(ang[2],delta)})
    at_center=[]
    for x in roots_at(center):
        _,ang,d,J=features(Q,O,V,x)
        three=Q.conj()@(O+x*V)@Q.conj().T
        masses=np.linalg.svd(three,compute_uv=False)
        solar=read(cfg["clocks"])["juno_NO_dm21"][0]
        at_center.append(dict(x=x,angles=ang,delta_deg=np.degrees(d),J=J,mbb_eV=abs(three[0,0])*math.sqrt(solar)/masses[1],
                              FC68_90=fc_read(ang[2],d) if good else None))
    rejection=None if good else compact_certificate(n,Q,O,V,lo,hi,angles)
    reactor_cert=None
    if k==(0,0) and psi==0:
        reactor_cert=compact_certificate(n,Q,O,V,lo,hi,angles,only_reactor=True)
    acc={}
    if good:
        pol=read(cfg["accelerator_FC"])["official_coordinates_rad"]
        zl,zh=(1-np.sqrt(1-angles[0]))/2
        # p_mu1 enforces cos(delta)>0 at all relevant points.
        assert abs(n[1])**2-angles[1,1]*(1-angles[2,0])-(1-angles[1,0])*angles[2,1]*zh>0
        dmin=(1-zh)*math.sqrt(zl*min(angles[1]*(1-angles[1]))*min(angles[2]*(1-angles[2])))
        assert min(np.array(pol["conf90_1"])[:,1])>np.pi/2
        for ix,label in enumerate(["68","90"]):
            if any(w["FC68_90_at_mid"][ix] for w in good):
                acc[label]="witness"
            else:
                lim=max(np.array(pol["conf"+label+"_0"])[:,1])
                acc[label]=compact_certificate(n,Q,O,V,lo,hi,angles,necessary_J_upper=dmin*math.sin(lim))
    return dict(k=k,psi=psi,rho_per_x=unit/weight,frequency_intervals=intervals,accepted=good,
                two_central_clock_roots=at_center,full_rejection_certificate=rejection,
                reactor_only_certificate=reactor_cert,FC_if_angles_pass=acc)

def phase_data(raw,A):
    n=np.mod(raw,3)
    signs=np.where((4*n-raw)%6==0,1,-1)
    coeff=A*np.array([1.,1.,4.])*32.**n
    return {"raw_phase_indices":raw.tolist(),"neutral_monomial_powers":n.tolist(),
            "real_signs":signs.tolist(),"absolute_up_Wilsons":coeff.tolist()}

def cp_reconstruction(n,W,A,low_b,root,step):
    transport=sim["FLOW"]
    flow,rep=transport["charged_background"](step)
    D=transport["unpack"](rep["down_texture_10PeV"])
    F=np.linalg.solve(D,W)
    flavor=read(cfg["scalar"])
    r=flavor["scalar_magnitude_over_cutoff"]
    phi=np.pi*flavor["scalar_angle_over_pi"]
    s=r*np.exp(1j*phi)
    X=flavor["real_X_over_cutoff"]
    I=X*s**4
    omega=np.exp(4j*phi)
    S=D[1,2].imag/omega.imag
    R=D[1,2].real-S*omega.real
    poles=np.array(read(cfg["quark_monomials"])["powers_A_B_C_E_F_R_D"])
    aa,bb,c,e,f,short,long=np.array([D[0,1].real,D[1,0].real,D[1,1].real,D[2,1].real,D[2,2].real,R,S])*np.exp(1j*poles*phi)
    rawD=np.array([[0,aa,0],[bb,c,short+long],[0,e,f]])
    P10=np.diag(np.exp(1j*phi*np.array([0,0,1])))
    P5=np.diag(np.exp(1j*phi*np.array([2,0,1])))
    np.testing.assert_allclose(rawD,P10@D@P5,atol=1e-15)
    FR=P5.conj()@F
    vr=FR@n
    br=P5@np.linalg.solve(F.T,low_b)
    heavy=flavor["cutoff_TeV"]*np.array([[0,s**3],[-s**3,1.]])
    bare=np.diag(A)@P10@P10
    revised=bare@np.diag([I/abs(I),I**2/abs(I)**2,1])
    np.testing.assert_allclose(revised@revised.conj().T,bare@bare.conj().T,atol=1e-16)
    top_phases=[np.angle(np.linalg.det(m)*np.linalg.det(rawD)*np.linalg.det(heavy)) for m in (bare,revised)]
    assert max(abs(np.array(top_phases)))<1e-12
    C=rawD.T@np.linalg.solve(revised,rawD)
    triplet=runpy.run_path(str(HERE/"independent_check.py"))
    pole=triplet["MASS"]
    scalar,leg=triplet["below_both_clusters"](flow,rep,pole)
    ew=read("../1631/inputs.json")["masses_GeV"]["W"]
    tW=(ew/pole)**2
    kin=math.sqrt(1-4*tW)*(1-4*tW+12*tW*tW)
    out=[]
    for row in root["two_central_clock_roots"]:
        rho=-row["x"]*root["rho_per_x"]
        low=FR.T@(C-np.outer(C@vr,C@vr)/(vr@C@vr)+rho*np.outer(br,br))@FR
        direct=sim["weak_observables"](low)
        np.testing.assert_allclose(direct["angles"],row["angles"],atol=1e-10)
        np.testing.assert_allclose(direct["J"],row["J"],atol=1e-10)
        scale=math.sqrt(read(cfg["clocks"])["juno_NO_dm21"][0])/np.linalg.svd(low,compute_uv=False)[1]
        m_b=abs(rho)*scalar*np.linalg.norm(leg*low_b)**2*scale
        product=(pole**2*m_b*1e-9/(8*np.pi*triplet["VEV"]**2))**2*kin
        out.append({"raw_real_b_coefficient":rho,"angles":row["angles"],"J":row["J"],
                    "delta_deg":row["delta_deg"],"mbb_eV_at_central_solar":row["mbb_eV"],
                    "conditional_b_mass_eV_at_2TeV":m_b,"conditional_triplet_width_product_GeV2":product,
                    "SK_T2K_FC68_90":row["FC68_90"]})
    return {"scalar_neutral_modulus":abs(I),"original_and_revised_full_tree_strong_angle":top_phases,
            "heavy_color_masses_TeV":np.linalg.svd(heavy,compute_uv=False),
            "original_absolute_up_Wilsons":A*np.array([1.,1.,4.]),
            "selected_absolute_up_Wilsons":A*np.array([32.,1024.,4.]),
            "up_charm_top_log_radial_slope_original":[0,0,2],
            "up_charm_top_log_radial_slope_selected":[4,8,2],
            "full_colored_determinant_angular_exponents":[12,24],"direct_forward_both_roots":out}

def calculate(step):
    n,Q,_,_,_,(W,A,b,_),_=sim["background"](step)
    rows=[]
    all_max=A[2]*4
    for k in itertools.product(range(6),repeat=2):
        for l in range(3):
            r=shape(n,Q,W,A,b[:,1],k,l*np.pi/3)
            strong=[]
            for first in (l,l+3):
                raw=(np.array([0,*k])+first)%6
                if sum(raw)%6==0:
                    strong.append(phase_data(raw,A))
            r["strong_compatible_raw_choices"]=strong
            rows.append(r)
    good=[r for r in rows if r["accepted"]]
    accepted=[r for r in good if r["FC_if_angles_pass"]["90"]=="witness"]
    assert len(good)==4 and len(accepted)==1
    selected=accepted[0]
    assert selected["k"]==(4,2) and abs(selected["psi"]-np.pi/3)<1e-13
    assert selected["strong_compatible_raw_choices"][0]["raw_phase_indices"]==[4,2,0]
    cert=[r["full_rejection_certificate"] for r in rows if r["full_rejection_certificate"]]
    base=rows[0]
    costs=[s for r in rows for s in r["strong_compatible_raw_choices"]]
    budget=sum(bool(max(s["absolute_up_Wilsons"])<=all_max+1e-12) for s in costs)
    assert len(costs)==36 and budget==12
    cp=cp_reconstruction(n,W,A,b[:,1],selected,step)
    def mixed(tau):
        theta=np.angle(1+tau*np.exp(-2j*np.pi/3))
        r=shape(n,Q,W,A,b[:,1],(-2*theta*3/np.pi,-theta*3/np.pi),theta)
        return {"tau":tau,"theta_deg":np.degrees(theta),"angle_witness":bool(r["accepted"]),
                "center_roots":r["two_central_clock_roots"],"full_angle_rejection":r["full_rejection_certificate"]}
    return {"frequency_ratio_min_center_max":ratio_window(),"registered_original":base,
            "equivalent_neutral_candidates":len(rows),"full_angle_failures":len(cert),
            "exclusion_visited_total":sum(c["cells"] for c in cert),"exclusion_max_depth":max(c["depth"] for c in cert),
            "strong_angle_preserving_raw_choices":len(costs),"at_original_max_single_Wilson":budget,
            "three_angles_accepted":good,"original_basis_verification":cp,
            "two_monomial_competition":[mixed(4.),mixed(16.)]}

if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--check",action="store_true")
    checking=args.parse_args().check
    result=calculate(cfg["log_step"])
    # Re-run the selected forward model only; the full torus exclusion is algebraic at fixed inputs.
    n,Q,_,_,_,(W,A,b,_),_=sim["background"](2*cfg["log_step"])
    alt=shape(n,Q,W,A,b[:,1],(4,2),np.pi/3)
    err=max(abs(a["delta_deg"]-bb["delta_deg"]) for a,bb in zip(alt["two_central_clock_roots"],result["original_basis_verification"]["direct_forward_both_roots"]))
    assert err<1e-6
    result["maximum_halfstep_root_delta_degree_change"]=err
    result=sim["FLOW"]["tidy"](result)
    dest=HERE/"results.json"
    if checking:
        assert result==json.loads(dest.read_bytes())
        print("1718 fixed CP basis, complete discrete search and raw strong-angle forward checks passed")
    else:
        dest.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
