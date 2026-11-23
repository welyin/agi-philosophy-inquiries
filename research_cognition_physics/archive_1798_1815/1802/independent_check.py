"""Triplet matching at its physical threshold and the gauge cost of an isolated triplet pair."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
LOOP=16*math.pi**2
CAP=4*math.pi


def read(name):
    return json.loads((HERE/name).read_bytes())


def low_universal(run):
    """Restore the common coefficient that cancels from the R1715 direction fit."""
    f=run.flow
    rep=run.rep
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(f["central"]("up"))
    y[1,:,:3]=f["low_ckm"](rep["gamma_MZ_deg"],f["angles"]("direct"))@np.diag(f["central"]("down"))
    lepton=rep["charged_left_MZ"]
    y[2,:,:3]=(np.array(lepton["real"])+1j*np.array(lepton["imag"]))@np.diag(f["central"]("charged_leptons"))
    quart=f["physics"]["quartic"]
    z=np.r_[y.ravel(),quart["lambda_paper"][0]/quart["lambda_paper_factor"],0j]
    for low,susy,dt,series in f["segments"](.1):
        y=z[:48].reshape(4,3,4)
        if low==f["MS"]:
            tan=f["TAN"]
            y/=np.array([tan,1,1,tan])[:,None,None]/math.hypot(tan,1)
        for j,pole in enumerate(f["fields"]["neutral_poles_GeV"]):
            if low==pole:
                y[3,:,j]=f["physics"]["Y0"][:,j]
        def rhs(z,g):
            y=z[:48].reshape(4,3,4)
            lam=z[-2].real
            grams=y@y.conj().transpose(0,2,1)
            tr=grams.trace(axis1=1,axis2=2).real
            if susy:
                dl=0
                scalar=6*tr[0]+2*tr[3]-6*g[0]/5-6*g[1]
            else:
                weights=np.array([3,3,1,1])
                total=weights@tr
                tr2=weights@(grams@grams).trace(axis1=1,axis2=2).real
                gy=.6*g[0]
                dl=24*lam*lam+(4*total-9*g[1]-3*gy)*lam-2*tr2+3/8*(2*g[1]**2+(g[1]+gy)**2)
                scalar=4*lam-3*g[1]+2*total
            return np.r_[f["beta"](y,g,susy).ravel(),dl,scalar]
        for i in range(0,len(series)-1,2):
            a=rhs(z,series[i])
            b=rhs(z+a*dt/2,series[i+1])
            c=rhs(z+b*dt/2,series[i+1])
            d=rhs(z+c*dt,series[i+2])
            z+=dt*(a+2*b+2*c+d)/6
    np.testing.assert_allclose(np.linalg.svd(z[:48].reshape(4,3,4)[0],compute_uv=False)[::-1],run.up,rtol=1e-8)
    return math.exp(z[-1].real)


def gauge_window():
    ew=read("../1369/inputs.json")
    cf=read("../1769/inputs.json")
    baseline=read("../1698/results.json")["shared_witness"]
    U=baseline["U_GeV"]
    strong=cf["strong_lattice_5f"]["ALPHA_2026"]
    top=cf["top_MS_matching_GeV"]
    corners=[]
    salvaged={}
    for a,b,c,d in itertools.product((-2,2),repeat=4):
        def shift(key,n):
            g,s=ew[key]
            return 4*math.pi*((g+n*s)**-2-g**-2)
        da1,da2=shift("g1_GUT",a),shift("g2",b)
        da3=1/(strong[0]+c*strong[1])+math.log((top[0]+d*top[1])/ew["mu_ref_GeV"])/(3*math.pi)-4*math.pi/ew["g3"][0]**2
        invari=2*math.pi*(5*da1-12*da2+7*da3)
        ell=-invari/30
        h=(2*math.pi*(da1-da2)+.4*ell)/5.6
        corners.append([ell,U*math.exp(h),U*math.exp(h-ell)])
    x=np.array(corners)
    # A triplet and a heavier sextet only; the (3,2,1/6) pair starts at the physical crossing.
    for label,da3 in (("original",0.),("ALPHA_central",1/strong[0]+math.log(top[0]/ew["mu_ref_GeV"])/(3*math.pi)-4*math.pi/ew["g3"][0]**2)):
        L=math.log(U/1e13)
        h,s=np.linalg.solve([[5.6-.4,3.2],[4+4,-5]],[.4*L,-2*math.pi*da3-4*L])
        ell=h+L
        aa=baseline["inverse_alpha_U"]-(8.6*h+3.6*ell+3.2*s)/(2*math.pi)
        reconstructed=aa-baseline["inverse_alpha_U"]+(np.array([8.6,3,-1])*h
                      +np.array([3.6,4.,0.])*ell+np.array([3.2,0.,5.])*s)/(2*math.pi)
        np.testing.assert_allclose(reconstructed,[0.,0.,da3],atol=1e-12)
        salvaged[label]={"triplet_GeV":1e13,"sextet_GeV":U*math.exp(h-s),
                         "new_crossing_GeV":U*math.exp(h),"inverse_alpha_U":aa}
    # logs are affine in the inverse-coupling shifts. Including ell<0 only relaxes the bound.
    return {"ell_interval_including_unphysical": [x[:,0].min(),x[:,0].max()],
            "U_interval_GeV":[x[:,1].min(),x[:,1].max()],"triplet_min_GeV":x[:,2].min(),
            "triplet_plus_sextet_only_witnesses":salvaged,
            "min_MT_over_U":math.exp(-x[:,0].max())}


def limits(run,mixed,x,universal,step=.08):
    seed,poles,cols,ii=run.mixed_seed(x)
    y0=seed[:48].reshape(4,3,4)
    k=seed[48:57].reshape(3,3)
    t=run.lower
    kev=mixed["neutral"](t.conj().T@k@t.conj())[1]["common_light_scale_eV"]
    vu=246.22*run.flow["TAN"]/math.hypot(1,run.flow["TAN"])/math.sqrt(2)
    poles=poles/universal
    k0=ii*kev*universal*1e-9/vu**2
    y=np.zeros((4,3,7),complex)
    y[:,:,:4]=y0
    z=np.r_[y.ravel(),np.eye(3).ravel(),0j]
    def rhs(z,time):
        y=z[:84].reshape(4,3,7)
        G=z[84:93].reshape(3,3)
        u,e,n=(y[i]@y[i].conj().T for i in (0,2,3))
        g=run.g2(time)
        L=e+n+(3*np.trace(u).real+np.trace(n).real-.6*g[0]-3*g[1])*np.eye(3)
        return np.r_[run.flow["beta"](y,g,True).ravel(),(L@G).ravel(),4*np.linalg.norm(n,2)]/LOOP
    def advance(z,time,dt):
        a=rhs(z,time)
        b=rhs(z+dt*a/2,time+dt/2)
        c=rhs(z+dt*b/2,time+dt/2)
        d=rhs(z+dt*c,time+dt)
        return z+dt*(a+2*b+2*c+d)/6
    def coefficient(z):
        G=z[84:93].reshape(3,3)
        return np.linalg.norm(G@k0@G.T,2)
    ratio0=np.linalg.norm(k0,2)
    roots={}
    last=run.start
    comparison=1e13
    singulars=[]
    for target in sorted([*poles,comparison,run.stop]):
        nseg=math.ceil(math.log(target/last)/step)
        dt=math.log(target/last)/nseg
        for s in range(nseg):
            time=math.log(last/run.start)+s*dt
            nextz=advance(z,time,dt)
            for label,cap in (("unit_product",1.),("cap_product",CAP)):
                if label not in roots and run.start*math.exp(time+dt)*coefficient(nextz)>=cap:
                    a,b=0.,dt
                    for _ in range(34):
                        m=(a+b)/2
                        if run.start*math.exp(time+m)*coefficient(advance(z,time,m))<cap:
                            a=m
                        else:
                            b=m
                    roots[label]=run.start*math.exp(time+(a+b)/2)
            z=nextz
        if target in poles:
            j=np.where(poles==target)[0][0]
            G=z[84:93].reshape(3,3)
            z[:84].reshape(4,3,7)[3,:,4+j]=G@cols[:,j]
            singulars.append(coefficient(z))
        if target==comparison:
            compar=comparison*coefficient(z)
            y=z[:84].reshape(4,3,7)
            gain=z[84:93].reshape(3,3)
            le=np.linalg.svd(y[2],full_matrices=False)[0][:,::-1]
            mat=le.conj().T@gain@k0@gain.T@le.conj()
            fractions=np.sum(abs(mat)**2,axis=1)/np.linalg.norm(mat)**2
        last=target
    # Each appended mass starts block diagonal. An upper bound on the old block grows by exp(4 int ||Y||² dt).
    # At the last threshold the earlier sector remains much smaller; this also works when recomputed without its small poles.
    assert 3*poles[1]*math.exp(z[-1].real)<poles[2]
    return {"tree_right_scales_GeV":poles,"II_equivalent_mass_spectral_eV":ratio0*vu**2*1e9,
            "kappa_at_10PeV_GeV_inverse":ratio0,"kappa_at_last_RH_GeV_inverse":singulars[-1],"triplet_mass_at_product_GeV":roots,
            "conditional_charged_single_tag_fraction_at_MT_1e13":fractions,
            "product_at_1e13":compar,"product_at_original_U":run.stop*coefficient(z),
            "vR_common_f_necessary_upper_if_MT_above_all_RH_GeV":math.sqrt(CAP)/singulars[-1],
            "all_RH_wavefunction_log_norm_upper":z[-1].real}


def calculate():
    base=runpy.run_path(str(HERE/"../1800/check.py"))
    run=base["Running"]()
    mixed=runpy.run_path(str(HERE/"../1799/check.py"))
    inputs=read("../1799/inputs.json")
    cases=dict(inputs["mixed_type_I_II"],**inputs["fixed_delta_minus70_mass_readout"])
    universal=low_universal(run)
    gauge=gauge_window()
    results={}
    maximum_error=0.
    for case,x in cases.items():
        results[case]={}
        for tag,f in (("old_absolute_matching",1.),("restored_common_running",universal)):
            row=limits(run,mixed,x,f)
            check=limits(run,mixed,x,f,step=.04)
            maximum_error=max(maximum_error,abs(check["triplet_mass_at_product_GeV"]["cap_product"]/
                                               row["triplet_mass_at_product_GeV"]["cap_product"]-1))
            results[case][tag]=row
    # Use the worst of the independent weak coupling errors, before the isolated contact triplet threshold.
    ew=read("../1369/inputs.json")
    Q=gauge["U_interval_GeV"][1]
    maxg=run.g2(math.log(Q/run.start))**-1
    for i,key in enumerate(("g1_GUT","g2")):
        g,s=ew[key]
        maxg[i]+=(g+2*s)**-2-g**-2
    upper_g12=1/maxg[:2]
    eta=(.6*upper_g12[0]+3*upper_g12[1])*math.log(Q/run.start)/LOOP
    gauge["max_g12_squared_before_T"]=upper_g12
    gauge["kappa_min_ratio_to_10PeV_any_Yukawa"]=math.exp(-2*eta)
    gauge["max_10PeV_II_kappa_for_isolated_triplet_GeV_inverse"]=CAP*math.exp(2*eta)/gauge["triplet_min_GeV"]
    gauge["degenerate_15_inverse_alpha_at_1e13_original_U"]=4*math.pi/run.g2(math.log(run.stop/run.start))[0]-7/(2*math.pi)*math.log(run.stop/1e13)
    return {"low_common_factor":universal,"five_cases":results,"two_marginal_gauge_box":gauge,"max_relative_step_error":maximum_error}


def clean(a):
    if isinstance(a,dict):
        return {k:clean(x) for k,x in a.items()}
    if isinstance(a,(list,np.ndarray)):
        return [clean(x) for x in a]
    if isinstance(a,(float,np.floating)):
        return float(format(a,".10g"))
    return a


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    data=json.dumps(clean(calculate()),indent=2,allow_nan=False)+"\n"
    target=HERE/"independent_results.json"
    if args.check:
        assert json.loads(target.read_bytes())==json.loads(data)
    else:
        target.write_text(data,encoding="utf-8")
    print("1802 independent",hashlib.sha256(data.encode()).hexdigest()[:12])
