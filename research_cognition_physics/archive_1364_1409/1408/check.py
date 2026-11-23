"""A normalized charged-current inequality on the rank-two TM1 heavy interface."""
import argparse
import json
from math import hypot, sqrt
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_text(encoding="utf8"))


def box(data,k):
    pi, stat, syst=data["pion"]["R_times_1e4"]
    theory, terr=data["pion"]["SM_times_1e4"]
    tau, t_err=data["tau"]["g_tau_over_e_leptonic"]
    spread=k*hypot(stat,syst)
    rmin=(pi-spread)/(theory+k*terr)
    rmax=(pi+spread)/(theory-k*terr)
    return [2-1/rmin-(tau+k*t_err)**2,2-1/rmax-(tau-k*t_err)**2]


def run():
    data=read("inputs.json")
    windows={str(k):box(data,k) for k in [0,1,2]}
    assert windows["0"][1]<0<windows["1"][1]
    lo,hi=0.,1.
    for _ in range(58):
        mid=(lo+hi)/2
        if box(data,mid)[1]>0:
            hi=mid
        else:
            lo=mid
    n=np.array(data["null_direction"])/sqrt(6)
    a=np.ones(3)/sqrt(3)
    b=np.array([0,1.,-1])/sqrt(2)
    # Physical full-rank real Dirac block with orthogonal weak profiles.
    s13sq=(1-sqrt(1-read(data["older_plane"])["daya_bay"]["sin_sq_2theta13"]))/2
    st=sqrt(3*s13sq)
    ct=sqrt(1-st**2)
    plane=np.column_stack([ct*a-st*b, st*a+ct*b])
    no=read(data["mass_shapes"])["NO"]
    solar=no["solar_gap"][0]
    masses=np.sqrt([solar,solar+no["atmospheric_gap"][0]])
    heavy=np.asarray(data["heavy_mass_illustration_GeV"])*1e9
    eigen=masses/(heavy+masses)
    X=(plane*eigen)@plane.T
    d=X.diagonal()
    H=(b@X@b)/(1-d[0])
    assert np.isclose(H,(d[1]+d[2]-2*d[0])/(1-d[0]),rtol=1e-11,atol=1e-27)
    assert 0<H<windows["2"][1] and np.linalg.norm(X@n)<1e-27
    assert np.isclose((1-eigen[1])*plane[0,1]**2,s13sq,atol=1e-13)
    # A permitted diagonal-Gram example in the one-width window; only kinematics is asserted here.
    u=np.array([1.,1.5,.5])/sqrt(3.5)
    v=np.cross(n,u)
    loss=.0046*np.outer(u,u)+1e-7*np.outer(v,v)
    de,dm,dt=loss.diagonal()
    r=(1-de)/(1-dm)
    grel=sqrt((1-dt)/(1-de))
    pi,stat,syst=data["pion"]["R_times_1e4"]
    std,err=data["pion"]["SM_times_1e4"]
    tg,te=data["tau"]["g_tau_over_e_leptonic"]
    assert (pi-hypot(stat,syst))/(std+err)<r<(pi+hypot(stat,syst))/(std-err) and tg-te<grel<tg+te
    true_r=(1-d[0])/(1-d[1])
    true_g=sqrt((1-d[2])/(1-d[0]))
    assert (pi-2*hypot(stat,syst))/(std+2*err)<true_r<(pi+2*hypot(stat,syst))/(std-2*err)
    assert tg-2*te<true_g<tg+2*te
    # Compress the measured diagonal combination to the allowed null plane.
    slopes={}
    for p in (.64,2/3,.70):
        orth=np.array([[sqrt(1-p),0],[-sqrt(p),0],[0,1.]])
        compressed=orth.T@np.diag([-2.,1.,1.])@orth
        assert np.allclose(np.linalg.eigvalsh(compressed),[3*p-2,1])
        slopes[str(p)]=float(np.trace(compressed@np.diag([.02,.0005])))
    assert slopes["0.64"]<0<slopes["0.7"]
    joint=[]
    for k in [0,1,2]:
        solar_max=.3036+k*.0064
        reactor_max=.0851+k*.0024
        wmin=(1-solar_max)*(1+sqrt(1-reactor_max))/2
        joint.append(dict(k=k,minimum_normalized_electron_weight=wmin,
                          largest_external_electron_loss=1-2/(3*wmin)))
    return dict(H_boxes=windows,first_common_width_factor=hi,
                directional_witness_numerators=slopes,reactor_transfer=joint,full_mass_example=dict(
        heavy_eV=heavy.tolist(),light_nonzero_eV=masses.tolist(),X_eigenvalues=eigen.tolist(),H=H,
        pion_R_over_SM=(1-d[0])/(1-d[1]),tau_over_e=sqrt((1-d[2])/(1-d[0]))),
        kinematic_example=dict(row_losses=loss.diagonal().tolist(),pion_R_over_SM=r,tau_over_e=grel))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result=run()
    if p.parse_args().check:
        def check(got,expected):
            if isinstance(got,dict):
                assert got.keys()==expected.keys()
                for key in got:
                    check(got[key],expected[key])
            elif isinstance(got,list):
                assert len(got)==len(expected)
                for x,y in zip(got,expected):
                    check(x,y)
            else:
                assert np.allclose(got,expected,rtol=1e-10,atol=1e-26)
        check(result,read("results.json"))
        print("1408 normalized weak-current inequality, windows and exact mass witness verified")
    else:
        print(json.dumps(result,indent=2))
