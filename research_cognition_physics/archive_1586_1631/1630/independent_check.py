"""Exact type-III current projection plus an independent threshold and calibration check."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DT = np.array([1,1,1,2,2])
ACTIVE = np.diag([1,1,1,0,0])
ZCH = np.diag(DT)
# SU(2) raising/lowering norm for m=0 to m=-1, j=1, and
# m=+1/2 to m=-1/2, j=1/2. Overall g/sqrt(2) is stripped.
CC = np.diag([math.sqrt(j*(j+1)-m*(m-1)) for j,m in [(0.5,0.5)]*3+[(1,0)]*2])


def randu(rng, n):
    return np.linalg.qr(rng.standard_normal((n,n))+1j*rng.standard_normal((n,n)))[0]


def full_test(rng):
    # The triplet mass scale is common to its charged/neutral components.
    # The electron mass matrix is sampled independently of the triplet Yukawa.
    b=(rng.normal(size=(3,2))+1j*rng.normal(size=(3,2)))*rng.uniform(.06,.45)
    t=randu(rng,2)
    bare=np.array([3.,rng.uniform(3.7,10)])
    heavy=t@np.diag(bare)@t.T
    neutral=np.block([[np.zeros((3,3)),b],[b.T,heavy]])
    charged=np.block([[np.diag([.0001,.03,.1]),math.sqrt(2)*b.conj()],
                      [np.zeros((2,3)),heavy.conj()]])
    _, masses, vh=np.linalg.svd(neutral)
    un=vh[::-1].conj().T
    masses=masses[::-1]
    phase=np.exp(-.5j*np.angle(np.diag(un.T@neutral@un)))
    un=un*phase
    uc=np.linalg.svd(charged)[0][:,::-1]
    v=uc.conj().T@CC@un
    zc=uc.conj().T@ZCH@uc
    direct=un.conj().T@ACTIVE@un
    from_w=2*np.eye(5)-v.conj().T@v
    assert np.max(np.abs(zc-v@v.conj().T))<1e-13
    assert np.max(np.abs(from_w-direct))<1e-13
    assert np.max(np.abs(v.conj().T@(2*np.eye(5)-zc)@v-direct))<1e-13
    c=direct[:3,:3]
    h=np.linalg.eigvalsh(np.eye(3)-c)
    assert abs(h[0])<1e-14
    assert abs(h[1]*h[2]/np.linalg.det(c) * (masses[3]*masses[4]) /
               (masses[1]*masses[2])-1) < 3e-10
    loss=3-np.trace(c@c).real
    assert loss >= 2*h[2]-h[2]**2-1e-13
    heavy_poles=masses[3:]
    assert all(heavy_poles >= bare-1e-13)
    assert all(bare >= (1-h[2])*heavy_poles - h[2]*masses[2]-1e-13)
    light_w=v[:3,:3]
    naive=np.linalg.eigvalsh(np.eye(3)-light_w@light_w.conj().T)
    relative=3-np.trace(c@c).real
    return bool(naive[0]<-1e-5), relative


def one_family(x):
    theta=.5*math.atan(2*x)
    neutral=[math.cos(theta),-math.sin(theta)]
    phi=math.atan(math.sqrt(2)*x)
    charged=[math.cos(phi),-math.sin(phi)]
    current=charged[0]*neutral[0]+math.sqrt(2)*charged[1]*neutral[1]
    return current, neutral[0]**2, 1+charged[1]**2


def dimension_six_check():
    x=.0004
    w, c, zc=one_family(x)
    # Light charged has T3=−Zc/2, neutral has T3=C/2.
    assert abs((w**2-1)/x**2-1)<1e-5
    assert abs((1-c)/x**2-1)<1e-5
    assert abs((zc-1)/(2*x**2)-1)<1e-5
    # Two diagonal flavors. Only left-handed W amplitudes, light neutrinos
    # treated massless, flavor-changing charged Z absent in this example.
    we=one_family(0.015)[0]
    wm=one_family(0.025)[0]
    gf_ratio=we*wm
    assert gf_ratio>1 and abs((gf_ratio-1)/((.015**2+.025**2)/2)-1)<.01
    # At O(epsilon^0) one massless chiral fermion channel has width
    # K |T3|^2, K=(g/cW)^2 mZ/(24 pi). For i!=j, include both
    # Z->i- j+ and Z->j- i+ in a charge-summed visible branching fraction.
    zc12=.02+.03j
    charged_sum=2*abs(-zc12/2)**2
    neutral_one=1/4
    assert abs(charged_sum/neutral_one-2*abs(zc12)**2)<1e-16
    return round(gf_ratio,12)


def threshold_bounds():
    old=json.loads((HERE.parent/"1611"/"results.json").read_bytes())
    mu_low,mu_high=old["whole_triplet_geometric_mean_GeV"]
    box=json.loads((HERE.parent/"1625"/"inputs.json").read_bytes())["box"]
    a0,a1=box["abs31_eV2"]
    r0,r1=box["solar_over_abs31"]
    max_light_eV=math.sqrt(a1*(1+r1))
    hmax=1-math.sqrt(.98)
    # If mu_i>=200 GeV, each heavy pole is <=mu_i times the factor.
    factor=(1+hmax*(max_light_eV/1e9)/200)/(1-hmax)
    poles=[mu_low**2,(mu_high*factor)**2]
    # Invariant pdet(I-C)/det C; both numerator and denominator independent
    # of an arbitrary unitary choice of light basis.
    windows={}
    for name,prod in [("NO",[a0*math.sqrt(r0),a1*math.sqrt(r1)]),
                      ("IO",[a0*math.sqrt(1+r0),a1*math.sqrt(1+r1)])]:
        windows[name]=[prod[0]*1e-18/poles[1],prod[1]*1e-18/poles[0]]
    return {"pole_product_GeV2":[float(format(k,'.13g')) for k in poles],
            "pdet_over_det_C_windows":{n:[float(format(x,'.13g')) for x in vals]
                                       for n,vals in windows.items()}}


def atlas_ratio():
    # The single-muon published width was obtained by multiplying that ratio
    # by a lepton-universal reference width. Use an upper cap for the multiplier.
    ratio = (510.5 - 2*math.hypot(2.5,13.5))/84.2
    beta = math.sqrt(1-4*(.107/90)**2)
    light_muon_width_floor = .5 * beta**3
    trace_lower = .98*ratio*light_muon_width_floor
    cmin = math.sqrt(trace_lower-2)
    h = 1-cmin
    # Valid also if the two individual threshold-to-bare factors vary independently.
    old=json.loads((HERE.parent/"1611"/"results.json").read_bytes())
    tlo,thi=old["whole_triplet_geometric_mean_GeV"]
    box=json.loads((HERE.parent/"1625"/"inputs.json").read_bytes())["box"]
    a0,a1=box["abs31_eV2"]
    r0,r1=box["solar_over_abs31"]
    max_light_GeV=(a1*(1+r1))**.5 * 1e-9
    pole_lo=(tlo/1.10)**2
    pole_hi=(thi*(1.10+h*max_light_GeV/200)/cmin)**2
    no_upper=a1*math.sqrt(r1)*1e-18/pole_lo
    io_lower=a0*math.sqrt(1+r0)*1e-18/pole_hi
    assert 2 < trace_lower <= 3 and no_upper < io_lower
    return dict(raw_ratio_lower=round(ratio,12), trace_C2_floor=round(trace_lower,12),
                h_upper=round(h,12), no_upper=float(format(no_upper,'.12g')),
                io_lower=float(format(io_lower,'.12g')))


def run():
    r=np.random.default_rng(1630)
    out=[full_test(r) for _ in range(72)]
    failures=sum(1 for negative,_ in out if negative)
    assert failures>=65
    w, c, _=one_family(.2)
    assert (1-w*w)<0<1-c
    return {"exact_current_and_singular_value_comparisons":len(out),
            "naive_light_only_deficit_negative":failures,
            "one_family_x_0p2":{"naive_deficit":round(1-w*w,11),"actual":round(1-c,11)},
            "pure_left_muon_fermi_ratio_diagonal_example":dimension_six_check(),
            "atlas_single_muon_ratio":atlas_ratio(),
            "matching_equal_threshold_condition":threshold_bounds()}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    file=HERE/"independent_results.json"
    if args.check:
        assert json.loads(file.read_bytes())==result
    else:
        content=json.dumps(result,indent=2)+"\n"
        with file.open("x",encoding="utf-8") as f:
            f.write(content)
    print("1630 independent SU(2), chiral current, pole, threshold and calibration checks passed")
