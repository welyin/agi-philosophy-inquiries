"""Angular inversion for two Bs acceptances and a separate exclusive pion readout."""
import argparse
import cmath
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def circle_interval(ell, n, k, mu_coeff, w, window):
    # x=cos(arg(N/L)); this quadratic is strictly decreasing on the physical interval here.
    a=(1+k)*ell
    assert a>(1-w)*n
    def value(x):
        return mu_coeff*((a-n*x)**2+w*n*n*(1-x*x))
    def root(rate):
        if rate>=value(-1):
            return -1.
        if rate<=value(1):
            return 1.
        low,high=-1.,1.
        for _ in range(65):
            x=(low+high)/2
            if value(x)>rate:
                low=x
            else:
                high=x
        return (low+high)/2
    assert window[0]<value(-1) and window[1]>value(1)
    return [root(window[1]),root(window[0])]


def calculate():
    with (HERE/"results.json").open(encoding="utf-8-sig") as file:
        authored=json.load(file)
    Hmu=(1.1663787e-5**2*.2307**2*5.36688*.1056583755**2
         *math.sqrt(1-4*(.1056583755/5.36688)**2)*1.615/(16*math.pi*6.582119569e-13))
    Htau=(1.1663787e-5**2*.1894**2*1.777**2*5.27938
          *(1-(1.777/5.27938)**2)**2*1.639/(8*math.pi*6.582119569e-13))
    n=.04185*1.469*math.sqrt(4.6)/(127.952*math.pi*.23121)
    k=.64
    w1=(1.43/1.61)*(1.577-.358*1.61)/(1.577-.358*1.43)
    w2=1.43/(1.61*1.092)
    external_pi=[.00393-.00041,.00393+.00041]
    rows={}
    for tag,branch in authored.items():
        left_B=complex(*branch["left_B_amplitude"])
        F=complex(*branch["Fsb"])
        ell=abs(left_B)
        x1=circle_interval(ell,n,k,Hmu*1e9,w1,[3.18,4.54])
        x2=circle_interval(ell,n,k,Hmu*1e9,w2,[2.88,3.77])
        low,high=max(x1[0],x2[0]),min(x1[1],x2[1])
        assert low<high
        angles=[abs(cmath.phase(n*complex(x,math.sqrt(1-x*x))-ell)) for x in (high,low)]
        actual_angle=abs(cmath.phase(F/left_B))
        left=complex(*branch["left_ub"])
        right=complex(*branch["right_ud_us_ub"][2])
        pi,axial=abs(left+right),abs(left-right)
        assert abs(pi*pi+axial*axial-2*.0041**2)<1.e-13
        assert abs(right-complex(*branch["right_ud_us_ub"][1])*F/branch["Fss"])<1.e-12
        assert abs(abs(F)**2/branch["Fss"]-branch["Fbb"])<1.e-12
        tau=Htau*axial**2*1.e4
        assert abs(tau-branch["charged_tau_times_1e4"])<1.e-8
        p=(k*left_B-F)*left_B.conjugate()/ell
        pred=[Hmu*1.e9*(p.real**2+w*p.imag**2) for w in (w1,w2)]
        assert max(abs(t-u) for t,u in zip(pred,branch["fitted_Bs_times_1e9_CMS_Run3"]))<1.e-7
        survives=angles[0]<=actual_angle<=angles[1]
        assert survives==(3.18<=pred[0]<=4.54 and 2.88<=pred[1]<=3.77)
        rows[tag]={"fixed_L_allowed_abs_arg_F_over_L_rad":angles,"realized_angle":actual_angle,
                   "Bs_survives_both":survives,"pion_vector_amplitude":pi,
                   "pion_within_one_width":external_pi[0]<=pi<=external_pi[1],
                   "tau_e4":tau,"tau_within_one_width":.64<=tau<=1.84}
    resolution=abs(rows["1"]["pion_vector_amplitude"]-rows["-1"]["pion_vector_amplitude"])/2
    return {"rows":rows,"exclusive_pi_one_width":external_pi,
            "half_split_between_opposite_beta_pion_amplitudes":resolution}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    value=calculate()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as f:
            json.dump(value,f,indent=2)
            f.write("\n")
    else:
        assert value==json.loads(path.read_text(encoding="utf-8"))
        print("1467 independent OK")
