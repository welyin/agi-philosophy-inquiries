"""A charm absolute leptonic rate test of two complex aligned scalar coefficients."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
MESONS = [0.1395706, 0.493677, 5.27934, 1.86966]  # pi+, K+, B+, D+ [GeV]


def read(n, name="results.json"):
    return json.loads((ARCHIVE/str(n)/name).read_text(encoding="utf-8"))


def q(x):
    return F(str(x))


class Box:
    """Closed rational interval; all coefficient enclosure endpoints are exact."""
    def __init__(self, low, high=None):
        self.lo = q(low)
        self.hi = self.lo if high is None else q(high)
        assert self.lo <= self.hi

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __add__(self, other):
        return Box(self.lo+other.lo, self.hi+other.hi)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        edges = [x*y for x in (self.lo,self.hi) for y in (other.lo,other.hi)]
        return Box(min(edges),max(edges))

    def __truediv__(self, other):
        assert other.lo*other.hi > 0
        return self * Box(1/other.hi,1/other.lo)

    def modulus(self):
        return float(max(abs(self.lo), abs(self.hi)))


def affine(point, b_row, d_row, pi_row):
    bx,by = [x-y for x,y in zip(b_row,pi_row)]
    dx,dy = [x-y for x,y in zip(d_row,pi_row)]
    kx,ky = [x-y for x,y in zip(point,pi_row)]
    det = bx*dy-by*dx
    b = (kx*dy-ky*dx)/det
    d = (bx*ky-by*kx)/det
    return [Box(1)-b-d,b,d]


def coefficients(yukawa, allowance, parts, pole_floor=None):
    """Split each of four independent up/down matching ratios into rational intervals."""
    import itertools
    grid = [F(-1)+2*F(k,parts) for k in range(parts+1)]
    rows = []
    for meson,(up,down) in zip(MESONS, yukawa):
        ratio = q(up)/q(down)
        square = Box(q(meson)**2)
        choices = []
        for lo,hi in zip(grid,grid[1:]):
            x = Box(ratio*(1+allowance*lo),ratio*(1+allowance*hi))
            down_weight = square / (Box(1)+x)
            # Use correlation coefficient_u = M_P^2 - coefficient_d.
            up_weight = square-down_weight
            if pole_floor is not None:
                pole = Box(1) / (Box(1)-square*Box(0,F(1,pole_floor**2)))
                down_weight,up_weight = down_weight*pole,up_weight*pole
            choices.append((down_weight,up_weight))
        rows.append(choices)
    for pi,k,b,d in itertools.product(*rows):
        yield (affine(k,b,d,pi), affine((Box(0),Box(0)),b,d,pi))


def num(x):
    return float(f"{x:.10g}")


def data(n, scale_allowance, wp):
    """B/π amplitude bound from Belle II and the registered inclusive pion width."""
    pi,_,bm,_ = MESONS
    mu,tau = .1056583755,1.77686
    fpi,fb = 130.50+n*.13,189.4-n*1.4
    bbranch = (1.24+n*math.hypot(.41,.19))*1e-4
    bwidth_us = bbranch/(1.637-n*.004)*1e6
    bbound = wp["weak_Vud"]/wp["weak_Vub"]*fpi/fb*mu/tau
    bbound *= (1-(mu/pi)**2)/(1-(tau/bm)**2)
    bbound *= math.sqrt(pi/bm*bwidth_us/(38.408-n*.007)/(1-scale_allowance))
    # The measured semileptonic moduli admit this lower bound for every CKM phase.
    inputs = read(1702,"inputs.json")["direct_semileptonic_moduli"]
    cb = inputs["Vcb"][0] + .93*math.hypot(*inputs["Vcb"][1:])
    s13 = wp["weak_Vub"]
    c13 = math.sqrt(1-s13**2)
    s12,c12,s23 = wp["weak_Vus"]/c13,wp["weak_Vud"]/c13,cb/c13
    vcd_min = s12*math.sqrt(1-s23**2)-c12*s23*s13
    dbound = (48.02+n*math.hypot(.48,.24,.12,.15))/(212.7-n*.6)/vcd_min
    dbound /= math.sqrt(1-scale_allowance)   # separate rate allowance beyond paper
    kmin = (.23135-n*.00046)/wp["weak_Vus_over_Vud"]
    return dict(Bpi=bbound,Dabs=dbound,Kpi=kmin,Vcd_min=vcd_min)


def certificate(k,g,t):
    b,d,kmin = t["Bpi"], t["Dabs"], t["Kpi"]
    b_g = g[0].modulus()+g[1].modulus()*b
    dg = g[2].modulus()
    assert d*dg < 1
    max_d_over_pi = d*b_g/(1-d*dg)
    k_base = k[0].modulus()+k[1].modulus()*b
    kd = k[2].modulus()
    k_upper = k_base+kd*max_d_over_pi
    assert kd > 0
    needed_d_over_pi = max(0.,(kmin-k_base)/kd)
    needed_dabs = needed_d_over_pi/(b_g+dg*needed_d_over_pi)
    # Minimum extra multiplicative factor on the K/pi rate implied by this bound.
    return k_upper, needed_dabs, max(1.,(kmin/k_upper)**2)


def compute():
    wp = read(1734,"independent_results.json")
    point = read(1732)["exploratory_k3"]["registered_with_lower_us"]
    down = point["down_strange_bottom_MZ"]
    up = read(1386,"inputs.json")["yukawa"]["up"]
    pairs = [(up[0][0],value) for value in down]+[(up[1][0],down[0])]
    kc,gc = next(coefficients(pairs,q(0),1))
    empirical = {}
    for n in [1,2]:
        for allowance in [.05,.1]:
            t = data(n,allowance,wp)
            empirical[f"{n}width_{int(100*allowance)}pct_rate"] = t
    output = {}
    for width,parts,floor in [(0,1,None),(.02,4,None),(.02,4,80)]:
        studies = {tag: [-math.inf,math.inf,math.inf] for tag in empirical}
        for k,g in coefficients(pairs,q(width),parts,floor):
            for tag,t in empirical.items():
                values = certificate(k,g,t)
                acc = studies[tag]
                acc[0] = max(acc[0], values[0])
                acc[1] = min(acc[1], values[1])
                acc[2] = min(acc[2], values[2])
        for tag,nums in studies.items():
            assert nums[0] < empirical[tag]["Kpi"] and nums[1] > empirical[tag]["Dabs"]
        label = f"mass_ratio_width_{width}" + (f"_exact_real_pole_above_{floor}GeV" if floor else "")
        output[label] = {
            tag:dict(zip(["Kpi_amplitude_upper","Dabs_min_to_match_Kpi_and_Bpi",
                          "Kpi_extra_rate_factor_lower"],map(num,values)))
            for tag,values in studies.items()}
    # Test both reduced complex amplitude identities with independently chosen complex coefficients.
    rows = np.array([[1,p*p*d/(u+d),p*p*u/(u+d)] for p,(u,d) in zip(MESONS,pairs)])
    test = rows@np.array([1.,.038-.112j,2.11+1.14j])
    base = test[[0,2,3]]
    ks = np.array([float(x.lo) for x in kc])
    gs = np.array([float(x.lo) for x in gc])
    assert abs(ks@base-test[1]) < 1e-13 and abs(gs@base-1.) < 1e-13
    return {
      "scope": "R1732 registered weak point; no new vector current; two flavor-aligned complex coefficients with common mu/tau factor, local amplitudes",
      "data": {"D_mu_fD_times_Vcd_effective_MeV":[48.02,.48,.24,.12,.15],
               "D_plus_decay_constant_MeV":[212.7,.6],"Kpi_RM123_axial":[.23135,.00046]},
      "affine_K_on_pi_B_D": list(map(num,ks)),
      "affine_1_on_pi_B_D": list(map(num,gs)),
      "experimental_sensitivity_windows": {k:{field:num(x) for field,x in t.items()} for k,t in empirical.items()},
      "certificates": output}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    dest = HERE/"independent_results.json"
    results = compute()
    if args.check:
        if json.loads(dest.read_text(encoding="utf-8")) != results:
            raise SystemExit("1735 independent results differ")
        print("1735 independent: BESIII charm and Belle/pion bounds exclude the two-coefficient repair at the frozen point")
    else:
        dest.write_text(json.dumps(results,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
