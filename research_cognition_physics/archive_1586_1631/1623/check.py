"""Finite axial Bhabha scattering and common-normalisation angular certificates."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT/"inputs.json").read_text(encoding="utf-8"))
BINS = DATA["lep"]["bins"]


class I:
    """Ordered bounds with array broadcasting for real operations."""
    def __init__(self, lo, hi=None):
        self.lo, self.hi = lo, lo if hi is None else hi

    def __add__(self, other):
        if not isinstance(other, I):
            other = I(other)
        return I(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, I):
            other = I(other)
        a, b, c, d = self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi
        return I(np.minimum(np.minimum(a,b),np.minimum(c,d)),
                 np.maximum(np.maximum(a,b),np.maximum(c,d)))

    __rmul__ = __mul__

    def inv(self):
        assert np.all(self.lo*self.hi > 0.)
        return I(1/self.hi, 1/self.lo)

    def square(self):
        a, b = self.lo**2, self.hi**2
        return I(np.where(self.lo*self.hi > 0., np.minimum(a,b),0.), np.maximum(a,b))


def components(energy, angle, r, width):
    """Each helicity component has a SM amplitude and coefficient of gF²/2."""
    cfg = DATA["couplings"]
    x, y = cfg["gY"]**2, cfg["g2"]**2
    sw, electric = x/(x+y), x*y/(x+y)
    left, right = math.sqrt(x+y)*(sw-.5), math.sqrt(x+y)*sw
    invs = energy.square().inv() if isinstance(energy,I) else I(1/energy**2)
    rz, wz = invs*cfg["mZ"]**2, invs*(cfg["mZ"]*cfg["GammaZ"])
    zp = I(1)-rz
    iz = (zp.square()+wz.square()).inv()
    zre, zim = zp*iz, -wz*iz
    k = (I(1)-angle)*.5
    photon_t = -k.inv()
    z_t = -(k+rz).inv()
    x_t = -(k+r).inv()
    pole = I(1)-r
    w = r*width
    norm_inv = (pole.square()+w.square()).inv()
    xr, xi = pole*norm_inv, -w*norm_inv
    terms = []
    for coef in (left**2, right**2):
        terms.append(((I(1)-k).square(), electric+coef*zre + electric*photon_t+coef*z_t,
                      coef*zim, xr+x_t, xi))
    terms.append((2*k.square(), electric+left*right*zre, left*right*zim, -xr, -xi))
    terms.append((I(2), electric*photon_t+left*right*z_t, I(0.), -x_t, I(0.)))
    return terms


def density(terms, g):
    value = I(0.)
    for weight, sr, si, ar, ai in terms:
        value += weight*((sr+ar*g).square()+(si+ai*g).square())
    return value


def kernel(terms,g):
    v=density(terms,g)
    return I(np.mean(v.lo, axis=-1), np.mean(v.hi, axis=-1))


def table(energy):
    values = np.array(DATA["lep"][str(energy)])
    edge = DATA["benchmark"]["box_widths"]*values[:,1]
    return (values[:,0]-edge)/values[:,2], (values[:,0]+edge)/values[:,2]


def fit(energy, ratio_lo, ratio_hi):
    dl, du = table(energy)
    mi = dl[:,None]/ratio_hi
    ma = du[:,None]/ratio_lo
    # Lower bound on max_i lower_i / min_j upper_j for each physical mass cell.
    return np.max(mi, axis=0)/np.min(ma, axis=0)


def angular_shapes(energy, masses, g, width=0.):
    r = I((np.asarray(masses)[:,None]/energy)**2)
    nodes, weights = np.polynomial.legendre.leggauss(32)
    answers = []
    for low, high in BINS:
        z = I(.5*(high+low)+.5*(high-low)*nodes)
        parts = components(energy, z, r, I(width))
        # kernel averages nodes uniformly. Include Gauss integration weights.
        parts = [(w*weights*16, sr, si, ar, ai) for w,sr,si,ar,ai in parts]
        answers.append(kernel(parts,I(g*g/2)).lo/kernel(parts,I(0.)).lo)
    return np.array(answers)


def certify(energy, g_band, bins=120, mass_steps=128):
    mass = np.linspace(0.,64.,mass_steps+1)/energy
    r = I(mass[:-1,None]**2, mass[1:,None]**2)
    width = I(0., DATA["benchmark"]["max_Gamma_over_m"])
    prepared, refs = [], []
    for low, high in BINS:
        grid = np.linspace(low,high,bins+1)
        terms = components(energy, I(grid[:-1],grid[1:]),r,width)
        prepared.append(terms)
        # Independent pure SM evaluation; bin integration encloses the exact integral.
        refs.append(kernel(terms,I(0.)))
    n = math.ceil((g_band[1]-g_band[0])/.002)
    gg = np.linspace(*g_band,n+1)
    best, absolute = float("inf"), float("inf")
    for a,b in zip(gg[:-1],gg[1:]):
        ratios = [kernel(term,I(a*a/2,b*b/2))*sm.inv()
                  for term,sm in zip(prepared, refs)]
        low = np.stack([v.lo for v in ratios])*(1-1e-9)
        high = np.stack([v.hi for v in ratios])*(1+1e-9)
        best = min(best,float(np.min(fit(energy,low,high))))
        absolute = min(absolute, float(np.min(low[0])))
    return {"min_relative_response_span":best,"max_symmetric_per_bin_error_before_overlap":
            (best-1)/(best+1),"min_absolute_backward_bin_over_SM":absolute}


def point(energy, m, g):
    ratios = angular_shapes(energy, [m], g)[:,0]
    lo, hi = table(energy)
    norm_lo, norm_hi = np.max(lo/ratios), np.min(hi/ratios)
    return {"ratios":ratios.tolist(), "normalization_limits":[float(norm_lo),float(norm_hi)],
            "critical_pair_zero_based":[int(np.argmax(lo/ratios)),int(np.argmin(hi/ratios))],
            "required_span":float(norm_lo/norm_hi)}


def one_sided_radiation():
    # Both lepton directions are measured, with |pi-theta1-theta2| < 10 degrees.
    z0,z1=BINS[0]
    acol=math.radians(10)
    theta=math.acos(-z0)
    beta=math.sin(acol)/(math.sin(theta)+math.sin(theta-acol))
    emin=200*math.sqrt((1-beta)/(1+beta))
    zmin,zmax=(z0-beta)/(1-z0*beta),(z1+beta)/(1+z1*beta)
    ms=np.linspace(0,64,65)[:,None]
    zs=np.linspace(zmin,zmax,65)
    es=np.linspace(emin,200,25)
    gs=np.linspace(.62,.65,7)
    lower=float('inf')
    for e1,e2 in zip(es[:-1],es[1:]):
        r=I((ms[:-1]/e2)**2,(ms[1:]/e1)**2)
        terms=components(I(e1,e2),I(zs[:-1],zs[1:]),r,I(0.,.10))
        sm=density(terms,I(0.)).hi
        for g1,g2 in zip(gs[:-1],gs[1:]):
            lower=min(lower,float(np.min(density(terms,I(g1*g1/2,g2*g2/2)).lo/sm)))
    lower *= 1-1e-9
    upper=table(200)[1][0]
    return {'max_abs_boost':beta,'min_energy_GeV':emin,'cm_cos_box':[zmin,zmax],
            'pointwise_rate_ratio_floor':lower,'two_width_abs_upper_over_SM':upper,
            'max_supported_fraction_times_absolute_transfer':upper/lower}


def run():
    weak = DATA["benchmark"]["gF"]
    result = {}
    for energy in (189,200):
        result[str(energy)] = {
            "g0":point(energy,0.,0.),
            "benchmark_g_0p632":{str(m):point(energy,m,.632) for m in (0,20,32,48,64)},
            "weak_band":certify(energy,weak),
        }
    result["189"]["additional_0p2_to_0p62"] = certify(189,DATA["benchmark"]["broader_gF"])
    result['200']['single_beam_collinear_ISR']=one_sided_radiation()
    assert result["189"]["g0"]["required_span"] < 1.
    assert result["189"]["weak_band"]["min_relative_response_span"] > 1.1
    assert result["200"]["weak_band"]["min_absolute_backward_bin_over_SM"] > 3.
    # QED crossing check: the spin-averaged squared pure photon amplitude divided by e^4.
    for z in (-.8,.1,.85):
        t, u = -(1-z)/2, -(1+z)/2
        helicity = 2*u*u*(1+1/t)**2 + 2*t*t+2/t**2
        textbook = 2*((u*u+t*t)+(u*u+1)/t**2+2*u*u/t)
        assert abs(helicity-textbook)<1e-10
    return result


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    return float(f"{x:.8g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=stable(run())
    out=ROOT/"results.json"
    if args.check:
        if json.loads(out.read_text(encoding="utf-8"))!=result:
            raise SystemExit("1623 finite propagator check mismatch")
    else:
        with out.open("x",encoding="utf-8") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("1623 LEP finite shape and rate enclosures verified")
