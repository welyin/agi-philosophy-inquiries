"""Exact tree-amplitude checks in units m_h=v=1 (restore dimensions in the report)."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path


def canonical_potential(raw3, raw4, z1, z2):
    # At the stationary point with V''=Z=1.
    return raw3-3*z1/2, raw4-3*z1*raw3+(19*z1*z1/4-2*z2)


def amplitude(z, coupling3, coupling4, energy, quadratic, cubic):
    """The t,u Higgs poles are kept explicitly; z is squared cosine."""
    tscale = (energy-4)/2
    offset = 1+tscale
    t_plus_u = -2*tscale
    t_times_u = tscale**2*(1-z)
    kinematic_squares = energy**2+t_plus_u**2-2*t_times_u
    kinematic_product = energy*t_times_u
    exchange = 1/(energy-1) - 2*offset/(offset**2-tscale**2*z)
    return (-coupling4-coupling3**2*exchange
            +quadratic*kinematic_squares+cubic*kinematic_product)


def certified_bins():
    # Exact outer bounds for the author's two S=25 bins. Symmetry allows w >= 0.
    # Per-cell dyadic bounds keep the sum small while preserving its direction.
    cells, w_cells, grid = 128, 8, 2**35
    half = F(1,20)
    moments = []
    for center in [F(28,3), F(9,2)]:
        low_sum = high_sum = square_high_sum = 0
        for i in range(cells):
            qlo = center-half+2*half*i/cells
            qhi = qlo+2*half/cells
            clo, chi = (26-qhi)/2, (26-qlo)/2
            for j in range(w_cells):
                wl, wh = half*j/w_cells, half*(j+1)/w_cells
                # 2c/(c²-w²) decreases in c, increases with |w|.
                lo = 1/(qhi-1)+2*chi/(chi**2-wl**2)
                hi = 1/(qlo-1)+2*clo/(clo**2-wh**2)
                low_sum += (lo*grid).__floor__()
                high_sum += (hi*grid).__ceil__()
                square_high_sum += (hi**2*grid).__ceil__()
        moments.append([F(n,grid*cells*w_cells) for n in [low_sum,high_sum,square_high_sum]])
    delta_lower = moments[1][0]-moments[0][1]
    assert delta_lower>0
    errors = []
    for b in [11,15]:
        rate_sum_upper = sum(b*b+50*b*hi+625*square for _,hi,square in moments)
        errors.append(rate_sum_upper/(5000*delta_lower))
    gap = 4-sum(errors)
    assert gap>F(1,50)
    return {"rate_inversion_gap_certified_lower":(gap*10**8).__floor__()/10**8}


def calculate():
    results = {}
    for name, g3, g4 in [("SM",F(3),F(3)), ("GW",F(5),F(11))]:
        a, b = F(3,7), F(-5,9)
        z1, z2 = 2*a, 2*(b+a*a)
        v3, v4 = g3+3*a, g4+6*a*g3+3*a*a+4*b
        assert canonical_potential(v3,v4,z1,z2)==(g3,g4)
        threshold = amplitude(F(0),g3,g4,F(4),F(0),F(0))
        for alpha,beta in [(F(0),F(0)),(F(19,7),F(-3,5))]:
            terms = [amplitude(z,g3,g4,F(6),alpha,beta) for z in [F(0),F(1,2),F(1)]]
            second = terms[0]-2*terms[1]+terms[2]
            assert second==g3**2/21
        # A local four-derivative operator changes threshold, but the angular combination cancels.
        def four_derivative(z, s):
            t = (s-4)/2
            return ((s-2)**2+2*(t+2)**2+2*t*t*z)/2
        assert four_derivative(F(0),F(4))==6
        assert sum(w*four_derivative(z,F(6)) for w,z in [(1,F(0)),(-2,F(1,2)),(1,F(1))])==0
        results[name]=dict(threshold=str(threshold),angular_combination=str(second),
                           unnormalized_coordinate_V3=str(v3),unnormalized_coordinate_V4=str(v4))
    assert results["SM"]["threshold"]=="12" and results["GW"]["threshold"]=="92/3"
    results.update(certified_bins())
    return results


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    opts=ap.parse_args()
    result=calculate()
    if opts.check:
        saved=json.loads(Path(__file__).with_name("independent_results.json").read_text(encoding="utf-8"))
        assert result==saved
        print("1390 exact radial metric, threshold and angular pole projection verified.")
    else:
        print(json.dumps(result,indent=2))
