"""1244: finite radar-to-normal-coordinate certificates.

Default execution recomputes and compares without writing.  --write creates the
initial results exclusively.  All certified inequalities use Fraction arithmetic.
The de Sitter calibration is independent of the generic budget example.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent
GRID = 10**45


def pack(x: Q):
    return {"exact": str(x), "decimal_diagnostic": float(x)}


def floor_grid(x: Q) -> Q:
    return Q((x.numerator * GRID) // x.denominator, GRID)


def ceil_grid(x: Q) -> Q:
    return -floor_grid(-x)


class Interval:
    """Closed rational interval; endpoints may be outward-rounded for output."""

    def __init__(self, lo, hi=None):
        self.lo = Q(lo)
        self.hi = Q(lo if hi is None else hi)
        assert self.lo <= self.hi

    @staticmethod
    def get(value):
        return value if isinstance(value, Interval) else Interval(value)

    def __add__(self, other):
        other = self.get(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.get(other))

    def __rsub__(self, other):
        return self.get(other) - self

    def __mul__(self, other):
        other = self.get(other)
        values = (self.lo * other.lo, self.lo * other.hi,
                  self.hi * other.lo, self.hi * other.hi)
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.get(other)
        assert not other.lo <= 0 <= other.hi
        return self * Interval(1 / other.hi, 1 / other.lo)

    def rounded(self):
        return Interval(floor_grid(self.lo), ceil_grid(self.hi))

    def square(self):
        if self.lo <= 0 <= self.hi:
            return Interval(0, max(self.lo**2, self.hi**2))
        return Interval(min(self.lo**2, self.hi**2),
                        max(self.lo**2, self.hi**2))

    def output(self):
        result = self.rounded()
        return {"lower": pack(result.lo), "upper": pack(result.hi),
                "width_upper": pack(result.hi - result.lo)}


def exp_positive(x: Q, degree=24) -> Interval:
    """Positive Taylor series and geometric majorant of its omitted terms."""
    assert x >= 0 and x < degree + 2
    term = Q(1)
    partial = Q(1)
    for k in range(1, degree + 1):
        term *= x / k
        partial += term
    first_omitted = term * x / (degree + 1)
    tail = first_omitted / (1 - x / (degree + 2))
    return Interval(partial, partial + tail).rounded()


def tanh_positive(x: Q) -> Interval:
    e = exp_positive(2 * x)
    # The rational map is increasing; retain its correlated numerator/denominator.
    return Interval((e.lo - 1) / (e.lo + 1),
                    (e.hi - 1) / (e.hi + 1)).rounded()


def atanh_positive(x: Interval, degree=16) -> Interval:
    assert 0 <= x.lo <= x.hi < 1
    def partial(q):
        return sum((q**(2*k + 1) / (2*k + 1)
                    for k in range(degree + 1)), Q(0))
    tail = x.hi**(2*degree + 3) / ((2*degree + 3) * (1 - x.hi**2))
    return Interval(partial(x.lo), partial(x.hi) + tail).rounded()


def log1p_positive(x: Q, degree=20) -> Interval:
    """Alternating-series bounds, used only for the accelerated-observer deletion."""
    assert 0 < x < 1 and degree % 2 == 0
    lower = sum(((-1)**(k+1) * x**k / k
                 for k in range(1, degree + 1)), Q(0))
    return Interval(lower, lower + x**(degree + 1) / (degree + 1)).rounded()


def matrix_multiply(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(len(b))), Q(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return [list(row) for row in zip(*a)]


def generic_budget():
    eps, G, T, L = Q(1, 10000), Q(1, 1000), Q(1, 2), Q(1, 10)
    eta_tau, eta_n, eta_model = Q(1, 100000), Q(1, 10000), Q(1, 100000)
    boost_norm = Q(2)
    assert 0 <= eps <= Q(1, 4) and 8*G*T < 1

    # Squared rational certificates replace numerical square-root comparisons.
    cplus_squared = (1 + eps) / (1 - eps)
    assert cplus_squared <= Q(9, 4)
    assert cplus_squared <= (1 + 2*eps)**2
    assert 1 / (1 - eps**2) <= 4  # c+ - c- = 2 eps / sqrt(1-eps^2)
    frame_bound = 8*G*T / (1 - 8*G*T)
    assert exp_positive(4*G*T).hi - 1 <= frame_bound
    d = L * (1 + 2*eps)
    # |v|<=3/2, sqrt(3)<=2 and |(1,v)|_1<=4 yield 56G.
    acceleration_bound = 56*G
    assert (Q(2) + Q(3, 2)) * 16*G == acceleration_bound
    terms = {"null_time_weight": 4*eps*L,
             "parallel_frame": d*frame_bound,
             "null_bending": acceleration_bound*d*d/2}
    geometric = sum(terms.values(), Q(0))
    read_terms = {"two_clock_readings": 2*eta_tau,
                  "direction_reading": d*eta_n,
                  "remaining_model_contract": eta_model}
    reading = sum(read_terms.values(), Q(0))
    one_observer = geometric + reading
    two_observer = (1 + boost_norm) * one_observer
    assert one_observer < Q(762, 10**6)
    target = Q(23, 10000)
    assert two_observer < target

    # Exact norm and metric identities for rapidity log(2), not sampled boosts.
    lorentz = [[Q(5,4), Q(-3,4)], [Q(-3,4), Q(5,4)]]
    eta = [[Q(1), Q(0)], [Q(0), Q(-1)]]
    assert matrix_multiply(matrix_multiply(transpose(lorentz), eta), lorentz) == eta
    gram = matrix_multiply(transpose(lorentz), lorentz)
    # 4I-Lambda^T Lambda=(15/8) [[1,1],[1,1]] is positive semidefinite.
    psd = [[4*(i == j) - gram[i][j] for j in range(2)] for i in range(2)]
    assert psd == [[Q(15,8), Q(15,8)], [Q(15,8), Q(15,8)]]

    # Deleting free motion: exact Minkowski metric, but a uniformly accelerated
    # observer has radar radius log(1+a*r)/a at the event (t=0,x=r).
    # This is a scope witness, not a newly authorized accelerated apparatus.
    acceleration, r = Q(1), Q(1,10)
    accelerated_radius = log1p_positive(acceleration*r) / acceleration
    deletion_gap = r - accelerated_radius
    assert deletion_gap.lo > Q(46,10000)
    assert accelerated_radius.hi < T
    assert deletion_gap.lo > two_observer

    return {
        "status": "PASS", "scientific_group": 1,
        "inputs": {k: pack(v) for k,v in {
            "epsilon": eps, "Gamma_component_bound": G, "time_tube": T,
            "coordinate_radius": L, "eta_tau": eta_tau, "eta_direction": eta_n,
            "eta_model": eta_model, "boost_operator_norm": boost_norm}.items()},
        "derived": {"cplus_squared": pack(cplus_squared),
                    "frame_error_bound": pack(frame_bound), "d_and_Rmax": pack(d),
                    "coordinate_null_acceleration_bound": pack(acceleration_bound)},
        "geometric_terms": {k: pack(v) for k,v in terms.items()},
        "reading_terms": {k: pack(v) for k,v in read_terms.items()},
        "single_observer_error_bound": pack(one_observer),
        "two_observer_error_bound": pack(two_observer),
        "conditional_target": pack(target), "conditional_margin": pack(target-two_observer),
        "deletion": {"removed_assumption": "geodesic observer",
                     "background": "exact Minkowski metric; epsilon=Gamma=0",
                     "acceleration": pack(acceleration), "event_radius": pack(r),
                     "radar_radius": accelerated_radius.output(),
                     "strict_coordinate_gap": deletion_gap.output(),
                     "claim": "dropping free motion invalidates the curvature-only bound"},
        "scope": "Conditional geometric/readout budget; no certification of actual precision, active device transport or full-resource SR permissions."
    }


def de_sitter_calibration():
    H, emission, D = Q(1,10), Q(1,10), Q(2)
    gamma, gamma_v = (D + 1/D)/2, (D - 1/D)/2
    u = tanh_positive(H*emission/2)
    assert (D*D*u).hi < 1  # all events lie inside the relevant static-patch window
    tau_b = (2/H) * atanh_positive(D*u)
    tau_return = (2/H) * atanh_positive(D*D*u)
    assert tau_b.hi < tau_return.lo

    # The de Sitter embedding null condition reduces to this exact polynomial.
    assert D*D - 2*gamma*D + 1 == 0
    assert gamma*gamma - gamma_v*gamma_v == 1
    # Riemann normal coordinates of q on B's geodesic: x_A=tau_B*u_B.
    normal_t, normal_x = gamma*tau_b, gamma_v*tau_b
    # Its normal coordinates in B's frame are exactly (tau_B,0).
    assert gamma*gamma - gamma_v*gamma_v == 1
    assert -gamma_v*gamma + gamma*gamma_v == 0
    radar_t = (emission + tau_return)/2
    radar_x = (tau_return - emission)/2
    delta_t, delta_x = radar_t - normal_t, radar_x - normal_x
    assert delta_t.lo > 0 and delta_x.lo > 0

    # B's own-axis readout is its local proper clock; no direction at R=0 needed.
    # Lambda applied to A's radar tuple differs from B's tuple by these scalars.
    transported_time_error = emission + tau_return/4 - tau_b
    transported_space_error = tau_return/4 - emission
    assert transported_time_error.lo > Q(7,10**6)
    assert transported_space_error.lo > Q(12,10**6)
    a_error_squared = delta_t.square() + delta_x.square()
    passive_error_squared = (transported_time_error.square()
                             + transported_space_error.square())
    assert passive_error_squared.hi < D*D*a_error_squared.lo
    # Strict curvature corrections to the tangent-space SR predictions.
    b_delay = tau_b - D*emission
    return_delay = tau_return - D*D*emission
    assert b_delay.lo > Q(5,10**6)
    assert return_delay.lo > Q(5,10**5)
    assert b_delay.hi < Q(6,10**6) and return_delay.hi < Q(6,10**5)
    assert tau_b.hi-tau_b.lo < Q(1,10**40)
    assert tau_return.hi-tau_return.lo < Q(1,10**40)

    return {
        "status": "PASS", "scientific_group": 2,
        "inputs": {"H": pack(H), "A_emission_proper_time": pack(emission),
                   "Doppler_factor": pack(D), "relative_speed": pack(gamma_v/gamma)},
        "tanh_half_emission": u.output(), "B_reception_proper_time": tau_b.output(),
        "A_return_proper_time": tau_return.output(),
        "normal_coordinates_in_A": {"t": normal_t.output(), "x": normal_x.output()},
        "radar_coordinates_in_A": {"T": radar_t.output(), "R": radar_x.output()},
        "A_radar_minus_normal": {"time": delta_t.output(), "space": delta_x.output()},
        "boosted_A_radar_minus_B_local_clock": {
            "time": transported_time_error.output(), "space": transported_space_error.output()},
        "A_radar_error_squared": a_error_squared.output(),
        "two_description_mismatch_squared": passive_error_squared.output(),
        "tangent_SR_clock_corrections": {"B": b_delay.output(), "A_return": return_delay.output()},
        "normal_coordinate_relation": "exact Lorentz transformation, same exponential map and base point",
        "series_certification": {"exp_degree": 24, "atanh_last_index": 16,
                                 "outward_grid_denominator": str(GRID)},
        "scope": "Independent exact curved-geometry calibration. It is NOT asserted to satisfy group 1's tube budget or to implement a complete sourced signal/clock apparatus."
    }


def source_hashes():
    paths = [
        "04_认知时空操作公理与洛伦兹运动学_阶段论文.md",
        "archive_1161_1200/research_note_1168.md",
        "archive_1161_1200/research_note_1169.md",
        "archive_1161_1200/research_note_1175.md",
    ]
    return {p: hashlib.sha256((RESEARCH/p).read_bytes()).hexdigest() for p in paths}


def run():
    return {"round": 1244, "status": "PASS", "scientific_groups": 2,
            "groups": [generic_budget(), de_sitter_calibration()],
            "source_hashes": source_hashes(),
            "scope": {"new_axioms_adopted": 0, "CT_adopted": False,
                      "active_SR_transport_proved": False,
                      "complete_apparatus_or_GR_source_proved": False,
                      "all_certified_arithmetic": "exact rational interval arithmetic"}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="create results.json exclusively")
    args = parser.parse_args()
    result = run()
    path = HERE / "results.json"
    if args.write:
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, sort_keys=True)
            stream.write("\n")
        print("WRITE_ONCE: 1244, 2 scientific groups PASS")
    else:
        with path.open(encoding="utf-8") as stream:
            saved = json.load(stream)
        assert result == saved, "saved result mismatch; read-only replay made no changes"
        print("READ_ONLY_MATCH: 1244, 2 scientific groups PASS")


if __name__ == "__main__":
    main()
