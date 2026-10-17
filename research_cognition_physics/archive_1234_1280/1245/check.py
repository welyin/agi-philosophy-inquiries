"""1245: exact-rational diagnostics for the joint classical ADM material.

The full-field constraint theorem is proved in proof.md.  Local jets here are
diagnostics, not a claim that a finite lattice represents the continuum algebra.
The 1242 scientific replay is inherited and is not a new scientific check group.
Default is read-only; --write creates results.json once using exclusive mode.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
from math import isqrt
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent


def dot(x, y):
    return sum((a*b for a, b in zip(x, y)), Q(0))


def matvec(a, x):
    return [dot(row, x) for row in a]


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    return [[dot(row, col) for col in transpose(b)] for row in a]


def contract(a, b):
    return sum((dot(row, col) for row, col in zip(a, b)), Q(0))


def inverse(a):
    n = len(a)
    aug = [list(row) + [Q(i == j) for j in range(n)]
           for i, row in enumerate(a)]
    for j in range(n):
        k = next(k for k in range(j, n) if aug[k][j])
        aug[j], aug[k] = aug[k], aug[j]
        pivot = aug[j][j]
        aug[j] = [x/pivot for x in aug[j]]
        for i in range(n):
            if i != j:
                factor = aug[i][j]
                aug[i] = [x-factor*y for x, y in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def sqrt_exact(x):
    x = Q(x)
    n, d = isqrt(x.numerator), isqrt(x.denominator)
    assert n*n == x.numerator and d*d == x.denominator
    return Q(n, d)


class Dual:
    """First derivative AD, all values and partial derivatives rational."""
    def __init__(self, value, gradient):
        self.value = Q(value)
        self.gradient = tuple(Q(x) for x in gradient)

    @classmethod
    def variable(cls, value, index, dimension):
        return cls(value, [int(i == index) for i in range(dimension)])

    def coerce(self, value):
        return value if isinstance(value, Dual) else Dual(value, [0]*len(self.gradient))

    def __add__(self, other):
        other = self.coerce(other)
        return Dual(self.value+other.value,
                    [a+b for a,b in zip(self.gradient, other.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-x for x in self.gradient])

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        return Dual(self.value*other.value,
                    [a*other.value+self.value*b
                     for a,b in zip(self.gradient, other.gradient)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        return Dual(self.value/other.value,
                    [(a*other.value-self.value*b)/(other.value**2)
                     for a,b in zip(self.gradient, other.gradient)])

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def __pow__(self, power):
        assert isinstance(power, int) and power >= 0
        out = self.coerce(1)
        for _ in range(power):
            out = out*self
        return out

    def sqrt(self):
        root = sqrt_exact(self.value)
        return Dual(root, [x/(2*root) for x in self.gradient])


H = [[Q(2), Q(1), Q(0)], [Q(1), Q(1), Q(0)],
     [Q(0), Q(0), Q(4)]]
HI = inverse(H)
SQRT_H = Q(2)


def field_density(q, momenta, gradients):
    F = 1+q[0]/5
    kinetic = momenta[0]**2/(2*SQRT_H)
    kinetic += sum((x**2/(2*SQRT_H*F) for x in momenta[1:]), Q(0))
    spatial = dot(gradients[0], matvec(HI, gradients[0])) + q[0]**2
    spatial += F*sum((dot(g, matvec(HI, g)) for g in gradients[1:]), Q(0))
    return kinetic+SQRT_H*spatial/2


def matter_group():
    assert matmul(H, HI) == [[Q(i == j) for j in range(3)] for i in range(3)]
    assert SQRT_H**2 == 4
    q = [Q(1,5), Q(1,7), Q(-1,9), Q(2,11)]
    momenta = [Q(2,7), Q(-3,8), Q(4,9), Q(-1,6)]
    gradients = [[Q(1,3), Q(-2,5), Q(1,7)],
                 [Q(-1,4), Q(2,9), Q(1,6)],
                 [Q(1,8), Q(1,5), Q(-1,3)],
                 [Q(-2,7), Q(1,11), Q(3,10)]]
    hessians = [[[Q((a+1)*(i+j+1), 31+i+j)
                  for j in range(3)] for i in range(3)] for a in range(4)]
    flat = q+momenta+[x for g in gradients for x in g]
    args = [Dual.variable(x, i, len(flat)) for i,x in enumerate(flat)]
    ad = field_density(args[:4], args[4:8],
                       [args[8+3*a:11+3*a] for a in range(4)])
    assert ad.value == field_density(q, momenta, gradients)
    F, Fp = 1+q[0]/5, Q(1,5)
    B = [Q(1), F, F, F]
    Bp = [Q(0), Fp, Fp, Fp]
    flux = [ad.gradient[8+3*a:11+3*a] for a in range(4)]
    div_flux = []
    for a in range(4):
        expected_flux = [SQRT_H*B[a]*x for x in matvec(HI, gradients[a])]
        assert tuple(expected_flux) == flux[a]
        div_flux.append(SQRT_H*sum((HI[i][j]*(B[a]*hessians[a][i][j]
                        +Bp[a]*gradients[0][i]*gradients[a][j])
                        for i in range(3) for j in range(3)), Q(0)))
        assert ad.gradient[4+a] == momenta[a]/(SQRT_H*B[a])
    field_interaction_derivative = ad.gradient[0]-SQRT_H*q[0]
    assert field_interaction_derivative != 0

    lapse_n, lapse_m = Q(5,4), Q(7,6)
    dn, dm = [Q(2,5), Q(-1,7), Q(3,11)], [Q(-1,3), Q(2,9), Q(1,8)]
    v = matvec(HI, [lapse_n*dm[i]-lapse_m*dn[i] for i in range(3)])
    variation_n = [lapse_n*(ad.gradient[a]-div_flux[a])-dot(dn, flux[a])
                   for a in range(4)]
    variation_m = [lapse_m*(ad.gradient[a]-div_flux[a])-dot(dm, flux[a])
                   for a in range(4)]
    field_pb = sum((variation_n[a]*lapse_m*ad.gradient[4+a]
                    -lapse_n*ad.gradient[4+a]*variation_m[a]
                    for a in range(4)), Q(0))
    dfield = [sum((momenta[a]*gradients[a][i] for a in range(4)), Q(0))
              for i in range(3)]
    assert field_pb == dot(v, dfield) and field_pb != 0
    field_legendre = dot(momenta, ad.gradient[4:8])-ad.value
    velocities = ad.gradient[4:8]
    field_lagrangian = SQRT_H/2*(velocities[0]**2
                       +F*dot(velocities[1:], velocities[1:])
                       -dot(gradients[0], matvec(HI, gradients[0]))-q[0]**2
                       -F*sum((dot(g, matvec(HI, g)) for g in gradients[1:]), Q(0)))
    assert field_legendre == field_lagrangian

    # Per-label physical momenta; fixed dmu=weight*dz occurs both in H and symplectic form.
    weight, action = Q(2,5), Q(3,2)
    mass = Q(1,10)+action*(1+q[0])/10
    particle_p = [3*mass/4, 3*mass/4, Q(0)]
    line_args = [Dual.variable(x,i,7) for i,x in enumerate([Q(0)]*3+particle_p+[action])]
    xad, pad, iad = line_args[:3], line_args[3:6], line_args[6]
    p_at_x = q[0]+dot(gradients[0], xad)
    mass_ad = Q(1,10)+iad*(1+p_at_x)/10
    energy_ad = (dot(pad, matvec(HI, pad))+mass_ad**2).sqrt()
    hn = weight*(lapse_n+dot(dn, xad))*energy_ad
    hm = weight*(lapse_m+dot(dm, xad))*energy_ad
    # PB in dz variables contains 1/weight; omitting it is a detectable error.
    line_pb = sum((hn.gradient[i]*hm.gradient[3+i]
                   -hn.gradient[3+i]*hm.gradient[i] for i in range(3)), Q(0))/weight
    dline = -weight*dot(v, particle_p)
    assert line_pb == dline and dline != 0
    wrong_unweighted_pb = weight*line_pb
    assert wrong_unweighted_pb != dline
    energy = energy_ad.value
    assert energy == 5*mass/4
    sigma = lapse_n*mass/energy
    vel_plus_shift = [lapse_n*x/energy for x in matvec(HI, particle_p)]
    assert sigma**2+dot(vel_plus_shift, matvec(H, vel_plus_shift)) == lapse_n**2
    beta = [Q(1,9), Q(-1,13), Q(1,15)]
    coordinate_vel = [vel_plus_shift[i]-beta[i] for i in range(3)]
    assert dot(particle_p, coordinate_vel)+mass*sigma == lapse_n*energy-dot(beta,particle_p)
    assert hn.gradient[6]/weight == sigma*(1+q[0])/10

    # Field-particle functional cross bracket.  The shared field PB has no label factor;
    # integrating the particle source retains exactly one factor weight.
    source_dp = weight*lapse_m*mass*(action/10)/energy
    cross_forward = -(lapse_n*momenta[0]/SQRT_H)*source_dp
    cross_reverse = +(lapse_m*momenta[0]/SQRT_H)*(weight*lapse_n*mass*(action/10)/energy)
    assert cross_forward != 0 and cross_forward+cross_reverse == 0
    total_material_pb = field_pb+line_pb+cross_forward+cross_reverse
    assert total_material_pb == dot(v,dfield)+dline

    # A one-way deletion cannot retain the old pointer rule as a Hamilton equation.
    # Rest motion isolates it without introducing an irrational comparison energy.
    proper_old_rate, proper_deleted_rate = (1+q[0])/10, Q(1,10)
    assert proper_old_rate-proper_deleted_rate == q[0]/10 != 0
    return {"status":"PASS", "scientific_group":1,
            "field_density":str(ad.value),
            "nonzero_Fp_interaction_derivative":str(field_interaction_derivative),
            "full_field_local_jet_PB":str(field_pb),
            "field_momentum_constraint_pairing":str(dot(v,dfield)),
            "label_density":str(weight), "worldline_energy":str(energy),
            "full_weighted_worldline_PB":str(line_pb),
            "worldline_momentum_constraint_pairing":str(dline),
            "wrong_PB_without_inverse_density":str(wrong_unweighted_pb),
            "nonzero_field_worldline_cross_terms":[str(cross_forward),str(cross_reverse)],
            "total_material_PB":str(total_material_pb),
            "field_Legendre_residual":"0", "worldline_Legendre_residual":"0",
            "proper_pointer_rate_gap_after_one_way_deletion":str(q[0]/10),
            "scope":"Exact local diagnostics. Full continuum identities follow from proof.md; no lattice closure or actual apparatus authorization."}


def gravity_legendre(h, hi, sqrt_h, k, a, b, lam, curvature, LambdaH):
    trace_k = contract(hi,k)
    kup = matmul(matmul(hi,k),hi)
    kval = contract(k,kup)
    lambda_l = lam/(3*lam-1)
    pi = [[sqrt_h/a*(kup[i][j]-lambda_l*trace_k*hi[i][j])
           for j in range(3)] for i in range(3)]
    pilow = matmul(matmul(h,pi),h)
    trace_pi = contract(h,pi)
    recovered_k = [[a/sqrt_h*(pilow[i][j]-lam*trace_pi*h[i][j])
                    for j in range(3)] for i in range(3)]
    assert recovered_k == k
    hkin = a/sqrt_h*(contract(pi,pilow)-lam*trace_pi**2)
    ham = hkin-b*sqrt_h*curvature+2*LambdaH*sqrt_h
    lag = 2*contract(pi,k)-ham
    expected = sqrt_h/a*(kval-lambda_l*trace_k**2)+b*sqrt_h*curvature-2*LambdaH*sqrt_h
    assert lag == expected
    return {"lambda_L":str(lambda_l), "kinetic_H":str(hkin),
            "Legendre_L":str(lag), "identity_residual":"0"}


def selected_group():
    k = [[Q(1,3),Q(1,7),Q(-1,11)], [Q(1,7),Q(-2,5),Q(1,13)],
         [Q(-1,11),Q(1,13),Q(3,8)]]
    generic = gravity_legendre(H,HI,SQRT_H,k,Q(7,3),Q(5,4),Q(2,5),Q(2,9),Q(1,17))
    selected = gravity_legendre(H,HI,SQRT_H,k,Q(6),Q(1,6),Q(1,2),Q(2,9),Q(1,17))
    assert generic["lambda_L"] == "2" and selected["lambda_L"] == "1"
    a,b,LambdaH = Q(6),Q(1,6),Q(0)
    kappa, Lambda_phys = a/2, a*LambdaH
    assert a*b == 1 and kappa == 3 and Lambda_phys == 0

    # Both bases of the SAME scalar+clock+record Hamiltonian; no square-root fitting.
    branches = []
    for storage, expected_rho in [(Q(1,5),Q(117,100)),(Q(0),Q(97,100))]:
        clock_mass = Q(7,10)+Q(1,10)+storage
        rotor_mass = Q(1,10)+Q(1,10)*Q(1)
        rho_clock, rho_rotor, rho_rods = clock_mass, Q(1,10)*rotor_mass, Q(3,20)
        rho = rho_clock+rho_rotor+rho_rods
        h2 = expected_rho
        assert rho == expected_rho
        hconstraint = -6*h2/a+2*LambdaH+rho
        assert hconstraint == 0
        branches.append({"storage_K":str(storage), "rho_clock":str(rho_clock),
                         "rho_recorder_including_internal_energy":str(rho_rotor),
                         "rho_rods":str(rho_rods), "rho_total_and_H_squared":str(rho),
                         "Hamiltonian_constraint":"0", "all_three_momentum_constraints":["0"]*3})
    # Nonzero-p off-shell source audit: the interaction must occur exactly once.
    p, action, number_density = Q(-1,1000), Q(1), Q(1,10)
    interaction_rho = number_density*action*p/10
    assert interaction_rho == Q(-1,100000)
    full_mass = Q(1,10)+action*(1+p)/10
    assert full_mass > 0 and (1+p)/10 > 0
    # 2 Lambda_H in the Hamiltonian corresponds to Lambda_phys=a Lambda_H.
    assert Q(1,6)*2*(Q(6)*Q(1,17)) == 2*Q(1,17)
    return {"status":"PASS", "scientific_group":2,
            "generic_before_selection":generic, "selected_tensor_Legendre":selected,
            "normalization":{"a":str(a),"b":str(b),"kappa":str(kappa),
                             "Lambda_H":str(LambdaH),"Lambda_phys":str(Lambda_phys)},
            "complete_initial_constraint_branches":branches,
            "nonzero_p_source_interaction_counted_once":str(interaction_rho),
            "scope":"Same full three-dimensional parent and the inherited 1242 bases. No inference of CT or of the low-order ansatz from current axioms."}


def inherited_check():
    code = STAGE/"1242"/"check.py"
    result_path = code.with_name("results.json")
    before = hashlib.sha256(result_path.read_bytes()).hexdigest()
    proc = subprocess.run([sys.executable,"-B","-X","utf8",str(code)],
                          cwd=RESEARCH.parent, capture_output=True,text=True,encoding="utf-8",check=True)
    assert hashlib.sha256(result_path.read_bytes()).hexdigest() == before
    old = json.loads(result_path.read_text(encoding="utf-8"))
    assert old["round"] == 1242 and old["status"] == "PASS"
    finite = old["group_finite_records"]
    gap = [Q(x) for x in finite["gap"]]
    remainder = [Q(x) for x in finite["remaining_separation"]]
    error = Q(finite["phase_error_per_history"])
    assert gap[0]-2*error == remainder[0] > 0
    assert gap[1]-2*error == remainder[1]
    assert old["group_deletions"]["fixed_geometry_record_gap_exact"] == "0"
    return {"status":"PASS", "new_scientific_groups":0,
            "replayed_round":1242, "stdout":proc.stdout.strip(),
            "results_sha256":before,
            "complete_readout_error_contract_per_history":str(error),
            "inherited_phase_gap":finite["gap"],
            "inherited_remaining_separation":finite["remaining_separation"],
            "inherited_B_separation":finite["B_separation"],
            "precision_remains_a_contract":True}


def hashes():
    relative = ["archive_342_369/research_note_351.md",
                "archive_342_369/research_note_352.md",
                "archive_1161_1200/1188/proof.md",
                "archive_1234_/1236/proof.md",
                "archive_1234_/1242/proof.md",
                "archive_1234_/1242/check.py",
                "archive_1234_/1242/results.json"]
    return {p:hashlib.sha256((RESEARCH/p).read_bytes()).hexdigest() for p in relative}


def run():
    return {"round":1245, "status":"PASS", "scientific_groups":2,
            "groups":[matter_group(),selected_group()],
            "inherited_replay":inherited_check(), "source_hashes":hashes(),
            "scope":{"new_axioms_adopted":0,"CT_adopted":False,
                     "conditional_joint_ADM_bridge":True,
                     "independent_horizon_entropy_proved":False,
                     "quantum_or_actual_readout_realization_proved":False,
                     "whole_GR_goal_proved":False,
                     "certified_arithmetic":"Fractions and exact first-derivative automatic differentiation"}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=run()
    destination=HERE/"results.json"
    if args.write:
        with destination.open("x",encoding="utf-8",newline="\n") as output:
            json.dump(result,output,ensure_ascii=False,indent=2,sort_keys=True)
            output.write("\n")
        print("WRITE_ONCE: 1245, 2 scientific groups PASS; inherited 1242 replay PASS")
    else:
        saved=json.loads(destination.read_text(encoding="utf-8"))
        assert saved == result,"read-only replay mismatch; no files changed"
        print("READ_ONLY_MATCH: 1245, 2 scientific groups PASS; inherited 1242 replay PASS")


if __name__ == "__main__":
    main()
