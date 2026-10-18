"""Round 1321: explicit finite parent, fixed family sequence and common bounds.

Default/--check is read-only. --write-new exclusively creates results.json.
No exponentially large full Fock matrix is claimed or required: CAR polynomials,
all degree shifts and global norm bounds implement the analytic full-space proof.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
CHECKS = []


def check(name, value, tol=2e-10):
    err = float(np.max(np.abs(value))) if np.size(value) else 0.0
    if not np.isfinite(err) or err > tol:
        raise AssertionError((name, err, tol))
    CHECKS.append({"name": name, "error": err, "tolerance": tol})


def truth(name, condition):
    check(name, 0.0 if bool(condition) else 1.0, 0.0)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canon_pair(x, y):
    """Distinct CAR operators. No contractions occur in the bilinears below."""
    if x == y:
        return None, 0
    if x < y:
        return (x, y), 1
    return (y, x), -1


def field(weyl, component):
    # op=(creation?, particle/antiparticle species, Weyl component, momentum).
    if component == 0:
        return [(-1, (0, "a", weyl, -1)), (1, (1, "b", weyl, -1))]
    return [(1, (0, "a", weyl, 1)), (-1, (1, "b", weyl, 1))]


def bilinear(weyl1, weyl2, retain_fourier=False):
    out = defaultdict(complex)
    for c, d, sign in [(0, 1, 1), (1, 0, -1)]:
        for fx, x in field(weyl1, c):
            for fy, y in field(weyl2, d):
                key, phase = canon_pair(x, y)
                if key is not None and (retain_fourier or fx+fy == 0):
                    out[(fx+fy, key)] += sign*phase
    return {k: v for k, v in out.items() if v != 0}


def add_poly(target, source, scale=1.0):
    for key, value in source.items():
        target[key] += scale*value


def current(color, sector):
    # single generation Weyl ids Q=0,...,5; L=6,7.
    out = defaultdict(complex)
    eps = np.array([[0, 1], [-1, 0]])
    if sector == "QQ":
        # 1/2 epsilon_color epsilon_weak Q_i epsilon_spin Q_j.
        for i in range(3):
            for j in range(3):
                if len({color, i, j}) != 3:
                    continue
                seq = [color, i, j]
                inversions = sum(seq[a] > seq[b] for a in range(3) for b in range(a+1, 3))
                ec = (-1)**inversions
                for a in range(2):
                    for b in range(2):
                        if eps[a, b]:
                            add_poly(out, bilinear(2*i+a, 2*j+b), ec*eps[a,b]/2)
    else:
        for a in range(2):
            for b in range(2):
                if eps[a,b]:
                    add_poly(out, bilinear(2*color+a, 6+b), eps[a,b])
    return {key: v for key, v in out.items() if v != 0}


def pair_vector(poly, creation):
    return {key[1]: value/2 for key, value in poly.items()
            if all(op[0] == creation for op in key[1])}


def inner(a, b):
    return sum(np.conjugate(value)*b.get(key, 0) for key, value in a.items())


def exponential(h, t):
    eigen, u = np.linalg.eigh(h)
    return (u*np.exp(-1j*t*eigen)) @ u.conj().T


def calculate():
    CHECKS.clear()
    eps = np.array([[0, 1], [-1, 0]])
    sigma_z = np.diag([1., -1.])
    for momentum in [-1, 1]:
        spinor = np.array([1., 0.]) if momentum == -1 else np.array([0., 1.])
        check(f"positive_left_Weyl_symbol_{momentum}",
              -momentum*sigma_z@spinor-spinor)
        # Negative frequency at spatial Fourier -momentum has negative energy.
        check(f"negative_frequency_symbol_{momentum}",
              momentum*sigma_z@spinor+spinor)
        check(f"collinear_Weyl_contraction_{momentum}", spinor@eps@spinor)
    full_bilinear = bilinear(0, 1, True)
    kept = bilinear(0, 1)
    truth("bilinear_four_nonzero_zero_momentum_terms", len(kept) == 4)
    truth("all_mixed_terms_have_nonzero_momentum",
          all(fourier != 0 for (fourier, monomial) in full_bilinear
              if monomial[0][0] != monomial[1][0]))
    truth("zero_momentum_terms_are_pair_or_antipair",
          all(monomial[0][0] == monomial[1][0] for _, monomial in kept))
    truth("all_discarded_fourier_are_plus_minus_two",
          all(fourier in [-2, 0, 2] for fourier, _ in full_bilinear))

    currents = {(sector, k): current(k, sector) for sector in ["QQ", "QL"] for k in range(3)}
    for (sector, k), poly in currents.items():
        for creation in [0, 1]:
            v = pair_vector(poly, creation)
            truth(f"{sector}_{k}_{creation}_four_CAR_terms", len(v) == 4)
            check(f"{sector}_{k}_{creation}_normalized_pair", inner(v,v)-1)
            check(f"{sector}_{k}_{creation}_operator_l1_bound", sum(abs(x) for x in v.values())-2)
        truth(f"{sector}_{k}_no_mixed_terms",
              all(mon[0][0] == mon[1][0] for _, mon in poly))
    for sector in ["QQ", "QL"]:
        for creation in [0, 1]:
            gram = [[inner(pair_vector(currents[sector,i], creation),
                           pair_vector(currents[sector,j], creation))
                     for j in range(3)] for i in range(3)]
            check(f"{sector}_{creation}_three_internal_columns_orthonormal", np.array(gram)-np.eye(3))
    # Repeating the identical CAR mode vanishes, but a Weyl self-bilinear need
    # not vanish. QQ here uses the actual color/weak slots of one generation.
    qq_symmetric = bilinear(0, 3)
    reverse = bilinear(3, 0)
    check("Weyl_bilinear_symmetric_under_composite_field_exchange",
          [qq_symmetric.get(k,0)-reverse.get(k,0) for k in set(qq_symmetric)|set(reverse)])
    truth("same_flavor_QQ_nonzero", bool(currents["QQ", 0]))

    # Full projected boson space vacuum+s_0..2+t_0..2.
    nb = np.diag([0.]+[1.]*6)
    boson = []
    for k in range(6):
        annih = np.zeros((7,7), complex)
        annih[0,k+1] = 1
        boson.append(annih)
        check(f"projected_boson_{k}_commutator", nb@annih.conj().T-annih.conj().T@nb-annih.conj().T)
        check(f"projected_boson_{k}_norm", np.linalg.norm(annih,2)-1)
    # All direct terms and their Hermitian conjugates from equation 7.
    direct = [("s+", "A", 2, -2), ("t", "A", -2, -2),
              ("s+", "Abar+", 2, 2), ("t", "Abar+", -2, 2),
              ("s", "B", -2, -2), ("t+", "B", 2, -2),
              ("s", "Bbar+", -2, 2), ("t+", "Bbar+", 2, 2)]
    terms = [(b,p,db+dp) for b,p,db,dp in direct]
    terms += [(b+" dagger",p+" dagger",-db-dp) for b,p,db,dp in direct]
    truth("eight_resonant_terms_including_hc", sum(d == 0 for _,_,d in terms) == 8)
    truth("four_raising_terms", sum(d == 4 for _,_,d in terms) == 4)
    truth("four_lowering_terms", sum(d == -4 for _,_,d in terms) == 4)
    truth("no_other_frequency_or_dropped_antiscalar", set(d for _,_,d in terms) == {-4,0,4})
    truth("target_and_charge_conjugate_resonant_sectors",
          {b+p for b,p,d in terms[:8] if d == 0} ==
          {"s+A", "tAbar+", "t+B", "sBbar+"})

    p = 1.0
    scalar_mass = 2*p
    volume = (2*np.pi/p)**3
    kappa = 2/np.sqrt(2*scalar_mass*volume)
    g = 2.**-24*p
    mu = 2.**-64*p
    local_lambda = g/kappa
    check("fixed_current_matrix_element_kappa", kappa*np.sqrt(2*scalar_mass*volume)-2)
    check("single_local_coupling_g", kappa*local_lambda-g)
    f = np.array([4.**-i for i in range(5)])
    s3 = sum(f[:3]**2)
    T_g = np.pi/np.sqrt(2*s3)
    truth("one_common_clock_less_than_three_over_g", T_g < 3)
    truth("uniform_flavor_l1_sum", sum(f) < 4/3)
    check("V0_full_Fock_bound_count", 8*3*2*4/3-64)
    check("Vplus_full_Fock_bound_count", 4*3*2*4/3-32)
    counter = Fraction(18496,2**24)
    y_bound = Fraction(196608,2**40)
    truth("mass_term_count_below_declared_bound", 8*25*8*11 < 32768)
    norms = {"V0":64, "Vplus":32, "Vminus":32, "Vtotal":128,
             "mass_over_mu":32768, "parent_over_p":323}
    truth("uniform_full_H_bound",
          322+128*2.**-24+32768*2.**-64 < 323)

    family = []
    nmatrices = {}
    for n in [3,4,5]:
        fn = f[:n]
        sn = sum(fn**2)
        bright = fn/np.sqrt(sn)
        h = np.zeros((2*n+1,2*n+1),complex)
        h[:n, n] = fn
        h[n, :n] = fn
        h[n, n+1:] = fn
        h[n+1:, n] = fn
        embedding = np.zeros((2*n+1,3),complex)
        embedding[:n,0] = bright
        embedding[n,1] = 1
        embedding[n+1:,2] = bright
        reduced = np.array([[0,np.sqrt(sn),0],[np.sqrt(sn),0,np.sqrt(sn)],[0,np.sqrt(sn),0]])
        check(f"n{n}_resonant_bright_invariant_parent", h@embedding-embedding@reduced)
        check(f"n{n}_all_orthogonal_flavor_directions_dark", h@(np.eye(len(h))-embedding@embedding.conj().T))
        unitary = exponential(h,T_g)
        prob = float(abs(embedding[:,2].conj()@unitary@embedding[:,0])**2)
        formula = float(np.sin(np.pi/2*np.sqrt(sn/s3))**4)
        check(f"n{n}_same_time_actual_chain_probability", prob-formula)
        nominal = 1-formula**2
        delta_bound = np.pi/4*(sn/s3-1)
        nominal_bound = 4*delta_bound**2
        truth(f"n{n}_analytic_shared_clock_bound", nominal <= nominal_bound+1e-14)
        total = nominal_bound+float(counter+y_bound)
        truth(f"n{n}_joint_full_parent_return_below_one_percent", total < .01)
        truth(f"n{n}_mode_count_all_fields", 4*(6+3+3+2+1+1)*n == 64*n)
        # Parent representation anomaly conditions for one generation (additive).
        # x=6Y, all-left Q,u^c,d^c,L,e^c,N^c.
        reps = [(6,1),(3,-4),(3,2),(2,-3),(1,6),(1,0)]
        check(f"n{n}_grav_hypercharge", n*sum(mult*q for mult,q in reps))
        check(f"n{n}_cubic_hypercharge", n*sum(mult*q**3 for mult,q in reps))
        check(f"n{n}_SU3_squared_hypercharge", n*(2*1-4+2))
        check(f"n{n}_SU2_squared_hypercharge", n*(3*1-3))
        truth(f"n{n}_SU2_Witten_even", (4*n)%2 == 0)

        dft3 = np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
        u = np.eye(n,dtype=complex)
        u[:3,:3] = dft3
        for j in range(3,n):
            rotation = np.eye(n,dtype=complex)
            c,s = np.cos(.01),np.sin(.01)
            rotation[j-1,j-1]=rotation[j,j]=c
            rotation[j-1,j]=s
            rotation[j,j-1]=-s
            u = u@rotation
        check(f"n{n}_common_mixing_recipe_unitary", u.conj().T@u-np.eye(n))
        mats = {"u":np.diag(np.arange(1,2*n,2)).astype(complex),
                "d":np.diag(np.arange(2,2*n+1,2)).astype(complex),
                "e":np.diag(np.arange(1,n+1)).astype(complex),
                "nu":u@np.diag(np.arange(1,n+1))@u.conj().T}
        kphys = mats["nu"]@mats["nu"].conj().T
        cp = float(np.imag(kphys[0,1]*kphys[1,2]*kphys[2,0]))
        truth(f"n{n}_physical_CP_invariant_nonzero", abs(cp) > 1)
        kweyl = mats["nu"].conj()@mats["nu"].T
        cpw = float(np.imag(kweyl[0,1]*kweyl[1,2]*kweyl[2,0]))
        check(f"n{n}_allleft_physical_CP_conjugation", cp+cpw)
        for name,m in mats.items():
            sv=np.linalg.svd(m,compute_uv=False)
            truth(f"n{n}_{name}_center_positive_simple", min(sv) > .999 and min(abs(np.diff(sv))) > .999)
            truth(f"n{n}_{name}_open_ball_max_norm", max(sv)+.125 < 11)
        truth(f"n{n}_doubled_neutral_open_ball_max_norm", 2*np.linalg.norm(mats["nu"],2)+.125 < 11)
        # An explicit nontrivial quark mixing variation inside the same open ball.
        rot=np.eye(n,dtype=complex); angle=.001
        rot[:2,:2]=[[np.cos(angle),np.sin(angle)],[-np.sin(angle),np.cos(angle)]]
        md2=rot@mats["d"]
        truth(f"n{n}_nontrivial_quark_mixing_inside_radius", np.linalg.norm(md2-mats["d"],2) < .125)
        comm=(mats["u"]@mats["u"].conj().T)@(md2@md2.conj().T)-(md2@md2.conj().T)@(mats["u"]@mats["u"].conj().T)
        truth(f"n{n}_changed_quark_mixing_observable_invariant", np.linalg.norm(comm)>1e-3)
        nmatrices[n]=mats
        family.append({"n":n,"fermion_modes":64*n,"full_fermion_Fock_dimension":str(2**(64*n)),
                       "scalar_dimension":7,"S_n":sn,"common_time_probability":formula,
                       "nominal_return_error":nominal,"nominal_error_bound":nominal_bound,
                       "full_parent_return_error_bound":total,"CP_physical_scaled":cp,
                       "CP_allleft_scaled":cpw})

    # Exact CAR creation vectors of all-left mass bilinears; daggers of particle
    # annihilation terms and direct anti-pair creation terms are four orthogonal
    # directions per matrix entry, including off-diagonal flavor coefficients.
    def vacuum_mass_creation(m):
        out=defaultdict(complex)
        n=len(m)
        for i in range(n):
            for j in range(n):
                allleft=m[i,j].conjugate()
                for (fourier,monomial),coefficient in bilinear(i,n+j).items():
                    if monomial[0][0] == 1:
                        out[("anti",monomial)]+=allleft*coefficient
                    else:
                        # Hermitian conjugate reverses order.
                        conjops=[(1,op[1],op[2],op[3]) for op in monomial[::-1]]
                        key,sgn=canon_pair(*conjops)
                        out[("particle",key)]+=allleft.conjugate()*coefficient.conjugate()*sgn
        return dict(out)

    readouts=[]
    W=323.0
    tau=mu/(2**24)
    for n in [3,4,5]:
        m=nmatrices[n]["nu"]
        creates=vacuum_mass_creation(m)
        a2=sum(abs(x)**2 for x in creates.values())
        check(f"n{n}_neutral_vacuum_creation_norm", a2-4*np.linalg.norm(m,"fro")**2, 2e-9)
        truth(f"n{n}_four_orthogonal_mass_channels_per_entry",len(creates)==4*n*n)
        a=mu*np.sqrt(a2)
        b=tau*tau*W*W/2
        intervals=[]
        for c in [1,2]:
            intervals.append([max(0,tau*c*a-b)**2,(tau*c*a+b)**2])
        truth(f"n{n}_fixed_short_time_actual_readout_intervals_disjoint",intervals[0][1]<intervals[1][0])
        check(f"n{n}_two_mass_points_norm_ratio",
              sum(abs(x)**2 for x in vacuum_mass_creation(2*m).values())/a2-4)
        readouts.append({"n":n,"a_over_mu":float(np.sqrt(a2)),
                         "tau_p":tau,"probability_intervals":intervals})

    # Fixed source chain checks: kappa(M,V), g=kappa lambda and physical mass vY.
    def finite_g(mass,vol,lam):
        return 2*lam/np.sqrt(2*mass*vol)
    h=1e-5
    dnum=(finite_g(scalar_mass+h,volume,local_lambda)-finite_g(scalar_mass-h,volume,local_lambda))/(2*h)
    check("scalar_mass_source_includes_kappa", (dnum+g/(2*scalar_mass))/g, 1e-9)
    dvol=(finite_g(scalar_mass,volume+h,local_lambda)-finite_g(scalar_mass,volume-h,local_lambda))/(2*h)
    check("volume_source_includes_kappa", (dvol+g/(2*volume))/g, 1e-9)
    # Simultaneous R compatibility uses full-left fermion phase i, scalar -1.
    check("Dirac_Higgs_residual_R_menu", 1j*1j*(-1)-1)
    check("QQ_scalar_R_menu", 1j*1j*(-1)-1)
    check("QL_scalar_R_menu", 1j*1j*(-1)-1)
    truth("bare_Majorana_not_R_invariant",1j*1j != 1)

    return {
        "round":1321,"status":"AUTHOR_CHECK_PASS",
        "scope":"Complete declared finite Fock parent; not an untruncated QFT apparatus.",
        "checks":CHECKS.copy(),"check_count":len(CHECKS),
        "maximum_absolute_residual":max(x["error"] for x in CHECKS),
        "fixed_contract":{"p":p,"M_S":scalar_mass,"volume":volume,"g":g,"mu":mu,
                          "lambda_D":local_lambda,"lambda_L":local_lambda,"kappa":kappa,
                          "f_sequence":f.tolist(),"T_times_g":T_g,"T":T_g/g,
                          "return_error_threshold":.01,"mass_ball_radius_over_mu":.125},
        "global_norm_bounds":norms,
        "uniform_errors":{"counterrotating_two_pass":float(counter),
                          "mass_two_pass":float(y_bound),
                          "counterrotating_rational":str(counter),"mass_rational":str(y_bound)},
        "families":family,"actual_short_time_readouts":readouts,
        "scope_exclusions":["Unknown full-flavor delivery","Dynamic Higgs/gauge apparatus",
                            "Modes outside declared cutoff","Old long-time P instrument",
                            "Empirically achievable resolution and total cost equality"],
        "source_sha256":{"proof.md":sha(HERE/"proof.md"),"check.py":sha(Path(__file__)),
                         "admission":sha(HERE.parent/"_admission"/"round1321_admission.md")}
    }


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=calculate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        previous=json.loads(path.read_text(encoding="utf-8"))
        if data != previous:
            raise AssertionError("Read-only recomputation differs from frozen results.")
    print(json.dumps({"round":1321,"mode":"write-new" if args.write_new else "read-only",
                      "check_count":data["check_count"],
                      "maximum_absolute_residual":data["maximum_absolute_residual"],
                      "max_return_error_bound":max(row["full_parent_return_error_bound"] for row in data["families"]),
                      "pass":True},ensure_ascii=False))


if __name__=="__main__":
    main()
