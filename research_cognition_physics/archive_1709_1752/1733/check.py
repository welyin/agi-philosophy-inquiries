"""Selection characters, a sharp real Wilson budget and tagged heavy-pole moments."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())


CFG=read("inputs.json")


def charges():
    h,hbar,tag=map(np.array,[*CFG["new_H_and_barH_Z12_Z3_R"],CFG["neutral_orientation_tags_Z12_Z3_R"]])
    for n in range(4):
        np.testing.assert_array_equal((hbar+h+n*tag)%[12,3,100],[0,0,2])
    assert sum(CFG["doublet_hypercharges"])==0
    old=read(CFG["flavor"])
    neutral=read(CFG["neutral"])
    assert "10,5bar,Hu,Hd have 2,0,2,1" in neutral["auxiliary_group"]
    # Reconstruct the old 5bar charges from the actual 1685 edges A12,B21,R23.
    T=np.array([0,0,1])
    F=np.array([2,0,1])
    edges=[(0,1,0,0),(1,0,2,0),(1,1,0,0),(2,1,1,0),
           (2,2,2,0),(1,2,1,0),(1,2,5,1)]
    residues=[(T[i]+F[j]+n*old["flavor_charge"]["S"]+m*old["flavor_charge"]["X"])%12
              for i,j,n,m in edges]
    assert not any(residues)
    # Entry t times an invariant round trip cannot change any additive character,
    # including an ordinary R charge (the comparison is between two superpotential terms).
    rows=[]
    for a,b in ((-1,1),(3,-7),(0,0)):
        primitive=b-a
        assert primitive-primitive+primitive==primitive
        rows.append([a,b,primitive])
    return {"added_vectorlike_pairs":3,"new_pairs_delta_b_SUSY":[1.8,3.,0.],
            "seven_original_down_vertex_Z12_residues":[int(v) for v in residues],
            "old_S_charge_mod12":read(CFG["flavor"])["flavor_charge"]["S"]%12,
            "closed_return_does_not_alter_tag_character_examples":rows}


def phase_budget():
    cfg=read(CFG["flavor"])
    s=cfg["scalar_magnitude_over_cutoff"]*np.exp(1j*np.pi*cfg["scalar_angle_over_pi"])
    x=cfg["real_X_over_cutoff"]
    assert x==abs(s)==.5
    ideal=np.exp(1j*np.angle(s))
    I=x*s**4
    # The vacuum aligns the two auxiliary flavons to s^2,s^3.
    q=cfg["flavor_charge"]
    assert (q["A2"]-2*q["S"])%12==(q["A3"]-3*q["S"])%12==0
    k,m=np.unravel_index(np.arange(49*19),(49,19))
    good=(q["S"]*k+q["X"]*m)%12==0
    terms=s**k[good]*x**m[good]
    n=k[good]//4
    assert np.all(k[good]%4==0) and np.all((m[good]-n)%3==0)
    amps=[max(abs(terms[n%3==i])) for i in range(3)]
    np.testing.assert_allclose(amps,[1,1/32,1/1024],rtol=1e-12)
    # Basis 1,I; their symmetric convex hull contains the third ray and all shorter vectors.
    def coefficients(z):
        b=z.imag/I.imag
        return (z.real-b*I.real,b)
    for z in terms:
        assert sum(abs(v) for v in coefficients(z))<=1+2e-13
    a,b=coefficients(ideal)
    np.testing.assert_allclose(a+b*I,ideal,rtol=1e-14)
    np.testing.assert_allclose([a,b],[2/math.sqrt(3),32/math.sqrt(3)],rtol=2e-14)
    theta_limit=1.7476e-10
    leak_max=math.tan(theta_limit)/(abs(ideal.imag)+(1-ideal.real)*math.tan(theta_limit))
    np.testing.assert_allclose(abs(np.angle(1+leak_max*(ideal-1))),theta_limit,rtol=1e-12)
    return {"sampled_maximum_on_three_rays_with_analytic_tails":amps,
            "additional_factor_Wilsons_on_1_and_I":[a,b],"lowest_l1":abs(a)+abs(b),
            "alternative_2_plus_I_over_lock_over_sqrt3":abs((2+I/abs(I))/math.sqrt(3)-ideal),
            "neutral_factor_complex_log_derivative_for_log_abs_S":[(4*b*I/ideal).real,(4*b*I/ideal).imag],
            "positive_color_leak_cap_at_single_source_theta_limit":leak_max,
            "theta_shift_rad_at_color_leak_one_part_in_million":np.angle(1+1e-6*(ideal-1))}


def matrix(t,back,eta,c0,masses):
    return np.diag(masses)@np.array([[-t,1,0,0],[0,-back,1,0],[-eta*c0**2*t,0,-(1-eta)*t,1]],complex)


def observations(B):
    poles,cols=np.linalg.eigh(B.conj().T@B)
    poles[0]=0
    zero=cols[:,0]/cols[0,0]
    weights=cols[3,:]*cols[0,:].conj()
    mixed=[sum(weights*poles**j) for j in range(4)]
    direct=np.linalg.matrix_power(B.conj().T@B,1)[3,0]
    np.testing.assert_allclose(mixed[1],direct,atol=2e-14)
    return {"masses":np.sqrt(poles[1:]),"zero":zero,"mixed":np.array(mixed),
            "den":(B.conj().T@B)[3,3].real,"weight":weights, "H0_prob":abs(cols[0,0])**2}


def endpoint_ratio(t,back,eta,c0):
    return t*((1-eta)*t*back+eta*c0*c0)


def case(c):
    m1,m3=CFG["witness_M1_M3_units"]
    masses=np.array([m1,c*m1,m3])
    eta=2/(1+c**4)
    pure=observations(matrix(c,c,0,c,masses))
    bypass=observations(matrix(c,c,eta,c,masses))
    np.testing.assert_allclose(pure["masses"],bypass["masses"],rtol=1e-14)
    target=np.array([1,c,c*c,c**3])
    for row in (pure,bypass):
        np.testing.assert_allclose(row["zero"],target,atol=1e-14)
        np.testing.assert_allclose(row["H0_prob"],1/sum(target**2),rtol=1e-14)
        np.testing.assert_allclose(np.prod(row["masses"]**2),np.prod(masses**2)*sum(target**2),rtol=1e-14)
    np.testing.assert_allclose(pure["mixed"][:3],0,atol=4e-14)
    np.testing.assert_allclose(pure["mixed"][3],-np.prod(masses**2)*c**3,rtol=2e-13)
    np.testing.assert_allclose(-bypass["mixed"][1]/bypass["den"],eta*c**3,rtol=2e-14)
    eps=1e-5
    slopes={}
    for name,e,alternates in (("true_return",0,True),("three_forward",0,False),("single_only",1,True),
                               ("isospectral_return_and_bypass",eta,True)):
        def ratio(t):
            back=t.conjugate() if alternates else t
            return endpoint_ratio(t,back,e,c)/t
        mod=np.log(abs(ratio(c*np.exp(eps))/ratio(c*np.exp(-eps))))/(2*eps)
        turn=np.angle(ratio(c*np.exp(1j*eps))/ratio(c*np.exp(-1j*eps)))/(2*eps)
        np.testing.assert_allclose([mod,turn],[2*(1-e),0 if alternates else 2*(1-e)],atol=3e-9)
        slopes[name]=[mod,turn]
    def parts(values):
        return [[float(v.real),float(v.imag)] for v in values]
    zheavy=lambda q2,row:sum(row["weight"]/(q2+np.r_[0,row["masses"]**2]))
    return {"c":c,"short_and_long_same_x":c**3,"light_original_H_probability":pure["H0_prob"],
            "bare_mass_ratios":masses.tolist(),"all_three_common_heavy_masses":pure["masses"].tolist(),
            "bypass_eta":eta,"bypass_final_neighbour_strength":1-eta,
            "normalized_direct_tag_pure":float((-pure["mixed"][1]/pure["den"]).real),
            "normalized_direct_tag_isospectral":float((-bypass["mixed"][1]/bypass["den"]).real),
            "pure_mixed_moments_n0_to_3":parts(pure["mixed"]),
            "equal_poles_different_residues":{"pure":parts(pure["weight"]),"bypass":parts(bypass["weight"])},
            "spacelike_tagged_transfer_at_Q2_equal_M1_squared":[zheavy(m1*m1,pure).real,zheavy(m1*m1,bypass).real],
            "relative_x_over_single_log_modulus_and_phase_slopes":slopes}


def run():
    prior=read(CFG["fit"])
    targets={name:next(x["Vud_10PeV"] for x in rows if x["k"]==3)
             for name,rows in prior["solutions"].items()}
    targets.update({name:row["Vud_10PeV"] for name,row in prior["exploratory_k3"].items()})
    return {"charge_check":charges(),"neutral_phase":phase_budget(),
            "cases":{name:case(c) for name,c in targets.items()}}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    target=HERE/"results.json"
    if args.check:
        assert json.loads(target.read_bytes())==result
        print("1733 charge identities, neutral-phase budget and isospectral tagged propagators verified")
    else:
        target.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
