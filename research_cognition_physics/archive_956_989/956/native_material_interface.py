"""956: an exact electronic dimer interface for old relational exchange tasks.

Hubbard physics is standard. This audit transports the SAME finite material
Hamiltonian, an actual charge effect and parameter sources to the spin sector.
It does not certify a full quantum-dot array, autonomous readout, QED or GR.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"native_material_interface_results.json"
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text("utf-8-sig"))
def norm(a):return float(np.linalg.norm(a,2))
def evolution(h,time):
    e,v=np.linalg.eigh(h);return (v*np.exp(-1j*time*e))@v.conj().T
def annihilate(mode):
    a=np.zeros((16,16))
    for state in range(16):
        if (state>>mode)&1:
            a[state^(1<<mode),state]=(-1)**((state&((1<<mode)-1)).bit_count())
    return a
def material(U,v):
    cs=[annihilate(i) for i in range(4)]
    ns=[c.T@c for c in cs]
    hopping=sum(cs[i].T@cs[i+2]+cs[i+2].T@cs[i] for i in (0,1))
    D=ns[0]@ns[1]+ns[2]@ns[3]
    sector=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    reduce=lambda a:sector.T@a@sector
    H=reduce(-v*hopping+U*D);D=reduce(D);hop=reduce(hopping)
    # Spin order up-up, up-down, down-up, down-down.
    single=sector.T@np.eye(16)[:,[5,9,6,10]]
    swap=np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]],float)
    singlet=np.array([0,1,-1,0])/np.sqrt(2)
    ps=np.outer(singlet,singlet)
    s=single@singlet
    # Fix double-occupancy phase from the original CAR hopping.
    d=-H@s/(2*v)
    assert abs(np.linalg.norm(d)-1)<1e-12 and norm(D@d-d)<1e-12
    om=math.sqrt(U*U+16*v*v);J=8*v*v/(om+U)
    occupancy=(1-U/om)/2
    c=math.sqrt(1-occupancy);b=math.sqrt(occupancy)
    W=single+(c-1)*np.outer(s,singlet)+b*np.outer(d,singlet)
    low=-J*ps
    return dict(H=H,D=D,hop=hop,W=W,low=low,ps=ps,swap=swap,
        single=single,s=s,d=d,omega=om,J=J,occupancy=occupancy)

def run():
    U=1.;v=.1;kappa=.7
    m=material(U,v);H=m["H"];D=m["D"];W=m["W"];ps=m["ps"]
    om=m["omega"];J=m["J"];d=m["occupancy"];low=m["low"]
    ev=np.linalg.eigvalsh(H)
    assert norm(W.T@W-np.eye(4))<1e-12
    intertwiner=norm(H@W-W@low)
    assert intertwiner<1e-12
    charge=norm(W.T@D@W-d*ps)
    sourceU=norm(W.T@D@W-((1-U/om)/2)*ps)
    sourcev=norm(W.T@(-m["hop"])@W+8*v/om*ps)
    assert max(charge,sourceU,sourcev)<1e-12
    J2=4*v*v/U;remainder=J2-J;remainder_bound=16*v**4/U**3
    rows=[]
    for angle in (.7,math.pi/2,math.pi):
        time=angle/J2
        exact_error=norm(evolution(H,time)@W-W@evolution(low,time))
        second_error=norm(evolution(H,time)@W-W@evolution(-J2*ps,time))
        bound=time*remainder_bound
        assert exact_error<1e-12 and second_error<=bound+1e-12
        rows.append(dict(exchange_phase=angle,time=time,
            exact_low_sector_error=exact_error,second_order_error=second_error,
            second_order_analytic_upper=bound))
    # Differentiate the SAME material H, and independently its low eigenvalue.
    def energy(a,b):return float(np.linalg.eigvalsh(material(a,b)["H"])[0])
    h=2e-4
    dU=(energy(U+h,v)-energy(U-h,v))/(2*h)
    dv=(energy(U,v+h)-energy(U,v-h))/(2*h)
    dUU=(energy(U+h,v)-2*energy(U,v)+energy(U-h,v))/(h*h)
    response_errors=dict(first_U=abs(dU-d),first_v=abs(dv+8*v/om),
        second_U=abs(dUU+8*v*v/om**3))
    assert response_errors["first_U"]<1e-7 and response_errors["first_v"]<1e-6
    assert response_errors["second_U"]<1e-7
    # Virtual high states, rather than P H'' P=0, supply static curvature.
    val,V=np.linalg.eigh(H);ground=V[:,0];curvature=0.
    for a in range(1,len(val)):
        curvature+=2*abs(V[:,a].conj()@D@ground)**2/(val[0]-val[a])
    assert abs(curvature+8*v*v/om**3)<1e-12
    # No spin reference: charge and H commute with all total SU(2) generators.
    cs=[annihilate(i) for i in range(4)]
    basis=np.eye(16)[:,[s for s in range(16) if s.bit_count()==2]]
    pauli=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
    spin_commutator=0.
    for p in pauli:
        spin=sum(.5*p[a,b]*cs[2*site+a].T@cs[2*site+b]
            for site in (0,1) for a in (0,1) for b in (0,1))
        spin=basis.T@spin@basis
        spin_commutator=max(spin_commutator,norm(D@spin-spin@D),norm(H@spin-spin@H))
    # Preserve complete six-dimensional poststates; charge readout can leave W.
    Q=np.eye(6)-W@W.T
    leakage=W.T@(D@Q@D+(np.eye(6)-D)@Q@(np.eye(6)-D))@W
    leakage_error=norm(leakage-2*d*(1-d)*ps)
    assert leakage_error<1e-12
    effects=(W.T@D@W,W.T@(np.eye(6)-D)@W)
    assert norm(sum(effects)-np.eye(4))<1e-12
    # Transport the OLD four-spin code: pair 12 singlet is logical zero.
    b16=np.eye(16)
    code=np.stack([(b16[5]-b16[6]-b16[9]+b16[10])/2,
        (2*b16[3]+2*b16[12]-b16[5]-b16[6]-b16[9]-b16[10])/np.sqrt(12)],axis=1)
    encoded_charge=code.T@np.kron(d*ps,np.eye(4))@code
    code_error=norm(encoded_charge-np.diag([d,0.]))
    assert code_error<1e-12
    # Actual 434 weights, positive at its frozen witness epsilon=1/512.
    eps=1/512;weights={}
    from itertools import combinations
    for block in ((0,1,2,3),(4,5,6,7)):
        for a,b in combinations(block,2):weights[(a,b)]=1.
    weights[(0,4)]=eps;weights[(1,5)]=eps
    weights[(0,1)]-=eps**2/12;weights[(4,5)]-=eps**2/4
    weights[(5,6)]+=eps**2/(6*np.sqrt(3))
    weights[(4,6)]-=eps**2/(6*np.sqrt(3))
    scale=.001
    hopping={f"{a+1}-{b+1}":math.sqrt(U*scale*w/2) for (a,b),w in weights.items()}
    matching_error=max(abs(2*x*x/U-scale*w) for x,w in zip(hopping.values(),weights.values()))
    assert min(weights.values())>0 and len(weights)==14 and matching_error<1e-14
    # Dropping the identity term changes material response despite same spin channel.
    Jprime=-8*kappa*v*v/om
    actual_force=Jprime*ps
    no_offset_force=-Jprime/2*m["swap"]
    assert norm(actual_force-no_offset_force-Jprime/2*np.eye(4))<1e-12
    files=[Path(__file__),STAGE.parent/"archive_429_466/434/encoded_exchange_response_audit.py",
        STAGE.parent/"archive_429_466/434/encoded_exchange_response_audit_results.json",
        STAGE/"955/joint_effective_window_results.json"]
    return dict(round=956,date="2026-10-07",all_scientific_checks_passed=True,
        parameters=dict(U=U,hopping=v,kappa=kappa,source_contract="static U,v; chosen v(R)=v0 exp[-kappa(R-R0)]"),
        material=dict(two_electron_dimension=6,low_spin_dimension=4,spectrum=ev.tolist(),
            exact_exchange_J=J,second_order_J=J2,
            second_order_remainder=remainder,second_order_remainder_upper=remainder_bound,
            charge_probability_singlet=d,charge_probability_triplet=0.,
            encoded_434_logical_charge_effect=encoded_charge.tolist(),
            unread_charge_instrument_leakage_singlet=2*d*(1-d),
            singlet_radial_force=Jprime,omitted_identity_source_error=norm(actual_force-no_offset_force),
            static_U_curvature=curvature,bare_compressed_U_curvature=0.),
        verification=dict(exact_H_intertwining_error=intertwiner,charge_effect_error=charge,
            U_source_error=sourceU,hopping_source_error=sourcev,SU2_commutator=spin_commutator,
            charge_poststate_leakage_formula_error=leakage_error,
            old_434_logical_readout_transport_error=code_error,finite_difference=response_errors,
            finite_time_rows=rows),
        old_434_match=dict(epsilon=eps,all_14_weights_positive=True,spin_energy_scale=scale,
            hopping_coefficients=hopping,max_hopping_over_U=max(hopping.values())/U,
            second_order_coefficient_error=matching_error,
            scalar_offset_must_be_kept=-scale*sum(weights.values()),
            full_array_high_order_and_long_time_transport_certified=False),
        scope=dict(native_electronic_exchange_charge_and_static_sources_common=True,
            old_exchange_universality_reused_not_reproved=True,
            elementary_new_neutral_portal_required_for_this_interface=False,
            original_955_bounds_transferred=False,
            spin_only_postmeasurement_closure_claimed=False,
            all_six_protocols_autonomously_implemented=False,
            material_to_full_SM_and_GR_error_certified=False,
            three_spatial_dimensions_derived=False,full_goal_completed=False),
        references=["https://arxiv.org/abs/cond-mat/9701055","https://arxiv.org/abs/1502.02194",
            "https://arxiv.org/abs/1105.0675"],
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    if a.write:assert not TARGET.exists()
    out=run()
    if a.write:
        with TARGET.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        old=read(TARGET)
        assert old["source_hashes"]==out["source_hashes"] and old["scope"]==out["scope"]
        for key in ("exact_exchange_J","charge_probability_singlet","static_U_curvature"):
            assert abs(old["material"][key]-out["material"][key])<1e-12
    print(json.dumps({k:v for k,v in out.items() if k!="source_hashes"},ensure_ascii=False,indent=2))
