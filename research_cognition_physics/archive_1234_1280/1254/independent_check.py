"""Independent round 1254 audit: exact constants and physical-R preparation.
No imports or calls to author/legacy compute functions. Only --write-once writes.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "independent_results.json"
AUTHOR = {
    "admission.md": HERE / "admission.md",
    "proof.md": HERE / "proof.md",
    "check.py": HERE / "check.py",
    "results.json": HERE / "results.json",
    "NEXT.md": HERE / "NEXT.md",
    "../research_note_1254.md": HERE.parent / "research_note_1254.md",
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rat(value):
    value = F(value)
    return {"numerator": value.numerator, "denominator": value.denominator}

def norm_check(value, label, limit=3e-10):
    norm = float(np.linalg.norm(value))
    assert norm < limit, (label, norm)
    return norm

def kron(*items):
    result = np.array([1], dtype=complex)
    for item in items:
        result = np.kron(result, item)
    return result

def classical_group():
    # Rebuild the force tube and PG clock bound from their ingredients.
    n_min, mass_upper = F(320), F(32)
    assert mass_upper / n_min == F(1,10)
    assert F(19) < F(9,2)**2  # m_max = 23+2 sqrt(19) < 32.
    x0, xmin = F(10), F(99,10)
    low_max, high_min, high_max = F(11,10), F(2), F(21,10)
    bracket = 1 + high_max / 2 * (1/xmin - 1/x0)
    bracket_cover = F(501,500)
    assert bracket == F(6607,6600) < bracket_cover
    force_cover = high_max * bracket_cover / (2*xmin*xmin)
    v_cover = F(11,1000)
    assert force_cover == F(1169,108900) < v_cover
    # A first exit before s=1 is impossible; no ODE solver is used here.
    tube_margin = x0 - force_cover/2 - xmin
    assert tube_margin > 0
    low_drop = low_max * bracket_cover / (4*xmin*xmin)
    high_drop = high_min / (4*x0*x0)
    assert low_drop == F(167,59400) and high_drop == F(1,200)

    def energy(b):
        return b - b*b/(2*x0)

    # E'(B)=1-B/x0 >0 throughout the closed source range.
    assert 1-high_max/x0 > 0
    u2_cover = F(2,5)
    u2_actual_cover = 2*energy(high_max)/xmin
    assert u2_actual_cover < u2_cover
    # v<=0 drops uv from the positive PG numerator; then square.
    clock2 = 1/(1-u2_cover) + v_cover*v_cover/(1-u2_cover)**2
    assert clock2 == F(600121,360000) < F(13,10)**2
    proper_per_boundary_min = F(10,13)
    radius_gap = high_drop*proper_per_boundary_min**2-low_drop
    assert radius_gap == F(1477,10038600) > F(1,7000)
    # Strict history inequality holds for t>0; both sides vanish at t=0.
    assert radius_gap * F(0)**2 == 0
    area_gap = 2*xmin*radius_gap
    assert area_gap == F(1477,507000)

    elo, ehi = energy(low_max), energy(high_min)
    energy_gap = ehi-elo
    boundary_radius = F(20)
    compactness_gap = 2*energy_gap/boundary_radius
    curvature_gap = 48*(ehi*ehi-elo*elo)/boundary_radius**6
    assert (elo,ehi,energy_gap)==(F(2079,2000),F(9,5),F(1521,2000))
    assert compactness_gap==F(1521,20000)
    assert curvature_gap==F(25913277,16000000000000)

    # Rational points diagnose the static dictionary independently.
    static_checks = []
    for u, lapse in [(F(3,5),F(4,5)),(F(5,13),F(12,13)),(F(7,25),F(24,25))]:
        assert u*u+lapse*lapse==1
        w=1-lapse; mu=x0*w
        cosh, sinh = 1/lapse, u/lapse
        unreduced_p, normal_energy = mu*u/lapse, mu/lapse
        h = mu-mu*mu/(2*x0)
        assert 2*h/x0==u*u
        assert 2*(1-cosh+u*sinh)-u*u==w*w
        assert sinh-u*cosh==0
        assert w/(cosh-1)==lapse
        assert h==normal_energy+mu*mu/(2*x0)-unreduced_p*u
        assert normal_energy**2-unreduced_p**2==mu**2
        static_checks.append({"u":rat(u),"proper_per_boundary":rat(lapse),
                              "mass":rat(mu),"energy":rat(h)})

    author=json.loads(AUTHOR["results.json"].read_text("utf-8-sig"))
    checks={
        "qmax":bracket,"amax":force_cover,"low_drop":low_drop,
        "clock2":clock2,"radius_gap":radius_gap,"area_gap":area_gap,
    }
    for key,value in checks.items():
        assert author["exact_classical"][key]==rat(value)
    for key,value in {"energy_gap":energy_gap,"compactness_gap":compactness_gap,
                      "scaled_K_gap":curvature_gap}.items():
        assert author["exact_independent_geometry"][key]==rat(value)
    return {
        "force_bracket":rat(bracket),"force_upper":rat(force_cover),
        "first_exit_margin":rat(tube_margin),"proper_speed_upper":rat(v_cover),
        "low_drop_coefficient":rat(low_drop),"high_drop_coefficient":rat(high_drop),
        "u_squared_upper_from_energy":rat(u2_actual_cover),
        "PG_clock_squared_upper_derived":rat(clock2),
        "radius_gap_coefficient":rat(radius_gap),"normalized_area_gap":rat(area_gap),
        "energy_gap":rat(energy_gap),"compactness_gap":rat(compactness_gap),
        "scaled_K_gap":rat(curvature_gap),"static_dictionary_diagnostics":static_checks,
        "strict_radius_inequality_domain":"0 < t <= 1; equality at t=0",
        "scope":"Exact rational constants for a classical finite tube; no numerical ODE or quantum geometry.",
    }

def preparation_group():
    i2,i3=np.eye(2),np.eye(3)
    vac=i3[:,0]; zero=i2[:,0]; one=i2[:,1]
    swap=np.zeros((9,9))
    for l in range(3):
        for r in range(3):
            swap[3*r+l,3*l+r]=1
    parity=[(np.eye(9)+swap)/2,(np.eye(9)-swap)/2]
    hc=1.5*np.eye(4)
    for k in range(3):
        hc[k+1,k]=hc[k,k+1]=math.sqrt((k+1)*(3-k))/2
    energies,hvectors=np.linalg.eigh(hc)
    spectrum_res=norm_check(energies-np.arange(4),"original FK energies")
    weight=np.abs(hvectors[0,:])**2
    weight_res=norm_check(weight-np.array([1,3,3,1])/8,"original J energy weights")
    mstar=17+2*math.sqrt(19)
    age=math.pi/4
    phi=hvectors@(np.exp(-1j*age*(mstar+energies))*(hvectors.conj().T@np.eye(4)[:,0]))
    code=np.column_stack([kron(phi,i3[:,q],vac) for q in [1,2]])

    # Rebuild the original gates. Dressing acts on D only, leaving true R.
    pauli_y=np.array([[0,-1j],[1j,0]])
    pauli_x=np.array([[0,1],[1,0]])
    b=-pauli_y
    ap,am=(i2+b)/2,(i2-b)/2
    ur=kron(ap,i2)-1j*kron(am,pauli_y)
    um=kron(np.outer(zero,zero),i2)-1j*kron(np.outer(one,one),pauli_x)
    gates=[kron(np.eye(8),swap),kron(ur,i2,i3,i3),kron(i2,um,i3,i3)]
    prefixes=[np.eye(72,dtype=complex)]
    for gate in gates:
        prefixes.append(gate@prefixes[-1])
    physical_r_effect=kron(np.eye(72),np.diag([0,1,1]))
    dressing_checks=[]
    for prefix in prefixes:
        dressed=kron(prefix,i3)
        norm_check(dressed.conj().T@dressed-np.eye(216),"prefix unitary")
        dressing_checks.append(norm_check(
            dressed.conj().T@physical_r_effect@dressed-physical_r_effect,
            "true physical R preserved by original dressing"))

    projectors=[]
    for n in range(6):
        p=np.zeros((36,36),complex)
        for k in range(4):
            for side in range(2):
                if k+2*side==n:
                    p+=kron(np.outer(hvectors[:,k],hvectors[:,k].conj()),parity[side])
        projectors.append(p)
    norm_check(sum(projectors)-np.eye(36),"six source projectors complete")
    for j,p in enumerate(projectors):
        for k,q in enumerate(projectors):
            norm_check(p@q-(p if j==k else 0),"source projector orthogonality")
    pieces=[p@code for p in projectors]
    r_effect=kron(np.eye(12),np.diag([0,1,1]))
    source_res=[]
    for j,number in enumerate([1,3,4,4,3,1]):
        source_res.append(norm_check(pieces[j].conj().T@pieces[j]-number/16*np.eye(2),
                                     "full-input source weight"))

    # Purely rational lower bounds, computed from b/x0 rather than supplied c^2.
    lower_bounds={}
    for source in [1,2]:
        w=F(source,10)
        d2=(2*w-w*w)/(1-w)**2
        lower=d2/(2*(1+d2))
        expected={1:(F(19,81),F(19,200)),2:(F(9,16),F(9,50))}[source]
        assert (d2,lower)==expected
        assert F(source)+F(32,320)<10
        lower_bounds[source]={"derivative_squared_lower":rat(d2),
                              "probability_lower_hbar1":rat(lower)}

    rows=[]
    for source in [1,2]:
        for resource in [320,1280,5120,81920]:
            momenta=[]
            for n in range(6):
                lam=source+(mstar+1+n)/resource
                w=lam/10
                u=math.sqrt(2*w-w*w)
                momenta.append(10*(math.atanh(u)-u))
            p=np.array(momenta)
            gram_abs=np.exp(-(resource*(p[:,None]-p[None,:]))**2/4)
            eta_min=sum(weight[k]*(1-gram_abs[k,k+2])/2 for k in range(4))
            assert eta_min>float({1:F(19,200),2:F(9,50)}[source])
            for phase_name,phase in [
                ("positive_overlap",np.zeros(6)),
                ("nonzero_spectral_phases",np.array([0,.2,-.4,.7,-.3,.9])),
            ]:
                # gram[j,k]=<phi_j|phi_k>. Keep all six-sector coherences.
                gram=gram_abs*np.exp(1j*(phase[None,:]-phase[:,None]))
                assert np.min(np.linalg.eigvalsh(gram))>-3e-11
                eta=sum(weight[k]*(1-gram[k,k+2].real)/2 for k in range(4))
                assert eta>=eta_min-3e-13
                effect=sum(gram[j,k]*(pieces[j].conj().T@r_effect@pieces[k])
                           for j in range(6) for k in range(6))
                scalar_res=norm_check(effect-eta*np.eye(2),"all-input true-R effect")
                isometry=sum(gram[j,k]*(pieces[j].conj().T@pieces[k])
                             for j in range(6) for k in range(6))
                norm_check(isometry-np.eye(2),"conditioned preparation isometry")
                choi=np.zeros((3,2,3,2),complex)
                matrix_res=[]
                for a in range(2):
                    for c in range(2):
                        unit=np.zeros((2,2)); unit[a,c]=1
                        full=sum(gram[k,j]*(pieces[j]@unit@pieces[k].conj().T)
                                 for j in range(6) for k in range(6))
                        reduced=np.trace(full.reshape(12,3,12,3),axis1=0,axis2=2)
                        target=np.zeros((3,3),complex)
                        target[1:,1:]=eta*unit
                        if a==c: target[0,0]=1-eta
                        matrix_res.append(norm_check(reduced-target,"true-R matrix unit channel"))
                        choi[:,a,:,c]=reduced/2
                choi=choi.reshape(6,6)
                blank=kron(np.outer(vac,vac),np.eye(2)/2)
                ideal=np.zeros(6,complex)
                ideal[2]=ideal[5]=1/math.sqrt(2)
                expected=eta*np.outer(ideal,ideal.conj())+(1-eta)*blank
                choi_res=norm_check(choi-expected,"passive-reference channel")
                ref_distance=float(np.sum(abs(np.linalg.eigvalsh(choi-blank)))/2)
                norm_check(np.array([ref_distance-eta]),"reference distance from original blank")
                rows.append({
                    "b":source,"N":resource,"phase_case":phase_name,
                    "P_R_nonvac":float(eta),"positive_overlap_probability":float(eta_min),
                    "all_input_effect_residual":scalar_res,
                    "matrix_unit_channel_residual":max(matrix_res),
                    "passive_reference_residual":choi_res,
                    "reference_half_trace_distance_from_blank":ref_distance,
                })
    return {
        "original_FK_spectrum_residual":spectrum_res,"original_J_weights_residual":weight_res,
        "physical_R_dressing_residuals":dressing_checks,"source_weight_residuals":source_res,
        "exact_uniform_bounds":{str(k):v for k,v in lower_bounds.items()},
        "matrix_diagnostics":rows,
        "scope":"Standard coherent packets on an explicitly extended real canonical chart; centers only are static. Physical R at preparation, not later transfer.",
        "not_claimed":["Exact zero-velocity quantum state","Half-line packet result",
                       "All squeezed packets excluded","Certified quantum geometry or actual geometric readout"],
    }

def calculate():
    return {
        "round":1254,"status":"PASS","independent_groups":2,
        "author_sha256":{key:sha(path) for key,path in AUTHOR.items()},
        "independent_code_sha256":sha(Path(__file__)),
        "group1_classical_geometry":classical_group(),
        "group2_quantum_preparation":preparation_group(),
    }

def compare(saved,current):
    if isinstance(current,dict):
        assert saved.keys()==current.keys(),(saved.keys(),current.keys())
        for key in current: compare(saved[key],current[key])
    elif isinstance(current,list):
        assert len(saved)==len(current)
        for a,b in zip(saved,current): compare(a,b)
    elif isinstance(current,float):
        assert math.isclose(saved,current,rel_tol=2e-9,abs_tol=3e-10),(saved,current)
    else:
        assert saved==current,(saved,current)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-once",action="store_true")
    args=parser.parse_args()
    result=calculate()
    if args.write_once:
        with OUTPUT.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        print("FIRST_WRITE_PASS_1254_INDEPENDENT: 2 groups")
    else:
        compare(json.loads(OUTPUT.read_text("utf-8-sig")),result)
        print("READ_ONLY_MATCH_1254_INDEPENDENT: 2 groups")
