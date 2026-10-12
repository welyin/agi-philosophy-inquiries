"""1134: induced positive-mass Poincare model with the same LLP-entangled task.

Default is read-only. --write-results writes only 1134/results.json.
Finite fibers and numerical Lorentz checks are diagnostics, not continuum proofs.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
FROZEN = {
    "archive_1120_/1129/proof.md": "c835b90f123547f7d224c72e1e0cbf89e1dd2fa377ae67ce90b1211922fcfefa",
    "archive_1120_/1129/results.json": "69ed793bd5847d31003eb05b18205ecf72af34c26f9ad2919ebe548e7a25a02c",
    "archive_1120_/1132/proof.md": "f05daf49972a0b7fe13f5d593a032998bf1c9008182c051f2e12bb6ea958df63",
    "archive_1120_/1132/results.json": "95dcbd80a2a9e5d3d7ba03e52f1d1aa7e931606c9735f0359623f7b4f969430c",
    "archive_935_955/research_note_941.md": "a12e2dd67212a77538f850908efc6057e450d4def00dccd59aaa7cfaea346b7f",
    "archive_1120_/1127/proof.md": "d1c8827cc12563b996c951da6039e041e08aa711de18e84a56bb240660660c27",
}


def norm(a):
    return float(np.linalg.norm(a, 2))


def exp_h(h, t):
    e, v = np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T


def sources():
    hashes = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert hashes == FROZEN
    old = json.loads((ROOT/"archive_1120_/1132/results.json").read_text(encoding="utf-8"))
    assert old["source_identity"]["complete_old_clock_dimension"] == 288
    assert old["exact_budgets"]["fractions"]["mass"] == str(10**80)
    return dict(source_sha256=hashes, complete_old_clock_dimension=288,
        source_and_support_preparations_changed=False,
        LLP_state_replaced_by_independent_packet=False,
        H_BT_is_new_dynamics=True, receiver_1133_imported=False,
        primary_full_texts_read=[
            dict(url="https://arxiv.org/html/1008.5215", sections="II, IV, V; equations 3, 11-14; strong versus weak cluster"),
            dict(url="https://arxiv.org/html/1008.5222", sections="I, equations 2-6; positive mass and rotation decomposition")],
        permission="M1134: new induced interacting Poincare candidate; complete experimental transport may be declared, but is not inherited old CO/boost endpoint identity.")


def budgets():
    M, A, sigma, t = F(10**80), F(10**34), F(1,10**30), F(5,4)
    C, lam, eps = F(33536000000000000000), F(11,7), F(1,10**5)
    v, Ceta = 16*lam*C, F(24780800000000000000000,49)
    mmin = M-v
    old = json.loads((ROOT/"archive_1120_/1132/results.json").read_text(encoding="utf-8"))["exact_budgets"]["fractions"]
    old_move = F(old["complete_moving_error"])
    old_handoff, old_causal = F(old["joint_handoff_error"]), F(old["finite_causal_error"])
    assert mmin > M/2
    assert F(1157516)+F(5,3)*8836**2 < 12000**2
    rows = {}
    for label,b1,b2 in [("code",268,28000),("finite_menu",600,400000)]:
        p1, p2 = b1+94*A,b2+188*b1*A+12000*A*A
        a1, a2, b2mix = 1/sigma+p1,1/sigma**2+2*p1/sigma+p2,p1/sigma+p2
        term1 = t*a2/(2*mmin)
        term2 = t*(a2+2*b2mix)/(2*M)
        term3 = Ceta*a1*t*t/(2*M)
        delta = term1+term2+term3
        assert delta < F(3000004,10**14)
        vals = dict(material_first_moment_norm=F(b1),material_second_moment_norm=F(b2),
            field_first_moment_norm=p1,field_second_moment_norm=p2,
            A1=a1,A2=a2,B2=b2mix,mass_shell_term=term1,
            original_fiber_term=term2,actual_original_force_term=term3,complete_BT_error=delta)
        rows[label] = dict(fractions={k:str(x) for k,x in vals.items()}, decimal_views={k:float(x) for k,x in vals.items()})
    code = F(rows["code"]["fractions"]["complete_BT_error"])
    menu = F(rows["finite_menu"]["fractions"]["complete_BT_error"])
    theta = F(11,7)
    dw1 = theta*4*40
    dw2 = theta*(4*7200+4*40**2)+theta**2*16*40**2
    assert dw1 < 252 and dw2 < 120000
    assert 268+252 < 600 and 28000+2*252*268+120000 < 400000
    assert old_handoff+code < F(1258186,10**9)
    assert old_causal+2*menu < F(120170,10**9)
    vals = dict(mass=M,amplitude=A,sigma=sigma,time=t,V_norm_upper=v,
        internal_mass_lower=mmin,source_cumulant_m2_upper=8836*A*A,
        source_cumulant_m4_upper=1157516*A*A,
        source_first_momentum_norm=94*A,source_second_momentum_norm=12000*A*A,
        finite_gate_first_derivative=dw1,finite_gate_second_derivative=dw2,
        gate_code_first_moment_upper=F(520),gate_code_second_moment_upper=F(283072),
        complete_handoff_error=old_handoff+code,finite_causal_error=old_causal+2*menu,
        old_clock_effect_error=old_move+code,
        order_one_boost_mean_K_lower=F(3,4)*mmin-F(5,4)*F(rows["code"]["fractions"]["A1"]))
    return dict(scope="Rest preparation and explicitly bounded two-branch menu only. All exact fractions preserve the original LLP-entangled preparation.",
        fractions={k:str(x) for k,x in vals.items()},decimal_views={k:float(x) for k,x in vals.items()},
        preparation_rows=rows,old_clock_reference_probability=float(np.sin((np.pi/2+1.25)/2)**6),
        arbitrary_normal_input_bound_claimed=False,nonzero_boost_inherits_rest_error=False,
        old_separate_force_formula_inherited_by_BT=False,spatial_translation_and_rotation_same=True)


def quad(fn,a,b,n=160):
    x,w = np.polynomial.legendre.leggauss(n)
    return float((b-a)/2*np.sum(w*fn((a+b)/2+(b-a)*x/2)))


def coherent_moments():
    # Isotropic six-mode coherent momentum distribution. Independent Poisson
    # occupation sums verify the fourth cumulant, rather than assuming Gaussian.
    rows = []
    ns = np.arange(38)
    diff = ns[:,None]-ns[None,:]
    for mu in [.2,1.,2.]:
        poisson = np.array([np.exp(-mu)*mu**int(n)/math.factorial(int(n)) for n in ns])
        probs = np.outer(poisson,poisson)
        axis2,axis4 = np.sum(probs*diff**2),np.sum(probs*diff**4)
        total2,total4 = 3*axis2,3*axis4+6*axis2**2
        m2,m4 = 6*mu,6*mu
        formula = m4+5*m2*m2/3
        assert abs(total2-m2)<1e-12 and abs(total4-formula)<2e-10
        assert abs(total4-5*m2*m2/3)>.9*m4
        rows.append(dict(per_direction_mean=mu,probability_mass=float(poisson.sum()),
            total_second_moment=float(total2),total_fourth_moment=float(total4),
            cumulant_formula=formula,missing_Gaussian_only_term=m4))
    b = lambda x: np.exp(-1/(x*(1-x)))
    ib = quad(b,0,1)
    rp = lambda x:b(x)/ib
    rpp = lambda x:b(x)*(1-2*x)/(x*x*(1-x)**2)/ib
    # Integrate beta derivative norms on the actual shell; radial R itself is
    # not needed for gradient/laplacian norms, and ||beta|| <= volume is analytic.
    grad2 = quad(lambda u:4*np.pi*(7+u)**2*rp(u)**2,0,1)
    lap2 = quad(lambda u:4*np.pi*(7+u)**2*(rpp(u)+2*rp(u)/(7+u))**2,0,1)
    assert np.sqrt(grad2)<376 and np.sqrt(lap2)<6157
    return dict(scope="Independent finite coherent Poisson occupation sums plus actual beta shell integrals. The continuous Fock moment proof uses its exact characteristic functional.",
        six_mode_rows=rows,actual_beta_gradient_norm=float(np.sqrt(grad2)),
        actual_beta_laplacian_norm=float(np.sqrt(lap2)),
        zero_mean_momentum_not_zero_fluctuation=True)


def kron_all(*a):
    ans = np.array([[1.]])
    for item in a:
        ans = np.kron(ans,item)
    return ans


def llp_fibers():
    # 1D translation subgroup diagnostic: retained field momentum, unknown Q,
    # a source register and complete internal two-level system. External K is
    # integrated continuously. No actual continuum CAR/CCR replacement claimed.
    labels = [(p,q,s,d) for p in [-1,0,1] for q in [0,1] for s in [0,1] for d in [0,1]]
    ps = np.array([x[0] for x in labels],float)
    ds = np.array([1+x[3] for x in labels],float)
    ss = np.array([x[2] for x in labels],float)
    pmat = np.diag(ps)
    sx = np.array([[0,1],[1,0.]])
    hop = np.diag(np.ones(2),1)+np.diag(np.ones(2),-1)
    v = .15*kron_all(hop,sx,sx,np.eye(2))
    hf = np.diag(np.sqrt(1+ps*ps)+.3*ss)
    force_norm = norm(1j*(v@pmat-pmat@v))
    initial = np.zeros((24,2),complex)
    pw = np.array([1,2,1])/np.sqrt(6)
    sw = np.array([np.sqrt(.7),1j*np.sqrt(.3)])
    dw = np.sqrt([.6,.4])
    for i,(p,q,s,d) in enumerate(labels):
        initial[i,q] = pw[p+1]*sw[s]*dw[d]
    assert norm(initial.conj().T@initial-np.eye(2))<1e-13
    sigma,t = .7,.8
    x,w = np.polynomial.legendre.leggauss(180)
    ks,weights = 9*x,9*w
    packet = lambda k:(2*sigma*sigma/np.pi)**.25*np.exp(-sigma*sigma*k*k)
    field_p2 = float(np.sum(pw**2*np.array([-1,0,1])**2))
    field_p4 = field_p2
    variance = 1/(4*sigma*sigma)
    A1 = np.sqrt(variance+field_p2)
    A2 = np.sqrt(3*variance**2+6*variance*field_p2+field_p4)
    B2 = np.sqrt(variance*field_p2+field_p4)
    gram = np.zeros((2,2),complex)
    wrong_gram = np.zeros((2,2),complex)
    internal_rho = np.zeros((24,24),complex)
    for k,weight in zip(ks,weights):
        state = packet(k-ps)[:,None]*initial
        independent = packet(k)*initial
        gram += weight*(state.conj().T@state)
        wrong_gram += weight*((state-independent).conj().T@(state-independent))
        internal_rho += weight*np.outer(state[:,0],state[:,0].conj())
    wrong = float(np.sqrt(np.linalg.eigvalsh(wrong_gram).max()))
    purity = float(np.real(np.trace(internal_rho@internal_rho)))
    assert norm(gram-np.eye(2))<1e-12 and wrong>.2 and purity<.95
    rows = []
    for mass in [8.,40.,200.]:
        mint = hf+np.diag(np.sqrt((mass+ds)**2+ps*ps))+v
        ev,evec = np.linalg.eigh(mint)
        error_gram = np.zeros((2,2),complex)
        output_gram = np.zeros((2,2),complex)
        for k,weight in zip(ks,weights):
            state = packet(k-ps)[:,None]*initial
            ho = hf+np.diag(np.sqrt((mass+ds)**2+(k-ps)**2))+v
            old_out = exp_h(ho,t)@state
            bt_out = (evec*np.exp(-1j*t*np.sqrt(ev*ev+k*k)))@evec.conj().T@state
            difference = bt_out-old_out
            error_gram += weight*(difference.conj().T@difference)
            output_gram += weight*(bt_out.conj().T@bt_out)
        actual = float(np.sqrt(max(0,np.linalg.eigvalsh(error_gram).max())))
        mmin = mass-norm(v)
        bound = t*A2/(2*mmin)+t*(A2+2*B2)/(2*mass)+force_norm*A1*t*t/(2*mass)
        assert actual < bound and norm(output_gram-np.eye(2))<2e-12
        rows.append(dict(diagnostic_mass=mass,whole_unknown_input_error=actual,
            analytic_moment_bound=float(bound),output_isometry_residual=norm(output_gram-np.eye(2))))
    assert rows[-1]["whole_unknown_input_error"]<rows[0]["whole_unknown_input_error"]
    return dict(scope="Finite internal, continuous one-momentum fiber diagnostic. Every field/source/internal output is kept and Gram norms cover unknown Q with arbitrary reference; not the actual QFT.",
        internal_dimension=24,unknown_dimension=2,quadrature_points=180,
        original_product_after_LLP_internal_purity=purity,
        erroneous_independent_packet_state_error=wrong,
        LLP_input_isometry_residual=norm(gram-np.eye(2)),
        exact_input_moments=dict(A1=float(A1),A2=float(A2),B2=float(B2)),rows=rows)


def boost_matrix(p,m):
    p = np.asarray(p,float)
    energy = np.sqrt(m*m+np.dot(p,p))
    b = np.eye(4)
    b[0,0] = energy/m
    b[0,1:] = b[1:,0] = p/m
    b[1:,1:] += np.outer(p,p)/(m*(energy+m))
    return b


def induced_representation():
    metric = np.diag([1.,-1.,-1.,-1.])
    l1 = boost_matrix([np.sinh(.4),0,0],1)
    l2 = boost_matrix([0,np.sinh(-.3),0],1)
    lorentz = l2@l1
    assert norm(lorentz.T@metric@lorentz-metric)<1e-13
    rows = []
    for m in [1.,3.,7.]:
        for p in [np.array([.2,-.7,.4]),np.array([-1.1,.3,.8])]:
            four = np.r_[np.sqrt(m*m+p@p),p]
            p1 = l1@four
            p2 = l2@p1
            w1 = np.linalg.inv(boost_matrix(p1[1:],m))@l1@boost_matrix(p,m)
            w2 = np.linalg.inv(boost_matrix(p2[1:],m))@l2@boost_matrix(p1[1:],m)
            w12 = np.linalg.inv(boost_matrix(p2[1:],m))@lorentz@boost_matrix(p,m)
            cocycle = norm(w12-w2@w1)
            rot = w12[1:,1:]
            assert cocycle<2e-13 and norm(rot.T@rot-np.eye(3))<2e-13
            assert abs(np.linalg.det(rot)-1)<1e-13
            h = 1e-5
            jac = np.zeros((3,3))
            for j in range(3):
                step = np.eye(3)[j]*h
                def shifted(z):
                    return (lorentz@np.r_[np.sqrt(m*m+z@z),z])[1:]
                jac[:,j] = (shifted(p+step)-shifted(p-step))/(2*h)
            jac_expected = p2[0]/four[0]
            jac_error = abs(np.linalg.det(jac)-jac_expected)
            assert jac_error<3e-9
            rows.append(dict(mass=m,mass_shell_residual=abs(p2@metric@p2-m*m),
                Wigner_cocycle_residual=cocycle,finite_difference_Jacobian_error=jac_error))
    # Independent target-momentum integral of a boosted normalized Gaussian,
    # using the inverse-boost wavefunction and its actual Jacobian amplitude.
    inv = np.linalg.inv(lorentz)
    norms = []
    for n in [64,88]:
        x,w = np.polynomial.legendre.leggauss(n)
        x,w = 10*x,10*w
        p = np.stack(np.meshgrid(x,x,x,indexing="ij"),axis=-1).reshape(-1,3)
        ww = (w[:,None,None]*w[None,:,None]*w[None,None,:]).reshape(-1)
        energy = np.sqrt(3**2+np.sum(p*p,axis=1))
        previous = np.c_[energy,p]@inv.T
        density = (previous[:,0]/energy)*np.exp(-np.sum(previous[:,1:]**2,axis=1))/np.pi**1.5
        val = float(ww@density)
        assert abs(val-1)<2e-8
        norms.append(dict(order=n,target_momentum_norm_squared=val))
    return dict(scope="Three-dimensional positive-mass canonical boosts: independent finite-difference Jacobian, Wigner rotation cocycle and transformed Gaussian norm. General spin/direct-integral representation is proved analytically.",
        rows=rows,boosted_wavefunction_norms=norms,
        common_old_spatial_rotations_and_translations=True,
        old_actual_boost_endpoint_identity_inferred=False)


def compare(actual,expected,path="results"):
    if isinstance(actual,dict):
        assert actual.keys()==expected.keys(),path
        for k in actual: compare(actual[k],expected[k],path+"."+k)
    elif isinstance(actual,list):
        assert len(actual)==len(expected),path
        for i,(a,b) in enumerate(zip(actual,expected)): compare(a,b,path+f"[{i}]")
    elif isinstance(actual,float):
        assert np.isfinite(actual) and np.isclose(actual,expected,atol=3e-10,rtol=3e-10),path
    else:
        assert actual==expected,path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write-results",action="store_true")
    args=p.parse_args()
    result=dict(round=1134,schema=1,status="PASS_INDUCED_POINCARE_WITH_SAME_NORMAL_REST_TASK",
        new_adopted_axioms=0,scientific_calibration=0,stage="active",numerical_check_groups=5,
        source_identity=sources(),exact_budgets=budgets(),normal_coherent_moments=coherent_moments(),
        LLP_entangled_fibers=llp_fibers(),induced_Poincare_diagnostic=induced_representation())
    target=HERE/"results.json"
    if args.write_results:
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("PASS: five groups; only 1134/results.json registered explicitly.")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
        print("PASS: five groups; six frozen source hashes and results match; no files written.")
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=="__main__":
    main()
