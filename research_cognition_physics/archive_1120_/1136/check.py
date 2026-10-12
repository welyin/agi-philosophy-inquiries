"""1136 finite source-deformation certificate. Default is read-only.

Only --write-results writes 1136/results.json. NumPy diagnostics do not replace
the continuum normal-Fock, unbounded-domain or BT-force proof.
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
    "archive_1120_/1131/proof.md": "0801684a9d4eb7a654a960cfd2a8fa916c85a9bb90dacc8638c2a71238e88023",
    "archive_1120_/1132/proof.md": "f05daf49972a0b7fe13f5d593a032998bf1c9008182c051f2e12bb6ea958df63",
    "archive_1120_/1132/results.json": "95dcbd80a2a9e5d3d7ba03e52f1d1aa7e931606c9735f0359623f7b4f969430c",
    "archive_1120_/1134/proof.md": "4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666",
    "archive_1120_/1134/results.json": "e5fb6abfac5ade09cbdb02f3a3719b2114464ebb6f97c4289e537e9175403894",
    "archive_1120_/1135/proof.md": "f3577777e8983a02d02bf66ef5f869b207ac4f94df1610202c1ce6e995e17068",
}


def sources():
    actual = {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FROZEN}
    assert actual == FROZEN
    return dict(sha256=actual, old_clock_dimension=288,
        source_and_Y_preparations_unchanged=True, LLP_entanglement_retained=True,
        Z_and_selector_are_new_normal_material=True,
        internal_kinetic_added_outside_old_support_square_root=True,
        whole_H_is_new_BT_dynamics=True, active_branch_claims_old_handoff=False,
        terminal_read_permission="New bounded cos(p/10) or (I+cos(p/10))/2 effective read; detector work is not included in propagation H.",
        no_global_charge_coupling_assumed=True, local_covariant_stress_claimed=False)


def pack(vals):
    return dict(fractions={k: str(v) for k,v in vals.items()},
                decimal_views={k: float(v) for k,v in vals.items()})


def budgets():
    M=I=mc=F(10**80); sig=F(1,10**30); t=F(5,4); A=F(10**34)
    eps=F(1,10**5); lam=F(11,7); C=F(33536000000000000000)
    v=16*lam*C; mmin=M+mc-v
    cx=F(10,9)*F(24780800000000000000000,49)
    cz=6*v; czz=F(2*10**25)
    dz=(175*t+10*cz*t*t+czz*t*t/2+cz*cz*t**3/3)/(2*I)
    dy=cx*(2*sig*t+t*t/(2*M*sig))
    dd=27*5*t/(8*M*M*sig*sig)
    ds=lam*C*512*3*t/(eps*A)
    assert dz<F(3,10**37)
    assert 5*F(10,9)+t+eps<7
    assert 2+t+eps<4*F(9,10)
    rows={}
    for label,b1,b2 in [("code",268,28000),("finite_menu",600,400000)]:
        p1=F(b1)+94*A; p2=F(b2)+188*b1*A+12000*A*A
        a1=1/sig+p1; a2=1/sig**2+2*p1/sig+p2; bmix=p1/sig+p2
        bt1=t*a2/(2*mmin); bt2=t*(a2+2*bmix)/(2*M)
        bt3=cx*a1*t*t/(2*M); bt=bt1+bt2+bt3
        proc=bt+dy+dd+dz+ds
        assert bt<F(4,10**8) and proc<F(106,10**8)
        rows[label]=pack(dict(p1=p1,p2=p2,A1=a1,A2=a2,B2=bmix,
            BT_mass_shell=bt1,BT_original_fiber=bt2,BT_force=bt3,
            delta_BT=bt,delta_process=proc))
    proc=F(rows["code"]["fractions"]["delta_process"])
    pmenu=F(rows["finite_menu"]["fractions"]["delta_process"])
    angle=F(43,100); cs=1-angle**2/2+angle**4/24
    loss=1-((1+cs)/2)**2; weight=F(21293,100000)
    raw=weight*loss; dm=40*lam*F(27,20)*eps
    gap=raw-8*dm-2*proc
    causal_area=lam*F(27,20)
    tail=causal_area**10/F(math.factorial(10))/(1-causal_area/11)
    causal=4*tail+2*pmenu
    handoff=80*lam*eps+proc
    clock=dy+dd+F(rows["code"]["fractions"]["delta_BT"])
    assert dm==F(297,350000)
    assert gap>F(1215,100000)>F(1,100)
    assert gap/2>F(1,200)
    assert handoff<F(1259,10**6) and causal<F(3,1000)
    bounds=dict(M=M,I=I,new_control_rest_mass=mc,m_min=mmin,v_star=v,
        C_x=cx,C_z=cz,C_zz=czz,delta_Z=dz,delta_Y=dy,delta_D=dd,delta_source=ds,
        mode_d_max=dm,central_weight_lower=weight,cos_upper=cs,
        fidelity_loss_lower=loss,raw_cos_gap=raw,mode_corrected_gap=raw-8*dm,
        actual_cos_gap=gap,actual_probability_gap=gap/2,
        baseline_handoff_error=handoff,baseline_bare_clock_error=clock,
        finite_causal_tail=tail,finite_causal_error=causal,
        control_mean_kinetic=F(1)/(8*I*F(1,20)**2),
        control_kinetic_norm_upper=F(175)/(2*I),
        added_control_and_support_rest_mass=M+mc)
    return dict(common=pack(bounds),branches=rows,
        cosine_baseline=math.exp(-.5),science_experiment_count=0)


def annihilators(n):
    out=[]
    for j in range(n):
        a=np.zeros((1<<n,1<<n),complex)
        for state in range(1<<n):
            if state>>j&1:
                a[state^(1<<j),state]=(-1)**((state&((1<<j)-1)).bit_count())
        out.append(a)
    return out


def norm(x):
    return float(np.linalg.norm(x,2))


def exp_h(h,t):
    e,v=np.linalg.eigh(h)
    return (v*np.exp(-1j*t*e))@v.conj().T


def car_checks():
    a=annihilators(4); eye=np.eye(16)
    car=max(norm(a[i]@a[j].conj().T+a[j].conj().T@a[i]-(eye if i==j else 0))
            for i in range(4) for j in range(4))
    # The four local sea modes have covariance I/2. Conditioning only A
    # leaves B maximally mixed; columns retain its four-dimensional purifier.
    bases=[]
    for rail in [0,1]:
        bases.append(np.eye(16)[:,[(1<<rail)|(b<<2) for b in range(4)]])
    gen=sum(a[i].conj().T@a[i+2]+a[i+2].conj().T@a[i] for i in range(2))
    max_fid=max_effect=max_gram=0.
    nb=a[2].conj().T@a[2]
    for alpha,beta in [(0,.2),(.6,1.1),(1.4,1.9),(-.2,2.2),(0,math.pi/2)]:
        ua,ub=exp_h(gen,alpha),exp_h(gen,beta)
        va,vb=ua@bases[0],ub@bases[0]
        fid=float(np.linalg.svd(va.conj().T@vb,compute_uv=False).sum()/4)
        target=((1+abs(math.cos(alpha-beta)))/2)**2
        max_fid=max(max_fid,abs(fid-target))
        # Complete unknown-Q isometry with retained sea purifier.
        iso=np.stack([(ua@b).reshape(-1)/2 for b in bases],axis=1)
        max_gram=max(max_gram,norm(iso.conj().T@iso-np.eye(2)))
        effect=iso.conj().T@np.kron(nb,np.eye(4))@iso
        wanted=.5*math.cos(alpha)**2*np.eye(2)+math.sin(alpha)**2*np.diag([1,0])
        max_effect=max(max_effect,norm(effect-wanted))
        bell=iso.reshape(64,2).reshape(-1)/math.sqrt(2)
        assert abs(np.vdot(bell,bell)-1)<1e-12
    assert car<1e-12 and max_fid<1e-12 and max_effect<1e-12 and max_gram<1e-12
    # Independent full positive-particle / hole realization of global charge.
    c=annihilators(8); vac=np.eye(256)[:,0]
    g=sum(c[i].conj().T@c[i+2]+c[i+2].conj().T@c[i]
          -c[i+4].conj().T@c[i+6]-c[i+6].conj().T@c[i+4] for i in range(2))
    jj=[]
    for i,j in [(0,1),(1,0)]:
        jj.append((c[i].conj().T@vac-c[i].conj().T@c[j].conj().T@c[j+4].conj().T@vac)/math.sqrt(2))
    J=np.stack(jj,axis=1); ge,gv=np.linalg.eigh(g); coords=gv.conj().T@J
    global_error=0.
    for th in [0,.13,.7,math.pi/2,2.3]:
        char=coords.conj().T@(np.exp(-1j*th*ge)[:,None]*coords)
        global_error=max(global_error,norm(char-(7*math.cos(th)+math.cos(3*th))/8*np.eye(2)))
    moments={str(n):norm(J.conj().T@np.linalg.matrix_power(g,n)@J-target*np.eye(2))
             for n,target in [(1,0),(2,2),(4,11)]}
    assert global_error<1e-12 and max(moments.values())<1e-11
    return dict(scope="Full 16-dimensional local CAR, retained four-dimensional sea purification and arbitrary-Q matrix units; separate 256-dimensional particle-hole global diagnostic.",
        CAR_residual=car,root_fidelity_residual=max_fid,complete_Q_Gram_residual=max_gram,
        B_effect_residual=max_effect,global_character_residual=global_error,
        global_moment_residuals=moments,global_charge_is_only_diagnostic=True,
        B_initial_occupation=.5,global_charge_variance=2,global_charge_fourth_moment=11)


def strain_grid(n):
    x=np.linspace(0,1,n)
    bump=np.zeros(n); interior=(x>0)&(x<1)
    bump[interior]=np.exp(-1/(x[interior]*(1-x[interior])))
    def cum(y):
        return np.r_[0,np.cumsum((y[:-1]+y[1:])*.5/(n-1))]
    R=cum(bump); R/=R[-1]; integral_R=cum(R)
    def strain(z):
        v=np.abs(z); u=np.clip((v-.075)/.025,0,1)
        middle=.075+.025*(u-np.interp(u,x,integral_R))
        return np.sign(z)*np.where(v<=.075,v,middle)
    return strain,float(.075+.025*(1-integral_R[-1]))


def gaussian_diagnostic():
    sig=.05; shift=.1; baseline=math.exp(-shift**2/(8*sig**2))
    rows=[]
    # The exact continuous Gaussian measure is evaluated by Hermite quadrature;
    # no finite-dimensional CCR or truncated Z preparation is substituted.
    for n,grid in [(192,16385),(256,32769)]:
        s,cap=strain_grid(grid); x,w=np.polynomial.hermite.hermgauss(n)
        z=-shift/2+math.sqrt(2)*sig*x
        delta=math.pi/2*(np.exp(-3*s(z))-np.exp(-3*s(z+shift)))
        fidelity=((1+np.abs(np.cos(delta)))/2)**2
        active=float(baseline*np.dot(w,fidelity)/math.sqrt(math.pi))
        rows.append(dict(nodes=n,strain_grid=grid,actual_saturation=cap,
            comparison_cosine=active,comparison_gap=baseline-active,
            normal_Z_variance=sig**2,p_variance=1/(4*sig**2)))
    convergence=abs(rows[0]["comparison_cosine"]-rows[1]["comparison_cosine"])
    assert convergence<2e-5 and rows[-1]["comparison_gap"]>.03
    min_angle=math.pi/2*(math.exp(.075)-math.exp(-.225))
    weight=baseline*math.erf(.5/math.sqrt(2))
    central=weight*(1-((1+math.cos(min_angle))/2)**2)
    assert central>.02148
    return dict(scope="Continuous normal Gaussian and finite-CAR comparison diagnostic. It does not assume the actual infinite sea attains maximal fidelity; the proof uses only an upper overlap bound.",
        baseline=baseline,rows=rows,quadrature_convergence=convergence,
        central_exact_numeric_weight=weight,central_numeric_min_angle=min_angle,
        central_numeric_gap=central,continuous_CCR_truncated=False)


def bt_force_diagnostic():
    rng=np.random.default_rng(1136); n=11
    u,_=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    base=(u*np.linspace(8,23,n))@u.conj().T
    x=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); v=(x+x.conj().T)/2
    v/=norm(v); z=.37; mass=base+math.sin(z)*v; derivative=math.cos(z)*v
    me,mv=np.linalg.eigh(mass); m=float(me.min()); C=norm(derivative)
    rows=[]
    def sqrt_mass(ma,k):
        e,w=np.linalg.eigh(ma)
        return (w*np.sqrt(e*e+k*k))@w.conj().T
    for k in [0.,1.,7.,40.]:
        he=np.sqrt(me*me+k*k); H=(mv*he)@mv.conj().T; gg=H-mass
        rhs=mass@derivative+derivative@mass
        D=mv@((mv.conj().T@rhs@mv)/(he[:,None]+he[None,:]))@mv.conj().T
        bounded_rhs=gg@derivative+derivative@gg
        E=-mv@((mv.conj().T@bounded_rhs@mv)/(he[:,None]+he[None,:]))@mv.conj().T
        eps=1e-5
        fd=(sqrt_mass(base+math.sin(z+eps)*v,k)-sqrt_mass(base+math.sin(z-eps)*v,k))/(2*eps)
        residual=norm(D-derivative-E); difference=norm(D-derivative)
        bound=(math.sqrt(m*m+k*k)-m)/math.sqrt(m*m+k*k)*C
        assert residual<2e-12 and norm(D)<=2*C+1e-12 and difference<=bound+1e-12
        assert norm(D-fd)<2e-8
        rows.append(dict(k=k,negative_Sylvester_residual=residual,
            finite_difference_residual=norm(D-fd),actual_BT_derivative_norm=norm(D),
            old_internal_derivative_norm=C,derivative_difference=difference,
            derivative_difference_bound=bound))
    return dict(scope="Noncommuting finite positive masses independently solve the Sylvester derivative; continuum force and domain are analytic in proof.",
        dimension=n,noncommuting_mass_derivative=norm(mass@derivative-derivative@mass),rows=rows,
        p_is_spatial_momentum=False,free_scalar_force_formula_inherited=False)


def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f"/{i}")
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-7,abs_tol=2e-10),(path,a,b)
    else: assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args()
    result=dict(round=1136,status="PASS",groups=5,new_axioms=0,science_experiments=0,stage="active",
        source_identity=sources(),exact_budgets=budgets(),full_CAR=car_checks(),
        normal_Gaussian=gaussian_diagnostic(),BT_force=bt_force_diagnostic())
    out=HERE/"results.json"
    if args.write_results: out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else: compare(result,json.loads(out.read_text(encoding="utf-8")))
    print("1136 PASS: 5 mathematical groups; 8 frozen sources; exact budgets; full CAR; normal Gaussian; noncommuting BT force.")


if __name__=="__main__": main()
