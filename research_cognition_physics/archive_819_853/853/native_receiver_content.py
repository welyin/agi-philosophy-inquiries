"""853: native mediated receiver calibration with original CAR content.

Three Gaussian oscillator mediators and three two-level receivers are an
explicit finite controlled diagnostic, NOT a finite-coupling solution of
the continuum model. Exact original 829 source probabilities are used.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'native_receiver_content_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1.,-1.]).astype(complex)

def run():
    gamma,_,_,_,comp,sx,sz=code.code_data()
    sy=[code.mul((0,0,3),code.mul(sz[b],sx[b])) for b in range(5)]
    sources=[code.mul((0,0,2),v) for v in (sz[0],sz[1],sy[3])]
    logical=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    eps=F(1,10);r0=F(0);r1=F(3,5);d=1024
    def phase(p):
        assert p in (0,2)
        return F(1 if p==0 else -1)
    def core(p):
        typ,ph=comp(p);return phase(ph) if typ=='scalar' else F(0)
    def ex(p,r):
        full=phase(p[2]) if p[:2]==(0,0) else F(0)
        return ((1-eps)-2*eps/(d-2))*core(p)+(1-eps)*r*core(code.mul(logical,p))+eps*d/(d-2)*full
    signs=list(itertools.product((-1,1),repeat=3))
    def joint(a,r):
        ans=F(0)
        for subset in itertools.product((0,1),repeat=3):
            word=code.product(sources[i] for i in range(3) if subset[i])
            ans+=F(np.prod([a[i] for i in range(3) if subset[i]],dtype=int))*ex(word,r)/8
        return ans
    inputs={str(r):[joint(a,r) for a in signs] for r in (r0,r1)}
    for r in (r0,r1):
        p=inputs[str(r)];assert sum(p)==1 and min(p)>0
        assert all(p[k]==(1+(1-eps)*r*np.prod(a))/8 for k,a in enumerate(signs))
    # The physical source is hbar*a_i*A_i. The two native diagnostic pulses:
    # U_C=exp(-i lambda_i sqrt(hbar) a_i A_i Q_i),
    # U_R=exp(+i eta_i sqrt(hbar) P_i Y_i/2).
    # U_C^* P_i U_C=P_i-lambda_i sqrt(hbar) a_i A_i.
    lambdas=(F(2,7),F(3,8),F(4,9));etas=(F(3,5),F(4,7),F(5,8))
    amplitudes=(F(1),F(1),F(5,6));variances=(F(7,10),F(5,8),F(9,11))
    hbar=F(1,3);hf=float(hbar)
    states={'plus_x':(I+X)/2,'plus_z':(I+Z)/2,'unpolarized':I/2}
    susceptibility={name:float(np.trace(rho@(1j*((-Y/2)@Z-Z@(-Y/2)))).real) for name,rho in states.items()}
    assert susceptibility=={'plus_x':1.,'plus_z':0.,'unpolarized':0.}
    def closed_means(name,h):
        bx=float(np.trace(states[name]@X).real);bz=float(np.trace(states[name]@Z).real)
        out={}
        for i in range(3):
            eta=float(etas[i]);lam=float(lambdas[i]);a=float(amplitudes[i]);nu=float(variances[i])
            damp=np.exp(-eta*eta*h*nu/2)
            for sign in (-1,1):out[i,sign]=damp*(bz*np.cos(eta*lam*h*a)-bx*sign*np.sin(eta*lam*h*a))
        return out
    max_quad=0.;max_two_orders=0.;unitarity=0.
    for name,rho in states.items():
        old=None
        for n in (32,64):
            nodes,weights=np.polynomial.hermite.hermgauss(n);weights=weights/np.sqrt(np.pi)
            means={}
            for i in range(3):
                eta=float(etas[i]);lam=float(lambdas[i]);amp=float(amplitudes[i]);nu=float(variances[i])
                for sign in (-1,1):
                    val=0.
                    for t,w in zip(nodes,weights):
                        p=np.sqrt(2*nu)*t-lam*np.sqrt(hf)*amp*sign
                        theta=eta*np.sqrt(hf)*p
                        u=np.cos(theta/2)*I+1j*np.sin(theta/2)*Y
                        unitarity=max(unitarity,float(np.linalg.norm(u.conj().T@u-I)))
                        val+=w*float(np.trace(rho@u.conj().T@Z@u).real)
                    means[i,sign]=val
                    max_quad=max(max_quad,abs(val-closed_means(name,hf)[i,sign]))
            if old is not None:max_two_orders=max(max_two_orders,max(abs(means[k]-old[k]) for k in means))
            old=means
    assert max_quad<2e-14 and max_two_orders<2e-14 and unitarity<2e-14
    def output(name,r,h):
        means=closed_means(name,h);out=[]
        for z in signs:
            out.append(sum(float(joint(a,r))*np.prod([(1+z[i]*means[i,a[i]])/2 for i in range(3)]) for a in signs))
        return out
    rows=[]
    prefactor=float((1-eps)*(r1-r0))
    for name in states:
        p0=output(name,r0,hf);p1=output(name,r1,hf);delta=np.array(p1)-p0
        assert min(p0+p1)>=0 and abs(sum(p0)-1)<1e-14 and abs(sum(p1)-1)<1e-14
        parity=float(sum(np.prod(z)*v for z,v in zip(signs,delta)))
        if name=='plus_x':
            analytic=-prefactor*np.prod([np.exp(-float(etas[i])**2*hf*float(variances[i])/2)*np.sin(float(etas[i]*lambdas[i]*amplitudes[i])*hf) for i in range(3)])
            assert abs(parity-analytic)<2e-15 and abs(parity)>1e-6
            assert max(abs(delta[k]-np.prod(z)*analytic/8) for k,z in enumerate(signs))<2e-15
        else:assert max(abs(delta))<2e-15
        rows.append(dict(preparation=name,susceptibility=susceptibility[name],min_probability=min(p0+p1),parity_content_difference=parity,max_bitstring_probability_difference=float(max(abs(delta)))))
    exact_leading=-(1-eps)*(r1-r0)*np.prod(lambdas)*np.prod(etas)*np.prod(amplitudes)
    assert isinstance(exact_leading,F)
    asymptotic=[]
    for h in (1.,.5,.25,.125):
        parity=-prefactor*np.prod([np.exp(-float(etas[i])**2*h*float(variances[i])/2)*np.sin(float(etas[i]*lambdas[i]*amplitudes[i])*h) for i in range(3)])
        asymptotic.append(dict(hbar=h,parity=parity,coefficient_after_dividing_hbar_cubed=parity/h**3))
    errs=[abs(row['coefficient_after_dividing_hbar_cubed']-float(exact_leading)) for row in asymptotic]
    assert all(a>b for a,b in zip(errs,errs[1:]))
    # Two pure one-particle receiver preparations can have equal free energy
    # and entropy but different read susceptibility. Flavor-blind H=E*N.
    H=np.diag([0.,2.,2.,4.]);rho_x=np.zeros((4,4),complex);rho_z=rho_x.copy()
    rho_x[1:3,1:3]=(I+X)/2;rho_z[1:3,1:3]=(I+Z)/2
    energy_x=float(np.trace(H@rho_x).real);energy_z=float(np.trace(H@rho_z).real)
    assert energy_x==energy_z==2. and np.allclose(rho_x@rho_x,rho_x) and np.allclose(rho_z@rho_z,rho_z)
    return dict(round=853,all_checks_passed=True,fresh_test_groups=1,
        diagnostic='original CAR source distribution plus exact controlled Gaussian mediator/receiver model; not the original continuum finite-coupling dynamics',
        original_code_joint_probabilities={k:[str(x) for x in v] for k,v in inputs.items()},
        source_amplitudes=[str(x) for x in amplitudes],lambdas=[str(x) for x in lambdas],etas=[str(x) for x in etas],probe_variances=[str(x) for x in variances],
        preparation_rows=rows,exact_hbar_cubed_leading_coefficient=str(exact_leading),asymptotic_rows=asymptotic,
        max_matrix_gaussian_vs_closed_error=max_quad,max_32_vs_64_quadrature_error=max_two_orders,max_receiver_unitarity_error=unitarity,
        equal_free_receiver_energies=[energy_x,energy_z],both_comparison_preparations_pure=True,
        new_continuum_receiver_species_and_profiles_are_explicit_inputs=True,
        full_curved_PDE_or_loop_coefficients_computed=False,original_branch_A_disproved=False,
        finite_coupling_continuum_probabilities_or_autonomous_apparatus_claimed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();result=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        old=json.loads(TARGET.read_text('utf-8'))
        assert result==old
    print(json.dumps(result,ensure_ascii=False,indent=2))
