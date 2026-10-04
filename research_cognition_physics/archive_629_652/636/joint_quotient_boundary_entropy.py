"""636: original quotient boundary records, entropy and electric sources.

Exact finite-support states of the ORIGINAL full Gauss Hilbert space.
They are not full-Hamiltonian equilibrium states or invariant sectors.
Hbar=1; all electric energies acquire hbar**2 in other units.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_region_energy_gluing as gluing

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_quotient_boundary_entropy_results.json'
MODULES = dict(Q=(1,0,1,1), u=(1,0,0,4), d=(1,0,0,-2),
               L=(0,0,1,-3), e=(0,0,0,-6), nu=(0,0,0,0))


def label_data(label):
    a,b,ell,q = label
    assert min(a,b,ell) >= 0
    return dict(label=list(label),
        residue=(2*(a+2*b)+3*ell+q) % 6,
        dimension=(a+1)*(b+1)*(a+b+2)*(ell+1)//2,
        casimirs=[float(Fraction(a*a+b*b+a*b+3*a+3*b,3)),
                  float(Fraction(ell*(ell+2),4)), float(q*q)])


def entropy(p):
    p = np.asarray(p, float)
    return float(-np.sum(p[p>0]*np.log(p[p>0])))


def entropy_parts(p, dimensions):
    center = entropy(p)
    edge = float(2*np.dot(p,np.log(dimensions)))
    return dict(center=center, representation_indices=edge, extended=center+edge)


def loop_holonomy(edges):
    bottom,left,right,top = edges
    return gluing.product(gluing.product(gluing.product(bottom,right),
                          gluing.inverse(top)),gluing.inverse(left))


def quotient_check():
    rng = np.random.default_rng(636)
    zeta = (np.exp(2j*np.pi/3)*np.eye(3),-np.eye(2),np.exp(1j*np.pi/3))
    data = [label_data(x) for x in MODULES.values()]
    errors = []
    for name,row in zip(MODULES,data):
        assert row['residue'] == 0
        R = gluing.rep(zeta,name)
        assert R.shape == (row['dimension'],)*2
        errors.append(float(np.max(abs(R-np.eye(len(R))))))
    # Actual original module characters, four links and four local gauge changes.
    endpoints = [(0,1),(0,2),(1,3),(2,3)]
    for _ in range(8):
        edges = [gluing.sample(rng) for _ in endpoints]
        gs = [gluing.sample(rng) for _ in range(4)]
        transformed = [gluing.product(gluing.product(gs[s],u),gluing.inverse(gs[t]))
                       for u,(s,t) in zip(edges,endpoints)]
        U,V = loop_holonomy(edges),loop_holonomy(transformed)
        for name in MODULES:
            errors.append(float(abs(np.trace(gluing.rep(U,name))-
                                    np.trace(gluing.rep(V,name)))))
    # Exact product of the three ORIGINAL label marginals, not a changed model.
    marginals = []
    labels = list(MODULES.values())
    for select in (lambda x:x[:2],lambda x:x[2],lambda x:x[3]):
        marginal = {}
        for label in labels:
            key = select(label)
            marginal[key] = marginal.get(key,Fraction(0))+Fraction(1,6)
        marginals.append(marginal)
    allowed = Fraction(0)
    for color,weak,q in itertools.product(*[list(x) for x in marginals]):
        w = marginals[0][color]*marginals[1][weak]*marginals[2][q]
        if label_data((*color,weak,q))['residue'] == 0:
            allowed += w
    assert allowed == Fraction(5,18)
    p = np.full(6,1/6)
    C = np.array([x['casimirs'] for x in data])
    mean = p@C
    covariance = (C.T*p)@C-np.outer(mean,mean)
    marginal_H = sum(entropy(list(m.values())) for m in marginals)
    assert np.allclose(mean,[2/3,1/4,11],rtol=0,atol=2e-14)
    assert abs(covariance[0,1])<1e-14
    assert abs(covariance[0,2]+8/3)<1e-14
    assert abs(covariance[1,2]+1.5)<1e-14
    assert max(errors)<4e-14 and marginal_H>entropy(p)+1
    return dict(modules={k:v for k,v in zip(MODULES,data)},
        original_center_and_loop_gauge_error=max(errors),
        actual_joint_entropy=entropy(p),
        sum_of_three_marginal_entropies=marginal_H,
        spurious_entropy=marginal_H-entropy(p),
        product_marginal_allowed_weight=str(allowed),
        product_marginal_forbidden_weight=str(1-allowed),
        casimir_means=mean.tolist(),casimir_covariance=covariance.tolist(),
        entropy=entropy_parts(p,[x['dimension'] for x in data]))


def schmidt_check():
    rows=[]
    for d in (1,2,3,6,8,10,27):
        # Coefficient matrix of chi_R(U_A U_B), in normalized PW bases.
        coefficient=np.zeros((d*d,d*d))
        for i,j in itertools.product(range(d),repeat=2):
            coefficient[i*d+j,j*d+i]=1/d
        reduced=coefficient@coefficient.T
        error=float(np.max(abs(reduced-np.eye(d*d)/(d*d))))
        matrix_entropy=entropy(np.diag(reduced))
        assert error<2e-15
        assert abs(np.trace(reduced)-1)<2e-15
        assert abs(matrix_entropy-2*np.log(d))<4e-14
        rows.append(dict(dimension=d,partial_trace_error=error,
                         entropy=matrix_entropy))
    # Mixture: trace the original six representation blocks separately.
    p=np.full(6,1/6);dims=[label_data(x)['dimension'] for x in MODULES.values()]
    spectrum=np.concatenate([np.full(d*d,w/(d*d)) for w,d in zip(p,dims)])
    expected=entropy_parts(p,dims)
    assert abs(entropy(spectrum)-expected['extended'])<2e-14
    return dict(rows=rows,original_module_mixture_spectrum=spectrum.tolist(),
                mixture_entropy=entropy(spectrum),decomposed_entropy=expected,
                two_cut_edges_same_representation_not_two_independent_labels=True)


def electric_check():
    # Local Gauss cancellation of the TWO outgoing active edges at vertex 0.
    T=gluing.gauge.generators(3);I=np.eye(3)
    singlet=I.reshape(-1)/np.sqrt(3)
    P=[np.kron(t,I) for t in T]
    Q=[-np.kron(I,t.T) for t in T]
    gauss=max(float(np.linalg.norm((p+q)@singlet)) for p,q in zip(P,Q))
    cross=sum(p@q for p,q in zip(P,Q))
    casimir=sum(p@p for p in P)
    cross_error=float(np.linalg.norm(cross@singlet+(4/3)*singlet))
    casimir_error=float(np.linalg.norm(casimir@singlet-(4/3)*singlet))
    assert max(gauss,cross_error,casimir_error)<2e-14
    # Original 617 non-diagonal geometry, with original COLOR coefficient.
    S=np.array([[.24,.31,-.12],[.31,-.07,.17],[-.12,.17,-.17]])
    eps=.73; b=float(gluing.geometry.PAR['b'][0])
    def coefficients(t,sigma):
        K=b/eps*np.exp(-2*sigma)*gluing.geometry.shape_exp(S,t)
        return K
    def contract(K):
        # Vertex 0: xx+yy-2xy. Vertex 1: yy. Vertex 2: xx.
        return float(2*K[0,0]+2*K[1,1]-2*K[0,1])
    t=.8;sigma=.12
    K=coefficients(t,sigma)
    kappa=contract(K)
    source=np.array([kappa,contract(K@S),-2*kappa])
    step=2e-6
    fd=np.array([(contract(coefficients(t+step,sigma))-
                  contract(coefficients(t-step,sigma)))/(2*step),
                 (contract(coefficients(t,sigma+step))-
                  contract(coefficients(t,sigma-step)))/(2*step)])
    fd_error=float(np.max(abs(fd-source[1:])))
    assert fd_error<3e-9
    irreps=[label_data(x) for x in ((0,0,0,0),(1,1,0,0),(3,0,0,0),(2,2,0,0))]
    assert all(x['residue']==0 for x in irreps)
    C=np.array([x['casimirs'][0] for x in irreps])
    dims=np.array([x['dimension'] for x in irreps])
    # Independent inverse-Cartan contraction checks every SU3 Casimir.
    inv_cartan=np.array([[2.,1.],[1.,2.]])/3
    for row in irreps:
        weight=np.array(row['label'][:2])
        assert abs(.5*weight@inv_cartan@(weight+2)-row['casimirs'][0])<1e-14
    delta=np.array([-5,16,-20,9])/160
    assert np.max(abs(np.stack([C**n for n in range(3)])@delta))<2e-14
    distributions=[np.full(4,.25)+delta,np.full(4,.25)-delta]
    outputs=[]
    for p in distributions:
        assert np.min(p)>0 and abs(sum(p)-1)<1e-14
        moments=np.array([p@C**n for n in range(1,4)])
        var=moments[1]-moments[0]**2
        outputs.append(dict(probabilities=p.tolist(),moments=moments.tolist(),
            electric_source_means=(moments[0]*source).tolist(),
            electric_source_covariance=(var*np.outer(source,source)).tolist(),
            entropy=entropy_parts(p,dims)))
    means_error=float(np.max(abs(np.array(outputs[0]['electric_source_means'])-
                                np.array(outputs[1]['electric_source_means']))))
    cov_error=float(np.max(abs(np.array(outputs[0]['electric_source_covariance'])-
                              np.array(outputs[1]['electric_source_covariance']))))
    difference={k:outputs[0]['entropy'][k]-outputs[1]['entropy'][k]
                for k in outputs[0]['entropy']}
    assert means_error<1e-13 and cov_error<1e-12
    assert abs(difference['extended'])>.1 and abs(difference['center'])>.001
    return dict(original_color_b=b,original_K=K.tolist(),
        kappa_and_geometry_derivatives=source.tolist(),
        omitted_off_diagonal_kappa_error=float(2*K[0,1]),
        source_finite_difference_error=fd_error,
        fundamental_gauss_error=gauss,cross_Casimir_error=cross_error,
        fundamental_Casimir_error=casimir_error,
        original_quotient_loop_irreps=irreps,states=outputs,
        source_mean_difference=means_error,source_covariance_difference=cov_error,
        entropy_difference_plus_minus=difference,
        full_Hamiltonian_equilibrium_or_source_equivalence_claimed=False)


def run():
    deps=('research_note_360.md','research_note_362.md','research_note_365.md',
          'research_note_589.md','research_note_603.md','research_note_617.md',
          'research_note_635.md','joint_region_energy_gluing.py',
          'joint_fermion_gauss_completion.py','joint_full_spatial_metric.py',
          'joint_quotient_gauge_completion.py','joint_curved_quantum_source.py',
          'joint_entropy_geometric_matching_results.json')
    return dict(round=636,tests_run=3,failures=0,errors=0,
        quotient_records=quotient_check(),boundary_schmidt=schmidt_check(),
        electric_geometry_and_entropy=electric_check(),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                           for name in deps},
        scope=dict(original_full_Gauss_state_space_used=True,
            original_nondiagonal_electric_geometry_retained=True,
            full_magnetic_Higgs_Yukawa_hopping_not_deleted=True,
            explicitly_chosen_finite_support_states_not_Gibbs=True,
            electric_source_moments_only_not_total_stress=True,
            no_new_particle_spectrum_from_loop_labels=True,
            entropy_algebra_and_trace_declared=True,
            continuous_gauge_entropy_and_Newton_matching_open=True,
            continuum_dimension_and_GR_still_open=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
