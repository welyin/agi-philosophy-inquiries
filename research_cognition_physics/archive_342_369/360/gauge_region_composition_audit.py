"""Round 360: a finite Z2 Gauss constraint separates regions from independent systems."""
import argparse
import itertools
import json
from pathlib import Path
import unittest
import numpy as np

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],dtype=complex)
Y=np.array([[0,-1j],[1j,0]],dtype=complex)
Z=np.diag([1,-1]).astype(complex)
PAULI=(I,X,Y,Z)
PLUS=np.array([1,1],dtype=complex)/np.sqrt(2)
MINUS=np.array([1,-1],dtype=complex)/np.sqrt(2)


def kron_all(parts):
    ans=np.array([1.0],dtype=complex)
    for part in parts:
        ans=np.kron(ans,part)
    return ans


def ket_density(v):
    return np.outer(v,v.conj())


def trace_distance(a,b):
    return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)


def span_rank(ops):
    return int(np.linalg.matrix_rank(np.stack([o.ravel() for o in ops],axis=1),tol=1e-10))


def cycle(n=4):
    eye=np.eye(2**n,dtype=complex)
    electric=[kron_all([X if j==i else I for j in range(n)]) for i in range(n)]
    gauss=[electric[(i-1)%n]@electric[i] for i in range(n)]
    projector=eye.copy()
    for gen in gauss:
        projector=projector@(eye+gen)/2
    embedding=np.column_stack([kron_all([PLUS]*n),kron_all([MINUS]*n)])
    loop=kron_all([Z]*n)
    return {'n':n,'eye':eye,'electric':electric,'gauss':gauss,
            'P':projector,'V':embedding,'W':loop}


def gauge_paulis(model,support):
    n=model['n']; ans=[]
    for labels in itertools.product(range(4),repeat=len(support)):
        parts=[I]*n
        for site,label in zip(support,labels):
            parts[site]=PAULI[label]
        op=kron_all(parts)
        if all(np.linalg.norm(op@g-g@op)<1e-10 for g in model['gauss']):
            ans.append(op)
    return ans


def compress(model,op):
    return model['V'].conj().T@op@model['V']


def partial_bipartite(rho,da,db,keep):
    four=rho.reshape(da,db,da,db)
    return np.trace(four,axis1=1,axis2=3) if keep==0 else np.trace(four,axis1=0,axis2=2)


def split_embedding():
    v=np.zeros((4,2),dtype=complex)
    v[0,0]=v[3,1]=1
    return v


def projected(rho,projector):
    raw=projector@rho@projector
    prob=float(np.trace(raw).real)
    if prob<=0:
        raise ValueError('zero success probability')
    return raw/prob,prob


def entropy(rho):
    vals=np.linalg.eigvalsh(rho)
    vals=vals[vals>1e-12]
    return float(-np.sum(vals*np.log2(vals)))


def report():
    m=cycle()
    a=[compress(m,o) for o in gauge_paulis(m,(0,1))]
    b=[compress(m,o) for o in gauge_paulis(m,(2,3))]
    full=[compress(m,o) for o in gauge_paulis(m,range(4))]
    rp=ket_density(PLUS); rm=ket_density(MINUS)
    v=split_embedding(); p=v@v.conj().T
    zero=np.array([1,0],complex)
    inp0=ket_density(np.kron(zero,zero))
    inp1=ket_density(np.kron(PLUS,PLUS))
    out0,s0=projected(inp0,p); out1,s1=projected(inp1,p)
    outmix,sm=projected((inp0+inp1)/2,p)
    return {
      'scope':'Specified finite Z2 gauge constraint and fixed regional electric algebras; not a derivation of gauge symmetry, gravity, general gravitational factorization, or failure of complex-quantum local tomography for independent systems.',
      'hypothesis_tested':'Do adjacent constrained regions automatically form two independently preparable full quantum types?',
      'cycle':{'links':4,'kinematic_dimension':16,'physical_dimension':2,
               'gauge_invariant_pauli_count':len(full),
               'physical_full_algebra_rank':span_rank(full),
               'region_A_algebra_rank':span_rank(a),'region_B_algebra_rank':span_rank(b),
               'shared_algebra_rank':span_rank(a)+span_rank(b)-span_rank(a+b),
               'product_regional_observable_rank':span_rank([x@y for x in a for y in b])},
      'phase_witness':{'physical_trace_distance':trace_distance(rp,rm),
                       'max_regional_product_expectation_difference':float(max(abs(np.trace((rp-rm)@x@y)) for x in a for y in b)),
                       'wilson_loop_expectations':[float(np.trace(r@compress(m,m['W'])).real) for r in (rp,rm)]},
      'boundary_extension':{'extended_dimension':4,'sewn_dimension':2,
                            'physical_projector':'(I+Z_A Z_B)/2',
                            'independent_maximally_mixed_probability_in_physical_sector':0.5,
                            'postselection_success':[s0,s1,sm],
                            'normalized_branch_nonaffinity_witness':trace_distance(outmix,(out0+out1)/2),
                            'normalized_branch_nonaffinity_exact':1/(6*np.sqrt(2))},
      'energy_account':{'penalty':'Delta/2 sum_v (I-G_v), Delta>0',
                        'smallest_positive_energy_over_Delta':2,
                        'open_string_endpoint_cost_over_Delta':2,
                        'closed_loop_endpoint_cost_over_Delta':0},
      'entropy_bits':[{'p':p0,'regional':entropy(np.diag([p0,1-p0])),
                       'pure_global':0.0,'dephased_global':entropy(np.diag([p0,1-p0]))}
                      for p0 in (0.2,0.5,0.8)],
      'proved':[
        'All gauge-invariant observables strictly supported on a proper part of a connected cycle reduce to the shared span of I and logical Z.',
        'The complete physical algebra is M2(C), including the nonlocal Wilson loop and complex phase.',
        'Independent regional preparation cannot give arbitrary flux pairs under the gluing constraint.',
        'An isometric boundary encoding preserves every unknown logical state with its reference; neither reduced region is a copy.',
        'Projecting independently prepared boundary states has a success cost and its normalized branch is not an affine channel.',
        'A full tensor-product ultraviolet space and a constrained nonfactorizing low-energy sector can coexist.'],
      'not_proved':['Cognitive origin of Gauss constraints','A preferred regional algebra',
                    'Gravitational boundary algebra or Einstein dynamics','Universal absence of all tensor decompositions'],
      'sources':['https://arxiv.org/pdf/1312.1183','https://arxiv.org/html/1601.04744v2']}


class Checks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.m=cycle()

    def test_01_gauss_projector_and_redundancy(self):
        m=self.m
        np.testing.assert_allclose(kron_all([I]*4),np.linalg.multi_dot(m['gauss']),atol=1e-12)
        np.testing.assert_allclose(m['P']@m['P'],m['P'],atol=1e-12)
        self.assertEqual(np.linalg.matrix_rank(m['P'],tol=1e-10),2)
        for g,h in itertools.product(m['gauss'],repeat=2):
            np.testing.assert_allclose(g@h,h@g,atol=1e-12)

    def test_02_encoding_and_complete_logical_algebra(self):
        m=self.m; v=m['V']
        np.testing.assert_allclose(v.conj().T@v,I,atol=1e-12)
        np.testing.assert_allclose(v@v.conj().T,m['P'],atol=1e-12)
        np.testing.assert_allclose(compress(m,m['electric'][0]),Z,atol=1e-12)
        np.testing.assert_allclose(compress(m,m['W']),X,atol=1e-12)
        np.testing.assert_allclose(compress(m,-1j*m['electric'][0]@m['W']),Y,atol=1e-12)
        full=gauge_paulis(m,range(4))
        self.assertEqual(len(full),32)
        self.assertEqual(span_rank([compress(m,o) for o in full]),4)

    def test_03_regional_algebras_and_shared_center(self):
        m=self.m
        a=[compress(m,o) for o in gauge_paulis(m,(0,1))]
        b=[compress(m,o) for o in gauge_paulis(m,(2,3))]
        self.assertEqual(len(a),4)
        self.assertEqual(span_rank(a),2)
        self.assertEqual(span_rank(b),2)
        self.assertEqual(span_rank(a+b),2)
        self.assertEqual(span_rank([x@y for x in a for y in b]),2)

    def test_04_complete_phase_witness(self):
        m=self.m; rp=ket_density(PLUS); rm=ket_density(MINUS)
        self.assertAlmostEqual(trace_distance(rp,rm),1)
        for a,b in itertools.product(gauge_paulis(m,(0,1)),gauge_paulis(m,(2,3))):
            self.assertLess(abs(np.trace((rp-rm)@compress(m,a@b))),1e-12)
        for rho,sign in ((rp,1),(rm,-1)):
            self.assertAlmostEqual(np.trace(rho@compress(m,m['W'])).real,sign)

    def test_05_regional_effect_probabilities(self):
        m=self.m
        states=[m['V']@ket_density(k)@m['V'].conj().T for k in (PLUS,MINUS)]
        for e0,e1,f0,f1 in itertools.product((0,.2,.7,1),repeat=4):
            ea=(e0+e1)*m['eye']/2+(e0-e1)*m['electric'][0]/2
            eb=(f0+f1)*m['eye']/2+(f0-f1)*m['electric'][2]/2
            probs=[np.trace(r@ea@eb).real for r in states]
            self.assertAlmostEqual(probs[0],probs[1])
            self.assertTrue(-1e-12<=probs[0]<=1+1e-12)

    def test_06_shared_flux_is_not_two_independent_bits(self):
        m=self.m
        for p0 in (0,.2,.5,.9,1):
            rho=m['V']@np.diag([p0,1-p0])@m['V'].conj().T
            self.assertAlmostEqual(np.trace(rho@m['electric'][0]@m['electric'][2]).real,1)
            self.assertAlmostEqual(np.trace(rho@m['electric'][0]).real,2*p0-1)
        # Intersection of independent full matrix factors contains only scalars.
        independent_A=[np.kron(x,I) for x in PAULI]
        independent_B=[np.kron(I,x) for x in PAULI]
        self.assertEqual(span_rank(independent_A)+span_rank(independent_B)-span_rank(independent_A+independent_B),1)

    def test_07_independent_preparation_violation(self):
        v=split_embedding(); p=v@v.conj().T
        for a,b in itertools.product(np.linspace(0,1,11),repeat=2):
            rho=np.kron(np.diag([a,1-a]),np.diag([b,1-b]))
            leak=np.trace(rho@(np.eye(4)-p)).real
            self.assertAlmostEqual(leak,a*(1-b)+(1-a)*b)
            if leak<1e-12:
                self.assertTrue((a==b==0) or (a==b==1))

    def test_08_actual_cut_and_boundary_embedding(self):
        v=split_embedding(); va=cycle(2)['V']; vb=cycle(2)['V']
        np.testing.assert_allclose(np.kron(va,vb)@v,self.m['V'],atol=1e-12)
        p=v@v.conj().T
        np.testing.assert_allclose(p,(np.eye(4)+np.kron(Z,Z))/2,atol=1e-12)
        rho=v@ket_density((np.array([1,1j]))/np.sqrt(2))@v.conj().T
        np.testing.assert_allclose(partial_bipartite(rho,2,2,0),I/2,atol=1e-12)
        self.assertAlmostEqual(np.trace(rho@np.kron(X,Y)).real,1)

    def test_09_encoding_preserves_unknown_reference_not_copies(self):
        v=split_embedding(); iso=np.kron(v,I)
        rng=np.random.default_rng(360)
        for _ in range(8):
            a=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
            rho=a@a.conj().T; rho/=np.trace(rho)
            encoded=iso@rho@iso.conj().T
            np.testing.assert_allclose(iso.conj().T@encoded@iso,rho,atol=1e-12)
        encoded=v@ket_density(PLUS)@v.conj().T
        self.assertAlmostEqual(trace_distance(partial_bipartite(encoded,2,2,0),ket_density(PLUS)),.5)

    def test_10_projective_sewing_is_not_free_deterministic_map(self):
        v=split_embedding(); p=v@v.conj().T; q=np.eye(4)-p
        np.testing.assert_allclose(p.conj().T@p+q.conj().T@q,np.eye(4))
        r0=np.diag([1,0,0,0]).astype(complex)
        r1=ket_density(np.kron(PLUS,PLUS))
        o0,s0=projected(r0,p); o1,s1=projected(r1,p)
        om,sm=projected((r0+r1)/2,p)
        self.assertAlmostEqual(s0,1); self.assertAlmostEqual(s1,.5); self.assertAlmostEqual(sm,.75)
        self.assertAlmostEqual(trace_distance(om,(o0+o1)/2),1/(6*np.sqrt(2)))
        np.testing.assert_allclose(p@((r0+r1)/2)@p,(p@r0@p+p@r1@p)/2)

    def test_11_local_endpoint_energy_and_closed_loop(self):
        m=self.m; penalty=sum((m['eye']-g)/2 for g in m['gauss'])
        eig=np.linalg.eigvalsh(penalty)
        self.assertAlmostEqual(min(x for x in eig if x>1e-8),2)
        psi=m['V']@PLUS
        for labels,cost in (([Z,I,I,I],2),([Z,Z,I,I],2),([Z,Z,Z,Z],0)):
            out=kron_all(labels)@psi
            self.assertAlmostEqual(np.vdot(out,penalty@out).real,cost)

    def test_12_local_entropy_does_not_certify_global_phase(self):
        v=split_embedding()
        for p0 in (.1,.2,.5,.9):
            ket=np.array([np.sqrt(p0),1j*np.sqrt(1-p0)])
            pure=v@ket_density(ket)@v.conj().T
            mixed=v@np.diag([p0,1-p0])@v.conj().T
            np.testing.assert_allclose(partial_bipartite(pure,2,2,0),partial_bipartite(mixed,2,2,0))
            self.assertAlmostEqual(entropy(pure),0)
            self.assertAlmostEqual(entropy(mixed),entropy(partial_bipartite(pure,2,2,0)))

    def test_13_all_proper_cycle_cuts_and_bounded_resources(self):
        for n in (3,4,5,6):
            m=cycle(n)
            np.testing.assert_allclose(m['P'],m['V']@m['V'].conj().T,atol=1e-12)
            for size in range(1,n):
                # Internal path constraints leave two flux states per region.
                va=np.column_stack([kron_all([PLUS]*size),kron_all([MINUS]*size)])
                vb=np.column_stack([kron_all([PLUS]*(n-size)),kron_all([MINUS]*(n-size))])
                np.testing.assert_allclose(np.kron(va,vb)@split_embedding(),m['V'],atol=1e-12)
            for electric in m['electric']:
                np.testing.assert_allclose(compress(m,electric),Z,atol=1e-12)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Checks)
    checked=unittest.TextTestRunner(verbosity=2).run(suite)
    if not checked.wasSuccessful():
        raise SystemExit(1)
    result=report()
    result['checks']={'run':checked.testsRun,'failures':0,'errors':0}
    if args.write_results:
        target=Path(__file__).with_name(Path(__file__).stem+'_results.json')
        if target.exists() and json.loads(target.read_text(encoding='utf-8'))!=result:
            raise RuntimeError('Refusing to overwrite a different saved result.')
        if not target.exists():
            target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
