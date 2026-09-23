"""Round 356: finite two-way mean dynamics from one extensive joint unitary.

Fixed finite reduced blocks converge with arbitrary reference. The interaction
is fully connected, not a spatial gravity model. No nonlinear single-copy map
is asserted on arbitrary preparations.
"""
import unittest
import numpy as np
from growing_stream_audit import main

I=np.eye(2,dtype=complex)
X=np.array([[0,1],[1,0]],complex)
Y=np.array([[0,-1j],[1j,0]],complex)
Z=np.diag([1.,-1.]).astype(complex)


def state(mu):
    return (I+np.sqrt(1-mu*mu)*X+mu*Z)/2


def power(matrix,n):
    out=np.array(1.,complex)
    for _ in range(n):
        out=np.kron(out,matrix)
    return out


def sums(n):
    return np.array([n-2*i.bit_count() for i in range(2**n)],float)


def selected_sums(k,l):
    return np.repeat(sums(k),2**l),np.tile(sums(l),2**k)


def kernel(n,k,l,mu_g,mu_m,theta):
    if not (0<=k<=n and 0<=l<=n and n>=1):
        raise ValueError('Selected blocks must lie in the two N-spin sectors.')
    sg,sm=selected_sums(k,l)
    dg=sg[:,None]-sg[None,:]
    dm=sm[:,None]-sm[None,:]
    energy=sg*sm
    direct=np.exp(-1j*theta*(energy[:,None]-energy[None,:])/n)
    from_m=(np.cos(theta*dg/n)-1j*mu_m*np.sin(theta*dg/n))**(n-l)
    from_g=(np.cos(theta*dm/n)-1j*mu_g*np.sin(theta*dm/n))**(n-k)
    return direct*from_m*from_g


def mean_kernel(k,l,mu_g,mu_m,theta):
    sg,sm=selected_sums(k,l)
    e=mu_m*sg+mu_g*sm
    return np.exp(-1j*theta*(e[:,None]-e[None,:]))


def apply(rho,schur,reference=1):
    d=len(schur)
    return (rho.reshape(d,reference,d,reference)*schur[:,None,:,None]).reshape(d*reference,d*reference)


def error_bound(n,k,l,mu_g,mu_m,theta):
    deterministic=3*abs(theta)*k*l/n
    fluctuation=theta**2*(l*l*(n-k)*(1-mu_g**2)+k*k*(n-l)*(1-mu_m**2))/n**2
    return min(1.,deterministic+fluctuation)


def partial(rho,dims,keep):
    out=rho.reshape(tuple(dims)*2)
    labels=list(range(len(dims)))
    for label in reversed(range(len(dims))):
        if label not in keep:
            index=labels.index(label)
            out=np.trace(out,axis1=index,axis2=index+len(labels))
            labels.remove(label)
    size=int(np.prod([dims[j] for j in labels]))
    return out.reshape(size,size)


def full_embedding(rho,n,k,l,mu_g,mu_m,reference=1):
    """Input ordering selected G, selected M, R; append independent environments."""
    env_g=power(state(mu_g),n-k)
    env_m=power(state(mu_m),n-l)
    initial=np.kron(np.kron(rho,env_g),env_m)
    dims=[2]*(k+l)+[reference]+[2]*(2*n-k-l)
    rg=k+l+1
    rm=rg+n-k
    order=list(range(k))+list(range(rg,rm))+list(range(k,k+l))+list(range(rm,len(dims)))+[k+l]
    moved=initial.reshape(dims+dims).transpose(order+[j+len(dims) for j in order])
    return moved.reshape(2**(2*n)*reference,2**(2*n)*reference)


def full_reduced(rho,n,k,l,mu_g,mu_m,theta,reference=1):
    full=full_embedding(rho,n,k,l,mu_g,mu_m,reference)
    sg,sm=selected_sums(n,n)
    phase=np.repeat(np.exp(-1j*theta*sg*sm/n),reference)
    full=phase[:,None]*full*phase.conj()[None,:]
    keep=list(range(k))+list(range(n,n+l))+[2*n]
    return partial(full,[2]*(2*n)+[reference],keep)


def distance(a,b):
    return float(np.sum(abs(np.linalg.eigvalsh(a-b)))/2)


def expectation(rho,op):
    return float(np.trace(rho@op).real)


def iid_reduced(n,k,l,mu_g,mu_m,theta):
    initial=np.kron(power(state(mu_g),k),power(state(mu_m),l))
    return apply(initial,kernel(n,k,l,mu_g,mu_m,theta))


def response(n,mu_self,mu_other,theta):
    gamma=(np.cos(2*theta/n)-1j*mu_other*np.sin(2*theta/n))**n
    return float(-np.sqrt(1-mu_self**2)*gamma.imag)


def variance_y(n,mu_g,mu_m,theta):
    one=iid_reduced(n,1,0,mu_g,mu_m,theta)
    mean=expectation(one,Y)
    if n==1:
        return 1-mean**2
    two=iid_reduced(n,2,0,mu_g,mu_m,theta)
    pair=expectation(two,np.kron(Y,Y))
    return float(1/n+(n-1)*pair/n-mean*mean)


def report():
    mu_g,mu_m,theta=.3,-.4,.6
    initial=np.kron(state(mu_g),state(mu_m))
    target=apply(initial,mean_kernel(1,1,mu_g,mu_m,theta))
    rows=[]
    for n in (10,100,1000,10000):
        out=iid_reduced(n,1,1,mu_g,mu_m,theta)
        rows.append({'N':n,'two_spin_trace_distance':distance(out,target),
                     'arbitrary_reference_upper':error_bound(n,1,1,mu_g,mu_m,theta),
                     'geometry_average_Y':response(n,mu_g,mu_m,theta),
                     'matter_average_Y':response(n,mu_m,mu_g,theta),
                     'geometry_average_Y_variance':variance_y(n,mu_g,mu_m,theta),
                     'variance_upper':(1+4*theta**2*(1-mu_g**2)*(1-mu_m**2))/n})
    cat_rows=[]
    for n in (10,100,1000):
        # Correlated M is half all-up, half all-down; same mu_M=0.
        cat_gamma=np.cos(2*theta)
        product_gamma=np.cos(2*theta/n)**n
        cat_rows.append({'N':n,'cat_single_G_error':float(abs(cat_gamma-1)/2),
                         'product_single_G_error':float(abs(product_gamma-1)/2)})
    return {'round':356,
            'scope':'Specified all-to-all two-sector mean interaction, fixed finite observed blocks and initially independent product environments; finite two-way classical mean feedback and concentration, not spatial geometry, HDA, global-state convergence or GR.',
            'parameters':{'mu_G':mu_g,'mu_M':mu_m,'theta':theta},
            'fixed_time_rows':rows,
            'finite_mean_response':{
                'geometry_Y':float(np.sqrt(1-mu_g**2)*np.sin(2*theta*mu_m)),
                'matter_Y':float(np.sqrt(1-mu_m**2)*np.sin(2*theta*mu_g)),
                'interaction_energy_per_N_in_units_g':mu_g*mu_m},
            'correlated_matter_counterexample':cat_rows,
            'resource_account':{'sector_qubits':'N and N','pairs':'N^2',
                                'coupling_per_pair':'g/N',
                                'operator_norm':'abs(g)*N',
                                'norm_per_total_qubit':'abs(g)/2',
                                'sum_of_incident_absolute_couplings_per_spin':'abs(g)',
                                'spatial_locality_proved':False,
                                'global_state_approximated':False},
            'source':'https://arxiv.org/abs/1804.00455 (related mean-interaction literature; bounds here derived directly)'}


class Checks(unittest.TestCase):
    def test_01_full_spectrum_and_energy_scaling(self):
        for n in (1,2,3):
            sg,sm=selected_sums(n,n)
            self.assertEqual(max(abs(sg*sm/n)),n)
            rho=np.kron(power(state(.3),n),power(state(-.4),n))
            self.assertAlmostEqual(float(np.diag(rho).real@(sg*sm/n))/n,-.12)

    def test_02_schur_kernel_is_a_channel(self):
        for n,k,l in ((2,1,1),(4,2,1),(10,2,2)):
            mat=kernel(n,k,l,.3,-.4,.6)
            np.testing.assert_allclose(mat,mat.conj().T,atol=1e-14)
            np.testing.assert_allclose(np.diag(mat),1.,atol=1e-14)
            self.assertGreater(np.linalg.eigvalsh(mat).min(),-1e-13)

    def test_03_full_joint_and_selected_formula(self):
        for n,k,l in ((1,1,1),(2,1,1),(3,2,1),(3,1,2),(3,1,0)):
            initial=np.kron(power(state(.3),k),power(state(-.4),l))
            direct=full_reduced(initial,n,k,l,.3,-.4,.6)
            predicted=apply(initial,kernel(n,k,l,.3,-.4,.6))
            np.testing.assert_allclose(direct,predicted,atol=2e-15)

    def test_04_entangled_reference_full_joint(self):
        k=l=1; reference=4
        v=np.eye(4).reshape(-1)/2
        rho=np.outer(v,v)
        actual=full_reduced(rho,2,k,l,.3,-.4,.7,reference)
        result=apply(rho,kernel(2,k,l,.3,-.4,.7),reference)
        np.testing.assert_allclose(actual,result,atol=3e-15)

    def test_05_arbitrary_reference_error_bound(self):
        rng=np.random.default_rng(356)
        for n,k,l in ((10,1,1),(30,2,1),(100,2,2),(100,0,2)):
            reference=3; d=2**(k+l)*reference
            for _ in range(3):
                v=rng.normal(size=d)+1j*rng.normal(size=d)
                v/=np.linalg.norm(v);rho=np.outer(v,v.conj())
                a=apply(rho,kernel(n,k,l,.3,-.4,.6),reference)
                b=apply(rho,mean_kernel(k,l,.3,-.4,.6),reference)
                self.assertLessEqual(distance(a,b),error_bound(n,k,l,.3,-.4,.6)+1e-13)

    def test_06_fixed_block_convergence(self):
        rows=report()['fixed_time_rows']
        for row in rows:
            self.assertLess(row['two_spin_trace_distance'],row['arbitrary_reference_upper'])
        for a,b in zip(rows,rows[1:]):
            self.assertGreater(a['two_spin_trace_distance']/b['two_spin_trace_distance'],9.)

    def test_07_both_mean_responses_survive(self):
        r=report()
        last=r['fixed_time_rows'][-1]
        for key,limkey in (('geometry_average_Y','geometry_Y'),('matter_average_Y','matter_Y')):
            self.assertAlmostEqual(last[key],r['finite_mean_response'][limkey],delta=4e-5)
            self.assertGreater(abs(last[key]),.3)

    def test_08_response_agrees_with_reduced_state(self):
        for n in (2,5,30):
            rho=iid_reduced(n,1,1,.3,-.4,.6)
            self.assertAlmostEqual(expectation(rho,np.kron(Y,I)),response(n,.3,-.4,.6))
            self.assertAlmostEqual(expectation(rho,np.kron(I,Y)),response(n,-.4,.3,.6))
            self.assertAlmostEqual(expectation(rho,np.kron(Z,I)),.3)
            self.assertAlmostEqual(expectation(rho,np.kron(I,Z)),-.4)

    def test_09_collective_concentration_with_finite_response(self):
        for n in (2,10,100,1000):
            var=variance_y(n,.3,-.4,.6)
            self.assertGreaterEqual(var,-1e-13)
            self.assertLessEqual(var,(1+4*.6**2*(1-.3**2)*(1-.4**2))/n+1e-13)

    def test_10_correlated_matter_blocks_mean_limit(self):
        rows=report()['correlated_matter_counterexample']
        for r in rows:
            self.assertGreater(r['cat_single_G_error'],.3)
        self.assertLess(rows[-1]['product_single_G_error'],.001)

    def test_11_mean_equations_and_energy(self):
        theta=.6;mu_g=.3;mu_m=-.4;e=1e-5
        def mean_state(mu,other,t):
            phase=np.diag(np.exp(-1j*t*other*np.array([1,-1])))
            return phase@state(mu)@phase.conj().T
        g=mean_state(mu_g,mu_m,theta);m=mean_state(mu_m,mu_g,theta)
        for rho,mu,other in ((g,mu_g,mu_m),(m,mu_m,mu_g)):
            finite=(mean_state(mu,other,theta+e)-mean_state(mu,other,theta-e))/(2*e)
            rhs=-1j*other*(Z@rho-rho@Z)
            np.testing.assert_allclose(finite,rhs,atol=2e-10)
        self.assertAlmostEqual(expectation(g,Z)*expectation(m,Z),mu_g*mu_m)

    def test_12_one_sector_reduction_and_no_nonlinear_unknown_map(self):
        k,l=1,0;n=13;theta=.7;mu=.4
        mat=kernel(n,k,l,.3,mu,theta)
        self.assertAlmostEqual(mat[0,1],(np.cos(2*theta/n)-1j*mu*np.sin(2*theta/n))**n)
        rho1=state(.3);rho2=state(-.8)
        np.testing.assert_allclose(apply(.4*rho1+.6*rho2,mat),
                                   .4*apply(rho1,mat)+.6*apply(rho2,mat),atol=1e-15)


if __name__=='__main__':
    main(__name__,'two_sector_mean_field_audit',report)
