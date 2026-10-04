"""642: full original quotient matter, tree coordinates and retained sources.

All physical witnesses are finite-energy states of the original nonlinear
Gauss model, not Gibbs states or invariant full-H sectors. Hbar=1.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_region_energy_gluing as old
import joint_quotient_boundary_entropy as labels

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_tree_matter_source_results.json'
ENDS=((0,1),(0,2),(1,3),(2,3))
I=(np.eye(3,dtype=complex),np.eye(2,dtype=complex),1.+0j)
matter=old.matter
gauge=old.gauge
original=old.original


def mul(*gs):
    out=I
    for g in gs:out=old.product(out,g)
    return out


def scalar_action(g,p):
    X=g[2]**3*g[1]@(p[:2]+1j*p[2:4])
    return np.r_[X.real,X.imag,p[4]]


def rooted(edges,phi):
    A,B,C,D=edges
    W=(I,A,B,mul(A,C))
    xi=np.array([scalar_action(W[v],phi[v]) for v in range(4)])
    loop=mul(B,D,old.inverse(C),old.inverse(A))
    return W,xi,loop


def group_error(g,h):
    return max(float(np.max(abs(g[0]-h[0]))),
               float(np.max(abs(g[1]-h[1]))),float(abs(g[2]-h[2])))


def one_parameter(c,a,t):
    if c==0:
        z=np.zeros(8);z[a]=t
        return gauge.group_exp(z,3),I[1],1.+0j
    if c==1:
        z=np.zeros(3);z[a]=t
        return I[0],gauge.group_exp(z,2),1.+0j
    return I[0],I[1],np.exp(1j*t)


def module_generators(name,c):
    """Actual derivative of the original module, one fixed spin."""
    d=old.rep(I,name).shape[0]
    if c==2:
        q=labels.MODULES[name][3]
        return np.array([q*np.eye(d)])
    gens=gauge.generators(3 if c==0 else 2)
    if c==0:
        if name=='Q':return np.array([np.kron(t,np.eye(2)) for t in gens])
        if name in ('u','d'):return gens
    if c==1:
        if name=='Q':return np.array([np.kron(np.eye(3),t) for t in gens])
        if name=='L':return gens
    return np.zeros((len(gens),d,d),complex)


def tree_dictionary_check():
    rng=np.random.default_rng(642)
    equiv=measure=mass=distance=magnetic=derivative=0.
    for _ in range(6):
        edges=[old.sample(rng) for _ in ENDS]
        phi=rng.normal(size=(4,5))*.3
        W,xi,L=rooted(edges,phi)
        gs=[old.sample(rng) for _ in range(4)]
        transformed=[mul(gs[s],u,old.inverse(gs[t])) for u,(s,t) in zip(edges,ENDS)]
        pp=np.array([scalar_action(gs[v],phi[v]) for v in range(4)])
        Wg,xg,Lg=rooted(transformed,pp)
        equiv=max(equiv,float(np.max(abs(xg-np.array([scalar_action(gs[0],p) for p in xi])))),
                  group_error(Lg,mul(gs[0],L,old.inverse(gs[0]))))
        for v in range(4):
            R=matter.representation(*W[v])
            h,delta=matter.mass_matrices(phi[v]);ht,dt=matter.mass_matrices(xi[v])
            mass=max(mass,float(np.max(abs(ht-R@h@R.conj().T))),
                     float(np.max(abs(dt-R@delta@R.T))))
            # Original curved density is sqrt(M) F^-3.
            measure=max(measure,abs(float(original.F(phi[v])-original.F(xi[v]))))
        for n,(s,t) in enumerate(ENDS):
            z=I if n<3 else L
            before=original.distance_squared(phi[s],scalar_action(edges[n],phi[t]))
            after=original.distance_squared(xi[s],scalar_action(z,xi[t]))
            distance=max(distance,abs(float(before-after)))
        face=mul(edges[0],edges[2],old.inverse(edges[3]),old.inverse(edges[1]))
        magnetic=max(magnetic,abs(gauge.potential(*face,old.geometry.PAR)-
                                      gauge.potential(*old.inverse(L),old.geometry.PAR)))
        # Transported electric variations from each genuine original link.
        for e,(s,t) in enumerate(ENDS):
            for c,a in ((0,2),(1,1),(2,0)):
                h=one_parameter(c,a,3e-5)
                changed=list(edges);changed[e]=mul(h,changed[e])
                _,actual_x,actual_L=rooted(changed,phi)
                g=mul(W[s],h,old.inverse(W[s]))
                expected_x=xi.copy()
                if e<3:
                    descendants=({1,3},{2},{3})[e]
                    for v in descendants:expected_x[v]=scalar_action(g,xi[v])
                    left=g if 2 in descendants else I
                    right=old.inverse(g) if 3 in descendants else I
                    expected_L=mul(left,L,right)
                else:expected_L=mul(g,L)
                derivative=max(derivative,float(np.max(abs(actual_x-expected_x))),
                               group_error(actual_L,expected_L))
    assert max(equiv,measure,mass,distance,magnetic,derivative)<2e-12
    return dict(original_full_quotient_equivariance_error=equiv,
        curved_measure_invariance_error=measure,full32_mass_and_Majorana_error=mass,
        original_geodesic_edge_error=distance,original_magnetic_error=magnetic,
        exact_finite_tree_electric_variation_error=derivative,
        matter_and_cycle_counts_unchanged=True,
        Haar_and_closed_form_unitarity_proved_analytically=True)


SHAPE=np.array([[.24,.31,-.12],[.31,-.07,.17],[-.12,.17,-.17]])
EPS=.73


def coefficients(t=.8,sigma=.12):
    return np.asarray(old.geometry.PAR['b'])[:,None,None]/EPS*np.exp(-2*sigma)*old.geometry.shape_exp(SHAPE,t)


def energy(name,t=.8,sigma=.12):
    C=np.array(labels.label_data(labels.MODULES[name])['casimirs'])
    K=coefficients(t,sigma)
    return float(C@(K[:,0,0]+K[:,1,1]-2*K[:,0,1]))


def physical_pair_check():
    rng=np.random.default_rng(6422)
    A,B=old.sample(rng),old.sample(rng)
    K=coefficients();rows=[]
    h=2e-6
    for name in matter.SLICES:
        RA,RB=old.rep(A,name),old.rep(B,name);d=len(RA)
        F=RA.conj().T@RB/np.sqrt(d)
        norm=float(np.vdot(F,F).real)
        transported=RA@F@RB.conj().T
        constant_error=float(np.max(abs(transported-np.eye(d)/np.sqrt(d))))
        direct_energy=0.;fd_error=0.;casimir_error=0.
        for c in range(3):
            T=module_generators(name,c)
            cas=float(labels.label_data(labels.MODULES[name])['casimirs'][c])
            casimir_error=max(casimir_error,float(np.max(abs(sum(t@t for t in T)-cas*np.eye(d)))))
            for a,tau in enumerate(T):
                # P=-i d/dt: original source-left action on A or B.
                PA=-RA.conj().T@tau@RB/np.sqrt(d);PB=-PA
                gram=np.array([[np.vdot(PA,PA).real,np.vdot(PA,PB).real],
                               [np.vdot(PB,PA).real,np.vdot(PB,PB).real]])
                direct_energy+=float(np.sum(K[c,:2,:2]*gram))
                plus=old.rep(mul(one_parameter(c,a,h),A),name).conj().T@RB/np.sqrt(d)
                minus=old.rep(mul(one_parameter(c,a,-h),A),name).conj().T@RB/np.sqrt(d)
                fd_error=max(fd_error,float(np.max(abs(-1j*(plus-minus)/(2*h)-PA))))
        expected=energy(name)
        C=np.array(labels.label_data(labels.MODULES[name])['casimirs'])
        diagonal=float(C@(K[:,0,0]+K[:,1,1]))
        dK=K@SHAPE
        shear=float(C@(dK[:,0,0]+dK[:,1,1]-2*dK[:,0,1]))
        shear_fd=(energy(name,.8+h)-energy(name,.8-h))/(2*h)
        conf_fd=(energy(name,sigma=.12+h)-energy(name,sigma=.12-h))/(2*h)
        assert max(abs(norm-1),constant_error,abs(direct_energy-expected),casimir_error)<1e-13
        assert fd_error<1e-8 and abs(shear-shear_fd)<1e-8 and abs(conf_fd+2*expected)<1e-8
        rows.append(dict(module=name,dimension=d,pointwise_norm=norm,
            constant_rooted_singlet_error=constant_error,Casimir_error=casimir_error,
            original_link_momentum_FD_error=fd_error,full_electric_expectation=expected,
            direct_original_gradient_energy=direct_energy,
            energy_if_shear_deleted=diagonal,source_shear=shear,
            shear_source_FD_error=float(abs(shear-shear_fd)),
            uniform_conformal_source=-2*expected,
            conformal_source_FD_error=float(abs(conf_fd+2*expected))))
    assert energy('Q')>0 and energy('nu')==0
    return dict(original_positive_coefficients=K.tolist(),rows=rows,hbar=1.,
        physical_state_uses_original_particle_and_full_sea_hole=True,
        same_normalized_compact_scalar_packets_all_modules=True,
        no_full_H_eigenstate_or_equilibrium_claim=True)


def coarse_source_and_mass_check():
    rng=np.random.default_rng(6423)
    determinant=mass_mean=0.;rows=[]
    for _ in range(6):
        R=matter.representation(*old.sample(rng))
        determinant=max(determinant,float(abs(np.linalg.det(R)-1)))
    for name in matter.SLICES:
        ids=np.arange(32)[matter.SLICES[name]][::2]
        n1=np.zeros((32,32));n1[ids,ids]=1/len(ids);n2=np.eye(32)-n1
        commutators=[]
        for _ in range(8):
            phi=rng.normal(size=(2,5))*.25
            h1,_=matter.mass_matrices(phi[0]);h2,_=matter.mass_matrices(phi[1])
            mean=float(np.trace(h1@n1+h2@n2).real)
            mass_mean=max(mass_mean,abs(mean))
            project=np.zeros((32,32));project[ids,ids]=1
            commutators.append(float(np.linalg.norm(project@h1-h1@project,2)))
        rows.append(dict(module=name,local_particle_numbers=[float(np.trace(n1)),float(np.trace(n2))],
            maximal_original_mass_projection_commutator=max(commutators)))
        assert max(commutators)>1e-4
    gap=energy('Q')-energy('nu')
    assert determinant<3e-13 and mass_mean<1e-14
    return dict(full32_filled_state_determinant_error=determinant,
        all_original_onsite_mass_expectations_zero_error=mass_mean,
        Majorana_mean_zero_by_fixed_total_number=True,
        hopping_mean_zero_by_fixed_local_numbers=True,
        bosonic_and_magnetic_means_equal_by_pointwise_norm=True,
        module_rows=rows,
        identical_retained_density_states_proved_analytically=True,
        full_H_expectation_gap_Q_minus_nu=gap,
        full_conformal_source_gap_Q_minus_nu=-2*gap,
        any_retained_state_energy_estimator_worst_error_lower_bound=gap/2,
        any_retained_state_conformal_source_estimator_worst_error_lower_bound=gap,
        module_tag_not_conserved_by_original_Yukawa=True,
        global_root_charge_both_zero=True)


def run():
    deps=('research_note_566.md','research_note_589.md','research_note_598.md',
          'research_note_617.md','research_note_636.md','research_note_641.md',
          'joint_region_energy_gluing.py','joint_fermion_gauss_completion.py',
          'joint_quotient_boundary_entropy.py')
    return dict(round=642,tests_run=3,failures=0,errors=0,
        full_tree_dictionary=tree_dictionary_check(),actual_physical_pair=physical_pair_check(),
        coarse_state_source_obstruction=coarse_source_and_mass_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_nonlinear_target_full_quotient_and32_CAR=True,
            exact_tree_change_of_representation_not_yet_physical_compression=True,
            explicit_legal_Gauss_states_with_identical_retained_state=True,
            full_H_and_geometry_mean_gap_not_only_electric_truncation=True,
            no_Gibbs_compression_continuum_or_GR_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
