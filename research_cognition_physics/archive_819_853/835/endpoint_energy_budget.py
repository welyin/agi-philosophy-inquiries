"""835 working: original-material energy variance of a controller endpoint.

Native sparse CAR checks validate local moment identities. The 51-node code
states are handled by the known degree-four compression, not huge matrices.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'endpoint_energy_budget_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime


def ensemble(nodes,occupied,polarized=None):
    mixed=[v for v in occupied if v!=polarized]
    for bits in itertools.product((0,1),repeat=len(mixed)):
        state=sum(1<<(32*v+30+b) for v,b in zip(mixed,bits))
        if polarized is not None:state|=1<<(32*polarized+30)
        yield {state:1.},2.**(-len(mixed))


def sparse_checks(old):
    rng=np.random.default_rng(835);rows=[]
    for nodes in (2,3,5):
        edges=[(v,v+1) for v in range(nodes-1)]
        initial=set(range(max(1,nodes//2)))
        added=set(range(nodes))-initial
        for trial in range(3):
            x=rng.normal(size=(nodes,5))*.6
            if trial==2:x[:,4]+=4
            h=np.zeros((32*nodes,32*nodes),complex);d=np.zeros_like(h)
            for v in range(nodes):
                sl=slice(32*v,32*(v+1));h[sl,sl],d[sl,sl]=old.mass_x(x[v])
            hop=np.zeros_like(h);edge_strength={}
            for v,w in edges:
                c=old.matter.gauge.group_exp(rng.normal(size=8)*.2,3)
                weak=old.matter.gauge.group_exp(rng.normal(size=3)*.2,2)
                rep=old.matter.representation(c,weak,np.exp(.13j))
                # Nontrivial spin unitary, retaining every original species.
                spin=old.matter.gauge.group_exp(rng.normal(size=3)*.2,2)
                coefficient=(.23-.11j)*(1+.1*v)
                edge=coefficient*rep@np.kron(np.eye(16),spin)
                a=slice(32*v,32*(v+1));b=slice(32*w,32*(w+1))
                hop[a,b]=edge;hop[b,a]=edge.conj().T
                edge_strength[(v,w)]=.5*float(np.linalg.norm(edge[30:32,30:32])**2)
            moments=[]
            for occ,polarized in ((initial,None),(set(range(nodes)),nodes-1)):
                mean=second=0.
                for state,weight in ensemble(nodes,occ,polarized):
                    image=old.car.quadratic(state,h+hop,d)
                    mean+=weight*float(old.dot(state,image).real)
                    second+=weight*float(old.dot(image,image).real)
                assert abs(mean)<1e-14
                moments.append(second)
            dirac=sum(abs(old.matter.Y['nu'])**2*float(x[v,:4]@x[v,:4]) for v in added)
            majorana=-sum(abs(old.matter.Y['s'])**2*x[v,4]**2 for v in added)
            hopping=sum(value for (v,w),value in edge_strength.items() if v in added and w in added)
            residual=abs(moments[1]-moments[0]-dirac-majorana-hopping)
            assert residual<3e-12
            rows.append(dict(nodes=nodes,full_CAR_modes=32*nodes,trial=trial,
                             native_mass_and_hopping_second_moment_difference=moments[1]-moments[0],
                             dirac_increment=dirac,majorana_increment=float(majorana),
                             hopping_increment=hopping,identity_residual=float(residual)))
    return rows


def packet(order,old):
    # Smooth compact strict-Gauss packet: r_H < .7, |x5-4| < .3.
    u,wu=np.polynomial.legendre.leggauss(order)
    r=.7*(u+1)/2;s=4+.3*u
    rr,ss=np.meshgrid(r,s,indexing='ij')
    w=(.35*wu[:,None])*(.3*wu[None,:])*rr**3/np.sqrt(1+(rr**2+ss**2)/6)
    w*=np.exp(-2/(1-(rr/.7)**2)-2/(1-((ss-4)/.3)**2))
    w/=w.sum()
    r2=float(np.sum(w*rr**2));s2=float(np.sum(w*ss**2))
    # 51-node chain; first 5 nodes initially carry the code, 46 are empty.
    # The final tag is spin-up at node 50; other output nodes have mixed
    # one-particle two-/four-point moments from the ten code blocks.
    edges_added=45;edge_strength=abs(.23-.11j)**2
    delta=46*(abs(old.matter.Y['nu'])**2*r2-abs(old.matter.Y['s'])**2*s2)+edges_added*edge_strength
    support_upper=46*(abs(old.matter.Y['nu'])**2*.7**2-abs(old.matter.Y['s'])**2*3.7**2)+edges_added*edge_strength
    ready_lower=46*abs(old.matter.Y['s'])**2*s2
    assert delta<support_upper<0
    return dict(Higgs_radius_second_moment=r2,singlet_second_moment=s2,
                native_Ynu_abs_squared=float(abs(old.matter.Y['nu'])**2),
                native_Ys_abs_squared=float(abs(old.matter.Y['s'])**2),
                common_mean_energy_difference=0,
                full_H_squared_endpoint_difference=float(delta),
                rigorous_support_upper_bound=float(support_upper),
                ready_resource_vacuum_loss_t2_Majorana_lower_coefficient=float(ready_lower),
                chain_nodes=51,initially_occupied_code_nodes=5,
                initially_empty_sterile_nodes=46,
                result_is_in_native_model_units=True)


def run():
    with ResearchRuntime(Layout()).installed():
        import joint_record_mass_feedback as old
        rows=sparse_checks(old)
        fine,coarse=packet(80,old),packet(48,old)
    error=max(abs(fine[k]-coarse[k]) for k in fine if isinstance(fine[k],float))
    assert error<1e-7
    return dict(round=835,status='working_not_formal',formal_test_groups_added=0,
                all_working_checks_passed=True,sparse_original_coefficient_checks=rows,
                organized_input_endpoint_packet=fine,quadrature_change_48_to_80=error,
                full_boson_H_retained_in_analytic_cancellation=True,
                strict_counterexample_requires_same_final_boson_and_other_material_state=True,
                old_round746_input_obstruction_not_reused_as_a_new_input_proof=True,
                all_compensated_or_autonomous_recovery_protocols_excluded=False,
                scope='A specified organized-code input with empty internal resource modes and the 834 zero-syndrome output can have identical mean energy but different full-H second moments when the rest of the state is held fixed. Original Majorana coupling also disturbs empty resources at order t^2. Changed environment, compensating resources and other endpoint contracts remain open.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(dict(status=r['status'],all_checks_passed=True,
        packet=r['organized_input_endpoint_packet'],quadrature_change=r['quadrature_change_48_to_80'],
        maximum_sparse_residual=max(z['identity_residual'] for z in r['sparse_original_coefficient_checks'])),ensure_ascii=False,indent=2))
