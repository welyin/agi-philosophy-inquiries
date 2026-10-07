"""858: inherited material-cell matching, Dirac reaction, actual-velocity gap.
The star uses the original radial target and node kinetic term. Other original
terms are retained analytically as a finite kappa-independent geometric source.
The small Dirac matrix is a separate algebra check, not the full graph model.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
import numpy as np
import material_cell_canonical_probe as inherited
HERE=Path(__file__).resolve().parent;TARGET=HERE/'material_cell_constraint_closure_results.json'

def arr(x):return np.array([[F(v) for v in row] for row in x],dtype=object)
def inverse(M):
    n=len(M);a=np.concatenate([M.copy(),np.eye(n,dtype=object)],axis=1)
    for i in range(n):
        k=next(j for j in range(i,n) if a[j,i])
        a[[i,k]]=a[[k,i]]
        a[i]=[F(v)/a[i,i] for v in a[i]]
        for j in range(n):
            if j!=i:a[j]=a[j]-a[j,i]*a[i]
    return a[:,n:]
def matstr(M):return [[str(v) for v in row] for row in M]
def canonical(n):
    J=np.zeros((2*n,2*n),dtype=object)
    for i in range(n):J[2*i,2*i+1]=F(1);J[2*i+1,2*i]=F(-1)
    return J

def dirac_check():
    # z=(q1,p1,q2,p2), then r1,r2,pi1,pi2.
    Jz=canonical(2);J=np.zeros((8,8),dtype=object);J[:4,:4]=Jz
    J[4:6,6:8]=np.eye(2,dtype=object);J[6:8,4:6]=-np.eye(2,dtype=object)
    A=arr([[2,1],[1,3]]);L=arr([[1,2,0,1],[0,1,3,-1]])
    Cgrad=np.zeros((4,8),dtype=object);Cgrad[:2,:4]=-L;Cgrad[:2,4:6]=A;Cgrad[2:,6:]=np.eye(2,dtype=object)
    C=Cgrad@J@Cgrad.T;B=L@Jz@L.T;Ai=inverse(A)
    Ci=np.block([[np.zeros((2,2),dtype=object),-Ai.T],[Ai,Ai@B@Ai.T]])
    assert np.array_equal(Ci,inverse(C)) and np.array_equal(C@Ci,np.eye(4,dtype=object))
    JD=J-J@Cgrad.T@Ci@Cgrad@J
    assert np.array_equal(JD[:4,:4],Jz) and np.count_nonzero(Cgrad@JD)==0
    R=Ai@L
    embed=np.zeros((8,4),dtype=object);embed[:4]=np.eye(4,dtype=object);embed[4:6]=R
    assert np.array_equal(embed.T@inverse(J)@embed,inverse(Jz))
    # A genuine nonseparable positive quadratic H for independent evaluation.
    S=np.eye(8,dtype=object)
    for i in range(8):S[i,i]=F(i+2)
    S[0,4]=S[4,0]=F(1,3);S[2,5]=S[5,2]=F(1,4)
    z=np.array([F(1,2),F(2,3),F(-1,4),F(3,5)],dtype=object);full=embed@z
    dH=S@full;dHred=embed.T@dH
    vector=(JD@dH)[:4];expected=Jz@dHred;frozen=Jz@dH[:4]
    assert np.array_equal(vector,expected) and np.any(vector!=frozen)
    assert np.array_equal(vector-frozen,Jz@R.T@dH[4:6])
    return dict(constraint_matter_bracket_B=matstr(B),constraint_inverse=matstr(Ci),
        reduced_symplectic_form_is_original=True,full_dirac_and_substituted_H_vectors_equal=True,
        omitted_reaction_vector=[str(v) for v in vector-frozen],constraint_tangency_residual='0')

def star(kappa,beta=F(0),scale=F(1)):
    d=F(1,4);w0=d**3;points=[(0,0,0),(d,0,0),(-d,0,0),(0,d,0),(0,-d,0),(0,0,d),(0,0,-d)]
    fields=[];T=F(0);dbeta=F(0);leading=F(0)
    for x,y,z in points:
        K=inherited.K(F(1),F(1,2)+x)
        v=[kappa+(1+beta)*y/kappa,2+z]
        p=[w0*a for a in inherited.mv(K,v)]
        vel=[a/(w0*scale) for a in inherited.mv(inherited.inv(K),p)]
        assert vel==[a/scale for a in v]
        kinetic=sum(p[i]*vel[i] for i in range(2))/2
        T+=kinetic
        dbeta+=w0*sum(K[0][j]*v[j] for j in range(2))*y/kappa/scale
        leading+=w0*K[0][0]/2
        # gamma=scale^(2/3) I. Its spatial s-gradient part is constant on this
        # stencil, so it cancels in every spatial B difference. Only the
        # exact momentum-dependent part is needed for dB.
        fields.append(dict(A=-vel[0]**2,Bvariable=-vel[1]**2,s=F(1,2)+x))
    sx=(fields[1]['s']-fields[2]['s'])/(2*d)
    Ay=(fields[3]['A']-fields[4]['A'])/(2*d)
    Bz=(fields[5]['Bvariable']-fields[6]['Bvariable'])/(2*d)
    density=8*sx*Ay*Bz/scale
    material=1/density;weight=w0*scale
    assert material==w0*scale**5/(1+beta)
    return dict(T=T,partial_beta_T=dbeta,leading_coefficient=leading,
        old_weight=weight,material_weight=material,ratio=material/weight)

def run():
    assert inherited.run()==json.loads(inherited.TARGET.read_text('utf-8'))
    rows=[]
    for k in (F(1),F(2),F(4),F(8)):
        d=star(k);T=d['T'];slope=d['partial_beta_T'];assert d['ratio']==1
        reduced=slope-T/4
        assert reduced<0 and slope>0
        rows.append(dict(kappa=str(k),seven_node_kinetic_energy=str(T),
            old_frozen_geometry_beta_response=str(slope),node_geometric_source=str(-T),
            reduced_node_beta_response=str(reduced),omitted_reaction=str(-T/4),
            rho_beta='1/4',G_velocity_difference_from_node_part=str(T/4)))
    roots=[]
    for scale in (F(7,8),F(1),F(9,8)):
        beta=scale**4-1;d=star(F(2),beta,scale);assert d['ratio']==1
        rho_beta=1/(4*(1+beta))
        roots.append(dict(scale=str(scale),beta=str(beta),old_weight=str(d['old_weight']),
            actual_old_flow_material_weight=str(d['material_weight']),
            matching_constraint_derivative_in_log_scale='4',rho_beta=str(rho_beta)))
    high=star(F(1024));lead=high['leading_coefficient']
    return dict(round=858,fresh_test_groups=1,all_checks_passed=True,
        previous_working_probe_reproduced=True,dirac_algebra=dirac_check(),
        original_radial_node_checks=rows,exact_implicit_matching_roots=roots,
        kappa_squared_kinetic_coefficient=str(lead),large_kappa_kinetic_divided_by_kappa_squared=float(high['T']/F(1024)**2),
        source_of_all_other_fixed_graph_terms='finite kappa-independent R_rest, retained symbolically; not set to zero in the full-model argument',
        argument_scope='A conditional second-class closure of old-flow material volumes exists where its geometric Jacobian is invertible. The full reduced Hamiltonian contains geometric-source reaction terms. Exact old radial data show that a one-cell matching constraint is locally regular but cannot preserve all old material velocities on its whole surface. This tests a specified proposed closure, not all internal geometry models.',
        multi_cell_constraint_rank_proved=False,self_consistent_new_velocity_reference_solved=False,
        quantum_reduction_or_Einstein_constraints_proved=False,original_graph_to_continuum_proved=False,full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:r[k] for k in ('round','all_checks_passed','kappa_squared_kinetic_coefficient','original_radial_node_checks')},ensure_ascii=False,indent=2))
