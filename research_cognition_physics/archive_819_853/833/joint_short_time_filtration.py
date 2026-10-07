"""833 working: joint source/record short-time filtration.

Exact Gaussian-integer differential/Clifford algebra tests the degree bounds.
A 16-dimensional correctable-error representation tests the Taylor recovery
bound. Neither calibration replaces the original graph Hamiltonian.
"""
from pathlib import Path
import argparse,itertools,json,math,sys
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_short_time_filtration_results.json'
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code
UNITS=(1,1j,-1,-1j)


def term(power,derivatives,pauli,coefficient=1):
    x,z,p=pauli
    return {(power,derivatives,x,z):coefficient*UNITS[p]}


def add(*operators):
    out={}
    for a in operators:
        for key,value in a.items():out[key]=out.get(key,0)+value
    return {key:value for key,value in out.items() if value!=0}


def times(a,b):
    out={}
    for (pa,da,xa,za),ca in a.items():
        for (pb,db,xb,zb),cb in b.items():
            sign=(-1)**((za&xb).bit_count())
            for j in range(min(da,pb)+1):
                factor=math.comb(da,j)*math.factorial(pb)//math.factorial(pb-j)
                key=(pa+pb-j,da+db-j,xa^xb,za^zb)
                out[key]=out.get(key,0)+ca*cb*sign*factor
    return {key:value for key,value in out.items() if value!=0}


def comm(a,b):return add(times(a,b),{key:-value for key,value in times(b,a).items()})


def filtration_checks():
    gamma,_,_,_,compress,_,_=code.code_data()
    es=[code.mul((0,0,1),code.mul(gamma[i],gamma[j])) for i,j in ((0,1),(4,5),(12,14))]
    degrees={}
    for k in range(4):
        for subset in itertools.combinations(es,k):
            p=code.product(subset);degrees[p[:2]]=2*k
    scalar=add(term(0,2,code.I,-1),term(4,0,code.I))
    cases=[('quadratic_source',add(scalar,term(1,0,es[1]),term(1,0,es[2])),term(2,0,es[0]),4,2),
           ('source_pair',add(scalar,term(1,0,es[2])),term(2,0,code.mul(es[0],es[1])),2,4),
           ('scalar_readout',add(scalar,*[term(1,0,e) for e in es]),term(3,0,code.I),6,0)]
    results=[]
    for name,h,obs,last,weight in cases:
        rows=[]
        for n in range(last+1):
            assert all(complex(v).real.is_integer() and complex(v).imag.is_integer() for v in obs.values())
            # Exact arithmetic: all real/imaginary coefficients are integers
            # comfortably below the exact binary64 integer limit.
            assert max([abs(v) for v in obs.values()]+[0])<2**40
            logic=[];maximum=0
            for (power,deriv,x,z),value in obs.items():
                degree=degrees[(x,z)];maximum=max(maximum,degree)
                assert degree+deriv<=weight+n
                if compress((x,z,0))[0]=='logical':
                    logic.append(dict(power=power,derivatives=deriv,
                                      real=int(complex(value).real),imaginary=int(complex(value).imag)))
            if n<last:assert not logic
            else:assert logic and maximum==6
            rows.append(dict(commutators=n,terms=len(obs),maximum_Majorana_degree=maximum,
                             logical_terms=logic))
            obs=comm(h,obs)
        results.append(dict(observable=name,initial_weight=weight,first_possible_logical_order=last,rows=rows))
    # Kinetic momentum coupled to a quadratic fermion violates the weight-two
    # Hamiltonian condition and can expose a source pair in one commutator.
    changed_h=term(0,1,es[2],-1j);pair=term(1,0,code.mul(es[0],es[1]))
    first=comm(changed_h,pair)
    assert any(compress((x,z,0))[0]=='logical' for (_,_,x,z) in first)
    return dict(cases=results,added_fermion_momentum_coupling_breaks_pair_first_order_blindness=True,
                coefficients_are_exact_Gaussian_integers=True,
                sharpness_examples_are_structural_not_original_Yukawa_coefficients=True)


def finite_recovery_check():
    i=np.eye(2);x=np.array([[0.,1.],[1.,0.]]);z=np.diag([1.,-1.])
    x1=np.kron(x,i);x2=np.kron(i,x)
    # Logical, syndrome(2 bits), environment. E1 E2 E3 = logical Z.
    h=np.kron(i,np.kron(x1,z))+np.kron(i,np.kron(x2,x))+.7*np.kron(z,np.kron(x1@x2,i))
    win=np.zeros((2,4,2,2))
    for logical in range(2):win[logical,0,0,logical]=1.
    win=win.reshape(16,2)
    h2norm=float(np.linalg.norm(h@h@win,2))
    eigen,vec=np.linalg.eigh(h)
    bell=np.eye(2).reshape(-1,order='F')/np.sqrt(2)
    rows=[]
    for t in (.03,.06,.12,.24):
        evolved=(vec@(np.exp(-1j*t*eigen)[:,None]*(vec.conj().T@win))).reshape(2,4,2,2)
        kraus=[]
        for syndrome in range(4):
            correction=z if syndrome==3 else i
            for env in range(2):kraus.append(correction@evolved[:,syndrome,env,:])
        j=sum(np.outer(k.reshape(-1,order='F'),k.reshape(-1,order='F').conj())/2 for k in kraus)
        defect=float(1-(bell.conj()@j@bell).real)
        bound=t**4*h2norm**2/4
        assert -1e-14<=defect<=bound+1e-13
        assert np.linalg.norm(sum(k.conj().T@k for k in kraus)-i)<1e-13
        rows.append(dict(time=t,entangled_input_infidelity=defect,Taylor_upper_bound=bound,
                         infidelity_over_t4=defect/t**4))
    return dict(H2_input_operator_norm=h2norm,rows=rows,
                recovery_is_ideal_syndrome_CP_not_native_decoder=True,
                finite_matrix_not_original_boson_dynamics=True)


def run():
    return dict(round=833,status='working_not_formal',formal_test_groups_added=0,
                all_working_checks_passed=True,source_filtration=filtration_checks(),
                actual_time_recovery_calibration=finite_recovery_check(),
                scope='For a scalar boson differential operator of order two plus coordinate-multiplication quadratic fermion couplings, commutation raises derivative-plus-Clifford degree by at most one. The code hides source means through order three, source pairs through order one, and scalar readouts through order five. Genuine time remainders require the declared source domains; graph-bounded parameter sources and bounded readouts admit them. Ideal recovery has a fourth-order infidelity bound, not an autonomous decoder or a fourth-order diamond bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
