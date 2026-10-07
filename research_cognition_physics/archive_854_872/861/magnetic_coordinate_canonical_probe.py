"""861 working: exact symplectic-potential test for a magnetic coordinate.
A periodic difference functional with noncommuting real SU(3) generators
checks the derivative momentum pullback. It is not a gauge-invariant lattice
GR model and does not solve continuum constraints.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'magnetic_coordinate_canonical_probe_results.json'

def zero(shape):return np.full(shape,F(0),dtype=object)
def comm(a,b):return a@b-b@a
def inner(a,b):return -2*sum((a@b)[i,i] for i in range(3))
def generators():
    out=[]
    for i,j in ((0,1),(1,2),(0,2)):
        a=zero((3,3));a[i,j]=F(1,2);a[j,i]=-F(1,2);out.append(a)
    return out

def run():
    n=5;gen=generators();A=zero((n,3,3,3));E=A.copy()
    for v in range(n):
        for i in range(3):
            A[v,i]=F(v+i+2,11)*gen[i]+F((v+2*i)%3-1,17)*gen[(i+1)%3]
            E[v,i]=F(v+2*i+1,19)*gen[(i+2)%3]
    def partial(z,i):return (np.roll(z,-1,axis=0)-np.roll(z,1,axis=0))/2 if i==0 else zero(z.shape)
    def covariant(z,i):return partial(z,i)+np.array([comm(A[v,i],z[v]) for v in range(n)],dtype=object)
    curvature=zero((n,3,3,3,3))
    for i in range(3):
        for j in range(3):curvature[:,i,j]=partial(A[:,j],i)-partial(A[:,i],j)+np.array([comm(A[v,i],A[v,j]) for v in range(n)],dtype=object)
    b=np.array([sum(inner(curvature[v,i,j],curvature[v,i,j]) for i in range(3) for j in range(i+1,3)) for v in range(n)],dtype=object)
    assert min(b)>0
    M=np.array([F(v+4,7) for v in range(n)],dtype=object)
    pi_r=np.array([F(v+1,11) for v in range(n)],dtype=object)
    pi_M=-F(3,4)*pi_r/M
    weight=pi_r/b
    shift=zero(E.shape)
    for j in range(3):
        for i in range(3):shift[:,j]-=F(3,2)*covariant(weight[:,None,None]*curvature[:,i,j],i)
    new_E=E+shift
    def pairing(a,z):return sum(inner(a[v,i],z[v,i]) for v in range(n) for i in range(3))
    rows=[]
    for trial in range(4):
        dA=zero(A.shape)
        for v in range(n):
            for i in range(3):dA[v,i]=F(v+i+trial-2,13)*gen[(i+trial)%3]+F(v-trial,23)*gen[(i+1)%3]
        db=zero((n,))
        for i in range(3):
            for j in range(i+1,3):
                dF=covariant(dA[:,j],i)-covariant(dA[:,i],j)
                db+=np.array([2*inner(curvature[v,i,j],dF[v]) for v in range(n)],dtype=object)
        dM=np.array([F(v-trial,29) for v in range(n)],dtype=object)
        dr=F(3,4)*(db/b-dM/M)
        theta_old=sum(pi_r*dr)+pairing(E,dA)
        theta_new=sum(pi_M*dM)+pairing(new_E,dA)
        omitted=sum(pi_M*dM)+pairing(E,dA)
        assert theta_old==theta_new and theta_old!=omitted
        rows.append(dict(trial=trial,symplectic_potential_residual='0',omitted_electric_momentum_shift_defect=str(theta_old-omitted)))
    energy=pairing(E,E)/2
    transported=pairing(new_E-shift,new_E-shift)/2
    wrong=pairing(new_E,new_E)/2
    assert energy==transported and energy!=wrong
    return dict(kind='round_861_working_canonical_probe',formal_reports=860,new_numbered_scientific_groups=0,
        all_checks_passed=True,noncommuting_color_matrices=True,positive_magnetic_norms=[str(x) for x in b],rows=rows,
        electric_energy_transport_residual='0',untransported_electric_energy_error=str(wrong-energy),
        metric_unimodular_variations_in_numerics=False,continuum_gauge_and_Einstein_constraints_solved=False,
        exact_quantum_or_graph_equivalence_proved=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
