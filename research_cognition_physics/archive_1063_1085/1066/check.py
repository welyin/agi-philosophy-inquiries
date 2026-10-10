"""Finite diagnostics for round 1066. Infinite topology is proved, not sampled."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),path
        for k in a: compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f'{path}[{i}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,abs_tol=1e-10,rel_tol=1e-9),(path,a,b)
    else:assert a==b,(path,a,b)

def principal_turn(z): return float(np.angle(z)/(2*math.pi))

def matrix_log(U):
    values,vectors=np.linalg.eig(U)
    return vectors@np.diag(1j*np.angle(values))@np.linalg.inv(vectors)

def matrix_exp(X):
    values,vectors=np.linalg.eig(X)
    return vectors@np.diag(np.exp(values))@np.linalg.inv(vectors)

def run():
    g=np.exp(2j*math.pi*.4)
    square_angle=principal_turn(g*g)
    chosen_root=np.exp(1j*math.pi*square_angle)
    assert abs(square_angle+.2)<1e-14
    assert abs(chosen_root-g)>1.99
    large={'g_turn':.4,'square_turn':square_angle,
           'root_of_square_turn':principal_turn(chosen_root),
           'chordal_reverse_defect':float(abs(chosen_root-g)),
           'original_U_reverse_consistency':False}

    quotients=[]; max_error=0.
    # Every projection is a known circle power map. A quotient may have a
    # wrapped endpoint while its root chain retains a nonzero lifted generator.
    for turn in (.4,-.31,.17):
        for power in (1,3,5,11):
            n0=0
            while abs(power*turn)/(2**n0)>=.24:n0+=1
            scaled=[]
            for n in range(n0,n0+5):
                source=np.exp(2j*math.pi*turn/(2**n))
                scaled.append((2**n)*principal_turn(source**power))
            err=max(abs(x-power*turn) for x in scaled)
            assert err<1e-10
            max_error=max(max_error,err)
            derror=0.
            for depth in (1,3,5):
                for numerator in (0,1,2**depth-1,2**depth):
                    t=numerator/(2**depth)
                    exp_path=np.exp(2j*math.pi*t*scaled[-1])
                    original_root=np.exp(2j*math.pi*turn/(2**depth))
                    from_root=original_root**(numerator*power)
                    derror=max(derror,float(abs(exp_path-from_root)))
            assert derror<1e-10
            max_error=max(max_error,derror)
            quotients.append({'turn':turn,'power':power,'chart_entry_n':n0,
                              'lifted_generator_turn':scaled[-1],
                              'principal_endpoint_turn':principal_turn(g**power) if turn==.4 else principal_turn(np.exp(2j*math.pi*turn*power)),
                              'scaled_log_error':err,'dyadic_path_error':derror})

    I=np.eye(2,dtype=complex)
    sigma=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    axis=np.array([2.,-1.,3.]);axis/=np.linalg.norm(axis)
    S=np.einsum('j,jab->ab',axis,sigma)
    def unitary(theta):return math.cos(theta)*I+1j*math.sin(theta)*S
    def adjoint(U):
        return np.array([[.5*np.trace(a@U@b@U.conj().T).real for b in sigma] for a in sigma])
    angle=.41
    paths=[]
    for n in (2,3,4,5,6):
        h=unitary(angle/(2**n))
        X=(2**n)*matrix_log(h)
        Y=(2**n)*matrix_log(adjoint(h))
        for t in (0.,.125,.375,.625,1.):
            left=adjoint(matrix_exp(t*X))
            right=matrix_exp(t*Y)
            err=float(np.linalg.norm(left-right))
            expected=float(np.linalg.norm(matrix_exp(t*X)-unitary(t*angle)))
            assert max(err,expected)<1e-10
            max_error=max(max_error,err,expected)
            paths.append({'n':n,'t':t,'adjoint_compatibility_error':err,
                          'source_path_error':expected})

    cuts=[]
    for j in (1,5,12):
        weight=2.**(-j)
        for eps in (1e-2,1e-4,1e-6):
            a,b=.5-eps,-.5+eps
            input_distance=min(abs(a-b),1-abs(a-b))*weight
            root_distance=min(abs(a/2-b/2),1-abs(a/2-b/2))*weight
            assert abs(input_distance-2*eps*weight)<1e-15
            assert abs(root_distance-(.5-eps)*weight)<1e-15
            assert abs(abs(a/2)*weight-abs(a)*weight/2)<1e-15
            cuts.append({'tail_index':j,'epsilon':eps,'input_cost_distance':input_distance,
                         'root_cost_distance':root_distance,'half_cost_exact':True})
    return {'round':1066,'scope':'Finite consistency checks only; LC, local contractibility and Lie conclusions are analytic.',
            'same_U_counterexample':large,'circle_quotient_paths':quotients,
            'su2_to_so3_compatibility':paths,'discontinuous_root_cut_witnesses':cuts,
            'max_matrix_or_phase_residual':max_error,
            'original_U_replay_claimed':False,'spatial_dimension_generated':False}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(result,json.loads(path.read_text(encoding='utf8')))
    print(json.dumps({'round':1066,'passed':True,'max_residual':result['max_matrix_or_phase_residual'],
                      'quotient_paths':len(result['circle_quotient_paths']),
                      'matrix_path_checks':len(result['su2_to_so3_compatibility'])},ensure_ascii=False))
