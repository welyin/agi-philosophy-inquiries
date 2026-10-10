"""Count-0 joint compatibility witness; finite checks do not prove global claims.
Default: deterministic JSON to stdout. --check: compare with saved results, read-only.
Requires only the project's existing Python and NumPy.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import numpy as np

I2=np.eye(2,dtype=complex)
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=complex)
ANISOTROPY=np.diag([1.0,0.75,0.5])
HERE=Path(__file__).resolve().parent

def tent(x,z):
    return max(0.0,1.0-float(np.linalg.norm(x-z)))

def task_distance(x,y):
    return min(1.0,float(np.linalg.norm(x-y)))

def clip_ball(v):
    return v/max(1.0,float(np.linalg.norm(v)))

def effect(x,o):
    return (I2+np.einsum('i,ijk->jk',clip_ball(x-o),PAULI))/2.0

def aligned_effect(x,o,O=None):
    # Primary noncovariant task; O is an explicit transported task parameter.
    if O is None: O=np.eye(3)
    v=ANISOTROPY@O@clip_ball(x-o)
    return (I2+np.einsum('i,ijk->jk',v,PAULI))/2.0

def psqrt(a):
    vals,vec=np.linalg.eigh((a+a.conj().T)/2.0)
    return (vec*np.sqrt(np.maximum(vals,0.0)))@vec.conj().T

def rodrigues(axis,angle):
    n=np.asarray(axis,dtype=float)
    n=n/np.linalg.norm(n)
    cross=np.array([[0,-n[2],n[1]],[n[2],0,-n[0]],[-n[1],n[0],0]])
    return math.cos(angle)*np.eye(3)+(1-math.cos(angle))*np.outer(n,n)+math.sin(angle)*cross

def spin_rotation(axis,angle):
    n=np.asarray(axis,dtype=float)
    n=n/np.linalg.norm(n)
    return math.cos(angle/2)*I2-1j*math.sin(angle/2)*np.einsum('i,ijk->jk',n,PAULI)

def rotation_plane(n,i,j,angle):
    out=np.eye(n)
    out[i,i]=out[j,j]=math.cos(angle)
    out[j,i]=math.sin(angle)
    out[i,j]=-math.sin(angle)
    return out

def ptrace_probe(qr):
    # Tensor order Q,R, each dimension 2.
    t=qr.reshape(2,2,2,2)
    return np.einsum('aras->rs',t)

def expectation(a,rho):
    return float(np.trace(a@rho).real)

def run():
    rng=np.random.default_rng(20261010)
    groups=[]
    all_residuals=[]
    checks=0
    def close(a,b,tol=4e-12):
        nonlocal checks
        residual=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
        all_residuals.append(residual)
        checks+=1
        assert residual<=tol,(residual,tol)
        return residual
    def le(a,b,tol=4e-12):
        nonlocal checks
        checks+=1
        assert a<=b+tol,(a,b,tol)
    def group(name,fn):
        before=checks
        start=len(all_residuals)
        detail=fn()
        groups.append({'name':name,'assertion_calls':checks-before,
            'max_equality_residual':max(all_residuals[start:],default=0.0),'detail':detail})
    points=rng.normal(size=(48,3))
    o=np.array([0.2,-0.3,0.1])
    r=0.5
    axis=np.array([1.0,2.0,-1.0])
    angle=0.73
    R=rodrigues(axis,angle)
    shift=np.array([0.17,-0.23,0.31])

    def metric_checks():
        for x,y in zip(points[::2],points[1::2]):
            d=task_distance(x,y)
            close(abs(tent(x,x)-tent(y,x)),d)
            for z in points[::8]:
                le(abs(tent(x,z)-tent(y,z)),d)
            for q_o in (o,np.zeros(3)):
                eig=np.linalg.eigvalsh(effect(x,q_o)-effect(y,q_o))
                le(max(abs(eig)),d)
            le(float(np.linalg.norm(clip_ball(x)-clip_ball(y))),float(np.linalg.norm(x-y)))
        for dist in (0.0,0.125,0.9,1.0,2.0):
            x=np.zeros(3); y=np.array([dist,0.0,0.0])
            close(abs(tent(x,x)-tent(y,x)),min(1.0,dist))
        return {'random_pairs':24,'tent_centers_per_pair':6,
                'distance_formula':'min(1, Euclidean distance)',
                'global_identity_status':'proved in model.md; finite samples only checked here',
                'declared_menu_only':True}
    group('declared_current_task_metric',metric_checks)

    def euclidean_checks():
        close(R.T@R,np.eye(3)); close(np.linalg.det(R),1.0)
        for x,y in zip(points[::2],points[1::2]):
            gx=R@x+shift; gy=R@y+shift
            close(task_distance(gx,gy),task_distance(x,y))
            close(tent(gx,R@o+shift),tent(x,o))
            close(R.T@(gx-shift),x)
        A=rodrigues([0,1,0],0.41)
        t=np.array([-0.2,0.1,0.6])
        for x in points[:6]:
            close(A@(R@x+shift)+t,(A@R)@x+(A@shift+t))
        return {'actual_action':'x -> R x + a; probe and passive reference unchanged',
                'stabilizer_of_anchor':'SO(3)',
                'implemented_here_as':'adopted classical pullback *-automorphisms, not derived autonomous control'}
    group('actual_euclidean_moves_and_test_transport',euclidean_checks)

    def meeting_checks():
        for x,y in zip(points[::2],points[1::2]):
            m=(x+y)/2; v=(x-y)/2
            close(m+v,x); close(m-v,y)
            close((y+x)/2,m); close((x+x)/2,x)
            close(((R@x+shift)+(R@y+shift))/2,R@m+shift)
            reflected=2*o-x
            close((x+reflected)/2,o)
            close(2*o-reflected,x)
            # The relative record changes sign under role exchange.
            close((y-x)/2,-v)
        return {'meeting':'m=(x+y)/2', 'kept_relative_record':'v=(x-y)/2',
                'full_reversible_map':'(x,y) <-> (m,v)',
                'actual_compensation':'(x,o) -> (2o-x,o)',
                'unknown_qubit_copied':False}
    group('meeting_with_retained_relative_record',meeting_checks)

    def coordination_checks():
        # Compact ball ||x-o||<=1; one actual positive-scaling step, not infinite execution.
        for direction in points[:16]:
            x=o+0.8*direction/np.linalg.norm(direction)
            rx=o+(x-o)/2
            close(np.linalg.norm(rx-o),np.linalg.norm(x-o)/2)
            le(np.linalg.norm(rx-o),np.linalg.norm(x-o)-0.1)
            for t in (0.0,0.25,0.75,1.0):
                b=o+(1-t/2)*(x-o)
                close(o+(b-o)/(1-t/2),x)
                le(np.linalg.norm(b-o),1.0)
            close(o+(1-0/2)*(x-o),x)
            close(o+(1-1/2)*(x-o),rx)
        close(o+(1-0.7/2)*(o-o),o)
        return {'compact_neighborhood_radius':1.0,'C':'||x-o||',
                'R':'o+(x-o)/2','B':'o+(1-t/2)(x-o)',
                'finite_cost_infinite_completion_claimed':False,
                'C_is_execution_energy':False}
    group('one_step_continuous_coordination',coordination_checks)

    def qubit_checks():
        U=spin_rotation(axis,angle)
        close(U.conj().T@U,I2)
        min_eig=1.0
        for x in points:
            E=effect(x,o)
            vals=np.linalg.eigvalsh(E)
            min_eig=min(min_eig,float(vals[0]))
            le(-float(vals[0]),0.0); le(float(vals[-1]),1.0)
            close(np.trace(E),1.0)
            kp=psqrt(E); km=psqrt(I2-E)
            close(kp.conj().T@kp+km.conj().T@km,I2)
            close(effect(o+R@(x-o),o),U@E@U.conj().T)
        gaps=[]
        for direction in points[:16]:
            u=direction/np.linalg.norm(direction)
            E=effect(o+r*u,o); F=effect(o-r*u,o)
            diff=E-F
            close(diff,r*np.einsum('i,ijk->jk',u,PAULI))
            gap=float(np.max(np.abs(np.linalg.eigvalsh(diff))))
            close(gap,r); gaps.append(gap)
        span=np.column_stack([np.array([np.trace(effect(o+r*e,o)).real]+[
            np.trace(effect(o+r*e,o)@s).real for s in PAULI])
            for e in (np.zeros(3),*np.eye(3))])
        close(np.linalg.det(span),r**3)
        return {'shell_radius':r,'antipodal_effect_gap':r,
                'effect_span_determinant':float(np.linalg.det(span)),
                'minimum_sample_effect_eigenvalue':min_eig,
                'repeatability_1073_required_or_claimed':False,
                'covariance':'isotropic control branch only; primary anisotropic branch deliberately fails it'}
    group('isotropic_qubit_control',qubit_checks)

    def anisotropic_checks():
        for x,y in zip(points[::2],points[1::2]):
            E=aligned_effect(x,o)
            vals=np.linalg.eigvalsh(E)
            le(-float(vals[0]),0.0); le(float(vals[-1]),1.0)
            close(np.trace(E),1.0)
            kp=psqrt(E); km=psqrt(I2-E)
            close(kp.conj().T@kp+km.conj().T@km,I2)
            le(max(abs(np.linalg.eigvalsh(E-aligned_effect(y,o)))),task_distance(x,y))
            # Joint test transport is not SU(2) conjugation of one fixed task.
            close(aligned_effect(R@x+shift,R@o+shift,R.T),E)
        for direction in points[:16]:
            u=direction/np.linalg.norm(direction)
            plus=aligned_effect(o+r*u,o); minus=aligned_effect(o-r*u,o)
            close(plus-minus,r*np.einsum('i,ijk->jk',ANISOTROPY@u,PAULI))
            gap=float(max(abs(np.linalg.eigvalsh(plus-minus))))
            le(r/2,gap); le(gap,r)
        Ey=rodrigues([0,1,0],math.pi/2)
        x=o+r*np.array([1.,0.,0.])
        before=aligned_effect(x,o); after=aligned_effect(o+Ey@(x-o),o)
        purity_before=float(np.trace(before@before).real)
        purity_after=float(np.trace(after@after).real)
        close(purity_before,5/8); close(purity_after,17/32)
        close(purity_before-purity_after,3/32)
        span=np.column_stack([np.array([np.trace(aligned_effect(o+r*e,o)).real]+[
            np.trace(aligned_effect(o+r*e,o)@s).real for s in PAULI])
            for e in (np.zeros(3),*np.eye(3))])
        close(np.linalg.det(span),r**3*3/8)
        return {'primary_effect':'(I + (diag(1,3/4,1/2) O s(x-o)) dot sigma)/2',
                'primary_fixed_task_O':'identity', 'shell_radius':r,
                'opposite_gap_uniform_lower':r/2, 'effect_span_determinant':float(np.linalg.det(span)),
                'same_fixed_task_rotation_purities':[purity_before,purity_after],
                'unitary_covariance_failure_exact_gap':'3/32',
                'moved_task_parameter':'O -> O R^-1; o -> R o + a',
                'repeatability_required_or_claimed':False,
                'actual_rotations_and_upper_interface_compatible_without_unitary_covariance':True}
    group('primary_anisotropic_interface_without_qubit_covariance',anisotropic_checks)

    def reference_checks():
        bell=np.array([1,0,0,1],dtype=complex)/math.sqrt(2)
        vec=rng.normal(size=4)+1j*rng.normal(size=4)
        vec=vec/np.linalg.norm(vec)
        states=[np.outer(bell,bell.conj()),np.outer(vec,vec.conj()),
                0.37*np.outer(vec,vec.conj())+0.63*np.eye(4)/4]
        for state in states:
            for x in (o+r*np.array([1.,0,0]),o+r*np.array([0,0,1.])):
                E=aligned_effect(x,o)
                kp=np.kron(psqrt(E),I2); km=np.kron(psqrt(I2-E),I2)
                plus=kp@state@kp.conj().T; minus=km@state@km.conj().T
                close(np.trace(plus+minus),1.0)
                close(ptrace_probe(plus+minus),ptrace_probe(state))
                close(np.trace(plus),np.trace(np.kron(E,I2)@state))
                close(np.trace(minus),np.trace(np.kron(I2-E,I2)@state))
                for branch in (plus,minus):
                    le(-float(np.min(np.linalg.eigvalsh(branch))),0.0)
                # An explicit coherent report isometry includes measurement back-action.
                V=np.vstack((kp,km))
                close(V.conj().T@V,np.eye(4))
                out=V@state@V.conj().T
                close(out[:4,:4],plus); close(out[4:,4:],minus)
        return {'reference_dimension':2,'input_states':3,'endpoints':2,
                'preserved':'unconditional passive reference marginal',
                'not_claimed':'probe state, individual conditional reference, or all later natural dynamics unchanged',
                'report_and_environment_retained':True}
    group('instrument_with_passive_quantum_reference',reference_checks)

    def noncommuting_checks():
        A=rodrigues([0,0,1],math.pi/2)
        B=rodrigues([1,0,0],math.pi/2)
        initial=r*np.array([0.,0.,1.])
        ab=o+A@B@initial; ba=o+B@A@initial
        close(ab,o+r*np.array([1.,0.,0.]))
        close(ba,o-r*np.array([0.,1.,0.]))
        plus_x=(I2+PAULI[0])/2
        p_ab=expectation(aligned_effect(ab,o),plus_x)
        p_ba=expectation(aligned_effect(ba,o),plus_x)
        close(p_ab,(1+r)/2); close(p_ba,0.5); close(p_ab-p_ba,r/2)
        close(task_distance(ab,ba),min(1.0,r*math.sqrt(2)))
        close(abs(tent(ab,ab)-tent(ba,ab)),task_distance(ab,ba))
        for t in (0.0,0.2,0.7,1.0):
            At=rodrigues([0,0,1],t*math.pi/2)
            Bt=rodrigues([1,0,0],t*math.pi/2)
            close(o+At@(o-o),o); close(o+Bt@(o-o),o)
        return {'actual_continuous_anchored_paths':['R_z(t*pi/2)','R_x(t*pi/2)'],
                'fixed_later_probe':'+X','same_reading_rule':True,
                'p_AB':p_ab,'p_BA':p_ba,'probability_gap':p_ab-p_ba,
                'tent_task_gap':task_distance(ab,ba),
                'geometric_noncommutation_already_visible_without_qubit_covariance':True}
    group('actual_anchored_order_difference',noncommuting_checks)

    def deletion_controls():
        a=rotation_plane(2,0,1,0.73); b=rotation_plane(2,0,1,-0.29)
        close(a@b,b@a)
        e2=np.array([0.3,0.4])
        eplus=(I2+e2[0]*PAULI[0]+e2[1]*PAULI[1])/2
        eminus=(I2-e2[0]*PAULI[0]-e2[1]*PAULI[1])/2
        close(np.max(abs(np.linalg.eigvalsh(eplus-eminus))),0.5)
        A4=rotation_plane(4,0,1,math.pi/2); B4=rotation_plane(4,1,2,math.pi/2)
        start=r*np.eye(4)[0]
        ab=A4@B4@start; ba=B4@A4@start
        le(0.5,float(np.linalg.norm(ab-ba)))
        hplus=r*np.eye(4)[3]; hminus=-hplus
        projected=lambda v:(I2+np.einsum('i,ijk->jk',clip_ball(v)[:3],PAULI))/2
        close(projected(hplus),I2/2); close(projected(hminus),I2/2)
        close(task_distance(hplus,hminus),1.0)
        return {'dimension_2':{'keeps':'tent tasks, midpoint, coordination, opposite qubit separation',
                                'deleted':'continuous anchored noncommuting Euclidean motions',
                                'SO2_commutation_residual':float(np.max(abs(a@b-b@a)))},
                'dimension_4':{'keeps':'classical actual tasks, midpoint, coordination, anchored noncommutation',
                                'deleted':'one complete continuous qubit task separating every opposite pair',
                                'opposite_e4_qubit_gap':0.0,'opposite_e4_tent_gap':1.0},
                'status':'explicit deletion examples only; global exclusion is not proved by these samples'}
    group('dimension_2_and_4_deleted_contract_controls',deletion_controls)
    return {'classification':'count-0 joint conditional compatibility witness',
            'selected_model':'C_0(R^3,M_2); primary anisotropic bounded effects in multiplier algebra',
            'dimension_generation_proved':False,
            'runtime':{'python':sys.version.split()[0],'numpy':np.__version__},
            'seed':20261010,'check_groups':len(groups),'assertion_calls':checks,
            'max_equality_residual':max(all_residuals,default=0.0),'groups':groups,
            'resource_ledger':{
                'source_registers':['actual endpoint x','optional second endpoint y','anchor o','unknown qubit Q','passive reference R'],
                'control_inputs':['test center z','finite Euclidean move parameters','finite coordination parameters','adopted meeting/compensation permission'],
                'records_retained':['meeting relative vector v','measurement report','measurement environment as required'],
                'role_fair_accounting':'same register and operation counts after swapping x,y; no erasure of their relative record',
                'precision':'exact continuum operations are ideal inputs; finite precision budgets are not derived',
                'execution_energy_or_autonomous_law_derived':False,
                'infinite_coordination_at_finite_resource_claimed':False},
            'scope':['global identities are analytic statements in model.md',
                     'finite NumPy checks validate formulas only',
                     'actual instrument/control permissions and the chosen R^3 base are inputs',
                     'declared single-call operational metric excludes arbitrary bounded effects and unrestricted repetition outside the task menu',
                     'no physical minimum scale or microscopic continuity is established',
                     'not a sole-qubit representation of all position states or translations'],
            'sha256_check':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}

if __name__=='__main__':
    result=run()
    if sys.argv[1:]==['--check']:
        saved=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        if saved!=result:
            raise SystemExit('saved results differ')
        print(json.dumps({'status':'pass','check_groups':result['check_groups'],
                          'assertion_calls':result['assertion_calls'],
                          'max_equality_residual':result['max_equality_residual']},ensure_ascii=False,sort_keys=True))
    elif sys.argv[1:]:
        raise SystemExit('usage: check.py [--check]')
    else:
        print(json.dumps(result,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False))
