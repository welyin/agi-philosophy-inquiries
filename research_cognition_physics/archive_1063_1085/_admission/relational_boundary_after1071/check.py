"""Source audit only: fixed consensus, vector and quadrupole boundary tasks.
Read-only stdout JSON; no installed packages other than existing NumPy.
"""
import json
import numpy as np
I=np.eye(2,dtype=complex)
pauli=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
x,y,z=pauli
k=np.kron
S=sum(k(t,t) for t in pauli)
w=1/4
omega=(np.eye(4)-w*S)/4
epsilon=3/32
checks=0
maximum=0.0

def eq(a,b,tol=2e-12):
    global checks,maximum
    e=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    maximum=max(maximum,e)
    assert e<tol,e
    checks+=1

def positive(a):
    global checks
    assert np.linalg.eigvalsh(a).min()>-2e-12
    eq(a,a.conj().T)
    eq(np.trace(a),1)
    checks+=1

def marginal(rho,keep):
    tensor=rho.reshape((2,)*6)
    axes=list(range(3))
    for a in reversed(range(3)):
        if a not in keep:
            pos=axes.index(a)
            tensor=np.trace(tensor,axis1=pos,axis2=pos+len(axes))
            axes.remove(a)
    return tensor.reshape(2**len(keep),2**len(keep))

def sigma(n):return np.einsum('i,ijk->jk',n,pauli)
def Q(n):return k(sigma(n),sigma(n))-S/3

def rho(n):return k(omega,I/2)+epsilon/2*k(Q(n),z)
def unitary(axis,theta):return np.cos(theta/2)*I-1j*np.sin(theta/2)*sigma(axis)
def rotated(U,n):return np.array([np.trace(t@U@sigma(n)@U.conj().T).real/2 for t in pauli])
def action(U,r):
    V=k(k(U,U),I)
    return V@r@V.conj().T

def probability(r,O):return float(np.trace(r@(np.eye(8)+O)/2).real)

axes=np.eye(3)
directions=[*axes, -axes[0], -axes[1], -axes[2],np.ones(3)/np.sqrt(3),np.array([2.,-1.,3.])/np.sqrt(14)]
rotations=[unitary(a,t) for a in [*axes,np.ones(3)/np.sqrt(3)] for t in [.23,.9,np.pi/2]]
for n in directions:
    positive(rho(n))
    eq(np.linalg.eigvalsh(Q(n)),[-4/3,0,2/3,2/3])
    eq(marginal(rho(n),[0,1]),omega)
    eq(marginal(rho(n),[0,2]),np.eye(4)/4)
    eq(marginal(rho(n),[1,2]),np.eye(4)/4)
    eq(rho(n),rho(-n))
    for m in directions:
        O=k(k(sigma(m),sigma(m)),z)
        eq(probability(rho(n),O),.5+2*epsilon*(float(m@n)**2-1/3))
    for U in rotations:
        eq(action(U,rho(n)),rho(rotated(U,n)))
        eq(k(U,U)@omega@k(U,U).conj().T,omega)

O_z=k(k(z,z),z)
p_z=probability(rho(axes[2]),O_z)
p_x=probability(rho(axes[0]),O_z)
eq(p_z,5/8);eq(p_x,7/16)
a=unitary(axes[2],np.pi/2);b=unitary(axes[0],np.pi/2)
r_ab=action(a@b,rho(axes[2]));r_ba=action(b@a,rho(axes[2]))
O_x=k(k(x,x),z)
p_ab=probability(r_ab,O_x);p_ba=probability(r_ba,O_x)
eq(p_ab,5/8);eq(p_ba,7/16)
eq(action(unitary(axes[0],np.pi),rho(axes[2])),rho(axes[2]))
eq(action(unitary(axes[2],.51),rho(axes[2])),rho(axes[2]))
# Fixed Bloch vectors of both stabilizer actions have zero intersection.
Rz=np.column_stack([rotated(unitary(axes[2],.51),n) for n in axes])
Rx=np.column_stack([rotated(unitary(axes[0],np.pi),n) for n in axes])
constraint=np.vstack([Rz-np.eye(3),Rx-np.eye(3)])
eq(np.linalg.matrix_rank(constraint,tol=1e-10),3)
# The contrasting vector extension: covariance on A with C retained.
t=1/32
def vector_state(n):return k(omega,I/2)+t/2*k(k(sigma(n),I),z)
for n in directions:
    positive(vector_state(n))
    eq(marginal(vector_state(n),[0,1]),omega)
    eq(marginal(vector_state(n),[1,2]),np.eye(4)/4)
    eq(marginal(vector_state(n),[0,2]),np.eye(4)/4+t*k(sigma(n),z))
    for U in rotations:
        eq(action(U,vector_state(n)),vector_state(rotated(U,n)))
    for m in axes:
        eq(probability(vector_state(n),k(k(sigma(m),I),z)),.5+2*t*float(m@n))
# Product extension has no directional dependence even after arbitrary common U.
product=k(omega,I/2)
for U in rotations:eq(action(U,product),product)
print(json.dumps({'audit':'relational_boundary_after1071','new_scientific_rounds':0,
 'checks':checks,'maximum_residual':maximum,'w':w,'quadrupole_epsilon':epsilon,
 'minimum_joint_eigenvalue':float(np.linalg.eigvalsh(rho(axes[2])).min()),
 'axis_probabilities':[p_z,p_x], 'order_probabilities':[p_ab,p_ba],
 'probability_gap':p_z-p_x,'exact_covariant_qubit_uniform_error_lower_bound':epsilon,
 'vector_amplitude':t,'spatial_dimension_derived':False},ensure_ascii=False,indent=2))
