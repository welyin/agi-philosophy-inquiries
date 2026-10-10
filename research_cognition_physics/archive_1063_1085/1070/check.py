"""Round1070: NumPy checks of analytic finite-window bound. Read-only."""
import json
import numpy as np
I=np.eye(2,dtype=complex)
S=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=complex)
checks=0
residuals=[]
def near(a,b):
    global checks
    v=float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    residuals.append(v); checks+=1
    assert v<2e-12,v
def rho(x):return (I+np.einsum('i,ijk->jk',x,S))/2
def bloch(m):return np.array([np.trace(m@s).real for s in S])
def channel(ks,m):return sum((k@m@k.conj().T for k in ks),np.zeros((2,2),complex))
def dist(a,b):return float(np.sum(np.abs(np.linalg.eigvalsh(a-b)))/2)
def cert(ks):
    global checks
    near(sum((k.conj().T@k for k in ks),np.zeros((2,2),complex)),I)
    vs=[k.reshape(-1,order='F') for k in ks]
    choi=sum((np.outer(v,v.conj()) for v in vs),np.zeros((4,4),complex))
    checks+=1; assert np.linalg.eigvalsh(choi)[0]>-2e-12
z=np.array([0.,0.,1.]); p=np.array([1.,0.],complex)
directions=[z,-z,np.array([1.,0,0]),np.array([0.,1,0]),np.array([2.,2.,1.])/3]
rows=[]
for r,a in [(0.5,0.25),(0.4,1/3),(0.75,0.25),(0.1,0.875),(0,1),(2/3,0)]:
    reset=[np.sqrt(1-a)*I]+[np.sqrt(a)*np.outer(p,I[:,j].conj()) for j in range(2)]
    ad=[np.diag([1.,np.sqrt(1-a)]),np.array([[0.,np.sqrt(a)],[0.,0.]])]
    cert(reset);cert(ad)
    errors=[[],[]]
    for n in directions:
        for f in [0.,0.5,1.]:
            x=r*f*n; target=rho(x+a*z)
            checks+=1; assert np.linalg.eigvalsh(target)[0]>-2e-12
            out=channel(reset,rho(x)); val=dist(out,target)
            near(bloch(out),(1-a)*x+a*z);near(val,a*np.linalg.norm(x)/2)
            errors[0].append(val)
            out=channel(ad,rho(x)); val=dist(out,target)
            k=1-np.sqrt(1-a)
            near(bloch(out),np.array([np.sqrt(1-a)*x[0],np.sqrt(1-a)*x[1],(1-a)*x[2]+a]))
            near(val,np.sqrt(k*k*(x[0]**2+x[1]**2)+a*a*x[2]**2)/2)
            checks+=1;assert val<=r*a/2+2e-12
            errors[1].append(val)
    for e in errors:near(max(e),r*a/2)
    rows.append({'r':r,'a':a,'analytic_optimum':r*a/2,'reset_max':max(errors[0]),'amplitude_damping_max':max(errors[1])})
# Four affine-independent legal mixed states span Hermitian matrices.
points=[np.zeros(3)]+[0.2*np.eye(3)[j] for j in range(3)]
B=np.array([[1.,*x] for x in points]);near(np.linalg.det(B),0.008)
U=(I-1j*S[1])/np.sqrt(2)
for x in points:near(U.conj().T@(U@rho(x)@U.conj().T)@U,rho(x))
# The ideal nonzero translation is not positive on all input states.
near(np.linalg.eigvalsh(rho((1+0.25)*z))[0],-0.125)
print(json.dumps({'round':1070,'passed':True,'checks':checks,'max_residual':max(residuals),
 'cases':rows,'fixed_example_exact':'1/16','diamond_bound_claimed':False,
 'actual_space_derived':False,'samples_replace_proof':False},ensure_ascii=False,indent=2,sort_keys=True))
