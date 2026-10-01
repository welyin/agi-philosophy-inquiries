"""Exploration only: original subgroup, one-site flat overlap, auxiliary pairing."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import joint_spinor_subgroup_mass as old
import joint_chiral_fibre_source as spin

_,ann=old.clifford()
gs=[]
for a in ann:gs.extend([a+a.conj().T,-1j*(a-a.conj().T)])
C=np.eye(32,dtype=complex)
for j in (1,3,5,7,9):C=C@gs[j]
T=[(C@g)[np.ix_(old.MASKS,old.MASKS)] for g in gs]
CD=spin.GAMMA[1]@spin.GAMMA[3]
B=1j*spin.G5@CD
e,V=np.linalg.eigh(spin.G5);vp=V[:,e>.5]
charges=np.array([sum([-2,-2,-2,3,3][j] for j in s) for s in old.STATES])
print('symmetric T',max(np.max(abs(t-t.T)) for t in T),'skew B',np.max(abs(B+B.T)))
print('spin basis pair',vp.T@B@vp)
for theta in (.17,.6,np.pi/3):
    us=[(np.cos(q*theta/2)*np.eye(4)+1j*np.sin(q*theta/2)*spin.GAMMA[3])@vp for q in charges]
    U=np.zeros((64,32),complex)
    for a,u in enumerate(us):U[4*a:4*a+4,2*a:2*a+2]=u
    E=np.arange(1,11,dtype=float);E/=np.linalg.norm(E)
    A=U.T@np.kron(sum(e*t for e,t in zip(E,T)),B)@U
    expected=np.kron(sum(e*t*(np.cos(theta) if a<6 else np.cos(1.5*theta)) for a,(e,t) in enumerate(zip(E,T))),vp.T@B@vp)
    print(theta,'pair residual',np.max(abs(A-expected)),'sv min',np.linalg.svd(A,compute_uv=False)[-1],
          'det magnitude',abs(np.linalg.det(A)))
