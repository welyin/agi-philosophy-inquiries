"""Round 1349. Direct physical tree-gauge propagation, without irrep sampling."""
import argparse
import json
from pathlib import Path
import numpy as np


def color(name):
    T=[]
    for i,j in [(0,1),(1,2),(2,0)]:
        e=np.zeros((3,3),complex);e[i,j]=1
        if name == "SU3":
            T.extend([(e+e.T)/2,1j*(e-e.T)/2])
        else:
            T.append(1j*(e-e.T))
    if name=="SU3":
        T.extend([np.diag([1.,-1.,0.])/2, np.diag([1.,1.,-2.])/np.sqrt(12)])
    C=np.trace(sum(x@x for x in T)).real/3
    Z=sum(np.linalg.matrix_power(np.kron(x,np.eye(3))+np.kron(np.eye(3),x),2) for x in T)/C
    for g in T:
        total=np.kron(g,np.eye(3))+np.kron(np.eye(3),g)
        assert np.max(abs(Z@total-total@Z))<1e-12
    return np.real_if_close(Z)


def evolution(Z,kappa,mu,t=.7,J=1., no_swap=False):
    exchange=np.eye(4).reshape(2,2,2,2).transpose(1,0,2,3).reshape(4,4)
    X=np.array([[0.,1.],[1.,0.]])
    # internal basis (chirality, flavour, bank, two active colours).
    I=np.eye(72)
    mass=t*np.kron(X,np.kron(mu*np.eye(4)+(0 if no_swap else exchange),np.eye(9)))
    y=np.kron(np.eye(8),kappa*(Z-np.eye(9)))+mass
    H=np.block([[0*I,J*I],[J*I,y]])
    eigen,basis=np.linalg.eigh(H)
    return eigen,basis


def observed(solver,T,bank):
    eigen,basis=solver
    initial=np.zeros((2,2,2,2,9,9),complex)
    # Both untouched conjugate colours are coherent spectator labels.
    initial[0,0,0,bank,:,:]=np.eye(9)/3
    final=(basis*np.exp(-1j*T*eigen))@(basis.conj().T@initial.reshape(144,9))
    p=np.sum(abs(final.reshape(2,2,2,2,9,9))**2,axis=(4,5))
    assert abs(p.sum()-1)<1e-12
    py=float(p[1].sum()); pflip=float(p[:,1].sum()); joint=float(p[1,1].sum())
    return dict(joint=joint, at_y=py, flipped=pflip,
                marginals_product=py*pflip, connected=joint-py*pflip,
                bank_changed=float(p[:,:,:,1-bank].sum()))


def exit_flip(delta,mass,T,J=1.):
    def transfer(u):
        freq=np.hypot(J,u/2)
        return -1j*J*np.exp(-1j*u*T/2)*np.sin(freq*T)/freq
    return abs(transfer(delta+mass)-transfer(delta-mass))**2/4


def run():
    groups={}
    for name in ["SU3","SO3"]:
        Z=color(name)
        eig=np.linalg.eigvalsh(Z)
        assert np.allclose(eig, [1]*3+[2.5]*6 if name=="SU3" else [0]+[1]*3+[3]*5)
        delta=1.3*(eig-1)
        rows={}
        for mu in [1.,2.]:
            full=evolution(Z,1.3,mu)
            flat=evolution(Z,0,mu)
            blind=evolution(Z,1.3,mu,no_swap=True)
            for T in [1.4,2.2]:
                tag=f"mu={mu:g},T={T:g}"
                data={}
                for b in [0,1]:
                    w=1-b/2
                    analytic=w*np.mean(exit_flip(delta,.7*(mu+1),T))
                    analytic+=(1-w)*np.mean(exit_flip(delta,.7*(mu-1),T))
                    output=observed(full,T,b)
                    assert np.isclose(analytic,output["joint"],atol=1e-12)
                    assert np.allclose(exit_flip(delta,.7*(mu+1),T),exit_flip(-delta,.7*(mu+1),T))
                    data["aligned" if b==0 else "orthogonal"]=output
                data["no_color_detuning_joint"]=[observed(flat,T,b)["joint"] for b in [0,1]]
                data["identity_mass_only_joint"]=[observed(blind,T,b)["joint"] for b in [0,1]]
                assert np.allclose(data["identity_mass_only_joint"][0],data["identity_mass_only_joint"][1])
                if mu==1:
                    assert np.isclose(data["orthogonal"]["joint"],.5*data["aligned"]["joint"])
                assert abs(data["orthogonal"]["connected"])>1e-4
                assert data["orthogonal"]["bank_changed"]>1e-5
                rows[tag]=data
        groups[name]=rows
    determinants={}
    for key in groups["SU3"]:
        matrix=np.array([[groups[name][key][b]["joint"] for b in ["aligned","orthogonal"]]
                         for name in ["SU3","SO3"]])
        value=float(np.linalg.det(matrix))
        if key.startswith("mu=1"):
            assert abs(value)<2e-14
        else:
            assert abs(value)>1e-7
        totals=float(matrix.sum())
        robust_error=2*abs(value)/(np.sqrt(totals**2+8*abs(value))+totals)
        same_group={}
        for name in ("SU3","SO3"):
            data=groups[name][key]
            flat=np.array(data["no_color_detuning_joint"])
            active=[data["aligned"]["joint"],data["orthogonal"]["joint"]]
            same_group[name]=float(np.linalg.det(np.array([active,flat])))
        determinants[key]=dict(value=value, entrywise_error_for_sign=robust_error,
                               within_group_active_vs_flat=same_group)
    return dict(full_basis_dimension=144,spectator_dimension=9,group=groups,
                color_by_bank_determinants=determinants)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    else:
        assert np.allclose(a,b,rtol=5e-11,atol=5e-12),(a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args(); result=run()
    path=Path(__file__).with_name("independent_results.json")
    if args.check:
        compare(result,json.loads(path.read_text(encoding="utf-8")))
    else:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2,allow_nan=False);out.write("\n")
    print("PASS",result["color_by_bank_determinants"])
