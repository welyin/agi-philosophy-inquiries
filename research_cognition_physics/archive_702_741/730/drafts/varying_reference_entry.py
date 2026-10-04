"""Non-numbered730 entry: inherited constant-mass diagonalization versus original fields."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_matter_ground_source as old
TARGET=HERE/'varying_reference_entry_results.json'

def run():
    N=8;q=old.matter.original.shared_source(N)
    hs=float(q['phi'][0,0,0,1]);ss=float(q['phi'][0,0,0,4])/1.06
    yn=old.matter.Y['nu'];ys=old.matter.Y['s'];a=abs(yn);b=abs(ys)
    pr=-np.angle(ys)/2;pl=-np.angle(yn)-pr
    phases=np.diag(np.exp(1j*np.array([pl,pr])))
    def fields(x):
        xx,yy,zz=x
        h=hs*(1+.05*np.sin(xx+yy));s=ss*(1+.06*np.cos(yy))
        dh=np.array([.05*hs*np.cos(xx+yy)]*2+[0.])
        ds=np.array([0.,-.06*ss*np.sin(yy),0.])
        return h,s,dh,ds
    def frame(x):
        h,s,dh,ds=fields(x);phi=np.array([0.,h,0.,0.,s])
        F=float(old.matter.original.F(phi));D=a*h/np.sqrt(F);MR=b*s/np.sqrt(F)
        theta=.5*np.arctan2(2*D,MR)
        grad=a*b*(s*dh-h*ds)/(b*b*s*s+4*a*a*h*h)
        O=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
        U=phases@O@np.diag([1j,1.])
        mass=np.array([[0.,yn*h],[yn*h,ys*s]])/np.sqrt(F)
        mm=(np.sqrt(MR*MR+4*D*D)-MR)/2;mp=mm+MR
        return theta,grad,U,mass,np.array([mm,mp])
    angles=[];errs=[];grams=[];maxgrad=0.
    for idx in np.ndindex((N,N,N)):
        x=2*np.pi*np.array(idx)/N;h,s,_,_=fields(x)
        errs.append(float(np.max(abs(q['phi'][idx]-np.array([0.,h,0.,0.,s])))))
        th,grad,U,M,mm=frame(x);angles.append(th);maxgrad=max(maxgrad,float(np.linalg.norm(grad)))
        errs.append(float(np.max(abs(U.T@M@U-np.diag(mm)))))
        pair,_=old.radial(h,s,old.IDS)
        eig=np.linalg.eigvalsh(old.bdg(*pair))
        errs.append(float(np.max(abs(np.sort(abs(eig))-np.sort(np.repeat(mm,4))))))
        grams.append(M.conj().T@M)
    i=int(np.argmin(angles));j=int(np.argmax(angles))
    comm=float(np.linalg.norm(grams[i]@grams[j]-grams[j]@grams[i],'fro'))
    x=np.array([.23,.57,.19]);theta,grad,U,M,mm=frame(x)
    step=1e-5;connection_errors=[]
    for axis in range(3):
        dx=np.eye(3)[axis]*step
        dU=(frame(x+dx)[2]-frame(x-dx)[2])/(2*step)
        connection=U.conj().T@dU
        expected=-1j*grad[axis]*np.array([[0,1],[1,0]])
        connection_errors.append(float(np.max(abs(connection-expected))))
    assert max(errs)<2e-13 and max(connection_errors)<2e-10 and comm>1e-7
    deps=('research_note_630.md','research_note_633.md','research_note_663.md','research_note_667.md',
          'research_note_729.md','joint_fermion_gauss_completion.py','joint_curved_quantum_source.py')
    return dict(entry_round=730,new_formal_round=False,original_grid=N,all_sampled_points=N**3,
        theta_range=[float(min(angles)),float(max(angles))],max_sampled_gradient=maxgrad,
        fixed_original_pair_indices=[i,j],squared_mass_commutator_Frobenius=comm,
        mass_and_original_background_error=max(errs),frame_connection_difference_error=max(connection_errors),
        nonconstant_frame_required=True,continuous_reference_state_not_yet_constructed=True,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        scope='Necessary mass-sector audit only. Constant rephasing retains original complex Y. This does not decouple the original electroweak theory; transformed gauge, kinetic, source and actual sterile record must all be retained. No new formal round, Hadamard construction, autonomous preparation or Einstein dynamics claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

