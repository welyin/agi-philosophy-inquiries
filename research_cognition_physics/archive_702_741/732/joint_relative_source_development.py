"""732: same-record smooth relative source and conditional linear development.

Numerical checks concern original730's complete finite local symbol, its source
identities and rank-four mode transport. The continuum coupled linear PDE
argument is in note732, not claimed as a numerical spacetime solution.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
import joint_dynamic_continuum_reference as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_relative_source_development_results.json'
ORDER=np.r_[np.arange(32),np.arange(64,96),np.arange(96,128),np.arange(32,64)]


def assemble(items):
    return old.old.block(items)[np.ix_(ORDER,ORDER)]


def maxabs(x):
    return float(np.max(np.abs(x)))


def frame_jet(g,dg):
    ev,V=np.linalg.eigh(g);assert ev.min()>0
    root=np.sqrt(ev);F=(V/root)@V.T
    divided=-1/(root[:,None]*root[None,:]*(root[:,None]+root[None,:]))
    dF=V@(divided*(V.T@dg@V))@V.T
    return F,dF


def mass_jet(phi,dphi):
    M=old.old.bdg(*old.old.matter.mass_matrices(phi))
    F=old.old.matter.original.F(phi)
    unit=np.eye(5)
    numerator=sum(dphi[a]*np.sqrt(old.old.matter.original.F(unit[a]))*
                  old.old.bdg(*old.old.matter.mass_matrices(unit[a])) for a in range(5))
    dM=numerator/np.sqrt(F)+M*(phi@dphi)/(6*F)
    return M,dM


def matrix_and_source(g,phi,a,a0,dg,dphi,da,da0):
    e,de=frame_jet(g,dg);M,dM=mass_jet(phi,dphi)
    gen=old.gauge_h(a,a0);dgen=old.gauge_h(da,da0)
    momentum=np.array([.31,-.27,.19]);Bs=[];Ds=[]
    for sign in (1,-1):
        k=sign*momentum
        K=sum(old.GAMMA[j]*sum(e[j,i]*k[i] for i in range(3)) for j in range(3))
        dK=sum(old.GAMMA[j]*sum(de[j,i]*k[i] for i in range(3)) for j in range(3))
        gauge=sum(old.GAMMA[j][:32,:32]@sum(e[j,i]*gen[i] for i in range(3)) for j in range(3))
        dgauge=sum(old.GAMMA[j][:32,:32]@sum(de[j,i]*gen[i]+e[j,i]*dgen[i]
                   for i in range(3)) for j in range(3))
        Bs.append(K+old.old.bdg(gauge,np.zeros_like(gauge))+M)
        Ds.append(dK+old.old.bdg(dgauge,np.zeros_like(dgauge))+dM)
    return assemble(Bs),assemble(Ds)


def background(t):
    p=old.point()
    return (p['g']+t*p['dg'],p['phi']+t*p['v'],p['a']+t*p['da'],p['a0']+t*p['da0'])


def Btime(t):
    bg=background(t)
    return matrix_and_source(*bg,*[np.zeros_like(x) for x in bg])[0]


def unitary(t,dt):
    B=Btime(t+dt/2);ev,V=np.linalg.eigh(B)
    return (V*np.exp(-1j*dt*ev))@V.conj().T


@lru_cache(maxsize=1)
def initial_modes():
    P,_,B,_=old.flow(.08,96)
    e=np.zeros(128,complex);e[30]=e[62]=1/np.sqrt(2)
    C=np.block([[np.zeros((64,64)),np.eye(64)],[np.eye(64),np.zeros((64,64))]])
    ce=C@e.conj();Q=np.outer(e,e.conj())+np.outer(ce,ce.conj());R=np.eye(128)-2*Q
    delta=(R@P@R-P)/2
    ev,V=np.linalg.eigh(delta);keep=abs(ev)>1e-10
    assert keep.sum()<=4
    values=ev[keep];modes=V[:,keep]
    assert maxabs((modes*values)@modes.conj().T-delta)<2e-13
    assert maxabs(B-Btime(0))<2e-13
    return P,delta,values,modes,C


@lru_cache(maxsize=1)
def transported():
    P,delta,values,modes,C=initial_modes()
    post=P+delta;time=0.;dt=.04/64
    modes=modes.copy();P=P.copy();post=post.copy()
    for _ in range(64):
        U=unitary(time,dt);time+=dt
        P=U@P@U.conj().T;post=U@post@U.conj().T;modes=U@modes
    change=(modes*values)@modes.conj().T
    return time,P,post,values,modes,change,C


def physical_generators():
    matter=old.old.matter;out=[]
    for color in matter.gauge.generators(3):
        H=np.zeros((32,32),complex)
        for name,spinweak in [('Q',4),('u',2),('d',2)]:
            sl=matter.SLICES[name];H[sl,sl]=np.kron(color,np.eye(spinweak))
        out.append(('color',None,H))
    for index,weak in enumerate(np.asarray(old.SIG)/2):
        H=np.zeros((32,32),complex)
        for name,colors in [('Q',3),('L',1)]:
            sl=matter.SLICES[name];H[sl,sl]=np.kron(np.eye(colors),np.kron(weak,np.eye(2)))
        out.append(('weak',index,H))
    H=np.zeros((32,32),complex)
    for name,sl in matter.SLICES.items():H[sl,sl]=old.old.chiral.CHARGES[name]*np.eye(sl.stop-sl.start)
    out.append(('circle',None,H))
    return out


def source(delta,D):
    return float(np.trace(D@delta).real/2)


def relative_mode_transport_check():
    time,P,post,values,modes,delta,C=transported();B=Btime(time)
    direct=post-P
    source_error=abs(source(direct,B)-float(sum(v*np.vdot(u,B@u).real/2 for v,u in zip(values,modes.T))))
    errors=[maxabs(delta-direct),maxabs(C@delta.conj()@C+delta),source_error,
            maxabs(modes.conj().T@modes-np.eye(len(values)))]
    assert max(errors)<3e-12
    assert np.linalg.eigvalsh(post).min()>-3e-12 and np.linalg.eigvalsh(post).max()<1+3e-12
    assert abs(source(delta,B))>.01
    return dict(physical_modes=64,Nambu_dimension=128,record_difference_rank=len(values),
        signed_mode_weights=values.tolist(),time=time,mode_vs_full_covariance_error=max(errors),
        actual_record_relative_energy=source(delta,B),
        finite_smooth_mode_reduction_not_a_continuum_discretization=True)


def gauge_and_mass_source_check():
    time,_,_,_,_,delta,_=transported();bg=background(time);g,phi,a,a0=bg
    B=Btime(time);dotdelta=-1j*(B@delta-delta@B)
    rows=[];operator_error=0.;balance_error=0.
    for kind,index,H in physical_generators():
        charge=assemble([old.old.bdg(H,np.zeros_like(H))]*2)
        dphi=np.zeros(5);da=np.zeros_like(a)
        X=phi[:2]+1j*phi[2:4]
        if kind=='weak':
            direction=np.eye(3)[index];change=1j*(old.SIG[index]/2)@X
            dphi[:4]=np.r_[change.real,change.imag];da=np.cross(a,direction)
        elif kind=='circle':
            change=3j*X;dphi[:4]=np.r_[change.real,change.imag]
        zero=[np.zeros_like(x) for x in bg]
        _,massD=matrix_and_source(*bg,zero[0],dphi,zero[2],zero[3])
        _,gaugeD=matrix_and_source(*bg,zero[0],zero[1],da,zero[3])
        D=massD+gaugeD
        operator_error=max(operator_error,maxabs(D-1j*(charge@B-B@charge)))
        rate=source(dotdelta,charge);force=source(delta,D)
        balance_error=max(balance_error,abs(rate+force))
        rows.append(dict(kind=kind,index=index,charge_rate=rate,
                         scalar_exchange=source(delta,massD),connection_exchange=source(delta,gaugeD)))
    assert operator_error<2e-13 and balance_error<2e-13
    assert max(abs(row['scalar_exchange']) for row in rows)>1e-7
    assert max(abs(row['connection_exchange']) for row in rows)>1e-8
    # Finite gauge-coordinate tangent is a check of the relative Ward algebra,
    # not a substitute for spacetime divergences or quantum gauge projection.
    return dict(generators=12,full_matrix_gauge_identity_error=operator_error,
        same_record_exchange_balance_error=balance_error,rows=rows,
        continuum_Ward_proof_separate_from_local_numerical_check=True)


def joint_background_work_check():
    time,_,_,_,_,delta,_=transported();bg=background(time);p=old.point()
    tangent=(p['dg'],p['v'],p['da'],p['da0']);zero=[np.zeros_like(x) for x in bg]
    B,totalD=matrix_and_source(*bg,*tangent)
    pieces={}
    for i,name in enumerate(('metric','all_five_scalars','weak_connection','circle_connection')):
        v=[x.copy() for x in zero];v[i]=tangent[i]
        _,D=matrix_and_source(*bg,*v);pieces[name]=source(delta,D)
    dotdelta=-1j*(B@delta-delta@B)
    no_state_work=source(dotdelta,B)
    direct=source(delta,totalD)
    h=2e-5;values=[]
    for dt in (h,-h):
        U=unitary(time,dt);dd=U@delta@U.conj().T
        values.append(source(dd,Btime(time+dt)))
    finite=(values[0]-values[1])/(2*h)
    fdD=(Btime(time+h)-Btime(time-h))/(2*h)
    assert abs(no_state_work)<1e-13 and abs(sum(pieces.values())-direct)<2e-13
    assert abs(finite-direct)<2e-8 and maxabs(fdD-totalD)<2e-8
    return dict(actual_original_background_work=pieces,complete_work_rate=direct,
        covariance_commutator_work=no_state_work,energy_time_difference=finite,
        work_identity_error=abs(finite-direct),matrix_jet_error=maxabs(fdD-totalD),
        background_is_original_first_jet_diagnostic_not_full_future_Einstein_solution=True)


def run():
    checks=('relative_mode_transport_check','gauge_and_mass_source_check','joint_background_work_check')
    results={name:globals()[name]() for name in checks}
    deps=('research_note_573.md','research_note_601.md','research_note_634.md',
          'joint_dynamic_continuum_reference.py','joint_dynamic_continuum_reference_results.json',
          'joint_source_constraint_response.py','joint_source_constraint_response_results.json',
          'round732_drafts/source_feedback_entry.py','round732_drafts/source_feedback_entry_results.json')
    return dict(round=732,tests_run=3,failures=0,errors=0,checks=list(checks),results=results,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope='Original730 same-background record difference gives smooth finite-mode sources satisfying joint relative Ward identities; with731 initial correction and573 fixed-gauge hyperbolicity it defines a short-time linear relative response. Numeric tests are full local matrix identities, not a continuum Einstein solve. No absolute reference, nonlinear semiclassical existence, causal detector or GR emergence is claimed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({key:result[key] for key in ('round','tests_run','failures','errors')}))
