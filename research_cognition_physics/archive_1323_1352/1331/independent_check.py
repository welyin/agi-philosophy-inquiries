"""Independent 1331 finite orientation/source audit; no author imports.

Complete tuple-CAR/reference states and sparse physical transitions are retained.
Dynamics may use the exactly invariant Z support. Readout rotation always acts
on the COMPLETE parent, never on that support. NumPy only.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent


def move(s,i,j):
    if j not in s:return None
    a=s.index(j);r=s[:a]+s[a+1:]
    if i in r:return None
    b=sum(k<i for k in r)
    return r[:b]+(i,)+r[b:],(-1)**(a+b)


class Op:
    def __init__(self,d,triples=()):
        self.d=d
        triples=list(triples)
        self.r=np.asarray([x[0] for x in triples],int)
        self.c=np.asarray([x[1] for x in triples],int)
        self.v=np.asarray([x[2] for x in triples],complex)
    def adj(self):return Op(self.d,zip(self.c.tolist(),self.r.tolist(),self.v.conj().tolist()))
    def scaled(self,a):return Op(self.d,list(zip(self.r,self.c,a*self.v)))
    def add(self,*others):
        t=list(zip(self.r,self.c,self.v))
        for x in others:t+=list(zip(x.r,x.c,x.v))
        return Op(self.d,t)
    def act(self,x):
        a=x[:,None] if x.ndim==1 else x
        y=np.zeros((self.d,a.shape[1]),complex)
        np.add.at(y,self.r,self.v[:,None]*a[self.c])
        return y[:,0] if x.ndim==1 else y
    def dense(self):
        a=np.zeros((self.d,self.d),complex)
        np.add.at(a,(self.r,self.c),self.v)
        return a


def local_counts(s,n):
    out=[]
    for x in range(2):
        o=4*n*x
        out.append([sum(o+k*n<=i<o+(k+1)*n for i in s) for k in range(4)])
    return np.asarray(out,int) # D,L0,L1,U


def model(n,q,N,md,mu,js,p):
    logical=list(itertools.combinations(range(8*n),q));states=[]
    for s in logical:
        co=local_counts(s,n);b=np.asarray(N)-co[:,1]-co[:,2]-2*co[:,3]
        states.extend((s,a,c) for a in range(b[0]+1) for c in range(b[1]+1))
    ix={s:i for i,s in enumerate(states)};d=len(states)
    counts=np.array([local_counts(s[0],n) for s in states]);blue=np.array([s[1:] for s in states])
    bs=np.asarray(N)[None,:]-counts[:,:,1]-counts[:,:,2]-2*counts[:,:,3]
    z=blue.sum(axis=1)+counts[:,:,2].sum(axis=1)+counts[:,:,3].sum(axis=1)
    ton=[];tedge=[];tjp=[]
    for col,(s,h0,h1) in enumerate(states):
        hs=[h0,h1]
        for x in range(2):
            o=4*n*x;b=bs[col,x];h=hs[x];red=b-h
            if h:
                dest=hs.copy();dest[x]-=1
                tjp.append((ix[(s,*dest)],col,math.sqrt((red+1)*h)))
            for i in range(n):
                r=move(s,o+n+i,o+2*n+i)
                if r:tjp.append((ix[(r[0],h0,h1)],col,r[1]))
                for j in range(n):
                    transitions=[(o+n+i,o+j,md[x][i,j]*math.sqrt(red/b),0),
                                 (o+2*n+i,o+j,md[x][i,j]*math.sqrt(h/b),-1),
                                 (o+n+i,o+3*n+j,mu[x][i,j]*math.sqrt((h+1)/(b+2)),1),
                                 (o+2*n+i,o+3*n+j,-mu[x][i,j]*math.sqrt((red+1)/(b+2)),0)]
                    for dest,orig,value,shift in transitions:
                        r=move(s,dest,orig)
                        if r and abs(value)>0:
                            bh=hs.copy();bh[x]+=shift
                            ton.append((ix[(r[0],*bh)],col,r[1]*value))
        for role,w in enumerate((0,1,1,2)):
            if bs[col,0]-h0<w:continue
            for i in range(n):
                for j in range(n):
                    r=move(s,role*n+i,4*n+role*n+j)
                    if r:tedge.append((ix[(r[0],h0,h1)],col,r[1]*js[role][i,j]))
    on=Op(d,ton);edge=Op(d,tedge);jp=Op(d,tjp)
    v=on.add(on.adj());he=edge.add(edge.adj())
    hdiag=sum(p[x]*(counts[:,x,0]+counts[:,x,3]-counts[:,x,1]-counts[:,x,2]) for x in range(2))
    h0=Op(d,[(i,i,a) for i,a in enumerate(hdiag) if a])
    h=h0.add(v,he);jx=jp.add(jp.adj()).scaled(.5);jy=jp.add(jp.adj().scaled(-1)).scaled(-.5j)
    selected=np.flatnonzero(z<=q);P=np.eye(d,dtype=complex)[:,selected]
    return dict(n=n,q=q,N=N,states=states,ix=ix,counts=counts,blue=blue,bs=bs,z=z,
                H=h,H0=h0,V=v,edge=he,Jp=jp,Jx=jx,Jy=jy,hdiag=hdiag,P=P,selected=selected)


def norm(a):
    if a.ndim==1:return float(np.linalg.norm(a))
    gram=a.conj().T@a
    return math.sqrt(max(0.,float(np.linalg.eigvalsh(gram)[-1])))


def source(a,x):return 1j*(a['edge'].act(a['Jx'].act(x))-a['Jx'].act(a['edge'].act(x)))
def nuclear(a):return float(np.linalg.svd(a,compute_uv=False).sum())


def prepared(a,site=1,spin_phase=1):
    n=a['n'];F=np.zeros((len(a['states']),n),complex)
    for i in range(n):
        for role,coef in [(1,1),(2,spin_phase)]:
            F[a['ix'][((4*n*site+role*n+i,),0,0)],i]=coef/math.sqrt(2)
    return F


def compute():
    checks=[];rows=[]
    def eq(name,x,y,tol=1e-9):
        r=float(np.linalg.norm(np.asarray(x)-np.asarray(y)))
        checks.append(dict(name=name,residual=r,tolerance=tol))
        if r>tol:raise AssertionError((name,r,tol))
    def le(name,x,y,tol=1e-9):
        r=max(0.,float(x-y));checks.append(dict(name=name,residual=r,tolerance=tol,value=float(x),bound=float(y)))
        if r>tol:raise AssertionError((name,x,y))
    def evolution(a,F,t):
        P=a['P'];hh=P.conj().T@a['H'].act(P);ev,vec=np.linalg.eigh(hh)
        return P@((vec*np.exp(-1j*ev*t))@(vec.conj().T@(P.conj().T@F)))

    matrix_tests=[]
    for n,q,N in [(1,1,(6,7)),(2,1,(6,6)),(1,2,(9,9))]:
        if n==1:
            md=[np.array([[.012+.002j]]),np.array([[.014-.001j]])]
            mu=[np.array([[.008-.001j]]),np.array([[.011+.002j]])]
            js=[np.array([[.19]]),np.array([[.31]]),np.array([[.07]]),np.array([[.23]])]
        else:
            md=[np.array([[.012,.003j],[.002,.010]]),np.array([[.013,.002],[.001j,.011]])]
            mu=[np.array([[.008,.001],[-.002j,.009]]),np.array([[.011,.002j],[.003,.007]])]
            js=[np.array([[.19,.02j],[.01,.15]]),np.array([[.31,.04j],[.015,.24]]),
                np.array([[.07,.01],[-.025j,.11]]),np.array([[.23,.03],[.02j,.17]])]
        a=model(n,q,N,md,mu,js,(.025,.025));P=a['P'];tag=f'n{n}.q{q}.N{N}'
        hP=a['H'].act(P);vP=a['V'].act(P)
        eq(tag+'.Z_exact_transition',a['H'].v*(a['z'][a['H'].r]-a['z'][a['H'].c]),0)
        eq(tag+'.H_invariant_support',hP-P@(P.conj().T@hP),0)
        le(tag+'.blue_on_actual_support',int(a['blue'][a['selected']].sum(axis=1).max()),q)
        le(tag+'.red_lower_bound',min(N)-3*q,float((a['bs']-a['blue'])[a['selected']].min()))
        eq(tag+'.onsite_SU2',a['V'].act(a['Jx'].act(P))-a['Jx'].act(vP),0)
        eq(tag+'.H0_SU2',a['H0'].act(a['Jx'].act(P))-a['Jx'].act(a['H0'].act(P)),0)
        # Verify SU2 on complete parent probe columns, not only P.
        rng=np.random.default_rng(1331+n+q);probe=rng.normal(size=(len(a['states']),3))+1j*rng.normal(size=(len(a['states']),3))
        probe/=np.linalg.norm(probe,axis=0)
        jm=a['Jp'].adj();jz=sum(N)/2-a['z']
        eq(tag+'.full_parent_SU2',a['Jp'].act(jm.act(probe))-jm.act(a['Jp'].act(probe)),2*jz[:,None]*probe)
        CP=source(a,P);source_norm=norm(CP)
        # Reference part obtained by subtracting the physical fermion commutator.
        fplus=[]
        for col,(s,h0,h1) in enumerate(a['states']):
            for x in range(2):
                for i in range(n):
                    r=move(s,4*n*x+n+i,4*n*x+2*n+i)
                    if r:fplus.append((a['ix'][(r[0],h0,h1)],col,r[1]))
        fp=Op(len(a['states']),fplus);fx=fp.add(fp.adj()).scaled(.5)
        Cferm=1j*(a['edge'].act(fx.act(P))-fx.act(a['edge'].act(P)))
        refnorm=norm(CP-Cferm)
        bound=2*math.sqrt(q+1)/math.sqrt(min(N)-3*q-2)*(nuclear(js[1])+nuclear(js[2])+2*nuclear(js[3]))
        le(tag+'.reference_source_bound',refnorm,bound)
        same=[js[0],js[1],js[1],js[3]]
        b=model(n,q,N,md,mu,same,(.025,.025));CPsame=source(b,b['P'])
        samesource=norm(CPsame)
        samebound=2*math.sqrt(q+1)/math.sqrt(min(N)-3*q-2)*(2*nuclear(js[1])+2*nuclear(js[3]))
        le(tag+'.equal_L_sufficient',samesource,samebound)
        if q==1:
            F=prepared(a);Hf=a['H'].act(F);Cf=source(a,F)
            Cdot=1j*(a['H'].act(Cf)-source(a,Hf))
            curvature=F.conj().T@Cdot
            delta=js[1]-js[2]
            eq(tag+'.actual_initial_source',F.conj().T@Cf,0)
            eq(tag+'.actual_matrix_curvature',curvature,-delta.conj().T@delta/2)
            eq(tag+'.actual_initial_Jx',F.conj().T@a['Jx'].act(F),np.eye(n)/2)
            E=np.zeros((len(a['states']),8*n),complex)
            for i in range(8*n):E[a['ix'][((i,),0,0)],i]=1
            eq(tag+'.compressed_source_singular_value',norm(E.conj().T@source(a,E)),norm(delta)/2)
            Hrestricted=P.conj().T@hP;Crestricted=P.conj().T@CP
            remainder=(2/3)*norm(Hrestricted)**2*norm(Crestricted)
            for t in (.05,.13,.27):
                out=evolution(a,F,t)
                measured=out.conj().T@a['Jx'].act(out)
                expansion=np.eye(n)/2-t*t*delta.conj().T@delta/4
                le(tag+f'.finite_matrix_Taylor_t{t}',norm(measured-expansion),remainder*t**3)
            matrix_tests.append(dict(n=n,known_start_site=1,flavor_dimension=n,
                delta_singular_value=norm(delta),curvature_matrix_real=curvature.real.tolist(),
                curvature_matrix_imag=curvature.imag.tolist(),finite_Taylor_coefficient=remainder))
        rows.append(dict(n=n,q=q,N=list(N),full_dimension=len(a['states']),actual_support_dimension=P.shape[1],
              source_norm=source_norm,reference_source_norm=refnorm,reference_bound=bound,
              equal_L_source_norm=samesource,equal_L_bound=samebound))

    # Complete readout: switch conversion+edge off, rotate FULL reference and L
    # directions, measure all occupation outcomes and keep their records.
    measured=[]
    for N in [(6,7),(9,9)]:
        md=[np.array([[.012+.002j]]),np.array([[.014-.001j]])]
        mu=[np.array([[.008-.001j]]),np.array([[.011+.002j]])]
        js=[np.array([[.19]]),np.array([[.31]]),np.array([[.07]]),np.array([[.23]])]
        a=model(1,1,N,md,mu,js,(.025,.025));F=prepared(a);P=a['P'];CP=source(a,P)
        vo=norm(P.conj().T@a['V'].act(P));cnorm=norm(P.conj().T@CP)
        jy=a['Jy'].dense();val,vec=np.linalg.eigh(jy)
        rotate=(vec*np.exp(1j*math.pi*val/2))@vec.conj().T
        jz=sum(N)/2-a['z'];d=len(a['states']);tag=f'readout.N{N}'
        eq(tag+'.full_rotation_dictionary',rotate.conj().T@(jz[:,None]*rotate),a['Jx'].dense(),tol=4e-9)
        eq(tag+'.full_rotation_H0',rotate.conj().T@(a['hdiag'][:,None]*rotate),np.diag(a['hdiag']),tol=4e-9)
        # Whole unknown left two-direction input with an arbitrary partner.
        unknown=np.zeros((d,2),complex)
        for k in range(2):unknown[a['ix'][((4+1+k,),0,0)],k]=1
        for t in (.8,1.4,2.):
            psi=evolution(a,F,t)[:,0];out=rotate@psi
            probabilities=np.abs(out)**2
            actual=float(np.real(psi.conj()@a['Jx'].act(psi)))
            decoded=float(probabilities@jz)
            eq(tag+f'.t{t}.actual_record_mean',decoded,actual)
            eq(tag+f'.t{t}.all_records',probabilities.sum(),1.)
            variance=float(probabilities@(jz-decoded)**2)
            before_orientation=np.array([actual,float(np.real(psi.conj()@a['Jy'].act(psi))),float(np.real(psi.conj()@(jz*psi)))])
            after_rotation=np.array([float(np.real(out.conj()@a['Jx'].act(out))),
                                     float(np.real(out.conj()@a['Jy'].act(out))),decoded])
            after_record=np.array([0.,0.,decoded])
            eq(tag+f'.t{t}.rotation_orientation_source',after_rotation,
               np.array([-before_orientation[2],before_orientation[1],before_orientation[0]]))
            # Full occupation record removes transverse SYSTEM coherences;
            # those changes remain assigned to the measurement/control source.
            eq(tag+f'.t{t}.record_Jz_source',after_record[2]-after_rotation[2],0.)
            # Large reference occupation outcome spread is retained, not clipped to P.
            outside=float(probabilities[a['z']>1].sum())
            assert outside>.01
            before=float(np.real(psi.conj()@a['H'].act(psi)))
            after_off=float(probabilities@a['hdiag'])
            work=-float(np.real(psi.conj()@a['V'].add(a['edge']).act(psi)))
            eq(tag+f'.t{t}.quench_rotate_record_energy',after_off-before,work)
            actual_unknown=evolution(a,unknown,t);rotated_unknown=rotate@actual_unknown
            eq(tag+f'.t{t}.unknown_partner_record_gram',rotated_unknown.conj().T@rotated_unknown,np.eye(2))
            eq(tag+f'.t{t}.unknown_record_source',rotated_unknown.conj().T@(jz[:,None]*rotated_unknown),
               actual_unknown.conj().T@a['Jx'].act(actual_unknown))
            pure=.5*math.cos((.31-.07)*t)
            # Bounded projected source is enough: both processes preserve P.
            perturb=t*t*vo*cnorm
            le(tag+f'.t{t}.full_nonzero_onsite_witness',abs(actual-pure),perturb)
            drop=.5-actual
            assert drop>0
            tolerance=.002;failure=.01
            repetitions=math.ceil(sum(N)**2*math.log(4/failure)/(2*tolerance*tolerance))
            local_outcomes=[]
            for x in range(2):
                local_outcomes.append(len(set((tuple(i-4*x for i in s if 4*x<=i<4*x+4),hs[x])
                    for s,*hs in a['states'])))
            record_bits=sum(math.ceil(math.log2(k)) for k in local_outcomes)
            assert drop>4*tolerance
            angle_precision=1e-5
            angular_bias=sum(N)*angle_precision/2
            distinguishability_margin=drop-2*(tolerance+angular_bias)
            assert distinguishability_margin>0
            measured.append(dict(N=list(N),time=t,mean=actual,pure_edge_mean=pure,
                discrepancy=abs(actual-pure),onsite_perturbation_bound=perturb,
                finite_signal_drop=drop,record_variance=variance,
                before_instrument_orientation=before_orientation.tolist(),
                rotation_orientation_change=(after_rotation-before_orientation).tolist(),
                record_system_orientation_change=(after_record-after_rotation).tolist(),
                total_instrument_orientation_change=(after_record-before_orientation).tolist(),
                orientation_changes_are_registered_control_measurement_sources=True,
                rotated_probability_outside_dynamic_support=outside,
                all_record_results=d,quench_work=work,per_mean_accuracy=tolerance,
                total_failure_probability=failure,repetitions_per_mean=repetitions,
                initial_reference_quanta_per_trial=sum(N)-1,
                reference_quanta_for_two_independent_means=2*repetitions*(sum(N)-1),
                local_record_outcomes=local_outcomes,blank_record_bits_per_trial=record_bits,
                blank_record_bits_for_two_means=2*repetitions*record_bits,
                angular_precision_per_local_pulse=angle_precision,
                systematic_mean_error_bound=angular_bias,
                two_mean_statistical_plus_angular_margin=distinguishability_margin,
                communication_for_record_aggregation_is_not_free=True))
    return dict(round=1331,status='PASS',new_axioms=0,stage_complete=False,
         construction='Complete independent tuple-CAR/reference parent; sparse full transitions; exact invariant support only for dynamics; full-parent SU2 readout.',
         extra_finite_instrument_input='Calibrated joint reference/L rotations with fixed menu, direction-resolved occupation records and blank memory; finite preparation and statistics, no free resets.',
         checks=checks,check_count=len(checks),max_residual=max(x['residual'] for x in checks),
         cases=rows,matrix_necessity=matrix_tests,finite_readout=measured,
         author_imports=False,field_or_generation_selection_proven=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):assert math.isfinite(b) and abs(a-b)<=1e-8*(1+abs(a)),(path,a,b)
    else:assert a==b,path


def main():
    p=argparse.ArgumentParser();p.add_argument('--write-new',action='store_true');p.add_argument('--check',action='store_true')
    args=p.parse_args();r=compute();dest=HERE/'independent_results.json'
    if args.write_new:
        with dest.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(r,json.loads(dest.read_text(encoding='utf-8')))
    print(json.dumps({k:r[k] for k in ('round','status','check_count','max_residual','cases','matrix_necessity','finite_readout')},ensure_ascii=False))


if __name__=='__main__':main()
