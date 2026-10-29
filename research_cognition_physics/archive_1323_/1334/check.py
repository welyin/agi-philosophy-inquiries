"""1334: fixed matrix clockwork, complete tree mass/weak/source checks.

No arbitrary portal matrix. NumPy only; no frozen implementation is imported.
All-n/all-M statements are analytic in proof.md; finite checks are diagnostics.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
S=np.diag(np.ones(2),-1).astype(complex)
TOL=3e-8


def norm(a): return float(np.linalg.norm(a,2))
def cp(a): return [float(np.real(a)),float(np.imag(a))]


def moments(m,r):
    out=np.zeros(6)
    for k in range(m+1):
        c=k*(k-1)/2
        out+=r**(2*k)*np.array([1,k,c,k*k,k*c,c*c])
    return out


def zmom(m,r,q):
    A,B,C,D,E,F=moments(m,r); x2=abs(q)**2
    return np.array([[A+D*x2+F*x2*x2,np.conj(q)*(B+E*x2),C*np.conj(q)**2],
                     [q*(B+E*x2),A+D*x2,B*np.conj(q)],
                     [C*q*q,B*q,A]],complex)


def zero_data(m,r,q):
    d=r*(np.eye(3)+q*S)
    powers=[np.eye(3,dtype=complex)]
    for k in range(m): powers.append(powers[-1]@d)
    x=np.vstack(powers)
    # SVD of stacked data avoids squaring its condition number for small Z eigenvalues.
    u,s,vh=np.linalg.svd(np.real_if_close(x),full_matrices=False)
    v=vh.conj().T
    z=x.conj().T@x
    inv=(v*(1/s**2))@v.conj().T
    root=(v*(1/s))@v.conj().T
    return dict(d=d,powers=powers,x=x,z=z,sv=s,v=v,inv=inv,root=root,
                v0=np.vstack(powers[::-1])@root)


def parent(m,r,q,g,M,h):
    z=zero_data(m,r,q); dim=3*(m+1)
    k=np.zeros((3*m,dim),complex)
    for j in range(m):
        k[3*j:3*j+3,3*j:3*j+3]=M*np.eye(3)
        k[3*j:3*j+3,3*j+3:3*j+6]=-M*z['d']
    p=np.eye(3,dim,k=3*m); pq=np.eye(3,dim)
    mass=np.vstack([h*g*p,k])
    ul,sv,vh=np.linalg.svd(np.real_if_close(mass),full_matrices=True)
    order=np.argsort(sv); sv=sv[order]; ul=ul[:,order]; vr=vh.conj().T[:,order]
    G=k@k.conj().T; Gi=np.linalg.inv(G); F=Gi@k
    c6=abs(g)**2*p@F.conj().T@F@p.T
    z.update(m=m,r=r,q=q,g=g,M=M,h=h,k=k,p=p,pq=pq,mass=mass,
             ul=ul,vr=vr,masses=sv,G=G,Gi=Gi,F=F,c6=c6,
             gamma2=float(np.linalg.eigvalsh(G)[0]))
    return z


def shape(z):
    t=np.trace(z).real
    return np.array([np.trace(z@z).real/t**2,np.linalg.det(z).real/t**3])


def structure_checks():
    rows=[]
    for m in (1,2,3,5):
        for r in (.7,1.,2.):
            for q in (.3+.4j,1.,2j):
                a=zero_data(m,r,q); z=zmom(m,r,q)
                Q=np.diag(np.exp(1j*np.angle(q)*np.arange(3)))
                errors=[norm(z-a['z'])/(1+norm(z)),
                        norm(z-Q@zmom(m,r,abs(q))@Q.conj().T)/(1+norm(z)),
                        norm(a['v0'].conj().T@a['v0']-np.eye(3))]
                rows.append(dict(m=m,r=r,q=cp(q),moment_relative_error=errors[0],
                    phase_relative_error=errors[1],zero_normalization_error=errors[2],
                    minimum_stacked_singular_value=float(a['sv'][-1]),
                    passed=bool(max(errors)<TOL)))
    step=1e-5
    def one(v):return shape(zmom(2,math.sqrt(v[0]),math.sqrt(v[1])))
    point=np.ones(2)
    jac=np.column_stack([(one(point+step*np.eye(2)[i])-one(point-step*np.eye(2)[i]))/(2*step)
                         for i in range(2)])
    exact=np.array([[151/1000,123/1000],[-1017/80000,-363/40000]])
    def six(v):
        tau,pu,pd,gu,gd,delta=v
        zu=zmom(2,math.sqrt(tau),math.sqrt(pu)); zd=zmom(2,math.sqrt(tau),math.sqrt(pd)*np.exp(1j*delta))
        su=shape(zu); sd=shape(zd)
        return np.r_[su,sd[0],gu**2*np.trace(np.linalg.inv(zu)).real,
                     gd**2*np.trace(np.linalg.inv(zd)).real,
                     (np.trace(zu@zd)/(np.trace(zu)*np.trace(zd))).real]
    point6=np.array([1.,1.,1.,.2,.3,math.pi/2])
    jac6=np.column_stack([(six(point6+step*np.eye(6)[i])-six(point6-step*np.eye(6)[i]))/(2*step)
                         for i in range(6)])
    zu=zmom(2,1,1); zd=zmom(2,1,1j); comm=zu@zd-zd@zu
    hu=np.linalg.inv(zu); hd=np.linalg.inv(zd); cc=hu@hd-hd@hu
    cpz=np.trace(comm@comm@comm); cph=np.trace(cc@cc@cc)
    certificate=dict(shape_jacobian=jac.tolist(),shape_jacobian_error=float(np.max(abs(jac-exact))),
        exact_shape_determinant=3093/16000000,finite_six_invariant_rank=int(np.linalg.matrix_rank(jac6,tol=1e-7)),
        smallest_six_jacobian_singular_value=float(np.linalg.svd(jac6,compute_uv=False)[-1]),
        CP_Z_trace=cp(cpz),CP_H_trace=cp(cph),CP_Z_error=float(abs(cpz-25104j)),
        CP_inverse_error=float(abs(cph-25104j/82**4)),
        passed=bool(np.max(abs(jac-exact))<1e-7 and np.linalg.matrix_rank(jac6,tol=1e-7)==6
                    and abs(cpz-25104j)<TOL and abs(cph-25104j/82**4)<TOL))
    return rows,certificate


def quad(ts):
    a,b,c=ts
    return np.array([(a-c)**2+b*b,2*b*(a+c),4*a*c],float)


def quad_min(co):
    xs=[-1.,1.]
    if co[2]>0 and -1 < -co[1]/(2*co[2]) < 1:xs.append(-co[1]/(2*co[2]))
    vals=[co[0]+co[1]*x+co[2]*x*x for x in xs]
    i=int(np.argmin(vals));return float(max(0.,vals[i])),xs[i]


def response(a,z):
    return abs(a['g'])**2*a['p']@a['k'].conj().T@a['Gi']@np.linalg.solve(a['G']-z*np.eye(len(a['G'])),a['k'])@a['p'].T


def weak_checks():
    # Published before computation; no refit across phase, weak-current, or source checks.
    u=parent(2,1,1,.2,10,1); d=parent(2,1,1,.3,13,1)
    au=(u['pq']@u['ul']).real; ad=(d['pq']@d['ul']).real
    Gu=au.T@au; Gd=ad.T@ad; rows=[]
    for delta in (0.,.4,math.pi/2,2.1,math.pi):
        Q=np.diag(np.exp(1j*delta*np.arange(3))); qa=np.kron(np.eye(3),Q)
        dc=parent(2,1,np.exp(1j*delta),.3,13,1)
        v=au.T@Q@ad; vl=v[:3,:3]
        poly_error=0.
        for i in range(9):
            for j in range(9):
                co=quad(au[:,i]*ad[:,j]);c=math.cos(delta)
                poly_error=max(poly_error,abs(abs(v[i,j])**2-(co[0]+co[1]*c+co[2]*c*c)))
        source_cov=norm(dc['mass']-qa@d['mass']@qa.conj().T)
        expected_ul=qa@d['ul']
        left_spectral=norm(dc['mass']@dc['mass'].conj().T@expected_ul-expected_ul@np.diag(d['masses']**2))
        projection=max(norm(v@v.conj().T-Gu),norm(v.conj().T@v-Gd),norm(Gu@Gu-Gu),norm(Gd@Gd-Gd))
        compensation=max(norm(Gu[:3,:3]-vl@vl.conj().T-v[:3,3:]@v[:3,3:].conj().T),
                         norm(Gd[:3,:3]-vl.conj().T@vl-v[3:,:3].conj().T@v[3:,:3]))
        quartets=[np.imag(vl[i,j]*vl[k,l]*np.conj(vl[i,l])*np.conj(vl[k,j]))
                  for i in range(3) for k in range(i+1,3) for j in range(3) for l in range(j+1,3)]
        cpq=float(max(abs(np.array(quartets))))
        errors=[source_cov,left_spectral,projection,compensation,poly_error]
        rows.append(dict(delta=delta,mass_covariance_error=source_cov,left_spectral_error=left_spectral,
            full_current_projection_error=projection,light_heavy_compensation_error=compensation,
            shared_phase_polynomial_error=poly_error,maximum_CP_quartet=cpq,
            light_current_abs=np.abs(vl).tolist(),neutral_u=Gu[:3,:3].tolist(),neutral_d=Gd[:3,:3].tolist(),
            light_heavy_strength_u=float(np.trace(v[:3,3:]@v[:3,3:].conj().T).real),
            passed=bool(max(errors)<TOL and (delta!=math.pi/2 or cpq>1e-8))))
    sectors=[]
    for label,a in [('u',u),('d',d)]:
        gamma2=a['gamma2']; lead2=abs(a['g'])**2/ a['sv']**2
        top=lead2[-1]; kappa=abs(a['g'])**2/(gamma2-top)
        residuals=[]
        for i,z in enumerate(a['masses'][:3]**2):
            f=abs(a['g'])**2*a['inv']-z*(np.eye(3)+response(a,z))
            residuals.append(norm(f@a['ul'][:3,i:i+1]))
        _,ks,vh=np.linalg.svd(a['k'],full_matrices=False)
        ve=a['p']@vh.conj().T
        spectral=abs(a['g'])**2*(ve*(1/ks**2))@ve.conj().T
        error=norm(spectral-a['c6'])
        actual2=a['masses'][:3]**2
        sectors.append(dict(department=label,complete_masses=a['masses'].tolist(),
            leading_light_masses=np.sqrt(lead2).tolist(),gamma_squared=gamma2,kappa=kappa,
            minimum_mass_gap=float(np.min(np.diff(a['masses']))),
            zero_kernel_error=norm(a['k']@a['v0']),schur_errors=residuals,
            C6_real=a['c6'].real.tolist(),C6_imag=a['c6'].imag.tolist(),
            C6_spectral_error=error,light_tree_interval_passed=bool(np.all(actual2<=lead2+TOL)
                and np.all(actual2+TOL>=lead2/(1+kappa))),
            passed=bool(top<gamma2 and min(np.linalg.eigvalsh(a['c6']))>-TOL and error<TOL
                        and max(residuals)<TOL and np.all(actual2<=lead2+TOL)
                        and np.all(actual2+TOL>=lead2/(1+kappa)))))
    return rows,sectors


def derivative(a,key):
    dm=np.zeros_like(a['mass']); dk=np.zeros_like(a['k']);dg=0j
    for j in range(a['m']):
        sl=(slice(3*j,3*j+3),slice(3*j+3,3*j+6))
        if key=='M':dk=a['k']/a['M'];break
        if key=='r':dk[sl]=-a['M']*(np.eye(3)+a['q']*S)
        if key=='q_real':dk[sl]=-a['M']*a['r']*S
        if key=='q_imag':dk[sl]=-1j*a['M']*a['r']*S
    if key=='g_real':dg=1
    if key=='g_imag':dg=1j
    dm[:3]=a['h']*dg*a['p']
    if key=='h':dm[:3]=a['g']*a['p']
    dm[3:]=dk
    dG=dk@a['k'].conj().T+a['k']@dk.conj().T
    dF=a['Gi']@(dk-dG@a['F'])
    dc=(2*np.real(np.conj(a['g'])*dg)*a['p']@a['F'].conj().T@a['F']@a['p'].T
        +abs(a['g'])**2*a['p']@(dF.conj().T@a['F']+a['F'].conj().T@dF)@a['p'].T)
    return dm,dc


def source_checks():
    rows=[]
    for label,args in [('u',[2,1,1,.2,10,1]),('d',[2,1,1j,.3,13,1])]:
        a=parent(*args)
        for key,idx,factor in [('r',1,1),('q_real',2,1),('q_imag',2,1j),
                               ('g_real',3,1),('g_imag',3,1j),('M',4,1),('h',5,1)]:
            step=1e-5; ap=list(args); am=list(args);ap[idx]+=step*factor;am[idx]-=step*factor
            plus=parent(*ap);minus=parent(*am);dm,dc=derivative(a,key)
            mass_error=norm((plus['mass']-minus['mass'])/(2*step)-dm)
            cerror=norm((plus['c6']-minus['c6'])/(2*step)-dc)
            hf=np.array([np.vdot(a['ul'][:,i],dm@a['vr'][:,i]).real for i in range(9)])
            hferr=float(np.max(abs((plus['masses']-minus['masses'])/(2*step)-hf)))
            rows.append(dict(department=label,source=key,parent_source_error=mass_error,
                C6_source_error=cerror,mass_Hellmann_Feynman_error=hferr,
                passed=bool(max(mass_error,cerror,hferr)<3e-7)))
    return rows


def portal_checks():
    rows=[]
    for M in (.1,.3,1.,10.,100.):
        a=parent(2,1,1+.3j,.2,M,1)
        upper=abs(a['g'])/a['sv'][-1]
        rows.append(dict(M=M,third_complete_mass=float(a['masses'][2]),
            exact_minmax_upper=upper,passed=bool(a['masses'][2]<=upper+TOL)))
    m=5;r=8;h=174
    examples=[];zero=[]
    for label,q,target in [('u',10,171.7),('d',8/5,2.89)]:
        a=zero_data(m,r,q);minimum=float(a['sv'][-1])
        frob=r**m/math.sqrt(3+2*m*m*q*q+(m*(m+1)/2)**2*q**4)
        g=target*minimum/h; conservative=target*frob/h
        masses=minimum/a['sv']*target
        examples.append(dict(department=label,m=m,r=r,q=q,h_GeV=h,target_highest_light_mass_GeV=target,
            sqrt_minimum_Z=minimum,all_M_necessary_g=g,Frobenius_necessary_g=conservative,
            ordinary_loop_diagnostic=g*g/(16*math.pi**2),Frobenius_loop_lower=conservative**2/(16*math.pi**2),
            normalized_leading_masses_GeV=masses.tolist(),
            leading_mass_curvature=float(masses[0]*masses[2]/masses[1]**2),
            ordinary_perturbative_interpretation_rejected=conservative>4*math.pi,
            nonperturbative_theories_excluded=False,passed=bool(minimum>=frob and conservative>4*math.pi)))
        zero.append(a)
    # Leading-only phase bounds; sample already fails the ordinary perturbative gate.
    ou=zero[0]['v'].real;od=zero[1]['v'].real
    mix={}
    for name,i,j in [('us',0,1),('ub',0,2),('cb',1,2)]:
        co=quad(ou[:,i]*od[:,j]);lower,c=quad_min(co)
        mix[name]=dict(coefficients=co.tolist(),minimum_abs=math.sqrt(lower),minimizing_cos_delta=c)
    return rows,dict(explored_before_admission=True,blind_test=False,examples=examples,
        leading_only_mixing_phase_bounds=mix,finite_controlled_SM_prediction=False,
        passed=bool(all(x['passed'] for x in examples)))


def limit_checks():
    rows=[]
    for m,r in [(2,.7),(2,1.),(3,1.),(3,2.)]:
        A,B,C,D,E,F=moments(m,r)
        gram=np.array([[F,E,C],[E,D,B],[C,B,A]])
        d2=D*F-E*E;d3=float(np.linalg.det(gram))
        target=d2**1.5/(F**1.5*math.sqrt(d3))
        a=zero_data(m,r,1024.);masses=1/a['sv']
        actual=float(masses[0]*masses[2]/masses[1]**2)
        rows.append(dict(m=m,r=r,x=1024.,gram_eigenvalues=np.linalg.eigvalsh(gram).tolist(),
            asymptotic_mass_curvature=target,finite_mass_curvature=actual,
            mass_curvature_error=abs(actual-target),
            passed=bool(d2>0 and d3>0 and abs(actual-target)<.002)))
    return rows


def compute():
    structure,local=structure_checks();weak,sectors=weak_checks();sources=source_checks()
    portals,explored=portal_checks();limits=limit_checks()
    rows=structure+[local]+weak+sectors+sources+portals+[explored]+limits
    provenance={str(p.relative_to(HERE.parent)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in [HERE/'check.py',HERE.parent/'_admission'/'round1334_admission.md']}
    return dict(round=1334,status='PASS' if all(x['passed'] for x in rows) else 'FAIL',
        provenance=provenance,structure_checks=structure,local_six_freedom_and_CP=local,
        full_parent_phase_checks=weak,controlled_sectors=sectors,source_checks=sources,
        all_M_minmax_checks=portals,explored_mass_sample=explored,fixed_r_limit_checks=limits,
        counts=dict(structure=len(structure),local=1,full_phase=len(weak),controlled_sectors=len(sectors),
                    sources=len(sources),portal=len(portals),explored=1,limits=len(limits)),
        new_axioms=0,n_three_is_input=True,autonomous_formation=False,
        exact_nonperturbative_poles=False,all_parameter_class_excluded=False,
        phase_common_to_all_readouts=True,stage_complete=False,next_round_admitted=False)


def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert set(a)==set(b),(path,'keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert len(a)==len(b),(path,'length')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=3e-8),(path,a,b)
    else:assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write-new',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();r=compute();assert r['status']=='PASS','Finite diagnostic failed; no result written.'
    p=HERE/'results.json'
    if args.write_new:
        with p.open('x',encoding='utf-8') as stream:json.dump(r,stream,ensure_ascii=False,indent=2);stream.write('\n')
    else:compare(r,json.loads(p.read_text(encoding='utf-8')))
    print(json.dumps(dict(round=1334,status='PASS',counts=r['counts'],results_written=args.write_new)))


if __name__=='__main__':main()
