"""1326 independent finite driven-reference, dressed-band and source checks.

No author implementation is imported. All matrices below are rebuilt from
occupation transitions; large reference spectra are independently assembled
from every fixed transverse-number Jacobi block.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
RESULT=HERE/'independent_results.json'
CHIRAL=np.array([1.,-1.,-1.,1.])
COST=np.array([0,1,1,2])


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def op(A): return float(np.linalg.norm(A,2))
def adj(A): return A.conj().T


class Checks:
    def __init__(self): self.rows=[]
    def near(self,name,a,b=0.,tol=2e-9):
        x=float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
        if not np.isfinite(x) or x>tol: raise AssertionError((name,x,tol))
        self.rows.append(dict(name=name,residual=x,tolerance=tol))
    def le(self,name,x,y,tol=2e-10):
        if not np.isfinite(x+y) or x>y+tol: raise AssertionError((name,x,y))
        self.rows.append(dict(name=name,value=float(x),upper=float(y),tolerance=tol))
    def yes(self,name,x):
        if not x: raise AssertionError(name)
        self.rows.append(dict(name=name,passed=True))


def reference(B,f):
    """Full reference spectrum: all N1 sectors, not only its ground sector."""
    spectra=[];zero=None
    for b in range(B+1):
        a=np.arange(B-b+1)
        S=np.diag(-np.sqrt(a[1:]),1)+np.diag(-np.sqrt(a[1:]),-1)
        R=np.diag(a+b)+f*S
        ev,Q=np.linalg.eigh(R)
        spectra.extend((float(e),b,j) for j,e in enumerate(ev))
        if b==0: zero=(ev,Q,S)
    spectra.sort()
    ev,Q,S=zero
    chi=Q[:,0].copy()
    if chi[0]<0: chi=-chi
    chip=Q[:,1:]@((Q[:,1:].T@S@chi)/(ev[0]-ev[1:]))
    return dict(ground=spectra[0][0],gap=spectra[1][0]-spectra[0][0],
                sector=spectra[0][1],chi=chi,chip=chip,
                derivative=float(chi@S@chi),eigenvalues=[x[0] for x in spectra],
                first_excited_sector=spectra[1][1])


def basis(B,restricted):
    return [(m,a,b) for m in range(4) for b in range(B+1)
            for a in range(B-b+1)
            if not restricted or b+int(m in (2,3)) in (0,1)]


def build(B,gd,gu,normalized,restricted):
    """Exact Nb<=B finite model, with boundary creation set to zero."""
    states=basis(B,restricted);idx={s:i for i,s in enumerate(states)}
    d=len(states);S=np.zeros((d,d));VD=np.zeros((d,d));VU=np.zeros((d,d))
    nb=np.array([a+b for m,a,b in states],float)
    pdiag=np.array([CHIRAL[m] for m,a,b in states])
    for j,(m,a,b) in enumerate(states):
        if a:
            i=idx[m,a-1,b];S[i,j]=S[j,i]=-math.sqrt(a)
        if m==0 and a+b:
            den=math.sqrt(a+b) if normalized else 1.
            for alpha,count in enumerate((a,b)):
                if count:
                    target=(alpha+1,a-int(alpha==0),b-int(alpha==1))
                    i=idx[target];VD[i,j]=VD[j,i]=math.sqrt(count)/den
        if m==3 and a+b<B:
            den=math.sqrt(a+b+2) if normalized else 1.
            for alpha,beta,sign in [(0,1,1),(1,0,-1)]:
                i=idx[alpha+1,a+int(beta==0),b+int(beta==1)]
                VU[i,j]=VU[j,i]=sign*math.sqrt((a,b)[beta]+1)/den
    V=gd*VD+gu*VU
    return dict(B=B,states=states,S=S,VD=VD,VU=VU,V=V,nb=nb,pdiag=pdiag,
                normalized=normalized,gd=gd,gu=gu,
                g=max(abs(gd),abs(gu))*(1. if normalized else math.sqrt(B+2)))


def embedding(model,coeff):
    E=np.zeros((len(model['states']),4))
    for row,(m,a,b) in enumerate(model['states']):
        if b==0:E[row,m]=coeff[a]
    return E


def inverse_sqrt(A):
    a,U=np.linalg.eigh(A)
    if a[0]<=0:raise AssertionError(('positive_encoding_Gram',a))
    return (U*(1/np.sqrt(a)))@adj(U),a,U


def band(model,f,p):
    ref=reference(model['B'],f)
    H=np.diag(model['nb']+p*model['pdiag'])+f*model['S']+model['V']
    ev,Q=np.linalg.eigh(H);low=Q[:,:4];high=Q[:,4:]
    E=embedding(model,ref['chi']);Ep=embedding(model,ref['chip'])
    Pi=low@low.T;G=E.T@Pi@E;R,ge,Gq=inverse_sqrt(G)
    J=Pi@E@R;hd=J.T@H@J
    # Kato derivative of the whole cluster, insensitive to its internal basis.
    bridge=(high.T@model['S']@low)/(ev[:4][None,:]-ev[4:,None])
    Pip=high@bridge@low.T+low@bridge.T@high.T
    Gp=Ep.T@Pi@E+E.T@Pip@E+E.T@Pi@Ep
    roots=np.sqrt(ge)
    divided=-1/(roots[:,None]*roots[None,:]*(roots[:,None]+roots[None,:]))
    Rp=Gq@(divided*(Gq.T@Gp@Gq))@Gq.T
    Jp=Pip@E@R+Pi@Ep@R+Pi@E@Rp
    hdp=Jp.T@H@J+J.T@model['S']@J+J.T@H@Jp
    source=J.T@model['S']@J;connection=J.T@Jp
    h0=E.T@H@E
    leak=model['V']@E-E@(E.T@model['V']@E)
    b=op(leak);delta=ref['gap']-2*(abs(p)+model['g'])
    if delta<=0:raise AssertionError(('band_gap',delta))
    return dict(model=model,f=f,p=p,ref=ref,H=H,eig=ev,vectors=Q,E=E,Ep=Ep,
                Pi=Pi,Pip=Pip,J=J,Jp=Jp,hd=hd,hdp=hdp,source=source,
                connection=connection,encoding_gram=G,h0=h0,
                low=ev[:4]-ref['ground'],
                leading=np.linalg.eigvalsh(h0)-ref['ground'],
                delta=delta,b=b,epsilon=b*b/delta,
                conservative_epsilon=model['g']**2/delta)


def evolve(H,t,X):
    ev,Q=np.linalg.eigh(H)
    return (Q*np.exp(-1j*t*ev))@(Q.T.conj()@X)


def inspect_band(c,label,x):
    model=x['model'];ref=x['ref'];E=x['E'];J=x['J'];H=x['H'];d=len(H)
    c.yes(label+'.reference_ground_sector',ref['sector']==0)
    c.near(label+'.E_isometry',E.T@E,np.eye(4))
    c.near(label+'.J_isometry',J.T@J,np.eye(4))
    c.yes(label+'.Gram_strictly_positive',np.linalg.eigvalsh(x['encoding_gram'])[0]>0)
    c.near(label+'.exact_invariant_low_band',H@J,J@x['hd'])
    c.near(label+'.cluster_projector_derivative',
           x['Pip']@x['Pi']+x['Pi']@x['Pip'],x['Pip'])
    c.near(label+'.antihermitian_connection',x['connection']+x['connection'].T)
    c.near(label+'.actual_source_identity',
           x['source'],x['hdp']-(x['hd']@x['connection']-x['connection']@x['hd']))
    for j in range(4):
        c.le(label+f'.ordered{j}.upper',x['low'][j],x['leading'][j])
        c.le(label+f'.ordered{j}.lower',x['leading'][j]-x['epsilon'],x['low'][j])
    c.le(label+'.actual_next_gap',x['delta'],x['eig'][4]-x['eig'][3])
    c.le(label+'.global_complement_lower',
         ref['ground']+ref['gap']-abs(x['p'])-model['g'],x['eig'][4])
    c.le(label+'.angle_bound',op((np.eye(d)-E@E.T)@x['Pi']),x['b']/x['delta'])
    c.le(label+'.canonical_encoding_distance',op(J-E),math.sqrt(2)*x['b']/x['delta'])
    t=.73;U=evolve(H,t,J);Ud=J@evolve(x['hd'],t,np.eye(4))
    c.near(label+'.full_exact_encoded_dynamics',U,Ud,tol=1e-8)
    rng=np.random.default_rng(1326)
    Z=rng.normal(size=(4,3))+1j*rng.normal(size=(4,3));Z/=np.linalg.norm(Z)
    c.near(label+'.unknown_partner_dynamics',U@Z,Ud@Z,tol=1e-8)
    c.near(label+'.partner_marginal',adj(U@Z)@(U@Z),adj(Z)@Z)
    bareU=evolve(H,t,E)
    c.le(label+'.product_preparation_error',op(bareU-E@evolve(x['hd'],t,np.eye(4))),
         2*op(J-E),tol=1e-8)
    # Independent source decomposition; no source is re-tuned on resonance.
    c.near(label+'.complete_source_energy',
           H,np.diag(model['nb'])+x['p']*np.diag(model['pdiag'])
           +x['f']*model['S']+model['gd']*model['VD']+model['gu']*model['VU'])
    K=np.array([a+b+COST[m] for m,a,b in model['states']],float)
    charge=np.array([.23+.37*COST[m]+.37*(a+b)
                     for m,a,b in model['states']])
    Kflow=1j*(H*K[None,:]-K[:,None]*H)
    supplied=1j*x['f']*(model['S']*K[None,:]-K[:,None]*model['S'])
    c.near(label+'.K_exchange_entirely_in_control_source',Kflow,supplied)
    c.near(label+'.charge_operator_dictionary',charge,.23+.37*K)
    c.near(label+'.charge_exchange_source',
           1j*(H*charge[None,:]-charge[:,None]*H),.37*supplied)
    c.near(label+'.conserved_total_energy',adj(U)@H@U,J.T@H@J)
    hbare=np.diag(model['nb']+x['p']*model['pdiag'])
    changes=[adj(U)@O@U-J.T@O@J for O in (hbare,x['f']*model['S'],model['V'])]
    c.near(label+'.all_energy_changes_close',sum(changes))
    # A nonconstant coordinate chart confirms that the connection cannot be
    # dropped merely because an instantaneous diagonal chart exists.
    G=np.array([[0,.3j,0,0],[-.3j,0,0,0],[0,0,.2,.1],[0,0,.1,-.2]],complex)
    R=evolve(G,x['f'],np.eye(4))
    ht=adj(R)@x['hd']@R;At=adj(R)@x['connection']@R-1j*G
    htprime=adj(R)@x['hdp']@R+1j*(G@ht-ht@G)
    c.near(label+'.source_covariant_chart',
           adj(R)@x['source']@R,htprime-(ht@At-At@ht))
    c.yes(label+'.nonzero_chart_connection_term',op(ht@At-At@ht)>1e-5)
    return dict(f=x['f'],dimension=d,reference_gap=ref['gap'],
                ground=ref['ground'],actual_band_gap=float(x['eig'][4]-x['eig'][3]),
                delta=x['delta'],global_norm_bound=model['g'],offdiagonal_norm=x['b'],
                epsilon=x['epsilon'],conservative_epsilon=x['conservative_epsilon'],
                exact_low=x['low'].tolist(),first_block=x['leading'].tolist(),
                encoding_distance=op(J-E),
                actual_source=np.diag(x['vectors'][:,:4].T@model['S']@x['vectors'][:,:4]).tolist(),
                reference_source=ref['derivative'])


def calibration(B,gd=.013,gu=.011):
    norm=build(B,gd,gu,True,True)
    r=reference(B,1.);E=embedding(norm,r['chi'])
    mean=E.T@norm['V']@E
    unit=build(B,1.,1.,False,True)
    unitmean=E.T@unit['V']@E
    cd=float(mean[1,0]/unitmean[1,0]);cu=float(mean[2,3]/unitmean[2,3])
    return norm,build(B,cd,cu,False,True)


def small_parent_checks(c):
    B=3;gd=.013;gu=.011;f=.8;p=.005
    full=build(B,gd,gu,True,False);red=build(B,gd,gu,True,True)
    idx={s:i for i,s in enumerate(full['states'])}
    Q=np.zeros((len(full['states']),len(red['states'])))
    for col,s in enumerate(red['states']):Q[idx[s],col]=1.
    C=np.array([b+int(m in (2,3)) for m,a,b in full['states']])
    K=np.array([a+b+COST[m] for m,a,b in full['states']])
    H=np.diag(full['nb']+p*full['pdiag'])+f*full['S']+full['V']
    Hr=np.diag(red['nb']+p*red['pdiag'])+f*red['S']+red['V']
    c.near('small.C_conserved',(C[:,None]-C[None,:])*H)
    c.yes('small.drive_breaks_old_K',op((K[:,None]-K[None,:])*H)>1.)
    c.near('small.complete_C01_embedding',H@Q,Q@Hr)
    c.near('small.source_embedding',full['S']@Q,Q@red['S'])
    # Ground reference spectrum computed two ways (full occupation and blocks).
    rstates=[(a,b) for b in range(B+1) for a in range(B-b+1)]
    ridx={s:i for i,s in enumerate(rstates)}
    R=np.diag([a+b for a,b in rstates]).astype(float)
    for j,(a,b) in enumerate(rstates):
        if a:R[ridx[a-1,b],j]=R[j,ridx[a-1,b]]=-f*math.sqrt(a)
    ref=reference(B,f)
    c.near('small.all_Jacobi_spectra',np.linalg.eigvalsh(R),ref['eigenvalues'])
    xf=band(full,f,p);xr=band(red,f,p)
    c.near('small.complete_low_spectra',xf['low'],xr['low'])
    c.near('small.complete_low_projector',xf['Pi'],Q@xr['Pi']@Q.T)
    c.near('small.complete_dressed_code',xf['J'],Q@xr['J'])
    c.near('small.complete_source_matrix',xf['source'],xr['source'])
    c.near('small.complete_reduced_source_connection',xf['connection'],xr['connection'])
    c.le('small.normalized_global_bound',op(full['V']),full['g'])
    # At zero drive, the old K<=B safety window is an exact invariant block.
    Hres=np.diag(1.+.4*K)+full['V']
    c.near('small.old_resonant_K',(K[:,None]-K[None,:])*Hres)
    safety=np.flatnonzero(K<=B)
    outside=np.flatnonzero(K>B)
    c.near('small.old_window_no_edge_leak',Hres[np.ix_(outside,safety)])
    row=inspect_band(c,'small',xr)
    return row


def derivative_and_instrument(c,label,x):
    model=x['model'];f=x['f'];p=x['p'];step=2e-4
    plus=band(model,f+step,p);minus=band(model,f-step,p)
    c.near(label+'.Jprime_finite_difference',
           (plus['J']-minus['J'])/(2*step),x['Jp'],tol=2e-7)
    c.near(label+'.hdprime_finite_difference',
           (plus['hd']-minus['hd'])/(2*step),x['hdp'],tol=2e-7)
    physical_derivative=np.diag(x['vectors'][:,:4].T@model['S']@x['vectors'][:,:4])
    c.near(label+'.full_Hellmann_Feynman',
           (plus['eig'][:4]-minus['eig'][:4])/(2*step),physical_derivative,tol=2e-7)
    c.near(label+'.excitation_Hellmann_Feynman',
           (plus['low']-minus['low'])/(2*step),
           physical_derivative-x['ref']['derivative'],tol=2e-7)
    # Full bounded source readout with both outcomes, partner and energy kept.
    scale=2*math.sqrt(model['B'])
    ev,Q=np.linalg.eigh(model['S'])
    c.le(label+'.source_operator_bound',float(max(abs(ev))),scale)
    Kp=(Q*np.sqrt((1+ev/scale)/2))@Q.T
    Km=(Q*np.sqrt((1-ev/scale)/2))@Q.T
    c.near(label+'.source_instrument_completeness',Kp.T@Kp+Km.T@Km,np.eye(len(ev)))
    rng=np.random.default_rng(926);z=rng.normal(size=(4,2))+1j*rng.normal(size=(4,2))
    z/=np.linalg.norm(z);psi=x['J']@z;outp=Kp@psi;outm=Km@psi
    pp=float(np.linalg.norm(outp)**2);pm=float(np.linalg.norm(outm)**2)
    source=float(np.trace(adj(psi)@model['S']@psi).real)
    c.near(label+'.source_probabilities',pp+pm,1.)
    c.near(label+'.actual_source_readout',(pp-pm)*scale,source)
    c.near(label+'.measurement_keeps_unknown_partner',
           adj(outp)@outp+adj(outm)@outm,adj(psi)@psi)
    before=float(np.trace(adj(psi)@x['H']@psi).real)
    after=float(np.trace(adj(outp)@x['H']@outp+adj(outm)@x['H']@outm).real)
    c.yes(label+'.measurement_energy_not_free',abs(after-before)>1e-8)
    return dict(f=f,finite_difference_step=step,source=source,probabilities=[pp,pm],
                source_bound=scale,system_measurement_energy_change=after-before,
                required_apparatus_energy_change=before-after,
                readout_permission_and_apparatus_are_inputs=True)


def compute():
    c=Checks();small=small_parent_checks(c)
    B=66;p=.005;norm,bare=calibration(B)
    models={'normalized':norm,'linear_once_calibrated':bare};rows=[];bands={}
    c.yes('large.complete_C01_dimension',len(norm['states'])==6*B+4)
    c.yes('large.includes_old_L64_window',B>=64+2)
    for label,model in models.items():
        c.le(label+'.actual_global_norm',op(model['V']),model['g'])
        for f in [1.,2.]:
            x=band(model,f,p);bands[label,f]=x
            rows.append(dict(menu=label,**inspect_band(c,label+f'.f{f:g}',x)))
    c.near('calibration.only_first_block_fixed',
           bands['normalized',1.]['h0'],bands['linear_once_calibrated',1.]['h0'])
    c.yes('calibration.full_parent_spectra_not_claimed_equal',
          op(np.diag(bands['normalized',1.]['low']-bands['linear_once_calibrated',1.]['low']))>1e-7)
    intervals={}
    for label in models:
        a=bands[label,1.];b=bands[label,2.]
        central=b['leading'][-1]-a['leading'][-1]
        interval=[float(central-b['epsilon']),float(central+a['epsilon'])]
        exact=float(b['low'][-1]-a['low'][-1])
        c.le(label+'.finite_response_lower',interval[0],exact)
        c.le(label+'.finite_response_upper',exact,interval[1])
        intervals[label]=dict(response_interval=interval,exact_response=exact,
                              heldout_energy_interval=[float(b['leading'][-1]-b['epsilon']),
                                                       float(b['leading'][-1])])
    N=intervals['normalized'];L=intervals['linear_once_calibrated']
    c.yes('heldout.energy_intervals_disjoint',N['heldout_energy_interval'][1]<L['heldout_energy_interval'][0])
    c.yes('heldout.finite_source_response_intervals_disjoint',N['response_interval'][1]<L['response_interval'][0])
    instruments=[]
    for label in models:
        instruments.append(dict(menu=label,**derivative_and_instrument(c,label+'.source',bands[label,2.])))
    files=[HERE/'independent_check.py',HERE/'proof.md',STAGE/'_admission'/'round1326_admission.md',
           STAGE/'1325'/'proof.md',STAGE/'_admission'/'driven_reference_math_screen_after1325.md']
    return dict(schema='round1326_independent_v1',status='PASS_INDEPENDENT_CHECKS_ONLY',
                scientific_acceptance=False,author_code_imported=False,
                checks=len(c.rows),max_absolute_residual=max(x['residual'] for x in c.rows if 'residual' in x),
                finite_model=dict(B=B,dimension=len(norm['states']),p=p,material_dimension=4,
                                  normalized_coefficients=[norm['gd'],norm['gu']],
                                  linear_fixed_coefficients=[bare['gd'],bare['gu']]),
                small_full_parent=small,low_band_checks=rows,
                heldout_comparison=intervals,actual_source_instruments=instruments,evidence=c.rows,
                scope=dict(reference_cutoff_is_exact_model_not_continuum_limit=True,
                           old_resonant_and_driven_low_band_settings_distinct=True,
                           Dirac_free_kernel_is_inherited_input=True,
                           dressed_preparation_permission_is_input=True,
                           phase_and_source_control_not_free=True,
                           Higgs_identity_not_certified=True,flavor_copies_not_scanned=True,new_axioms=0),
                source_bindings={p.relative_to(STAGE).as_posix():digest(p) if p.exists() else None for p in files})


def compare(a,b,path='result'):
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys():raise AssertionError(path+': keys')
        for k in a:compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):raise AssertionError(path+': list')
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        if isinstance(b,bool) or not isinstance(b,(int,float)) or not np.isfinite(b) or not math.isclose(a,b,rel_tol=5e-10,abs_tol=5e-11):
            raise AssertionError((path,a,b))
    elif type(a) is not type(b) or a!=b:raise AssertionError((path,a,b))


def main():
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group()
    g.add_argument('--check',action='store_true');g.add_argument('--preview',action='store_true')
    g.add_argument('--write-new',action='store_true');a=p.parse_args()
    if a.write_new and RESULT.exists():raise FileExistsError(RESULT)
    if not a.write_new and not a.preview and not RESULT.is_file():raise FileNotFoundError(RESULT)
    r=compute()
    if a.write_new:
        if any(x is None for x in r['source_bindings'].values()):raise AssertionError('unfinished source binding')
        with RESULT.open('x',encoding='utf8',newline='\n') as f:
            json.dump(r,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    elif not a.preview:compare(r,json.loads(RESULT.read_text(encoding='utf8')))
    print(json.dumps(dict(status='PASS',mode='preview' if a.preview else 'write-new' if a.write_new else 'check',
                         checks=r['checks'],max_residual=r['max_absolute_residual'],
                         heldout=r['heldout_comparison']),indent=2))


if __name__=='__main__':main()
