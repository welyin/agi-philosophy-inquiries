"""625: one finite CP compression for records and source variations.

The infinite full-Gauss convergence result is analytic in the note.
Numerics use an explicitly declared finite-difference diagnostic of the
original two-radial node quadratic form and its original singlet effects.
This is not a computation of the original full graph/CAR spectrum.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_source_preserving_compression_results.json'
BETA=.12
TIMES=(.12,.19)
HBAR=.7
VOLUME=.8


def diagnostic():
    n=8
    hs=np.linspace(.15,1.55,n+2)[1:-1]
    ss=np.linspace(-1.1,1.1,n+2)[1:-1]
    h,s=np.meshgrid(hs,ss,indexing='ij');h=h.ravel();s=s.ravel()
    f=original.M-(h*h+s*s)/6
    measure=np.sqrt(original.M)*h**3/f**3
    root_inv=1/np.sqrt(measure)
    def derivative(step):
        return (np.diag(np.ones(n-1),1)-np.diag(np.ones(n-1),-1))/(2*step)
    dh=np.kron(derivative(hs[1]-hs[0]),np.eye(n))*root_inv[None,:]
    ds=np.kron(np.eye(n),derivative(ss[1]-ss[0]))*root_inv[None,:]
    ghh=f*(1-h*h/(6*original.M));gss=f*(1-s*s/(6*original.M))
    ghs=-f*h*s/(6*original.M)
    T=HBAR**2/(2*VOLUME)*(dh.T@((measure*ghh)[:,None]*dh)
       +ds.T@((measure*gss)[:,None]*ds)+dh.T@((measure*ghs)[:,None]*ds)
       +ds.T@((measure*ghs)[:,None]*dh))
    phi=np.zeros((n*n,5));phi[:,1]=h;phi[:,4]=s
    W=np.diag(VOLUME*original.node_potential(phi))
    E,V=np.linalg.eigh(T+W)
    L=[V.T@(np.sqrt(.5+sign*np.sin(s)/4)[:,None]*V) for sign in (1,-1)]
    TT=V.T@T@V;WW=V.T@W@V
    probabilities=np.exp(-BETA*(E-E[0]));probabilities/=probabilities.sum()
    assert np.linalg.eigvalsh(T).min()>0
    assert np.max(abs(sum(x@x for x in L)-np.eye(n*n)))<1e-13
    return dict(T=TT,W=WW,E=E,L=L,prob=probabilities,dimension=n*n,
                box=dict(h=[.15,1.55],s=[-1.1,1.1],nodes_per_axis=n),
                minimum_target_F=float(f.min()),minimum_kinetic_eigenvalue=float(np.linalg.eigvalsh(T).min()))


def prepare(d,rank,rethermalize=False):
    q=d['prob'][:rank].copy()
    if rethermalize:q/=q.sum()
    else:q[0]+=d['prob'][rank:].sum()
    L=[x[:rank,:rank] for x in d['L']]
    defects=[(x@x)[:rank,:rank]-small@small for x,small in zip(d['L'],L)]
    return q,L,defects


def histories(d,rank,epsilon,plus=1.,minus=-.4,reset=True,contact=True,rethermalize=False):
    q,L,defects=prepare(d,rank,rethermalize)
    T=d['T'][:rank,:rank];W=d['W'][:rank,:rank]
    def unitary(u,duration):
        H=np.exp(-6*u)*T+np.exp(6*u)*W if contact else T+W+6*u*(W-T)
        e,v=np.linalg.eigh(H)
        return (v*np.exp(-1j*duration*e/HBAR))@v.conj().T
    output=[]
    us=[(unitary(epsilon*plus,t),unitary(epsilon*minus,t)) for t in TIMES]
    for history in itertools.product((0,1),repeat=len(TIMES)):
        x=np.diag(q).astype(complex)
        for r,(up,um) in zip(history,us):
            x=up@x@um.conj().T
            tail=np.trace(defects[r]@x)
            x=L[r]@x@L[r]
            if reset:x[0,0]+=tail
        output.append(np.trace(x))
    return np.array(output)


def coefficients(d,rank,step=.001,**options):
    z0=histories(d,rank,0.,**options)
    zp=histories(d,rank,step,**options)
    zm=histories(d,rank,-step,**options)
    lp=np.log(zp.sum());lm=np.log(zm.sum());l0=np.log(z0.sum())
    pplus=histories(d,rank,step,plus=1.,minus=1.,**options).real
    pminus=histories(d,rank,-step,plus=1.,minus=1.,**options).real
    return dict(probability=z0.real,score=(pplus-pminus)/(2*step*z0.real),
                log_linear=(lp-lm)/(2*step),log_quadratic=(lp+lm-2*l0)/(2*step**2),
                branch_quadratic=(zp+zm-2*z0)/(2*step**2),normalization=z0.sum())


def scalar_data(c):
    return dict(probabilities=c['probability'].tolist(),scores=c['score'].tolist(),
                log_linear=[float(c['log_linear'].real),float(c['log_linear'].imag)],
                log_quadratic=[float(c['log_quadratic'].real),float(c['log_quadratic'].imag)],
                normalization_error=float(abs(c['normalization']-1)))


def cp_check(d):
    rows=[]
    for rank in (6,12,24,40,64):
        q,L,defects=prepare(d,rank)
        cp=min(float(np.linalg.eigvalsh((x+x.T)/2).min()) for x in defects)
        tp=float(np.max(abs(sum(x@x+e for x,e in zip(L,defects))-np.eye(rank))))
        full_moment=float(d['prob']@((d['E']-d['E'][0]+1)**2))
        cut_moment=float(q@((d['E'][:rank]-d['E'][0]+1)**2))
        norm=histories(d,rank,.003,plus=1.,minus=1.).sum()
        assert cp>-2e-13 and tp<2e-13 and abs(norm-1)<2e-12
        assert cut_moment<=full_moment+2e-12
        rows.append(dict(rank=rank,defect_minimum_eigenvalue=cp,trace_preservation_error=tp,
                         total_history_normalization_error=float(abs(norm-1)),
                         shifted_A2=cut_moment,full_shifted_A2=full_moment,
                         initial_removed_probability=float(d['prob'][rank:].sum())))
    return dict(rows=rows,reset_keeps_passive_reference_in_the_analytic_construction=True)


def joint_coefficients_check(d):
    full=coefficients(d,d['dimension']);fine=coefficients(d,d['dimension'],step=.0005)
    rows=[]
    for rank in (6,12,24,40,64):
        c=coefficients(d,rank)
        row=scalar_data(c);row.update(rank=rank,
            probability_error=float(np.max(abs(c['probability']-full['probability']))),
            score_error=float(np.max(abs(c['score']-full['score']))),
            quadratic_error=float(abs(c['log_quadratic']-full['log_quadratic'])),
            branch_quadratic_error=float(np.max(abs(c['branch_quadratic']-full['branch_quadratic']))))
        assert row['normalization_error']<2e-12
        rows.append(row)
    assert rows[-1]['quadratic_error']<1e-12 and rows[-1]['score_error']<1e-12
    assert full['log_quadratic'].real<0
    return dict(rows=rows,full_finite_diagnostic=scalar_data(full),
                half_step_quadratic_difference=float(abs(fine['log_quadratic']-full['log_quadratic'])),
                finite_difference_coefficients_not_interval_certificates=True,
                fixed_64_dimensional_check_not_full_Gauss_cutoff_rate=True)


def wrong_connections_check(d):
    rank=12
    proper=coefficients(d,rank)
    rethermalized=coefficients(d,rank,rethermalize=True)
    omitted_contact=coefficients(d,rank,contact=False)
    noreset=histories(d,rank,0.,reset=False).sum()
    gaps=dict(rethermalized_linear_source_difference=float(abs(proper['log_linear']-rethermalized['log_linear'])),
              omitted_contact_quadratic_difference=float(abs(proper['log_quadratic']-omitted_contact['log_quadratic'])),
              bare_projected_instrument_probability_loss=float((1-noreset).real))
    assert gaps['rethermalized_linear_source_difference']>1e-5
    assert gaps['omitted_contact_quadratic_difference']>1e-5
    assert gaps['bare_projected_instrument_probability_loss']>1e-7
    return dict(rank=rank,**gaps,original_state_and_Hessian_not_refitted=True)


def run():
    d=diagnostic()
    names=('joint_curved_quantum_source.py','research_note_592.md','research_note_612.md',
           'research_note_622.md','research_note_623.md','research_note_624.md')
    return dict(round=625,tests_run=3,failures=0,errors=0,
        diagnostic_definition={k:d[k] for k in ('dimension','box','minimum_target_F','minimum_kinetic_eigenvalue')},
        beta=BETA,times=list(TIMES),hbar=HBAR,volume=VOLUME,
        cp_instrument=cp_check(d),joint_source_coefficients=joint_coefficients_check(d),
        incompatible_connections=wrong_connections_check(d),
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names},
        scope='Full original finite-graph theorem: one finite CP spectral approximation converges jointly in finite-history source coefficients through order two under common D and finite A2. Numerics only check an explicitly declared 64-dimensional radial finite-difference diagnostic, not full graph/CAR spectra, continuum, local EFT, or a physical free reset.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))
