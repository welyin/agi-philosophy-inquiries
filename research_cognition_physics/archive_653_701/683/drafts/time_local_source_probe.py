"""683 executed entry: strict physical-time source versus678 nearest-time source.
Only an auxiliary-regulator source replacement; no positive bulk construction.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_rational_physical_limit as prior
base=prior.base
TARGET=HERE/'time_local_source_probe_results.json'
def err(a):return float(np.max(abs(a),initial=0))

def run():
    links,e,phis=prior.fixture();mat=base.fixed_matrices(e,phis)
    jm,jp,M,Mbar,pair=mat;n=len(M);r=n//2
    g5=jp@jp.T-jm@jm.T
    u,v,_,h,gap=base.kernel(links);x=g5@h;xnorm=float(np.linalg.norm(x,2))
    assert xnorm<=3+1e-12
    T=np.vstack((.5*jm.T,.25*jp.T@M))
    # Endpoint-only column order psi,chi,z0,zL; intermediate z layers do not
    # appear in either observation. This avoids falsely claiming a large-body run.
    strict=np.zeros((n,4*n),complex)
    strict[:r,:n]=.5*jm.T
    strict[r:,:n]=-.75*jp.T@M
    strict[r:,n:2*n]=jp.T
    strict[:,2*n:3*n]=2*T;strict[:,3*n:]=-2*T
    positive_rows=np.r_[np.arange(64,128),np.arange(192,256)]
    negative_columns=np.concatenate([np.arange(4*64//2)+n*i for i in range(4)])
    assert err(strict[np.ix_(positive_rows,negative_columns)])==0
    exact=prior.soft(v@v.conj().T-u@u.conj().T,mat,lam=.37)
    source=np.eye(n)[:,[66,93]]
    rows=[]
    golden=(1+np.sqrt(5))/2
    for a in (.4,.2,.1,.05):
        layers=int(np.ceil(4/a**2))
        eps,_,_=prior.regulate(h,g5,a,layers)
        b=2*np.eye(n)+a*x
        f0=np.linalg.solve(b,(np.eye(n)+eps)/2)
        fl=np.linalg.solve(b,(np.eye(n)-eps)/2)
        pull=np.zeros((4*n,2*n),complex)
        pull[:2*n]=np.eye(2*n);pull[2*n:3*n,:n]=f0;pull[3*n:,:n]=fl
        modified=strict@pull
        original=prior.soft(eps,mat,lam=.37)
        nearest=strict.copy()
        nearest[:,2*n:3*n]=T@b;nearest[:,3*n:]=-T@b
        assert err(nearest@pull-original['Phi'])<2e-13
        formula=np.zeros_like(modified)
        formula[:,:n]=-a*T@np.linalg.solve(b,x@eps)
        assert err(modified-original['Phi']-formula)<2e-13
        actual=float(np.linalg.norm(modified-original['Phi'],2))
        bound=a*xnorm/(2*(2-a))
        assert actual<=bound+2e-13
        assert np.linalg.norm(np.vstack((f0,fl)),2)<=1/(2-a)+2e-13
        n0=prior.soft(eps,mat,lam=0)['N']
        modifiedN=n0+.37*modified.T@pair@modified
        mass_error=float(np.linalg.norm(modifiedN-original['N'],2))
        mass_bound=.37*np.linalg.norm(pair,2)*bound*(2*golden+bound)
        assert mass_error<=mass_bound+2e-13
        coef=[]
        for count in (0,2):
            z=source[:,:count]
            def coeff(N,Phi):
                c=Phi.T@z
                return base.pf(np.block([[N,c],[-c.T,np.zeros((count,count))]]))
            old=coeff(original['N'],original['Phi'])
            new=coeff(modifiedN,modified)
            target=coeff(exact['N'],exact['Phi'])
            coef.append(dict(sources=count,finite_original=base.old.cpair(old),
                strict_time=base.old.cpair(new),sign_limit=base.old.cpair(target),
                absolute_change=float(abs(new-old))))
        rows.append(dict(a=a,layers=layers,a_times_layers=a*layers,
            original_cross_time_source_norm=float(np.linalg.norm(nearest[np.ix_(positive_rows,negative_columns)],2)),
            strict_cross_time_source_norm=0.,effective_source_change=actual,
            uniform_in_layers_source_bound=float(bound),
            complete_mass_change=mass_error,complete_mass_bound=float(mass_bound),
            original_source_to_sign_error=float(np.linalg.norm(original['Phi']-exact['Phi'],2)),
            strict_source_to_sign_error=float(np.linalg.norm(modified-exact['Phi'],2)),
            physical_coefficients=coef))
    return dict(date='2026-10-02',entry_round=683,not_formal_round=True,
        original_nonflat_full_16_channels=True,Wilson_gap=gap,Wilson_X_norm=xnorm,
        uniform_source_bound_derived_not_fitted=True,
        regulator_a_not_physical_time_spacing=True,
        original_large_auxiliary_matrices_not_constructed=True,rows=rows,
        strict_time_observable_in_auxiliary_variables_only=True,
        full_bulk_reflection_positivity_still_unproved=True,
        complete_Gauss_average_bound_is_analytic_not_numerically_integrated=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in ('research_note_669.md','research_note_673.md','research_note_677.md',
                      'research_note_678.md','research_note_682.md','joint_rational_physical_limit.py')})

if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=683,all_checks_passed=True,
        source_errors=[r['effective_source_change'] for r in result['rows']],
        upper_bounds=[r['uniform_in_layers_source_bound'] for r in result['rows']])))

