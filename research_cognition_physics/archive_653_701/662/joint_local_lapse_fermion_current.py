"""662: original full CAR contribution to local lapse brackets.

Finite conditional coefficient checks, not a full Gauss thermal integral.
Common-domain and original-record results are proved in the numbered note.
The Fourier test uses one actual compact gauge direction and all 64 CAR modes;
its Nambu vectors test coefficients, not physical many-body states.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_full_graph_transfer_sources as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_local_lapse_fermion_current_results.json'
N=np.array([1.2,.7]); M=np.array([.6,1.3])


def local_blocks(data):
    d=data['d']; d1=np.zeros_like(d); d2=np.zeros_like(d)
    d1[:32,:32]=d[:32,:32]; d2[32:,32:]=d[32:,32:]
    return old.nambu(data['h']/2,d1),old.nambu(data['h']/2,d2)


def comm(a,b):return (a@b-b@a)/1j


def local_gauge_and_fock():
    data=old.graph_data();a,b=local_blocks(data);j=comm(a,b)
    assert old.norm(a+b-data['target_BdG'])<2e-13
    full=comm(N[0]*a+N[1]*b,M[0]*a+M[1]*b)
    coefficient=float(N[0]*M[1]-N[1]*M[0])
    identity=old.norm(full-coefficient*j)
    assert identity<2e-13 and old.norm(j)>.05
    rng=np.random.default_rng(662);rr=[];pp=[]
    for p in data['phi']:
        c=old.old.gauge.group_exp(.2*rng.normal(size=8),3)
        w=old.old.gauge.group_exp(.2*rng.normal(size=3),2)
        z=np.exp(.15j*rng.normal())
        rr.append(old.old.representation(c,w,z))
        x=z**3*w@(p[:2]+1j*p[2:4]);pp.append(np.r_[x.real,x.imag,p[4]])
    rh=old.blockdiag(*rr); rn=old.blockdiag(rh,rh.conj())
    transform=data['transform']; rt=transform@rn@transform.conj().T
    link=rr[0]@old.old.representation(*data['link'])@rr[1].conj().T
    h=np.zeros((64,64),complex);h[:32,32:]=.23*link;h[32:,:32]=.23*link.conj().T
    nh=transform@old.nambu(h,np.zeros_like(h))@transform.conj().T
    h0,d0=old.old.mass_matrices(pp[0]);h1,d1=old.old.mass_matrices(pp[1])
    bm=transform@old.nambu(old.blockdiag(h0,h1),old.blockdiag(d0,d1))@transform.conj().T
    at,bt=local_blocks(dict(h=nh[:64,:64],d=bm[:64,64:]))
    covariance=[old.norm(at-rt@a@rt.conj().T),old.norm(bt-rt@b@rt.conj().T),
                old.norm(comm(at,bt)-rt@j@rt.conj().T)]
    assert max(covariance)<2e-12
    # Independent many-body check only in the original invariant neutral sector.
    aa,bb=local_blocks(old.graph_data(neutral=True))
    modes=[30,31,62,63];ids=modes+[x+64 for x in modes]
    aa=aa[np.ix_(ids,ids)];bb=bb[np.ix_(ids,ids)]
    def lift(x):
        return old.car.fock(x[:4,:4],x[:4,4:])-.5*np.trace(x[:4,:4])*np.eye(16)
    af=lift(aa);bf=lift(bb);jf=comm(af,bf)
    fock_error=old.norm(jf-lift(comm(aa,bb)))
    rho=old.car.exp_h(af+bf,1.3);rho/=np.trace(rho)
    mean=np.trace(rho@jf);variance=np.trace(rho@jf@jf)-mean*mean
    assert fock_error<2e-13 and abs(mean)<2e-13 and variance.real>1e-5
    return dict(original_CAR_modes=64,local_matrix_current_norm=old.norm(j),
        positive_lapses=[N.tolist(),M.tolist()],lapse_coefficient=coefficient,
        lapse_identity_error=identity,independent_endpoint_gauge_errors=covariance,
        neutral_exact_Fock_identity_error=fock_error,
        neutral_stationary_mean=[float(mean.real),float(mean.imag)],
        neutral_stationary_variance=float(variance.real),
        conditional_matrix_test_not_full_Gauss_state=True)


def electric_mixed_bracket():
    data=old.graph_data();phi=data['phi']
    h0,d0=old.old.mass_matrices(phi[0]);h1,d1=old.old.mass_matrices(phi[1])
    zero=np.zeros((32,32),complex)
    b0=old.nambu(old.blockdiag(h0,zero),old.blockdiag(d0,zero))
    b1=old.nambu(old.blockdiag(zero,h1),old.blockdiag(zero,d1))
    charges=np.r_[np.ones(12),np.full(6,4),np.full(6,-2),np.full(4,-3),[-6,-6,0,0]]
    base=old.old.representation(*data['link']);grid=96
    t=np.arange(grid)*2*np.pi/grid;kappa=.37
    ks=[];dk=[];ddk=[]
    for x in t:
        link=np.exp(1j*x*charges)[:,None]*base
        pieces=[]
        for order in range(3):
            entry=.23*(1j*charges[:,None])**order*link
            hop=np.zeros((64,64),complex)
            hop[:32,32:]=entry;hop[32:,:32]=entry.conj().T
            pieces.append(old.nambu(hop,np.zeros_like(hop)))
        ks.append(pieces[0]);dk.append(pieces[1]);ddk.append(pieces[2])
    ks=np.array(ks);dk=np.array(dk);ddk=np.array(ddk)
    freq=np.fft.fftfreq(grid,d=1/grid)
    def deriv(v,order):
        return np.fft.ifft((1j*freq[:,None])**order*np.fft.fft(v,axis=0),axis=0)
    def act(b,v):return np.einsum('tij,tj->ti',b,v,optimize=True)
    def l2(v):return float(np.sqrt(np.mean(np.sum(abs(v)**2,axis=1))))
    rng=np.random.default_rng(6621)
    amplitudes=rng.normal(size=(5,128))+1j*rng.normal(size=(5,128))
    psi=np.exp(1j*t[:,None]*np.arange(-2,3)[None,:])@amplitudes
    psi/=l2(psi);dpsi=deriv(psi,1)
    nh=float(N.mean());mh=float(M.mean())
    bn=N[0]*b0+N[1]*b1+nh*ks
    bm=M[0]*b0+M[1]*b1+mh*ks
    car=act(bn,act(bm,psi))-act(bm,act(bn,psi));car/=1j
    tk_over_i=1j*kappa*(act(ddk,psi)+2*act(dk,dpsi))
    rows=[]
    for allocation in (np.array([1.,0.]),np.array([.5,.5])):
        nt=float(allocation@N);mt=float(allocation@M)
        def hn(v):return -nt*kappa*deriv(v,2)+act(bn,v)
        def hm(v):return -mt*kappa*deriv(v,2)+act(bm,v)
        actual=(hn(hm(psi))-hm(hn(psi)))/1j
        mixed=(nt*mh-mt*nh)*tk_over_i
        residual=l2(actual-mixed-car)
        assert residual<2e-10
        rows.append(dict(electric_endpoint_weights=allocation.tolist(),
            mixed_term_norm=l2(mixed),matrix_term_norm=l2(car),
            complete_bracket_residual=residual,
            omitted_mixed_term_error=l2(actual-car),
            omitted_matrix_term_error=l2(actual-mixed)))
    assert rows[0]['mixed_term_norm']>.05 and rows[0]['matrix_term_norm']>.01
    assert rows[1]['mixed_term_norm']<1e-14
    # Original hypercharge direction, not an invented scalar proxy for the link.
    source_check=max(old.norm(old.old.representation(np.eye(3),np.eye(2),np.exp(1j*x))-
                              np.diag(np.exp(1j*x*charges))) for x in (.13,.37))
    assert source_check<1e-13
    return dict(grid=grid,original_physical_hypercharges=charges.tolist(),
        wave_frequencies=[-2,-1,0,1,2],electric_coefficient=kappa,
        original_subgroup_identity_error=source_check,rows=rows,
        exact_finite_Fourier_polynomials_no_continuum_extrapolation=True,
        Nambu_coefficient_test_not_physical_Fock_wavefunction=True,
        other_gauge_directions_and_full_geometry_covered_analytically=True)


def run():
    deps=('research_note_585.md','research_note_587.md','research_note_588.md',
        'research_note_589.md','research_note_598.md','research_note_623.md','research_note_624.md',
        'joint_full_graph_transfer_sources.py','joint_fermion_gauss_completion.py',
        'joint_gauss_fermion_influence.py','round662_drafts/local_fermion_lapse_probe_results.json')
    return dict(date='2026-10-02',round=662,tests_run=2,failures=0,errors=0,
        local_gauge_and_Fock=local_gauge_and_fock(),electric_mixed_bracket=electric_mixed_bracket(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original fixed finite graph, declared positive local lapse allocation, full CAR mass and hopping. Common domain, Gauss and original scalar records analytically compatible; coefficient diagnostics retain original64 modes. No unique lapse allocation, full graph numerical Gauss integral, ADM closure, continuum chiral mapping or quantum GR.',
        all_checks_passed=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
