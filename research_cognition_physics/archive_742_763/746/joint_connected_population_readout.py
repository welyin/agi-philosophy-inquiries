"""746: original graph fourth-order CAR cancellations and certified native T readout.

Graph locality is an analytic commutator-word proof in the note. Finite
coefficient tests retain full original species and nontrivial gauge/spin links;
they are not bosonic graph time evolution. The interval certificate is exact.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round746_drafts'))
import higgs_readout_certificate as certificate
import joint_record_mass_feedback as inherited
matter=inherited.matter
TARGET=HERE/'joint_connected_population_readout_results.json'
def comm(a,b):return a@b-b@a
def nambu(h,d):return np.block([[h,d],[-d.conj(),-h.conj()]])

def coefficient_check():
    rng=np.random.default_rng(746);N=3;n=32*N;node=0
    P=np.zeros((n,n),complex);P[30,30]=P[31,31]=1
    ell=lambda K:np.trace(P@K[:n,:n])
    sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
    rows=[]
    for trial in range(6):
        x=rng.normal(size=(N,5))*.6
        h=np.zeros((n,n),complex);d=np.zeros_like(h)
        for v in range(N):
            hs,ds=inherited.mass_x(x[v]);sl=slice(32*v,32*v+32)
            h[sl,sl]=hs;d[sl,sl]=ds
        hopping=np.zeros_like(h)
        for v,w in ((0,1),(1,2),(2,0)):
            C=matter.gauge.group_exp(rng.normal(size=8),3)
            W=matter.gauge.group_exp(rng.normal(size=3),2)
            z=np.exp(1j*rng.normal());R=matter.representation(C,W,z)
            spin=sum(a*s for a,s in zip(rng.normal(size=3),(sx,sy,sz)))
            edge=(.3+.17j)*R@np.kron(np.eye(16),spin)
            sv=slice(32*v,32*v+32);sw=slice(32*w,32*w+32)
            hopping[sv,sw]=edge;hopping[sw,sv]=edge.conj().T
        B=nambu(h+hopping,d);C=nambu(hopping,np.zeros_like(hopping))
        xx=x[node].copy();xx[4]=0
        hd,dd=inherited.mass_x(xx);loc=np.zeros_like(h);loc[:32,:32]=hd
        D=nambu(loc,np.zeros_like(loc))
        # First force for A=f(T): scalar multiple of local Dirac block.
        third=abs(ell(comm(B,D)))
        cubic=abs(ell(comm(B,comm(B,D))))
        mixed=[]
        for i in range(5):
            hi,di=inherited.LINEAR_MASS[i]
            hh=np.zeros_like(h);dd=np.zeros_like(d);hh[:32,:32]=hi;dd[:32,:32]=di
            Qi=nambu(hh,dd)
            mixed.extend([abs(ell(Qi)),abs(ell(comm(C,Qi)))])
        assert max([third,cubic]+mixed)<2e-12
        # Native target identity grad_K T = (2*x_H/u,0) in x coordinates.
        xx=x[node];u=1+xx@xx/6;a=xx[:4]@xx[:4]
        grad=-a*xx/(3*u*u);grad[:4]+=2*xx[:4]/u
        force=grad+xx*(xx@grad)/6
        expected=np.r_[2*xx[:4]/u,0.]
        error=float(np.max(abs(force-expected)));assert error<2e-15
        rows.append(dict(full_modes=n,edges=3,nonzero_original_species_and_links=True,
            third_occupation_word=float(third),three_mass_or_hop_word=float(cubic),
            external_edge_two_factor_words=float(max(mixed)),target_force_identity_error=error))
    weights=np.array([.4,1.,2.5]);eps=.7;psi=(weights/eps**3)**(1/6)
    edge_direct=eps*((psi+np.roll(psi,-1))/2)**2
    edge_from_weights=((weights**(1/6)+np.roll(weights,-1)**(1/6))/2)**2
    relation=float(np.max(abs(edge_direct-edge_from_weights)));assert relation<2e-15
    # Bounded readout cost follows old instrument identity and this exact maximum.
    t=np.linspace(0,1,10001)
    gradient_squared=24*t*(1-t)**2
    assert max(gradient_squared)<=32/9+1e-14
    return dict(rows=rows,common_geometry_edge_volume_error=relation,
        exact_target_gradient_bound='K^-1(dT,dT) <= 32/9 for M=2',
        native_single_read_injection_bound='0 <= D_T <= hbar²/(9 w_v)',
        field_time_evolution_not_simulated=True)

def run():
    a=coefficient_check();b=certificate.run()
    assert b==json.loads((HERE/'round746_drafts/higgs_readout_certificate_results.json').read_text('utf8'))
    deps=('research_note_598.md','research_note_623.md','research_note_723.md','research_note_745.md',
          'joint_fermion_gauss_completion.py','round745_drafts/exact_history_density_results.json',
          'round745_drafts/certify_history_integral.py','round746_drafts/higgs_readout_certificate.py')
    return dict(round=746,tests_run=2,failures=0,errors=0,
        checks=['original_graph_CAR_word_and_common_geometry_checks','native_Higgs_readout_fourth_sign'],
        coefficient_checks=a,certificate=b,
        dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
        scope='In the original fixed finite-graph, fixed positive geometry model with species-preserving original gauge-covariant hopping, the existing T readout has a local fourth-order population contrast for a declared normal Gauss ready state. Analytic word cancellations remove arbitrary scalar graph potentials and hopping from that coefficient; a strict interval certifies its sign. No weakening of edges, graph boson time simulation, ideal occupation instrument, autonomous detector or dynamic gravity is claimed.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=746,tests_run=2,strict_negative=r['certificate']['strict_negative'])))
