"""723: full original reference histories, Gauss witness, energy-weighted transport.

The infinite-dimensional and full-Gauss conclusions are analytic. The first
group uses the actual one-node flow-coordinate resolvent on radial Gauss
wavefunctions; the other groups use the inherited64D radial diagnostic.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_squared_reference_readout as old
from round723_drafts import joint_reference_entry as entry
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_reference_process_transport_results.json'


def trace_norm(a):
    return float(np.sum(abs(np.linalg.eigvalsh((a+a.conj().T)/2))))


def radial_gauss(n):
    # Two real Hermite profiles have Gaussian tails; both are in original E.
    L=10.;xi=np.linspace(-L,L,n,endpoint=False);dx=2*L/n
    sigma=.65;center=-.5;q=(xi-center)/sigma
    a=np.exp(-q*q/2)/(np.pi**.25*np.sqrt(sigma));b=np.sqrt(2)*q*a
    aprime=-q*a/sigma;bprime=np.sqrt(2)*(1-q*q)*a/sigma
    z,wz=old.gl(48,-.5,.5);Z=z[:,None]
    bz=(1-4*z*z)**2;norm=np.sqrt(np.dot(wz,bz*bz));bz/=norm
    bzprime=-16*z*(1-4*z*z)/norm
    lo=np.full((len(z),n),-40.);hi=np.full_like(lo,5.)
    for _ in range(90):
        m=(lo+hi)/2
        v=((1+Z*Z/6)*m+np.exp(2*m)/12)/old.M
        lo=np.where(v<xi,m,lo);hi=np.where(v>=xi,m,hi)
    R=np.exp((lo+hi)/2);D=1+(R*R+Z*Z)/6
    T=old.M*R*R/(2*D);field_effect=.5+.25*np.sin(T)
    weights=wz[:,None]*dx
    lam=64.;shift=1.;hbar=old.old.HBAR
    momentum=2*np.pi*np.fft.fftfreq(n,d=dx)
    multiplier=lam/(lam+1j*(shift+hbar*hbar*momentum**2))
    ka=np.fft.ifft(multiplier*np.fft.fft(a))
    kb=np.fft.ifft(multiplier*np.fft.fft(b))
    ka_prime=np.fft.ifft(1j*momentum*multiplier*np.fft.fft(a))
    kb_prime=np.fft.ifft(1j*momentum*multiplier*np.fft.fft(b))
    phi=np.zeros(R.shape+(5,));phi[...,0]=R*np.sqrt(old.M/D)
    phi[...,4]=Z*np.sqrt(old.M/D)
    potential=old.old.original.node_potential(phi)
    def energy(profile,derivative):
        U=bz[:,None]*profile
        Uxi=bz[:,None]*derivative;Uz=bzprime[:,None]*profile
        AR=Uxi*((1+Z*Z/6)/R+R/6)/old.M-.5*U*(4/R-R/(2*D))
        Az=Uz+Uxi*Z*np.log(R)/(3*old.M)+.25*U*Z/D
        kinetic=abs(AR)**2*(1+R*R/6)+abs(Az)**2*(1+Z*Z/6)
        kinetic+=np.real(AR.conj()*Az)*R*Z/3
        return float(np.sum(weights*(hbar*hbar/2*kinetic+potential*abs(U)**2)))
    probs=[];initial=[];outputs=[];all_probs=[]
    for sign in (1,-1):
        psi=(a+sign*1j*b)/np.sqrt(2)
        dpsi=(aprime+sign*1j*bprime)/np.sqrt(2)
        kpsi=(ka+sign*1j*kb)/np.sqrt(2)
        kprime=(ka_prime+sign*1j*kb_prime)/np.sqrt(2)
        probs.append(float(np.sum(weights*bz[:,None]**2*abs(kpsi)**2*field_effect)))
        initial.append(energy(psi,dpsi))
        outputs.append(energy(kpsi,kprime)+energy(psi-kpsi,dpsi-kprime))
        record=[]
        for wave in (kpsi,psi-kpsi):
            for eff in (field_effect,1-field_effect):
                record.append(float(np.sum(weights*bz[:,None]**2*abs(wave)**2*eff)))
        all_probs.append(record);assert abs(sum(record)-1)<2e-12
    kappa=np.sqrt(shift-1j*lam)/hbar
    assert kappa.real>old.M and abs(initial[0]-initial[1])<1e-12
    assert abs(probs[0]-probs[1])>1e-4
    return dict(grid=n,actual_history_probabilities=probs,difference=probs[0]-probs[1],
                all_four_record_probabilities=all_probs,original_initial_energy=initial,
                post_Q_total_original_energy=outputs,exponential_kernel_decay=float(kappa.real),
                required_negative_tail_rate=old.M,
                norm_error=float(abs(dx*np.dot(a,a)-1)),
                real_pair_overlap=float(dx*np.dot(a,b)))


def gauss_witness_check():
    assert entry.run()==json.loads(entry.TARGET.read_text('utf8'))
    rows=[radial_gauss(n) for n in (1024,2048)]
    differences=[abs(x-y) for x,y in zip(rows[0]['actual_history_probabilities'],rows[1]['actual_history_probabilities'])]
    assert max(differences)<1e-11
    # Full-law equality is proved by conjugation, not a finite moment test.
    return dict(rows=rows,record_grid_difference=max(differences),
                full_original_H_and32_CAR_retained=True,
                spectral_law_equality_and_energy_domain_proved_analytically=True,
                FFT_quadrature_not_a_full_graph_spectrum=True)


def matrix_data(g,basis):
    t=old.fixture(g,'T');s=old.fixture(g,'s')
    transform=lambda a:basis.conj().T@a@basis
    menu=[];derivatives=[]
    for d in (t,s):
        ev,u=np.linalg.eigh(d['f'])
        field=[(u*np.sqrt(.5+sign*.25*np.sin(ev)))@u.conj().T for sign in (1,-1)]
        menu.append([transform(k) for k in field])
        derivatives.append([np.zeros_like(menu[-1][0]) for _ in field])
        menu.append([transform(k) for k in d['Ks']])
        derivatives.append([transform(k) for k in d['Kprimes']])
    refs=[transform(t['f']),transform(s['f']),transform(t['b']*t['I']-t['R']),transform(s['b']*s['I']-s['R'])]
    return dict(H=transform(t['H']),G=transform(t['G']),menu=menu,derivatives=derivatives,refs=refs)


def gibbs_with_derivative(h,g):
    beta=old.old.old625.BETA;e,u=np.linalg.eigh(h)
    values=np.exp(-beta*(e-e[0]));Z=sum(values)
    diff=e[:,None]-e[None,:]
    numerator=values[:,None]-values[None,:]
    divided=np.empty_like(diff)
    mask=abs(diff)>1e-9
    np.divide(numerator,diff,out=divided,where=mask)
    divided[~mask]=(-beta*(values[:,None]+values[None,:])/2)[~mask]
    unprime=u@(divided*(u.conj().T@g@u))@u.conj().T
    rho=(u*(values/Z))@u.conj().T
    prime=unprime/Z-rho*np.trace(unprime).real/Z
    return rho,prime


def process(data,rank,waiting,derivative=False,omit_tail=False):
    h=data['H'][:rank,:rank];g=data['G'][:rank,:rank]
    rho,rhop=gibbs_with_derivative(h,g)
    initial=rho.copy();initial_prime=rhop.copy()
    leaves=[(rho,rhop)]
    e,u=np.linalg.eigh(h);wait=(u*np.exp(-.19j*e/old.old.HBAR))@u.conj().T
    reset=np.zeros_like(rho);reset[0,0]=1
    for index,(menu,derivs) in enumerate(zip(data['menu'],data['derivatives'])):
        new=[]
        for y,yp in leaves:
            if waiting and index:
                y=wait@y@wait.conj().T
                yp=wait@yp@wait.conj().T  # derivatives not claimed with moving waits
            for K,Kp in zip(menu,derivs):
                B=K[:rank,:rank];Bp=Kp[:rank,:rank]
                tail=K[rank:,:rank];tailp=Kp[rank:,:rank]
                D=tail.conj().T@tail;Dp=tailp.conj().T@tail+tail.conj().T@tailp
                if omit_tail:D=np.zeros_like(D);Dp=np.zeros_like(Dp)
                out=B@y@B.conj().T+np.trace(D@y)*reset
                outp=(Bp@y@B.conj().T+B@yp@B.conj().T+B@y@Bp.conj().T
                      +(np.trace(Dp@y)+np.trace(D@yp))*reset)
                new.append((out,outp))
        leaves=new
    def embed(a):
        out=np.zeros_like(data['H'],dtype=complex);out[:rank,:rank]=a;return out
    outputs=[embed(y) for y,_ in leaves];primes=[embed(yp) for _,yp in leaves]
    assert omit_tail or abs(sum(np.trace(y).real for y in outputs)-1)<2e-12
    if derivative:assert abs(sum(np.trace(y).real for y in primes))<1e-10
    return outputs,primes,embed(initial),embed(initial_prime)


def weighted_error(a,b,root):
    return sum(trace_norm(root@(x-y)@root) for x,y in zip(a,b))


def matrix_fixture():
    g=.017;t=old.fixture(g,'T');e,basis=np.linalg.eigh(t['H'])
    data=matrix_data(g,basis);A=np.diag(e+1-e[0]);root=np.diag(np.sqrt(np.diag(A)))
    ranks=[]
    for rank in (16,32,48,64):
        while rank<64 and abs(e[rank]-e[rank-1])<1e-10:rank+=1
        if rank not in ranks:ranks.append(rank)
    return g,e,basis,data,A,root,ranks


def process_transport_check():
    g,e,basis,data,A,root,ranks=matrix_fixture()
    full=process(data,64,True)[0];rows=[]
    inv=np.diag(1/np.diag(root));budget_factor=1.
    for menu in data['menu']:
        average=sum(k.conj().T@A@k for k in menu)
        budget_factor*=float(np.linalg.eigvalsh(inv@average@inv).max())
    for rank in ranks:
        outputs,_,initial,_=process(data,rank,True)
        probs=[float(np.trace(y).real) for y in outputs]
        energy=float(np.trace(A@sum(outputs)).real)
        budget=budget_factor*float(np.trace(A@initial).real)
        assert energy<=budget+1e-9 and min(probs)>0
        error=weighted_error(outputs,full,root)
        reference_error=max(abs(np.trace(f@(sum(outputs)-sum(full)))) for f in data['refs'])
        source_error=float(abs(np.trace(data['G']@(sum(outputs)-sum(full)))))
        wrong=process(data,rank,True,omit_tail=True)[0]
        loss=1-sum(np.trace(y).real for y in wrong)
        rows.append(dict(rank=rank,energy_weighted_full_record_error=error,
                         four_reference_mean_error=float(reference_error),terminal_source_error=source_error,
                         original_shifted_energy=energy,finite_matrix_budget=budget,
                         probability_loss_if_tail_dropped=float(loss),all16_record_probabilities=probs))
    assert rows[-1]['energy_weighted_full_record_error']<1e-11
    assert rows[0]['probability_loss_if_tail_dropped']>1e-3
    assert rows[-2]['energy_weighted_full_record_error']<rows[0]['energy_weighted_full_record_error']
    return dict(rows=rows,same_original_waiting=True,complete_CP_tail_kept=True,
                fixed_graph_finite_diagnostic_not_spatial_refinement=True)


def source_transport_check():
    g,e,basis,data,A,root,ranks=matrix_fixture();step=2e-6
    plus=matrix_data(g+step,basis);minus=matrix_data(g-step,basis)
    full,full_prime,_,_=process(data,64,False,True);rows=[]
    for rank in ranks:
        outputs,primes,initial,initial_prime=process(data,rank,False,True)
        yp,_,ip,_=process(plus,rank,False);ym,_,im,_=process(minus,rank,False)
        fd=[(p-m)/(2*step) for p,m in zip(yp,ym)]
        fd_error=weighted_error(fd,primes,root)
        derivative_error=weighted_error(primes,full_prime,root)
        change=sum(outputs)-initial;dchange=sum(primes)-initial_prime
        analytic=float(np.trace(data['G']@change+data['H']@dchange).real)
        Ep=np.trace(plus['H']@(sum(yp)-ip)).real
        Em=np.trace(minus['H']@(sum(ym)-im)).real
        energy_fd=float((Ep-Em)/(2*step))
        assert fd_error<2e-5 and abs(energy_fd-analytic)<3e-6
        rows.append(dict(rank=rank,energy_weighted_derivative_error=derivative_error,
                         full_record_derivative_vs_fd_error=fd_error,
                         joint_energy_change_derivative=analytic,independent_energy_derivative=energy_fd))
    assert rows[-1]['energy_weighted_derivative_error']<1e-11
    assert rows[-2]['energy_weighted_derivative_error']<rows[0]['energy_weighted_derivative_error']
    return dict(rows=rows,own_compressed_Gibbs_preparation=True,
                fixed_original_projection=True,immediate_four_reference_history=True,
                no_moving_wait_dynamic_second_derivative_claim=True)


def run():
    results=dict(original_Gauss_witness=gauss_witness_check(),
                 common_reference_process=process_transport_check(),
                 immediate_joint_source=source_transport_check())
    names=('joint_squared_reference_readout.py','research_note_722.md','research_note_652.md',
           'research_note_625.md','research_note_704.md','research_note_707.md','research_note_708.md',
           'round723_drafts/joint_reference_entry.py','round723_drafts/joint_reference_entry_results.json')
    return dict(round=723,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in names},
                scope='Full original fixed-graph Gauss model: equal four single-reference spectral laws need not fix ordered records. New full-reference instruments admit common energy-weighted finite CP histories with fixed original waiting and terminal form sources; immediate joint records and own Gibbs preparation admit first geometry derivative transport. No autonomous apparatus, moving-wait second derivative, uniform spatial limit, relational coordinate inverse or generated gravity.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
