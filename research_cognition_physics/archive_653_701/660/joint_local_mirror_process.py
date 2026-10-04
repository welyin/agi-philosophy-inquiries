"""660: actual onsite mirror completion, original weight and observation map.

Reflection positivity is an analytic free-overlap extension, not inferred from
the numerical reflection symmetry. No interacting-gauge or original-H claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_spatial_auxiliary_geometry as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_local_mirror_process_results.json'


def err(a): return float(np.max(np.abs(a)))
def norm(a): return float(np.linalg.norm(a))
def cpair(z): return [float(z.real),float(z.imag)]
def diag(a,b):
    return np.block([[a,np.zeros((len(a),len(b)),complex)],
                     [np.zeros((len(b),len(a)),complex),b]])


def data(kind='spatial',value=1.,e=None,rephase=False,nx=2,nt=2):
    if kind=='spatial':
        f=old.frame(value,nx,nt);sites=f['sites'];count=len(sites)
        u=np.kron(f['V'],np.eye(16))
        v=np.kron(f['vec'][:,f['ev']>0],np.eye(16))
        jm=np.kron(np.kron(np.eye(count),old.old.VM),np.eye(16))
        jp=np.kron(np.kron(np.eye(count),old.old.VP),np.eye(16))
        sign=(f['vec']*np.sign(f['ev']))@f['vec'].conj().T
        g=np.kron(np.eye(count),old.old.spin.G5)
        d=np.kron((np.eye(4*count)+g@sign)/2,np.eye(16))
        if e is None:
            e=np.random.default_rng(66011).normal(size=(count,10))*.08
            e[:,0]=1.;e/=np.linalg.norm(e,axis=1)[:,None]
        m=np.zeros_like(d);mb=np.zeros_like(d)
        for i,ei in enumerate(e):
            t=sum(a*b for a,b in zip(ei,old.old.T));sl=slice(64*i,64*i+64)
            m[sl,sl]=np.kron(old.old.B,t)
            mb[sl,sl]=np.kron(old.old.B,t.conj().T)
    else:
        u,v,d,_=old.old.frames(value);sites=[(0,0)]
        jm=np.kron(np.eye(16),old.old.VM);jp=np.kron(np.eye(16),old.old.VP)
        if e is None:
            e=np.random.default_rng(66012).normal(size=10);e/=np.linalg.norm(e)
        t=sum(a*b for a,b in zip(e,old.old.T))
        m=np.kron(t,old.old.B);mb=np.kron(t.conj().T,old.old.B)
    if rephase:
        u=u.copy();v=v.copy()
        u[:,0]*=np.exp(7j*value);v[:,0]*=np.exp(11j*value)
    qp=jp@jp.conj().T;qm=jm@jm.conj().T
    nambu=np.block([[m@qp,-d.T],[d,mb@qm]])
    s=diag(np.column_stack((u,v)),np.column_stack((jm.conj(),jp.conj())))
    nc=s.T@nambu@s;r=u.shape[1]
    a=u.T@m@u;bb=jm.conj().T@mb@jm.conj()
    k=jm.conj().T@d@u;kl=jp.conj().T@d@v
    c=u.T@m@qp@v;f=v.T@m@qp@v
    transform=np.eye(4*r,dtype=complex)
    transform[2*r:3*r,:r]=np.linalg.solve(bb,k)
    transform[3*r:,:r]=np.linalg.solve(kl.T,c.T)
    transform[3*r:,r:2*r]=-.5*np.linalg.solve(kl.T,f)
    expected=np.zeros_like(nc)
    expected[:r,:r]=a;expected[2*r:3*r,2*r:3*r]=bb
    expected[r:2*r,3*r:]=-kl.T;expected[3*r:,r:2*r]=kl
    return dict(N=nambu,S=s,Nc=nc,T=transform,expected=expected,A=a,bb=bb,
                K=k,Kl=kl,C=c,F=f,u=u,v=v,D=d,Qp=qp,Qm=qm,M=m,Mb=mb,
                e=e,sites=sites,r=r)


def exact_completion_check():
    rows=[]
    for kind,value in (('spatial',1.),('holonomy',.17),('holonomy',.31)):
        q=data(kind,value);ti=np.linalg.inv(q['T'])
        identities=max(err(q['D']@q['u']-q['Qm']@q['u']),
                       err(q['D']@q['v']-q['Qp']@q['v']),
                       err(q['K'].T@np.linalg.solve(q['bb'],q['K'])-
                           q['u'].T@q['M']@q['Qm']@q['u']),
                       err(ti.T@q['Nc']@ti-q['expected']))
        actual=old.old.pfaffian(q['N'])
        predicted=np.linalg.det(q['Kl'])*old.old.pfaffian(q['A'])/np.linalg.det(q['S'])
        relative=float(abs(actual-predicted)/abs(predicted))
        assert identities<3e-12 and relative<3e-11
        assert abs(old.old.pfaffian(q['bb'])-1)<2e-12
        assert abs(np.linalg.det(q['T'])-1)<1e-12
        rows.append(dict(branch=kind,parameter=value,Nambu_dimension=len(q['N']),
            exact_congruence_and_chiral_error=identities,full_signed_Pfaffian=cpair(actual),
            original_factorized_Pfaffian=cpair(predicted),relative_Pfaffian_error=relative))
    return dict(rows=rows,phase_not_replaced_by_absolute_value=True,
        pairing_convention='exp(+xi.T N xi / 2), kinetic coefficient 1, each pair coefficient 1/2',
        no_fifth_direction_limit_required=True)


def observation_check():
    q=data();r=q['r'];g=-np.linalg.inv(q['Nc']);eta=q['T']@g@q['T'].T
    target=-np.linalg.inv(q['expected'])
    full_error=err(eta-target)
    b_error=err(g[:r,:r]+np.linalg.inv(q['A']))
    bare_bar=norm(g[3*r:,3*r:]);dressed_bar=norm(eta[3*r:,3*r:])
    rng=np.random.default_rng(66013);j=rng.normal(size=(4*r,4))
    tj=np.linalg.solve(q['T'].T,j)
    source_error=err(j.T@g@j-tj.T@target@tj)
    assert max(full_error,b_error,source_error)<2e-10
    assert bare_bar>1 and dressed_bar<2e-11
    # The same identity also fixes four-point Wick insertions (not only Z).
    cov=j.T@g@j;tcov=tj.T@target@tj
    def wick4(a):return a[0,1]*a[2,3]-a[0,2]*a[1,3]+a[0,3]*a[1,2]
    four_error=float(abs(wick4(cov)-wick4(tcov)))
    assert four_error<2e-7
    return dict(all_transformed_Wick_kernel_error=full_error,b_kernel_error=b_error,
        bare_left_bar_pair_norm=bare_bar,dressed_left_bar_pair_norm=dressed_bar,
        transported_source_two_point_error=source_error,transported_source_four_point_error=four_error,
        bare_left_bar_observables_not_identified=True,
        physical_time_support_of_dressed_observables_not_established=True,
        original_independent_bar_auxiliary_observables_not_identified=True)


def reflection_check():
    rows=[]
    for nx,nt in ((2,2),(2,4)):
        q=data(nx=nx,nt=nt);sites=q['sites'];perm=[sites.index((nt-1-t,x)) for t,x in sites]
        rr=np.eye(len(sites))[perm]
        gamma=np.kron(np.kron(rr,old.old.spin.GAMMA[3]),np.eye(16))
        zeros=np.zeros_like(gamma);refl=np.block([[zeros,gamma.T],[gamma,zeros]])
        eq=data(nx=nx,nt=nt,e=q['e'][perm])
        reflection_error=err(q['N']+refl.T@eq['N'].conj()@refl)
        kinetic=np.block([[zeros,-q['D'].T],[q['D'],zeros]])
        interaction=q['N']-kinetic
        positive=np.repeat([t>=nt//2 for t,x in sites],64)
        positive=np.concatenate((positive,positive));negative=~positive
        cross_error=err(interaction[np.ix_(positive,negative)])
        plus=np.zeros_like(interaction);plus[np.ix_(positive,positive)]=interaction[np.ix_(positive,positive)]
        # E is reflected along with coefficients when applying theta to V_plus.
        vref=eq['N']-kinetic;pref=np.zeros_like(vref)
        pref[np.ix_(positive,positive)]=vref[np.ix_(positive,positive)]
        split_error=err(interaction-plus+refl.T@pref.conj()@refl)
        assert max(reflection_error,cross_error,split_error)<3e-12
        rows.append(dict(nx=nx,nt=nt,reflection_error=reflection_error,
            interaction_cross_half_error=cross_error,onsite_theta_split_error=split_error))
    return dict(rows=rows,
        numerical_symmetry_not_used_as_a_proof_of_reflection_positivity=True,
        analytic_RP_input='Free overlap link-reflection theorem, m0=1, unit gauge links, even AP time; onsite V_plus+theta(V_plus), positive product S9 measure.',
        general_gauge_RP_not_claimed=True,common_Hamiltonian_and_all_time_limits_not_claimed=True)


def source_check():
    theta=.23;step=2e-5;nodes,weights=np.polynomial.legendre.leggauss(6)
    radial=(nodes+1)/2;weights=weights/2*12*radial**2*(1-radial)
    actual=[];target=[];untransported=[];maxpoint=0.
    for angle in (theta-step,theta,theta+step):
        total=0j;fact=0j;wrong=0j
        for rho,w in zip(radial,weights):
            e=np.zeros(10);e[0]=np.sqrt(rho);e[6]=np.sqrt(1-rho)
            q=data('holonomy',angle,e=e,rephase=True)
            full=old.old.pfaffian(q['N'])
            numerator=np.linalg.det(q['Kl'])*old.old.pfaffian(q['A'])
            proper=numerator/np.linalg.det(q['S'])
            maxpoint=max(maxpoint,float(abs(full-proper)))
            total+=w*full;fact+=w*proper;wrong+=w*numerator
        actual.append(total);target.append(fact);untransported.append(wrong)
        assert abs(total-old.old.full(angle))<3e-12
    # Smooth complex logarithmic derivative, with no phase branch crossing.
    derivative=(np.log(actual[2]/actual[1])-np.log(actual[0]/actual[1]))/(2*step)
    wrong=(np.log(untransported[2]/untransported[1])-np.log(untransported[0]/untransported[1]))/(2*step)
    z,dz=old.old.measure(theta)
    predicted=-float(np.sum(old.old.Q*np.tan(old.old.Q*theta/2)))+dz/z
    assert maxpoint<3e-12 and abs(derivative-predicted)<2e-7
    assert abs((wrong-derivative)-18j)<2e-8
    return dict(theta=theta,full_original_S9_integral=cpair(actual[1]),
        exact_radial_quadrature_nodes=6,pointwise_factorization_error=maxpoint,
        full_source_log_derivative=cpair(derivative),inherited_analytic_source=predicted,
        source_error=float(abs(derivative-predicted)),
        discarded_frame_Jacobian_spurious_response=cpair(wrong-derivative),
        arbitrary_frame_phase_is_not_a_physical_response=True,
        gauge_holonomy_checked_algebraically_not_as_interacting_RP=True)


def run():
    deps=('joint_spatial_auxiliary_geometry.py','joint_subgroup_measure_source.py',
          'research_note_612.md','research_note_658.md','research_note_659.md',
          'round660_drafts/mirror_schur_probe_results.json','round660_drafts/source_receipt.json')
    return dict(date='2026-10-02',round=660,tests_run=4,failures=0,errors=0,
        exact_completion=exact_completion_check(),observation_map=observation_check(),
        free_reflection_interface=reflection_check(),common_source=source_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Exact finite original chiral weight and dressed surviving observations through onsite mirror completion; free finite-box RP uses the mature overlap theorem plus an explicit onsite extension. General gauge RP, original fixed-graph CAR/H, common infinite time process, continuum GR/SM are open.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=660,tests=result['tests_run'],passed=result['all_checks_passed'])))
