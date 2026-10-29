"""Round 1318 independent field enumeration, complete spectra and common sources.

No author code is imported or executed. Default/--check is read-only;
--write-new exclusively creates independent_results.json.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'independent_results.json'
ADMISSION = HERE.parent / '_admission' / 'round1318_admission.md'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def norm(a):
    return float(np.linalg.norm(np.atleast_2d(a), 2))


class Audit:
    def __init__(self):
        self.errors = {}
        self.count = 0

    def truth(self, label, test):
        self.count += 1
        self.errors.setdefault(label, 0.)
        if not bool(test):
            raise AssertionError(label)

    def close(self, label, a, b, tol=2e-9):
        err = norm(np.asarray(a)-np.asarray(b))/max(1., norm(b))
        self.count += 1
        self.errors[label] = max(self.errors.get(label, 0.), err)
        if err > tol:
            raise AssertionError((label, err, tol))


def fields(n, a):
    # count, colour fundamental sign (0 singlet), SU(3)L sign, x=3X, r=i^s
    b = n-a
    return {
        'L': (n, 0, 1, -1, 1), 'Qp': (a, 1, 1, 1, 1),
        'Qm': (b, 1, -1, 0, 1), 'uc': (n, -1, 0, -2, 1),
        'dc': (n, -1, 0, 1, 1), 'ec': (n, 0, 0, 3, 1),
        'Nc': (n, 0, 0, 0, 1), 'Fc': (n, 0, 0, 0, -1),
        'Uc': (a, -1, 0, -2, -1), 'Dc': (b, -1, 0, 1, -1)}


def anomaly_sums(n, a):
    ans = dict(C3=0, L3=0, C2x=0, L2x=0, x3=0, gravx=0,
               rC=0, rL=0, rgrav=0, rx2=0, r2x=0)
    for count, c, l, x, s in fields(n, a).values():
        dc, dl = (3 if c else 1), (3 if l else 1)
        ans['C3'] += count*c*dl
        ans['L3'] += count*l*dc
        ans['C2x'] += count*bool(c)*dl*x
        ans['L2x'] += count*bool(l)*dc*x
        ans['x3'] += count*dc*dl*x**3
        ans['gravx'] += count*dc*dl*x
        ans['rC'] += count*bool(c)*dl*s
        ans['rL'] += count*bool(l)*dc*s
        ans['rgrav'] += count*dc*dl*s
        ans['rx2'] += count*dc*dl*s*x*x
        ans['r2x'] += count*dc*dl*s*s*x
    return ans


def field_tests(audit):
    cancelled, selected = [], []
    for n in range(1, 25):
        for a in range(n+1):
            an = anomaly_sums(n, a)
            expected = dict(C3=0, L3=2*(3*a-n), C2x=0,
                            L2x=3*a-n, x3=6*(n-3*a), gravx=0,
                            rC=4*n, rL=4*n, rgrav=16*n,
                            rx2=24*n, r2x=0)
            audit.truth('exact_field_anomaly_counts', an == expected)
            clean = all(an[k] == 0 for k in ['C3','L3','C2x','L2x','x3','gravx'])
            audit.truth('integer_classification_equivalence', clean == (n == 3*a))
            if clean:
                cancelled.append([n,a,n-a])
                if 11-2*n > 0:
                    selected.append([n,a,n-a])
    audit.truth('positive_full_threshold_selection', selected == [[3,1,2]])
    audit.truth('n6_survives_without_positive_sign', anomaly_sums(6,2)['L3'] == 0 and 11-12 < 0)
    audit.truth('n6_low_threshold_positive', 11-4*6/3 == 3)
    # Leading RG solution: a negative b0 still permits an ordinary finite window.
    g0, L, b0 = .1, 10., -1.
    gend = 1/np.sqrt(g0**-2+b0*L/(8*np.pi**2))
    audit.truth('finite_window_does_not_force_asymptotic_freedom', g0 < gend < .101)
    rows = fields(3,1)
    scalar = {'eta': (1,-1,2), 'rho': (1,2,2), 'chi': (1,-1,0)}
    terms, bare = [], []
    # Enumerate all two-Weyl one-scalar invariants for this field content,
    # using the 3*3*3 epsilon and 3*bar3 delta invariant tensors only.
    for xname, yname in itertools.combinations_with_replacement(rows, 2):
        _,c1,l1,x1,s1 = rows[xname]
        _,c2,l2,x2,s2 = rows[yname]
        if not ((c1 == c2 == 0) or (c1 == -c2 and c1 != 0)):
            continue
        if l1 == l2 == 0 and x1+x2 == 0 and (s1+s2)%4 == 0:
            bare.append(xname+' '+yname)
        for sn,(ls,xs,ss) in scalar.items():
            for star in [False,True]:
                t = -1 if star else 1
                tri = sorted([l1,l2,t*ls])
                inv = tri in [[-1,0,1],[1,1,1],[-1,-1,-1]]
                if inv and x1+x2+t*xs == 0 and (s1+s2+ss)%4 == 0:
                    terms.append(xname+' '+yname+' '+sn+('*' if star else ''))
    expected = ['L L rho','L ec rho*','L Nc eta*','L Fc chi*',
                'Qp uc eta*','Qp dc rho*','Qp Uc chi*',
                'Qm uc rho','Qm dc eta','Qm Dc chi']
    audit.truth('all_allowed_renormalizable_yukawas', sorted(terms) == sorted(expected))
    audit.truth('only_bare_neutral_cross', bare == ['Nc Fc'])
    audit.truth('same_kind_majorana_forbidden_by_r', (1+1)%4 != 0 and (-1-1)%4 != 0)
    audit.truth('same_kind_majorana_allowed_without_r', rows['Nc'][3] == rows['Fc'][3] == 0)
    # r_eff on vacuum and all-left neutral variables; a symmetric Weyl mass
    # can connect the two opposite phase spaces only.
    rg = np.diag([-1.,-1.,1.])
    for index, sn in enumerate(['eta','rho','chi']):
        audit.close('vacuum_preserves_residual', (1j**scalar[sn][2])*rg[:,index], np.eye(3)[:,index])
    rn = np.array([-1j]*3+[-1j]*3+[1j]*3+[1j]*3)
    allowed = np.isclose(rn[:,None]*rn[None,:],1.)
    audit.truth('residual_all_quadratic_blocks', np.all(~allowed[:6,:6]) and np.all(~allowed[6:,6:]) and np.all(allowed[:6,6:]))
    # Full all-fields common centre kernel, not a global anomaly certificate.
    kernel=[]
    for zc,zl,zx in itertools.product(range(3),repeat=3):
        ok=all((c*zc+l*zl+x*zx)%3 == 0 for _,c,l,x,_ in rows.values())
        ok=ok and all((l*zl+x*zx)%3 == 0 for l,x,_ in scalar.values())
        if ok: kernel.append([zc,zl,zx])
    audit.truth('faithful_centre_kernel', kernel == [[0,0,0],[1,1,1],[2,2,2]])
    alpha=np.pi/3
    lowembed=np.diag([-1.,-1.,1.])*np.exp(1j*np.array([-alpha,-alpha,2*alpha]))
    audit.close('low_G6_centre_embedding',lowembed,np.exp(2j*np.pi/3)*np.eye(3))
    return dict(integer_scan_n_max=24,local_anomaly_free=cancelled,
                positive_full_threshold=selected,yukawas=terms,bare_masses=bare,
                centre_kernel=kernel,n6_finite_window_g=gend,
                warning='Ordinary local anomalies and listed discrete tests only; no joint global-background certification')


def dft(n):
    j=np.arange(n)
    return np.exp(2j*np.pi*j[:,None]*j[None,:]/n)/np.sqrt(n)


def neutral_blocks(hscale=1., muscale=1.):
    u=dft(3)
    A=(u*np.array([.01,.02,.03]))@u.conj().T
    B=hscale*1e-4*np.array([[0,1,0],[-1,0,1],[0,-1,0]],complex)
    C=muscale*1e-4*np.eye(3)
    D=np.diag([10.,12.,15.])
    return A,B,C,D


def charge_conjugated_dictionary(audit):
    rng=np.random.default_rng(131819)
    def matrix():
        return rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
    Yn,YF,mu=matrix(),matrix(),matrix()
    hw=matrix();hw=hw-hw.T
    # Original all-left order (nu,F_L,Nc,Fc), with physical portal
    # bar(F_L) YF F_R and original all-left Nc^T mu Fc.
    S=np.zeros((12,12),complex)
    for i,j,block in [(0,2,Yn.conj()),(1,3,YF.conj()),(0,1,hw),(2,3,mu)]:
        S[3*i:3*i+3,3*j:3*j+3]=block
        S[3*j:3*j+3,3*i:3*i+3]=block.T
    perm=np.r_[np.arange(3),np.arange(9,12),np.arange(6,9),np.arange(3,6)]
    Sp=S[np.ix_(perm,perm)]
    physical=Sp[:6,6:].conj()
    expected=np.block([[Yn,hw.conj()],[mu.conj().T,YF.T]])
    audit.close('generic_charge_conjugated_Dirac_dictionary',physical,expected)
    audit.close('Weyl_reorder_no_same_phase_mass',Sp[:6,:6],np.zeros((6,6)))
    audit.close('Weyl_reorder_no_same_phase_mass',Sp[6:,6:],np.zeros((6,6)))
    audit.truth('generic_heavy_portal_transpose_is_necessary',norm(YF-YF.T)>1.)
    return dict(original_portal='bar(F_L) YF F_R',new_Dirac_heavy_block='YF.T',
                original_bare='Nc.T mu Fc',new_Dirac_lower_left='mu.conj().T',
                comment='The displayed diagonal point hides this transpose; generic matrices do not.')


def neutral_test(audit,hscale=1.,muscale=1.):
    A,B,C,D=neutral_blocks(hscale,muscale)
    M=np.block([[A,B],[C,D]])  # The author's displayed point is physical.
    Dw=M.conj()  # Transport to all-left coefficient; not an arbitrary relabel.
    U,s,Vh=np.linalg.svd(M)
    order=np.argsort(s)
    s=s[order]; U=U[:,order]; Vh=Vh[order]
    audit.close('physical_SVD',U@np.diag(s)@Vh,M)
    audit.close('Weyl_physical_conjugation',M@M.conj().T,(Dw@Dw.conj().T).conj())
    audit.close('antisymmetric_LL_rho',B.T,-B)
    Sw=np.block([[np.zeros_like(Dw),Dw],[Dw.T,np.zeros_like(Dw)]])
    sw=np.sort(np.linalg.svd(Sw,compute_uv=False))
    audit.close('symmetric_Weyl_double_spectrum',sw,np.repeat(s,2))
    phase=np.r_[-1j*np.ones(6),1j*np.ones(6)]
    audit.close('full_neutral_residual_invariance',phase[:,None]*Sw*phase[None,:],Sw)
    delta=np.block([[np.zeros((3,3)),B],[C,np.zeros((3,3))]])
    radius=norm(delta)
    base=np.array([.01,.02,.03,10.,12.,15.])
    audit.truth('full_rank_gap_perturbation',s[0]>=.01-radius and s[2]<.03+radius+1e-12 and s[3]>10-radius)
    H=M@M.conj().T
    P=U[:,:3]@U[:,:3].conj().T
    K=U[:3,:3]@np.diag(s[:3]**2)@U[:3,:3].conj().T
    cp=float(np.imag(K[0,1]*K[1,2]*K[2,0]))
    audit.truth('physical_nonzero_CP',abs(cp)>1e-12)
    Uw,sw2,_=np.linalg.svd(Dw)
    ix=np.argsort(sw2)[:3]
    Kw=Uw[:3,ix]@np.diag(sw2[ix]**2)@Uw[:3,ix].conj().T
    cpw=float(np.imag(Kw[0,1]*Kw[1,2]*Kw[2,0]))
    audit.close('CP_chart_sign',cpw,-cp,tol=1e-18)
    audit.close('full_mass_squared_projection',K,(H@P)[:3,:3])
    audit.truth('three_light_projection_not_exactly_unitary',norm(U[:3,:3]@U[:3,:3].conj().T-np.eye(3))>1e-13)
    # Exact finite-z elimination of the heavy rows of the Hermitian mass square.
    Hll,Hlh,Hhh=H[:3,:3],H[:3,3:],H[3:,3:]
    gamma=float(np.linalg.eigvalsh(Hhh)[0])
    def schur(z):
        return Hll-z*np.eye(3)-Hlh@np.linalg.solve(Hhh-z*np.eye(3),Hlh.conj().T)
    finite=[]
    for z in [.00005,.0002,.001]:
        S=schur(z)
        inverse=np.linalg.inv(H-z*np.eye(6))[:3,:3]
        audit.close('exact_resolvent_block',inverse,np.linalg.inv(S),tol=2e-8)
        diff=norm(S-schur(0)+z*np.eye(3))
        bound=z*norm(Hlh)**2/(gamma*(gamma-z))
        audit.truth('finite_z_heavy_elimination_bound',diff<=bound*(1+1e-5)+1e-14)
        finite.append(dict(z=z,error=diff,bound=bound))
    naive=A-B@np.linalg.solve(D,C)
    naive_m=np.sort(np.linalg.svd(naive,compute_uv=False))
    audit.truth('zero_momentum_Schur_not_physical_spectrum',norm(naive_m-s[:3])>1e-14)
    # Full physical source for log(w_eta): A changes; B,C,D do not.
    source=np.zeros_like(M); source[:3,:3]=A
    Mp=source
    Hp=Mp@M.conj().T+M@Mp.conj().T
    eigprime=np.real(np.diag(U.conj().T@Hp@U))
    eps=1e-4
    def eigs_at(t):
        Dt=M.copy(); Dt[:3,:3]=np.exp(t)*A
        return np.sort(np.linalg.svd(Dt,compute_uv=False))**2
    fd=(-eigs_at(2*eps)+8*eigs_at(eps)-8*eigs_at(-eps)+eigs_at(-2*eps))/(12*eps)
    audit.close('complete_neutral_source_Hellmann_Feynman',eigprime,fd,tol=2e-8)
    E,F,G=Hll,Hlh,Hhh
    Ginv=np.linalg.inv(G)
    K0=E-F@Ginv@F.conj().T
    C0=F@Ginv@Ginv@F.conj().T
    zstar=norm(E); eta=zstar/gamma
    def generalized(kernel,metric):
        vals,vecs=np.linalg.eigh(metric)
        invroot=(vecs/np.sqrt(vals))@vecs.conj().T
        return np.linalg.eigvalsh(invroot@kernel@invroot)
    lower=generalized(K0,np.eye(3)+C0/(1-eta))
    upper=generalized(K0,np.eye(3)+C0)
    audit.truth('exact_count_low_and_high_poles',np.all(s[:3]**2<=zstar+1e-12) and np.all(s[3:]**2>=gamma-1e-10))
    audit.truth('generalized_lower_pole_bound',np.all(lower<=s[:3]**2+1e-14))
    audit.truth('generalized_upper_pole_bound',np.all(s[:3]**2<=upper+1e-14))
    # Check the source after elimination, including both mixing vertices.
    Ep,Fp,Gp=Hp[:3,:3],Hp[:3,3:],Hp[3:,3:]
    z=.0002
    R=np.linalg.inv(G-z*np.eye(3))
    Sp=Ep-Fp@R@F.conj().T-F@R@Fp.conj().T+F@R@Gp@R@F.conj().T
    def eliminated_at(t):
        Mt=M.copy();Mt[:3,:3]=np.exp(t)*A
        Ht=Mt@Mt.conj().T
        return Ht[:3,:3]-z*np.eye(3)-Ht[:3,3:]@np.linalg.solve(Ht[3:,3:]-z*np.eye(3),Ht[3:,:3])
    fdS=(-eliminated_at(2*eps)+8*eliminated_at(eps)-8*eliminated_at(-eps)+eliminated_at(-2*eps))/(12*eps)
    audit.close('complete_eliminated_neutral_source',Sp,fdS,tol=2e-9)
    return dict(hscale=hscale,muscale=muscale,masses=s.tolist(),
                physical_CP=cp,Weyl_CP=cpw,light_projection_defect=norm(U[:3,:3]@U[:3,:3].conj().T-np.eye(3)),
                naive_Schur_mass_error=norm(naive_m-s[:3]),finite_elimination=finite,
                perturbation_norm=radius,source_derivatives=eigprime.tolist(),
                low_squared_pole_bounds=dict(lower=lower.tolist(),upper=upper.tolist(),eta=eta))


def visible_pole_enclosure(audit):
    U=dft(3)
    A=(U*np.array([.05,.08,.11]))@U.conj().T
    B=.16*np.array([[0,1,1j],[-1,0,.5],[-1j,-.5,0]])
    C=.12*np.array([[1,.2,0],[0,1,.1j],[.2,0,1]])
    D=np.diag([1.,1.5,2.])
    M=np.block([[A,B],[C,D]])
    H=M@M.conj().T
    E,F,G=H[:3,:3],H[:3,3:],H[3:,3:]
    gi=np.linalg.inv(G)
    K0=E-F@gi@F.conj().T;C0=F@gi@gi@F.conj().T
    zstar=norm(E);gamma=np.linalg.eigvalsh(G)[0];eta=zstar/gamma
    def generalized(metric):
        eigen,W=np.linalg.eigh(metric)
        inverse_root=(W/np.sqrt(eigen))@W.conj().T
        return np.linalg.eigvalsh(inverse_root@K0@inverse_root)
    lower=generalized(np.eye(3)+C0/(1-eta))
    upper=generalized(np.eye(3)+C0)
    poles=np.sort(np.linalg.svd(M,compute_uv=False))**2
    audit.truth('visible_enclosure_contract',eta<1 and poles[0]>0 and np.max(upper-lower)>1e-6)
    audit.truth('visible_enclosure_lower',np.all(lower<=poles[:3]+1e-14))
    audit.truth('visible_enclosure_upper',np.all(poles[:3]<=upper+1e-14))
    for z in poles[:3]:
        S=E-z*np.eye(3)-F@np.linalg.solve(G-z*np.eye(3),F.conj().T)
        audit.truth('visible_exact_low_Schur_root',np.min(np.abs(np.linalg.eigvalsh(S)))<1e-12)
    return dict(eta=float(eta),exact_low_squared_poles=poles[:3].tolist(),
                lower=lower.tolist(),upper=upper.tolist(),
                comment='A separate allowed-matrix diagnostic makes the enclosure wider than floating-point noise; not the displayed common point.')


def su3generators():
    out=[]
    for i,j in [(0,1),(0,2),(1,2)]:
        x=np.zeros((3,3),complex); x[i,j]=x[j,i]=.5
        y=np.zeros((3,3),complex); y[i,j]=-.5j; y[j,i]=.5j
        out.extend([x,y])
    t3=np.diag([.5,-.5,0.])
    t8=np.diag([1.,1.,-2.])/(2*np.sqrt(3))
    return out+[t3,t8]


def gauge_model(w,g,gx):
    w=np.asarray(w,float)
    t3=np.array([.5,-.5,0.]);t8=np.array([1.,1.,-2.])/(2*np.sqrt(3))
    X=np.array([-1.,2.,-1.])/3
    v=np.column_stack((g*t3,g*t8,gx*X))
    M2=2*v.T@np.diag(w*w)@v
    d=np.array([0.,g/np.sqrt(3),gx/3])
    lam,E=np.linalg.eigh(M2)
    Em=E[:,1:]; lm=lam[1:]
    photon=np.array([1/g,-1/(np.sqrt(3)*g),1/gx])
    photon/=np.linalg.norm(photon)
    return v,M2,d,lm,Em,photon


def gauge_kernel(w,g,gx,z=0.):
    v,M2,d,lm,Em,p=gauge_model(w,g,gx)
    R=(Em/(lm-z))@Em.T
    return float(d@R@d),R,(v,M2,d,lm,Em,p)


def gauge_test(audit,w,g,gx):
    w=np.asarray(w,float)
    F0,R,data=gauge_kernel(w,g,gx)
    v,M2,d,lm,Em,p=data
    audit.close('neutral_charge_sum',v.sum(axis=0),np.zeros(3))
    audit.close('nonuniversal_charge_common_difference',d,-v[2])
    audit.close('photon_null',M2@p,np.zeros(3))
    audit.close('nonuniversal_current_has_no_photon',d@p,0.)
    a,b,c=w*w
    target=1/(2*(c+a*b/(a+b)))
    audit.close('full_neutral_exchange_sum_rule',F0,target,tol=2e-12)
    # Full 9x9 gauge mass from all Hermitian generators; the two neutral
    # massive eigenvalues are a subspace of the same action's entire spectrum.
    gs=su3generators()
    full=np.zeros((9,9))
    X=np.array([-1.,2.,-1.])/3
    for k in range(3):
        phi=np.eye(3)[:,k]*w[k]
        columns=np.column_stack([g*T@phi for T in gs]+[gx*X[k]*phi])
        full+=2*np.real(columns.conj().T@columns)
    fullspec=np.linalg.eigvalsh(full)
    audit.truth('eight_gauge_masses_one_photon',np.sum(fullspec>1e-8)==8)
    audit.close('neutral_subblock_same_action',full[np.ix_([6,7,8],[6,7,8])],M2)
    audit.close('Cartan_block_no_offdiagonal_gauge_mixing',full[np.ix_([6,7,8],range(6))],np.zeros((3,6)))
    # The T13 pair is also electromagnetically neutral, but connects an
    # ordinary quark to a heavy one; it is not sourced by ordinary-ordinary
    # bilinears in this r_eff preserving vacuum. Do not omit it from the theory.
    Qactive=np.diag([1/3.,-2/3.,1/3.])
    residual=np.diag([-1.,-1.,1.])
    for k in [2,3]:
        audit.close('extra_neutral_X0_charge',Qactive@gs[k]-gs[k]@Qactive,np.zeros((3,3)))
        audit.close('extra_neutral_X0_ordinary_compression',gs[k][:2,:2],np.zeros((2,2)))
        audit.close('extra_neutral_X0_residual_odd',residual@gs[k]@residual,-gs[k])
        audit.close('extra_neutral_X0_mass',full[k,k],g*g*(w[0]**2+w[2]**2)/2)
    # Physical left-quark representations: anti-triplet upper/lower ordering
    # is (d,-u,D), so the ordinary u has positive T3 in both families.
    up_plus=np.array([g/2,g/(2*np.sqrt(3)),gx/3])
    up_minus=np.array([g/2,-g/(2*np.sqrt(3)),0.])
    dn_plus=np.array([-g/2,g/(2*np.sqrt(3)),gx/3])
    dn_minus=np.array([-g/2,-g/(2*np.sqrt(3)),0.])
    audit.close('up_down_same_nonuniversal_vector',up_plus-up_minus,d)
    audit.close('up_down_same_nonuniversal_vector',dn_plus-dn_minus,d)
    Uu=dft(3)@np.diag(np.exp(1j*np.array([.2,-.3,.4])))
    Ud=dft(3).conj()@np.diag(np.exp(1j*np.array([-.1,.3,.6])))
    P=np.diag([1.,0.,0.]);V=Uu.conj().T@Ud
    Pu=Uu.conj().T@P@Uu;Pd=Ud.conj().T@P@Ud
    for col in range(2):
        e=Em[:,col];kappa=float(d@e)
        Gu=Uu.conj().T@np.diag([up_plus@e,up_minus@e,up_minus@e])@Uu
        Gd=Ud.conj().T@np.diag([dn_plus@e,dn_minus@e,dn_minus@e])@Ud
        Nu=Gu-(up_minus@e)*np.eye(3);Nd=Gd-(dn_minus@e)*np.eye(3)
        audit.close('neutral_flow_mass_eigenstate_dictionary',Nu,kappa*Pu)
        audit.close('two_sector_flow_similarity',Nd,V.conj().T@Nu@V)
        audit.close('neutral_flow_rank_one_polynomial',Nu@Nu,kappa*Nu)
        audit.close('neutral_flow_rank_one_polynomial',Nd@Nd,kappa*Nd)
    # A joint parameter path changes the background and both Yukawa left
    # rotations. The four-current tensor derivative must include all three.
    Au=np.array([[0,.2,.1j],[-.2,0,.3], [.1j,-.3,0]],complex)
    Ad=np.array([[0,.1j,.3],[.1j,0,-.2],[-.3,.2,0]],complex)
    def unitary(A,t):
        ev,Q=np.linalg.eigh(1j*A)
        return (Q*np.exp(-1j*t*ev))@Q.conj().T
    Pup=Uu.conj().T@(P@Au-Au@P)@Uu
    Pdp=Ud.conj().T@(P@Ad-Ad@P)@Ud
    Fprime=-4*w[2]**2*F0**2
    tensorprime=Fprime*np.kron(Pu,Pd)+F0*(np.kron(Pup,Pd)+np.kron(Pu,Pdp))
    def tensor_at(t):
        ww=w.copy();ww[2]*=np.exp(t)
        uu=unitary(Au,t)@Uu;dd=unitary(Ad,t)@Ud
        ff=gauge_kernel(ww,g,gx)[0]
        return ff*np.kron(uu.conj().T@P@uu,dd.conj().T@P@dd)
    step=1e-4
    td=(-tensor_at(2*step)+8*tensor_at(step)-8*tensor_at(-step)+tensor_at(-2*step))/(12*step)
    audit.close('both_Yukawa_rotations_and_background_source',tensorprime,td,tol=2e-9)
    finite=[]
    for fraction in [.01,.2,.7]:
        z=fraction*lm[0]
        F,Rz,_=gauge_kernel(w,g,gx,z)
        audit.truth('finite_neutral_window_bound',F>=F0*(1-1e-12) and F<=F0/(1-fraction)*(1+1e-12))
        # Full parameter derivative includes both external charge/current and
        # the complete mixed propagator. Differentiating one heavy mass is not enough.
        source=[]
        eps=1e-4
        for index in range(5):
            vp=np.zeros_like(v);dp=np.zeros_like(d)
            if index<3:
                Mprime=4*w[index]**2*np.outer(v[index],v[index])
            else:
                coords=[0,1] if index==3 else [2]
                vp[:,coords]=v[:,coords];dp[coords]=d[coords]
                Mprime=2*(vp.T@np.diag(w*w)@v+v.T@np.diag(w*w)@vp)
            analytical=float(2*dp@Rz@d-d@Rz@Mprime@Rz@d)
            def at(t):
                ww=w.copy();gg=g;xx=gx
                if index<3: ww[index]*=np.exp(t)
                elif index==3: gg*=np.exp(t)
                else: xx*=np.exp(t)
                return gauge_kernel(ww,gg,xx,z)[0]
            finite_difference=(-at(2*eps)+8*at(eps)-8*at(-eps)+at(-2*eps))/(12*eps)
            audit.close('full_mixed_neutral_source',analytical,finite_difference,tol=2e-8)
            source.append(analytical)
        finite.append(dict(z=z,F=F,relative_increment=(F-F0)/F0,source=source))
    # The zero-energy sum is independent of nonzero gauge couplings, while
    # its single-heavy-eigenstate contribution generally is not this sum.
    heavy=float((d@Em[:,1])**2/lm[1])
    audit.truth('single_heavy_not_full_neutral_sum',abs(heavy-F0)>1e-12*F0)
    for coupling in [0,1]:
        vp=v.copy()*0;dp=d.copy()*0
        ix=[0,1] if coupling==0 else [2]
        vp[:,ix]=v[:,ix];dp[ix]=d[ix]
        Mp=2*(vp.T@np.diag(w*w)@v+v.T@np.diag(w*w)@vp)
        audit.close('zero_energy_gauge_source_cancellation',2*dp@R@d-d@R@Mp@R@d,0.,tol=2e-12)
    return dict(w=w.tolist(),g=g,gX=gx,Cartan_neutral_mass_squared=lm.tolist(),
                full_gauge_mass_squared=fullspec.tolist(),F0=F0,closed_F0=target,
                heavy_only_F0=heavy,finite=finite)


def scalar_vacuum_test(audit):
    # Independent route: exact quadratic expansion in the 18 ambient real
    # fields, not a numerical Hessian in the author's Gram coordinates.
    w=np.array([1.,1.,100.]);lam=.2;k=.3;f=.01;eps=.02
    W=float(np.prod(w));background=np.diag(w).astype(complex)
    def delta(x):
        return x[:9].reshape(3,3)+1j*x[9:].reshape(3,3)
    def quadratic(x):
        d=delta(x)
        value=4*lam*np.sum(w*w*np.real(np.diag(d))**2)
        for i,j in itertools.combinations(range(3),2):
            zij=w[i]*d[i,j]+w[j]*d[j,i].conjugate()
            value+=k*abs(zij)**2
            if (i,j)==(0,2): value+=2*eps*np.real(zij**2)
        value+=f*W*sum(np.vdot(d[:,i],d[:,i]).real/w[i]**2 for i in range(3))
        A=d/w[:,None]
        value-=f*W*np.real(np.trace(A)**2-np.trace(A@A))
        return float(value)
    eye=np.eye(18)
    H=np.zeros((18,18))
    for i in range(18):
        H[i,i]=2*quadratic(eye[i])
        for j in range(i):
            H[i,j]=H[j,i]=quadratic(eye[i]+eye[j])-quadratic(eye[i])-quadratic(eye[j])
    eigen=np.linalg.eigvalsh(H)
    audit.truth('scalar_ten_positive_eight_gauge_flat',np.sum(eigen>1e-7)==10 and np.max(np.abs(eigen[:8]))<1e-9)
    # Gauge orbit directions are independently generated from actual fields.
    orbit=[]
    for T in su3generators():
        d=1j*T@background
        orbit.append(np.r_[d.real.ravel(),d.imag.ravel()])
    X=np.diag([-1.,2.,-1.])/3
    d=1j*background@X
    orbit.append(np.r_[d.real.ravel(),d.imag.ravel()])
    orbit=np.column_stack(orbit)
    audit.truth('eight_broken_gauge_directions',np.linalg.matrix_rank(orbit,tol=1e-9)==8)
    audit.close('scalar_Hessian_annihilates_gauge_orbit',H@orbit,np.zeros_like(orbit),tol=2e-9)
    def potential(Phi):
        G=Phi.conj().T@Phi
        x=G.diagonal().real
        value=lam*np.sum((x-w*w)**2)
        value+=k*sum(abs(G[i,j])**2 for i,j in itertools.combinations(range(3),2))
        value+=2*eps*np.real(G[0,2]**2)
        value+=np.sum(f*W/(w*w)*(x-w*w))-2*f*(np.linalg.det(Phi).real-W)
        return float(value)
    rng=np.random.default_rng(131818)
    for _ in range(12):
        x=rng.normal(size=18);x/=np.linalg.norm(x)
        step=1e-3
        central=(potential(background+step*delta(x))+potential(background-step*delta(x)))/2
        audit.close('ambient_scalar_quadratic_against_exact_potential',central/step**2,quadratic(x),tol=2e-6)
        linear=(potential(background+step*delta(x))-potential(background-step*delta(x)))/(2*step)
        # Odd cubic terms need not vanish at finite step; the first variation does.
        audit.truth('stationary_scalar_first_variation_finite_error',abs(linear)<1e-4)
    # r_eff fixes the background but acts nontrivially on some fluctuations.
    rg=np.diag([-1.,-1.,1.]);rs=np.diag([-1.,-1.,1.])
    audit.close('residual_symmetry_full_scalar_background',rg@background@rs,background)
    return dict(parameters=dict(w=w.tolist(),lambda_all=lam,k_all=k,f=f,epsilon=eps),
                ambient_real_Hessian_eigenvalues=eigen.tolist(),broken_gauge_rank=8,
                smallest_positive_eigenvalue=float(eigen[8]),
                scope='Strict local minimum transverse to the eight broken directions; residual electromagnetism still acts on some fluctuations. Not unique global minimum or single-Higgs decoupling.')


def run():
    audit=Audit()
    counts=field_tests(audit)
    dictionary=charge_conjugated_dictionary(audit)
    spectra=[neutral_test(audit,*scales) for scales in [(1,1),(20,30),(0,30)]]
    enclosure=visible_pole_enclosure(audit)
    audit.truth('same_r_different_permitted_neutral_spectrum',norm(np.array(spectra[0]['masses'])-np.array(spectra[1]['masses']))>1e-9)
    gauges=[gauge_test(audit,*point) for point in [([1.,1.,100.],.65,.4),([1.2,.8,5.],.7,.3),([2.,1.,7.],.4,.8)]]
    vacuum=scalar_vacuum_test(audit)
    return dict(status='PASS_WITH_SCOPE',round=1318,independent_of_author_code=True,
                assertion_count=audit.count,check_categories=len(audit.errors),
                max_normalized_error=max(audit.errors.values()),errors=audit.errors,
                exact_field_enumeration=counts,charge_conjugated_dictionary=dictionary,
                neutral_complete_spectra=spectra,visible_pole_enclosure=enclosure,
                common_neutral_gauge_kernel=gauges,scalar_vacuum=vacuum,
                bindings={'independent_check.py':sha(__file__),'admission':sha(ADMISSION)},
                scope='Declared 331 representations and background: local anomaly/leading running, exact tree quadratic spectra and finite neutral response. No cognitive uniqueness, full global anomaly, single-Higgs decoupling or experimental certificate.')


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    result=run()
    if args.write_new:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write('\n')
    elif OUT.exists():
        if result != json.loads(OUT.read_text(encoding='utf-8')):
            raise AssertionError('Independent saved result differs from recomputation')
    print(json.dumps({k:result[k] for k in ['status','assertion_count','check_categories','max_normalized_error']},ensure_ascii=False))


if __name__=='__main__':
    main()
