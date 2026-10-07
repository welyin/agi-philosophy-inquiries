"""982: SAME native material, two electromagnetic modes, one TT mode.

Finite matrices audit algebra and a witness. The report proves a full-Fock,
finite-time comparison bound separately; matrix convergence is not that proof.
This is a selected matched-mode EFT test, not a full SM/GR realization.
"""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'native_radiative_channels_results.json'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def kron(*args):
    a=np.array([[1.]])
    for b in args:a=np.kron(a,b)
    return a
def evolution_on(h,t,v):
    e,w=np.linalg.eigh(h)
    return w@(np.exp(-1j*t*e)[:,None]*(w.conj().T@v))
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=3e-6,abs_tol=2e-8),(a,b)
    else:assert a==b,(a,b)

def run():
    old=load('native956',STAGE/'956/native_material_interface.py')
    m=old.material(1.,.1);h6=m['H'];J=m['J'];d=m['occupancy'];om=m['omega'];D=1+J
    cs=[old.annihilate(i) for i in range(4)];ns=[x.T@x for x in cs]
    sec=np.eye(16)[:,[i for i in range(16) if i.bit_count()==2]]
    Q6=sec.T@(ns[0]+ns[1]-ns[2]-ns[3])@sec/2
    ground=math.sqrt(1-d)*m['s']+math.sqrt(d)*m['d']
    middle=Q6@m['d'];upper=-math.sqrt(d)*m['s']+math.sqrt(1-d)*m['d']
    emb=np.stack([ground,middle,upper],axis=1)
    hm=np.diag([0.,D,om]);Q=np.array([[0,math.sqrt(d),0],
        [math.sqrt(d),0,math.sqrt(1-d)],[0,math.sqrt(1-d),0]])
    X=np.array([[0.,0,1],[0,0,0],[1,0,0]])
    errors=dict(isometry=norm(emb.T@emb-np.eye(3)),
        native_h=norm((h6+J*np.eye(6))@emb-emb@hm),
        native_Q=norm(Q6@emb-emb@Q),gap_addition=abs(om-D-J),
        dipole_gap_identity=abs(J*(1-d)-D*d))
    assert max(errors.values())<1e-13
    parity=np.diag([1.,-1,1.])
    assert norm(parity@X-X@parity)<1e-14 and norm(parity@Q+Q@parity)<1e-14
    g=J/40000;T=math.pi/(2*g);r=math.sqrt(d/(1-d))
    cl=g/math.sqrt(1-d);ch=g/math.sqrt(d)
    assert abs(cl/math.sqrt(J)-ch/math.sqrt(D))<1e-16
    zcoef=(cl/g)**2/J+(ch/g)**2/D
    assert abs(zcoef-50*om)<1e-11
    phi=(1+math.sqrt(5))/2;x=g*T
    grav_rwa=(math.sin(phi*x)+math.sin(x/phi))**2/5
    cascade_rwa=(math.cos(phi*x)-math.cos(x/phi))**2/5
    ma=1/(2*J)+r+r/om+1/(2*D)+1/r+1/(r*om)+1/(2*om)
    assert ma<27 and 3+r+1/r<10 and zcoef<54
    bound=g*(95+math.pi/2*(2312+54))
    assert bound<.00372 and .38<grav_rwa<.4
    rows=[]
    for n in (4,5):
        eye=np.eye(n);ann=np.diag(np.sqrt(np.arange(1,n)),1)
        im=np.eye(3);fieldid=np.eye(n**3)
        als=[kron(ann,eye,eye),kron(eye,ann,eye),kron(eye,eye,ann)]
        nsf=[a.T@a for a in als];As=[kron(im,a) for a in als]
        h0=kron(hm,fieldid)+sum(w*kron(im,z) for w,z in zip((J,D,om),nsf))
        vl=kron(Q,als[0]+als[0].T)/math.sqrt(1-d)
        vh=kron(Q,als[1]+als[1].T)/math.sqrt(d)
        vg=kron(X,als[2]+als[2].T);V=vl+vh+vg
        l_l=np.zeros((3,3));l_l[1,2]=1
        l_h=np.zeros((3,3));l_h[0,1]=1
        l_g=np.zeros((3,3));l_g[0,2]=1
        W=sum(kron(l,a.T)+kron(l.T,a) for l,a in zip((l_l,l_h,l_g),als))
        Z=zcoef*kron(Q@Q,fieldid);H=h0+g*V+g*g*Z
        # Basis: material(g,e,u), low photon, high photon, TT graviton.
        index=lambda s,l,h,k:((s*n+l)*n+h)*n+k
        basis=np.eye(3*n**3)
        indices=[index(2,0,0,0),index(1,1,0,0),index(0,1,1,0),index(0,0,0,1),index(0,0,0,0)]
        P=basis[:,indices];Pin=P[:,[4,0]]
        wp=P.T@W@P
        expected=np.zeros((5,5))
        for a,b in ((0,1),(1,2),(0,3)):expected[a,b]=expected[b,a]=1
        assert norm(wp-expected)<1e-13 and norm(W@P-P@wp)<1e-13
        # Solve commutator only for the explicit NON-resonant linear terms.
        E=np.diag(h0);diff=E[:,None]-E[None,:];rem=V-W
        nz=abs(rem)>1e-12
        assert np.min(abs(diff[nz]))>2*J-1e-12
        A=np.zeros_like(H);A[nz]=rem[nz]/diff[nz]
        comm_error=norm(h0@A-A@h0-rem)
        assert comm_error<1e-11 and norm(A.T+A)<1e-11
        AP=norm(A@P);APin=norm(A@Pin);VA_AW=norm((V@A-A@W)@P)
        assert AP<68 and APin<27 and VA_AW<2312
        exact=evolution_on(H,T,Pin)
        # RWA propagation in its exact finite invariant subspace, not a Fock truncation.
        rw=P@evolution_on(P.T@h0@P+g*wp,T,P.T@Pin)
        isometry_error=norm(exact-rw)
        assert isometry_error<bound
        grav_mask=np.diag(kron(im,nsf[2]))>.5
        pg=float(np.sum(abs(exact[grav_mask,1])**2))
        pg_rwa=float(np.sum(abs(rw[grav_mask,1])**2))
        cascade_index=index(0,1,1,0)
        assert abs(pg_rwa-grav_rwa)<1e-8
        assert abs(abs(rw[cascade_index,1])**2-cascade_rwa)<1e-8
        energy_before=Pin[:,1]@H@Pin[:,1]
        energy_after=np.vdot(exact[:,1],H@exact[:,1]).real
        energy_error=abs(energy_after-energy_before)
        assert energy_error<1e-10
        material_energy=float(np.vdot(exact[:,1],kron(hm,fieldid)@exact[:,1]).real)
        field_energies=[float(np.vdot(exact[:,1],w*kron(im,z)@exact[:,1]).real)
            for w,z in zip((J,D,om),nsf)]
        interaction_energy=float(np.vdot(exact[:,1],(g*V+g*g*Z)@exact[:,1]).real)
        q=(As[2]+As[2].T)/math.sqrt(2*om)
        p=-1j*math.sqrt(om/2)*(As[2]-As[2].T)
        # Test on a finite-support core, where truncation has no boundary term.
        source_residual=norm((1j*(H@p-p@H)+om*om*q+g*math.sqrt(2*om)*kron(X,fieldid))@P)
        velocity_residual=norm((1j*(H@q-q@H)-p)@P)
        assert max(source_residual,velocity_residual)<1e-12
        rows.append(dict(occupation_dimension_per_mode=n,full_dimension=len(H),
            commutator_error=comm_error,A_P_norm=AP,A_Pin_norm=APin,
            VA_minus_AW_P_norm=VA_AW,isometry_error=isometry_error,
            graviton_probability=pg,rwa_graviton_probability=pg_rwa,
            material_excitation_energy=material_energy,field_energies=field_energies,
            interaction_and_selfpolarization_energy=interaction_energy,
            conserved_energy_before=float(energy_before),conserved_energy_after=float(energy_after),
            energy_error=float(energy_error),source_equation_residual=source_residual,
            velocity_equation_residual=velocity_residual))
    assert abs(rows[0]['graviton_probability']-rows[1]['graviton_probability'])<1e-7
    # Quadratic TT normalization and source convention are explicit, not inferred from Q.
    pol=np.diag([1.,-1,0])/math.sqrt(2);k=np.array([0.,0,om])
    assert abs(np.trace(pol))<1e-14 and abs(np.trace(pol.T@pol)-1)<1e-14
    assert norm(k@pol)<1e-14
    sources=[STAGE/'956/native_material_interface.py',STAGE/'research_note_940.md',
        STAGE/'research_note_965.md',STAGE/'research_note_971.md',
        STAGE/'981/drafts/common_parent_contract_v1.md',STAGE/'research_note_981.md']
    return dict(round=982,all_scientific_checks_passed=True,
        parameters=dict(U=1.,v=.1,eta=0.,native_dimension=6,exact_active_dimension=3,
            J=J,Delta=D,Omega=om,occupancy=d,g_EM=g,g_TT=g,time=T,
            electric_mode_coefficients=[cl,ch],selfpolarization_coefficient_over_g2=zcoef),
        native_intertwiner_errors=errors,
        analytic=dict(A_coefficient_sum=ma,full_Fock_isometry_bound=bound,
            rational_isometry_upper=.00372,rwa_graviton_probability=grav_rwa,
            rwa_two_photon_probability=cascade_rwa,rwa_graviton_all_time_upper=.8,
            full_graviton_probability_lower=.38-.00372,
            full_graviton_probability_upper=.4+.00372,
            isolated_RWA_prediction=1.,isolated_RWA_error_lower=.59628,
            Hamiltonian_lower_bound_relative_to_m0_minus_J=-g*g/om),
        matrix_witnesses=rows,
        scope=dict(same_native_h_and_Q=True,TT_physical_normalization_specified=True,
            independent_quadrupole_matching_input=True,
            all_unknown_inputs_in_g_u_and_passive_reference_bound=True,
            full_Fock_bound_is_analytic_not_cutoff_convergence=True,
            actual_mode_environment_matching_certified=False,
            qstar_derived_from_flat_material=False,full_SM_matching_completed=False,
            nonlinear_Einstein_dynamics_completed=False,arbitrary_quantum_source_to_classical_geometry=False,
            real_graviton_detection_claim=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in sources})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes',)},ensure_ascii=False,indent=2))
