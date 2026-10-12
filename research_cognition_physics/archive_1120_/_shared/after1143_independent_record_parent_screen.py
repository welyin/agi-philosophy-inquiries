"""Non-numbered independent-record finite-cell screen; default readonly.
Finite oscillator truncation is a diagnostic, never a CCR/normality proof.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

FROZEN = {'research_cognition_physics/archive_1120_/research_note_1128.md': '2c424638a38e87bd42234c3e6c3b77af860f2d81109be0dafa1e4c36c4f8cb44', 'research_cognition_physics/archive_1120_/research_note_1129.md': '2a3c198113b36ca17f7ebd19819ca9b0d14ceea7d5ad8b5b9c6ff3fad81d2fe0', 'research_cognition_physics/archive_1120_/research_note_1133.md': '65b43a519fb3d8eff90600a6cfcfa9fb2f929b85d1cd9a3253eb497c8a88778c', 'research_cognition_physics/archive_742_763/research_note_742.md': '01382af9384d07bb0847bfaf6690d4cd051aafce9436c25b9e96e87263403329', 'research_cognition_physics/archive_742_763/research_note_743.md': '9d31e9ec2ad8f9d02ec74cd99b18ba27584ee0aaba8cbfdd9bf0d7a084bac83e', 'research_cognition_physics/archive_1120_/_shared/after1143_local_parent_preflight.md': '388fe038bbf640f5affce5411103640d9b2a4a451511d4f75745ed9de1ebe4bd'}
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).with_suffix('.json')


def annihilator(k,n=4):
    out=np.zeros((2**n,2**n),complex)
    for x in range(2**n):
        if (x>>k)&1:
            out[x^(1<<k),x]=(-1)**((x&((1<<k)-1)).bit_count())
    return out


def exp_h(h,t):
    lam,vec=np.linalg.eigh(h)
    return (vec*np.exp(-1j*t*lam))@vec.conj().T


def apply_controlled(s,x,t,state):
    """state shape material, source, columns; full material eigenspaces."""
    lam,vec=np.linalg.eigh(s)
    result=np.zeros_like(state)
    for value in sorted(set(np.round(lam,10))):
        ids=np.abs(lam-value)<1e-9
        proj=vec[:,ids]@vec[:,ids].conj().T
        part=np.einsum('ac,cbj->abj',proj,state)
        result+=np.einsum('bc,acj->abj',exp_h(x,t*value),part)
    return result


def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):return [clean(v) for v in x]
    if isinstance(x,float):return float(format(x,'.12g'))
    if isinstance(x,np.floating):return float(format(float(x),'.12g'))
    return x


def calculate():
    for name,sha in FROZEN.items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha,name
    a0,a1,p0,p1=[annihilator(k) for k in range(4)]
    ident=np.eye(16)
    az=a1.conj().T@a1
    ap=(a0+a1)/np.sqrt(2)
    ax=ap.conj().T@ap
    b=p0.conj().T@p1+p1.conj().T@p0
    ya=-1j*a0.conj().T@a1+1j*a1.conj().T@a0
    na=a0.conj().T@a0+a1.conj().T@a1
    nptr=p0.conj().T@p0+p1.conj().T@p1
    ar=np.sqrt(np.pi/2)
    kappa=np.pi/4
    qdim=40
    osc=np.zeros((qdim,qdim),complex)
    for n in range(1,qdim):osc[n-1,n]=np.sqrt(n)
    q=(osc+osc.conj().T)/np.sqrt(2)
    p=(osc-osc.conj().T)/(1j*np.sqrt(2))
    number=np.diag(np.arange(qdim))
    source=np.zeros((16,qdim,16),complex)
    source[:,0,:]=ident
    mat_results={}
    diagnostics={}
    for name,a in [('Z',az),('X',ax),('off',np.zeros((16,16)))]:
        mid=apply_controlled(b,p,ar,source)
        mid=apply_controlled(a,q,ar,mid)
        loop=apply_controlled(b,p,-ar,mid)
        loop=apply_controlled(a,q,-ar,loop)
        target=exp_h(a@b,ar*ar)
        expected=np.zeros_like(source);expected[:,0,:]=target
        residual=float(np.linalg.norm((loop-expected).reshape(16*qdim,16),2))
        assert residual<2e-10,(name,residual)
        cont=ident if name=='off' else exp_h(ya,kappa)
        total=cont@target
        joint=np.einsum('ac,cbj->abj',cont,loop)
        expected2=np.zeros_like(source);expected2[:,0,:]=total
        assert np.linalg.norm((joint-expected2).reshape(16*qdim,16),2)<2e-10
        assert np.linalg.norm(total@na-na@total)<1e-12
        assert np.linalg.norm(total@nptr-nptr@total)<1e-12
        mat_results[name]={
            'projector_residual':float(np.linalg.norm(a@a-a)),
            'commuting_independent_materials':float(np.linalg.norm(a@b-b@a)),
            'full_CAR_record_unitarity':float(np.linalg.norm(total.conj().T@total-ident)),
            'full_16_column_source_return_residual':residual,
            'data_number_conserved':float(np.linalg.norm(total@na-na@total)),
            'pointer_number_conserved':float(np.linalg.norm(total@nptr-nptr@total)),
        }
        if name=='Z':
            pl,pv=np.linalg.eigh(p)
            sine=(pv*np.sin(pl/ar))@pv.conj().T
            values=[];energies=[]
            for column in (5,6):
                v=mid[:,:,column]
                values.append(float(np.einsum('ab,bc,ac->',v.conj(),sine,v).real))
                energies.append(float(np.einsum('ab,bc,ac->',v.conj(),number,v).real))
            exact_sine=np.exp(-1/(2*np.pi))*np.sin(1.)
            assert abs(values[0])<2e-10
            assert abs(values[1]+exact_sine)<2e-10
            assert abs(energies[0]-np.pi/4)<2e-10
            assert abs(energies[1]-np.pi/2)<2e-10
            diagnostics={'source_sine_read_A0_A1':values,'source_effect_probability_gap':abs(values[1]-values[0])/2,
                'intermediate_source_N_A0_A1':energies,'input_dependent_source_energy_difference':energies[1]-energies[0],
                'exact_normal_source_effect_gap_formula':float(exact_sine/2),
                'source_truncation_dimension':qdim,
                'full_truncated_CCR_operator_defect':float(np.linalg.norm(q@p-p@q-1j*np.eye(qdim),2)),
                'truncation_used_as_normal_source_proof':False}
    assert abs(np.linalg.norm(az@ax-ax@az,2)-.5)<1e-12
    sigma=F(1,10**7);mass=F(10**16);L=F(1400)
    def delta(t):
        t=F(t);return L*(t*sigma+t*t/(4*mass*sigma))
    d5=delta(5);d2=delta(2)
    gap=F(25,72)-2*d2
    deriv_bus=F(9,7)*2*(512+16)
    deriv_cont=F(11,14)*512
    positive_shift=F(144,7)**2+F(1,4)+13
    assert max(deriv_bus,deriv_cont)<L
    assert positive_shift<500
    assert d5<F(1,1000) and gap>F(17,50)
    budget={
        'sigma':str(sigma),'control_inertia_M':str(mass),'time_T':5,
        'declared_vector_Lipschitz_bound':int(L),
        'bus_derivative_upper':str(deriv_bus),'continue_derivative_upper':str(deriv_cont),
        'autonomous_joint_isometry_bound_T5':str(d5),'autonomous_joint_isometry_bound_T5_decimal':float(d5),
        'autonomous_joint_isometry_bound_T2':str(d2),'autonomous_joint_isometry_bound_T2_decimal':float(d2),
        'actual_bounded_source_probability_gap_strict_lower':str(gap),
        'actual_bounded_source_probability_gap_lower_decimal':float(gap),
        'normal_initial_control_energy':str(mass/2+1/(8*mass*sigma*sigma)),
        'selfadjoint_positive_shift_bound':str(positive_shift),
        'energy_difference_not_inferred_from_trace_distance':True,
        'spatial_microcausality_or_four_momentum_proved':False,
    }
    return clean({'status':'PASS_PREFLIGHT_ONLY','new_numbered_round':None,'scientific_checks_added':0,
       'mathematical_diagnostic_groups':3,'frozen_sources':FROZEN,
       'full_CAR_instruments':mat_results,'noncommuting_questions_norm':float(np.linalg.norm(az@ax-ax@az,2)),
       'finite_source_diagnostic':diagnostics,'normal_autonomous_analytic_budget':budget})


def main():
    args=argparse.ArgumentParser();args.add_argument('--write-results',action='store_true');opts=args.parse_args()
    result=calculate()
    if opts.write_results:
        assert not OUT.exists(),'Refuse overwriting saved qualification results'
        OUT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        print('WROTE',OUT)
    else:
        saved=json.loads(OUT.read_text(encoding='utf-8'))
        assert result==saved,'Readonly recomputation mismatch'
        print('PASS_PREFLIGHT_ONLY: 3 groups, normal source proof analytic, no numbered round')
    print(json.dumps({'T5_error':result['normal_autonomous_analytic_budget']['autonomous_joint_isometry_bound_T5_decimal'],
         'actual_source_probability_gap_lower':result['normal_autonomous_analytic_budget']['actual_bounded_source_probability_gap_lower_decimal'],
         'worst_finite_source_return_residual':max(v['full_16_column_source_return_residual'] for v in result['full_CAR_instruments'].values())},ensure_ascii=False))

if __name__=='__main__':main()
