"""915 working: actual finite-quadrature sensitivity of the Ward-reduced source.
These are coefficients of a conditional algebraic error bound, not estimates
of the unknown physical errors, and not a continuum certificate.
"""
from pathlib import Path
import sys,json,hashlib,argparse,time
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent
sys.path.insert(0,str(STAGE/'914'))
import weighted_ward_pairing as prior
ex=prior.ex;weighted=prior.weighted;ms=ex.ms
TARGET=HERE/'weighted_source_budget_results.json'

def coefficients(src,w,dw):
    rows=[];g=src.jets['g'];Mg=src.data['metric'];B=src.data['B']
    for s,C,k in zip(src.sources,src.currents,src.ks):
        q={}
        q['metric_value']=w[:,None,None]*s['g']
        q['scalar_value']=w[:,None]*s['phi']+np.einsum('bm,bmi->bi',dw,s['dphi'])
        q['connection_value']=w[:,None,None]*s['A']+np.einsum('br,brma->bma',dw,s['dA'])
        q['scalar_first_jet']=w[:,None,None]*s['dphi']
        q['connection_first_jet']=w[:,None,None,None]*s['dA']
        q['xi']=k[:,3,None]*(np.einsum('bmn,bm,bnr->br',Mg,dw,g)+np.einsum('bmn,bn,bmr->br',Mg,dw,g))
        q['eta']=np.einsum('bma,bm->ba',C,dw)
        q['Deta']=np.zeros_like(src.jets['A'],dtype=complex)
        q['Deta'][...,:8]=k[:,3,None,None]*np.einsum('bmni,bn->bmi',B-B.swapaxes(1,2),dw)
        rows.append(q)
    return rows

def evaluate(coeff,inputs):return sum(np.sum(c*inputs[k]) for k,c in coeff.items())

def independent_source(src,channel,w,dw,inputs):
    v=dict(g=inputs['metric_value'],phi=inputs['scalar_value'],A=inputs['connection_value'],dphi=inputs['scalar_first_jet'],dA=inputs['connection_first_jet'])
    wx=ex.oldref.multiply(v,w,dw)
    out=np.sum(ms.pair(src.sources[channel],wx))
    xi=inputs['xi'];eta=inputs['eta'];Deta=inputs['Deta']
    low=np.einsum('bmn,bn->bm',src.jets['g'],xi)
    qg=dw[:,:,None]*low[:,None,:]+dw[:,None,:]*low[:,:,None]
    qa=dw[:,:,None]*eta[:,None,:]
    qf=dw[:,None,:,None]*Deta[:,:,None,:]-dw[:,:,None,None]*Deta[:,None,:,:]
    qM=np.sum(src.data['metric']*qg,axis=(1,2))+np.sum(src.data['B']*qf[...,:8],axis=(1,2,3))
    return out+np.sum(src.currents[channel]*qa)+np.sum(src.ks[channel][:,3]*qM)

def run():
    start=time.time();bg=ex.mb.Background(17);src=ex.loop.LoopSource(bg,2,8);src.kernels()
    u,z,mode=weighted.receiver_modes();w,dw,stats=weighted.weight_jets(u,z,src.points)
    old=json.loads((STAGE/'914/weighted_ward_pairing_results.json').read_text('utf-8'))
    assert abs(stats['maximum_weight_abs']-old['weight']['maximum_weight_abs'])<1e-12
    coeff=coefficients(src,w,dw);names=list(coeff[0]);rows=[];rng=np.random.default_rng(915)
    for channel,c in enumerate(coeff):
        factors={k:float(np.sum(abs(a))) for k,a in c.items()}
        # Diagnostic errors in the finite input array, not a chosen physical
        # tolerance. Complex perturbations test the safe complex-linear bound.
        eps={k:(i+1)*1e-8 for i,k in enumerate(names)}
        random={k:(rng.uniform(-1,1,a.shape)+1j*rng.uniform(-1,1,a.shape))*eps[k]/np.sqrt(2) for k,a in c.items()}
        calc=evaluate(c,random);direct=independent_source(src,channel,w,dw,random)
        bound=sum(factors[k]*eps[k] for k in names)
        worst={k:np.divide(a.conj(),abs(a),out=np.zeros_like(a),where=abs(a)>0)*eps[k] for k,a in c.items()}
        adversarial=evaluate(c,worst)
        assert abs(calc-direct)<1e-16 and abs(calc)<=bound*(1+1e-12)
        assert abs(adversarial-bound)<max(1e-16,bound*2e-13)
        rows.append(dict(channel=channel,finite_input_error_amplification=factors,diagnostic_bound=bound,
           diagnostic_random_error=float(abs(calc)),independent_formula_error=float(abs(calc-direct)),
           diagnostic_aligned_complex_error=float(abs(adversarial)),unknown_actual_physical_input_errors={k:None for k in names}))
    const=coefficients(src,np.ones_like(w),np.zeros_like(dw))
    constant_extra=max(ex.maximum(c[k]) for c in const for k in ('xi','eta','Deta'))
    assert constant_extra==0
    result=dict(round=915,status='working',date='2026-10-06',N=17,source_samples=len(src.points),weight=stats,rows=rows,
      constant_weight_extra_coefficients_max=constant_extra,
      input_groups_not_independent_physical_fields=True,
      error_scope='Fixed finite source samples, kernels, background and weight. Errors in dependent extractor inputs must be bounded consistently; unknown background, modes, support and quadrature errors are additional.',
      full_continuous_operator_norm_certified=False,actual_physical_input_errors_bounded=False,
      actual_finite_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
      elapsed_seconds=round(time.time()-start,3),
      source_hashes={str(p.relative_to(STAGE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),STAGE/'914/weighted_ward_pairing.py',STAGE/'914/weighted_receiver_pairing.py',STAGE/'908/relational_loop_source.py')})
    print(json.dumps(result,ensure_ascii=False,indent=2),flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
