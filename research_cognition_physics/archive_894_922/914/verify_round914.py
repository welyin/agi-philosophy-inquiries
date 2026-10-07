"""914 delivery: archive evidence plus fresh original-background product test."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
import numpy as np
import weighted_ward_pairing as science
ex=science.ex;HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_914_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fresh_product_check():
    bg,field,_=ex.previous.rebuilt.build_pair()
    src=ex.loop.LoopSource(bg,2,8);x=src.points[::64];p=ex.second_jets(bg,x);B=len(x)
    rng=np.random.default_rng(1914)
    xi=rng.normal(size=(B,4))*.02;dx=rng.normal(size=(B,4,4))*.03
    ddx=rng.normal(size=(B,4,4,4))*.04;ddx=(ddx+ddx.swapaxes(1,2))/2
    a=rng.normal(size=(B,12))*.05;da=rng.normal(size=(B,4,12))*.04
    dda=rng.normal(size=(B,4,4,12))*.03;dda=(dda+dda.swapaxes(1,2))/2
    w=rng.normal(size=B)+.2j*rng.normal(size=B);dw=rng.normal(size=(B,4))+.2j*rng.normal(size=(B,4))
    ddw=rng.normal(size=(B,4,4))+.2j*rng.normal(size=(B,4,4));ddw=(ddw+ddw.swapaxes(1,2))/2
    gauge=ex.gauge_jets(p,xi,dx,ddx,a,da,dda)
    wx=xi*w[:,None];wdx=dx*w[:,None,None]+dw[:,:,None]*xi[:,None,:]
    wddx=ddx*w[:,None,None,None]+dw[:,:,None,None]*dx[:,None,:,:]+dw[:,None,:,None]*dx[:,:,None,:]+ddw[:,:,:,None]*xi[:,None,None,:]
    wa=a*w[:,None];wda=da*w[:,None,None]+dw[:,:,None]*a[:,None,:]
    wdda=dda*w[:,None,None,None]+dw[:,:,None,None]*da[:,None,:,:]+dw[:,None,:,None]*da[:,:,None,:]+ddw[:,:,:,None]*a[:,None,None,:]
    scaled=ex.oldref.multiply(gauge,w,dw)
    total=ex.gauge_jets(p,wx,wdx,wddx,wa,wda,wdda)
    q={k:total[k]-scaled[k] for k in total}
    low=np.einsum('bmn,bn->bm',p['g'],xi);eta=np.einsum('br,bra->ba',xi,p['A'])-a
    expectedg=dw[:,:,None]*low[:,None,:]+dw[:,None,:]*low[:,:,None]
    expectedA=dw[:,:,None]*eta[:,None,:]
    deta=np.einsum('bmr,bra->bma',dx,p['A'])+np.einsum('br,bmra->bma',xi,p['dA'])-da
    Deta=deta+ex.ms.bracket(p['A'],eta[:,None,:])
    expectedF=dw[:,None,:,None]*Deta[:,:,None,:]-dw[:,:,None,None]*Deta[:,None,:,:]
    actualF=q['dA']-q['dA'].swapaxes(1,2)+ex.ms.bracket(q['A'][:,:,None,:],p['A'][:,None,:,:])+ex.ms.bracket(p['A'][:,:,None,:],q['A'][:,None,:,:])
    errs=dict(metric=ex.maximum(q['g']-expectedg),connection=ex.maximum(q['A']-expectedA),scalar=ex.maximum(q['phi']),scalar_jet=ex.maximum(q['dphi']),curvature=ex.maximum(actualF-expectedF))
    assert max(errs.values())<2e-13
    # Independent direct old reference derivative and alpha extraction at fresh
    # actual source points, with the chosen natural feature pairing.
    _,rho,gram,Q,extract=ex.internal_data(p)
    dz=ex.feature_delta(p,gauge);lie=ex.lie_feature(p,xi,dx)
    recovered=extract({k:dz[k]-lie[k] for k in dz})
    errs['internal_extraction']=ex.maximum(recovered-a)
    errs['diff_extraction']=ex.maximum(np.linalg.solve(ex.jacobian_from_jets(p),ex.oldref.reference_delta(p,gauge)[...,None])[...,0]-xi)
    assert errs['internal_extraction']<1e-11 and errs['diff_extraction']<1e-9
    return errs

def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():assert k not in frozen or frozen[k]==v,k;frozen[k]=v
    for n in range(776,914):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(HERE/'drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    results=[read(HERE/f'{n}_results.json') for n in ('material_extractor_naturality','weighted_receiver_pairing','weighted_ward_pairing')]
    for d in results:
        for rel,digest in d['source_hashes'].items():assert sha(STAGE/rel)==digest,rel
        assert not d['full_goal_completed'] and not d['full872_response_computed']
    mat,raw,ward=results;m=mat['rows'][0]
    assert m['feature_naturality_error']<1e-12 and m['omitted_tensor_index_terms_error']>.01
    assert mat['derivative_order_bound']==dict(Ldiff=1,Lint=2,Pi=3)
    assert ward['weighted_Ward_decomposition_error']<1e-15
    assert ward['constant_weight_original_source_error']<1e-16
    assert abs(ward['constant_weight_quadrature_defect']['real'][1])>1e-7
    assert ward['previous_direct_calculation_reproduced_error']<1e-14
    assert not ward['finite_source_error_certified']
    fresh=fresh_product_check()
    note=STAGE/'research_note_914.md';s=note.read_text('utf-8');assert s.count('$$')==12
    assert re.findall(r'\\tag\{(\d+)\}',s)==[str(i) for i in range(1,7)]
    docs=[note,STAGE/'915/drafts/STATUS.md']+[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for p in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',p.name)
        if m and p.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=914)==list(range(1,915))
    oldfiles=[STAGE/'913/research_round_913_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/finite_scope_decision_checks.json']+[STAGE/f'research_note_{n}.md' for n in (773,785,803,869,870,872,899,908,911,912,913)]
    newfiles=[note,Path(__file__),STAGE/'915/drafts/STATUS.md']+[HERE/(n+s) for n in ('material_extractor_naturality','weighted_receiver_pairing','weighted_ward_pairing') for s in ('.py','_results.json')]
    for p in newfiles:
        if p.suffix=='.py':ast.parse(p.read_text('utf-8'))
    d=dict(round=914,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,cumulative_numbered_test_groups_from_913=3699,formal_reports=914,display_equations=6,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      fresh_original_background_generator_checks=fresh,
      original_nonconstant_mode_and_source_calculated=True,pure_gauge_quadrature_defect_independently_separated=True,
      frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles},
      uniform_reference_support_certified=False,actual_finite_observable_error_certified=False,full872_response_computed=False,full_goal_completed=False,
      visual_checks_performed=False,app_goal_changed=False)
    if not writing:
        old=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert old[k]==d[k],k
    return d
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
