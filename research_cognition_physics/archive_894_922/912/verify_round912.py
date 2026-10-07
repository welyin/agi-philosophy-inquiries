"""Verify912 same-source interface, actual defect identities and frozen history."""
from pathlib import Path
import argparse,ast,hashlib,json,re,sys
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.dont_write_bytecode=True;sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
import receiver_forced_response as science
TARGET=HERE/'research_round_912_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,912):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for p in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json',
        '887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json',
        '897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json',
        '904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json'):
        add(read(STAGE/p)['files'])
    audit=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(audit['preserved_files']);add(audit['evidence_and_audit_hashes'])
    add(read(STAGE/'911/drafts/finite_scope_checkpoint_checks.json')['evidence_and_audit_hashes'])
    add(read(HERE/'drafts/effective_scope_reaudit_checks.json')['evidence_and_audit_hashes'])
    for k,v in frozen.items():assert sha(ROOT/k)==v,k
    layout=Layout().verify()
    saved=read(science.TARGET)
    for k,v in saved['source_hashes'].items():assert sha(STAGE/k)==v,k
    assert saved['full104_forced_tangent_computed'] and saved['same_semidiscrete_source_used_throughout']
    assert not saved['actual_finite_observable_error_certified'] and not saved['full872_numerical_feedback_completed']
    assert [(x['N'],x['steps'],x['Tend']) for x in saved['rows']]==[(17,2,.00025),(17,4,.00025),(17,2,-.00025),(25,2,.00025)]
    for x in saved['rows']:
        assert x['initial']['additional_matter_counterflow_max']==0
        assert x['initial']['momentum_residual_minus_zero_mode']<1e-12
        assert not x['physical_Pi_projection_implemented'] and not x['continuum_error_certified']
        assert x['final_total_constraints']['H']>1e-4 and x['final_total_constraints']['M']>1e-4
        assert x['final_response_max']['E']>1e-5 and x['final_response_max']['phi']>1e-6
        for key in ('initial_defect','final_defect'):
            d=x[key]
            for k in ('time_jet_identity_error','source_identity_error','Dirac_identity_error'):assert d[k]<1e-12
            assert d['finite_grid_source_correction_upper_max']>=d['source_correction_max']
    fd=saved['rows'][0]['independent_finite_amplitude_checks'];ratio=fd[0]['maximum_error']/fd[1]['maximum_error']
    assert 3.9<ratio<4.1 and fd[-1]['maximum_error']<1e-9
    # Fresh same-source and defect checks on the actual original17^3 initial data.
    r=science.r
    with r.ev.ResearchRuntime(r.ev.Layout()).installed():
        import joint_reference_constraint_strata as old
        model=r.ev.Model(17,old);analytic=science.tangent.AnalyticModel(17,old)
        y,u,A,init=science.initial(model,analytic,old)
        current=science.product_defect(model,y,u)
        c=science.total_constraints(model,analytic,y,u,A)
        original,z,force,load=r.source(model,y,u)
        actual,af,al,*_=science.source_data(model,y,u)
        source_difference=max(science.norm(original[k]-actual[k]) for k in original)
        force_difference=max(science.norm(force[k]-af[k]) for k in force)
        load_difference=max(science.norm(load[k]-al[k]) for k in load)
        assert max(source_difference,force_difference,load_difference)<1e-13
        for k in c:assert abs(c[k]-saved['rows'][0]['initial_total_constraints'][k])<1e-12,k
        for k in ('time_jet_identity_error','source_identity_error','Dirac_identity_error','initial_spin_curl_decomposition_error'):assert current[k]<1e-12,k
    note=STAGE/'research_note_912.md';body=note.read_text('utf-8')
    assert body.count('$$')==16 and '\t' not in body
    assert re.findall(r'\\tag\{(\d+)\}',body)==[str(i) for i in range(1,9)]
    docs=[note,STAGE/'913/drafts/STATUS.md']+list((HERE/'drafts').glob('*.md'))
    docs += [STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]
    docs += [STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    links=0
    for doc in docs:
        prose=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',prose):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for file in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',file.name)
        if m and file.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=912)==list(range(1,913))
    oldfiles=[STAGE/'911/research_round_911_checks.json',HERE/'drafts/STATUS.md',HERE/'drafts/working_checks.json',HERE/'drafts/effective_scope_reaudit_checks.json']
    oldfiles += [STAGE/f'research_note_{n}.md' for n in (731,754,872,899,902,904,905,906,911) if (STAGE/f'research_note_{n}.md').exists()]
    newfiles=[note,Path(__file__),science.TARGET,HERE/'receiver_forced_response.py',HERE/'drafts/source_representation_probe.json',STAGE/'913/drafts/STATUS.md']
    for file in newfiles:
        if file.suffix=='.py':ast.parse(file.read_text('utf-8'))
    result=dict(round=912,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_911=3697,formal_reports=912,display_equations=8,
        frozen_inputs={str(file.relative_to(ROOT)):sha(file) for file in oldfiles},
        new_scientific_and_entry_files={str(file.relative_to(ROOT)):sha(file) for file in newfiles},
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
        same906_source_difference=source_difference,same905_force_difference=force_difference,same_constraint_load_difference=load_difference,
        fresh_original17_grid_defect_identities_checked=True,independent_amplitude_error_ratio=ratio,
        full104_forced_tangent_implemented=True,physical_A_and_continuum_error_certified=False,
        actual_finite_observable_error_certified=False,full872_numerical_feedback_completed=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False)
    if not writing:
        old=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert old[k]==result[k],k
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    if args.write:assert not TARGET.exists()
    result=run(args.write)
    if args.write:TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
