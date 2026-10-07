"""931 delivery: original charged matter, polarization and local matching audit."""
from pathlib import Path
import sys,json,hashlib,re,ast,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout
TARGET=HERE/'research_round_931_checks.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(writing=False):
    frozen={}
    def add(d):
        for k,v in d.items():
            assert k not in frozen or frozen[k]==v,k
            frozen[k]=v
    for n in range(776,931):
        old=read(STAGE/f'{n}/research_round_{n}_checks.json')
        for k in ('frozen_inputs','new_scientific_and_entry_files'):add(old[k])
    for path in ('875/drafts/clock_transport_working_checks.json','878/drafts/working_checks.json','884/drafts/working_checks.json','887/drafts/working_checks.json','894/drafts/working_checks.json','896/drafts/working_checks.json','897/drafts/working_checks.json','899/drafts/working_checks.json','900/drafts/working_checks.json','904/drafts/working_checks.json','907/drafts/working_checks.json','909/drafts/working_checks.json','912/drafts/working_checks.json','913/drafts/working_checks.json','915/drafts/working_checks.json','917/drafts/working_checks.json'):add(read(STAGE/path)['files'])
    au=read(STAGE/'909/drafts/scope_reaudit_checks.json');add(au['preserved_files']);add(au['evidence_and_audit_hashes'])
    for path in ('911/drafts/finite_scope_checkpoint_checks.json','912/drafts/effective_scope_reaudit_checks.json'):add(read(STAGE/path)['evidence_and_audit_hashes'])
    au=read(STAGE/'914/drafts/finite_scope_decision_checks.json');add(au['evidence_hashes']);add(au['new_document_and_verifier_hashes'])
    au=read(STAGE/'919/drafts/effective_scope_after_918_checks.json');add(au['frozen_inputs_verified']);add(au['new_document_and_verifier_hashes'])
    for rel,digest in frozen.items():assert sha(ROOT/rel)==digest,rel
    layout=Layout().verify()
    result=read(HERE/'charged_matter_response_results.json')
    for rel,digest in result['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    assert result['round']==931 and result['all_scientific_checks_passed']
    assert result['q_star']==.27
    assert abs(result['physical_Dirac_charge_squared_sum']-8/3)<1e-14
    for key in ('frozen_full_mass_spectrum_error','current_projector_spectral_error','Feynman_cut_spectral_error','independent_dispersion_error','radial_scaling_error'):assert result[key]<1e-13,key
    assert result['radial_source_derivative_error']<1e-9
    assert len(result['finite_momentum_rows'])==6
    for row in result['finite_momentum_rows']:
        assert row['remainder']<=row['remainder_bound']*(1+1e-11)
        assert row['finite_difference_error']<1e-9
    assert abs(result['C1']-.7043879163161342)<1e-12
    assert .00046<result['subtracted_response_difference']<.00047
    assert result['spectral_rows'][0]['density']==0 and result['spectral_rows'][1]['density']>0
    assert result['finite_matching_slope_shift']==.07
    assert abs(result['finite_matching_change']-.0091)<1e-14
    for key in ('same_local_matching_at_q_star_but_different_radial_source','finite_local_polynomials_cannot_cancel_pair_cut','original_matter_used_no_extra_dilaton_species','scalar_source_and_EM_kernel_from_same_determinant','spacelike_subtraction_removes_only_zero_momentum_local_coefficient','no_radial_axion_generated_with_fixed_mass_phases'):assert result[key]
    for key in ('bare_gauge_normalization_selected','all_higher_derivative_matching_coefficients_fixed','full_interacting_SM_or_hadron_response_computed','graph_continuum_mapping_or_actual_reference_instrument_constructed','free_on_shell_photon_speed_modified_claimed','all_six_protocols_realized_in_continuum_candidate','full_backreaction_or_GR_generated','whole_stage_completed','full_goal_completed'):assert not result[key]
    note=STAGE/'research_note_931.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==22 and re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,12)]
    for phrase in ('总响应尚未唯一','零动量标量插入','局部匹配反例','不含QCD相互作用','成熟','停止','整体目标保持，未结项'):
        assert phrase in prose,phrase
    working=HERE/'drafts/matter_response_working.md'
    assert '一次有限检验及停止条件' in working.read_text('utf-8')
    newdocs=[note,working,STAGE/'932/drafts/STATUS.md']
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
    ledger=nav[-1].read_text('utf-8-sig')
    mapping=ledger.split('## 六条共同协议：全局缺口对应与检验优先级（截至931）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    assert '当前认知假说补充：共同指认的六条协议' in nav[1].read_text('utf-8-sig')
    assert '已采用六条共同认知协议' in nav[2].read_text('utf-8-sig')
    links=0
    for doc in newdocs+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for q in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',q.name)
        if m and q.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=931)==list(range(1,932))
    assert '001—931轮共931份' in (STAGE.parent/'README.md').read_text('utf-8-sig')
    assert '231—931的701份' in (STAGE/'阶段成果总览.md').read_text('utf-8-sig')
    for q in nav:
        s=q.read_text('utf-8-sig');assert '931：原带荷物质约束电磁响应' in s and '932/drafts/STATUS.md' in s
    oldfiles=[STAGE/'930/research_round_930_checks.json',HERE/'drafts/STATUS.md']
    newfiles=newdocs+[Path(__file__),HERE/'charged_matter_response.py',HERE/'charged_matter_response_results.json']
    for q in newfiles:
        if q.suffix=='.py':ast.parse(q.read_text('utf-8'))
    previous=read(STAGE/'930/research_round_930_checks.json')
    out=dict(round=931,date='2026-10-06',all_delivery_checks_passed=True,fresh_test_groups=1,formal_reports=931,
      cumulative_numbered_test_groups_from_930=previous['cumulative_numbered_test_groups_from_929']+1,
      historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,local_links_checked=links,
      original_mass_and_charge_data_reused=True,
      current_spectrum_parameter_integral_and_dispersion_cross_checked=True,
      same_matter_soft_scalar_source_connected=True,
      finite_momentum_expansion_bound_analytical_not_grid_inference=True,
      local_matching_freedom_explicit_and_not_declared_program_failure=True,
      no_new_fundamental_response_species_required_for_this_contribution=True,
      full_on_shell_photon_dispersion_not_claimed=True,
      full_interacting_SM_and_GR_not_claimed=True,
      EM_technical_branch_stopped=True,
      full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
      frozen_inputs={str(q.relative_to(ROOT)):sha(q) for q in oldfiles},new_scientific_and_entry_files={str(q.relative_to(ROOT)):sha(q) for q in newfiles})
    if not writing:
        before=read(TARGET)
        for k in ('frozen_inputs','new_scientific_and_entry_files'):assert before[k]==out[k],k
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    if a.write:assert not TARGET.exists()
    d=run(a.write)
    if a.write:TARGET.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in d.items() if k not in ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
