"""Publish finite state-dependent strong-field response and mechanism adoption."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1004'


def main():
    old=(STAGE/'1003/verify_round1003.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('mechanism synthesis 1003','finite strong-field response 1004')
    prefix=prefix.replace('research_round_1003_checks.json','research_round_1004_checks.json')
    prefix=prefix.replace('range(776,1003)','range(776,1004)')
    marker='    for path,keys in extra.items():'
    prefix=prefix.replace(marker,
        "    extra['1004/drafts/resource_adoption_checks.json']=('frozen_evidence_hashes','audit_file_hashes')\n"+marker,1)
    checks=r'''
    import importlib.util
    import numpy as np
    spec=importlib.util.spec_from_file_location('horizon1004',HERE/'finite_horizon_response.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'finite_horizon_response_results.json');core.compare(core.calculate(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1004
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('rigorous_quadrature_error_certified','numerical_calibration_is_four_dimensional',
        'finite_window_certifies_finite_HH_reservoir','exact_finite_window_KMS_ratio_claimed',
        'same_background_is_two_self_consistent_backreaction_solutions',
        'autonomous_physical_detector_constructed','hawking_spectrum_derived_from_cognitive_principles',
        'full_goal_completed'):assert not saved[key],key
    p=saved['parameters'];L=p['dimensionless_half_width'];Q=p['dimensionless_frequency_cut']
    assert abs(p['alpha']/p['detector_proper_acceleration']-16)<1e-13
    assert abs(p['alpha']*p['support_half_width']-8)<1e-14
    b=2*math.pi
    tail=4*L*L/math.pi*math.exp(-b*Q)/(1-math.exp(-b*Q))*(Q/b+1/b**2)
    assert math.isclose(tail,saved['spectral_H_minus_B_tail_upper_bound'],rel_tol=1e-14)
    assert 0<tail<6.80e-26
    assert math.isclose(saved['smooth_difference_at_zero'],p['alpha']**2/(24*math.pi),rel_tol=1e-14)
    # A separate scalar power-series evaluation checks the coincidence branch.
    for z in (0.,.001,.01):
        series=(1/12-z*z/240+z**4/6048-z**6/172800)/(2*math.pi)
        assert abs(float(core.smooth_kernel(np.array(z)))-series)<1e-17
    for row in saved['rows']:
        latest=row['evaluations'][-1]
        assert latest['absolute_algorithm_difference']<2.3e-13
        assert row['refinement_difference']<3.5e-13
        assert abs(2*row['H_minus_U']-latest['time_H_minus_B'])<1e-14
        assert row['H_minus_U']>0
    note=STAGE/'research_note_1004.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','累计3786','不是一般R处的固有加速度','测试场比较','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'strong_field_adoption_v1.md',HERE/'finite_horizon_response.py',
        HERE/'finite_horizon_response_results.json',Path(__file__),HERE/'drafts/horizon_adoption_decision.md',
        HERE/'drafts/horizon_review_notes.md',HERE/'drafts/publish1004.py',STAGE/'1005/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1004）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [p for p in newfiles if p.suffix=='.md']+nav:
        content=re.sub(r'\$\$.*?\$\$','',doc.read_text('utf-8-sig'),flags=re.S)
        for link in re.findall(r'\]\(([^)]+)\)',content):
            if re.match(r'^[a-zA-Z]+://',link) or link.startswith('#'):continue
            target=(doc.parent/link.split('#')[0].strip('<>')).resolve()
            assert target.exists() or (writing and target==TARGET.resolve()),(doc,link)
            links+=1
    nums=[]
    for pth in STAGE.parent.rglob('research_note_*.md'):
        m=re.fullmatch(r'research_note_(\d+).md',pth.name)
        if m and pth.parent.name.startswith('archive_'):nums.append(int(m[1]))
    assert sorted(n for n in nums if n<=1004)==list(range(1,1005))
    assert '001—1004轮共1004份' in nav[0].read_text('utf-8-sig')
    assert '231—1004的774份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1004：几何、场态与有限强场记录' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1003/research_round_1003_checks.json',HERE/'drafts/STATUS.md',
              HERE/'drafts/resource_adoption_checks.json',HERE/'drafts/NEXT.md']
    prev=read(oldfiles[0])
    out=dict(round=1004,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1004,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_1003=prev['cumulative_numbered_test_groups_from_1002']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,one_plus_one_algorithms_agree=True,
        analytic_frequency_tail_bound_verified=True,
        four_dimensional_radial_spectrum_computed=False,
        physical_backreaction_certified=False,full_goal_completed=False,
        visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(p.relative_to(ROOT)):sha(p) for p in oldfiles},
        new_scientific_and_entry_files={str(p.relative_to(ROOT)):sha(p) for p in newfiles})
    if not writing:
        before=read(TARGET)
        for key in ('frozen_inputs','new_scientific_and_entry_files'):assert before[key]==out[key],key
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run(args.write)
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    print(json.dumps({k:v for k,v in out.items() if k not in
        ('frozen_inputs','new_scientific_and_entry_files')},ensure_ascii=False,indent=2))
'''
    paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals={p:p.read_bytes() for p in paths};updates={}
    for p,raw in originals.items():
        s=raw.decode('utf-8-sig').replace('\r\n','\n')
        title='## 1004：几何、场态与有限强场记录';mark='## 1004工作审计：任务资源与关系参考'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[1004报告]({pre}research_note_1004.md)补同一传播几何下真实场态与仪器的记录差；'
            +'采用4D有限响应判别，1+1双算法校准一致，有限波包版本保新增应力。'
            +f'[结果]({pre}1004/finite_horizon_response_results.json) · '
            +f'[核验]({pre}1004/research_round_1004_checks.json)。正式1004／累计3786，整体未完成。\n\n'
            +f'[强场机制补充]({pre}1004/strong_field_adoption_v1.md)分清几何访问、场态、探测与来源；'
            +'未签收4D灰体数值、共同几何反作用或完整蒸发。停止强场细化，'
            +f'接[1005]({pre}1005/drafts/STATUS.md)回辐射输运、再吸收与组织形成反馈。'
            +'应用目标、957及条件性空间接口保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—1003轮共1003份')==1
            s=s.replace('001—1003轮共1003份','001—1004轮共1004份',1)
        if p==paths[4]:
            before='当前正式1003／累计3785，1003已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式1004／累计3786，1004已结项',1)
        if p==paths[5]:
            assert s.count('231—1003的773份')==1
            s=s.replace('# 231—1003轮阶段成果总览','# 231—1004轮阶段成果总览',1)
            s=s.replace('231—1003的773份','231—1004的774份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至1003）';assert s.count(before)==1
            s=s.replace(before,before.replace('1003','1004'),1)
            # Retain the original C titles and earlier evidence instead of redefining the ledger.
            updates_by_id={
                'C13':('1004采用实际场态与有限探测的强场接口；几何访问未选定量子关联',
                    '1+1仅校准，微观熵／完整蒸发未处理；同背景不同态非共同反作用解'),
                'C19':('供给按任务分类；1004强场态、有限波包及仪器准备共同登记',
                    'HH无限浴非有限制备；真实参考、接触与仪器仍有输入'),
                'C20':('整体v2保共同任务；1004区分4D解析接口、1+1校准及能源应力',
                    '测试场结果未认证共同反作用；尾界不等于物理总误差')}
            for cid,(state,boundary) in updates_by_id.items():
                pattern=r'^\|'+cid+r' [^\n]*$';matches=re.findall(pattern,s,re.M);assert len(matches)>=1
                cells=matches[0].split('|');assert len(cells)==6
                cells[2]+='；1004';cells[3]=state;cells[4]=boundary
                s,count=re.subn(pattern,lambda _:'|'.join(cells),s,count=1,flags=re.M);assert count==1
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式1004／累计3786。强场采用保真实场态和仪器，'
                '同一因果结构不单独决定记录。有限响应校准不签收全部几何反馈；'
                '下一项回物质—辐射交换与形成反馈，不继续扩展黑洞计算。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1003整体v2](../../1003/overall_operation_hypothesis_v2.md)、'
                '[1004供给审计](../../research_note_1004_working.md)及'
                '[1004强场补充](../../1004/strong_field_adoption_v1.md)按实际状态、任务和来源接续。'
                '下一项补发射、吸收、输运和组织形成的共同机制。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with (HERE/'verify_round1004.py').open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 1004 finite strong-field response and mechanism.')


if __name__=='__main__':main()
