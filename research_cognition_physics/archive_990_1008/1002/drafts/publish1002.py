"""Publish a bounded thermal-record reuse interface and lifecycle adoption."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1002'


def main():
    old=(STAGE/'1001/verify_round1001.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('nuclear formation interface 1001','thermal record reuse interface 1002')
    prefix=prefix.replace('research_round_1001_checks.json','research_round_1002_checks.json')
    prefix=prefix.replace('range(776,1001)','range(776,1002)')
    checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('reuse1002',HERE/'thermal_record_reuse.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'thermal_record_reuse_results.json');core.compare(core.calculate(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1002
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('original_charge_interaction_implements_writer','all_resources_prepared_from_thermal_states',
        'population_effect_is_qnd_under_full_active_H','infinite_retention_certified',
        'cyclic_auxiliary_restoration_certified','physical_readout_instrument_constructed',
        'joint_clock_communication_lifecycle_certified','full_goal_completed'):assert not saved[key],key
    old=read(STAGE/'980/finite_thermal_records_results.json')
    p=np.array(saved['parameters']['thermal_populations']);e=np.array(saved['parameters']['energies'])
    assert saved['parameters']['thermal_populations']==old['parameters']['thermal_populations']
    assert saved['parameters']['energies']==old['parameters']['energies']
    q=float(p[2]+p[3]);C=float(p[0]+p[1]-p[2]-p[3]);g=saved['parameters']['g']
    assert abs(q-saved['peak_equal_prior_error'])<1e-14
    assert abs(C-saved['population_contrast'])<1e-14
    for row in saved['samples']:
        err=.5*(1-C*math.sin(g*row['time'])**2)
        assert abs(err-row['equal_prior_error'])<2e-12
    intervals=old['eigensystem_certificate']['coarse_analytic_certificate']['thermal_intervals']
    q_hi=sum(F(str(intervals[i][1])) for i in (2,3))
    assert q_hi<F(188903,1000000)
    bound=F(99,100)*F(188903,1000000)+F(1,200)+F(13,14000)+F(1,10000)
    assert bound==F(saved['certified_after_reset_error_bound']['exact'])<F(193043,1000000)
    aux=np.array(saved['parameters']['auxiliary_energies'])
    # Enumerate endpoint joint probabilities, independently of the H exponential.
    joint=np.zeros((4,5))
    for i in range(4):joint[3-i,i+1]=p[i]
    marginal_m=joint.sum(axis=1);marginal_b=joint.sum(axis=0)
    dM=float(e@(marginal_m-p));dB=float(aux@marginal_b-aux[0])
    assert abs(dM-saved['conditional_1_memory_energy_gain'])<1e-13
    assert abs(dB-saved['conditional_1_auxiliary_energy_change'])<1e-12
    assert abs(dM+dB)<1e-12
    def ent(arr):
        arr=np.asarray(arr);arr=arr[arr>0]
        return float(-np.dot(arr,np.log(arr)))
    mutual=ent(marginal_m)+ent(marginal_b)-ent(joint.ravel())
    assert abs(mutual-saved['conditional_1_memory_auxiliary_mutual_information'])<1e-12
    assert marginal_b[0]==0 and math.isclose(sum(marginal_b[1:]),1,abs_tol=1e-15)
    assert saved['unknown_source_dephased_locally']
    assert saved['auxiliary_initial_pure_energy_state_is_input']
    assert saved['engineered_energy_matched_interaction_is_input']
    note=STAGE/'research_note_1002.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for term in ('整体目标未完成','累计3785','不是持续全过程的QND读取','完整记录了来源标签','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'organization_lifecycle_adoption_v1.md',HERE/'thermal_record_reuse.py',
        HERE/'thermal_record_reuse_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1002.py',STAGE/'1003/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1002）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1002)==list(range(1,1003))
    assert '001—1002轮共1002份' in nav[0].read_text('utf-8-sig')
    assert '231—1002的772份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1002：热材料再写与功能分工' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1001/research_round_1001_checks.json',HERE/'drafts/STATUS.md',
              STAGE/'980/finite_thermal_records_results.json']
    prev=read(oldfiles[0])
    out=dict(round=1002,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1002,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_1001=prev['cumulative_numbered_test_groups_from_1000']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,old_thermal_state_reused=True,
        finite_window_rational_bound_independently_verified=True,
        auxiliary_energy_entropy_and_record_all_accounted=True,
        joint_clock_communication_lifecycle_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
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
        title='## 1002：热材料再写与功能分工';mark='## 1001：核组成与材料资源'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[1002报告]({pre}research_note_1002.md)将980实际热输出接到有限精度人口再写；'
            +'新增有限辅助体同时供能与记录，整窗叠旧误差后错率<0.193043。'
            +f'[结果]({pre}1002/thermal_record_reuse_results.json) · '
            +f'[核验]({pre}1002/research_round_1002_checks.json)。正式1002／累计3785，整体未完成。\n\n'
            +f'[机制补充]({pre}1002/organization_lifecycle_adoption_v1.md)按变量、时段与任务质量分工；'
            +'工程化耦合及纯辅助准备明示，未签收相位循环或全功能装置。停止器件优化，'
            +f'接[1003]({pre}1003/drafts/STATUS.md)整合整体假说新版与剩余机制缺口。应用目标及957保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—1001轮共1001份')==1
            s=s.replace('001—1001轮共1001份','001—1002轮共1002份',1)
        if p==paths[4]:
            before='当前正式1001／累计3784，1001已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式1002／累计3785，1002已结项',1)
        if p==paths[5]:
            assert s.count('231—1001的771份')==1
            s=s.replace('# 231—1001轮阶段成果总览','# 231—1002轮阶段成果总览',1)
            s=s.replace('231—1001的771份','231—1002的772份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至1001）';assert s.count(before)==1
            s=s.replace(before,before.replace('1001','1002'),1)
            rows={
                'C19':('|C19 真空、热态、非平衡准备|1、3、6、H2；645、955—980、994—1002|'
                    '1002同一980热态可作有限差错再写；辅助体能源／熵／记录共同入账|'
                    '纯辅助与工程化作用是新增准备；实际资源不能由比较热态冒充|'),
                'C20':('|C20 跨尺度、误差与共同来源|4、6、H3；854、898、899、936、950、957、971、981—1002|'
                    '1002按实际状态接续旧复位误差，固定材料谱／人口及整个读窗保留|'
                    '不能运输976相位合同；有限40维示例未完成SM／GR物理匹配|'),
                'C23':('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961—1002|'
                    '功能分工保留；1002证明非纯热材料仍可有限精度记录，来源相干与辅助资料保全局|'
                    '未证明免费纯准备、永久保持或循环辅助恢复；宇宙初始箭头未推得|')}
            for cid,row in rows.items():
                s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M)
                assert count==1,cid
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式1002／累计3785。热材料可按有限任务精度再用，'
                '辅助体的低熵准备、能源和资料不得遗漏。停止器件优化，下一项整合整体假说新版。'
                '完整共同生命周期及整体仍未完成；以下旧取舍保留当时范围。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1002补充](../../1002/organization_lifecycle_adoption_v1.md)'
                '把物质—资源链接到按任务质量划分的功能生命周期，旧热输出只接所声明人口记录。'
                '下一项汇总991—1002的整体解释与真实缺口，不继续修理本轮写入器。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with (HERE/'verify_round1002.py').open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 1002 thermal record reuse and lifecycle mechanism.')


if __name__=='__main__':main()
