"""Publish the bounded nuclear/material mechanism; preserve all prior rounds."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1001'


def main():
    old=(STAGE/'1000/verify_round1000.py').read_text('utf-8')
    prefix=old.split('    import importlib.util',1)[0]
    prefix=prefix.replace('dynamic thermoelastic interface 1000','nuclear formation interface 1001')
    prefix=prefix.replace('research_round_1000_checks.json','research_round_1001_checks.json')
    prefix=prefix.replace('range(776,1000)','range(776,1001)')
    checks=r'''
    import importlib.util
    from fractions import Fraction as F
    spec=importlib.util.spec_from_file_location('nuclear1001',HERE/'nuclear_formation_screen.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'nuclear_formation_results.json');core.compare(core.calculate(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1001
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('capture_rate_derived','cosmic_history_simulated','physical_abundances_predicted',
        'mass_uncertainty_propagated','work_extraction_certified','joint_record_lifecycle_certified',
        'cognitive_derivation_of_nuclear_force','full_goal_completed'):assert not saved[key],key
    masses={key:F(value['exact']) for key,value in saved['mass_inputs_MeV'].items()}
    mp,mn,md=[masses[key] for key in ('proton','neutron','deuteron')]
    W=mp+mn;B=W-md
    eg=F(saved['capture_photon_MeV']['exact']);kd=F(saved['deuteron_recoil_MeV']['exact'])
    threshold=F(saved['photodissociation_threshold_MeV']['exact'])
    assert eg+kd==B and (W-eg)**2-eg**2==md*md
    assert md*md+2*md*threshold==W*W
    assert F(saved['matched_energy_residual_MeV']['exact'])==0
    assert F(saved['double_count_residual_MeV']['exact'])==B
    assert math.isclose(math.hypot(float(md),float(B))+float(B)-float(W),
        saved['photon_equals_binding_on_shell_energy_excess_MeV'],rel_tol=2e-9)
    # Independently reconstruct densities and solve the equilibrium condition
    # by monotone bisection, not the generating program's quadratic formula.
    inp=saved['thermal_inputs'];a=inp['neutron_constituent_fraction'];b=inp['proton_constituent_fraction']
    for row in saved['equilibrium_rows']:
        T=row['temperature_MeV'];y=row['Yd']
        nb=inp['eta']*2*inp['zeta3']/math.pi**2*T**3
        K=.75*(2*math.pi/T)**1.5*float(md/(mp*mn))**1.5*math.exp(float(B)/T)
        assert math.isclose(nb*K,row['R'],rel_tol=2e-14)
        lo,hi=0.0,min(a,b)
        for _ in range(180):
            mid=(lo+hi)/2
            if mid-nb*K*(a-mid)*(b-mid)>0:hi=mid
            else:lo=mid
        assert math.isclose((lo+hi)/2,y,rel_tol=5e-13)
        assert math.isclose(row['Yn']+row['Yp']+2*y,1,abs_tol=2e-15)
        assert math.isclose(row['Yp']+y,b,abs_tol=2e-15)
        assert math.isclose(row['Yn']+y,a,abs_tol=2e-15)
        assert row['max_thermal_fugacity']<1e-8
    note=STAGE/'research_note_1001.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==10
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,6)]
    for term in ('整体目标未完成','累计3784','四份受限平衡状态','直接复用971','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'nuclear_material_resource_adoption_v1.md',HERE/'nuclear_formation_screen.py',
        HERE/'nuclear_formation_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1001.py',STAGE/'1002/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1001）',1)[1].split('### 当前取舍',1)[0]
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
    assert sorted(n for n in nums if n<=1001)==list(range(1,1002))
    assert '001—1001轮共1001份' in nav[0].read_text('utf-8-sig')
    assert '231—1001的771份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1001：核组成与材料资源' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'1000/research_round_1000_checks.json',HERE/'drafts/STATUS.md',STAGE/'research_note_971.md']
    prev=read(oldfiles[0])
    out=dict(round=1001,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1001,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_1000=prev['cumulative_numbered_test_groups_from_999']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,old_recoil_interface_reused=True,
        restricted_equilibrium_independently_solved=True,
        full_nuclear_rates_or_yields_verified=False,joint_record_lifecycle_certified=False,
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
        title='## 1001：核组成与材料资源';mark='## 1000：宏观相与热—力输运'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[1001报告]({pre}research_note_1001.md)接核束缚、反应组成、环境保存与有限资源；'
            +'反冲直接复用971，受限核平衡显示正结合能不足以保证环境中大量保存。'
            +f'[结果]({pre}1001/nuclear_formation_results.json) · '
            +f'[核验]({pre}1001/research_round_1001_checks.json)。正式1001／累计3784，整体未完成。\n\n'
            +f'[机制补充]({pre}1001/nuclear_material_resource_adoption_v1.md)区分存在、形成、保存与运行，'
            +'采用多来源元素与冷却机制；未预测现实丰度或提取功。停止核率／恒星细化，'
            +f'接[1002]({pre}1002/drafts/STATUS.md)补同一材料／环境中的功能分工与共同生命周期。应用目标及957保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—1000轮共1000份')==1
            s=s.replace('001—1000轮共1000份','001—1001轮共1001份',1)
        if p==paths[4]:
            before='当前正式1000／累计3783，1000已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式1001／累计3784，1001已结项',1)
        if p==paths[5]:
            assert s.count('231—1000的770份')==1
            s=s.replace('# 231—1000轮阶段成果总览','# 231—1001轮阶段成果总览',1)
            s=s.replace('231—1000的770份','231—1001的771份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至1000）';assert s.count(before)==1
            s=s.replace(before,before.replace('1000','1001'),1)
            rows={
                'C17':('|C17 质量、耦合、稳定物质结构|2、3、6；434、941—958、971、981、998—1001|'
                    '统计稳定与宏观相保留；1001补核组成、形成与环境保存，反冲复用971|'
                    '核质量／作用／状态为输入；未算完整核谱、材料相与现实丰度|'),
                'C19':('|C19 真空、热态、非平衡准备|1、3、6、H2；645、955—980、994—1001|'
                    '有限接收者、热—力账保持；1001区分平衡组成、反应时间及释能去向|'
                    '四份核平衡不是共同热史；真实速率、接热与可用功仍需对应过程|'),
                'C21':('|C21 完整物质与相互作用|1、4、5；531—553、735、939、951、957、971、981、1001|'
                    'P981输入保留；1001采用核有效作用／反应网络及多来源元素机制|'
                    '未从认知推出强力、核谱或核率；核界面筛选不替代SM全部匹配|'),
                'C23':('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961—1001|'
                    '相与输运保持；1001将组成历史痕迹、有限燃料及保存环境纳入解释|'
                    '材料痕迹不保证全部历史可读；共同记录生命周期与初始箭头未推得|')}
            for cid,row in rows.items():
                s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M)
                assert count==1,cid
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式1001／累计3784。核组成、环境保存与材料资源已补机制；'
                '核率和现实产额不冒领。下一项回同一环境中的功能分工和共同生命周期，停止核候选细化。'
                '整体仍未完成；以下旧取舍保留当时范围。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1001补充](../../1001/nuclear_material_resource_adoption_v1.md)'
                '把核组成、冷却与有限供给接入整体v1，复用971而不重证反冲。下一项补资源进入同一组织后的'
                '保持／读写／钟尺／通信／更新分工，不把完整器件制造作为前置门槛。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    target=HERE/'verify_round1001.py'
    with target.open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 1001 nuclear/material interface to eight navigation files.')


if __name__=='__main__':main()
