"""Publish the phase/transport mechanism and its bounded compatibility screen."""
from pathlib import Path
import re

STAGE=Path(__file__).resolve().parents[2]
RESEARCH=STAGE.parent
HERE=STAGE/'1000'


def main():
    previous=(STAGE/'999/verify_round999.py').read_text('utf-8')
    prefix=previous.split('    import importlib.util',1)[0]
    prefix=prefix.replace('Verify mechanism synthesis 999','Verify dynamic thermoelastic interface 1000')
    prefix=prefix.replace('research_round_999_checks.json','research_round_1000_checks.json')
    prefix=prefix.replace('range(776,999)','range(776,1000)')
    checks=r'''
    import importlib.util
    from fractions import Fraction as F
    import numpy as np
    spec=importlib.util.spec_from_file_location('core1000',HERE/'thermoelastic_common_account.py')
    core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
    saved=read(HERE/'thermoelastic_common_account_results.json');core.compare(core.run(),saved)
    assert saved['all_scientific_checks_passed'] and saved['round']==1000
    for rel,digest in saved['source_hashes'].items():assert sha(ROOT/rel)==digest,rel
    for key in ('long_time_ode_simulated','microscopic_phase_derived','transport_coefficients_derived',
                'noise_or_memory_lifetime_certified','all_discrete_nonuniform_modes_damped',
                'einstein_geometry_derived','quantum_channel_constructed','full_goal_completed'):
        assert not saved[key],key
    assert saved['discrete_checkerboard_null_mode_preserved']
    p={key:F(value['exact']) for key,value in saved['parameters'].items()}
    rho,K,b,c,T0,kappa=[p[key] for key in ('rho','K','b','c','T0','kappa')]
    state={key:[F(x['exact']) for x in value] for key,value in saved['finite_state'].items()}
    eps,v,T=state['strain'],state['velocity'],state['temperature']
    assert len(T)==4 and all(t>0 for t in T)
    # Independently differentiate E, S and A with the delivered actual rates.
    def exact(row,key):return F(row[key]['exact'])
    for key,expected in (('matched',F(0)),('omitted_thermal_feedback',F(1,25))):
        row=saved[key]
        de,dv,dt=[[F(x['exact']) for x in row[n]] for n in
                  ('strain_rate','velocity_rate','temperature_rate')]
        dE=sum(rho*vel*acc+(K*ep+b*T0)*dx+c*dtemp for ep,vel,acc,dx,dtemp in zip(eps,v,dv,de,dt))
        dS=sum(b*dx+c*dtemp/t for dx,dtemp,t in zip(de,dt,T))
        dA=sum(rho*vel*acc+K*ep*dx+c*(1-T0/t)*dtemp for ep,vel,acc,dx,dtemp,t in zip(eps,v,dv,de,dt,T))
        assert dE==expected==exact(row,'energy_rate')
        assert dS==F(1,12)==exact(row,'entropy_rate')
        assert dA==dE-T0*dS==exact(row,'availability_rate')
    edge=kappa*sum((T[i]-T[(i+1)%4])**2/(T[i]*T[(i+1)%4]) for i in range(4))
    assert edge==F(1,12)==F(saved['edge_entropy_production']['exact'])
    mode=saved['continuous_linear_symbol'];k=F(mode['k']['exact'])
    expected=[rho*c,rho*kappa*k**2,(K*c+b*b*T0)*k**2,K*kappa*k**4]
    assert [F(x['exact']) for x in mode['polynomial_coefficients']]==expected
    assert expected[1]*expected[2]-expected[0]*expected[3]==F(1,2)==F(mode['routh_margin']['exact'])
    for z in mode['roots']:
        root=complex(z['real'],z['imag'])
        assert root.real<0 and abs(np.polyval([float(x) for x in expected],root))<1e-12
    assert F(mode['isothermal_speed_squared']['exact'])==4
    assert F(mode['adiabatic_speed_squared']['exact'])==F(9,2)
    note=STAGE/'research_note_1000.md';prose=note.read_text('utf-8')
    assert prose.count('$$')==16
    assert re.findall(r'\\tag\{(\d+)\}',prose)==[str(i) for i in range(1,9)]
    for term in ('整体目标未完成','累计3783','棋盘格','未构造全部微观量子通道','386或425'):
        assert term in prose,term
    newfiles=[note,HERE/'material_phase_transport_adoption_v1.md',HERE/'thermoelastic_common_account.py',
        HERE/'thermoelastic_common_account_results.json',Path(__file__),HERE/'drafts/adoption_decision.md',
        HERE/'drafts/review_notes.md',HERE/'drafts/publish1000.py',STAGE/'1001/drafts/STATUS.md']
    for pth in newfiles:
        if pth.suffix=='.py':ast.parse(pth.read_text('utf-8'))
    nav=[STAGE.parent/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    mapping=nav[-1].read_text('utf-8-sig').split(
        '## 六条共同协议：全局缺口对应与检验优先级（截至1000）',1)[1].split('### 当前取舍',1)[0]
    assert re.findall(r'^\|(C\d\d) ',mapping,re.M)==[f'C{i:02d}' for i in range(1,28)]
    links=0
    for doc in [pth for pth in newfiles if pth.suffix=='.md']+nav:
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
    assert sorted(n for n in nums if n<=1000)==list(range(1,1001))
    assert '001—1000轮共1000份' in nav[0].read_text('utf-8-sig')
    assert '231—1000的770份' in nav[5].read_text('utf-8-sig')
    for pth in nav:assert '1000：宏观相与热—力输运' in pth.read_text('utf-8-sig')
    oldfiles=[STAGE/'999/research_round_999_checks.json',HERE/'drafts/STATUS.md']
    prev=read(oldfiles[0])
    out=dict(round=1000,date='2026-10-07',all_delivery_checks_passed=True,
        formal_reports=1000,fresh_test_groups=1,
        cumulative_numbered_test_groups_from_999=prev['cumulative_numbered_test_groups_from_998']+1,
        historical_unique_files_verified=len(frozen),historical_manifest_evidence=layout,
        local_links_checked=links,energy_entropy_availability_separately_verified=True,
        positive_entropy_alone_insufficient_in_stated_countermodel=True,
        microscopic_phase_derived=False,joint_record_lifecycle_certified=False,
        full_goal_completed=False,visual_checks_performed=False,app_goal_changed=False,
        frozen_inputs={str(pth.relative_to(ROOT)):sha(pth) for pth in oldfiles},
        new_scientific_and_entry_files={str(pth.relative_to(ROOT)):sha(pth) for pth in newfiles})
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
    target=HERE/'verify_round1000.py';assert not target.exists()
    paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[
        STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md',
                         '_shared/notes/unified_physics_condition_ledger_current.md')]
    originals={p:p.read_bytes() for p in paths};updates={}
    for p,raw in originals.items():
        s=raw.decode('utf-8-sig').replace('\r\n','\n')
        title='## 1000：宏观相与热—力输运';mark='## 999：选择性组织与整体假说'
        assert title not in s and s.count(mark)==1,p
        pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
        block=(title+'\n\n'
            +f'[1000报告]({pre}research_note_1000.md)复用645静态共源，补同一材料的集体运动与导热。'
            +'给定热弹性支保总能源、增热熵、降相对可用量；漏反向反馈即使熵产为正也会失去能源闭合。'
            +f'[结果]({pre}1000/thermoelastic_common_account_results.json) · '
            +f'[核验]({pre}1000/research_round_1000_checks.json)。正式1000／累计3783，整体未完成。\n\n'
            +f'[机制补充]({pre}1000/material_phase_transport_adoption_v1.md)采用相中慢变量及相容输运，'
            +'保必要记录与涨落；未生成微观材料相、真实系数或记忆寿命。停止该本构细化；'
            +f'接[1001]({pre}1001/drafts/STATUS.md)补核反应、元素组成与材料资源链。应用目标及957保持。\n\n')
        s=s.replace(mark,block+mark,1)
        if p==paths[0]:
            assert s.count('001—999轮共999份')==1
            s=s.replace('001—999轮共999份','001—1000轮共1000份',1)
        if p==paths[4]:
            before='当前正式999／累计3782，999已结项';assert s.count(before)==1
            s=s.replace(before,'当前正式1000／累计3783，1000已结项',1)
        if p==paths[5]:
            assert s.count('231—999的769份')==1
            s=s.replace('# 231—999轮阶段成果总览','# 231—1000轮阶段成果总览',1)
            s=s.replace('231—999的769份','231—1000的770份',1)
        if p==paths[-1]:
            before='## 六条共同协议：全局缺口对应与检验优先级（截至999）';assert s.count(before)==1
            s=s.replace(before,before.replace('999','1000'),1)
            rows={
                'C17':('|C17 质量、耦合、稳定物质结构|2、3、6；434、941—958、981、998—1000|'
                       '999统计稳定性保留；1000明示材料相与慢变量，给有限热—力共源过程|'
                       '实际相、成键、反应和寿命仍需输入／匹配；未从SM微观生成材料|'),
                'C19':('|C19 真空、热态、非平衡准备|1、3、6、H2；645、955—980、994—1000|'
                       '有限接收者保持；1000区分总能源、宏观熵与相对可用量，热—力交换共账|'
                       '正温局部热力学、材料相及输运为输入；非整体量子熵增长或已提取功|'),
                'C23':('|C23 统计、退相干与时间箭头|1、4、6；305、354、522—523、929、947、961—999、1000|'
                       '旧有限资源机制复用；1000将相中慢变量、耗散及涨落相容纳入解释|'
                       '确定性热例不认证噪声或记录寿命；宇宙初始边界与完整时间箭头未推得|'),
                'C26':('|C26 宏观经典、流体等|1、3、4、6；350、354—355、645、758—763、973—1000|'
                       '1000采用宏观相与集体输运机制，热弹性给共同能源／熵及漏反馈反例|'
                       '自由能不独自决定输运率；未完成全部相、电输运或微观匹配，声学不等于GR|')}
            for cid,row in rows.items():
                s,count=re.subn(r'^\|'+cid+r' [^\n]*$',lambda _:row,s,count=1,flags=re.M)
                assert count==1,cid
            anchor='### 当前取舍\n';assert s.count(anchor)==1
            s=s.replace(anchor,anchor+'\n正式1000／累计3783。已补相中慢变量、共同热—力输运及能源／熵区分；'
                '停止本构与器件细化，下一项核反应／元素组成／材料资源机制。整体仍未完成；'
                '以下旧取舍保留当时范围。\n',1)
            replacement=('4. 957数学假说v0.2保持；[1000补充](../../1000/material_phase_transport_adoption_v1.md)'
                '将整体v1接到宏观相、集体运动和输运，复用645而不重证静态共源。下一项补核束缚／反应'
                '如何改变组成和资源，不先计算全部核物理、恒星或宇宙过去。')
            s,count=re.subn(r'^4\. 957数学假说v0\.2保持；[^\n]*$',lambda _:replacement,s,count=1,flags=re.M)
            assert count==1
        if b'\r\n' in raw:s=s.replace('\n','\r\n')
        updates[p]=(b'\xef\xbb\xbf' if raw.startswith(b'\xef\xbb\xbf') else b'')+s.encode('utf-8')
    assert all(p.read_bytes()==raw for p,raw in originals.items()),'Concurrent navigation edit'
    with target.open('x',encoding='utf-8') as dest:dest.write(prefix+checks)
    for p,data in updates.items():p.write_bytes(data)
    print('Published 1000 phase and transport interface to eight navigation files.')


if __name__=='__main__':main()
