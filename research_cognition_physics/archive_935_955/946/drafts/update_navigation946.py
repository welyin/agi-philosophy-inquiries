"""One checked update to living navigation; frozen evidence stays unchanged."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 946：规范不变门户中的记录、径向物质与弱场引力'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 945：记录与Newton引力相位的共同有限过程'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[946报告]({pre}research_note_946.md)将规范不变Higgs门户的二次部门接入同一记录—Newton过程，三项对易效果有联合分布；有限质量下径向内容差>0.0034096、记录差>0.0212445、引力对照差>0.0216928。'
        +f'[结果]({pre}946/portal_common_process_results.json) · [核验]({pre}946/research_round_946_checks.json)。正式946／累计3731。'
        +'门户及物理表示是输入，完整非线性SM／Einstein匹配和六协议未签收。另限定945：路径0与相干路径不对易，旧表单项概率保持，不合称一次经典联合。'
        +f'停止门户细化，接[947]({pre}947/drafts/STATUS.md)按全局交付选择操作访问与能力运输。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—945轮共945份')==1
        s=s.replace('001—945轮共945份','001—946轮共946份')
        old='945在同一有限过程内接通记录与领先Newton引力相位，联合任务均有正下界；原E_rec及全部规范物质／六协议未因此并入。[覆盖表v0.2](archive_764_/945/drafts/common_model_scope_v0_2.md)按存在性重审验收，不先要求唯一参数、面积熵或UV完成。停止本例优化，接[946](archive_764_/946/drafts/STATUS.md)选择共同规范物质与操作接口。完整阶段尚未完成。'
        new='946以明确规范不变门户，将Higgs径向二次场与记录、Newton项接到同一有效H；三份对易效果有联合概率和有限质量正下界。945路径标签／相干读数的联合措辞已限定；原数值保持。[覆盖表v0.2](archive_764_/945/drafts/common_model_scope_v0_2.md)的存在性口径继续，完整SM／Einstein匹配与六协议仍未签收。停止门户优化，接[947](archive_764_/947/drafts/STATUS.md)选择共同操作访问与能力运输。完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        new='[946](archive_764_/research_note_946.md)已在同一记录—Newton有效H中加入规范不变Higgs门户的二次部门，三类有界任务均有正下界；并限定945的非对易路径读数表述。门户与物理表示为输入，完整非线性SM／Einstein及六协议未签收。沿[存在性验收范围](archive_764_/945/drafts/common_model_scope_v0_2.md)，不把唯一参数、面积熵或UV完成作为普遍门槛。停止门户细化，接[947](archive_764_/947/drafts/STATUS.md)回到共同操作访问与能力运输。'
        s,count=re.subn(r'^\[945\]\(archive_764_/research_note_945\.md\)已将[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1,p
    if p==paths[5]:
        assert s.count('# 231—945轮阶段成果总览')==1 and s.count('231—945的715份')==1
        s=s.replace('# 231—945轮阶段成果总览','# 231—946轮阶段成果总览').replace('231—945的715份','231—946的716份')
    if p==paths[-1]:
        old='六条共同协议：全局缺口对应与检验优先级（截至945）';assert s.count(old)==1
        s=s.replace(old,'六条共同协议：全局缺口对应与检验优先级（截至946）')
        rows={
          'C03':'|C03 事件身份、关系记录、访问|1、2、3、4；391、923、929、932—934、946|946给路径、关系记录、径向场三份对易效果的共同分布，并限定945非对易读数措辞|不同测量设置的合法概率不等于可合成一个经典联合；实际访问合同仍需实现|',
          'C14':'|C14 内部规范群、全局形式、连接|1、4、5；375、930—934、946|946的H†H门户保持选定SM规范群并接入实际记录来源|规范群是输入；本轮二次径向块不等于完整规范动力学或全部实验恢复|',
          'C17':'|C17 质量、Higgs、Yukawa等|1、3、6；925—927、931、941—946|946同一正Hessian给混合质量、记录传播和Higgs径向响应|门户、真空与参数为明确物理输入；非线性及真实参数匹配未认证，唯一性留后续|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、944—946|946重算两质量有限运动界，三项同一菜单任务有正下界|附加匹配误差小于0.001足以保本例正信号，但实际父理论误差尚未证明；不自动控制无界应力|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、939、943—946|记录内容可进入同一门户的Higgs径向二次部门；非零联合信号已核|有限准备、路径重合、径向读口仍为有效权限输入；完整六协议与物理材料访问待共同运输|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—946|946两场混合、噪声与交互能源共用同一H；总能量／动量守恒|父势及全部场应力应共同变分；完整非线性／曲背景有限匹配未签收，不继续候选高阶修补|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        new='4. 946使记录、Newton路径和Higgs径向二次部门共用一份过程，三项对易效果有正下界；停止门户及仪器优化。按[存在性范围](../../945/drafts/common_model_scope_v0_2.md)接[947](../../947/drafts/STATUS.md)优先选择共同操作访问与能力运输。完整物理匹配与六协议仍开放；一个候选的全部内部修复不是纲领前提。'
        s,count=re.subn(r'^4\. 945把领先Newton[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents for round 946.')
