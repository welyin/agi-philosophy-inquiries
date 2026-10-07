"""Update living navigation and the existence-stage acceptance boundary."""
from pathlib import Path
import re
STAGE=Path(__file__).resolve().parents[2];RESEARCH=STAGE.parent
paths=[RESEARCH/n for n in ('README.md','research_direction.md','RESEARCH_STATE.md')]+[STAGE/n for n in ('README.md','文件索引.md','阶段成果总览.md','跨阶段主题索引.md','_shared/notes/unified_physics_condition_ledger_current.md')]
original={p:p.read_bytes() for p in paths};updates={}
heading='## 945：记录与Newton引力相位的共同有限过程'
for p,raw in original.items():
    bom=raw.startswith(b'\xef\xbb\xbf');crlf=b'\r\n' in raw
    s=raw.decode('utf-8-sig').replace('\r\n','\n')
    oldheading='## 944：有限场通信与同一反作用的受控连接'
    assert heading not in s and s.count(oldheading)==1,p
    pre='archive_764_/' if p.parent==RESEARCH else '../../' if p==paths[-1] else ''
    block=(heading+'\n\n'+f'[945报告]({pre}research_note_945.md)将940的领先Newton约束势接入944同一量子过程，实际质量谱、记录和场噪声共用字典；有限质量下引力路径差>0.0220301、联合记录差>0.0155322。'
        +f'[结果]({pre}945/joint_record_newton_results.json) · [核验]({pre}945/research_round_945_checks.json)。正式945／累计3730。'
        +f'[范围表v0.2]({pre}945/drafts/common_model_scope_v0_2.md)区分已证连接与输入；完整应力、相对论／SM及六协议未签收。停止本例优化，'
        +f'接[946]({pre}946/drafts/STATUS.md)选择规范物质与操作接口。\n\n')
    s=s.replace(oldheading,block+oldheading,1)
    if p==paths[0]:
        assert s.count('001—944轮共944份')==1
        s=s.replace('001—944轮共944份','001—945轮共945份')
        old='939整体恢复表及940—943接口保持原范围；944在明确有效材料中连接有限记录、场噪声和机械交换，并给有限质量误差界。停止材料优化，接[945](archive_764_/945/drafts/STATUS.md)审计同一物理对象的覆盖与匹配，不能把不同模型各自成功相加签收。完整阶段尚未完成。'
        new='945在同一有限过程内接通记录与领先Newton引力相位，联合任务均有正下界；原E_rec及全部规范物质／六协议未因此并入。[覆盖表v0.2](archive_764_/945/drafts/common_model_scope_v0_2.md)按存在性重审验收，不先要求唯一参数、面积熵或UV完成。停止本例优化，接[946](archive_764_/946/drafts/STATUS.md)选择共同规范物质与操作接口。完整阶段尚未完成。'
        assert old in s;s=s.replace(old,new,1)
    if p in paths[1:3]:
        new='[945](archive_764_/research_note_945.md)已将领先Newton约束势与944记录置于同一H，使用实际内部质量；有限移动材料中两类有界任务均有正下界。[覆盖表v0.2](archive_764_/945/drafts/common_model_scope_v0_2.md)明确这尚非完整Einstein／SM／六协议共同实现；唯一参数选择、面积熵及UV完成不作为当前普遍门槛。停止软核、路径与器件优化，接[946](archive_764_/946/drafts/STATUS.md)由全局覆盖选择规范物质与操作接口。'
        s,count=re.subn(r'^\[944\]\(archive_764_/research_note_944\.md\)给重粒子[^ \n]*[^\n]*$',lambda m:new,s,count=1,flags=re.M);assert count==1,p
    if p==paths[5]:
        assert s.count('# 231—944轮阶段成果总览')==1 and s.count('231—944的714份')==1
        s=s.replace('# 231—944轮阶段成果总览','# 231—945轮阶段成果总览').replace('231—944的714份','231—945的715份')
    if p==paths[-1]:
        old='六条共同协议：全局缺口对应与检验优先级（截至944）';assert s.count(old)==1
        s=s.replace(old,'六条共同协议：全局缺口对应与检验优先级（截至945）')
        rows={
          'C10':'|C10 共同传播几何与普适耦合|1、4、5；334—340、924、930、940—945|945同一实际内部质量与位置进入Newton约束及量子记录过程|共同锥、等效耦合为物理输入；领先慢速恢复不等于完整普适应力匹配|',
          'C12':'|C12 Einstein项及有效修正|1、4、6；已有同作用接口、940、945|明示Einstein物理输入后有约束阶共同恢复，完整有效覆盖仍开放|当前先验同对象恢复及误差；选出几何系数的生成性要求留后续，不作存在性门槛|',
          'C13':'|C13 熵、面积项、引力耦合|1、4、6；313—315及热路线|部分同源系数已有；参考和平衡合同依路线而定|采用热力学引力路线或声称恢复相应任务时核验；不要求每份有效共同实现先推出面积熵|',
          'C17':'|C17 质量、Higgs、Yukawa等|1、3、6；925—927、931、941—945|同一质量可连接更新、阈值及领先引力，但原规范物质匹配未完成|可选相容质量机制和参数作为物理输入；当前验共同恢复，唯一数值与机制起源留后续|',
          'C20':'|C20 跨尺度、误差与共同来源|4、6、H3；854、898、944—945|945扩展有限质量界，同时运输记录与Newton路径有界任务|仍不自动控制无界应力或完整父理论匹配；只在声明任务域合并预算|',
          'C21':'|C21 主体、探针、记忆与通信材料|2、3、6、H2；853、939、943—945|同一有效过程内，关系记录与引力路径相位共同可辨|相干位置准备及重合读口为输入；完整六协议与规范材料未签收，不继续器件优化|',
          'C22':'|C22 应力、守恒与反作用|1、3、4、6、H2；935—945|945同一H保能量／动量并给Newton来源、力和相位；944场噪声完整保留|全部p场应力、后Newton和曲背景反馈未匹配；当前结果不冒充完整Einstein源|',
          'C27':'|C27 共同数据与可区分预测|1、3、5、6；现有同字典要求、945|945同准备省略Newton项产生有限可辨差，属于理论内部对照|核所声称覆盖的既有物理与误差；不要求当前已有独特新预测或拟合全部宇宙参数|'}
        for key,row in rows.items():
            s,count=re.subn(r'^\|'+key+r' [^\n]*$',lambda m:row,s,count=1,flags=re.M);assert count==1
        old='4. 944已取得有限场通信与机械反作用的共同连接，有限质量概率差有严格正下界。停止本有效材料优化；接[945](../../945/drafts/STATUS.md)以整体覆盖选择物理匹配。939局部形式分支、942材料工具与944有效H各有范围，未构成已完成的统一模型；任何单一候选的全部修复均非纲领普遍门槛。'
        new='4. 945把领先Newton约束与记录接入同一有效H，两类有界任务均有正下界；停止软核、路径和器件优化。按[覆盖表v0.2](../../945/drafts/common_model_scope_v0_2.md)接[946](../../946/drafts/STATUS.md)选择共同规范物质与操作接口。独立参数的唯一性、面积熵和UV完成不作当前普遍门槛；有限实际匹配仍需证明。'
        assert old in s;s=s.replace(old,new,1)
    if crlf:s=s.replace('\n','\r\n')
    updates[p]=(b'\xef\xbb\xbf' if bom else b'')+s.encode('utf-8')
assert all(p.read_bytes()==raw for p,raw in original.items())
for p,raw in updates.items():p.write_bytes(raw)
print('Updated eight living navigation documents and the stage acceptance boundary.')
