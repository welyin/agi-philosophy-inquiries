"""Finish the verified migration in the user's 15 sibling archive directories.

The source can be the untouched old path or the partially migrated intermediate
path. Every scientific byte is checked against the original ZIP before moving.
"""
import json
import os
from pathlib import Path
import re
import shutil
import zipfile
import organize_research_231_775 as m

ROOT,RESEARCH,OLD,INTERMEDIATE=m.ROOT,m.RESEARCH,m.OLD,m.NEW
PHASES=m.PHASES
LAST=RESEARCH/'archive_764_775'
FINAL_MIG=LAST/'_migration'


def phase_dir(p):
    return f"archive_{p['start']}_{p['end']}"


def atomic_write(p,raw):
    temporary=p.with_name(p.name+'.migration_tmp')
    assert not temporary.exists()
    p.parent.mkdir(parents=True,exist_ok=True)
    with temporary.open('xb') as f:
        f.write(raw)
    os.replace(temporary,p)


def write(p,text):
    assert not p.exists(),p
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(text,encoding='utf8',newline='\n')


def rel(p,base):
    return m.relative(p,base)


def main():
    source_mig=INTERMEDIATE/'_migration'
    oldplan=json.loads((source_mig/'plan.json').read_text('utf8'))
    snapshot=json.loads((source_mig/'snapshot_checks.json').read_text('utf8'))
    archive=source_mig/snapshot['file']
    assert m.digest(archive)==snapshot['sha256']
    assert not (FINAL_MIG/'manifest.json').exists()
    entries=[]
    prefix={m.folder(p):phase_dir(p) for p in PHASES}
    for e in oldplan['entries']:
        parts=Path(e['destination']).parts
        if parts[0] in prefix:
            destination=Path(prefix[parts[0]],*parts[1:])
        else:
            destination=Path('archive_764_775',*parts)
        entries.append(dict(**{k:v for k,v in e.items() if k!='destination'},
                            destination=destination.as_posix(),intermediate=e['destination']))
    mapping={(OLD/e['original']).resolve():(RESEARCH/e['destination']).resolve() for e in entries}
    assert len(set(mapping.values()))==8920
    assert all(p.is_relative_to(RESEARCH.resolve()) for p in mapping.values())
    # A link to the former flat directory should reach the new research index.
    m.NEW=RESEARCH
    records=[]
    with zipfile.ZipFile(archive) as z:
        for i,e in enumerate(entries):
            raw=z.read('research_cognition_physics/archive_231_/'+e['original'])
            assert m.sha(raw)==e['original_sha256']
            old=OLD/e['original']
            intermediate=INTERMEDIATE/e['intermediate']
            existing=[p for p in (old,intermediate) if p.is_file()]
            assert len(existing)==1,(e['original'],existing)
            src=existing[0]
            current=src.read_bytes()
            if e['rewrite_markdown_links']:
                # Intermediate edits were link-only; the final text is regenerated
                # directly from the verified original, never by editing edited prose.
                enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
                def without_targets(data):
                    text=data.decode(enc)
                    for a,b,_,_ in reversed(list(m.links(text))):
                        text=text[:a]+'<LOCAL_LINK>'+text[b:]
                    return text
                assert without_targets(current)==without_targets(raw),e['original']
            else:
                assert current==raw,e['original']
            dst=RESEARCH/e['destination']
            assert src.resolve().is_relative_to(RESEARCH.resolve())
            assert dst.resolve().is_relative_to(RESEARCH.resolve()) and not dst.exists()
            changed,edits=m.rewrite_links(raw,old,dst,mapping) if e['rewrite_markdown_links'] else (raw,[])
            dst.parent.mkdir(parents=True,exist_ok=True)
            src.replace(dst)
            if changed!=current:
                atomic_write(dst,changed)
            records.append(dict(**e,current_sha256=m.sha(changed),current_bytes=len(changed),link_edits=edits))
            if (i+1)%2000==0:
                print(f'Final layout: {i+1}/8920 files',flush=True)
    FINAL_MIG.mkdir(parents=True,exist_ok=True)
    for p in list(source_mig.iterdir()):
        assert p.is_file()
        name='plan_intermediate.json' if p.name=='plan.json' else p.name
        target=FINAL_MIG/name
        assert not target.exists() and target.resolve().is_relative_to(RESEARCH.resolve())
        p.replace(target)
    for base in (OLD,INTERMEDIATE):
        for d in sorted((p for p in base.rglob('*') if p.is_dir()),key=lambda p:len(p.parts),reverse=True):
            assert d.resolve().is_relative_to(base.resolve())
            d.rmdir()
        assert not list(base.iterdir())
        base.rmdir()
    document(records,mapping,snapshot)
    manifest=dict(version=2,date='2026-10-04',rounds=[231,775],phases=PHASES,
                  destination_base='research_cognition_physics',original_file_count=len(records),
                  entries=records,scientific_code_and_json_unchanged=True,
                  mathematical_prose_unchanged=True,goal_unchanged_and_paused=True,
                  new_scientific_rounds=0,user_selected_sibling_archive_layout=True)
    m.dump(FINAL_MIG/'manifest.json',manifest)
    print('Finished 15 sibling archives; goal remains paused.',flush=True)


def document(records,mapping,snapshot):
    def note(n,base):
        return rel(mapping[OLD/f'research_note_{n}.md'],base)
    rows=[]
    for phase in PHASES:
        base=RESEARCH/phase_dir(phase)
        text=f"# {phase['start']}—{phase['end']}轮：{phase['slug']}\n\n[返回研究总目录](../README.md) · [文件索引](文件索引.md) · [复算说明](../archive_764_775/_migration/README.md)\n\n"
        text+=f"## 研究问题\n\n{phase['question']}\n\n## 阶段成果\n\n"+'\n'.join('- '+x for x in phase['results'])+'\n\n'
        text+=f"## 适用范围与未完成项\n\n{phase['limits']}\n\n这是研究过程的主题分期，不代表本阶段全部开放问题已关闭。认知动机、附加建模输入、解析证明和数值校准按原报告区分。\n\n## 关键报告\n\n"
        text+='、'.join(f'[{n}]({note(n,base)})' for n in phase['keys'])+'\n\n'
        text+=f"## 与后续阶段的连接\n\n{phase['transition']}\n\n## 全部轮次\n\n|轮次|报告|\n|---|---|\n"
        for n in range(phase['start'],phase['end']+1):
            title=mapping[OLD/f'research_note_{n}.md'].read_text('utf-8-sig').splitlines()[0].lstrip('# ').replace('|','／')
            text+=f'|{n}|[{title}]({note(n,base)})|\n'
        write(base/'README.md',text)
        index='# 文件索引\n\n[阶段README](README.md)。原文件名保留；跨阶段脚本使用[统一复算入口](../archive_764_775/_migration/README.md)。历史核验中的路径/哈希按原快照解释。\n\n|类型|文件|原路径|\n|---|---|---|\n'
        for e in records:
            if e['destination'].startswith(phase_dir(phase)+'/'):
                index+=f"|{e['kind']}|[{Path(e['destination']).name}]({rel(RESEARCH/e['destination'],base)})|`{e['original']}`|\n"
        write(base/'文件索引.md',index)
        rows.append(f"|{phase['id']}|{phase['start']}—{phase['end']}|[{phase['slug']}]({phase_dir(phase)}/README.md)|")
    overview='''# 231—775轮阶段成果总览

## 当前阶段定位

**已经有阶段性成果，完整统一目标尚未完成。** 545份编号研究报告形成三类主要资产：空间和操作之间的条件性定理、限制特定候选的严格反例，以及共同承载物质、参考、记录和几何反作用的联合模型构造。

001—230的有限维量子结论继续作为基础。本次按问题转折整理231—775，不把报告数量、数学框架的相似或有限数值结果当成全部物理已被推导的证据。

## 一、可以复用的空间定理

382—386建立完整反向接口的三维上界、重定向下界及实际位移邻域与有限尺度证书。384已经消去旧下界中的额外Lipschitz条件；425给一致半幅与成本收缩到光滑坐标的条件连接。522—523将热参考与实际方向仪器放到同一实现中。

这些结论的价值在于指出何种实际接口足够。仍未证明FUCP或任意认知过程都必须满足它们的所有前提；qubit内部的三个分量不自动等于三维物理空间。

## 二、真正排除候选的边界

342—344给既有操作、传播和共同几何条件下的不唯一性。604排除指定中心差分物种接法；649排除原独立几何动能直接下降到指定面积匹配商；699严格排除声明的辅助候选对全部时间参数具有物理正性；759排除对全部输入统一成立的单确定背景轨道合同。

各反例的前提和量词分别保留。它们限制后续设计，但不是对所有认知模型、所有经典近似或所有量子引力路线的反证。

## 三、同一材料共同承担多个角色

531以后联合表示、反常、质量、真空、尺度与几何。548—573及753给同一物质提供关系参考并与动态经典几何相容的构造；574、598及后续连接有限图量子过程。625、702—704、723等将实际记录、热态和来源放入同一近似族。

这比把不同模型的结果并列更强，但四维、群/物种、经典作用和部分参数仍为输入，不能写成标准模型或Einstein作用已经由认知唯一生成。

## 四、最新共同量子构造

730—741接通声明费米部门、实际来源与有限阶反作用。764—767在原共同背景上建立完整线性物理代数及正Hadamard态；768—772保持同一物理态并修复完整首阶局部来源。773—775连接局部形式修复、共同插入与实际来源首项。

775关闭首项相容接口。高次共同规范化N1/N2、实际相互作用准备和记录、严格UV/全局拼接、原图到连续的共同映射仍开放。阶段性归档并不等于完整量子引力结项。

## 五、分期与恢复

15个阶段按研究对象和方法转折划分，原文件各归属一次，用主题索引表达跨阶段复用。目录整理不增加研究轮次。目标保持暂停；用户恢复后从已保存的776入口继续，先核高次共同规范化，不重做773—775首项或旧空间证明。

本总览概括已存报告的证据等级，不冒充逐行独立审稿。本次逐项核对迁移及代表性复算，没有重跑全部545轮或做图像检验。

[15阶段目录](../README.md) · [跨阶段主题](跨阶段主题索引.md) · [776入口](next_round_776/STATUS.md) · [迁移核验](_migration/README.md)
'''
    write(LAST/'阶段成果总览.md',overview)
    topics=[('空间与坐标',[372,382,383,384,386,425,522,523]),('引力前提与非唯一性',[301,304,324,342,344,351,358]),('内生交互与传播',[429,451,459,465,466,504]),('共同物质与几何',[531,548,549,551,572,573,753]),('记录、热态与来源',[574,598,625,702,704,723,730,735,741]),('限定失败',[604,649,699,759]),('关联与当前量子入口',[756,758,762,763,767,771,773,774,775])]
    text='# 跨阶段主题索引\n\n[研究总目录](../README.md)。原文件只归属一个阶段；此处按复用问题连接。\n\n'
    for title,ns in topics:
        text+='## '+title+'\n\n'+'、'.join(f'[{n}]({note(n,LAST)})' for n in ns)+'\n\n'
    write(LAST/'跨阶段主题索引.md',text)
    write(LAST/'_shared/README.md','# 跨阶段公共材料\n\n没有唯一编号归属的专题和辅助材料集中保留于此，不构成新增研究阶段。[研究总目录](../../README.md) · [迁移和复算](../_migration/README.md)。完整文件列表见[末阶段文件索引](../文件索引.md)。\n')
    write(LAST/'_history/README.md','# 历史导航与快照\n\n原长篇索引在`indexes`中，链接已经迁移；历次导航快照在`navigation_snapshots`中，字节完全不改。快照中的旧相对路径是历史证据，不能作为当前导航；原布局由[复算入口](../_migration/README.md)恢复。当前以[研究总目录](../../README.md)和[状态](../../RESEARCH_STATE.md)为准。\n')
    migration_readme='''# 迁移与复算说明

原archive_231_的8920份文件已实际移动。545份正式报告按15个并列archive目录归属；每阶段包含notes、code、results和必要的drafts。跨阶段公共材料及历史导航集中保存在archive_764_775的_shared与_history中，不表示它们全部产生于764轮以后。

## 冻结证据

科学代码、JSON结果和历史核验保持原字节；Markdown阅读版只更改链接目标，数学及文字正文不改。[迁移映射](manifest.json)记录每个旧/新路径、哈希及精确链接改动。

[原工作区快照](original_workspace_231_775.zip)保存迁移前全部原件及复算上下文。[快照核验](snapshot_checks.json)已逐文件验证。历史核验中的相对路径和哈希指原件，不应与链接迁移后的阅读版混用。

## 统一复算入口

在本目录用既有Python运行时执行：

```powershell
& 'C:\\Users\\admin\\.cache\\codex-runtimes\\codex-primary-runtime\\dependencies\\python\\python.exe' -B -X utf8 .\\replay.py --script verify_round775.py
```

脚本名可用旧相对路径，也可用manifest中的新相对路径；额外参数放在`--`后。运行器验证当前科学代码与原件一致，在隔离临时目录恢复原布局并调用同一Python/NumPy，结束后只清理自己创建的临时目录。旧脚本直接在各code目录运行仍可能缺少跨阶段依赖，因此统一使用此入口；没有为改目录而改1858份Python算法。

## 核验与适用范围

[迁移核验](migration_checks.json)检查全部8920文件、允许的链接变化、15个目录和545份报告。[代表性复算](replay_checks.json)逐项列实际运行结果。未重跑全部545轮，未重审全部证明，未增加科学计数，也未做图像检查。原本已有的历史断链单独记录，不隐去。

[研究总目录](../../README.md) · [阶段成果](../阶段成果总览.md)
'''
    write(FINAL_MIG/'README.md',migration_readme)
    # Preserve external prior bytes in the already verified workspace ZIP.
    nav_changes=[]
    for p in [ROOT/'README.md',RESEARCH/'archive_223_230/README.md']:
        original=p.read_bytes()
        new,edits=m.rewrite_links(original,p,p,mapping)
        if p==ROOT/'README.md':
            txt=new.decode('utf-8-sig').replace('\r\n','\n')
            head,rest=txt.split('\n\n',1)
            txt=head+'\n\n**231—775轮已按15个并列archive目录归档，研究目标暂停。** [阶段目录与成果](research_cognition_physics/README.md)。报告、代码和结果已实际迁移；恢复研究从776入口继续。下方保留历史进展。\n\n'+rest
            new=txt.encode('utf8')
        atomic_write(p,new)
        nav_changes.append(dict(path=p.relative_to(ROOT).as_posix(),before_sha256=m.sha(original),after_sha256=m.sha(new),link_edits=edits))
    readme='''# 认知物理研究：阶段论文与研究档案

## 当前状态

目标由用户暂停。231—775轮共545份编号报告，已实际分拆到以下15个并列目录。已有条件性空间定理、限定反例和联合模型构造；完整统一目标尚未完成。

[阶段成果总览](archive_764_775/阶段成果总览.md) · [跨阶段主题索引](archive_764_775/跨阶段主题索引.md) · [当前状态](RESEARCH_STATE.md) · [研究方向](research_direction.md)

## 231—775阶段目录

|阶段|轮次|主题|
|---|---|---|
'''
    readme+='\n'.join(rows)+'''

## 前两阶段论文

- [可组合认知结构与复量子状态空间](可组合认知结构与复量子状态空间_阶段论文.md)：[001—222档案](archive_001_222/README.md)。
- [可组合认知结构与有限维量子理论](可组合认知结构与有限维量子理论_阶段论文.md)：[223—230档案](archive_223_230/README.md)。

## 使用方式

每个archive的README说明问题、成果、假设边界、关键报告和全部轮次；原始报告、代码及结果在其子目录中。旧稿和历史快照全部保留，迁移不新增科学轮次。[复算与原始快照](archive_764_775/_migration/README.md)保留旧固定路径的可执行环境。

恢复研究从[776已有入口](archive_764_775/next_round_776/STATUS.md)开始；本次整理不会自动恢复目标。
'''
    direction='''# 研究方向：认知本体论与现代物理的共同模型

## 目标保持

本次只做231—775分期归档，不改写应用中的目标。继续寻找同一结构在明确范围内共同承载量子理论、3+1时空、广义相对论、标准模型及测量记录。允许认知到物理与物理到认知的双向构建；认知观察提出假说，数学检查连接，物理预测检验后果。

## 已有材料

[阶段总目录](README.md)、[成果总览](archive_764_775/阶段成果总览.md)与[主题索引](archive_764_775/跨阶段主题索引.md)分别组织历史与跨阶段复用。原全部方向记录已保存在迁移快照，不把历史“下一步”当作当前指令。

## 恢复后的顺序

1. 从[776入口](archive_764_775/next_round_776/STATUS.md)核原混合体系的高次共同因果规范化N1/N2。
2. 复用773—775的局部修复与首项匹配，不再以完整in-in有效作用作为该首项的必要前置。
3. 接实际相互作用准备与记录、原图到连续和共同跨尺度映射，压缩独立输入并形成可区分后果。

空间接口回用382—386、425及522—523；384已消去的额外Lipschitz条件不得重新列为缺口。604、649、699、759的限定失败按原量词保留。

## 执行约束

目标当前暂停，待用户恢复。新科学轮次仍交付编号报告、推导/反例及必要的可复算代码结果，整理不新增编号。新增研究另建工作目录，不覆盖本次归档。无新应用任务、定时任务或图像检查。区分认知动机、附加输入、成熟定理、解析证明、数值验证与物理解释；预置Einstein作用的模型验证不是独立生成引力。

历史复算使用[统一入口](archive_764_775/_migration/README.md)。
'''
    state='''# 研究状态

## 2026-10-04：目标暂停，231—775分期归档

- 最新完成：775；下一编号：776，未开始新科学轮次。
- 231—775共545份正式报告，实际分拆至15个并列archive目录。
- 最新历史核验记载累计3510项科学检查、1637份编号科学文件、3829份保护证据；迁移检查另计，不代表重跑全部历史实验。
- [阶段目录](README.md)、[阶段成果](archive_764_775/阶段成果总览.md)、[主题索引](archive_764_775/跨阶段主题索引.md)。

## 已接通

原共同背景的完整线性物理代数与正Hadamard态；同一自由态完整首阶来源的局部Ward修复；775进一步直接连接该来源与局部BV反常首项。条件性三维、光滑坐标和先前同源记录成果继续复用。

## 尚未接通

高次共同规范化N1/N2；实际相互作用准备、正性和记录；严格UV/全局拼接；原有限图到连续模型的共同映射。维数、群/物种、经典作用及若干参数仍按输入记账。

## 恢复入口

[776冻结入口](archive_764_775/next_round_776/STATUS.md)。用户恢复目标前，只完成本次整理，不自动开启776。新增科学工作应另建目录，引用新路径，不恢复向旧archive_231_写入。

## 证据

所有历史文件与迁移前导航已保存原字节快照。日常阅读版仅调整链接，代码/结果不变；通过[统一复算入口](archive_764_775/_migration/README.md)重建旧布局。旧状态中的历史安排不覆盖此页。
'''
    for name,text in [('README.md',readme),('research_direction.md',direction),('RESEARCH_STATE.md',state)]:
        p=RESEARCH/name
        original=p.read_bytes()
        atomic_write(p,text.encode('utf8'))
        nav_changes.append(dict(path=p.relative_to(ROOT).as_posix(),before_sha256=m.sha(original),after_sha256=m.digest(p),navigation_compacted=True))
    m.dump(FINAL_MIG/'navigation_changes.json',nav_changes)


if __name__=='__main__':
    main()
