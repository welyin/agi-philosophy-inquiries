from pathlib import Path
import hashlib,json,re,sys

ROOT=Path.cwd()
HERE=ROOT/'research_cognition_physics/archive_1106_/_shared/axiom_wording_revision'
sys.path.insert(0,str(ROOT/'scripts'))
from organize_research_1044_1062 import text_edits
from research_layout import relocation_original_bytes

def sha(b):return hashlib.sha256(b).hexdigest()
originals={}
updates={}
def get(name):
    raw=(ROOT/name).read_bytes()
    originals[name]=raw
    return raw.decode('utf-8')
def put(name,text):updates[name]=text.encode('utf-8')
def replace(text,before,after):
    assert text.count(before)==1,(before[:100],text.count(before))
    return text.replace(before,after,1)
def math(text):return re.findall(r'\$\$.*?\$\$|(?<!\$)\$(?!\$)[^$\r\n]*\$(?!\$)',text,re.S)

layer_name='research_cognition_physics/_migration/motion_active_1106/layer.json'
layer=json.loads(get(layer_name))
old_entries={e['original']:e for e in layer['entries']}
nav_names=['research_cognition_physics/'+x for x in ['README.md','research_direction.md','RESEARCH_STATE.md','ROADMAP.md']]
paper_name='research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md'

# Only the live heading blocks are revised. Everything after the historical marker stays byte-identical.
for name in nav_names:
    text=get(name)
    marker='## 历史导航（截至1105，保留原研究时点）'
    front,tail=text.split(marker,1)
    front=replace(front,'1086—1105已封存，阶段论文及历史研究不改。','1086—1105洛伦兹运动学阶段已经结项；历史研究保持原文，当前论文的体系说明随公理整理维护。')
    front=replace(front,'继承F＋U＋C＋P四原则、CO1—CO5认知原子操作五公理与已采用的SR-T／A／R／O／F；SR-S只在明确启用的可选定标分支使用，生成完备G仍未采用。继续区别原子公理、组合权限与宏观读出，不将新输入追认为旧原则的推论。','**当前认知公理体系由F＋U＋C＋P四原则、CO1—CO5认知原子操作五公理和SR-T／A／R／O／F时空操作五公理共同组成。** 四原则承担量子系统与组合的基础职责，CO规定原子关系操作，SR规定完整时空任务及参考结构；它们分层共同使用，SR不是单原子门的列表。体系随研究逐阶段扩充，不预设最初四原则已经完备，后续持续整理与收敛。SR-S仍为明示启用的可选定标扩展，生成完备G未采用。')
    front=replace(front,'旧四原则与CO五公理继续共同使用，G未采用。','四原则、CO与SR共同构成当前认知公理体系，G未采用。')
    front=replace(front,'新增条款作为理论输入采用，不追认为旧认知原则的推论。实际钟、全部影响、物质交互以及新旧全部条件的共同实现仍待后继完成。','这些是当前认知公理体系内已经完成的运动学结论。下一阶段继续连接实际钟、全部纳入研究的影响及物质交互，并核验整个体系在同一模型中的共同实现。')
    put(name,front+marker+tail)

text=get(paper_name)
paper_before=text
pairs=[
('**采用状态。** 依据本次明确授权，本文正式采用SR-T、SR-A、SR-R、SR-O、SR-F五项时空操作公理，并正式定义加入SR-S的可选定标扩展版本。基础五项不包含SR-S。既有认知四原则F＋U＋C＋P及认知原子操作五公理CO1—CO5继续作为共同基础；本文没有把新增条款追认为旧原则的推论，也未采用全宇宙生成完备猜想G。历史报告中的“候选”“未采用”是当时状态，保留原文。',
 '**体系定位。** 洛伦兹运动学阶段已经完成。当前认知公理体系由F＋U＋C＋P四原则、CO1—CO5认知原子操作五公理，以及本文正式采用的SR-T、SR-A、SR-R、SR-O、SR-F时空操作五公理共同组成。SR是该体系的正式组成部分，SR-S为明示选用的定标扩展，基础版本不包含SR-S；全宇宙生成完备猜想G仍未采用。认知公理体系随研究逐阶段扩充，不预设最初四原则已经完备，后续将持续检查相容性、整理依赖并压缩冗余。公理的采用与定理的证明分别登记；历史报告中的“候选”“未采用”保留其原研究时点。'),
('哪些附加操作条件能够确定惯性参考之间的运动学。','认知时空操作公理怎样确定惯性参考之间的运动学。'),
('全域平直字典、实际参考权限及跨全部准备的期限边界是明示输入。','全域平直字典、实际参考权限及跨全部准备的期限边界是本体系时空公理的明确内容。'),
('以及新增公理和旧基础的完整联合实现，仍是后续证明责任。','以及当前认知公理体系的完整共同实现，仍是后续证明责任。'),
('本文SR条款约束完整任务及其事件描述；其中SR-A明确属于平直运动学的描述条件。它们不是五条单原子门的列表。将来若给出其原子来源，须证明原子协变、组合继承和实际读口相容，不能只将宏观条款改名为“原子性质”。',
 '当前认知公理体系具有分层职责：四原则规定量子系统与组合的基础结构，CO五公理规定基础原子关系操作，SR五公理规定完整时空任务、事件描述及惯性参考。其中SR-A承担全域平直运动学的描述职责；SR整体不是五条单原子门的列表。后续整理可研究各层的派生关系及更精简的共同基础；若将某条SR公理进一步联系到原子机制，须给出原子协变、组合继承和实际读口相容的连接。'),
('给出旧认知基础与具体时空运动学之间的输入账。','说明当前认知公理体系中各层公理对时空运动学定理承担的职责。'),
('该强量词现作为新增假设正式承担，不再伪装成“认知有限”的直接推论。','该跨全部准备的共同边界由SR-F明确规定；固定组织的能力边界与这一共同边界在体系中分别登记。'),
('本文不把采用新公理写成已经验证全部联合相容性。','整个认知公理体系在同一物理模型中的联合相容性，仍须按声明的任务和有效域继续验证。'),
('这给SR-R的原子—组合来源路线，不要求不同任意宏观装置能力相同。但原子协变、任务生成覆盖、实际读口及准备族闭合仍要证明或明确新增为条件，(25)本身不会无条件生成它们。',
 '这给正式公理SR-R与原子—组合机制之间的一条连接路线，不要求任意两个宏观装置能力相同。SR-R已经属于当前体系；进一步通过(25)压缩其独立地位时，须验证原子协变、任务生成覆盖、实际读口和准备族闭合，逐项记录由哪些公理或机制承担。'),
('在同一模型同时满足新旧基础和真实准备资格','在同一模型同时满足当前认知公理体系和真实准备资格'),
('也不把预设Lorentz物质模型的复现登记成其认知来源证明。','后续将分别登记公理内的演绎结论、具体机制的建模与验证，以及尚待进一步派生的结构。'),
('四原则、五条原子操作公理及新增时空公理的完整共同实现尚待验证；全部实际钟、信号、物质与相互作用的狭义相对论连接尚未结项。正式采用公理改变了研究的出发条件，不将待证来源追认为既有定理。原强目标仍未完成，本次收尾范围是1086—1105的公理化惯性运动学阶段。',
 '四原则、CO五公理与SR五公理共同构成当前认知公理体系，洛伦兹运动学已在该体系内完成推导。本次收尾范围是1086—1105的公理化惯性运动学阶段。下一阶段继续核验该体系的完整共同实现，以及实际钟、全部纳入研究的信号、物质与相互作用之间的狭义相对论连接；原较强目标中的这些连接尚未完成。研究将沿扩充、检验、整理与收敛的路径继续推进。')]
for a,b in pairs:text=replace(text,a,b)
assert math(text)==math(paper_before),'Formula content changed'
a=paper_before.index('## 4. 有限任务速度选择定理');b=paper_before.index('## 8. 相容性见证')
c=text.index('## 4. 有限任务速度选择定理');d=text.index('## 8. 相容性见证')
assert paper_before[a:b]==text[c:d],'Proof sections changed'
put(paper_name,text)

name='research_cognition_physics/archive_1106_/README.md'
text=get(name)
text=replace(text,'继承四原则F＋U＋C＋P、CO1—CO5与正式采用的SR-T/A/R/O/F。SR-S为明示采用的可选定标分支；生成完备G仍未采用。已证明的完整Lorentz子群不重列为缺口；实际钟律、所有影响的共同锥、材料与交互并未自动由它给出。',
 '本阶段以**当前认知公理体系**为共同基础：F＋U＋C＋P四原则、CO1—CO5原子关系操作五公理，以及SR-T/A/R/O/F时空操作五公理。四原则不被预设为已经完备的最终清单；体系将随研究扩充，并在各阶段整理与收敛。三组公理承担量子结构、原子关系操作和完整时空任务的不同职责，SR不被缩成单原子门。SR-S为明示启用的可选定标扩展，G未采用。洛伦兹运动学阶段已完成；本阶段继续实际钟、全部纳入研究的影响、材料与交互的共同连接。')
text += '\r\n[认知公理体系措辞与兼容性修订](_shared/axiom_wording_revision/README.md)。\r\n'
put(name,text)
name='research_cognition_physics/archive_1106_/_shared/goal.md'
text=get(name)
goal_quote_before=text.split('## 目标原文',1)[1].split('## 起点与验收',1)[0]
text=replace(text,'## 起点与验收','## 当前体系说明\r\n\r\n上面的目标原文及工具回执按启动时内容保存。当前研究统一使用“认知公理体系”：四原则、CO原子操作五公理及SR时空操作五公理均为其正式组成部分；SR-S仍是可选定标扩展，G未采用。体系逐阶段扩充并整理收敛，不预设最初四原则已经穷尽认知的公理结构。这里调整公理身份的说明，不改工具目标、数学条款或阶段验收范围。\r\n\r\n## 起点与验收')
text=replace(text,'1086—1105阶段已结项。继承[时空阶段论文](../../认知时空操作公理与洛伦兹运动学_阶段论文.md)中正式采用的五项基础SR公理，以及四原则、CO五公理；不改其量词和适用范围。SR-S仅属于可选定标扩展，G未采用。公理采用不是证明它们来自旧原则，也不是证明现实满足它们。',
 '1086—1105洛伦兹运动学阶段已结项。本阶段承接[时空阶段论文](../../认知时空操作公理与洛伦兹运动学_阶段论文.md)及当前完整认知公理体系，保留各项公理的量词、层次和适用范围。SR-S仅属于可选定标扩展，G未采用。已完成的运动学定理直接复用；共同物理实现与现实预测按本阶段任务继续核验。')
assert goal_quote_before==text.split('## 目标原文',1)[1].split('## 当前体系说明',1)[0]
put(name,text)
name='research_cognition_physics/archive_1106_/_shared/dependency_ledger.json'
ledger=json.loads(get(name))
ledger['foundation']['name']='current_extensible_cognitive_axiom_system'
ledger['foundation']['F_U_C_P']='正式组成部分：量子系统与组合的基础原则；不预设该初始四原则已构成最终完备公理集'
ledger['foundation']['CO1_to_CO5']='正式组成部分：基础原子关系操作及其组合，保留原量词与适用域'
ledger['foundation']['SR_T_A_R_O_F']='正式组成部分：完整时空任务、事件字典及惯性参考结构；洛伦兹运动学已在当前体系内推导'
ledger['foundation']['SR_S']='明确启用时采用的可选定标扩展'
ledger['foundation']['development']='逐阶段扩充、核验相容性、整理依赖并收敛；公理采用与演绎证明分开登记'
items=ledger.pop('inherited_adopted_inputs')
for item in items:
    item['axiom_content']=item.pop('extra_input')
    item['status']='optional_calibration_extension' if item['id']=='SR-S' else 'formal_component_of_current_cognitive_axiom_system'
ledger['inherited_cognitive_axioms']=items
for item in ledger['open_connections']:
    if item['id']=='MATTER':item['must_not_assume']='须在同一模型核验钟、控制、记录与交互的共同预测、守恒和反作用，区分公理后果与具体机制的实现'
    if item['id']=='JOINT':item['question']='当前完整认知公理体系及各物理连接是否在同一模型和共同有效域内联合满足？'
put(name,json.dumps(ledger,ensure_ascii=False,indent=2)+'\n')
name='research_cognition_physics/archive_1106_/文件索引.md'
text=get(name)
text=replace(text,'[继承输入与待证连接]','[当前认知公理体系与待证连接]')
text += '\r\n- [认知公理体系措辞与兼容性修订](_shared/axiom_wording_revision/README.md)\r\n'
put(name,text)

# Extend the already existing outer layer, with no new wrapper or frozen-manifest change.
entries=[]
for name in nav_names+[paper_name]:
    raw=originals[name]
    prior=relocation_original_bytes(raw,old_entries[name]) if name in old_entries else raw
    revised=updates[name]
    entry=dict(original=name,destination=name,original_sha256=sha(prior),current_sha256=sha(revised),original_bytes=len(prior),current_bytes=len(revised),edits=text_edits(prior.decode('utf-8'),revised.decode('utf-8')))
    assert relocation_original_bytes(revised,entry)==prior
    entries.append(entry)
assert layer['unchanged_support_sha256'].pop(paper_name)==sha(originals[paper_name])
layer['entries']=entries
layer['scope']='Four live root navigation files and one current stage paper; reversible wording-only revision of cognitive axiom-system identity'
layer['wording_revision']='research_cognition_physics/archive_1106_/_shared/axiom_wording_revision/README.md'
put(layer_name,json.dumps(layer,ensure_ascii=False,indent=2)+'\n')
name='scripts/run_research_motion_active.py'
text=get(name)
text=replace(text,'Only four live root navigation files have a new reversible outer layer. Frozen','Four live root navigation files and the current SR paper share this reversible\nouter layer. Frozen')
text=replace(text,"('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}\n", "('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}\nPAPER = 'research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md'\nLIVE_MARKDOWN = NAVIGATION | {PAPER}\n") if '\r\n' not in text else replace(text,"('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}\r\n", "('README.md', 'research_direction.md', 'RESEARCH_STATE.md', 'ROADMAP.md')}\r\nPAPER = 'research_cognition_physics/认知时空操作公理与洛伦兹运动学_阶段论文.md'\r\nLIVE_MARKDOWN = NAVIGATION | {PAPER}\r\n")
text=replace(text,'exact four changed readers','listed changed Markdown readers')
text=replace(text,'if key in NAVIGATION else path','if key in store.motion_entries else path')
text=replace(text,'if name not in NAVIGATION or','if name not in LIVE_MARKDOWN or')
text=replace(text,'Only the four distinct root navigation files may change','Only the four root navigation files and current SR paper may change')
text=replace(text,"if set(self.motion_entries) != NAVIGATION:","if not NAVIGATION <= set(self.motion_entries) <= LIVE_MARKDOWN:")
text=replace(text,'return dict(files=4, manifest_sha256=self.motion_layer_sha256,','return dict(files=len(self.motion_entries), manifest_sha256=self.motion_layer_sha256,')
text=replace(text,'old_archive_and_paper_modified=False, new_scientific_rounds=0,','old_archive_modified=False, current_paper_revised=PAPER in self.motion_entries,\n                    published_paper_recoverable=True, new_scientific_rounds=0,')
text=replace(text,'newest four-document inverse layer','newest live-Markdown inverse layer')
compile(text,name,'exec')
put(name,text)
name='research_cognition_physics/_migration/motion_active_1106/README.md'
text=get(name)
text=replace(text,'本层只给研究根目录的README、research_direction、RESEARCH_STATE和ROADMAP添加当前目标栏。1086—1105档案、阶段论文、旧运行器与旧清单完全保留。','本层维护研究根目录README、research_direction、RESEARCH_STATE、ROADMAP的当前导航，以及当前洛伦兹阶段论文的认知公理体系措辞。1086—1105档案、旧运行器、旧清单和原论文发布回执保持原字节；历史复算通过同一层的可逆差分恢复发布时论文。')
text=replace(text,'`layer.json`保存4项可逆文本差分及前后哈希，不存整份文档副本；','`layer.json`目前保存5项可逆Markdown差分及前后哈希，不存整份文档副本；')
text=replace(text,'`checks.json`保存本次全部旧层兼容性校验与1105默认复算的实际输出。','`checks.json`保留启动时的旧层兼容性校验；本次措辞修订及当前校验见[修订记录](../../archive_1106_/_shared/axiom_wording_revision/README.md)。')
text=replace(text,'新入口先验证真实导航文件并在内存中恢复结项时字节','新入口先验证真实活动文档并在内存中恢复结项时字节')
text=replace(text,'以后更新这4份活动导航时','以后更新列入本层的活动Markdown时')
put(name,text)

# Check before any write; do not overwrite a concurrently changed document.
assert not (HERE/'changes.json').exists(),'Revision already applied'
for name,raw in originals.items():assert (ROOT/name).read_bytes()==raw,('Concurrent edit',name)
record=dict(schema='cognitive_axiom_wording_revision_v1',scope='Current explanatory prose only; no new mathematical axiom or research round',files={name:dict(before_sha256=sha(originals[name]),after_sha256=sha(raw)) for name,raw in updates.items()},paper_math_count=len(math(paper_before)),paper_math_sha256=sha(json.dumps(math(paper_before),ensure_ascii=False).encode()),paper_proof_sections_4_to_7_sha256=sha(paper_before[a:b].encode()),goal_objective_quote_sha256=sha(goal_quote_before.encode()),goal_receipt_sha256=sha((ROOT/'research_cognition_physics/archive_1106_/_shared/goal_start.json').read_bytes()),old_history_after_navigation_marker_unchanged=True,new_axioms=0,new_rounds=0,scientific_increment=0)
for name,raw in updates.items():(ROOT/name).write_bytes(raw)
(HERE/'changes.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
(HERE/'README.md').write_text('''# 当前认知公理体系措辞修订

依据用户最新要求，当前文档统一将F＋U＋C＋P四原则、CO1—CO5认知原子操作五公理、SR-T/A/R/O/F时空操作五公理表述为同一可扩展认知公理体系。四原则不预设为最终完备清单；体系逐阶段扩充、核验、整理与收敛。SR-S仍为可选定标扩展，G未采用。

四原则、原子关系操作和完整时空任务分层承担职责。洛伦兹运动学阶段已经结束；实际钟、全部纳入研究的影响及物质交互之间的共同连接仍按当前阶段继续研究。没有新增数学条款、研究轮次或科学计数，没有把公理采用改写成某项演绎证明。

## 修改及保留

- 更新当前洛伦兹论文、四根导航当前区及1106阶段说明、目标说明、依赖账和索引。
- 论文公式、定理证明及量词边界保持；目标原文和goal_start原始回执保持。
- 历史报告、旧采用回执、冻结运行器与归档清单保持。四根导航“历史导航”之后保持。
- 复用已有motion_active_1106活动差分层，增加当前论文的可逆措辞差分。最外层运行器仅扩展列名与读接口，未新建包装层、目录副本或ZIP。原论文发布版可在内存中按哈希恢复。

## 记录与复核

[逐文件变更哈希](changes.json) · [实际复核结果](checks.json)。revise_wording.py是本次一次性维护程序，已有变更记录时拒绝重复执行。启动时的旧checks、原始采用与发布记录保留其历史版本；本记录只说明当前措辞修订。

历史复算入口仍为scripts/run_research_motion_active.py。复核包括正文公式与证明相等、当前本地链接和格式、全部旧层兼容性及1105默认结果重现。没有进行图像检验。
''',encoding='utf-8')
print(json.dumps(dict(status='APPLIED_WORDING_ONLY',files=len(updates),paper_math_unchanged=True,proof_sections_unchanged=True),ensure_ascii=False))
