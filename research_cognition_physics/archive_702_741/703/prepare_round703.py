"""Prepare703; old scientific files and entry evidence remain frozen."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


def replace(text,mapping):
    pats=[r'(?<!\d)'+re.escape(k)+r'(?!\d)' if k.isdecimal() else re.escape(k)
          for k in sorted(mapping,key=len,reverse=True)]
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


names=('unified_physics_condition_ledger_703.md','round703_drafts/research_note_703_draft.md',
       'round703_drafts/final_review.txt','round704_drafts/STATUS.md',
       'round703_drafts/bounded_source_entry.py','round703_drafts/bounded_source_entry_results.json',
       'round703_drafts/bounded_source_entry.md','round703_drafts/check_and_publish_entry.py',
       'round703_drafts/entry_checks.json','round703_drafts/prepare_entry_publication.py',
       'round703_drafts/literature_scope_audit.json','round703_drafts/initial_geometry_diagnostic.py',
       'round703_drafts/initial_geometry_diagnostic_results.json')
protected=2823+3+len(names);assert protected==2839
ledger=(HERE/'unified_physics_condition_ledger_702.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：702原完整过程的共同热态、演化与记录',
                      '# 联合条件总账：703共同几何来源、热响应与原记录')
ledger=ledger.replace('接[701全账](unified_physics_condition_ledger_701.md)，回填[702报告](research_note_702.md)。[结果](joint_gibbs_preserving_transfer_results.json)、[核验](research_round_702_checks.json)。',
    '接[702全账](unified_physics_condition_ledger_702.md)，回填[703报告](research_note_703.md)。[结果](joint_geometry_thermal_limit_results.json)、[核验](research_round_703_checks.json)。')
ledger=ledger.replace('## 当前共同对象及仍存在的分支',
    '**703当前增量：** 原正几何域内可取固定束缚参照，保留原全部相互作用与几何变化。统一扇形半群及热迹界给一般几何一二阶热来源、Hessian接触、归一态与立即原记录的共同极限。623原响应已经存在，本轮补近似族的来源极限；含等待的真实来源导数、动态量子反馈及空间细化仍开放。\n\n## 当前共同对象及仍存在的分支')
updates={
    'C01 量子对象':'703固定参照正过程的热态及一般几何来源导数共用原Gauss空间',
    'C03 事件记录':'703固定原Kraus立即读取的联合后态及一二阶几何响应共同收敛',
    'C04 内部演化':'703各背景继承702真实时间零阶极限；有序正热时来源导数共同收敛',
    'C09 参考系统':'703自身热态归一对一般几何参数也可共同求导',
    'C20 尺度映射':'703固定图时间极限与几何一二阶导数交换；无空间一致速率',
    'C22 来源反作用':'703原一般几何热响应和Hessian接触与正近似族匹配，不另拟合应力'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    if row.startswith('|原固定图完整量子过程|'):
        cells=row.split('|');cells[2]+='；703一般几何热来源、接触及立即记录共同极限';rows[i]='|'.join(cells)
    for key,value in updates.items():
        if row.startswith('|'+key+'|'):
            cells=row.split('|');cells[2]+='；'+value
            if key=='C22 来源反作用':
                cells[3]=cells[3].replace('动态量子反馈、一般几何求导与共同重整化',
                    '动态量子反馈、原辅助／连续分支的一般几何来源匹配与共同重整化；原H_F热来源近似极限已由703接通')
            rows[i]='|'.join(cells);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n';at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''## 703一般几何来源与同一热过程共同收敛

- 623已给原共同域、来源算符界和真实二阶响应；不重复。703从正几何范围的参考动能及原位势取固定A=eta R，D(x)=H(x)−A仍含正动能和束缚，原总H完全不变。
- 有限系数的小复邻域中，D(z)是共同形式域的闭扇形族，Re D(z)≥delta R−c。Vogt–Voigt定理的输入逐项核验；无界几何来源不被当有界扰动。
- 固定A允许以Schatten Holder给T_n(z)统一迹界exp(c beta)Tr exp(−beta A)。实轴702迹收敛与Vitali/Cauchy共同给来源导数的迹范数极限。
- 任意原C²有限系数路径用链式法则接入，不新增解析历史公理。原Hessian接触、变字典导数和指定局部lapse分配保持来源，不独立拟合。
- 归一态、log Z来源及固定原Kraus立即读取共同收敛；有序正热时插入同样受控。含真实等待的近似来源导数仍需另接，624原连续响应存在并非空白。
- 原体积量子诊断中log Z二阶为10.8984146941，接触贡献−26.8366720040；漏接触会得37.7350866981。代码用非对易jet、独立Duhamel及差分三份核对。
- 原单次读口在对称诊断恒为1/2，正式校准使用两次真实读口及完整cq后态；初版诊断保留，不把零响应计成有效证据。解析全图结论不依赖有限径向矩阵。
- 四分支、699反例及全部空间边界不变；382—386、425、522—524逐项继承。C07／C08旧空间输入不因热来源定理消失。

## 本轮合并与下一项

C01／C03／C04／C09／C20和C22在原固定图中共享一族保几何热来源及接触项的正过程。新增表示参数不改变极限，但未选背景、温度、参数、维数、局部lapse分配或Einstein作用。

接[704](round704_drafts/STATUS.md)：核同一近似自身热态与真实等待的来源响应，回查621—624／643／655；不能改用独立实时分割或固定原态而省略身份。停止本轮热导数、有限格点及常数优化，目标不变。
'''
write('unified_physics_condition_ledger_703.md',ledger)
mapping={'joint_gibbs_preserving_transfer':'joint_geometry_thermal_limit',
         '701':'702','702':'703','703':'704','3309':'3311','3311':'3313',
         '1415':'1418','1418':'1421','2809':'2823','2823':'2839'}
verify=replace((HERE/'verify_round702.py').read_text('utf8'),mapping)
verify=verify.replace('import joint_flow_resolution_limit as prior_model',
                      'import joint_gibbs_preserving_transfer as prior_model')
verify=verify.replace('original full-model Gibbs, real-time, records and energy/entropy joint limit.',
                      'original geometric thermal-source jets, Hessian contacts and actual records.')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+'    names='+repr(names)+'\n'+verify[b:]
write('verify_round703.py',verify)
pub=replace((HERE/'publish_round702.py').read_text('utf8'),mapping|{'## 348.':'## 349.','## 253.':'## 254.'})
a=pub.index('summary=');b=pub.index('planned={}',a)
pub=pub[:a]+'''summary=('**第703轮完成：** [共同几何来源、热响应与原记录]({p}research_note_703.md)'
         '原完整固定图的几何一二阶热来源、Hessian接触、归一态和立即原记录有同一正近似极限；'
         '固定参照只改表示，C²路径不需新增解析公理。'
         '两组、十六式通过，最新703／3313，1421份编号科学文件、2839份保护证据。'
         '[核验]({p}research_round_703_checks.json)、[全条件账]({p}unified_physics_condition_ledger_703.md)。'
         '含等待的近似来源导数及空间连续仍开放；旧空间与目标不变。')
order=('**当前执行顺序（703后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[704共同热源与真实记录响应]({p}round704_drafts/STATUS.md)，'
       '核同一近似族的真实等待、来源及热规整，不重复623—624已存在的原响应；'
       '保留四分支和699反例范围，停止静态热导数优化。')
'''+pub[b:]
pub=pub.replace('原完整过程的共同热态、演化与记录','共同几何来源、热响应与原记录')
pub=pub.replace('原完整过程的热参考、记录和演化共用一个极限','原几何来源及热记录共用一个正极限')
pub=pub.replace('旧空间接口继承，热尾紧性不替代空间细化','旧空间接口继承，几何热来源不替代时空生成')
write('publish_round703.py',pub)
write('postcheck_round703.py',replace((HERE/'postcheck_round702.py').read_text('utf8'),mapping))
with (HERE/'round703_drafts/research_note_703_draft.md').open('xb') as f:
    f.write((HERE/'research_note_703.md').read_bytes())
review='''703 primary-agent mathematical/code/scope review; no independent agent.
Previous goal turn completed702 and703 entry: progress.
Detailed process inspection denied; ordinary Python process query found none.
Original623 domain/source response already solved; new gap is approximation/source derivative interchange.
Fixed eta R removes a small existing kinetic/onsite portion; remainder keeps ellipticity and confinement.
No original interaction removed. New finite-step representation is not silently identified with702.
R includes compact-group electric kinetic, all noncompact scalar directions and finite Fock identity.
Common sectorial form hypotheses and semigroup bound supplied before invoking Vogt/Voigt.
Complex parameter holomorphy uses independent real/imaginary coupling coordinates, not conjugation.
Fixed Gibbs anchor plus Schatten Holder bounds all n in trace norm.
Real trace convergence inherited702; Banach Vitali and Cauchy give genuine derivative limits.
Finite-dimensional coefficient analyticity plus C2 chain rule covers original C2 paths and Hessian contacts.
Normalization via original positive Z and locally uniform convergence, no chosen independent states.
Finite-step geometric source is derivative of entire split, not blindly inserted bare G.
Ordered positive thermal-time sources covered; real-time parameter derivatives explicitly still open.
Immediate original Kraus cq outputs, not arbitrary source-dependent waiting dynamics.
Volume diagnostic includes original kinetic -6 and onsite +6 scaling, canonical local mass volume0.
All jet cross terms retained; target Hessian independently checked by Duhamel including contact.
Single plus event is symmetry-constant; replaced by actual two-read event and cq jets before publication.
Initial diagnostic code/results preserved as historical draft, not authoritative final reproducer.
General shape coefficient bounds analytic over compact box; samples and radial matrix only calibrate.
No spatial-uniform, quantum geometry, chiral dictionary or GR claim.
384 removed Lipschitz stays removed;386/425 alternatives;523/524 reused; E not s.
Two groups,16 equations; next704 targets real-time/source interface without repeating old originals.
'''
for name in ('research_note_703.md','joint_geometry_thermal_limit.py','joint_geometry_thermal_limit_results.json',
             'unified_physics_condition_ledger_703.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round703_drafts/final_review.txt',review)
print('703 prepared; protected evidence',protected)
