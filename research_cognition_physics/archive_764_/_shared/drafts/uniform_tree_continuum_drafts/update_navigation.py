"""Append a checked, unnumbered theorem bridge with versioned navigation."""
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parents[1]
RESEARCH=HERE.parent
ROOT=RESEARCH.parent
report=json.loads((HERE/'uniform_tree_continuum_checks.json').read_text('utf8'))
assert report['all_reported_checks_passed']
assert report['total_protected_including_this_review']==1004
assert report['unchanged_numbered_scientific_tests']==2546
paths=[ROOT/'README.md',RESEARCH/'README.md',RESEARCH/'research_direction.md',
    RESEARCH/'RESEARCH_STATE.md',HERE/'README.md',HERE/'spatial_premise_closure_audit.md',
    HERE/'three_dimensional_four_conditions_review.md']
baseline={p:p.read_bytes() for p in paths}
snapshot=HERE/'uniform_tree_continuum_drafts/navigation_before_review'
snapshot.mkdir(exist_ok=False)
planned={}
for p,raw in baseline.items():
    enc='utf-8-sig' if raw.startswith(b'\xef\xbb\xbf') else 'utf8'
    nl='\r\n' if b'\r\n' in raw else '\n'
    body=raw.decode(enc).replace('\r\n','\n')
    prefix=('research_cognition_physics/archive_231_/' if p.parent==ROOT else
        'archive_231_/' if p.parent==RESEARCH else '')
    banner=(f'**519后宏观来源接续（不增轮次）：** [均匀关系树的CRT极限]({prefix}uniform_tree_continuum_review.md)'
        '精确对接517来源、条件二叉分枝树及成熟随机树定理；奇数规模、归一化及来源误差已核。'
        '在声明的瞬时树路径度量下，非退化紧致缩放为n⁻¹ᐟ²，极限实树的Hausdorff维数为2。'
        '这不确定动态传播空间或模型的物理维数；467的既有路线限制未变，不再做放大均匀树的三维拟合。'
        '3项独立诊断单列，编号仍519／2546；原1000份证据保持，新4份单列保护，总计1004份。')
    first,rest=body.split('\n',1)
    body=first+'\n\n'+banner+'\n'+rest
    if p==RESEARCH/'research_direction.md':
        body+='''

**519后成熟定理接续：** [均匀图来源的宏观范围](archive_231_/uniform_tree_continuum_review.md)及[核验](archive_231_/uniform_tree_continuum_checks.json)已完成。瞬时单位边树几何的均匀来源有CRT极限，非退化紧致缩放指数确定为1／2；这不生成实际测距或自主增长，也不能用Hausdorff维数2替代物理空间维数。保持编号519／2546，后继须继承1004份保护证据。

**下一项去重约束：** 不再通过增加I、改变统一幂次缩放、拟合CRT指数或静态随机游走来修复三维。继续核实际动态图中的位置／路径报告及有来源的复合端点，保留经典历史条件态和量子关系；任意新的几何选择机制须明确新增输入。单凭成熟普适计算或一个形式log U也不签收宏观空间，347—348及既有QCA审计直接复用。
'''
    if p==RESEARCH/'RESEARCH_STATE.md':
        body+='''

**519后未编号宏观来源审查：** [正文](archive_231_/uniform_tree_continuum_review.md)、[代码](archive_231_/uniform_tree_continuum_probe.py)、[结果](archive_231_/uniform_tree_continuum_probe_results.json)、[核验](archive_231_/uniform_tree_continuum_checks.json)、[只读复算](archive_231_/verify_uniform_tree_continuum_review.py)。3项诊断不并入2546，868份编号科学文件不变；1004份保护证据。三维、GR和阶段论文收尾仍未完成。
'''
    if p==HERE/'spatial_premise_closure_audit.md':
        body+='''

### 157.2 519后已有宏观随机树定理接续，不开520

上一目标回合完成519，真实来源的返回信号有小于信号的误差证书，属于有效进展，不累计阻碍。本回合按导航及最新结果检查，未发现520／521编号占用；没有重启或中止未知进程。

为避免在更大图上继续盲目拟合三维，查明517均匀全标号二叉树来源与条件Galton–Watson满二叉树精确等律。删指定叶后加左右序仅作数学辅助，每个平面形状有I!(I＋1)!个标号；原无序形状不等权。Aldous原Theorem23的gcd条件不能直接套{0,2}，改按Le Gall明确的奇数子序列说明，配合GH归一化及已有Hausdorff维数定理。

[审查](uniform_tree_continuum_review.md)确认：声明的瞬时路径度量在n⁻¹ᐟ²缩放下趋向Brownian CRT，其Hausdorff维数为2；更强幂次缩放塌缩，更弱者直径不紧，均不是三维修复。517有限源误差η_I→0可传递有界报告与弱分布极限；均值须另核幅度或一致可积性。既没有条件化包在根，也没有把条件历史都当作均匀图。

该接续补明特定来源的宏观极限，未减少实际距离读者、动态图传播等同关系或单个系统增大规模的输入。467已排除瞬时树近欧氏路线，故不为这份成熟分类新开520。三项诊断、6式与独立完整复算通过，正式稿与终审稿字节一致；新增4份证据，全部1000旧证据保持，总计1004。编号519／2546与868份编号科学文件不变。

后继仍回到实际定位菜单及复合端点的来源，不重复静态树放大、树随机游走或更多尺度指数。微观标签不必是空间点，但采用新的平均／历史／关系编码时，必须用同一实际报告核验它，不能仅更名以绕开旧反例。此结果不反驳整个宏观空间纲领或FUCP的所有可能加强，目标继续开放。
'''
    if p==HERE/'three_dimensional_four_conditions_review.md':
        body+='''

### 62.1 均匀瞬时关系树的极限不改变四项条件

[519后未编号接续](uniform_tree_continuum_review.md)把517来源精确接到Brownian CRT极限；它给候选图度量的缩放与分布，没有给实际位置群、半幅、重定向或方向接口。Hausdorff维数2也不表示获得了光滑二维空间。图增大不能自动修复467已证的瞬时树路线限制；动态传播与历史端点路线另核。编号519／2546，总保护证据1004，三维与GR仍开放。
'''
    planned[p]=body.replace('\n',nl).encode(enc)
    label=('project_README.md' if p==ROOT/'README.md' else
        'research_'+p.name if p.parent==RESEARCH else 'archive_'+p.name)
    with (snapshot/label).open('xb') as f:f.write(raw)
for p,raw in baseline.items():assert p.read_bytes()==raw,str(p)
for p,data in planned.items():
    assert p.read_bytes()==baseline[p],str(p)
    p.write_bytes(data)
with (snapshot/'manifest.json').open('x',encoding='utf8') as f:
    json.dump({str(p.relative_to(ROOT)):{'before':hashlib.sha256(baseline[p]).hexdigest(),
        'after':hashlib.sha256(data).hexdigest()} for p,data in planned.items()},f,indent=2,ensure_ascii=False)
print(json.dumps({'navigation_files_updated':len(planned),'snapshots_preserved':len(baseline)}))
