"""Reproduce integrated v1.1 from the byte-preserved v1.0 paper."""
from pathlib import Path
import hashlib

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'paper_versions/finite_quantum_v1.0.md'
TARGET=HERE.parent/'可组合认知结构与有限维量子理论_阶段论文.md'


def build():
    assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()=='51525ce376f665b9b21488a6b23c19ee55e7679fd7181302404faacb3266d0d7'
    text=SOURCE.read_text(encoding='utf-8')
    def replace(old,new):
        nonlocal text
        assert old in text,old[:80]
        text=text.replace(old,new)
    def section(start,end,replacement):
        nonlocal text
        a,b=text.index(start),text.index(end)
        text=text[:a]+replacement+'\n\n'+text[b:]

    replace('**第二阶段论文 v1.0 · 2026年9月20日**','**第二阶段论文 v1.1（整合版）· 2026年9月21日**')
    replace('研究范围：第223—228轮','研究范围：第223—230轮')
    replace('若明确要求操作路径连续，或另取可逆群闭性，则得到精确全酉权限。',
            '只要某一个现有对象存在一条非恒定连续可逆操作路径，就得到所有现有对象的精确全酉权限；非平凡Time已提供这一条件。')
    replace('本阶段以此结项；条件的最小性、具体哈密顿量、物理时钟、无限维极限、场论与引力留作后续问题。六轮共有55项可复算检查，',
            '有限自适应实验的总误差由各步骤误差之和控制，后选择需记入成功率。本阶段以此结项；原合同是否自动保证连续可逆种子、具体哈密顿量、物理时钟、无限维极限、场论与引力留作后续问题。八轮共有70项可复算检查，')
    section('### 2.2 精确实现的正则条件','### 2.3 连续可逆动力学的时间合同',r'''### 2.2 精确实现的连续种子条件

令$G_A$为实际可逆通道群，$G_A^0$为恒等连通分支，$H_A$为恒等路径分支。本文区分状态路径、实际操作路径及群的闭性。

**条件$R_{\mathrm{seed}}$。** 某个现有对象D上存在一条非恒定连续可逆操作路径：

$$
\gamma:[0,1]\longrightarrow G_D,\qquad\gamma\text{ 连续且非恒定}.
$$

非恒定指通道确实变化；不可观测整体相位不算。乘上起点的逆可将起点移到恒等。无需预先要求该路径遍及纯态集，或为每个对象分别提供路径。

v1.0的Reg表示两个可选条件：每个对象的恒等路径分支在纯态上传递（$C_{\mathrm{path}}$），或者每个实际可逆群均闭。§6将证明，在至少含一个非平凡对象的完整合同内，种子、全酉权限及两种Reg条件等价。种子减少了预先核验的内容，不宣称它是逻辑上严格更弱的独立公理。

原C的连通分支不能未经证明改成路径分支。本文未从F＋U＋C＋P推出种子必然存在，也未构造符合全部合同的无路径反模型。原合同下的逼近与带种子条件的精确定理分别成立。''')
    replace('它与C的纯态间可达性、Reg的控制正则性分别记账。连续可控不意味着已经指定自然演化。',
            '它与C的纯态间可达性和连续种子的存在性分别记账。非平凡Time提供种子，恒等时间群不提供种子；连续可控不意味着自然演化已选定。')
    replace('若另有Reg，则每个已有对象上的全部酉通道',r'若另有$R_{\mathrm{seed}}$（或v1.0的Reg），则每个已有对象上的全部酉通道')
    replace('K只确定到实标量恒等项。','K只确定到实标量恒等项。\n5. **有限实验逼近。** 固定有限深度、有限结果的自适应实验，各步带记录仪器的diamond误差为ε_j时，最终含参考输出的迹范数误差不超过Σ_jε_j，读出总变差不超过其一半。')
    replace(r'\mathrm{Reg}',r'R_{\mathrm{seed}}')
    replace('精确分支再用Reg','精确分支使用连续种子及简单Lie代数')
    replace('| 交换子及薛定谔形式 |','| 有限自适应实验 | 仪器稠密性、有限次替换及CPTP收缩 | 保留记录与参考的统计界 |\n| 交换子及薛定谔形式 |')
    section('### 6.3 从稠密到精确','## 7 全部测量、通道和有限仪器',r'''### 6.3 稠密群中的一条可逆路径

原合同下，选择AA中的实际门逼近$U\otimes I$；第二份A固定准备后丢弃，就逼近A上的任意酉。这些诱导近似通道未必可逆，不能直接算作$G_A$的元素。

设$G\subseteq PU(N)$稠密，N≥2，且含非恒定连续路径。其恒等路径分支H是G的非平凡正规子群。[Yamabe定理](https://doi.org/10.18910/6228)使H成为解析Lie子群，允许浸入而不要求闭。其Lie代数h非零，正规性使Ad(G)保持h；有限维线性空间的闭性和G的稠密性使Ad(PU(N))也保持h。因此

$$
0\ne\mathfrak h\triangleleft\mathfrak{su}(N)
\quad\Longrightarrow\quad
\mathfrak h=\mathfrak{su}(N),\qquad H=G=PU(N).
$$

最后一步使用简单Lie代数及恒等邻域：h为全体使H包含开邻域，而连通PU没有真开子群。闭包只用于证明线性空间的不变性，没有把极限操作预先算作实际操作。

### 6.4 一个对象的种子推广到所有对象

设D满足$R_{\mathrm{seed}}$，A任意。取S＝AD及Q＝SS。式(11)给$G_Q$稠密；在一个D因子上施加原路径，其余因子恒等，使$G_Q$含非恒定路径。§6.3给$G_Q=PU(d_Q)$。执行$U_A\otimes I_{\rm aux}$，固定准备并丢弃辅助，便精确实现A上的U及其逆；该等式保持任意外部纠缠参考，未知输入无需读取或复制。

所有现有对象因此获得全酉权限。下一节使全酉权限等价于全部有限仪器权限；全酉权限又允许在任一非平凡对象上选非标量K，构成非平凡连续可逆一参数群。因此在完整合同内，种子存在、某个非平凡连续可逆一参数群存在、全酉权限及全仪器权限等价。全PU自动满足$C_{\mathrm{path}}$和闭群；两种旧Reg的充分性已由227轮证明。

仍开放的分支是：是否可能所有对象均无非恒定可逆连续路径，却满足U、C及全部组合合同？本文未证明该分支存在或不可能。

### 6.5 有限乘积与逆函数定理

若种子为$\exp(tX)$，其中非零$X\in\mathfrak{su}(N)$，其实际共轭方向的实线性包在Ad(G)和Ad(PU(N))下不变，简单性使它为全部su(N)。可选$m=N^2-1$个实际$g_j$，使$X_j=\operatorname{Ad}_{g_j}X$成基。有限乘积

$$
f(t_1,\ldots,t_m)=\prod_{j=1}^{m}\exp(t_jX_j),
\qquad df_0(v)=\sum_{j=1}^{m}v_jX_j
$$

的导数可逆，逆函数定理使其像覆盖恒等邻域，实际群即全PU。每个共轭脉冲由实际门、种子和逆门有限组合，未用无限Trotter极限赋予权限。

该构造对共轭基的小误差稳定：Hermitian种子K的Hilbert–Schmidt范数归一时，每个共轭酉的算子误差≤ε，使基坐标矩阵的算子误差≤$2\sqrt m\,\varepsilon$。若理想基最小奇异值为σ，取$\varepsilon<\sigma/(2\sqrt m)$仍满秩。稠密门只需达到有限正容差，再用其自身的精确共轭方向应用逆函数定理。

结论给存在性，不给寻找门、设定精确参数或控制耗时的效率上界。随机混合恒等与Z共轭通常把纯态变成混态，不能用普通通道的连续性代替可逆种子。''')
    replace('Reg下该门精确可用',r'$R_{\mathrm{seed}}$下该门精确可用')
    replace('因此Reg下式(3)',r'因此$R_{\mathrm{seed}}$下式(3)')
    replace('没有Reg时，令S=ABE','仅用原合同时，令S=ABE')
    replace('## 8 连续可逆演化的生成元',r'''### 7.4 固定有限自适应实验的误差界

设最大深度为n，每层每个历史节点的带记录仪器以diamond误差ε_j逼近。节点有限，228轮允许分别选择实际仪器。保留经典历史与外部参考，按历史块选择后续步骤。第j层误差至多为ε_j乘以各块概率之和。逐步替换及CPTP收缩给

$$
\|\omega-\omega'\|_1\leq E=\sum_{j=1}^{n}\varepsilon_j,
\qquad d_{\rm TV}(p,q)\leq\min\{1,E/2\}.
$$

故对每个固定有限量子实验及任意正容差，原合同允许一个满足容差的实际实验。不交换“每个实验存在实现”与“一个实现适用所有实验”的量词，不保证相同成本或无限次实验等价；两个固定且不同的仪器仍可被区分。

等先验区分最终量子输出的最优成功率为$1/2+\|\omega-\omega'\|_1/4$，见[Watrous第3章](https://cs.uwaterloo.ca/~watrous/TQI/TQI.double.3.pdf)的Holevo–Helstrom定理和通道范数性质。本节将标准工具与权限稠密性衔接，不把统计界当作新的物理定律。

### 7.5 后选择需保留成功率

对事件未归一化输出A、B，概率p、q，完整误差E<p时有q>0，且

$$
\frac12\left\|\frac A p-\frac B q\right\|_1
\leq\frac{\|A-B\|_1+|p-q|}{2p}\leq E/p.
$$

分母不能省略。两个忽略输入的替换仪器均以概率p成功、以1−p失败；失败输出相同，成功输出为正交态。完整带记录diamond距离为2p，成功后的条件态迹距离却为1。平均等待成功需1/p次独立尝试，r次至少一次成功的概率不超过rp。固定p>0时仍可使E远小于p；不能在p趋零时要求不带分母的统一界。

## 8 连续可逆演化的生成元''')
    replace('Time不要求非零K。','Time不要求非零投影生成元；标量K使通道恒等。非标量K提供连续可逆种子，因此同时给精确操作完备性。')
    replace('六轮脚本使用既有Python 3.12.14、NumPy 2.3.5。55项检查覆盖以下内容；结项复核的实际执行情况、文件哈希和公式审核单独记录在[结项检查](archive_223_/stage2_closure_checks.json)。',
            '八轮脚本使用既有Python 3.12.14、NumPy 2.3.5。70项检查覆盖以下内容；全阶段回归见[最终检查](archive_223_/stage2_completion_checks.json)，整合版与目录迁移见[当前核验](archive_223_/integration_checks.json)。')
    replace('| **合计** | **55** | **对应六轮实现；不包含第一阶段旧测试的重复累计** |',
            '| 229 | 8 | 共轭方向满秩与扰动稳定性、局部逆映射、局部控制与凸随机化边界 |\n| 230 | 7 | 含参考的有限反馈实验、全部记录、后选择及最优区分界 |\n| **合计** | **70** | **对应八轮实现；不包含第一阶段旧测试的重复累计** |')
    replace('这些数值是实现精度，不是定理的容差前提。','229轮在N＝2、3、4、6、8核验共轭方向秩为N²−1；N＝2、4的近恒等目标反演残差低于3×10⁻¹⁵。230轮保留Bell参考及全部历史，核对1、2、4层反馈实验的加法误差界，并给稀有事件反例。样本不证明所有整数维对象存在。\n\n这些数值是实现精度，不是定理的容差前提。')
    replace('原合同是否足以给精确全酉权限，Reg能否削弱？本文不作必需性断言，条件性结论不依赖解决此问题。',
            '种子已经减少预先核验条件；原合同是否自动排除所有对象均无非恒定可逆路径的分支仍开放。没有完整反模型，不作必需性或最小性断言。')
    replace('下一阶段优先明示多体局域性和时间输入，再研究相应生成元及受控极限。',
            '新阶段按用户的物理生成猜想，先审计事件、影响关系与因果次序，再研究几何、粒子、质量、规范及引力，不把目标物理结构作为隐含前提。')
    section('## 11 结论','## 参考文献',r'''## 11 结论

第223—230轮将第一阶段的复状态结构推进为有限维量子操作理论：原F＋U＋C＋P给相容量子表示、全部有限仪器及固定有限自适应实验的任意精度逼近；某个对象上的非恒定连续可逆操作路径足以给所有现有对象的精确操作完备性。指定Time后得到薛定谔形式，非平凡Time又提供连续种子。

**第二阶段于第230轮完成收尾。**

v1.1将229—230轮直接整合进定义、综合定理、证明、实验范围和证据表；v1.0与补充稿作为历史版本保留。结项范围是有限维量子操作与连续可逆动力学形式的条件性重建。具体自然动力学、时空、场论、规范群与引力留待后续生成研究；原合同是否自动提供连续可逆种子仍开放。''')
    replace('## 附录A：证据与复算入口','10. J. Watrous. *The Theory of Quantum Information*, Chapter 3. [作者公开章节](https://cs.uwaterloo.ca/~watrous/TQI/TQI.double.3.pdf)。用于区分与通道范数工具。\n\n## 附录A：证据与复算入口')
    replace('阶段目录沿用用户指定名称`archive_223_`，本次仅将其状态标为223—228轮已结项，不迁移或覆盖既有研究文件。导航中的本次阶段决定优先于历史笔记末尾的“下一步”。',
            '阶段目录按用户要求重命名为`archive_223_230`。32份逐轮研究文件保持原字节。v1.0论文在[历史版本](archive_223_/paper_versions/finite_quantum_v1.0.md)保留，[原补充稿](archive_223_/STAGE2_ADDENDUM.md)保留当时语境；本整合版为当前论文。旧版中的路径按原写作位置解释，当前使用本附录与档案索引。')
    needle='| 228 | [仪器完备性](archive_223_/research_note_228.md) | [Python](archive_223_/instrument_completion_bridge.py) | [JSON](archive_223_/instrument_completion_bridge_results.json) |'
    replace(needle,needle+'\n| 229 | [连续种子](archive_223_/research_note_229.md) | [Python](archive_223_/continuous_seed_bridge.py) | [JSON](archive_223_/continuous_seed_bridge_results.json) |\n| 230 | [有限实验与收尾](archive_223_/research_note_230.md) | [Python](archive_223_/finite_protocol_closure.py) | [JSON](archive_223_/finite_protocol_closure_results.json) |')
    section('每轮检查记录见[阶段索引]','**编辑说明。**','''每轮检查见[阶段索引](archive_223_/README.md)。[原清单](archive_223_/STAGE2_MANIFEST.json)和[补充清单](archive_223_/STAGE2_CLOSURE_ADDENDUM.json)继续锁定两批证据；[当前核验入口](archive_223_/verify_archive.py)处理目录与版本映射，可显式复跑70项检查。旧核验脚本保留历史合同，当前使用以下入口：

```powershell
python -B -X utf8 research_cognition_physics/archive_223_230/verify_archive.py
python -B -X utf8 research_cognition_physics/archive_223_230/verify_archive.py --run-tests
```''')
    replace('本次为综合与结项，不新增第229轮。','本次整合不新增研究编号；229—230轮的假设、结论与边界已合入正文。')
    return text.replace('archive_223_/','archive_223_230/')


if __name__=='__main__':
    output=build()
    if TARGET.exists() and TARGET.read_bytes()!=SOURCE.read_bytes() and TARGET.read_text(encoding='utf-8')!=output:
        raise RuntimeError('Current paper differs from source and reproducible output; inspect before replacing.')
    TARGET.write_text(output,encoding='utf-8',newline='\n')
    print('Integrated v1.1 written; v1.0 remains byte-identical.')
