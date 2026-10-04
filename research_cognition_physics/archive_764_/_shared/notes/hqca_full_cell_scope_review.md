# 第516轮后未编号审计：数据通用性不自动成为完整胞元的操作能力

日期：2026-09-28。最新编号科学状态仍为516／2531。本次应用既有对称性工具，给418原模型一项新的接口判定；不另登记517。代码：[hqca_full_cell_scope_probe.py](../code/hqca_full_cell_scope_probe.py)；结果：[JSON](../results/hqca_full_cell_scope_probe_results.json)。3项独立诊断，不加入编号科学检查。

## 1. 当前目标与被检验的具体命题

上一目标回合修正了“不同接入口必须有完全相同内部未来”的过强要求，改变了后继验收，记为进展。本回合返回认知原则到实际物理过程的连接：342—343已经证明原形式操作合同不选择自然传播或指定维数，但尚不能把一份通用自治计算装置直接补成后来内部实现要求的完整反模型。

本次检验的候选是：**418允许任意加长经典程序并减小工作数据通道误差，是否因此也能把每个完整12维物理胞元当成具有全量子操作能力的对象？**

答案在明确输入范围内为否，而且下界不随链长或等待时间减小。这既不是原认知原则的反例，也不是要求所有认知必须按原胞元划分。把完整胞元选作对象、保持其固定物理输出因子，是本次被审计的具体识别；编码的工作qubit仍保留418的结论。

### 去重与成熟工具

- [345—346](../../../archive_342_369/research_note_345.md)已有编程与程序返回范围，但没有本次418规则的全胞元误差地板。
- [404](../../../archive_370_428/research_note_404.md)针对另一QCA的固定正常程序和长时返回，不能替代本次全链长、任意时间、精度可变的结论。
- [418](../../../archive_370_428/research_note_418.md)已经明示未知输入是工作数据；本次不撤销其通用性，只审向全部物理胞元的升级。
- [421§4](../../../archive_370_428/research_note_421.md)、[447](../../../archive_429_466/research_note_447.md)已经包含守恒量、相干辅助和相对信息边界。对称性不能免费产生所缺参考，是成熟工具；本次只核原规则和指定接口，不声称发现新超选择定理。
- [Bartlett、Rudolph、Spekkens综述](https://arxiv.org/abs/quant-ph/0610030)讨论参考系、超选择限制及关系编码。以下有限维证明自含，不引用该综述替代本模型计算。

## 2. 原规则实际守恒什么

每胞元为六维程序乘数据qubit，程序字母为点、箭头、I、W、S、T。严格复用418代码和式(2)，不重新选择自然H：

$$
\mathcal H_j=\mathbb C^6_{p_j}\otimes\mathbb C^2_{d_j},\quad
R=\sum_{A\in\{I,W,S,T\}}\bigl(
|A\cdot\rangle\langle\cdot A|\otimes I+
|A\triangleright\rangle\langle\triangleright A|\otimes A\bigr),
\quad h=-(R+R^\dagger),\quad H_L=\sum_jh_{j,j+1}.
\tag{1}
$$

每项只交换两个程序字母的位置，不改变其数量。对任一字母s定义其本地占据q_s及全链计数Q_s，得到

$$
q_s=|s\rangle\langle s|_p\otimes I_d,\quad
Q_s=\sum_jq_s^{(j)},\qquad
[h,q_s\otimes I+I\otimes q_s]=0,
\quad[H_L,Q_s]=0.
\tag{2}
$$

原两胞元144维矩阵的每个非零跃迁均逐项保持六种整数计数，程序核验不依赖浮点判定本征值为零。局部恒等式求和证明所有有限链长；不靠扫描链长外推。

## 3. 全胞元操作的有限精度障碍

固定旧胞元j作为输入和输出对象A。其余所有可用程序、读者、控制器及参考资源合称辅助E，初态sigma与未知A／被动参考独立。只需假定它对T计数不变，不要求各辅助独立，也不要求它在全部程序基上对角。

完整实现由同一H的等待和不引入该对称性破坏的读取／反馈构成；更一般地，允许任意相应协变CPTP过程Lambda，因而以下障碍比仅等待H更宽。令q=q_T^(j)，Q_E=Q_T-q，有

$$
[\sigma_E,Q_E]=0,\qquad
\Phi_\sigma(\rho)=\operatorname{Tr}_{\bar A}\Lambda(\rho\otimes\sigma_E),\qquad
\Phi_\sigma(e^{i\theta q}\rho e^{-i\theta q})
=e^{i\theta q}\Phi_\sigma(\rho)e^{-i\theta q}.
\tag{3}
$$

后式来自把输入变换扩成整体变换、使用sigma不变性与Lambda协变性，再作偏迹。输出旁观者及完整记录仍在整体内部，偏迹只定义本次胞元接口。若另授予不协变程序门或读数，已改变合同，不能偷加后继续引用本界。

取胞元酉V在程序I、T二维子空间作Hadamard，其余程序字母恒等，并对数据恒等。它是“完整胞元允许全部酉”必须包含的一个具体门。对任意数据态rho_d，有

$$
V|I\rangle_p=|+\rangle_p:=\frac{|I\rangle_p+|T\rangle_p}{\sqrt2},\quad
X_{IT}=(|I\rangle\langle T|+|T\rangle\langle I|)_p\otimes I_d,\quad
\|X_{IT}\|=1,\quad \rho_0=|I\rangle\langle I|_p\otimes\rho_d.
\tag{4}
$$

rho_0对q不变，因此实际输出对q不变，I与T间的块严格为零。目标输出的X_IT期望为1。以D为半迹距、d_diamond为半diamond距离，迹范数对偶性给

$$
\operatorname{Tr}(X_{IT}\Phi_\sigma(\rho_0))=0,\quad
\operatorname{Tr}(X_{IT}V\rho_0V^\dagger)=1,
\qquad
\boxed{d_\diamond(\Phi_\sigma,\operatorname{Ad}_V)
\ge D(\Phi_\sigma(\rho_0),V\rho_0V^\dagger)\ge\tfrac12.}
\tag{5}
$$

该界覆盖任意有限辅助数量、经典程序混合、等待时间以及满足式(3)的有限完整记录反馈；取这些实现的通道范数闭包也不能绕过。它不要求精确操作、机械释放、永久输出或恢复裸H_A。

X_IT是检验“完整矩阵对象”主张的数学见证，不宣称受限装置已能自行测它。若只把对称不变效果作为实际权限，原胞元本就不具备所主张的完整复矩阵操作接口；关系编码可以是另一种对象识别，未被排除。若讨论条件成功分支，还必须逐分支协变或明确结果标签的群作用；非选择通道协变不足以自动保证每个后选择分支都满足式(5)。本次主结论不依赖后选择。

## 4. 相干辅助改变输入条件，但不能免费生成

允许任意sigma_E，定义对其总T计数的twirl。用同一个实现把sigma替换为twirl后的状态，式(5)适用；对初始辅助的迹距收缩和三角不等式给

$$
\mathcal T_E(\sigma)=\frac1{2\pi}\int_0^{2\pi}
e^{i\theta Q_E}\sigma e^{-i\theta Q_E}\,d\theta,\qquad
a_E=D(\sigma,\mathcal T_E\sigma),\qquad
d_\diamond(\Phi_\sigma,\operatorname{Ad}_V)\le\epsilon
\ \Longrightarrow\ a_E\ge\tfrac12-\epsilon.
\tag{6}
$$

这是必要条件，不是相干量达到该数就足以实现V。全部能参与执行的相位参考和控制器都必须包含在E内；不能先用未入账的参考生成程序相干。对称初始资源加本规则和对称仪器不能生成违反总计数对称性的资源，这是式(3)的直接后果。它不禁止模型一开始就拥有明示的相干资源。

一个同H的边界例说明不能把式(5)误写成“该规则永远不能得到相干程序”。取三胞元，数据全0；所有门都保持00，此子空间对原H严格不变。一个空位及两个门字母形成三点路径。在时间t=pi/sqrt(2)有

$$
e^{-itH_3}|\cdot,x,y\rangle_p|000\rangle_d
=-|x,y,\cdot\rangle_p|000\rangle_d,
\quad t=\frac\pi{\sqrt2},\quad x,y\in\{I,T\};
\qquad |\cdot,I,+\rangle\longmapsto-|I,+,\cdot\rangle.
\tag{7}
$$

证明只用原局部跃迁：空位的三个位置构成等权开放三链，生成元为负邻接矩阵，端点在该时刻完美转移，门序列不变。线性性允许y为相干叠加。相干辅助在原固定中间胞元准备出plus；把相同辅助人口改为I、T经典均匀混合，输出也为混合，与plus的半迹距恰为1/2。

但这不是任意未知胞元的Hadamard实现。若旧中间程序与被动R纠缠，式(7)实际上给

$$
|\cdot\rangle_1\otimes|\psi\rangle_{2R}\otimes|+\rangle_3
\longmapsto-|\psi\rangle_{1R}\otimes|+\rangle_2\otimes|\cdot\rangle_3,
\qquad \operatorname{supp}\psi_2\subseteq\operatorname{span}\{|I\rangle,|T\rangle\}.
\tag{8}
$$

旧未知信息被移到左胞元；固定中间输出是替换通道，而不是V。程序对完整两维输入和参考核验了该转移，避免只看已知I输入就宣布操作能力补齐。这里没有生成相干来源，也没有证明完整C。

## 5. 复算与对目标的影响

3项检查分别核：原144维局部规则的六项整数计数；24组完整12维输入Choi协变及特定状态误差；原三胞元零数据部门的相干来源边界、经典混合对照和未知输入／参考转移。全部链长、时间和任意输入的结论来自式(2)—(6)，有限Choi距离仅为诊断，不充作diamond上界。没有重跑旧实验或图像检查。

运行已有Python与NumPy：

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics/archive_231_/hqca_full_cell_scope_probe.py'
```

结果只允许首次独占写入；修改前的初版代码保留为[hqca_full_cell_scope_probe_before_reference_check.txt](hqca_full_cell_scope_probe_before_reference_check.txt)。本稿通过独立审阅后逐字节复制为正式未编号审计，旧975份保护证据保持。

**认知动机：** 内部执行不能只在一部分数据上成立，就把程序、控制和记忆的同类能力视为已经完成。

**额外建模输入：** 418固定模型、固定完整胞元接口、辅助对称性与协变操作菜单。它们不是全体认知必须采用的要求。

**解析结论：** 这一具体升级失败，半diamond误差至少1/2；相干资源或关系编码等改变合同的路线仍开放。

**物理范围：** 既未从认知推出三维，也未反证宏观有效三维或GR。原生支撑为一维不自动排除三维有效部门；不能靠可编译任何报告来证明空间必须具有某个维数。

**后继：** 停止仅增加经典程序来修复完整胞元的候选。若继续该不足性分支，只接受有来源的相干辅助和完整未知输入接口，或另一份明确对象识别；普通相干电池与一般编码工具已在421、447、486，不重复开轮。正向空间分支仍须给实际定位任务及跨尺度选择，不能把本审计当作又一项所有空间都必须通过的新门槛。
