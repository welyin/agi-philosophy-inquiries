# 第777轮：辅助接触的BV正则分裂与同态自由约化

日期：2026-10-05。接[776](research_note_776.md)、[冻结入口](777/drafts/STATUS.md)。[代码](777/joint_auxiliary_bv_reduction.py)、[结果](777/joint_auxiliary_bv_reduction_results.json)、[核验](777/research_round_777_checks.json)、[当前联合条件账](../_shared/notes/unified_physics_condition_ledger_current.md)。三组精确检查，无图像检查、独立代理或应用目标变更。

## 1. 新连接与验收边界

**原770规范固定作用中的16个辅助乘子，可以连同反场作局部BV正则变换，精确分裂成传播部门与一个可缩辅助对。776的反场二次源接触，恰是这个传播部门保持经典主方程所需的二次反场项。** 它不是一个可任意删掉的归一常数，也不需要新的物理物种。

同一变换还保原自由二点对象，给自由星代数及自由BV复形的约化；物理自由态不另选。**但时间序乘积不能直接逐项投影，局部接触仍非零。** 因此本轮关闭“经典BV消元及同态自由代数能否共同实现”的接口，没有签收原全部高次N1/N2或连续量子BV推送。

|层次|内容|
|---|---|
|认知动机|同一系统的不同描述应保留操作与一致性关系，不能只保状态变量的数目|
|原输入|753联合背景、原四维／群／物种／作用／参数、770固定背景的线性规范条件与可逆N|
|继承证明|768／775同一自由核及BV导子；773局部收缩；774保首项提升；776高次基本辅助源与Wick输送|
|本轮对象映射|原辅助源、反场配对、经典主方程和同一自由星乘可采用同一局部字典|
|未增加的内容|无新增认知公理、粒子、参考态、经典动力学或量子反常抵消假设|
|未完成|共同高阶时间序规范化、实际圈反常及相互作用正态、实际仪器、全局／跨尺度映射|

上轮776产生正式报告及可复算结果，为有效进展。本轮先读导航、776和777入口并检查进程，未发现运行中的Python研究；再回查770、773—775及全阶段辅助消元关键词。773已讨论非最小可缩对，本轮不重算该一般事实；新增的是**与776同一个Schur接触、反场字典及原自由态的共同实现**。空间382—386、425、522—523及限定反例维持原范围。

## 2. 固定符号，防止反场号混淆

采用右导数作用于第一个输入、左导数作用于第二个输入的BV括号，s作用在左：

$$
 (F,G)=\sum_i\left(
 F\frac{\overleftarrow\delta}{\delta\phi^i}
       \frac{\overrightarrow\delta G}{\delta\phi_i^*}
 -F\frac{\overleftarrow\delta}{\delta\phi_i^*}
       \frac{\overrightarrow\delta G}{\delta\phi^i}\right),
 \qquad sF=(S,F).
 \tag{1}
$$

连续式包含积分及分部积分。物理玻色涨落记v，原ghost为c，anti-ghost为bar c，辅助场为b。设a为bar c的**规范共轭反场**，rho为b的反场。v、b、a为偶，bar c、rho为奇。本约定中使s bar c=b的非最小项为−a b。因此，776中把+b的外部系数记为a_src时，**a_src=−a**；二次项a_src N⁻¹a_src／2仍为+aN⁻¹a／2。这只是明示符号字典，不回写冻结报告中的源约定。

令Y=sv=R(B₀+v)c。使用770恢复普通规范参数后的闭规范表示：原最小BV作用S_min满足经典主方程，并在v*上仿射，其v*项为v*Y。原费米、ghost反场和需要的frame表示仍包含在S_min中。此处只变换原16个b部门；不凭此签收额外frame部门的量子规范化。

G、N只依赖固定B₀，N=Nᵗ点态可逆，U=N⁻¹G。原规范固定作用写为

$$
 S_{\mathrm{gf}}=S_{\min}
 +\langle b,Gv\rangle-\tfrac12\langle b,Nb\rangle
 -\langle\bar c,GY\rangle-\langle a,b\rangle .
 \tag{2}
$$

这是原线性规范固定，不要求Y对v线性。将量子场依赖的N或一般高阶规范固定换进来会增加项，不属于本轮证明。

## 3. 两个必须同时进行的正则变换

在776的β=b−Uv之外，定义

$$
 \widetilde v^*=v^*+U^{\mathrm t}\rho,\qquad
 e=\beta+N^{-1}a,\qquad
 \eta=\bar c+N^{-1}\rho,
 \qquad (v,c,a,\rho)\ \text{其余不变}.
 \tag{3}
$$

第一步是Schur字段变换的余切提升；Uᵗ是**完整形式转置**，含系数与分部积分。第二步同时移动辅助字段和anti-ghost，而非只给β补一个值。直接核基本括号得到(v,tilde v*)、(e,rho)、(eta,a)仍为规范配对，其余交叉为零。例如(e,eta)中，(β,N⁻¹rho)与(N⁻¹a,bar c)恰好抵消。

变换和逆变换均为有限阶局部微分映射；不需要Green逆或增加边界条件。工作区仍在既有正则片内，测试紧支；有实际边界时不能略去分部积分边界项。U含背景系数，其转置不能误用为U本身。

将逆字典b=e+Uv−N⁻¹a、bar c=eta−N⁻¹rho、v*=tilde v*−Uᵗrho代回(2)，得到**完整非线性作用的精确分裂**：

$$
 \begin{split}
 S_{\mathrm{gf}}&=S_{\mathrm{red}}-\tfrac12\langle e,Ne\rangle,\\
 S_{\mathrm{red}}&=S_{\min}(v,\widetilde v^*,\ldots)
 +\tfrac12\langle Gv,N^{-1}Gv\rangle
 -\langle\eta,GY\rangle-\langle a,Uv\rangle
 +\tfrac12\langle a,N^{-1}a\rangle .
 \end{split}
 \tag{4}
$$

证明中不能漏掉两个抵消：S_min的反场移动给−rho UY，anti-ghost移动给+rho UY；辅助平方与−ab合并给−eNe／2＋aN⁻¹a／2。S_red不含e或rho。

由于(3)保BV括号，原主方程及(4)直接给

$$
 (S_{\mathrm{red}},S_{\mathrm{red}})_{\mathrm{red}}=0,
 \quad se=0,\quad s\rho=-Ne,
 \quad s\eta=Uv-N^{-1}a,\quad sa=GY.
 \tag{5}
$$

特别地，s²eta=UY−N⁻¹GY=0。若删去aN⁻¹a／2，只剩s eta=Uv，平方会留下UY；这个非零项是离壳的，不能借ghost自由方程令它消失。776源求和得到的二次项与这里的BV要求是同一项：两项接口由一个原作用共同决定。

## 4. 允许消元的理想，以及不允许声称的同构

直接令β=0不构成离壳BV约化，因为sβ=−UY一般非零。正确的辅助对是e、rho；其微分理想在s下稳定。令pi将e、rho及其全部jet置零，i为剩余代数的包含。取奇导子kappa(e)=−N⁻¹rho、kappa(rho)=0，并向jet延伸使其与全导数交换，令N_aux数这对变量及其jet：

$$
 s\kappa+\kappa s=N_{\mathrm{aux}},\qquad
 h=\kappa N_{\mathrm{aux}}^{-1}\ (N_{\mathrm{aux}}>0),
 \qquad sh+hs=1-i\pi .
 \tag{6}
$$

N随时空变化时，κ(∂e)=−∂(N⁻¹rho)，所以逆系数的全部导数仍在，未违背776的接触jet要求。每个固定阶项有限局部，紧支性和模散度的商保持。因此该配对不改变局部BRST同调H(s)或H(s|d_H)。

成熟依据是[Barnich—Brandt—Henneaux，hep-th/9405109，§15、Theorems 15.1—15.2](https://arxiv.org/pdf/hep-th/9405109)。这里没有把辅助消元定理当原创；本轮给出了满足其局部代数可解前提的**原对象、精确变换和收缩**。

一个重要边界是：pi为链映射，不是整个离壳代数的BV括号同态。例如

$$
 \pi(e^A(x),\rho_B(y))=\delta^A_B\delta(x,y),\qquad
 (\pi e^A(x),\pi\rho_B(y))_{\mathrm{red}}=0.
 \tag{7}
$$

剩余变量使用其自身非退化BV括号；不能将这次消元描述成“对任意输入直接取商仍保全部反括号”。真正不变的是上述局部同调及由包含实现的等价，量子接触还须单独输送。

## 5. 同一自由量子代数能跟随，但时间序乘积不能裸投影

反场不收缩，775给Wββ=Wβv=0，其他辅助交叉也为零。在(3)后，We•=Wrho•=0，剩余传播变量的W与原自由物理准备相同。故对适用的微因果泛函和形式反场多项式，普通自由Wick星乘满足

$$
 \pi(F\star_W G)=\pi F\star_{W_{\mathrm{red}}}\pi G,
 \qquad \pi D_0=D_{0,\mathrm{red}}\pi .
 \tag{8}
$$

第一式逐收缩成立，因为W没有指向待投影变量的腿；第二式由(4)的二次部分分裂、775已有自由Wick/BV交织及线性正则字典成立。原自由物理商和状态的识别随之保持；未构造新的相互作用Hilbert表示或荷的定义域。

**(8)不适用于原始时间序乘积。** e的自由先进／推迟逆仍为−N⁻¹δ，故在适配的辅助零Wick减除下

$$
 \pi T_2(e(x),e(y))=-i\hbar N^{-1}(x)\delta(x,y),
 \qquad T_{2,\mathrm{red}}(\pi e(x),\pi e(y))=0.
 \tag{9}
$$

原770减除还需776的Q_C输送。这不是重复宣告二点接触存在；新的结论是：**同一个pi可保自由星乘与BV链，却不能保时间序族。** 因而经典辅助等价和同一自由态并不足以签收连续量子BV推送或全N1/N2。

原a=j=反场=0时，e=β、eta=bar c、tilde v*=v*；传播自由算符仍为原D₁、D₀及原费米算符。原局部关系观测O(v)、771预定的线性修复C₁(v)也不变。变换不增加原粒子内容或删去任何物理自由度。

## 6. 原基本源子菜单确实不闭合

不能仅用776“J不依赖量子涨落”的子菜单代替整个插入代数。原规范固定费米泛函已有bar c·b；在原坐标，

$$
 s(\bar c\,b)=b^2,
 \qquad s(\bar c\,GY)=b\,GY
 \quad(sY=0).
 \tag{10}
$$

第二项的辅助系数含动态ghost和通常的物理涨落，不是外部固定J。第一项显示一般BV伙伴已包含辅助复合插入；这不是等到猜测某个未知圈反常后才可能出现的情形。

式(4)提供了有意义的替代：可以在仅含传播字段、ghost、原费米及其反场的**约化完整BV作用**上推进规范化，再检查回到原插入的局部接触字典。它不是改换物理作用来回避困难，也不提前宣称这种量子回提已经成立。

773的局部形式同调与774的修复可通过经典局部等价转述；线性字典保总次数，特别有

$$
 \pi C_1(v)=C_1(v),\qquad iC_1(v)=C_1(v).
 \tag{11}
$$

因此775已经核实的首项相容不需要再次附加。真正仍缺的是**同一套传播时间序乘积、基本场方程与局部QAP，以及复合插入的量子回提**，不是再证明一个辅助二点矩阵可逆。

对接[Fröb，1803.10235，式(137)—(140)、Theorems 5、8—9](https://arxiv.org/pdf/1803.10235)时需保留a²：该文允许自由二次BV生成元含二次反场，且区分不传播的反场与动态字段。它未把任意未核验的T族自动变成共同规范化；下一轮要用当前约化对象逐项验收。

## 7. 精确诊断与负对照

采用既有Python，无新增依赖；有限超多项式用有理系数并严格保Grassmann符号。独立规范对照取sq=c、sr=q c、sc=0，故r−q²／2不变，

$$
 S_{\min}^{\mathrm{test}}=\tfrac12(r-q^2/2)^2+p\,c+t\,q c,
 \qquad Gv=\tfrac23q-\tfrac15r,\qquad N=\tfrac32.
 \tag{12}
$$

它是有非线性变换与非恒定Faddeev系数的诊断，不是原完整场论的替代模型。

三组通过：

1. 全部100个生成元括号精确保持；原／约化主方程及(4)精确成立。删除a²项后，s²eta=4c／9−2qc／15；漏掉反场移动留下−4c·rho／9＋2qc·rho／15。直接β=0的理想同样不稳定。
2. 231个非零超单项式核s²=0、(6)和链映射；(10)给精确b²。另核(7)的反括号投影缺陷为1，防止将链等价过度表述为任意BV乘积等价。
3. 两点差分模板、N=diag(2,3)核局部导数字典的转置顺序。用U代替Uᵗ使基本配对最大缺陷为7／6。模板只检查分部积分字典，不被当成新的空间模型。

```powershell
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics\archive_764_796\777\joint_auxiliary_bv_reduction.py' --check
```

连续、任意有限jet及原全部字段的结论由上述解析映射承担，不由有限单项式数目推广。未计算原圈反常、没有独立论文级审稿。

## 8. 对联合目标的增量与下一步

本轮将“辅助时间序接触”“原BV闭合”“同一自由物理态”“保首项局部修复”连接到同一个作用字典。没有新增四维、规范群或物种的推导，也没有消去原经典作用等独立输入。

当前账见[联合条件账](../_shared/notes/unified_physics_condition_ledger_current.md)：C01、C16、C19、C22中的一个共同描述接口得到补齐；全部量子过程、记录与反作用共同实现尚未完成。旧条件总账继续保留，不复制整份历史作为每轮新账。

下一项[778入口](778/drafts/STATUS.md)：在(4)的原约化传播复形上核高次因果规范化和QAP，特别审查二次反场项、同一T₁和基本场方程能否共同实现，并列出量子回提真正需要的接触。若只能在约化表示证明，须明确其与原物理关系观测的同一性及未覆盖的插入，不改写全N1/N2的验收口径。
