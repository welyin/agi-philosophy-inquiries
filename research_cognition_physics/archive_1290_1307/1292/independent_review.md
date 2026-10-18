# 1292 独立审阅：残余群、同源质量与有限响应

2026-10-18。审阅者：/root/spatial_parent_next。结论：**PASS_WITH_SCOPE**。作者最后补明完整零 jet 背景后，未发现阻断本轮准入目标的数学或范围问题。此结论不是长程目标验收，也不新增公理。

## 1. 实际审阅与独立方法

全文审阅本轮[证明](proof.md)、[作者核验](check.py)及[结果](results.json)、[准入](../_shared/round1292_admission.md)、[主报告](../research_note_1292.md)和[NEXT](NEXT.md)。核对1035原荷表、1060共同背景反常运输和1289稳定子／局部复合字典的继承终点。另读主报告引用的[共同因果锥范围筛查](../_admission/massless_channel_cognitive_scope_screen.md)，未把其有界历史审计当作全公理形式独立性证明。

[独立代码](independent_check.py)不导入作者或历史计算函数。全部数值核验使用 Fraction：

- 以原表示的有理相位检查中心下降、Yukawa及 SNN；独立枚举残余角／离散对并作轨道商。
- 改用 $(B,W^3)$ 排序，由动能在基向量及其和上的取值恢复 Hessian；独立核流矩阵、精确逆矩阵响应及删除凝聚的对照。
- 样本取 $k=-2,-1,0,1,3$，另取一组不同的正有理参数。样本不承担全整数及连续动量窗量词；这些由以下解析复核承担。
- 仅做同一标量动能与常时间势能源二次型的交叉检查，不将有限代数样本冒充完整量子 Ward、重整化或仪器来源证书。

## 2. 全局残余与 Spin 识别

从覆盖角直接求得 $W=\operatorname{diag}(e^{3i\theta},e^{-3i\theta})$。原中心三次幂将 $\theta$ 与 $\theta+\pi$ 识别，故 $\varphi=2\theta$ 才是普通 $2\pi$ 周期角；保留的色商为 $(\omega_3,e^{2\pi i/3})$。因此

$$
G_H=(SU(3)\times U(1)_{3Q})/Z_3,\qquad
G_{HS}^{\rm gauge}=(SU(3)\times Z_{18|k|})/Z_3
\quad(k\ne0),
$$

后者分量群为 $Z_{6|k|}$，不是 $Z_{18|k|}$。

$r^a$ 必须同时用弱中心 $(-I_2)^a$ 补偿；只运输超荷角会漏掉该条件。独立得到

$$
K_k=\{(\varphi,a):e^{-18ik\varphi}(-1)^a=1\},
\quad
x=(e^{i\pi/(18k)},1),\quad z=(1,2),
\quad
K_k\simeq Z_{36|k|}\times Z_2^F .
$$

色识别由 $x^{12k}$ 给出；含 $r$ 的内部群分量为 $Z_{12|k|}\times Z_2^F$。再识别 $(-1_{\rm Spin},z)$，由于采用的 Spin 分支连通，完整分量群成为 $Z_{12|k|}$。这不删除费米子宇称。负 $k$ 只使这份生成元的角向取向反转；不得据此宣称两符号完整模型等价。

$k=0$ 保留 $U(1)$，$K_0=U(1)\times Z_2^F$，完整 Spin 识别后连通。作者三种群的计数各属其明确对象，主报告没有混用。

## 3. 质量、局部插值和有限差分

同一标量动能给两个荷向量的正二次型，独立恢复

$$
\det M_k^2=9k^2g^2g'^2v^2w^2,\qquad
M_k^2=C_k^T\operatorname{diag}(v^{-2},w^{-2})C_k .
$$

因而非零 $k$ 的质量及局部横向流秩均为2；$k=0$ 为1。正参数下非零非对角元保证两个中性特征值不同，满秩 $C_k$ 保证各领先极点的残数不消失。这是声明展开阶的关联函数字典，不是全量子粒子谱或未知输入码的认证。

在固定例中，令 $n=k^2\ge1$：

$$
\delta(n)=\frac{18n}{72n+1},\qquad
\delta(n)-\frac{18}{73}
=\frac{18(n-1)}{73(72n+1)}\ge0 .
$$

这独立证明全部非零整数 $k$ 的统一质量下界。对 $s\in[9/584,9/292]$，

$$
R_0(s)\ge\frac12,\qquad
R_k(s)\le\frac{s}{s+18/73}\le\frac19<\frac18,
$$

所以作者保守的 $R_0-R_k\ge3/8$ 成立。上端 $R_0=82/155$。独立逆矩阵还给固定参数下的严格正差，记 $a_0=1/4$、$b=36k^2$：

$$
R_0-R_k=
\frac{b(s+a_0)^2}
{(s+2a_0)\,[s^2+s(2a_0+b)+a_0b]}>0.
$$

每支归一响应总误差至多 $1/64$ 时，剩余差至少 $11/32$；这是条件容差，不是已经取得的高阶或实测预算。$w=0$ 时非零 $k$ 的中性质量矩阵退回零 $k$，明确说明固定尺度下界不可删去。有限精度也不能自行认证精确零质量。

对局部费米插值，纯内部剩余角 $\varphi=2\pi/(18|k|)$ 保 $H,S$ 而使 $N\mapsto-N$，故背景处单 $N$ 及其导数的线性系数为零。此推论依赖作者最终明确的完整背景：规范势、所有普通导数 jets、费米背景均为零。不能扩到有额外带荷背景或参考的菜单，也不禁止成对过程、非局域着装或 Majorana 二次项。

## 4. 物理输入与来源边界

以下限定是通过结论的一部分：

1. 新复 $S$、其荷、径向势、非零相态和正耦合是公开物理输入。$k=0$ 仍是复字段；所选势的相位模、动能及来源保留。树级角向零 Hessian 不能被悄悄当作已消去的自由度。
2. 1060的声明 Spin–$Z_4$ 与 $G_6$ 费米反常结论直接继承；新增普通标量不重写该结论。全背景反常平凡不保证各主丛都有这一无缺陷 Higgs 截面，也不是所有边界／Pin／高形式背景或全 UV 证书。
3. 完整经典来源从同一作用变分，包含全部规范、标量、势、费米和 Yukawa 部门。作者式(21)符号正确；总守恒使用共同场方程，部门能流不必各自守恒。来源定义、经典守恒和实际无界源测试误差不能混为一项。
4. $R_k(s)$ 是平移不变参考附近的领先欧氏／类空场强核。有限正动量窗及其加权积分没有自动给出紧时空探针、统计样本或真实记录。误差若要用于 O 层，须另有匹配；本轮没有将它赋零。
5. 保留无质量 Abelian 通道只在这一固定 SNN 分支选择 $k=0$。原 Dirac 分支、更改字段或质量机制仍属外层比较；共同因果锥本身没有被重解释成光子要求。
6. NEXT正确区分“各模型满足同一角色”与“各模型必须同预测”。本轮没有证明全部分支均具共同操作资格，也不要求一个后继自由见证对所有原观测不可见。

因此本轮可验收为同一有效作用中的结构—响应条件连接。新增的实际仪器、完整热几何共同实例及标准模型唯一性均未签署；它们也不因未在本轮调用就成为本轮必须扩建的技术任务。

## 5. 执行回执

一组独立实现，未重算1060反常或历史仪器。

|执行|回执|结果|
|---|---|---|
|独立 calculate 内存精确预检|1154ae|exit 0，PASS_WITH_SCOPE|
|首次 --write 排他保存|c3fe90|exit 0，CREATED|
|不同工作目录 --check|a452cc|exit 0，READ_ONLY_MATCH|

独立科学 payload SHA256：`91237d84ef965ee762224d29e5d07a514f76ac931eb787f31fb56490fb8ce636`。

作者结果的全部 source_bindings 已逐项对磁盘字节核对一致。未修改任何作者文件、导航或旧证据。下面绑定实际读过的最终版本；尚未生成的 acceptance 不参与循环绑定。

## 6. 最终文件绑定

|根相对路径|SHA256|
|---|---|
|research_cognition_physics/archive_1290_1307/1292/proof.md|`03d79f08a37949a7200b42af0d1091f2c8441dcfe3a4e0e8dc9dcfb9efc2c99a`|
|research_cognition_physics/archive_1290_1307/1292/check.py|`4acff6d5b8688251540ab163833b35b175733f3c6a663c3043ba11a8b233c586`|
|research_cognition_physics/archive_1290_1307/1292/results.json|`f8b3508de5810e9fa9e2d53839ec61871140965147de76d0b0960c3f4462607b`|
|research_cognition_physics/archive_1290_1307/_shared/round1292_admission.md|`6495e115078a6fca4a04a8f9df9a3103afb0e5644e0eaf5d72ea8dd5430b6b3a`|
|research_cognition_physics/archive_1290_1307/research_note_1292.md|`0b4ebe663784f0d988cc2a4c81fc31792b73ad999e87b5490600098992b255fa`|
|research_cognition_physics/archive_1290_1307/1292/NEXT.md|`bf0c84932c2933614299de6f1de9ab3991d24005bd51c08ab53b081db52503d4`|
|research_cognition_physics/archive_1290_1307/1292/independent_check.py|`b10e074132110bbb7955385629d2387d5f1f58c44c2bc2a950bcdce486694d3f`|
|research_cognition_physics/archive_1290_1307/1292/independent_results.json|`c1c247a87ec6733b2218f61c0c8f685a0bbf4076be101cc18df4e710c5a499ea`|
|research_cognition_physics/archive_1009_1043/1035/drafts/joint_selection_derivation.md|`2ef6676648f8960c23ec24b69d0e396ab6e47cd9c20a216898b1d00ceea225ad`|
|research_cognition_physics/archive_1046_1062/1060/proof.md|`a973063186de08b100e72c8a49a824854b0132ffa9d173ed0300865dd2911805`|
|research_cognition_physics/archive_1288_1289/1289/proof.md|`e796fa7da1834a6fbdd7fed29892b888d53c4dd52c1b39fd7f180e30eb3a3b0d`|
|research_cognition_physics/archive_1290_1307/_admission/massless_channel_cognitive_scope_screen.md|`b5148a3451a018ef485bf70c7eb0178263a666f177d0b12f9b5ab92fe9cd7765`|
