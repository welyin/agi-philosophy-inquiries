# 有限正能观察者熵引理：独立数学与范围审阅

2026-10-16。审阅者：`curved_contact_scope`。

**结论：PASS_FINITE_INNER_ACTION_ENTROPY_ONLY。** 主稿的精确恒等式、非交换态的参考系统处理、Gamma 误差界及硬能窗版本均通过独立解析核对。适合作为本项目的非编号准入引理保存；本次不预约 1271，不增加科学组，不授予 JSS 局域场、实际观察者或几何熵资格。

## 1. 独立解析核验

### 1.1 参考迹与正能 corner

令 $w(y)=(\beta/Z)e^{-\beta y}$。有限材料的简并度已经由矩阵迹及 $Z$ 计入，直接积分得 $\tau(P)=1$。对于材料参考 $\sigma$ 及独立正常波包 $f_\beta(q)=\sqrt\beta e^{-\beta q/2}\boldsymbol1_{q\ge0}$，确有

$$
f_\beta(y-H)\sigma f_\beta(y-H)^*=w(y)P(y).
$$

因此所选迹确实是一个事前指定参考态在该代数上的限制，未依待解释工作态重新拟合。它仍是带中心有限代数上的一项参考选择；$\tau(P)=1$ 不使这份迹成为抽象代数唯一允许的迹。归一有限迹下 $S_\tau\le0$ 的符号正确。

### 1.2 任意非交换材料态与完整纯化

置 $r(y)=\operatorname{Tr}D(y)$、$\rho_y=D(y)/r(y)$。独立展开密度对数给

$$
S_\tau(d)=h(Y)+\int r(y)S(\rho_y)dy+\log\beta-\log Z
-\beta\big(\langle q\rangle+\operatorname{Tr}\rho H\big).
$$

完整纯化使每个 $y$ 条件下的材料—参考向量为纯态，且参考平均态保持原边缘，故

$$
\int r(y)S(\rho_y)dy=S(\rho)-I(Y:R).
$$

能量标签的 cq 延拓满足 $Y\leftarrow J\rightarrow R$；这对简并能层也成立。因此

$$
h(Y)-h(q)=I(Y:J),\qquad
\epsilon=I(Y:J)-I(Y:R)=I(Y:J\mid R)\ge0.
$$

这证明主稿精确恒等式及 $\epsilon\le I(Y:J)$，没有先对 $\rho$ 去相干。对易态的 $\epsilon=0$、纯态的 $\epsilon=I(Y:J)$ 均正确。这里的 $J$ 是证明使用的 cq 标签，未实施一次免费能量测量或生成免费记录。原伙伴若不是完整纯化，不能直接用其边缘取代公式中的 $R$。

### 1.3 Gamma 与硬能窗

对非负能量移位 $a$，直接密度比给

$$
D(p(\cdot-a)\Vert p)=a/\vartheta-(k-1)\mathbb E\log(1+a/Q).
$$

对数二阶下界及 Gamma 两个逆矩使线性项恰好抵消，剩余项为 $a^2/[2(k-2)\vartheta^2]$。再由混合分布的相对熵分解获得主稿的状态相关二阶矩界。该证明也给出谱宽 $\Delta$ 下对全部材料态统一的充分上界 $\Delta^2/[2(k-2)\vartheta^2]$；不用按未知态逐一改变准备。

硬窗版本保持原代数和原迹，仅把状态条件化到 $Q\le B$。以未截断 Gamma 密度作共同比较，额外常数恰为 $-\log P_B$；上尾截断使 $\mathbb E_cQ^{-1}\ge\mathbb EQ^{-1}$，线性项非正，而 $\mathbb E_cQ^{-2}\le\mathbb EQ^{-2}/P_B$。主稿的硬窗界因此成立。若另改代数为双边能窗 corner，其迹应为 $1-e^{-\beta B}$；主稿已正确避免混用。

原平滑波包满足 $\|f'\|_2^2=1/[4(k-2)\vartheta^2]$。硬窗波包一般不在 $H^1$，本引理只需其 $L^2$ 正常性；它不提供有限共轭动量方差或平滑制备权限。

## 2. 实际复算与数值范围

已全文审读证明和两个程序，未导入作者数值结果作为解析证明，也未另造第二份数值实现。

- 主程序默认只读复算：启动回执 `a3952c`，完成回执 `2ce1ab`，exit 0，30 例 `READ_ONLY_MATCH`。两阶求积差 $1.6296317593977871\times10^{-13}$，恒等式残差 $4.496986008399595\times10^{-14}$。最终主程序只更新参考纯化的注释，数值内容未变。
- 硬窗程序最终默认只读复算：`260dcf`，exit 0，12 例 `READ_ONLY_MATCH`。两阶求积差 $3.694822225952521\times10^{-13}$，恒等式残差 $1.441788312003811\times10^{-14}$。
- 最终五文件读回与哈希回执：`0eec0e`，exit 0。

这些是数值稳定性诊断。有限节点、有限积分上端及二进制浮点比较不构成严格求积、尾积分或舍入误差证书；一般量词由第 1 节解析推导承担。

## 3. 净增与停止线

本地查重未发现同一“完整纯化的精确熵余项＋正能 Gamma 资源界＋硬能窗”结果。JSS 的慢观察者近似、旧内作用/参考共轭审计及有限编码熵恒等式是前置工具，不能将它们等同于本式的显式余项。本审阅不宣称全体文献中的首创。

新增内容是一项有限材料代数中的直接熵误差控制。它避免以小迹距替代熵连续性，也表明所选模型的任意指定有限误差可以用有限观察者能量预算取得。资源预算是充分条件，不是最优必要下界；硬窗态的能量上界也不是整套制备和控制装置的资源上界。

以下身份仍未由此计算得到：

- III 型局域场交叉积或因子结构；原 1270 齐次量子块到局域场的运输。
- 该参考及观察者的实际制备、内部固有时、完整来源与因果访问。
- 同一局部引力约束所给的面积/Wald 项，以及该项与这里代数熵的物理身份。

初始观察者与材料—参考独立、所选代数及参考迹是明确模型输入。dressing 是所用表示/限制的构造，不能当作免费物理门。当前结果不采用 TR 或 CT 新角色，不改变既有目标。下一步应核一份具体同作用场/观察者字典，而不是继续仅优化 Gamma 常数或凭这份有限代数公式宣布几何热身份成立。

## 4. 最终输入绑定

| 输入 | SHA256 |
| --- | --- |
| [finite_observer_entropy_candidate_after1270.md](finite_observer_entropy_candidate_after1270.md) | `c0b64a581211ef289e30caaf7a78da4290248b8e4ab56ef2030d5192a9b8dd1b` |
| [finite_observer_entropy_check_after1270.py](finite_observer_entropy_check_after1270.py) | `dcda3baa7681f5c37bd1fda6524fdab78656fc366a2ddc51a7c6d5f32d6eee01` |
| [finite_observer_entropy_check_after1270.json](finite_observer_entropy_check_after1270.json) | `f93f7abe9aa6ffe3b3d01c518a1b576d02b0eb3c359f5c686ef47a8ec755e45d` |
| [finite_observer_hard_window_after1270.py](finite_observer_hard_window_after1270.py) | `89b8c94d5d291b2398a2cf4952a05771454463f18ea77300c11740791fd4b276` |
| [finite_observer_hard_window_after1270.json](finite_observer_hard_window_after1270.json) | `00518136f8d4f95f15d80c4029c8cdb0727065c6f10bf4c885713447342e236e` |
| [README.md](README.md)（工作稿名 `positive_observer_README_after1270.md`，发布时原样改名） | `c9ab19220ccbfa18d209afd7d58fd3f5ebd6c4df60b8c2358801d32eba70872b` |

## 5. 最终文稿确认

**PASS：最终主稿及概览的表述、数值和下一步范围通过。** `4e2bf9` 全文核对概览与微分熵补充，`6dc484`、`74c9fc` 确认最终措辞及哈希。正文现明确 $h(q)$ 是能源分布的微分熵，不能当作纯观察者波包的 von Neumann 熵；$J$ 只是证明标签，$Y$ 是所选代数中的组合量，未获实际读取资格。概览三个数值上界均向上取整，实际余项和均能与已复算结果相符。

硬能窗断言适用于工作观察者准备；实现所选迹的指数参考是另一个指定的正常参考态，本文没有把它也认证为同一硬能窗内的准备。此区分不改变任何熵恒等式。

概览正确保留 CT 选后的一致性用途，同时没有把它冒充 Einstein 的独立生成；保留 1258 局部场部门及 1270 有限共同过程，但明确零模不得双计、FLRW 不能只换参考名便成为精确 de Sitter。停止继续优化 Gamma 参数，转向具体同作用对象的有限字典，符合本次准入范围。数值源码和结果未再改变，本次只追加文稿确认与六文件最终绑定，不重复运行、不增科学组。
