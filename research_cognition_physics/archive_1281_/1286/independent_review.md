# 1286 独立审阅：等能未知编码、有限量子记录与共同平均来源

2026-10-18。**PASS_WITH_SCOPE。** 独立核验1组；没有必修项。仅签本轮声明的齐次量子有效类，不将建模输入升级为认知公理，不授予完整4D局域量子仪器或实际末读取权限。

已全文核对[证明](proof.md)、[作者代码](check.py)、[结果](results.json)、[报告](../research_note_1286.md)与[后继](NEXT.md)，并按[公理审计](../_shared/axiom_audit_20261018.md)及[角色分层](../_shared/effective_role_layers_20261018.md)区分认知角色、有效类别、定量输入与条件结论。

## 1. 独立方法与等能编码

[独立代码](independent_check.py)不导入作者函数。作者使用径向幂单项式微分和Gaussian矩；本核验使用正交广义Laguerre径向振子基及三对角Jacobi作用。对各部门维数 $d=3,2,1$，令 $\alpha=d/2-1$，基为 $L_n^{(\alpha)}(v/2)$ 乘Gaussian，则

$$
vL_n^{(\alpha)}=2(2n+\alpha+1)L_n^{(\alpha)}-2(n+1)L_{n+1}^{(\alpha)}-2(n+\alpha)L_{n-1}^{(\alpha)},
$$

$$
\|L_n^{(\alpha)}g_d\|^2=\frac{(d/2)_n}{n!},\qquad
T_d+\frac v8=n+\frac d4.
$$

由此逐次作用完整 $H$；范数平方是正权系数平方和。这里没有压缩传播、数值积分或有限维截断Hamiltonian。每个固定次 $H^k$ 的精确有限展开不等于假定整条传播留在有限基中。

独立得到

$$
\mathcal J^\dagger\mathcal J=I,\qquad
\mathcal J^\dagger H\mathcal J=\frac{4561}{78}I,\qquad
\mathcal J^\dagger Y\mathcal J=\operatorname{diag}(14,266/25).
$$

$H$ 对径向指标的带宽不超过2，$Y$ 不超过1，故相隔3及以上的相关交叉矩元为零。独立对角式为

$$
E_n=\frac{357n^2+1683n+860}{156};\qquad E_3=\frac{16E_0+9E_6}{25}=\frac{4561}{78}.
$$

全部能源分项分别核对，包含中性场动能 $1/8$ 及交叉势的 $7$、$133/25$。等能不是扣除记录资源后的结果；删除交叉耦合时总能及背景需重算，作者已明确这一点。

## 2. 四阶余项与全窗正差

两码从 $H^0$ 至 $H^4$ 的全部范数平方与作者结果逐项精确相同。正交基所得范数整数上界为

$$
(1,79,12463,2925311,968191741),\qquad
(1,124,44613,20742456,12107230522).
$$

对有界 $E_r$，将四阶时间导数直接展开在左右演化态上并用Cauchy–Schwarz，确得

$$
|p_j^{(4)}(\eta)|\le\sum_{k=0}^4\binom4k\|H^k\Psi_j\|\,\|H^{4-k}\Psi_j\|.
$$

$\Psi_j\in D(H^4)$ 足以支撑此式；它不要求无界对易子的全空间算符范数或全时间Taylor级数收敛。两个上界为 $4717138248$ 与 $56732896010$，总和 $M=61450034258$。

独立检查二阶项的符号：$-\partial_rV\,E_r'=Yr\sin r/2$，而两码共享初始记录Gaussian及其自由项，所以

$$
(p_0-p_1)''(0)=\frac{42}{25}e^{-1/2}.
$$

实Hamiltonian、效果和准备使两概率均为偶函数。使用 $e^{-1/2}>1/2$，对整个 $T\le\eta\le2T$ 必须按 $\eta^2$ 因子化控制，而非把两个不同端点的粗界任意拼合：

$$
p_0(\eta)-p_1(\eta)>
\eta^2\left[\frac{21}{50}-\frac{M(2T)^2}{24}\right],\qquad
\frac{21}{50}-\frac{M(2T)^2}{24}
=\frac{1229274982871}{3000000000000}>\frac25.
$$

因此 $T=10^{-6}$ 时全窗严格差 $>1/2500000000000=4\times10^{-13}$ 成立。这同时下界两个实际记录字段边缘的半迹距离，但没有认证达到该精度的实际终末仪器。删掉交叉势后，记录与其余材料的传播因子化，记录边缘对两码恒同；这是有效的删条件对照。

## 3. 未知伙伴与共同来源的量词

完整 $U_\eta\mathcal J\otimes I_R$ 是对任意未知码及伙伴的等距，保持全联合量子信息；没有声称原材料边缘不受扰、未知态被复制或码空间在传播中不变。压缩能源是标量，且 $U_\eta$ 保 $H$，所以该全部输入域共享 $E_*$。所有编码态及其交叉项均保Gauss及规范场三分量旋转，足以保持本有效域的均匀各向同性平均来源。

独立用四份有理度规jet核中性共形项与显示边界项的抵消、混合势尺度变换；再对有效作用的 $N,a$ 分别变分，取得

$$
\rho=E_*/a^4,\qquad P=\rho/3.
$$

来源包含中性自由项、曲率改进及混合势。独立核lapse方程、尺度方程和连续性后，$b^2=4561/1287$、$a=1+b\eta$、$\tau=\eta+b\eta^2/2$ 相容。几何不是由记录概率反定义，任意未知码也没有被分别配上另一背景。

这仅是同一声明有效作用的非选择平均来源。它不使能源涨落为零，不为条件记录分支分配几何，也不证明完整4D量子截断、空间子区域的局域接续或微因果仪器。$r=as$ 的末效果使用共同已定义的尺度字典；若调用实际读取，其控制、误差及来源须另列。本轮没有调用它，故不能把读出费用默认为已包含，也不必为当前量子关联结论另建未调用装置。

## 4. 执行记录与判定

- 独立结果首次排他创建：回执 `3a8202`，exit 0，`EXCLUSIVE_FIRST_WRITE`。
- 作者默认只读复算，经现有matter/gauge runner：回执 `efcc6e`，exit 0，`READ_ONLY_MATCH`。
- 独立默认只读复算，从 `C:/Users/admin` 经同一runner调用：回执 `4bcbb6`，exit 0，`READ_ONLY_MATCH`。

[独立结果](independent_results.json)保留精确范数、能源分账、全窗余项、独立来源变分及限定范围。没有重复旧大轮实验，没有增加模型、参数扫描或新公理。准许按 **PASS_WITH_SCOPE** 签收本轮共同量子记录／均匀来源实例；不将其计为整个阶段已完成。

## 5. 最终字节绑定

以下七件为本审阅实际读取、复算及绑定的版本；作者五件未修改。

|文件|SHA256|
|---|---|
|[proof.md](proof.md)|`a81bb431bd1a4d6e9e96936939536272d3ed9810fdcfebdaf6f011069777c8a4`|
|[check.py](check.py)|`46c4caffa87eb4f64bb475205ca7bdf3c3b8eec1c82dea68fff529b44e51933f`|
|[results.json](results.json)|`3b0aa05c7ff4ca7e085591834b8e5115aa82aaf043d2b862ebf2a28f9ac7cc67`|
|[NEXT.md](NEXT.md)|`69f45b525cd83d64ed8b74b040f9de623294d9582e880d80270d44b172c6e19b`|
|[research_note_1286.md](../research_note_1286.md)|`f4a9bd29432aa4d98c592113a6cc7145ffe9fb29b4be369ca5d7f4e318f153f3`|
|[independent_check.py](independent_check.py)|`47ff1d30f21225f5a9dfd3a5d1093607e253123cec6ac39477578fcb6e866ee4`|
|[independent_results.json](independent_results.json)|`27af768341c83f84fac9927aee207bfbef7c3a8b57d824bd68e7213e7669f73d`|
