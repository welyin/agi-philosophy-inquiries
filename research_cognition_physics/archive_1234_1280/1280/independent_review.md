# 1280独立审阅：自然共形参考的绝对来源与有限匹配

2026-10-17。结论：**PASS_WITH_SCOPE**。独立核验1组；仅写本轮独立代码、结果和本审阅。未导入作者或1279的实现。

## 1. 审阅内容与成熟来源

已通读[固定规格](fixed_spec.md)、[准入](admission.md)、[作者证明](proof.md)、作者代码／结果、[主报告](../research_note_1280.md)及[NEXT](NEXT.md)，并核对后继观察者筛查仅作候选。当前CT与HG采用状态保持，本轮是推导后场方程的有效实例，不是新的引力生成证明。

实际外部来源核查如下：

- [Pinamonti 0806.0803v2](https://arxiv.org/html/0806.0803)：§2 Lemma 2.4的Hadamard运输，§3.3式(9)及Theorem 3.1的Wick平方共形选择。
- [Brown–Cassidy D16,1712](https://doi.org/10.1103/PhysRevD.16.1712)：出版方摘要；本次未获取该原PDF全文。具体变换式通过[Awad–Johnson hep-th/0008211v3](https://arxiv.org/pdf/hep-th/0008211)第7页式(17)、第9页式(24)独立核对。其参考文献37实际列D15,2810，作者稿已正确区分。
- [Halpern–Neiman 1502.04106v2](https://arxiv.org/html/1502.04106)：§II式(2)—(8)、§IV的null正则场和半面参考；保留无焦散与未处理零频模范围。

成熟场论输入不因本轮数值诊断而成为新定理。以下曲率及有限误差核验针对已固定对象。

## 2. 共同背景、bulk参考与最小来源

**自然坐标。** 在过去线性段 $r=b(u-u_c)$、$b>0$，可直接积分

$$
\eta(u)=\eta(-2)+\frac1b\left(\frac1{r(-2)}-\frac1{r(u)}\right).
$$

因而 $u\downarrow u_c$ 时 $\eta\to-\infty$；未来 $r=1$、$\eta=u$。该最大正支确实共形到整个Minkowski空间，光滑正共形因子保持全局双曲性及Hadamard资格。它没有越过Rosen零点，且不是1279的artanh参考或周期横向模型。

**曲率。** 从度规直接求联络有 $\Gamma^u{}_{ab}=0$，故 $\ell=du$平行；$R_{ab}=2q\ell_a\ell_b$。由于 $\ell^2=0$，$R=0$、$R_a{}^cR_{cb}=0$、$R_{ab}R^{ab}=0$。由 $\nabla_cR_{ab}=2q'\ell_c\ell_a\ell_b$再求导，得到 $\Box R_{ab}=0$。共形平坦还给 $C_{abcd}=0$。因此相关 $H^{(1)},H^{(3)}$、标准曲率平方变分及Bach项均为零；不能只用“标量不变量全零”替代这项张量检查。

**绝对来源。** 共形运输的真空及预先固定的Wick处方给 $\langle\varphi^2\rangle_{\rm ren}=0$；两端的 $R$本身也为零。共形应力的平直项与局部异常项均零，再由

$$
T^{\rm conf}_{ab}-T^{\rm min}_{ab}
=\frac16(G_{ab}+g_{ab}\Box-\nabla_a\nabla_b)\varphi^2
$$

得到原最小参考应力零。这不是从最小／共形KG算符相同直接推出Hilbert源相同。工作相干态仍使用最小作用，故其完整声明标量平均源为 $f'^2\ell_a\ell_b$，与原 $r$ 方程相容。工作态的改进项一般非零，未被省去。

该结论依赖已固定兼容处方以及重整化 $G,\Lambda$。它不许可事后增减 $cG_{ab}$，不推广至任意新增非局部作用；参考本身作为同背景比较态也不被要求单独满足这份有脉冲的半经典方程。

## 3. 局部正常性与有限横向量词

平面位移由光滑实KG解定义，给正代数相干Hadamard态。任一相对紧区域的相应Cauchy数据可以紧支延拓，因而局部正常。这不提供整个平面的有限能量Fock向量。

作者§4的有限 $A_\perp$说法可接受，限于声明的lightfront导数／CCR代数：正频单粒子空间的横向部分为 $L^2(\mathbb R^2)$，分解成 $L^2(A_\perp)\oplus L^2(A_\perp^c)$；真空及半面模作用分解。由于 $1_{A_\perp}\in L^2$且纵向位移光滑紧支，该子代数上的位移正常，响应为面积乘密度。

这只是原态的代数限制及表示：没有把 $1_{A_\perp}$当作一份新的bulk波包，没有周期紧化、实体墙或实际末读，没有将有限null响应等同完整有限资源过程。因此本轮从密度写有限代数响应是严格的，但不能借此领取全局有限总资源。

## 4. 新模权与几何余额的独立核验

对脉冲 $u\in[-2,-1]$，$r(u)\le r(t)\le1$给

$$
r(u)^2\le\frac{-\eta(u)r(u)^2}{-u}\le1,
\qquad \delta_\eta=1-(499/500)^2=\frac{999}{250000}.
$$

$B$是最小工作源的物理boost密度，$M_\eta$是自然参考的模密度，二者独立定义。1279的分部积分式 $C=J-B\ge0$可继承，旧artanh参考却不可继承。精确有理核验为

$$
\frac{|M_\eta-B|}{B}
\le\frac{375}{498002}+\frac{999}{250000}\left(1+\frac{375}{498002}\right)
=\frac{591628623}{124500500000}<\frac1{200}.
$$

严格余量 $30873877/124500500000>0$。有限几何片仍取 $[-2,0]$，$\Delta a-E_g=B/18000$、$s=36000\pi$，故上述同一误差控制模响应与独立几何余额。固有钟曲线及 $1\le d\tau_m/d\tau\le250000/249001$成立；这不是实际观察者KMS或探针误差定理。

## 5. 独立实现与实际回执

[独立代码](independent_check.py)以度规一、二阶jet自行生成联络、Riemann、Ricci及全部Weyl分量，在三个精确Fraction非平坦jet上诊断；没有从作者假定的Ricci形式回填Weyl。有限jet不是通式证明，通式由§2解析链承担。代码另独立实现 $y=2u+3$脉冲、自然 $\eta$积分、两档RK4及三档Gauss求积。

- 首次排他结果：`a54496`，exit 0，`CREATED_ONCE`。
- 异工作目录默认只读：`e5f5e1`，exit 0，`READ_ONLY_MATCH`。
- 细档诊断：$E=6.052923538596655$，$B=9.078700802079652$，$M_\eta=9.078235145658025$，$(M_\eta-B)/B=-5.1291085781776214\times10^{-5}$。
- $J-B-C$残差约 $4.38\times10^{-15}$，有限Raychaudhuri余额残差约 $2.45\times10^{-15}$；$\eta(-2)=-2.0000613483535683$。

作者和独立积分一致。不同采样端点使记录的模权比最小诊断略有不同；二者都不用于证明连续权界。浮点残差、网格差不是严格积分误差、有限仪器误差或完整量子涨落控制。

## 6. 验收范围与停止线

**无必修问题，限定通过。** 本轮补上的是真正bulk参考、该处方下绝对最小参考源及同一工作态平均源自洽，再以有限界接到null模响应和独立几何余额。不是单纯重复相减记账，也不是另一次Einstein选择证明。

主报告／NEXT正确保留了实际读口、加入材料的来源、噪声和反作用合同尚未完成；局部平均源结论未升级为稳定性或涨落定理。[观察者后继筛查](../_shared/jacobson_observer_completion_screen_after1280.md)的未来平坦段核身份只被作为候选，未纳入本轮已认证的实际仪器。全局能量、参考恢复和参数唯一性未领取，也没有被升为所有有限有效实例的新前置门槛；比较参考无须另造制备工厂。

本独立核验为1组；成熟参考理论与旧几何恒等式继承，不重复计数。

## 7. 最终版本绑定

作者和独立结果的已有文件绑定均已逐项读取核对；以下为本次实际审阅版本，后生成的主线验收文件不预先绑定。

| 工作区根相对文件 | SHA256 |
|---|---|
| `research_cognition_physics/archive_1234_/1280/admission.md` | `b4269abd82b794d5c02ea613f7d7259203441a4ce9a2eb640436962fab683754` |
| `research_cognition_physics/archive_1234_/1280/fixed_spec.md` | `8329836937b059188a61825a4bf9ff96c9e706c3addffc854f191122e7cffe7b` |
| `research_cognition_physics/archive_1234_/1280/proof.md` | `760534ff5559be3ab95dc121054b18e9272f534f7d6b393cd6a271b7d8b3e8ec` |
| `research_cognition_physics/archive_1234_/1280/check.py` | `ad2f58a9acbefc0f057c04ebeadcdf290730ff3df5d88e42043e8481d6e42cd4` |
| `research_cognition_physics/archive_1234_/1280/results.json` | `69157a63c980f8ecdc430a3e4154e76a738939352602fe052ddccb5dabaf7f85` |
| `research_cognition_physics/archive_1234_/1280/NEXT.md` | `84607b3787fcff5a0567435cff293266345ea2b2ba55548e90d44db9653a6fd2` |
| `research_cognition_physics/archive_1234_/1280/independent_check.py` | `2b9c58473b7274d423d6ff354fe8b063f2966cefdf5bcea63e71e7991dc5616d` |
| `research_cognition_physics/archive_1234_/1280/independent_results.json` | `ebdfbb810bf3aa28885962a450c7f4b545931438f7397213ee4b6df13adb7d5f` |
| `research_cognition_physics/archive_1234_/research_note_1280.md` | `e827501aeb7bc990e1fbe1c5c4c80efe828b06239c90db890d27a89bf28b0e86` |
| `research_cognition_physics/archive_1234_/1279/proof.md` | `e19d39d32677f49bdd54121b4dbbc908ac5fe917de489918bc3722b918e00523` |
| `research_cognition_physics/archive_1234_/1279/results.json` | `21b80fbc89d1686427faa42bdec7d77020bb3e31203b2159fd451d37e71ff629` |
| `research_cognition_physics/archive_1234_/_shared/ct_hg_adoption_after1278.md` | `b3d92a66dba200c3e5ebad79f9805597a687b316b4b808d6bd166be9e7c33721` |
| `research_cognition_physics/archive_1234_/_shared/jacobson_observer_completion_screen_after1280.md` | `128a8610dff85ca2629fd796f1495c1c6fc5481e24d8ef13f1ddcf91254c8112` |
