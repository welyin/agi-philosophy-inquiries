# 1279独立审阅：同背景类光参考与物理流有限匹配

2026-10-17。结论：**PASS_WITH_SCOPE**。独立核验1组。未导入作者代码；连续结论由解析证明与精确有理数承担，浮点积分仅作交叉诊断。

## 1. 实际审阅对象与来源

已全文核对本轮[准入](admission.md)、[固定规格](fixed_spec.md)、[证明](proof.md)、作者代码／结果、[主报告](../research_note_1279.md)及[NEXT](NEXT.md)。CT和HG采用状态依照既有正式采用记录，不改历史前提。本轮是已经推导的场方程上的共同接口实例，不另计引力生成。

成熟来源实际读取范围为[Halpern–Neiman，arXiv:1502.04106v2](https://arxiv.org/html/1502.04106) §II式(2)—(8)、§IV与§VI。原文支持无焦散类光面上的场重标定、所选单调参数的自由标量参考及半面热性；零频模处理和bulk延拓均有边界。本轮采用这些输入，新增的是同一曲域物理来源的有限比较。没有从该来源领取实际热探针、工作态绝对应力或任意bulk区域资格。

## 2. 独立解析核验

**背景与域。** 对 $ds^2=-2du\,dv+r^2d\mathbf{x}^2$ 和最小无质量标量 $f(u)$，KG方程、$T_{uu}=f'^2$、$R_{uu}=-2r''/r$ 及其余Einstein分量相容；$u$是该类光生成线的仿射参数。横向有限环面和有限 $u$ 域均保留。Volterra首零点论证给全域 $r\ge1/2$，脉冲上 $r\ge499/500$；没有将解扩张为无限无焦散过去。

独立地，用 $y=2(u+2)-1$、$t=y/\sqrt{1-y^2}$ 得

$$
E=8\int_{\mathbb R}t^2(1+t^2)^{3/2}e^{-2t^2}dt
\le\frac{55}{8}\sqrt{\frac\pi2}<\frac{165}{16}<36.
$$

两半区Cauchy界给 $E\ge4$。此处独立上界仅交叉核验，正式预算仍用预声明的 $E\le36$。

**相干参考与物理源。** $\psi=r\varphi$ 的辛形式标准化；$w=500\operatorname{artanh}(u/500)$ 将整段无焦散域映到实线。$rf$ 在 $w$ 中光滑紧支，因此所用无零频类光Fock位移范数有限。负半面完整模响应可因该位移支撑缩写为脉冲积分，不能据此将任意有限小段的模Hamiltonian认作同一表达式。

对 $B,J,M$ 使用作者固定定义，独立展开平方并分部积分得到

$$
J-B=\int[rr'+(-u)qr^2]f^2du\ge0.
$$

紧支端点项消失，且 $\sup f^2\le E/4$。因脉冲长度1、$-u\le2$、$r'\le A_q$，有 $C\le3A_qE/4$。与 $B\ge(499/500)^2E>0$ 合并可合法除以 $B$。另一方面，积分形式的 $\operatorname{artanh}$ 上下界给 $0\le1-\xi/(-u)\le1/62500$。因此

$$
\frac{|M-B|}{B}
\le\frac{375}{498002}+\frac1{62500}\left(1+\frac{375}{498002}\right)
=\frac{23935877}{31125125000}<\frac1{1250}.
$$

严格余量是 $964223/31125125000$。这同时计入曲域场重标定与有限参考区间两项误差，不是直接设置模流等于物理能流。

**几何与钟。** 在几何片 $[-2,0]$（不同于脉冲支持 $[-2,-1]$）独立复核有限Raychaudhuri余额：$\Delta a-E_g=B/18000$。由于 $s=36000\pi$，$s(8\pi G)=2\pi$，模响应与几何余额使用同一个相对界。给定纵向双曲线确有固有参数 $\tau$ 和几何加速度 $\kappa$；参数率满足 $1\le d\tau_m/d\tau\le62500/62499$。这不是将整个类光模自同构等同为已实现的时样观察者热演化。

## 3. 实际独立复算及回执

[独立代码](independent_check.py)只使用标准库、NumPy及声明函数，自行以 $y$ 变量实现脉冲、反向RK4和Gauss积分；未调用或导入作者实现。严格分数部分另重建常数。两档RK与三档Gauss均为预声明诊断。

- 首次排他保存：`25285c`，exit 0，`CREATED_ONCE`。
- 从工作区 `scripts` 子目录默认只读重放：`dfac20`，exit 0，`READ_ONLY_MATCH`。
- 8192步诊断：$E=6.052923538596655$，$B=9.078700802079652$，$M=9.078739913053282$。
- $(M-B)/B=4.307992352888211\times10^{-6}$；$J-B-C$ 残差约 $4.38\times10^{-15}$，几何余额残差约 $2.45\times10^{-15}$。

这些与作者独立结果一致，但浮点残差、网格差均未被当作严格积分误差、完整量子误差或实际观测容差。严格证书不依赖其小数值。

## 4. 通过范围与未领取内容

通过的是同一固定背景、同一所选类光参考、同一正常相干扰动的**参考相减来源—模响应—独立几何余额**有限匹配。来源相减使用同处方；没有把不同背景的参考相减拼接，也没有把参考的绝对应力设零。

主报告与NEXT准确保留：参考作为比较态不要求另建准备工厂；几何固有时间线已经存在，不能再说“没有钟”。尚需具体有效合同的是工作态完整来源的基线及反作用、同一有限观察者菜单的响应—模流匹配。这里未认证实际仪器、全部工作态来源或整体共同目标；也未增设全球bulk构造／全UV理论作为普遍必要门槛。

没有发现必修数学或范围问题。本轮数量为一组作者核验和一组本独立核验；成熟类光定理及旧有限Raychaudhuri身份按继承处理。

## 5. 审阅版本绑定

以下哈希由实际工作区文件读取生成；作者和独立结果中的已有绑定也已逐项核对。根端后续生成的归档／验收文件不在本次预先绑定之内。

| 工作区根相对文件 | SHA256 |
|---|---|
| `research_cognition_physics/archive_1234_/1279/admission.md` | `ffa5c5374df39df8a0c78bc06d5238f523c9591be4f7379111132dde0b9bbd69` |
| `research_cognition_physics/archive_1234_/1279/fixed_spec.md` | `c69e5171068701c2889dc1afd5288a3b1e037c75218c9c5b0ef522143b799987` |
| `research_cognition_physics/archive_1234_/1279/proof.md` | `e19d39d32677f49bdd54121b4dbbc908ac5fe917de489918bc3722b918e00523` |
| `research_cognition_physics/archive_1234_/1279/check.py` | `f9b5768a62e1000f68d27ac1f53ddf70683e9e7eab96d1b4bf73b0e8c223b2fa` |
| `research_cognition_physics/archive_1234_/1279/results.json` | `21b80fbc89d1686427faa42bdec7d77020bb3e31203b2159fd451d37e71ff629` |
| `research_cognition_physics/archive_1234_/1279/independent_check.py` | `2d60bad98a504233b71691634835476ff9bd3c6100352e4b140bcc0df509805a` |
| `research_cognition_physics/archive_1234_/1279/independent_results.json` | `f6fca6626271b6690bfc8e497d0b87ce324576eefc5e41f0a2ce14349ac11698` |
| `research_cognition_physics/archive_1234_/1279/NEXT.md` | `d406b7b5859f98be1edd61113806f60e4ee51866385538cf052a22ac8998d6ac` |
| `research_cognition_physics/archive_1234_/research_note_1279.md` | `0135ded12f379db52abf3cfaf99e237a38841ab316a5bb117e43c439e76de4c9` |
| `research_cognition_physics/archive_1234_/_shared/ct_hg_adoption_after1278.md` | `b3d92a66dba200c3e5ebad79f9805597a687b316b4b808d6bd166be9e7c33721` |
