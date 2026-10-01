# 第581轮：共同几何与物质有效项的条件压缩

日期：2026-10-01。接[580量子有效项](research_note_580.md)与[条件账](unified_physics_condition_ledger_580.md)。链接按archive_231_解析，完成以独立终审及冻结为准。

## 1. 实际推进：增加算符不等于增加同样多的独立条件

580已经证明固定场变量的原二导数作用不足，但保留了场重定义的可能。本轮沿“先合并”的要求，检验允许同一几何参与变换时，能否把当前标量一圈局部系数压缩为更少的算符结构。

**可以在明确的一阶EFT范围做到。** 当前b4可由显式度规与标量重定义，转写为Euler项、两个四梯度组合、二导数的US项与一份共同势修正。给出的是同一作用与观察量的表示转换，不是删掉所有新物理条件；两个剩余四梯度项也未完成全部冗余分类。

本轮继承原K、U，额外明确采用动态Euclidean Einstein—标量领先作用作为冗余判据。引力作用、四维、归一和EFT展开是输入；不把这一步当引力生成。只处理580的零规范曲率标量部门，不能替代572的非零规范共同来源。

## 2. 成熟对接、模型与阶数

[Criado与Pérez-Victoria，1811.09413，第2—5节](https://arxiv.org/html/1811.09413)说明局部微扰场重定义及源的共同变换，并区分一阶EOM消元与高阶等价。本文给当前b4与原物质的具体系数、变量替换及余项检查，不将成熟等价定理当新发现。

|层次|输入或结论|
|---|---|
|认知动机|物质、共同几何与记录应有一致表示，不能逐项虚增认知原则|
|领先作用|四维Euclidean动态度规、Einstein归一κ=1、574原K和U；不量子化引力|
|量子输入|580已算的标量一圈局部系数，不含规范／费米子／引力环|
|本轮增量|完整局部系数的离壳因式分解、明确一阶场重定义、原作用实际替换与长度泛函拉回|
|范围|有限匹配系数、局部有界光滑资料、F远离0、同一EFT阶；全导数有明确边界处理|
|未完成|全部独立物理算符、非零规范源、全框架量子匹配、原图连续极限及实验预测|

记λ为微扰阶的记账量，使用有限重整化局部系数或形式级数，不将数值λ直接取作发散的1/ε：

$$
I[g,\phi]=I_0[g,\phi]+\lambda\int\sqrt g\,b_4+O(\lambda^2),\qquad
I_0=\int\sqrt g\left[-\frac R2+\frac S2+U\right],\quad
S=K_{ab}\partial_\mu\phi^a\partial^\mu\phi^b,\quad Q=B_{\mu\nu}B^{\mu\nu}.
\tag{1}
$$

这不是580冻结度规标量行列式忽然包含引力环，而是在所得局部有效作用中指定动态几何的领先方程，检验哪些表示冗余可合并。

## 3. 全部曲率项的精确离壳分解

定义度规误差张量及单倍Euler导数：

$$
E_{\mu\nu}=R_{\mu\nu}-B_{\mu\nu}-Ug_{\mu\nu},\quad e=R-S-4U,\qquad
H_{\mu\nu}:=\frac1{\sqrt g}\frac{\delta I_0}{\delta g^{\mu\nu}}
=-\frac12\left(E_{\mu\nu}-\frac e2g_{\mu\nu}\right).
\tag{2}
$$

H不是580的正二倍响应。以下所有等式不要求E=0；用E留下完整离壳差。四维Euler密度保留为独立拓扑项：

$$
\mathcal E_4=\mathrm{Riem}^2-4\mathrm{Ric}^2+R^2,\qquad
5\left[\frac{\mathrm{Riem}^2-\mathrm{Ric}^2}{180}+\frac{R^2}{72}\right]
=\frac{\mathcal E_4}{36}+\frac{\mathrm{Ric}^2}{12}+\frac{R^2}{24}.
\tag{3}
$$

不据四维拓扑性质删除它的量子反常、高阶维数正规化作用或有边界贡献。将580 b4拆成一份候选约化密度：

$$
b_{4,\mathrm{met}}=\frac{\mathcal E_4}{36}-\frac{7S^2}{216}+\frac{11Q}{108}
+\frac{US}{18}+U^2-\frac23U\operatorname{tr}A
+\frac12\operatorname{tr}A^2-\frac16\operatorname{tr}(AM_X).
\tag{4}
$$

这里A、M_X仍是580的目标协变Hessian及梯度矩阵。定义

$$
D_{\mu\nu}=\frac{R_{\mu\nu}+B_{\mu\nu}+Ug_{\mu\nu}}{12}
+g_{\mu\nu}\left[\frac{R+S+4U}{24}-\frac{\operatorname{tr}A}{6}-\frac S9\right],\qquad
T_{\mu\nu}=-2D_{\mu\nu}+(\operatorname{tr}_gD)g_{\mu\nu}.
\tag{5}
$$

逐项平方差给精确恒等式：

$$
b_4-b_{4,\mathrm{met}}=E_{\mu\nu}D^{\mu\nu}=H_{\mu\nu}T^{\mu\nu}.
\tag{6}
$$

证明：Ric²差为E:(Ric+B+Ug)，R²差为e(R+S+4U)；其余R相关差为−e(trA/6+S/9)。再由式(2)有E=−2H+(trH)g，得到第二个等号。式(4)中−7/216、11/108与1/18来自同一组系数，不能任意分别选择。

## 4. 物质EOM和总导数也须共同保留

设τ(φ)为目标协变张力，即拉回连接定义的∇^μ∂_μφ；Eφ=grad_K U−τ是领先标量Euler导数抬指标形式。原势的链式法则为

$$
\Box U=\operatorname{tr}(AM_X)+\langle\operatorname{grad}U,\tau\rangle_K,
\qquad
-\frac16\operatorname{tr}(AM_X)
=\frac{|\operatorname{grad}U|_K^2}{6}
-\frac{\langle\operatorname{grad}U,E_\phi\rangle_K}{6}-\frac{\Box U}{6}.
\tag{7}
$$

于是可以定义完整约化系数：

$$
\boxed{b_{4,\mathrm{red}}=\frac{\mathcal E_4}{36}-\frac{7S^2}{216}
+\frac{11Q}{108}+\frac{US}{18}
+\left[U^2-\frac23U\operatorname{tr}A+\frac12\operatorname{tr}A^2
+\frac{|\operatorname{grad}U|_K^2}{6}\right],}
\quad
b_4=b_{4,\mathrm{red}}+H:T-\frac{\langle\operatorname{grad}U,E_\phi\rangle_K}{6}-\frac{\Box U}{6}.
\tag{8}
$$

方括号只依赖φ，可并入共同势修正；US/18可由δK_ab=U K_ab/9表示为动能修正。它们受同一原U约束，不是重新挑选的几份物质来源。一般有限匹配仍可有额外自由参数，本文只压缩所列局部系数的表示。

## 5. 实际变量替换及可逆范围

令右侧全部使用新变量计算，做

$$
g_{\mathrm{old}}^{\mu\nu}=g_{\mathrm{new}}^{\mu\nu}-\lambda T^{\mu\nu},\qquad
\phi_{\mathrm{old}}^a=\phi_{\mathrm{new}}^a+\frac\lambda6(\operatorname{grad}_K U)^a.
\tag{9}
$$

I0的一阶变化为−λ∫√g H:T+λ∫√g〈Eφ,gradU〉/6，恰消式(8)的EOM项。对无边界周期域或恰当边界处理，得到

$$
I[g_{\mathrm{old}},\phi_{\mathrm{old}}]
=I_0[g_{\mathrm{new}},\phi_{\mathrm{new}}]
+\lambda\int\sqrt{g_{\mathrm{new}}}\,b_{4,\mathrm{red}}+O(\lambda^2).
\tag{10}
$$

这是微扰一阶等价；不能只用领先EOM反复替换而宣称全阶相同。T含Ricci，因此逆是导数展开中的局部形式逆，不是完整高阶PDE的全局可逆定理。若在紧域保持F≥Fmin>0、正度规一致界、所需有限阶导数有界，且λT和λgradU足够小，则作用的Taylor余项有依赖这些界的Cλ²控制。它不覆盖无界细化频率、F边界、强曲率或任意场族的一致极限。

## 6. 几何与记录也必须拉回

任意旧观察量须改写为O_new=O_old∘变换，外源耦合、准备与边界资料亦需同样映射。不能删完作用项后依然把未经变换的新度规读数当成原读数。局部非奇异重定义的量子Jacobian及源处理遵循所用正规化；本轮未由此证明两种部分量子化的完整等价。

例：一维坐标闭路上旧度规线长，在新背景g_new=δ时为

$$
L_{\mathrm{old}}=\int_0^{2\pi}\sqrt{g_{\mathrm{old},11}}\,dx
=2\pi+\frac\lambda2\int_0^{2\pi}T_{11}\,dx+O(\lambda^2).
\tag{11}
$$

代码的原势周期背景给一阶系数−.462515160632。忽略观察量拉回就会丢掉一阶效应；这只是清楚的几何泛函见证，不声称已造出完整量子钟尺仪器。

剩余四梯度组合不能仅看某个系数符号判稳定性。Euclidean正B有S²/4≤Q≤S²：

$$
-\frac{7S^2}{216}+\frac{11Q}{108}
=\begin{cases}5S^2/72,&Q=S^2,\\-S^2/144,&Q=S^2/4.\end{cases}
\tag{12}
$$

符号变化不证明原理论失稳；表示、领先约束、其余作用项与观察量均已改变，也不能把截断高导数作用当全阶非扰动理论。两个结构是否还可进一步冗余，仍须完整基底／振幅或观测分析。

## 7. 同一原势的实际复算

[代码](joint_effective_operator_reduction.py)汇总两份入口、[结果](joint_effective_operator_reduction_results.json)保存完整证据，默认全部复算并比较JSON。五组检查首次纳入编号计数；未重算旧轮作为新发现。

实际替换不以线性预期生成新作用值：直接计算φ_old的原U、K与新导数，以及g_old的曲率。对仅依赖x¹的正度规ds²=e^{2a}dx₁²+e^{2b}(dx₂²+dx₃²+dx₄²)，代码使用

$$
R=-6e^{-2a}[b''+2(b')^2-a'b'],\qquad
\Delta_K U=K^{ab}\partial_a\partial_bU-\frac F{3M}\phi^a\partial_aU.
\tag{13}
$$

后一式由574测度√detK=√M F⁻³及K逆的散度展开，另经独立有限差分散度核验；导数不是只用一维势替代五分量目标。

|检查|结果|
|---|---|
|完整b4与离壳因式分解|三个一般局部张量及一个满足领先度规关系的张量，误差≤1.12×10⁻¹⁶|
|原U目标链式关系|原势径向测地线周期族，分部积分与全目标梯度／Eφ收缩差约2.24×10⁻¹¹|
|RS块系数移动|随机不变量多项式误差≤2.23×10⁻¹⁶，单项去除会联动S²与US|
|实际联合重定义|λ=.04、.02、.01、.005时，一阶扣除余项约−5.6050×10⁻⁵、−1.4149×10⁻⁵、−3.5544×10⁻⁶、−8.9077×10⁻⁷；相邻比3.9615、3.9806、3.9903，符合二阶控制|
|长度观察量拉回|一阶系数−.462515160632；λ=.005时拉回长度6.28087448，而未变换裸长度为2π；二阶余项约1.75×10⁻⁶|

这些具体样本是解析Taylor和张量恒等式的交叉核验，不承担任意场族的证明。

## 8. 本轮减少什么，下一步接哪里

[581条件账](unified_physics_condition_ledger_581.md)将580中一批曲率和混合项合并为同一EFT表示转换。没有增加一套新认知公理，也没有把算符数减少等同于已决定全部有限参数。580固定变量的离壳结论和本轮允许改变变量的条件性约化并不冲突。

下一步优先回到共同来源：572保留非零规范应力，此时Ric=B+Ug不再是完整领先度规关系，580的零规范曲率标量系数也不再是全部标量环。需要同时接入规范背景的目标连接与应力，核压缩是否产生不能忽略的物质—规范交叉项。不能在纯标量部门把条件收束后，就称原非零规范共同源已经统一。

553与574的量子变量匹配、连续尺度、物理预测及统一目标继续开放。本轮未解决578—579的非微扰图能源障碍。
