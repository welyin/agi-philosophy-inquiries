# 第1293轮独立数学与计算审阅

2026-10-18。结论：**PASS_WITH_SCOPE**。最终作者稿及主线报告的数学链和范围一致，无必修问题。此审阅仅验收本轮明确的二次有效场实例，不代替主线 acceptance，不新增公理，也不宣告长程目标完成。

## 1. 实际审阅与独立性

全文审阅[作者证明](proof.md)、[作者代码](check.py)、[作者结果](results.json)、[事前准入](../_shared/round1293_admission.md)、[本轮报告](../research_note_1293.md)及[后继停止线](NEXT.md)。作者结果中9份来源绑定逐项与当前字节相符。旧730／798／813／818的继承只沿已读回的原范围；1287仅继承 bump 与时样核估计方法，不搬其字段、信号或标签。

[独立代码](independent_check.py)不导入作者或旧轮计算。它以不同CAR位序、相关的忠实Gaussian参考、完整算符块与两相反Fourier单元复核准备、来源和能量。严格连续量词由解析证明承担；有限矩阵和径向积分只作诊断。

- 首次排他保存[独立结果](independent_results.json)：`2c4a1d`，exit 0。
- 异工作目录 `D:/workspace` 执行 `--check`：`ab4593`，exit 0，`READ_ONLY_MATCH`。
- 独立科学 payload SHA256：`5319d26adf9e721bd959097cca4793ef03775ed95be16acc5403da574202ec25`。
- 最大CAR／来源算符残差为 $3.889171018395362\times10^{-15}$；Fourier恒等式最大残差为 $6.892309123117358\times10^{-16}$。这些数值不是连续场近似误差。

## 2. 连续传播和严格常数

由同一 $h_D$ 的平方，跨手性块为 $-im\sin(t\sqrt{-\Delta+m^2})/\sqrt{-\Delta+m^2}$；乘三维KG初速度尾，严格时样核为 $im^2J_1(ms)/(4\pi s)$。该符号与质量平方因子正确。

固定空间／时间窗确给

$$
 r\in[3/16,5/16],\qquad s^2\in[3/32,55/256],\qquad
 \frac{J_1(s)}s\ge\frac{1993}{4096}>\frac{31}{64}.
$$

核在两个正 bump 的全部积分域同相，不能以振荡抵消消掉信号。独立有理核验使用 $e^{4/3}<4$、$e>2$、交错项比值不超过 $55/2048$。它重新核CAR的L2归一而未套KG的L1归一。作者消去 $\pi$ 的较强幅度界成立；独立程序另以 $3<\pi<4$ 得同一公开结论

$$
 |\alpha(t)|>\frac{31}{2^{27}},\qquad
 p_0(t)-p_1(t)>\frac{961}{2^{54}}
 \quad(t\in[7/16,1/2]).
$$

此前因果外零差、删除质量后的跨手性零块，以及事前非负归一时间权平均的继承均成立。删除模型采用自身正常参考，不是对原真空免费质量 quench。该计算不排除其它质量机制或其它传输任务。

## 3. 未知伙伴、真实码态与接收量

固定单占据的两味码处于同一奇宇称扇区，密度矩阵为偶；同参考角的分级乘积因而有效。保的是原参考在补空间的整个限制，原有限模式与补空间的旧关联确被替换，未冒称保留。正常性来自原表示中的有限CAR替换；紧支对象是全Cauchy模式，不是正频投影后仍紧支的单粒子。

独立16维Fock检查从非因子化参考开始，验证替换映射完备性、共同角限制、未知输入与伙伴的保留。接收CAR算符的码块满足同一对角块及零交叉块；这个算符身份给任意伙伴的初始独立性，不靠单个Bell样本概括。混合单占据码确有 $\langle n_1n_2\rangle=0$，而其准自由替代给 $1/4$，因此只对二次源使用协方差的限定必要且已遵守。

共同演化后占据概率差恰为单粒子重叠模平方；接收偶双线性的输入相关部分同样运输相干矩阵元。它是接收代数中的记录关联，不声称未知qubit已无噪转移，也不把数学效果认作已制造的仪器。

## 4. 同一完整二次来源与能量

作用先对 $g,A,H$ 变分再取背景。两选定代际的 $Y_\nu=I_2,Y_e=2I_2$、其余不混入及共同参考不变，使所有已登记背景变分保持味单位结构。物理 $L,\nu_R,e_R$ 取 $U$，全左手记号 $N,e^c$ 取 $U^*$ 的最终表述正确。旧离散 $r$ 给两手插值量同为 $-i$，未知码没有混合不同离散荷。

在来源形式域中，Schur身份与直接压缩均给 $\mathcal E_t^*(J_\zeta)=j_\zeta(t)I_2$。带电弱流／Higgs混合腿的零压缩与其余部门输入无关基线没有被删去；中性Higgs相位源仍在。一般背景Ward须保 $j\,\delta A+\mathcal F\,\delta H$，本背景的平移守恒不推广为各项荷各自守恒。已继承全字段的反常范围，未以两中性领先分量重新冒领全反常证明。

独立有限CAR检查任意味单态二次算符的码块；另用真实复Yukawa变分 $\delta Y=i\operatorname{diag}(1,-1)$ 得非标量压缩 $-\sin(2t)\operatorname{diag}(1,-1)$，而共同Higgs相位方向给 $-\sin(2t)I$。这具体核定：固定耦合背景来源菜单可以标量，任意参数spurion不能被一并纳入。

能量独立展开为

$$
 \operatorname{Tr}(h_D\Delta C_\rho)
 =-2\operatorname{Tr}(Fh_DP_-)
 =\operatorname{Tr}(F|h_D|)
 =2\langle g_A,\sqrt{-\Delta+1}\,g_A\rangle.
$$

其中 $V^\dagger h_DV=0$ 同时消去 $FPF$ 项和未知码项，不是把参考移除成本设零。由 $\max_{z\ge1}z^2e^{-z}=4/e^2<1$ 得 $|\nabla b|<2/a$；内半球下界给 $\|\nabla g\|^2<524288$，故 $0<E_{\rm rel}<1450$。独立两相反Fourier单元的三个未知码密度矩阵复核同一恒等式，最大残差 $1.3322676295501878\times10^{-15}$；这是代数诊断，不是连续谱截断。

## 5. 通过的精确范围

1. 成立于原字段的明确领先二次有效对象、登记初态和有限局域偶效果；Higgs着装完整复合量的全阶CAR身份未取得。
2. 来源码上标量只指每份固定平滑紧支背景测试的期望及其与被动伙伴的联合期望。来源高阶矩、条件化末读、味敏感装置、任意Yukawa变分均未被该结论覆盖。
3. 共同参考／其他字段绝对来源保持为公开基线；相对正常序能源界不等于制备成本、局部正能、弱引力界或平直自引力解。
4. $n=3,4$ 分别承接这份合同，不要求两个模型预测或全部来源相同，也不形成整个相互作用SM的代际自由反模型。
5. 完整SM相互作用、Higgs复合高阶及真实读口误差未赋零。以后若两输出误差为 $\varepsilon_0,\varepsilon_1$，有界效果差只能继承 $961/2^{54}-\varepsilon_0-\varepsilon_1$；无界来源须另有其域预算。

报告和NEXT保留这些边界，并明确停止优化小概率／制造仪器支线，回到有结构选择收益的比较。未见需作者修订的数学或范围阻断。

## 6. 最终版本绑定

以下为实际审阅和复算时的字节哈希，不绑定尚待主线签署的 acceptance，避免循环依赖。

|文件（工作区相对路径）|SHA256|
|---|---|
|`research_cognition_physics/archive_1290_1307/1293/proof.md`|`ed3281d145225f2512d58a664fa8d1ad6160cdab5c1bc4885e9a1daa5f8331db`|
|`research_cognition_physics/archive_1290_1307/1293/check.py`|`d33d18f99dd30aa14b4892b3e21323af4f8ef32a671312bb6158b0f7586c9547`|
|`research_cognition_physics/archive_1290_1307/1293/results.json`|`604010929e2c9b02e9206addd01e2cb4501e4a84b1e07fd4921e92c6088a4070`|
|`research_cognition_physics/archive_1290_1307/_shared/round1293_admission.md`|`3bb3fc23a5393cdea60654728bbb16b4ae45b7622102747c56de9a837d50fa40`|
|`research_cognition_physics/archive_1290_1307/research_note_1293.md`|`bd9008778b0b4d7b77ba3e47c4416a9bfd0cfe77c7300d558da3a1fa7351e654`|
|`research_cognition_physics/archive_1290_1307/1293/NEXT.md`|`a10cd4e7999b7b41e86ed6c18cdc3070c42460b247673f04d843430d7d029b38`|
|`research_cognition_physics/archive_1290_1307/1293/independent_check.py`|`87b27d005d6b72f540ba33abda2536e026eca11caafec30425445310f3f996fa`|
|`research_cognition_physics/archive_1290_1307/1293/independent_results.json`|`e543207f82d45bd8663fa0e50e95cd8b8928b8e9f31561004646ab9e070ff16e`|
