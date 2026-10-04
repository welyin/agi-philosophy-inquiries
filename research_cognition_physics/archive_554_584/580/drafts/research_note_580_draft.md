# 第580轮：同一曲目标的量子有效项与几何响应

日期：2026-10-01。接[579残余能源](../../research_note_579.md)、[574原曲目标](../../research_note_574.md)及[553曲率运行](../../../archive_531_553/research_note_553.md)。稿件链接按archive_231_解析；正式完成以独立终审及冻结核验为准。

## 1. 本轮解决的共同接口

578—579已经限制原有限图量子处方的直接细化，单一能源减项不足。本轮转向成熟的背景场方法，检查**同一H⁵物质在同一固定Einstein几何中的标量量子修正，需要什么局部作用项，以及这些项怎样响应几何变化**。

得到两个限定结论：

1. 在指定四维Euclidean标量一圈处方及固定场变量中，仅重定义原势U、二导数目标度量K及纯几何项，不能吸收所算标量行列式。原目标曲率固定一个非零四梯度组合，并出现曲率—梯度混合项。
2. 混合项在平坦背景的取值为零，度规变分却一般非零；丢掉它会丢失共同几何来源的一部分。固定背景的能源或作用数值相等不够。

这是当前共同候选的有效项与几何响应连接，不是一般非线性sigma模型重整化的新发现。没有把本计算当原图的重整化极限，也没有证明全模型所有环合计后仍有相同净系数。

## 2. 输入、成熟来源与历史去重

|层次|地位|
|---|---|
|认知动机|物质、记录及几何须使用同一量子处方和资源来源|
|继承模型|574的五实分量目标K、原势U、F>0分支；不加χ探针|
|附加处方|给定四维Euclidean背景，一圈、协变背景场分裂与标量测度；只积分五标量，冻结Einstein度规与规范背景|
|本轮范围|规范曲率为零的局部部门，无边界或全导数可丢弃的条件；适当IR处方下的局部UV系数|
|解析增量|当前目标的完整标量系数、周期离壳非吸收见证、度规响应与冻结变量的兼容边界|
|尚未完成|全规范／费米子／引力环、在壳独立算符基、有限匹配、原图连续极限、Lorentz观测及量子Einstein解|

[Alonso–Jenkins–Manohar，1511.00724，第III节式(73)、(78)—(79)](https://arxiv.org/html/1511.00724)给协变sigma模型涨落方法；[Vassilevich，hep-th/0306138，式(4.28)及(2.31)](https://arxiv.org/pdf/hep-th/0306138)给Laplace型热核与Euclidean行列式约定。本轮只做当前K、U与背景处方的映射和收缩。

553已研究纯时空曲率反项，故本轮不重新宣布R²或W²为新结构。553固定Jordan度规；574及本轮固定Einstein度规。两者的场依赖变换须在量子变量中明确，见第7节。离壳背景不必使D全谱为正；IR零模／负模影响完整行列式定义，不能据局部UV系数声称全部背景稳定。

## 3. 同一目标的Jacobi算符

记φ为原五实分量，r²=h²+s²。采用574的M=2、F=M−r²/6及原正定L、真空参数u：

$$
K_{ab}=\frac{\delta_{ab}}F+\frac{\phi_a\phi_b}{6F^2},\qquad
U=\frac{d^\top Ld}{4F^2},\quad d=(h^2-u_h,s^2-u_s)^\top,\qquad
I_0=\int\sqrt g\left[\frac12K_{ab}\partial_\mu\phi^a\partial^\mu\phi^b+U\right].
\tag{1}
$$

目标是五维H⁵、截面曲率k=−1/6；这是场空间维数，时空四维另为输入。在目标及Euclidean时空正交标架记X_mu为背景梯度、A为U的目标协变Hessian：

$$
M_X=\sum_\mu X_\mu X_\mu^\top,\quad S=\operatorname{tr}M_X,\quad Q=\operatorname{tr}M_X^2,
\qquad P=A-k(SI-M_X),\qquad
\Omega_{\mu\nu}=k(X_\mu X_\nu^\top-X_\nu X_\mu^\top),\qquad D=-\nabla^2+P.
\tag{2}
$$

这里M_X不是式(1)的常数M，秩至多4。采用正曲率球的Riemann符号；负k给非负Jacobi动能势，提供符号检查。对目标维数m直接收缩：

$$
\operatorname{tr}P=\operatorname{tr}A-k(m-1)S,\qquad
\operatorname{tr}P^2=\operatorname{tr}A^2-2k(S\operatorname{tr}A-\operatorname{tr}AM_X)
+k^2[(m-2)S^2+Q],\qquad
\sum_{\mu,\nu}\operatorname{tr}\Omega_{\mu\nu}^2=2k^2(Q-S^2).
\tag{3}
$$

最后一项是反对称矩阵的平方迹，通常非正，不是Frobenius范数平方。μ、ν双重求和已含两个方向。

## 4. 热核密度与极点分别定义

定义去掉整体(4π)⁻²的局部系数b4；略去总导数时：

$$
\operatorname{Tr}e^{-tD}\sim(4\pi t)^{-2}\int\sqrt g\,[5+t b_2+t^2b_4+\cdots],\qquad
b_4=\frac12\operatorname{tr}P^2+\frac1{12}\operatorname{tr}\Omega_{\mu\nu}^2
-\frac R6\operatorname{tr}P
+5\left[\frac{\mathrm{Riem}^2-\mathrm{Ric}^2}{180}+\frac{R^2}{72}\right].
\tag{4}
$$

源文D=−(∇²+E)，所以E=−P。代m=5、k=−1/6得到

$$
\boxed{b_4=\frac12\operatorname{tr}A^2+
\frac{S\operatorname{tr}A-\operatorname{tr}(AM_X)}6
+\frac{S^2}{27}+\frac Q{54}
-\frac R6\operatorname{tr}A-\frac{RS}{9}
+5\left[\frac{\mathrm{Riem}^2-\mathrm{Ric}^2}{180}+\frac{R^2}{72}\right].}
\tag{5}
$$

R是时空曲率，不是目标R_K=−10/3。前两项分别修正势及二导数部分；R trA可要求场依赖非最小曲率项，亦不能擅自当常数Einstein系数。即使允许这些扩展，四梯度与RS仍超出原二导数算符类。

采用Euclidean Γ1=½Tr logD、d=4−2ε及proper-time约定，实玻色因子给

$$
\Gamma_{1,\mathrm{div}}=-\frac1{2(4\pi)^2\epsilon}\int\sqrt g\,b_4,
\qquad I_{\mathrm{ct}}=+\frac1{2(4\pi)^2\epsilon}\int\sqrt g\,b_4.
\tag{6}
$$

单个平直自由实标量P=m_s²有b4=m_s⁴/2，故极点−m_s⁴/(64π²ε)，可核符号及二分之一。文献Lorentz记号不直接复制为Euclidean同符号。式(6)没有确定有限反项、非局部作用或某物理态的真空能。

## 5. 固定变量离壳算符类的不闭合

在平直周期盒中取目标一小段单位速测地线γ，以q(x)=a sin(nx¹)参数化φ=γ(q)，固定a和整数n，横向体积为1、x¹周期2π。于是Q=S²：

$$
\int\left(\frac{S^2}{27}+\frac Q{54}\right)
=\frac1{18}\int_0^{2\pi}a^4n^4\cos^4(nx^1)\,dx^1
=\frac{\pi a^4 n^4}{24}>0.
\tag{7}
$$

任意固定的δU(φ)沿该族积分与n无关；任意δK_ab(φ)的二导数项为n²乘同一常数。纯曲率及场依赖f(φ)R在平直盒为零，全导数的周期积分为零。因此

$$
C_0+C_2n^2=\frac{\pi a^4}{24}n^4\quad\text{不可能对全部正整数 }n\text{成立}.
\tag{8}
$$

这已经是不能只改U、K及无标量梯度几何项的离壳反例，亦排除了把四梯度项称为全导数。有限数值只检查代表函数；全域结论由频率幂次证明。此为测试函数族，不预言任意高频物理仍在EFT适用范围。

允许含导数的场重定义、最低阶EOM或改变度规变量时，算符可在基底间移动；本轮不声称S²、Q、RS各自都是独立可观测耦合。完整其它环也可能改变净结果。后继应压缩这些冗余，而非把每个热核单项提升为新认知原则。

## 6. 同一有效项的几何响应

定义B_mu_nu=K_ab ∂muφa∂nuφb，S=g^mu_nu B_mu_nu，Q=B_mu_nu B^mu_nu。固定φ，取α=1/27、β=1/54。对I4=∫√g L4、L4=αS²+βQ：

$$
\mathcal V^{(4)}_{\mu\nu}:=\frac2{\sqrt g}\frac{\delta I_4}{\delta g^{\mu\nu}}
=-g_{\mu\nu}L_4+4\alpha S B_{\mu\nu}+4\beta B_{\mu\rho}B_\nu{}^\rho,
\qquad g^{\mu\nu}\mathcal V^{(4)}_{\mu\nu}=0\quad(d=4).
\tag{9}
$$

迹为零是该局部四导数项的四维尺度性质，不等于完整量子理论无迹反常。对IRS=∫√g RS，R变分与S变分必须同时保留：

$$
\mathcal W_{\mu\nu}:=\frac2{\sqrt g}\frac{\delta I_{RS}}{\delta g^{\mu\nu}}
=2\left[S G_{\mu\nu}+(g_{\mu\nu}\Box-\nabla_\mu\nabla_\nu)S+R B_{\mu\nu}\right],
\qquad g^{\mu\nu}\mathcal W_{\mu\nu}=6\Box S\quad(d=4).
\tag{10}
$$

若漏R B项，迹错误地多出−2RS。所有公式均为Euclidean正二倍逆度规变分；Lorentz物理T常用负二倍变分，还涉及Wick及作用符号，不能直接以本表称物理压力。

平坦盒上的同一测地线背景给S=a²n²cos²(nx¹)，故

$$
R=0,\quad I_{RS}=0,\qquad
\mathcal W_{11}=0,\quad
\mathcal W_{jj}=-4a^2n^4\cos(2nx^1)\ (j=2,3,4).
\tag{11}
$$

因此b4中的−RS/9仍有非零几何响应。它不是能从平坦背景作用值判断可删的项。另一独立核验采用g_mu_nu=e^{2σ(x¹)}δ_mu_nu，仍令B_11=(q′)²，其余零：

$$
\sqrt g RS=-6(q')^2[\sigma''+(\sigma')^2],\qquad
\delta_{\sigma=f}I_{RS}=-6\int(q')^2(f''+2\sigma'f')
=-\int\sqrt g\,f\,g^{\mu\nu}\mathcal W_{\mu\nu}.
\tag{12}
$$

式(9)—(12)接通指定局部有效项与几何变分；只有把它们连同共同物质变化放回同一有效作用，才能继续建立同一来源。尚未得到完整应力的态期望、守恒的全部量子来源或自洽Einstein解。

## 7. 553与574的量子变量不能直接拼接

经典关系gE=F(φ)gJ依赖所积分的标量。于是

$$
\delta g_E=F\delta g_J+\delta F\,g_J,\qquad
\delta g_E=0\ \Longrightarrow\ \delta g_J=-\frac{\delta F}{F}g_J,
\qquad
\delta g_J=0\ \Longrightarrow\ \delta g_E=\delta F\,g_J.
\tag{13}
$$

两个“只积分标量、冻结各自度规”的程序沿不同的配置空间方向积分。不能把553的部分β函数直接视为当前目标的同一量子处方。δF一阶偶为零亦不自动解决二阶变分、背景分裂、测度及源的匹配。

此为部分量子化的区别，**不是**完整Jordan／Einstein量子理论不等价。要合并需声明全部积分变量、规范固定、Jacobian、反项及同尺度匹配。574有限图Hamiltonian到本四维协变路径积分的映射亦未证明。

## 8. 可复算证据与核验范围

[代码](../joint_scalar_effective_geometry.py)、[结果](../joint_scalar_effective_geometry_results.json)默认复算并逐项比较完整JSON，不覆盖旧证据。四组已保存入口首次晋升正式检查，加三组独立验证，共七组。

|检查|保存结果|
|---|---|
|同一原势|五分量真空Hessian三角向值近零，非零值约.04683851145、.13805672781，连接原两径向谱；不把角向零值认定物理粒子|
|直接曲率收缩|矩阵法与不变量展开差5.55×10⁻¹⁷，正交标架变换保持|
|梯度／曲率缩放|固定种子四梯度系数.11054149035437；R—梯度交叉−.16876816042289=−S/9|
|冻结变量映射|δF=−.00448799944466，两个背景冻结方向不同，有限差分误差2.77×10⁻¹²|
|周期非吸收|a=.24时n=1、2、3、5，I4约.000434294、.00694870、.0351778、.271434；用n=1、2拟合C0+C2n²后n=3残差.01737175|
|四梯度度规变分|三个一般正定度规与非对角变化，有限差分误差≤9.74×10⁻¹¹，四维迹为数值零|
|曲率梯度响应|曲背景直接／张量迹变分11.4002114213466，有限差分差4.09×10⁻¹³；漏RB时为11.3923127034333。平坦背景IRS=0但所选共形变化响应54.2867210540|

最后一组初版区分错误公式的断言要求差值>.01，实际差约.00789872；解析式及两种正确算法一致，因此按远大于10⁻¹⁰数值误差的>.001重新核验。没有调整方程或保存失败结果。

## 9. 条件压缩与下一步

[580条件账](../unified_physics_condition_ledger_580.md)增加的是共同量子处方的必要说明，未新增基本认知公理。当前已排除“只保原二导数标量项即可吸收这一标量行列式”的固定变量离壳接法；保留完整EFT、有限截止及其它量子处方。

下一步优先检验：在明确允许动态几何及受控场重定义后，哪些局部项只是同一物理结构的不同表示，哪些确需新的独立匹配数据。至少要把物质、度规及记录的变换共同追踪，不能删项后仍沿用原观察量。580出现的更多算符不自动等于同样多的新自由参数。

统一目标、维数与手征物质来源、观测相容性仍开放；本轮不是连续量子引力或全物理结项。
