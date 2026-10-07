# 第840轮：原区域参考有限性的模谱充要条件

日期：2026-10-05。接[839](research_note_839.md)和[冻结入口](840/drafts/STATUS.md)。

配套：[代码](840/regional_modular_domain.py) · [结果](840/regional_modular_domain_results.json) · [核验](840/research_round_840_checks.json) · [复算](840/verify_round840.py)。

## 1. 主结论与适用范围

**在829—839指定的组织准备下，其相对原自由区域参考的相对熵有限，当且仅当有限编码模式属于原区域一粒子模生成元的绝对一阶矩域。** 这将839整族有限性问题变成同一原区域协方差的一份谱判据；没有证明此前任意选取的光滑紧支模式都满足它。

|层次|本轮地位|
|---|---|
|认知动机|同一组织的记录内容与工作参考要能共同比较，而非逐内容另换参考|
|候选要求与替代|若目标应用需要有限参考相对熵，必须核本轮域条件；也可研究其它准备/参考，不强加为认知公理|
|继承输入|730原Gaussian自由参考、829有限编码、838二点重置、839同一区域正常态及极限|
|解析增量|有限模能量的双边界与条件熵界，导出原区域有限性的充要判据|
|数值验证|固定同一协方差的嵌套压缩、全部交叉项、模谱端点；非原PDE计算|
|开放|原既选模式的实际域条件、动力能源连接、应力识别、相互作用及面积响应|

相对熵的单位是nat；模谱是无量纲参考谱，不能直接称为真实物理能源或准备功。这里仍是原自由费米区域，不包括已完成的全引力区域量子化。

## 2. 原区域对象及需检查的矩

沿839，H_B为同一区域的Cauchy空间，E_B是其嵌入全空间的投影。P为730原纯协方差，但其区域压缩一般不是投影：
$$
A=E_BPE_B|_{H_B},\qquad
h=\log\frac{I-A}{A},\qquad
L=-\log A-\log(I-A),\qquad CAC=I-A .
\tag{1}
$$

端点0、1通过非负谱形式处理；若编码模式含端点谱质量，下面的绝对矩为∞，不把零本征值加底噪。Q是20维C不变编码投影，q=20、d=2^{q/2}=1024。取任意C实正交基u₁,…,u_q：
$$
Z_B=\operatorname{Tr}_Q|h|
:=\sum_{j=1}^{q}\langle u_j,|h|u_j\rangle,\qquad
Z_B<\infty\ \Longleftrightarrow\
\operatorname{ran}Q\subset\operatorname{Dom}|h|^{1/2}.
\tag{2}
$$

这是绝对一阶矩域，不要求每个u都在Dom h中；Q上的迹与基无关。原空间与真实动力H没有在式(1)中被重新定义。

## 3. 同一有限限制的模能量和条件熵

沿839真正嵌套的H_n⊃ran Q，记A_n为A的压缩，R_n=I_n−Q。Gaussian化态g_n的协方差是
$$
\widetilde A_n=Q/2+R_nA_nR_n,\qquad
h_n=\log\frac{I_n-A_n}{A_n},\qquad
\Delta K_n=\frac12\operatorname{Tr}_{H_n}
 h_n(\widetilde A_n-A_n).
\tag{3}
$$

1/2来自自对偶计数，保配对和跨块项。γ_n是原参考限制，ν_n是其余部限制。先在支持相容的有限层计算；非相容时b_n=∞。重置熵增满足
$$
s_n=S(g_n)-S(\gamma_n)
=\log d+S(\nu_n)-S(\gamma_n),\qquad
0\le s_n\le2\log d,\qquad
b_n:=D(g_n\Vert\gamma_n)=\Delta K_n-s_n .
\tag{4}
$$

上下界分别来自次可加性与Araki—Lieb不等式；并非从无限区域普通熵相减得到。有限层增加余部时，[Lieb—Ruskai的强次可加性](https://www.numdam.org/item/RCP25_1973__19__A5_0/)使s_n单调增加；839使b_n单调增加。因此ΔK_n也增加，但不据此断言|h_n|的压缩单调。这些成熟熵不等式不计为本轮发现。

## 4. 有限模能量的双边估计

令A_Q=QA_nQ、H_Q=Qh_nQ、Z_n=Tr_Q|h_n|。自对偶性给Tr_Q h_n=0。将838重置差代入式(3)：
$$
\Delta K_n
=\frac12\left[
 \operatorname{Tr}_Q(A_QH_Q)
 -2\operatorname{Tr}_Q(A_nh_n)\right].
\tag{5}
$$

所有乘积均由原完整A_n生成，未将跨块关联删除后再计算。由于∥A_Q−I_Q/2∥≤1/2、∥H_Q∥₁≤Z_n，且A_n、h_n对易，
$$
0\le\Delta K_n
\le\frac12\left(\frac12Z_n+Z_n\right)
=\frac34Z_n .
\tag{6}
$$

此处Tr_Q(A_nh_n)=Tr_Q((A_n−I/2)h_n)；后者非正，因为h为A的递减log-odds。标量t=|h_n|给|A_n−1/2|=(1/2)tanh(t/2)。再用t/(1+e^t)≤1/e：
$$
\begin{aligned}
\Delta K_n
&\ge\operatorname{Tr}_Q|A_n-I/2|\,|h_n|-\frac14 Z_n\\
&=\frac14 Z_n-\operatorname{Tr}_Q
 \frac{|h_n|}{1+e^{|h_n|}}
\ge\frac14 Z_n-\frac{q}{e}.
\end{aligned}
\tag{7}
$$

结合式(4)：
$$
\frac14Z_n-\frac qe-2\log d
\le b_n\le\frac34Z_n .
\tag{8}
$$

常数只需一致，不宣称最优。纯模式若与Q有重叠，重置使参考支持不再包含后态，b_n=∞；若与Q正交，纯模式在两态中共同固定，可删除后使用同一证明。这和839支持处理一致。

## 5. 从压缩到同一原区域谱

对−log应用[Hansen—Pedersen算符Jensen不等式](https://arxiv.org/pdf/math/0204049)的等距压缩形式，得到二次型意义的
$$
|h_n|\le L_n:=-\log A_n-\log(I_n-A_n)
\le E_n L E_n|_{H_n},\qquad
|h|\le L\le |h|+2\log2\,I .
\tag{9}
$$

可先用−log(x+δ)证明有界正则化版本，再令δ↓0取非负形式极限；这只是证明手段，不更改物理参考。于是Z_B有限时，所有Z_n≤Z_B+2qlog2，式(8)给统一b_n上界。

反向不假定Z_n收敛。把A_n在H_n外扩成半单位，其强极限为A。对连续有界函数f_M(x)=min{M,|log((1−x)/x)|}（端点定义为M）：
$$
\operatorname{Tr}_Qf_M(A_n)\longrightarrow
\operatorname{Tr}_Qf_M(A),\qquad
Z_B\le\liminf_{n\to\infty} Z_n .
\tag{10}
$$

第二式先固定M取极限，再令M增加。若原谱在端点有质量，结论仍按∞成立。故不能靠不同截断处的抵消，把无穷原模矩藏起来。

## 6. 原区域充要条件和误差范围

由839的b_n↗b_B及式(8)—(10)，得到
$$
\boxed{\ b_B<\infty\ \Longleftrightarrow\ Z_B<\infty\ },\qquad
\max\left\{0,\frac{Z_B}{4}-\frac qe-2\log d\right\}
\le b_B
\le\frac34\bigl(Z_B+2q\log2\bigr).
\tag{11}
$$

式(11)在∞情形亦表达相应发散。再由839，整族ω_{ε,r}相对同一原参考有限，当且仅当这一个域条件成立；有限记录项c直接相加。无需逐内容增加新参考条件。

这里的区域极限与正常态前提已由839完成；新要检查的是**既选紧支模式在具体A中的谱权重**。本轮没有用48维诊断替原无限背景验收这项。若条件成立，839的背景相对熵差才可安全相减；若不成立，内容彼此仍可比较，完整目标也未被反证。

## 7. 动力能源不能未经证明替代模域

[637—645](../archive_629_652/README.md)在指定Gibbs参考下能用真实H控制相对熵；本轮的局部γ_n一般不是那个Gibbs态。813/829的光滑资料和动力能源矩，不能仅凭符号相似代替式(2)。

一个抽象CAR例子说明缺少谱比较的危险。取一粒子“能源”E_j=j，模式系数a_j=√(e²−1)e^{-j}，参考占据p_j=(1+exp(e^{3j}))^{-1}，按自对偶伙伴补全。它是合法Gaussian协方差，且可作纯Gaussian扩展，但未被证明是原时空Hadamard参考：
$$
\sum_{j\ge1}|a_j|^2j^{2m}<\infty
\quad\hbox{每个有限 }m,\qquad
2\sum_{j\ge1}|a_j|^2
 \left|\log\frac{1-p_j}{p_j}\right|
=2(e^2-1)\sum_{j\ge1}e^j=\infty .
\tag{12}
$$

所以“所有多项式能源矩有限”作为抽象信息仍不足。这个例子不是原730背景的反例，也未排除原Dirac局域结构提供更强模谱估计。真正下一步是证明那项连接，而不是继续增加任意多项式矩。

## 8. 验证与下一项

一组联合检查，累计3624：

- 同一48维自对偶协方差，固定20维编码Q，按20、24、28、36、48维嵌套压缩，逐级保原配对和交叉项。b、s、ΔK均单调，式(3)—(9)通过。
- 双边界及独立Gaussian相对熵公式残差小于4.7×10⁻¹⁵；代码中“full reference”仅指这一固定有限诊断，不是原730 PDE解。
- 一模式log-odds升到120时，b≈59.30685，核模谱大值的真实增长；抽象序列部分和用于校准式(12)，发散由解析几何级数证明。
- 未计算原A_B本征谱、未签收所选模式的Z_B有限性；没有面积、Newton系数或自治准备新结论。

接[841入口](841/drafts/STATUS.md)：回到原模式的共同实现，核局部性、物理正则域和模域能否由同一已用背景/参考条件保证。若需对接局部模流或成熟有限相对熵定理，逐项映射其前提，避免把平直真空或任意模平滑模式替换成原对象。
