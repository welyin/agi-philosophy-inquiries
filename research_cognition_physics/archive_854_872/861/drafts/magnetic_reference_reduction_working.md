# 861工作稿：磁标签的配置图册、规范固定与完整动量运输

日期：2026-10-06。正式860／累计3645。接[入口](STATUS.md)、[860](../../research_note_860.md)及[858](../../research_note_858.md)。[精确校准](../magnetic_coordinate_canonical_probe.py) · [结果](../magnetic_coordinate_canonical_probe_results.json)。本稿还未完成材料胞元与约化Hamiltonian/尺度来源的全部连接，不发表861。

## 1. 范围和旧成果区别

先研究860的两导数S_base在经典零费米部门的ADM相空间，保既有p和原非阿贝尔颜色材料。860的导数通信作用仍按局部EFT阶处理；不能从其形式量子存在性直接推断已给无穷阶精确Hamiltonian。

858给的是在590独立控制模型中**额外施加物理约束**，且原材料参考依赖旧速度。当前改查实际Einstein经典约束下的**坐标规范固定**：h作为时间，Y=(s,p,M_h)作为空间标签。规范固定是同一物理解的呈现，不等于另加w=V_R、Π=0。检索磁标量/配置型参考主题只命中859—860，未发现已有完整此映射。

## 2. h切片上有真正配置型的磁标签

在dh类时、n·dh>0的片中，以h=τ选时间，则D_i h=0，M_h限制为

$$
M_\Sigma=\frac12\gamma^{ik}\gamma^{jl}
\langle F^c_{ij},F^c_{kl}\rangle .
$$

它只含空间γ和A及空间导数。可先在任意ADM相片定义M_Σ这份配置函数，再以h=τ限制；离开该规范面时不冒称M_Σ等于协变M_h。

χ=(h−τ,s−y¹,p−y²,M_Σ−y³)全部是配置型函数，彼此Poisson括号为零。内部Gauss与它们对易。859四参考满秩意味着在h等值面上J^a_i=D_iY^a可逆；n·dh给正的实际v_h。

## 3. 规范矩阵是微分算符，不能漏掉lapse梯度

设C_perp[N]+C_i[β^i]为原Einstein约束生成元；空间生成元取能对内部规范不变量实现Lie导数的版本。规范面上

$$
\delta\chi^0=v_hN,\qquad
\delta\chi^a=\mathcal B^a[N]+J^a_i\beta^i,
\qquad
\mathsf M=\begin{pmatrix}v_h&0\\\mathcal B&J\end{pmatrix}.
$$

不能将所有B都当作乘法速度。磁标签含空间连接导数，若v_Aj=δC_perp[1]/δΠ_Aj，则其正常变化有

$$
\mathcal B^M[N]=N\mathcal B^M[1]
+2\langle F^{ij},v_{Aj}\rangle\partial_iN .
$$

度规动量的正常变化另计在B^M[1]；内部纯规范项因M规范不变而消失。正是这份梯度项阻止“只代入单位lapse的旧速度”。

在同一局部正则片，矩阵虽非纯代数，其逆仍明确：给定f=(f⁰,f^a)，

$$
N=\frac{f^0}{v_h},\qquad
\beta^i=(J^{-1})^i_a
\left(f^a-\mathcal B^a\!\left[\frac{f^0}{v_h}\right]\right).
$$

对紧支测试仍保支集，不需要额外椭圆求逆。保持h=τ、Y=y时，取f⁰=1、f^a=0，从实际相点共同算lapse与shift。只给局部规范保持，未给全局无Gribov图册。

原约束和Gauss都成立时，其第一类括号弱为零；χ与C的完整第二类矩阵有上述可逆混合块。对任意两份配置泛函F(q)、G(q)，其Dirac括号为零，因为χχ括号零且F、G与χ对易。这是**在该经典局部规范面上的配置读数**陈述，不代表体积与全部物质动量对易，也不代表原量子图代数已经映射成功。

## 4. 原磁参考还允许显式的配置/动量图册

这一步补上“配置形式”不等于“完整canonical运输”的风险。选固定局部坐标体积，分解

$$
\gamma_{ij}=e^{2r/3}\bar\gamma_{ij},\quad\det\bar\gamma=1,\quad
\sqrt\gamma=e^r,\quad
M_\Sigma=e^{-4r/3}\bar M,\quad
\bar M=\tfrac12\bar\gamma^{ik}\bar\gamma^{jl}\langle F_{ij},F_{kl}\rangle .
$$

753/859色磁场给M>0；小片上可用M代替共形r，而保留原A和五份单位行列式度规配置。逆为r=(3/4)log(\bar M/M)。r与\barγ本身按固定坐标密度定义，不把\barγ冒称普通空间张量。

写旧辛势中的π_r δr、P_bδb，b包含\barγ和A。对紧支变化，完整拉回为

$$
\pi_M=-\frac{3\pi_r}{4M},\qquad
P_b^{\rm new}=P_b+\frac34(D_b\bar M)^*
\left(\frac{\pi_r}{\bar M}\right),\qquad
\Theta_{\rm old}=\int\left(\pi_M\delta M+P_b^{\rm new}\delta b+\cdots\right).
$$

此处D_b是Fréchet导数，星号含空间分部积分。单位行列式度规的变化限制在无迹切空间。对颜色连接，记\bar F^{ij}=\barγ^{ik}\barγ^{jl}F_kl，有

$$
\Pi_A^{{\rm new},j}=
\Pi_A^j-\frac32D_i\left(\frac{\pi_r}{\bar M}\bar F^{ij}\right)
=\Pi_A^j+2D_i\left(\frac{M\pi_M}{\bar M}\bar F^{ij}\right).
$$

因此把色磁标量当配置坐标会产生确定的电动量平移。原电场、完整Hamiltonian、约束、记录和来源都必须作逆变换；不能保新M名称却把新Π_A当原电场。这个变换在声明非零磁片中可逆，纯配置逆含空间导数，但无新时间导数。

连续内部Gauss的保持可由原/新辛势及生成函数的规范不变性直接核对；不能拿有限差分版本当精确非阿贝尔格点规范理论。

## 5. 体积接口已简化，但原图匹配未完成

在h切片上，856的材料密度形式现在给

$$
n=\eta\frac{|\det D_i(s,p,M_\Sigma)|}{\sqrt\gamma},\qquad
V_R=\frac{\sqrt\gamma}{\eta|\det D_i(s,p,M_\Sigma)|} .
$$

其当前表达只含配置。转到Y标签坐标后Jacobian为1，可写V_R=\sqrt{\gamma_Y}/η。这里η和胞元离散方式仍属声明的材料测度/调节字典，不能从坐标变换直接推导其数值或原578的量子谱下界。

有两条必须继续核清的边界：

1. 配置对易不等于与全部原物质代数对易。被动canonical变换保原全部Poisson关系，必须连电场的上述平移一起看；随后Einstein约束约化也不是原590独立控制模型的等价证明。
2. 要把V_R代入尺度减除，需要完整约化H及其几何/材料来源。当前仅得到参考图册及其辛运输，没有证明原量子图、热态和记录过程等价，也没有取消579的剩余发散。

## 6. 精确校准与下一项

代码在五点周期差分泛函上使用不对易的实SU(3)反Hermitian生成元，保非零磁范数，计算真实曲率变化及其转置。四份独立配置变化的旧/新辛势差严格为零；漏电动量平移则四份都非零。完整逆运输保电能，直接将新动量代作旧电场则电能不同。

这是canonical导数/分部积分校准；尚未数值验证无迹度规变化，未求Einstein约束，也不是精确格点规范模型。当前不将它单独计为861科学完成。

下一步把第3—4节共同代入原约束的局部约化，核物理Hamiltonian、完整来源及材料体积反项；涉及860高导数通信时明确EFT阶与约化顺序。若必须新增动力学条件，应列为新增输入，不能借图册改名跳过。
