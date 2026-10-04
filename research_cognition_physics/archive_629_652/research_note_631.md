# 第631轮：原整代物质的张量应力谱与共同曲率系数

日期：2026-10-01。接[630](research_note_630.md)，复用[553](../archive_531_553/research_note_553.md)、[599](../archive_585_628/research_note_599.md)、[601](../archive_585_628/research_note_601.md)和[620](../archive_585_628/research_note_620.md)。[代码](631/joint_tensor_stress_spectrum.py)、[结果](631/joint_tensor_stress_spectrum_results.json)、[核验](631/research_round_631_checks.json)、[条件账](631/unified_physics_condition_ledger_631.md)。三组检查、十四式；主代理审查，无新增独立代理审查。

## 1. 从一个几何方向扩展到同一物质的完整二点切割

630把原整代质量、实时谱与共形来源相连，但共形变化只读取应力迹。本轮保持同一原质量、真空和连续动能，检验其余度规来源。**结果：同一原物质给一个迹通道和五个横向无迹通道；后者的正谱复现553已有费米Weyl平方运行系数。质量、量子噪声和曲率对数因此再共享一份数据。**

本轮新增原物质的明确张量谱、质量来源的联合Gram矩阵、旧热核的正确归一及色散所需额外局部资料。一般张量分解和正谱原理是成熟成果，不计为新发现。

|层次|内容|
|---|---|
|认知动机|同一几何的不同变化方向与同一物质／态相容|
|继承输入|630原一代16-Weyl质量、半权计数、3+1平直背景、canonical费米子、真空|
|新增检查对象|Hilbert／Belinfante对称应力的全部非接触二点切割；一般类时外动量的Lorentz重建|
|解析结果|同一原谱的迹／自旋2分解、联合正性及Weyl平方对数匹配|
|数值检验|原所有质量的72方向Dirac矩阵积分、协变投影、旧热核有理系数、色散余项|
|未完成|一般曲背景的全量子Ward身份、动态引力、实际真空自洽性、图到连续映射、区域记录|

“自旋2”在这里标识**类时外动量下应力响应的张量表示**，不是发现新引力子。五个分量不是无质量引力子的五种在壳偏振。固定平直背景也尚未被证明满足原完整量子几何方程。

## 2. 成熟框架与项目对象

[Karateev，2012.08538，§3及§5](https://arxiv.org/html/2012.08538)给应力二点函数的迹／自旋2谱分解、正性和UV归一的关系；其定义会把动量幂和2π放在不同位置，不能直接照搬常数。下面从630相同的canonical矩阵和相空间测度固定归一。

[Hu与Verdaguer，0802.0658，§3](https://arxiv.org/html/0802.0658)说明应力噪声是正的双点分布，局部c-number重整化不改变其连通涨落；其随机引力方程预先采用半经典几何方程。这里只对接噪声／响应接口，没有用该方程证明Einstein作用来源，也不把二点匹配升级为完整量子过程。

553已存在Weyl平方项、费米粒子计数和运行结果。本轮复算相同费米部门与真实切割之间的映射，不再宣布首次发现高阶曲率项。630的原质量、参考态和冻结证据保持不变。

## 3. 从canonical应力得到两个谱

令ρ_m为630的质量双线性谱。在总空间动量零、ω>0的真空成对通道，守恒使应力的时间分量无切割；空间应力的迹与质量来源满足：

$$
\rho_{0\mu,\alpha\beta}=0,\qquad
T_{\mathrm{tr}}=\frac1{\sqrt3}\sum_iT_{ii},\qquad
(T_{\mathrm{tr}})_{-+}=-\frac1{\sqrt3}O_{-+}.
\tag{1}
$$

(1)是非接触跃迁关系；真空一次期望、时间排序接触项及局部反常不由它删除。符号与630的(−+++)和正Hamiltonian质量来源一致。

给任意实对称空间矩阵e，记 \(T_e=e_{ij}T_{ij}\)。canonical一粒子顶点为

$$
V_e=(e\mathbf k)\cdot\boldsymbol\alpha,\qquad
H_{\mathbf k}=\boldsymbol\alpha\cdot\mathbf k+\beta m,\quad
P_\pm=\frac12(1\pm H_{\mathbf k}/E),\quad
V_O=m\beta .
\tag{2}
$$

Clifford迹和球面二阶／四阶矩给

$$
\operatorname{tr}(P_-V_eP_+V_e)
 =2\left(|e\mathbf k|^2-\frac{(\mathbf k\cdot e\mathbf k)^2}{E^2}\right),\quad
\left\langle\operatorname{tr}(P_-V_eP_+V_e)\right\rangle_\Omega
 =\frac{2k^2}{3}\operatorname{tr}e^2
 -\frac{2k^4}{15E^2}\left[(\operatorname{tr}e)^2+2\operatorname{tr}e^2\right].
\tag{3}
$$

当tr e=0、tr e²=1时，乘630相同的壳测度kE/(2π)，再取原d_a=1/2，得到

$$
\rho_2(\omega)=\sum_a\frac{d_a\omega^4}{240\pi}
 \left(1-\frac{4m_a^2}{\omega^2}\right)^{3/2}
 \left(3+\frac{8m_a^2}{\omega^2}\right)\theta(\omega-2m_a),
\qquad \rho_0(\omega)=\frac{\rho_m(\omega)}3 .
\tag{4}
$$

两者非负且来自同一批物质。对于本真空自由场，谱支持在未来／过去类时成对区域；不存在额外单粒子引力极点。其余动量可由已输入的Lorentz协变重建。设 \(s=-p^2>0\)：

$$
\pi_{\mu\nu}=\eta_{\mu\nu}+\frac{p_\mu p_\nu}{s},\qquad
P^{(0)}_{\mu\nu,\rho\sigma}=\frac13\pi_{\mu\nu}\pi_{\rho\sigma},\qquad
P^{(2)}_{\mu\nu,\rho\sigma}
 =\frac12(\pi_{\mu\rho}\pi_{\nu\sigma}+\pi_{\mu\sigma}\pi_{\nu\rho})
 -P^{(0)}_{\mu\nu,\rho\sigma}.
\tag{5}
$$

对正能类时动量，

$$
\rho_{\mu\nu,\rho\sigma}(p)
 =\rho_2(\sqrt s)P^{(2)}_{\mu\nu,\rho\sigma}
 +\frac{\rho_m(\sqrt s)}3P^{(0)}_{\mu\nu,\rho\sigma},\qquad
p^\mu\rho_{\mu\nu,\rho\sigma}=0 .
\tag{6}
$$

这是所声明自由费米部门的完整非接触应力二点切割；不是全模型含标量、规范、引力圈的完整应力，更不是非线性几何方程的证明。真实完整二次度规变分还含来源自身的二阶变分及接触项。

## 4. 与630共同来源合并

选择五个正交归一的空间无迹e_A。在来源约定与630一致时，径向q、共形σ和五个归一剪切来源的切割矩阵为

$$
\rho_{\mathrm{joint}}=
\begin{pmatrix}
\rho_m vv^\mathsf T&0\\
0&\rho_2 I_5
\end{pmatrix},\qquad
v=(b,1),\quad b=3.71869156775,\qquad
N_{\mathrm{joint}}=\frac12\rho_{\mathrm{joint}}\quad(\omega>0).
\tag{7}
$$

迹与无迹部门的交叉为零来自旋转对称真空和角积分；不能推广到任意状态、各向异性背景或非线性来源。q与σ之间的非零交叉必须保留。两谱严格正时(7)秩6；只能用一个标量随机来源复现全部二点涨落的指定接法失败。

在630同一频率ω=1.19727400923：

$$
\rho_m=0.0262795552418,\qquad
\rho_2=0.0583060690399,\qquad
N_{\rm shear}=0.0291530345200 .
\tag{8}
$$

从原16个质量逐项直接做7×7 Gram矩阵的72方向角积分，最大解析／矩阵差2.78×10⁻¹⁶；20个类时动量的协变投影收缩误差1.37×10⁻¹³。Gram矩阵同时检查O、归一应力迹与五个剪切方向的相关性，未通过把五个谱手工拷贝来替代计算。

作为边界检验，保持同样费米分量数而令质量为零：

$$
\rho_m=\rho_0=0,\qquad
\rho_2(\omega)=\frac{D}{80\pi}\omega^4,\quad D=\sum_ad_a=8 .
\tag{9}
$$

所以没有质量迹涨落不代表没有几何涨落。此处O是**相对质量缩放来源**；原绝对标量q来源在零背景处须重新求导，未被宣布为零。局部迹反常也是独立的c-number／接触资料，(9)没有否认它。

## 5. 同一正谱复现旧Weyl平方系数

高频时(4)给

$$
\rho_2(\omega)=\frac{D}{80\pi}\omega^4
 -\frac{S_2}{24\pi}\omega^2-\frac{S_4}{8\pi}+O(\omega^{-2}).
\tag{10}
$$

因此该应力迟致核的z⁴对数发散系数是D/(80π²)。独立从599采用的Dirac平方算符、spin联络曲率迹及[Vassilevich热核式4.28](https://arxiv.org/pdf/hep-th/0306138)得到每个Dirac的纯几何系数：

$$
b_{4,D}^{\rm geom}
 =\frac{-7\,\mathrm{Riem}^2-8\,\mathrm{Ric}^2+5R^2}{360}
 =-\frac1{20}W^2+\frac{11}{360}E_4 ,
\qquad \text{忽略已声明的体全导数}.
\tag{11}
$$

该成熟系数在553已使用。本轮用有理数复算其基底变换，并核对应力谱的归一。取 \(h_{ij}=2\epsilon e_{ij}\)、tr e=0、tr e²=1、只依赖时间；在Euclidean频率Q下，Weyl平方作用的二阶展开与Hessian为：

$$
\left.\int W^2\right|_{\epsilon^2}=2\int(\epsilon'')^2,\qquad
\frac{\delta^2}{\delta\epsilon^2}
 \left(b_W\int W^2\right)=4b_WQ^4,\qquad
4\beta_{b_W}^{\rm f}=\frac{D}{80\pi^2},
\quad \beta_{b_W}^{\rm f}=\frac{D/20}{16\pi^2}.
\tag{12}
$$

最后一式与(10)完全相符。符号须保：费米裸Euclidean圈的W²对数系数为−(D/20)/(16π²)，反项及553的β约定为正；\(\Gamma_E^{(2)}=-\Pi(iQ)\) 与630一致。不能以正谱直接断言任意局部四阶引力作用无ghost或全尺度稳定。

本轮是一代16-Weyl，c_f=D/20=0.4、β_bW^f=0.00253302959106；不是553三代加五实标量与十二矢量的完整c=293/120。Λ=500m_max时，谱端点给z⁴对数系数0.0101320565298，解析／热核值0.0101321183642，差异为预期质量幂次尾。

质量较低次的几何来源、真空能源及Einstein项还需要接触项和一次期望共同匹配。本轮只用最高四导数项作无歧义连接，没有仅凭切割就选出Newton常数、Λ或所有有限反项。

## 6. 几何核增加哪些共同条件

因为ρ₂增长为ω⁴，剪切核需要三次减除；630质量核的两次减除不足以原样覆盖它：

$$
\Pi_{2,\mathrm{ren}}(z)=a_0+a_2z^2+a_4z^4+\mathcal R_2(z),\qquad
\mathcal R_2(z)=\frac{z^6}{\pi}\int_0^\infty
 \frac{\rho_2(\omega)}{\omega^5(\omega^2-z^2)}\,d\omega .
\tag{13}
$$

a₀、a₂、a₄在这个单独通道中是有限匹配数据；合入完整协变作用后不能对每个度规分量独立选择，须同真空、Einstein、曲率反项及相应接触项共用。单通道三常数不意味着全理论有三条新增认知公理。

原正质量阈值ω₀不变，且

$$
K_6=\frac1\pi\int_0^\infty\frac{\rho_2(\omega)}{\omega^7}\,d\omega
 =\frac1{1344\pi^2}\sum_a\frac{d_a}{m_a^2}
 =0.0539861936287,\qquad
0\le\mathcal R_2(z)\le
 \frac{z^6K_6}{1-z^2/\omega_0^2}\quad(|z|<\omega_0).
\tag{14}
$$

这来自与630相同的x=2m/ω代换；相应积分为5/7。192／384点复算的最大变化6.72×10⁻¹⁹。此界仅约束所定义的非局部色散余项，不能控制任意有限局部项、全部曲背景或严格有限时间来源的全部频谱尾。

可以用一个标量和五个无迹经典Gaussian变量匹配这些二阶协方差，这是随机引力接口的局部二点表述；自由费米应力本身不因此成为Gaussian量子过程，其高阶关联仍在。引入随机变量不会替代原量子物质和实际记录。

## 7. 成果、剩余输入及下一步

本轮把630的迹方向扩为原同一物质的全部平直非接触应力二点结构，并将新的张量通道同553／599曲率系数联立。明确排除了“共形／单标量响应足以替代全部几何涨落”以及“质量核两次减除自动覆盖完整应力”的指定接法。

仍保留：3+1、给定真空与连续动能、有限反项、实际相互作用、背景自洽性、原图到连续理论映射。反常关闭和切割守恒都不自动构造全部接触Ward身份；这里也没有解Einstein方程或推导引力子。

[632入口](632/drafts/STATUS.md)应优先核同一态、背景方程和来源接触项的共同条件：630／631诊断背景尚非原全理论的自洽解，不能径直作为随机引力的背景。先回查旧真空与同阶变分工作，再把新增实时／张量资料接进去，避免重复固定点或单系数扫描。继续整合条件，认知设计后置。

