# 第771轮：同一完整一圈来源的局部联合Ward修复

日期：2026-10-04。接[770](../../research_note_770.md)及[771工作报告](research_note_771_working.md)。[代码](../joint_local_ward_repair.py)、[结果](../joint_local_ward_repair_results.json)、[核验](../research_round_771_checks.json)、[条件账](../unified_physics_condition_ledger_771.md)。三组检查、十九式；主代理终审，无独立代理或图像检验。

## 1. 本轮实际关闭的接口

**在原753在壳光滑背景、767／768同一自由物理Hadamard态，以及下述明确的局部Wick规范化类别内，完整首阶tadpole来源可以共同满足无穷小微分同胚和内部规范Ward身份。** 735原费米来源直接加入；本轮补此前缺失的耦合玻色、引力和ghost部门。没有用某个参考真空把实际来源置零。

这是固定背景上一圈、一次来源插入的局部规范化存在性与构造。**不证明全部相互作用QME、规范固定独立性、全阶背景独立、有限量子强度自洽，或原Q过程到连续理论的等价。** 各有限物理重整化常数仍需另定；形式首阶响应不能当作完整量子引力。

上一目标轮完成770并保存771工作入口，属于有效进展。本次核导航、最新报告、结果和进程；未发现运行中的Python。回用734、735、765、768—770，不重证空间382—386、425、522—523；604和649／699范围保持。没有改目标、建任务或定时任务。

|层次|采用内容|
|---|---|
|认知动机|同一状态、记录与几何/物质反馈须共同满足交换关系|
|继承模型|原四维作用、群/表示/参数、F>0、短时全局双曲背景及紧初片|
|继承量子输入|原完整自由物理态、BRST扩展及735费米规范化|
|本轮处方选择|固定共同尺度；采用标准局部Hadamard有限jet；允许实、局部协变、状态无关的来源Wick接触重定义|
|解析增量|完整辅助场消元的局部接触、原两个耦合算符的余项jet、共同体积接触及联合Ward|
|数值角色|完整主部的接触号校准、非交换矩阵势的真实Hadamard输运jet、联合权重及分部积分检验|

允许类别只使用原算符、配对和有限阶背景jet的局部表达，光滑依赖F>0补丁；保原实性与尺度阶次，不引入状态相关系数、非局部逆、物理质量的倒数或新参考态。这是明确的局部复合来源处方，尚不主张它已来自全相互作用时间序乘积的共同规范化。

## 2. 将辅助块消元，但保全部接触

以下J_L、J_d及J_op专指玻色与ghost来源，原费米贡献在式(15)另加一次。量子阶数按ℏ的首阶系数记，不在每式重复写ℏ。

仍用765／768的V、W′配对及星号：

$$
 D_1=P+KK^\star,\quad D_0=K^\star K,\quad
 \mathbb L=\begin{pmatrix}P&K&0&0\\K^\star&-I&0&0\\
 0&0&0&D_0\\0&0&D_0&0\end{pmatrix},\qquad
 D_1K=KD_0.
 \tag{1}
$$

V含全部几何、规范和五标量混合，不能替成独立场之和。完成辅助平方，令β=b−K* v，旧字段由新字段经T给出：

$$
 T=\begin{pmatrix}I&0&0&0\\K^\star&I&0&0\\0&0&I&0\\0&0&0&I\end{pmatrix},
 \qquad T^\star\mathbb L T=\mathbb L_d
 =\begin{pmatrix}D_1&0&0&0\\0&-I&0&0\\0&0&0&D_0\\0&0&D_0&0\end{pmatrix}.
 \tag{2}
$$

T与T⁻¹都是局部有限阶微分变换，不取任何Green逆。由768实际态的强交织和两腿方程，变换后的实际核的v-β、β-v和β-β块全为零；v块仍是原λ₁，ghost仍是原λ₀。没有删除实际物理模式。

把770减除同样输送得H_t=T⁻¹𝕳T⁻*，并与对角处方H_d（v块H₁，β块0，ghost两个非零块H₀）比较。按算符指标记S_o=H_t−H_d，则仅辅助相关块非零：

$$
 (S_o)_{v\beta}=-q',\quad(S_o)_{\beta v}=-q,\quad
 (S_o)_{\beta\beta}=K^\star H_1K-D_0H_0-H_0D_0,
 \quad q=K^\star H_1-H_0K^\star,\quad q'=H_1K-KH_0.
 \tag{3}
$$

q、q′光滑，最后一块可写K*q′−H₀D₀，也光滑。它们由同一局部核构造，所需重合jet局部且状态无关。此处用局部Hadamard减除，不用非局部实际λ代替H。

T随背景变化，不能只保(2)而遗漏它的导数。转回密度Hessian，定义Ωη=T⁻¹δηT，得到

$$
 \widehat{\mathbb L}_d=T^{\rm st}\widehat{\mathbb L}T,\qquad
 T^{\rm st}\bar\delta_\eta\widehat{\mathbb L}T
 =\bar\delta_\eta\widehat{\mathbb L}_d
   -\Omega_\eta^{\rm st}\widehat{\mathbb L}_d
   -\widehat{\mathbb L}_d\Omega_\eta.
 \tag{4}
$$

st是密度/字段双指标的正式转置；它与(2)的纤维星号不可混用。令𝖧_t、𝖧_d、𝖲分别为H_t、H_d、S_o转成两个字段上指标后的核。由光滑C及两腿方程测试分部积分，原来源J_L与对角密度插入J_d满足

$$
 \begin{aligned}
 J_L(\eta)&=J_d(\eta)+\ell_T(\eta),\\
 \ell_T(\eta)&=-\tfrac12\operatorname{STr}_\Delta
       [\mathsf S\bar\delta_\eta\widehat{\mathbb L}_d]
 +\tfrac12\operatorname{STr}_\Delta
       [ (\widehat{\mathbb L}_d\mathsf H_t)\Omega_\eta^{\rm st}
            +(\mathsf H_t\widehat{\mathbb L}_d)\Omega_\eta].
 \end{aligned}
 \tag{5}
$$

所有括号内的方程余项均光滑；因此ℓ_T是可由原局部核计算的有限接触。普通规范参数与765协变参数之间的背景字典，必须随T、Ω、配对共同变化。没有声称“辅助场不传播，所以其来源接触自动为零”。

## 3. 保原配对，将来源变成两个正常双曲算符插入

采用标准实场半权和一对ghost的Grassmann权，令w₁=1/2、w₀=−1。M_i为两个束的原零阶密度配对；H_i、C_i按算符指标，C_i=λ_i−H_i。重复工作报告已证的密度乘积法则，

$$
 \theta_{i,\eta}=M_i^{-1}\bar\delta_\eta M_i,\qquad
 \ell_M(\eta)=\sum_{i=0,1}w_i\operatorname{Tr}_\Delta
       [(D_iH_i)\theta_{i,\eta}],\qquad
 J_{\rm op}=J_d+\ell_M.
 \tag{6}
$$

β的实际核和H_d块都为零，在J_d内无插入；其原接触已在(5)，未被遗忘。于是

$$
 J_{\rm op}(\eta)=\sum_{i=0,1}w_i\operatorname{Tr}_\Delta
      [C_i\bar\delta_\eta D_i],\qquad
 D_iC_i=-D_iH_i,\quad C_iD_i=-H_iD_i.
 \tag{7}
$$

这不是两个无混合真空算符：D₁、D₀依旧含原所有低阶混合。这里只把一次插入的余项计算化到正常双曲对象；原111辅助系统本身不被称为正常双曲算符。

## 4. 原矩阵算符的实际局部余项通式

每个D_i唯一写成□_{∇i}+E_i。由其原正式自伴性，∇i保持相应非退化纤维配对，E_i自伴；无需纤维配对正定。令σ为半测地间隔，ν=1/(8π²)，细化770未固定的局部有限部分为标准减除jet：

$$
 H_i=\nu\left(\frac{u_i}{\sigma_+}
       +\sum_{n\geq0}v_{n,i}\sigma^n\log(\sigma_+/\mu^2)\right),
 \qquad A_i=[v_{1,i}].
 \tag{8}
$$

这里只取足够阶的渐近展开；不假设无穷级数收敛。核的密度因子按M_i恢复。μ为固定共同长度尺度，实际态保持。其他局部有限部分的选择只改变明确的局部接触。

[Kamiński，1904.03708，Theorem 1及§I](https://arxiv.org/html/1904.03708)证明上述矩阵系数的双点伴随对称，允许非正定配对、任意签名及附加连接；原D₁／D₀满足这些前提。不是把标量Moretti定理乘字段数。

为使通式可计算，A_i等于−D_i局部热核二次系数的四分之一。以Ω_i为∇i曲率，有

$$
 A_i=\frac14\left\{\frac12E_i^2+\frac16RE_i
  +\frac16\Box_{\nabla_i}E_i+\frac1{12}\Omega_{i,\mu\nu}\Omega_i^{\mu\nu}
  +\frac{5R^2-2R_{\mu\nu}R^{\mu\nu}
     +2R_{\mu\nu\rho\sigma}R^{\mu\nu\rho\sigma}+12\Box R}{360}I_i\right\}.
 \tag{9}
$$

这里使用[Vassilevich，hep-th/0306138，§2.1及式(4.28)](https://arxiv.org/html/hep-th/0306138)的局部系数约定；其a₄是本文热参数二次项。只用局部恒等式，不引入Euclidean量子态。E_i和Ω_i由原完整算符确定，所有矩阵乘积保次序。

下面给所需的双点余项，而不仅给其迹。由(8)的输运递推，在对角线的一阶jet上D_iH_i等于ν[(□σ+2)v₁+2σ^{;μ}∇_μv₁]。利用[□σ]=4、[∇σ]=0和[∇′∇σ]=−g，得

$$
 [D_iH_i]=[H_iD_i]=6\nu A_i,\qquad
 [\nabla'(D_iH_i)]=\nu(6[\nabla'v_{1,i}]-2[\nabla v_{1,i}]),
 \qquad [\nabla'(H_iD_i)]=8\nu[\nabla'v_{1,i}].
 \tag{10}
$$

左、右递推由双点伴随对称配对。该计算也可对照[Zahn，1407.1994，命题4.3与4.8](https://arxiv.org/html/1407.1994)：这里只复用其正常双曲输运余项计算；上面已按本项目D=+□+E和σ约定重推，不调用该文的手征反常结论。结合Synge规则，Q_i=D_iH_i−H_iD_i满足

$$
 [Q_i]=0,\qquad [\nabla'Q_i]=-2\nu\nabla A_i.
 \tag{11}
$$

这是一份矩阵恒等式，连接曲率、非交换势与原低阶混合都留在A_i；并非用对角热核迹猜测所有双点项。Wightman处方没有Feynman逆核的delta；所列方程余项光滑，因果反对称部分不改变这些局部jet。

## 5. 用同一个局部体积接触关闭剩余Ward

对普通无穷小背景规范参数（含微分同胚ξ及内部参数），D_i按束的自然背景作用变换。其生成元形式为

$$
 U_i=\xi^\mu\nabla_{i,\mu}+u_i^{(0)},\qquad
 \bar\delta_R D_i=[U_i,D_i],\qquad
 b=\tfrac12\operatorname{tr}A_1-\operatorname{tr}A_0.
 \tag{12}
$$

u_i^(0)保留张量指标、内部规范、参数字典及参数导数的全部零阶项。由[Q_i]=0，它们在本次交换余项收缩中消失；这是计算结果，不是预先删掉内部规范或物质。用光滑C_i先分部积分，再用(11)，

$$
 J_{\rm op}(R)=-\sum_i w_i\operatorname{Tr}_\Delta(Q_iU_i)
       =-2\nu\int d{\rm vol}_g\,\xi^\mu\nabla_\mu b.
 \tag{13}
$$

纯内部规范变换ξ=0时该缺陷已经为零。对微分同胚，定义同一局部来源一形式

$$
 \ell_{\rm vol}(\eta)=-2\nu\int d{\rm vol}_g\,b\,
        \bar\delta_\eta\log\sqrt{|g|},\qquad
 \ell_{\rm vol}(R)=-2\nu\int d{\rm vol}_g\,b\nabla_\mu\xi^\mu
       =2\nu\int d{\rm vol}_g\,\xi^\mu\nabla_\mu b.
 \tag{14}
$$

δlog√|g|表示1/2 g^{μν}δηg_{μν}，本身是变分标量，未引入坐标依赖背景密度。紧支测试保证边界项为零。由(5)—(7)、(13)—(14)，得到所需完整修复：

$$
 J_{\rm ren}^{b+gh}=J_L+\ell_M-\ell_T+\ell_{\rm vol},\qquad
 J_{\rm ren}=J_{\rm ren}^{b+gh}+J_F^{735},\qquad
 J_{\rm ren}(R)=0.
 \tag{15}
$$

故770式(14)有明确局部解ℓ=ℓ_M−ℓ_T+ℓ_vol。这个来源按原作用同时变分全部物理背景字段，不是只修T而另给物质任意源。式(15)与任意同背景合法物理态差相容：

$$
 J_{\rm ren}[\lambda']-J_{\rm ren}[\lambda]
   =J_L[\lambda']-J_L[\lambda].
 \tag{16}
$$

局部接触不会抹掉实际准备或既有来源差。ℓ_vol是来源Wick接触，**不是**把b当作固定函数后声称δ∫b√|g|仅有体积项。一般作用变分还含δb；本轮未证明全部来源一形式的有效作用可积性。735已允许这种局部复合插入修复；不能与“加入不变作用的Euler导数无法修非守恒缺陷”的旧结论混淆。

同一接触可同时赋给770相等的实际tadpole与背景二次插入；二者之差仍为零，已证自由规范固定插入不被重新打开。是否能把这套处方扩展到所有插入与相互作用BRST身份，是另一个尚未完成的条件。

## 6. 归一、复算与适用边界

标量单束时w=1/2，M的体积接触为3νA δlog√|g|，(14)为−νA δlog√|g|，总和对应应力修正

$$
 \Delta T_{\mu\nu}=2\nu A g_{\mu\nu}
   =\tfrac13g_{\mu\nu}\langle:\!\varphi(-D)\varphi\!:\rangle,
 \qquad \langle:\!\varphi(-D)\varphi\!:\rangle=6\nu A.
 \tag{17}
$$

这与[Moretti Theorem 2.1](https://arxiv.org/html/gr-qc/0109048)的四维系数相符；它是符号/归一回检，完整构造来自上面的原耦合算符与辅助接触。

三组独立用途的复算如下。

1. 原111主部上，独立构造T、H_t及源变分，检查(2)—(5)的全部矩阵次序与两类接触。校准中的零核仅是齐次矩阵解，不是替换实际物理态；数据也不是原背景的Hadamard数值。
2. 平直Lorentz背景上的非交换线性矩阵势E(x)，实际积分输运递推：

$$
 v_0(x,y)=-\tfrac14(E(x)+E(y)),\qquad
 v_1(x,y)=-\tfrac12\int_0^1 t\,E(y+t(x-y))v_0(y+t(x-y),y)dt.
 \tag{18}
$$

高斯积分与独立复步导数核输运、伴随对称、左右jet及(11)。E矩阵不对易，不能以独立质量标量替代。

3. 用两矩阵束的周期局部jet验证w₁=1/2、w₀=−1及(13)—(14)的积分抵消，明确展示遗漏ghost或反转修复符号的非零残差。此周期参数校准不是紧时间物理宇宙模型。

第一组消元接触为0.0158506763，身份误差≤2.23×10⁻¹⁶。第二组势交换子范数0.0696，输运与jet残差≤2.09×10⁻¹⁷。第三组未修复Ward为−2.65042716×10⁻⁵，体积接触给相反值，合计残差≤2.72×10⁻²⁰；遗漏ghost时残差为−0.000565491。

本轮不数值求解完整753背景的重整化来源，也不由有限矩阵宣称连续结论。连续结论由局部递推、原算符映射与上述明确处方证明；有限常数可再加入同时守恒的许可项。适用范围限原F>0光滑补丁、紧支无穷小规范变换及当前一圈来源；不新处理全局规范反常或F=0。

## 7. 接回主线的下一步

735费米来源与本轮来源处方共同给出原在壳背景上的有限光滑联合来源，可作为754原完整约束/形式首阶响应的输入：

$$
 P\,w_1+j_{\rm ren}=0,\qquad K^\star j_{\rm ren}=0.
 \tag{19}
$$

式(19)只标出下一接口；原背景、初始修正、实际记录和同态来源变化须共同跟踪。接[772](../../772/drafts/STATUS.md)：先核这份绝对来源下的形式共同响应及其剩余局部常数，再区分它与自洽量子状态/实际过程之间的缺口。不重复局部余项或继续矩阵精度。目标保持。
