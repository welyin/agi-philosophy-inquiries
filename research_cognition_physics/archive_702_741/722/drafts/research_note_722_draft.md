# 第722轮：原平方参考的能源障碍与可容许读取

日期：2026-10-04。接[721](../../research_note_721.md)、[722入口](squared_reference_entry.md)。[代码](../joint_squared_reference_readout.py)、[结果](../joint_squared_reference_readout_results.json)、[核验](../research_round_722_checks.json)、[条件账](../unified_physics_condition_ledger_722.md)。三组复算；主代理审查，无新增独立代理复核。

## 1. 本轮关闭什么

652已经把原导数平方参考Q_f=W_f−V_f²放在完整有限图量子过程内，721已经证明原速度读口保有限能源。两者不能直接拼成“任何Q_f谱函数读口都保能源”。本轮完成两项相互约束的结果：

1. 一个具体平方相位读口，在原模型允许的一节点分支，把正规、Gauss不变、有限完整能源的态送出能源形式域。
2. 对包含原空间差分W_f的完整Q_f，构造另一种两结果读口，证明它在足够宽的响应参数下保持原完整能源域，并有能源域内的几何一阶导数。

第二项没有恢复第一种instrument的效果或后态；它回答的是同一个Q_f是否仍有实际、非平凡、有限能源的读取。因此不能把第一项扩大成所有平方参考不可读取，更不能扩大成空间或统一计划的反证。

|层次|本轮地位|
|---|---|
|认知动机|参考读数、记录后态与资源须属于同一过程|
|继承对象|598／623完整固定图H、Gauss、全部CAR与两类质量；652的T、s、V_f及Q_f|
|新增操作选择|指定有理Kraus与响应参数；不是新增认知公理，也不是自主装置|
|解析证明|原五维反例提升；完整Q_f的能源域保持、实际来源与有限历史|
|数值范围|原五维径向态积分；旧64维径向诊断加条件化邻居差分；不模拟完整Gauss热谱|
|边界|给定正外部几何、固定有限图；额外量子θ、空间一致极限和量子引力未覆盖|

554的非对易参考、560／652的平方压缩及704的热准备导数直接复用，不重复计为发现。空间382—386、425、522—524的既有条件及699的限定结论保持。

## 2. 把入口尾部提升到原五维Gauss态

采用623嵌入坐标x=φ/√F，写x=(x_H,z)、R=|x_H|、M=M₀²=2。原目标的坐标与测度直接拉回为

$$
D=1+(R^2+z^2)/6,\quad F=M/D,\quad
g_x=I-\frac{xx^{\mathsf T}}{6+|x|^2},\qquad
d\mu=D^{-1/2}\,d^5x .
\tag{1}
$$

固定w>0。721给原T速度场Y_T=(F/w)x_H·∂_{x_H}，保持z和角度。在R>0上

$$
\xi=\frac wM\left[A(z)\log R+\frac{R^2}{12}\right],\quad A=1+z^2/6,
\qquad d\mu=m\,d\xi\,dz\,d\Omega_3,\quad
m=\frac{MR^4}{w(A+R^2/6)^{3/2}},\qquad u=\sqrt m\,\psi .
\tag{2}
$$

每条非零径向轨道的ξ范围是整条实线，Y_Tξ=1，半密度变换后V_T=−iℏ∂ξ。此处使用652已经确定的自伴流，没有选择新的边界扩展。R=0是原目标内部的零测度集合；ξ坐标的负无穷端不等于删除Higgs零点或F=0边界。

取一节点无边分支，因而W_T=0；保全部32模式CAR，选其空真空作为态，不删去Hamiltonian中的Dirac／Majorana算符。令Y₀是归一常数S³谐波，b为支撑在|z|<1/2内的归一光滑函数，取

$$
u_0(\xi,z,\Omega)=\sqrt{3/2}\,(1-|\xi|)_+\,b(z)Y_0(\Omega)
 \otimes|0_{\rm CAR}\rangle .
\tag{3}
$$

它归一且Gauss不变。支撑中R有严格正下界和有限上界，ψ₀具有平方可积的一阶弱导数；原势、质量在支撑上有界。所以ψ₀属于原完整H的能源形式域，不要求它在H²域。空真空不需要是原H的基态或不变态。数值用紧支撑H¹多项式b替代光滑b，同样满足本证明。

## 3. 原能源对负ξ尾部有指数要求

四维Higgs径向的局部Hardy工具是成熟结论；可对照[Frank—Seiringer引言式(1.1)](https://arxiv.org/pdf/0803.0503)。这里N=4、p=2，直接由下式复核，不使用该文的分数阶主定理：

$$
\int_{\mathbb R^4}|\nabla v|^2-\int_{\mathbb R^4}\frac{|v|^2}{R^2}
=\int_{\mathbb R^4}\left|\nabla v+\frac{x_H}{R^2}v\right|^2\ge0 .
\tag{4}
$$

先在紧支撑光滑函数上积分，再按H¹闭包延拓，逐Fock分量和z积分。在R<2、|z|<1内，原度量／测度与Euclidean的比较常数有限。给ψ乘径向与z截止，截止在R≤1、|z|≤1/2上等于1；结合623完整H的形式等价，取A_H=H+c_H≥1，得到

$$
\int_{R\le1,\ |z|\le1/2}\frac{|\psi|^2}{R^2}\,d\mu
\le C\big(q_{\rm kin}[\psi]+\|\psi\|^2\big)
\le C' q_{A_H}[\psi] .
\tag{5}
$$

全部质量保留在形式等价中，不用纯玻色H替代完整H。该工具不同于579的全H⁵谱底／Hardy界，此次用途是内部Higgs零点附近的四维横向控制。

在ξ≤0、|z|≤1/2上R≤1，式(2)给

$$
R^{-2}=\exp\left[-\frac{2M\xi}{wA(z)}+\frac{R^2}{6A(z)}\right]
\ge e^{-2c_*\xi},\qquad c_*=\frac{M}{w(1+1/24)}>0 .
\tag{6}
$$

因此有限原能源要求u负ξ尾部的这一指数加权范数有限。

入口已经解析证明tent的平方流尾部。令U=exp(−iV_T²/ℏ²)，归一tent在ξ中的传播相当于单位自由Schrödinger参数，记传播结果为a₁。沿用入口的两次分段积分及余项界，得到

$$
|a_1(\xi)|^2\ge\frac{6\sin^2(1/4)}{\pi|\xi|^4},\qquad
|\xi|\ge \frac{33}{2\sin(1/4)} .
\tag{7}
$$

原tent是实函数，反向传播有相同模。式(5)—(7)使Uψ₀和U†ψ₀都不属于D(A_H^{1/2})。发散是原局部动能要求与严格多项式尾部不相容，不是数值截断拟合。选exp(−iV²/ℏ²)仅固定一个合法读口参数，不改变原H的ℏ。

对应的真实两结果instrument为

$$
J_\pm=(U\pm iU^\dagger)/2,\qquad
\sum_\pm J_\pm^\dagger J_\pm=I,\qquad
\sum_\pm J_\pm\rho J_\pm^\dagger
=\tfrac12(U\rho U^\dagger+U^\dagger\rho U) .
\tag{8}
$$

取ρ=|ψ₀〉〈ψ₀|，非选择后态的正移位能源无限，至少一个非零概率分支离开原能源域。可先以有界谱截断A_H验证能源等式，再用单调收敛；若两分支都在形式域，其线性组合Uψ₀也应在域，亦得矛盾。**有界、完备、Gauss合法的Kraus仍可不满足有限能源任务。**

## 4. 完整平方参考的另一种读取

现在回到任意固定有限图，f=T或s，保原652差分W_f≥0，并选固定上界b≥sup W_f。设能源Hilbert空间𝓔=D(A_H^{1/2})，范数为移位能源平方根。T、s及其目标梯度全域有界；有限图、正外部几何使W_f及全部目标梯度有界。乘法的Leibniz规则及623形式等价给

$$
R_f=b-Q_f=V_f^2+D_f\ge0,\qquad D_f=b-W_f,\qquad
d=\|D_f\|_{\mathcal B(\mathcal E)}<\infty .
\tag{9}
$$

这里b仅是对同一Q_f的可逆标号平移，不是丢弃配置的阈值；W_f不与V_f对易也没关系。正几何紧集内可统一选b和d。额外量子θ若存在，其导数条件另核，本文仍沿721排除该扩展。

721的双向能源界可吸收多项式因子，给某些固定图常数C_*≥1、ω>0：

$$
\|e^{-i\tau V_f/\hbar}\|_{\mathcal B(\mathcal E)}
\le C_*e^{\omega|\tau|} .
\tag{10}
$$

光滑核稠密、局部一致有界和原光滑流给𝓔上的强连续群。其Laplace积分直接给

$$
\|(V_f-\alpha)^{-1}\|_{\mathcal B(\mathcal E)}
\le\frac{C_*}{|\operatorname{Im}\alpha|-\hbar\omega},
\qquad |\operatorname{Im}\alpha|>\hbar\omega .
\tag{11}
$$

取α²=iλ、λ>0，用两个一阶预解式相减，得到

$$
\|(V_f^2-i\lambda)^{-1}\|_{\mathcal B(\mathcal E)}
\le\eta_\lambda=
\frac{C_*}{\sqrt\lambda(\sqrt{\lambda/2}-\hbar\omega)}=O(\lambda^{-1}) .
\tag{12}
$$

当dη_λ<1，在𝓔上作Neumann级数，Z₀=(V_f²−iλ)⁻¹满足

$$
(R_f-i\lambda)^{-1}=Z_0(I+D_fZ_0)^{-1},\qquad
\|(R_f-i\lambda)^{-1}\|_{\mathcal B(\mathcal E)}
\le\frac{\eta_\lambda}{1-d\eta_\lambda} .
\tag{13}
$$

这一逆确实等于原物理Hilbert空间中的预解式：对𝓔输入，所得向量由Z₀落入D(V_f²)，满足原(R_f−iλ)方程；原自伴算符的逆唯一。没有在形式层面擅自交换D_f和V_f。一组便于检查、无最优性主张的充分条件为

$$
\lambda>\max\{8\hbar^2\omega^2,\ 4\sqrt2 C_*d\} .
\tag{14}
$$

定义同一Q_f的实际Kraus，而不仅是效果：

$$
K_0=\lambda(\lambda+iR_f)^{-1},\qquad
K_1=iR_f(\lambda+iR_f)^{-1}=I-K_0 .
\tag{15}
$$

二者保Gauss、CAR宇称；它们的物理Hilbert范数不超过1，而且

$$
E_0=K_0^\dagger K_0=\frac{\lambda^2}{\lambda^2+R_f^2},\qquad
E_1=\frac{R_f^2}{\lambda^2+R_f^2},\qquad E_0+E_1=I .
\tag{16}
$$

它们不是Lüders平方根Kraus，不可删除相位。由式(13)得N_λ=λη_λ/(1−dη_λ)，对于每个正常有限能源态

$$
\sum_{j=0}^1\operatorname{Tr}A_HK_j\rho K_j^\dagger
\le\{N_\lambda^2+(1+N_\lambda)^2\}\operatorname{Tr}A_H\rho<\infty .
\tag{17}
$$

故所有非零概率的条件后态都仍有限能源；低概率条件分支没有统一的归一后能源上界。本轮没有借大λ使状态截断，也没有从有限矩阵全有界性推出式(17)。

## 5. 这仍在读取原Q，而非另换一个参考

K_j与R_f的所有谱投影对易，所以非选择读后保原Q的整个单变量谱分布。记录0概率非恒定；若可获得所有λ>λ_min的精确概率，则

$$
\operatorname{Tr}\!\left[P_{R_f}(B)\sum_jK_j\rho K_j^\dagger\right]
=\operatorname{Tr}[P_{R_f}(B)\rho],\qquad
\frac{p_0(\lambda)}{\lambda^2}
=\int_{[0,\infty)}\frac{d\mu_{R_f^2,\rho}(x)}{\lambda^2+x} .
\tag{18}
$$

右边为有限正测度的Stieltjes变换：在负实谱之外解析，一个开区间上的精确值通过恒等定理及边界反演唯一决定测度。R_f≥0又使R_f²与Q_f=b−R_f的谱标号一一对应。这只证明该整族概率不必丢失Q的单变量分布信息；**不是有限菜单、有限样本或稳定反演结论**。单个二值读数当然不能重建任意分布。

这也不意味着读取无扰动：Q之外的相干、原H、其它参考以及后续等待一般改变。完整量子态、不对易参考的联合任务仍不能用单变量分布替代。

## 6. 同一原几何来源可以保留到总导数

令外部几何参数为γ，在固定Hilbert识别和正参数紧集中工作。沿652，V_f=w(γ)⁻¹𝒱_f，𝒱_f不依赖γ；W_f及其目标梯度对γ连续可微。固定b、λ，记ζ=∂γlog w，则

$$
\partial_\gamma R_f=-2\zeta V_f^2-\partial_\gamma W_f .
\tag{19}
$$

定义Z=(λ+iR_f)⁻¹。由于R_fZ=(I−λZ)/i，且V_f²Z=R_fZ−D_fZ，右边是𝓔有界算符。故(∂γR_f)Z也在𝓔有界。原共同D(V_f²)、预解式恒等式和紧集上一致界给算符范数可微性，不仅是光滑核上的形式微分：

$$
\partial_\gamma Z=-iZ(\partial_\gamma R_f)Z,\qquad
\partial_\gamma K_0=-i\lambda Z(\partial_\gamma R_f)Z,
\qquad \partial_\gamma K_1=-\partial_\gamma K_0
\quad\hbox{in }\mathcal B(\mathcal E) .
\tag{20}
$$

原G=∂γH是𝓔有界形式。因此对固定有限能源ρ，或按704在能源夹权迹范数中可微的ρ(γ)，读后完整能源差的总导数有意义：

$$
\partial_\gamma\Delta E
=\operatorname{Tr}\rho\left(\sum_jK_j^\dagger G K_j-G\right)
+2\operatorname{Re}\sum_j\operatorname{Tr}\rho K_j^\dagger H(\partial_\gamma K_j)
+\operatorname{Tr}(\partial_\gamma\rho)\left(\sum_jK_j^\dagger HK_j-H\right).
\tag{21}
$$

含H的乘积都按𝓔上的双线性形式解释，不能要求每个K_jψ都在H的算符域。原Gibbs准备的该一阶条件由704复用。此处补齐的是**式(15)新读口**的域，不反向宣称721相位读口对任意有限能源态也有能源总导数。

有限次读取与同一原H的真实等待组合，沿721相同的结果树论证，若第a次读口的式(17)常数为C_a，则

$$
\sum_{\boldsymbol j}\operatorname{Tr}A_HK_{\boldsymbol j}\rho K_{\boldsymbol j}^\dagger
\le\left(\prod_a C_a\right)\operatorname{Tr}A_H\rho .
\tag{22}
$$

不丢罕见结果、不重置参考态。式(21)当前签收立即读口；含变化等待的总来源继续按704的既有域合同，不能只凭式(22)扩展所有历史导数。

## 7. 可复算结果与范围审查

第一组在原五维测度做两套积分：流坐标积分与独立径向积分的最大差3.56×10⁻¹⁵。w=1、沿旧数值ℏ=0.7，式(3)多项式b版本归一，原动能4.534043758、原正势0.150029276，均有限。空CAR真空的质量期望为零，但质量算符及形式约束仍保留。尾部发散由式(5)—(7)证明，有限区间的对数下界只作校准。

第二、三组沿652的64维径向表示，给f固定邻居值0.2、单位格距以诊断非零W=(f−0.2)²，并在共形比较下令W(γ)=e^(−4γ)W(0)。这是**条件化邻居参考诊断**，不是完整双节点Gauss模型；完整图结论靠式(9)—(17)。H仍是旧原径向H，未把参考W当作新相互作用加进H。

在γ=0.017、λ=16时：

|读口|记录0概率|漏W造成的概率差|实际加能|完整几何导数|其中仪器变化项|
|---|---:|---:|---:|---:|---:|
|Q_T|0.805082676|0.006694141|0.623868487|−16.426608011|−10.008322706|
|Q_s|0.878806006|0.010799974|0.035371340|−1.046522183|−0.677940890|

Kraus完备与谱分布身份误差≤9.04×10⁻¹⁶；独立几何差分与三来源之和误差≤7.67×10⁻⁸。两次读取、原等待和全部四条记录分别验算归一与本有限矩阵自身能源预算。λ=16是诊断参数，不被冒称已数值求得原全图式(14)的充分常数。

## 8. 减少的缺项与下一步

当前减少的是“原导数平方参考能否有保资源的实际读取”这一独立缺项：答案在固定图上肯定，同时具体相位读取存在严格反例。读取参数、原差分处方、固定几何和准备仍为输入。没有生成原H、三维、规范群、Einstein作用或自主操作权限。

下一项接[723](../../723/drafts/STATUS.md)：停止更多谱函数、尾部和常数优化，回查554、647—652、704及707—708，检查**同一完整参考任务**所需的T、s、Q_T、Q_s实际记录能否沿既有共同态／区域映射输送，以及哪些额外量词才足以复用旧条件性空间结论。单参考合法不等于四参考无扰动联合锐读，也不等于坐标重建；成熟结论的直接应用只记入口，不另加轮次。

保持384已消去Lipschitz、386与425的替代连接、522—524热参考及实际方向仪器的精确范围。699限定失败与其余共同分支不变，统一目标继续开放。
