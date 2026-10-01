# 第656轮：原空间手征测度与共同传播投影子的连接

日期：2026-10-02。接[655](research_note_655.md)、[原空间入口](round656_drafts/actual_spatial_measure_entry.md)。[代码](joint_spatial_auxiliary_geometry.py)、[结果](joint_spatial_auxiliary_geometry_results.json)、[核验](research_round_656_checks.json)、[条件账](unified_physics_condition_ledger_656.md)。

## 1. 本轮解决什么

655把原完整有限图动力学、固定原态记录和常时间来源接到同一个正分割，但没有证明该分割等于原空间overlap手征测度。656回到606／615的原矩阵，恢复真实空间依赖，不另选一个好用的转移核。

**实际新连接：**在明示自由规范背景下，原辅助Pfaffian的全部S⁹配置之绝对值，等于同一手征投影子压缩一个内部酉矩阵所产生的行列式。其非线性变化有由同一投影子固定的全局界；九个切向方向共享一个正图Laplacian。恢复空间传播会同时改变时间耦合，不能把旧时间链和独立空间项随意拼接。

共同二维内部平面上还得到带符号Pfaffian的精确非负恒等式，包含零点；一般S⁹只证明绝对值身份，未证明全配置相位或反射正性。

上一目标回合只核对并说明研究顺序，没有改变科学证据，按无新增进展处理。本轮重新核验655全部2218份保护证据并复算655，未发现运行中的Python。回查615、653、655和历史关键词；旧Slater重叠、矩阵对数和overlap局域性是成熟工具，不另算新定理。本轮新增的是它们与原空间辅助测度、原时间链及空间来源的具体连接。

|层次|本轮内容|
|---|---|
|认知动机|同一过程应共同承担传播、内部变量及其响应|
|继承输入|四维Euclidean Wilson结构、m₀=1、原16内部表示、B、十个T矩阵、辅助S⁹处方|
|适用分支|单位规范链路、零物理标量、有限周期空间与反周期时间；不推广到一般规范场|
|额外诊断参数|空间Wilson系数κ，κ=1是原值，其他值只用于追踪共同来源|
|解析结论|压缩恒等式、全S⁹模长界、平面精确正性、Hessian及局域尾界|
|数值作用|独立原Pfaffian、压缩行列式、矩阵导数及解析小格结果互相核验|
|未签收|全S⁹积分、实际物理态映射、正时间拼接、空间连续极限和量子引力|

本轮对接[Kikukawa原文§4.2—4.5](https://arxiv.org/html/1710.11618v3#S4.SS2)。该文讨论原辅助处方；其中自由背景全配置正性的论证还使用实性和零集假定及数值观察。本轮不把它升级为无条件定理，也不以模长代替原相位。未做图像检查。

## 2. 原投影子与内部矩阵

令X为606的自由Wilson矩阵，H为其Hermitian化，V为负谱正交帧。暂把内部16维指标剥离；P作用于格点与四分量spin：

$$
H=\gamma_5X=H^\dagger,\qquad
P=1_{(-\infty,0)}(H)=VV^\dagger,\qquad V^\dagger V=I.
\tag{1}
$$

实移位矩阵包含反周期边界符号。615的B为反对称酉矩阵，十个T_a为对称矩阵。取T₀作参考，K₀=I，K_a=T₀†T_a，a=1,…,9。原Clifford关系给

$$
B_{\rm lat}^\dagger P^*B_{\rm lat}=P,\qquad
K_a^\dagger=-K_a,\quad \{K_a,K_b\}=-2\delta_{ab}I_{16},\quad
\operatorname{tr}K_a=0,\quad \operatorname{tr}(K_aK_b)=-16\delta_{ab}.
\tag{2}
$$

第一式直接由B†γ_μ*B=γ_μ与B†γ₅*B=γ₅及谱函数得到，B_lat=I_sites⊗B。没有要求P逐格对角。令E_x∈S⁹，e₀=(1,0,…,0)，u=V⊗I₁₆：

$$
A(E)=u^T\operatorname{diag}_x\!\left(B\otimes\sum_{a=0}^{9}E_x^aT_a\right)u,
\quad A_0=A(e_0),\quad f(E)=\frac{\operatorname{Pf}A(E)}{\operatorname{Pf}A_0}.
\tag{3}
$$

使用同一帧的比值消除任意固定基底相位。式(2)使VᵀB_lat V酉，故A₀酉且无零点。设Π_x是格点投影，F_x=V†Π_xV，Σ_xF_x=I。直接乘A₀†给

$$
D(E):=A_0^{-1}A(E)=\sum_xF_x\otimes K(E_x)
=\mathcal V^\dagger\mathcal U(E)\mathcal V,
\quad K(E)=E^0I+\sum_{a=1}^9E^aK_a,\quad
\mathcal V=V\otimes I_{16},\quad \mathcal U=\operatorname{diag}_x(I_4\otimes K(E_x)).
\tag{4}
$$

K(E)†K(E)=I由式(2)及|E|=1得到，因此𝒰是酉矩阵。这是原辅助矩阵的具体身份，未更换测度或增加传播规则。

## 3. 全部S⁹配置：模长、漏出与全局界

记𝒫=P⊗I₁₆，R=(I−𝒫)𝒰𝒱，L=R†R。酉压缩给D†D=I−L，同时f²=detD。于是

$$
|f(E)|=\det(I-L)^{1/4}\le1,\qquad 0\le L\le I.
\tag{5}
$$

这里L是从原负谱子空间被内部变换带出之量的Gram矩阵，不是物理系统真实丢失概率。若||L||<1，矩阵对数逐本征值给

$$
\frac14\operatorname{Tr}L\le-\log|f(E)|
\le\frac{\operatorname{Tr}L}{4(1-\|L\|)}.
\tag{6}
$$

有零点时下界仍成立，−log0=+∞，上界不作有限断言。内部迹满足tr(K(E_x)†K(E_y))=16E_x·E_y。利用P²=P，对所有配置得到

$$
\frac14\operatorname{Tr}L
=\frac12\sum_{x<y}w_{xy}|E_x-E_y|^2,
\qquad w_{xy}=8\|P_{xy}\|_F^2\quad(x\ne y).
\tag{7}
$$

因此同一传播P对整个S⁹被积式的模长给出非线性约束，并非只在二维平面或小扰动附近有效。不过不能用|f|替代f去做原路径积分；相位还须另核。

## 4. 共同二维内部平面：精确正性

取E_x=(cosθ_x,sinθ_x,0,…)，D_θ=diag_x(e^{iθ_x}I₄)，M_θ=V†D_θV。K₁有±i本征值，各八重。式(4)按这两支分解，得到

$$
f(\theta)=|\det M_\theta|^8
=\det(M_\theta^\dagger M_\theta)^4\ge0.
\tag{8}
$$

证明不能只说“从参考态连续延拓所以不会变号”：这会漏掉零点。严谨做法是把两边看成角变量的有限Laurent多项式；两边的平方都等于detD。Laurent多项式环无零因子，故它们全局相等或全局相反。参考θ=0时均为1，确定全局正号，包括零点。

满模条件也无需小角假设。有限维酉压缩所有奇异值≤1，行列式模为1当且仅当原子空间被保持，等价于𝒰与𝒫对易：

$$
|f(E)|=1\quad\Longleftrightarrow\quad
P_{xy}\otimes[K(E_y)-K(E_x)]=0\ \text{对所有 }x,y
\quad\Longleftrightarrow\quad E_x=E_y\ \text{在每条 }w_{xy}>0\text{ 的边上}.
\tag{9}
$$

最后一步用K的实线性单射性。若该图连通，则仅整体共同方向达到满模；在平面上等价于角度模2π相同。它不证明原S⁹积分处于有序相：配置熵、体积极限和相位尚未处理。平面配置在一般S⁹配置空间内测度为零，式(8)不能代替全S⁹积分。

式(8)也可表示成填满RanP的Slater态与D_θ作用后的重叠模之八次方。这个Slater态使用Euclidean格点／spin空间，**不是**598原物理32CAR态；维数和迹公式相似不构成物理态识别。

## 5. 九个切向方向共享同一个Laplacian

在共同参考方向附近用E_x=(sqrt(1−|y_x|²),y_x)参数化，y_x∈R⁹。Pfaffian无零的邻域内，

$$
d\log f=\tfrac12\operatorname{Tr}(D^{-1}dD),\qquad
d^2\log f=\tfrac12\operatorname{Tr}(D^{-1}d^2D-D^{-1}dD\,D^{-1}dD).
\tag{10}
$$

代入D(0)=I、式(2)、ΣF_x=I，得到一阶为零和以下实Hessian。tr(F_xF_y)=||P_xy||²_F，故

$$
\partial_{y_x^a}\partial_{y_y^b}\log|f|\big|_0
=\delta_{ab}\begin{cases}
8\|P_{xy}\|_F^2,&x\ne y,\\
8[\operatorname{tr}(P_{xx}^2)-\operatorname{tr}P_{xx}],&x=y
\end{cases}
=-\delta_{ab}(L_w)_{xy},\quad L_w=\operatorname{diag}(w\mathbf1)-w.
\tag{11}
$$

对角恒等式是投影条件Σ_yP_xyP_yx=P_xx的结果，不能独立调节。于是

$$
-\log|f(E)|=\frac12\sum_{x<y}w_{xy}|y_x-y_y|^2+O(\|y\|^3),
\qquad \ker(L_w\otimes I_9)\text{每连通分量有九个共同方向零模}.
\tag{12}
$$

有限图上余项是局部Taylor余项，未给无限体积一致控制。这是一个辅助权重的局部二次作用，不是新量子理论、光锥或三维空间的推导。

## 6. 恢复空间会同时改写时间部分

取三个周期空间格、两个反周期时间格，另两空间方向各一个格。空间Wilson项乘κ。k=0,±2π/3，ω=±π/2：

$$
H(k,\omega)=\gamma_5[\kappa(1-\cos k)I+i\kappa\sin k\,\gamma_1+i\sin\omega\,\gamma_4],
\qquad r(k)^2=1+2\kappa^2(1-\cos k).
\tag{13}
$$

P(k,ω)=(I−H/r)/2。直接有限Fourier求和，r=sqrt(1+3κ²)，三类非对角权重为

$$
w_s=\frac{4\kappa^2}{1+3\kappa^2},\qquad
w_t=\frac89(1+2/r)^2,\qquad
w_\times=\frac89(1-1/r)^2.
\tag{14}
$$

s为同时间不同空间；t为异时间同空间；×为异时间异空间。κ=0给三条独立时间链，每条权重8；κ=1给w_s=1、w_t=32/9、w_×=2/9，每行度数6：

$$
\operatorname{spec}L_w\big|_{\kappa=1}=
\{0,11/3,11/3,8,31/3,31/3\}.
\tag{15}
$$

同一空间系数的导数不是只作用空间权重：

$$
\partial_\kappa(w_s,w_t,w_\times)\big|_1=(1/2,-8/3,1/3).
\tag{16}
$$

一般有隙有限矩阵也可不选择可微特征帧，直接对谱投影求导。H=Σ_iλ_i|i〉〈i|，p_i=1_{λ_i<0}，则

$$
(P')_{ij}=\frac{p_i-p_j}{\lambda_i-\lambda_j}(H')_{ij}\quad(p_i\ne p_j),
\qquad (P')_{ij}=0\quad(p_i=p_j),\qquad
w'_{xy}=16\operatorname{Re}\operatorname{tr}(P_{xy}^\dagger P'_{xy}).
\tag{17}
$$

前两式在H本征基书写，末式回到位置基。负／正谱之间有隙，同侧简并不导致除零。这条来源响应由原H决定；κ不是物理空间度规，不将它的导数直接命名为完整应力。

现在令每一时间片内部空间方向完全相同，但两个时间片的角度相差δ。任何归一为“空间均匀时等于1”的纯同时间空间乘子均无法改变旧时间链权重。然而原Pfaffian给

$$
-\log f=\begin{cases}
12\delta^2+O(\delta^4),&\kappa=0,\\
6\delta^2+O(\delta^4),&\kappa=1.
\end{cases}
\tag{18}
$$

因此“旧时间链乘一个上述纯空间修正”这类接法被排除，连这个空间均匀切片都不匹配。此结论强于入口单点不等，但仍不排除带混合时空项、额外状态或记忆的正转移。被积式差异不直接证明积分后的所有观测都不同。

## 7. 耦合尾由同一原谱隙控制

在原各向同性四维自由格上，a_μ=1−cosk_μ∈[0,2]。Clifford关系直接给

$$
H(k)^2=\left[1+2\sum_{\mu<\nu}a_\mu a_\nu\right]I,
\qquad 1\le |H|\le7.
\tag{19}
$$

由(I−H²/25)的范数≤24/25，把signH写成H/5乘(1−z)^(−1/2)的二项级数；第n系数binom(2n,n)/4^n≤1。截断到m阶的范围≤2m+1条原格边，算符误差≤(7/5)(24/25)^(m+1)/(1−24/25)。对距离d≥2的不同格点取m=floor((d−2)/2)，并用四分量块的Frobenius范数≤2倍算符范数，得到

$$
w_{xy}\le 9800\left(\frac{24}{25}\right)^{2\lfloor d(x,y)/2\rfloor},\qquad d(x,y)\ge2.
\tag{20}
$$

这是很松但与图大小无关的上界；它把原谱隙传到同一个辅助Hessian。不是最近邻截断、物理因果界或维数选择。成熟的overlap局域性不作为新的独立轮次；这里需要的是与式(7)、(11)的明确对象映射。

## 8. 复算及验收

命令为现有Python运行 `joint_spatial_auxiliary_geometry.py`；已保存结果后默认只复算比较，`--write-results`只允许独占新建，避免覆盖。

四组检查：

1. 原Pfaffian与平面行列式身份，在κ=0、0.5、1及两种时间／空间格长上共12个一般配置核验；另验零点和统一方向。相位直接算Pfaffian，不从行列式随便取平方根。
2. 全部九个内部Clifford关系；三个κ的原完整矩阵Hessian；一般九方向组合的原Pfaffian有限差分。
3. 五个κ的三类权重、κ=1的解析谱、谱投影导数和有限差分；空间均匀时间切片给12与6的二次系数。
4. 平面及三组全S⁹配置的压缩漏出恒等式和上下界；独立四维Fourier符号核验。全S⁹三例−log|f|约0.525273、7.267241、18.817182，对应下界0.521719、6.474637、13.862153。全S⁹模长身份误差≤1.43×10⁻¹⁴。

初版数值仅含平面漏出；后续解析推导把它扩展到全S⁹模长，故保留[初版结果](round656_drafts/initial_planar_results.json)，正式结果覆盖新增身份但不覆盖旧文件。本轮没有靠改参数掩盖失败，也没有执行图像检验或独立代理审查。

**条件压缩：**在这个明确分支上，辅助变量之间的空间权重、时间权重、混合权重及其κ来源不再是四套独立输入，都由同一个P固定。C16测度因此与C04传播和C22来源发生实质连接。

**仍未消去：**原四维格、群和物种、辅助处方、物理参考态、S⁹积分相位、反射正性及空间连续映射。不得把一个配置的最大模当成物理真空，不得把Slater压缩空间当成原CAR物理空间，不得把自由辅助一致性当成标准模型和GR已经统一。

接[657完整积分与时间拼接](round657_drafts/STATUS.md)，先检验原多空间测度是否支持同一正过程；不继续优化κ或辅助装置，认知设计后置，总目标不变。
