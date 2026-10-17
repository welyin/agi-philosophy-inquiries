# 1268材料分项：同齐次读取、异空间主部

2026-10-16。按[准入](admission.md)交付一个作者科学组。CT未采用；不修改1267或旧模型。本文与[读取分项](readout_proof.md)分别验收，不能把一种材料的来源和另一种材料的量子读取拼成共同证书。

**核心命题：**两种明确协变延拓可以给相同的齐次尺度Hamiltonian、密度和压力，却给不同的杆—数流主部。原几何依赖粒子质量产生纵向杆—无压流混合；只把读量改为线性杆标量，不能消去该混合。直接调杆动能系数则避开这项数密度混合，但必须公开其新的控制/记录材料身份。无截断主部异常不自动成为原有限频带任务失败。

## 1. 原作用和诊断量词

[1245式(2)、(6)、(12)](../1245/proof.md)的原杆和第一流来自

$$
\mathcal L_{\rm rod}=-\frac{\gamma}{2}\mathcal B,\qquad
\gamma=F(p)=1+p/5>0,\qquad
\mathcal B=\sum_Ag^{\mu\nu}\partial_\mu\phi^A\partial_\nu\phi^A,
$$

$$
S_c=\int d\mu_c\int(\Theta_c-M\,d\tau_c).
\tag{1}
$$

1247加入的标量、振子、程序和记忆保留。[1249](../1249/canonical.md)的齐次约化不替代全空间材料资格。保持窗原脉冲关闭，第一流内部质量为正，原三杆空间梯度非零。

在一个局部共动正交标架的**冻结系数jet**取

$$
g_{\mu\nu}=\operatorname{diag}(-1,1,1,1),\quad
U=(1,\mathbf0),\quad n>0,\quad
\partial_i\phi^A=s\delta_i^A,\quad
\partial_t\phi^A=0,\quad \mathcal B_0=3s^2 .
\tag{2}
$$

$s=k/a$是当地杆梯度，后文$q$是扰动波数，不与引力耦合$\kappa$混同。冻结局部系数是计算主符号的方法，不是假定整个宇宙、自由几何或实际读取窗口静止。

原候选质量项为

$$
M=M_{\rm base}+c\,f(\mathcal B),\qquad
c=\varepsilon\alpha g(\theta)p_*,\quad p_*=\pm1,\quad
\bar M=M_{\rm base}+cf(\mathcal B_0)>0.
\tag{3}
$$

$p_*$是保持的指针$\sigma_z$部门，不把$|+x\rangle$的零平均替代两个部门。这里$f$是$\mathcal B$的函数，与前筛查的$q_{\rm rod}$函数经链式法则相连；原直接读$a$时$f(\mathcal B)=\sqrt{3k^2/\mathcal B}$。在平台内求导，完成区另需自己的资格。

## 2. 固定数流下的杆块

保持$n,c,\gamma$不变，杆方程为

$$
\partial_\mu\bigl[(\gamma+2nc f_{\mathcal B})
  \partial^\mu\phi^A\bigr]=0.
$$

其二阶主矩阵是

$$
\mathcal P_{AB}^{\mu\nu}
=A g^{\mu\nu}\delta_{AB}
+4nc f_{\mathcal B\mathcal B}
 \partial^\mu\phi^A\partial^\nu\phi^B,\qquad
A=\gamma+2nc f_{\mathcal B}.
\tag{4}
$$

沿波矢分解两个横向杆模和一个纵向杆模。横向速度平方为1；纵向为

$$
v_{\rm rod,L}^2=\frac{D}{A},\qquad
D=A+4nc f_{\mathcal B\mathcal B}s^2 ,
\tag{5}
$$

前提$A>0,D>0$。对$f=\sqrt{3k^2/\mathcal B}$，

$$
f_{\mathcal B}=-\frac{f}{2\mathcal B},\qquad
f_{\mathcal B\mathcal B}=\frac{3f}{4\mathcal B^2},
\qquad
A=\gamma-\frac{ncf}{\mathcal B},\qquad D=\gamma.
\tag{6}
$$

所以$c>0$且$A>0$时固定数流纵向杆块快于原度规；小$c$不使精确速度重新等于1。对线性$f=\mathcal B/(3k^2)$，$f_{\mathcal B\mathcal B}=0$，该**固定数流杆块**确实回到原锥。下一节说明这还不是完整材料的结论。

## 3. 保原无压流：不能删除密度—杆混合

用守恒标签的实际位移$\boldsymbol\xi$表示第一流：

$$
X(t,\mathbf z)=\mathbf z+\boldsymbol\xi(t,\mathbf z),\qquad
\phi^A(t,\mathbf x)=s x^A+\varphi^A(t,\mathbf x).
$$

从同一世界线作用
$-\int n\,d^3z\,M(\mathcal B(t,X))\sqrt{1-|\dot X|^2}\,dt$
展开二阶主项。标签测度$n\,d^3z$固定，Eulerian数密度的变化没有冻结。用

$$
\delta_1\mathcal B=2s\,\nabla\!\cdot\!\boldsymbol\varphi,\qquad
\delta_2\mathcal B=|\nabla\boldsymbol\varphi|^2-
|\dot{\boldsymbol\varphi}|^2
$$

以及位置取值的项
$\boldsymbol\xi\cdot\nabla\delta_1\mathcal B$。对紧支或周期扰动分部积分后，

$$
\boxed{
\mathcal L^{(2)}_{\rm prin}
=\frac A2\bigl(|\dot{\boldsymbol\varphi}|^2
-|\nabla\boldsymbol\varphi|^2\bigr)
-2nc f_{\mathcal B\mathcal B}s^2
 (\nabla\!\cdot\!\boldsymbol\varphi)^2
+\frac K2|\dot{\boldsymbol\xi}|^2
+C(\nabla\!\cdot\!\boldsymbol\xi)
  (\nabla\!\cdot\!\boldsymbol\varphi),
}
\tag{7}
$$

$$
K=n\bar M>0,\qquad C=2nc f_{\mathcal B}s.
$$

交叉符号可直接核：
$-2nc f_{\mathcal B}s\,\boldsymbol\xi\cdot
\nabla(\nabla\cdot\boldsymbol\varphi)$积分后为式(7)的正号。无压流没有凭空获得正的$(\nabla\cdot\xi)^2$恢复项。为修复而添加压缩能或其它材料会改变候选，不是此推导允许的省略项。

取纵向平面波$e^{i(qx-\omega t)}$、$y=\omega^2/q^2$，由式(7)得到

$$
\begin{pmatrix}D-Ay&-C\\-C&-Ky\end{pmatrix}
\binom{\varphi_L}{\xi_L}=0,
\qquad
AKy^2-DKy-C^2=0.
\tag{8}
$$

若$A,K,D>0$且$C\ne0$，两根

$$
y_\pm=\frac{D\pm\sqrt{D^2+4AC^2/K}}{2A}
$$

满足$y_-<0<y_+$。因此不是只有固定流背景下的轻微锥偏移：这个实际无压流—杆主块还含负的纵向速度平方。对线性$f$，$D=A$而多项式在$y=1$处为$-C^2$，故

$$
y_+>1,\qquad y_-<0
\quad\text{只要 }c f_{\mathcal B}\ne0 .
\tag{9}
$$

这排除“改成线性$\mathcal B$就自动保留完整原共同锥”的推断。

### 3.1 为什么没有以冻掉物质制造异常

本节保留原第一流位移及其质量惯性，连续方程由标签表示自动满足。原内部钟满足$\dot\theta=1$，指针部门保持；在这一均匀共动jet的非零频率主模可取$\delta\theta=\delta P=0$。$g'$引起的$J$方程、原内部振子及记忆方程是低阶传输/有限频率项，不能作为额外的正密度梯度能加入式(7)。

原1245§3已要求位置变分固定标签内部动作；不能把不同标签的$\nabla J$另加为机械力，以取消式(7)。其它标量与杆的$F_p$交叉项在该二阶位移/场表示中至多含一阶扰动导数。所写作用不含度规导数，谐和规范下度规二阶主块与式(7)之间的混合也至多一阶；原场/几何不会提供一个漏掉的同阶$\xi$—$\varphi$抵销。

这里给的是明确背景和扰动类的材料主符号证据，不是所有非均匀背景的全PDE分类、非线性存在性定理或实际量子通道因果反证。若对流体重新施加刚性固定、增加弹性/压强或采用不同材料作用，须重算，不能沿用(8)或沿用原无压模型名称。

## 4. 有限波数/窗口：UV异常不是有限任务失败

对常系数块，负根的增长率为

$$
\Gamma(q)=|q|\sqrt{
\frac{\sqrt{D^2+4AC^2/K}-D}{2A}}
=|q|\sqrt{
\frac{2C^2}{K(\sqrt{D^2+4AC^2/K}+D)}}.
\tag{10}
$$

由$D>0$，

$$
\Gamma(q)\le \frac{|q|\,|C|}{\sqrt{KD}}.
\tag{11}
$$

因此无截断$q\to\infty$的常系数方程有无界快增长，不能直接取得原强双曲/共同锥资格。然而在事前固定$|q|\le q_{\max}$、持续时间$T$内，这一块的指数增长因子至多

$$
\exp\!\left[\frac{q_{\max}T|C|}{\sqrt{KD}}\right].
\tag{12}
$$

$c=\varepsilon\alpha gp_*$时，局部率按$O(|q|\varepsilon|\alpha g|)$缩小。式(12)只控制这个**常系数增长因子**，不是完整传播算符范数、准备用量或实际记录偏差。真实时变系数还要计基底变化、其它项及有限带准备/传播；不能把$\int\Gamma\,dt$直接称该完整问题的证书。正根导致的有限提前量、边界泄漏也须另与实际允许事件预算比较。

代码给一个正常正惯性、有限$q_{\max},T$的有理系数诊断，使(12)的指数小于$1/100$，同时$y_-<0$仍严格成立。该例仅显示逻辑范围不同；它不是原选定极小$\varepsilon$任务已失败或已实现的证据。当前只停止未经检查继承原无截断材料资格，不因高波数异常否定有限有效候选或阶段目标。

## 5. 一个更小的替代：直接调杆系数

与质量路线不同，公开选择一个体积型有效作用项

$$
\mathcal L_{\rm int}^{\rm vol}
=-d(\theta,P)\mathcal B,\qquad
d=\frac{\varepsilon\alpha}{3k^2}g(\theta)P.
\tag{13}
$$

这里不再乘$n$。在$P,\theta$当作有明确材料身份、且变分时固定的独立标量参数时，杆的主系数仅为

$$
\gamma_{\rm eff}=\gamma+2d,\qquad
\mathcal P^{\mu\nu}_{AB}=\gamma_{\rm eff}
g^{\mu\nu}\delta_{AB}.
\tag{14}
$$

一个同时覆盖$\pm1$部门的非退化条件是

$$
\frac{2\varepsilon|\alpha|\|g\|_\infty}{3k^2}
<\gamma_{\min}.
\tag{15}
$$

在归一旧$C^2$脉冲、$k^2=1/10$时，左边为
$175\varepsilon|\alpha|/(12\tau)$。这是明确有限不等式，不能把它与已经证明的完整因果实现混同。

**独立变分来源**为

$$
T_{\mu\nu}^{\rm vol}
=2d\sum_A\partial_\mu\phi^A\partial_\nu\phi^A
-d\mathcal B\,g_{\mu\nu}.
\tag{16}
$$

其作用没有度规导数，但控制/记录变量不是外部无源系数。合并杆后，杆方程给出的来源交换包含

$$
\nabla_\mu T_{\rm rod}^{\mu\nu}
=-\frac{\mathcal B}{2}\nabla^\nu\gamma
-\mathcal B\nabla^\nu d.
\tag{17}
$$

原$p$来源只抵销第一项；实际$\theta/P$和控制者必须抵销新增第二项，不能因为$P$简并就漏掉它。

如果$\theta/P$是原流的Eulerian内部值，原内部辛项按数密度加权，而式(13)不按$n$加权。还原为“每标签质量”会得到$d\mathcal B/n$，显含数密度；这超出1245的旧$M(\text{场值,内部量})$直接Legendre模板。若采用新的局部场控制或有限整体有效指针，同样须公开其辛项、来源、准备及空间字典。式(14)只核固定控制背景的杆块；在均匀控制jet下没有式(7)的密度混合，但尚未证明一般控制扰动的完整主部或CC。

## 6. 核心删除证据：同约化、同齐次来源，异空间响应

令$C_0=3k^2$，沿原齐次族

$$
n=a^{-3},\qquad \mathcal B=C_0a^{-2},\qquad
f(\mathcal B)=\sqrt{C_0/\mathcal B}=a,\qquad d=c/C_0 .
$$

质量型和体积型的坐标体积Hamiltonian密度恰相同：

$$
a^3ncf=ca,\qquad a^3d\mathcal B=ca.
\tag{18}
$$

独立来源亦在该齐次族相同：

$$
\rho_{\rm mass}=\rho_{\rm vol}=c/a^2,\qquad
P_{\rm mass}=P_{\rm vol}=-c/(3a^2).
\tag{19}
$$

即两者可给同一七对附加$\varepsilon\alpha g aP$门及同一齐次压力，却仍有第3、5节不同的空间延拓。只检验$a$读取和均匀能源账不能选择其中一份材料资格。

差别可直接由同一变分展示，不用新扫描。在固定$\mathcal B,g,\theta,P$下改变数密度，交互密度响应为

$$
\partial_n\rho_{\rm mass}=cf,\qquad
\partial_n\rho_{\rm vol}=0.
\tag{20}
$$

原基线材料对$\delta n$的共同响应在此相减后保留，未被删除。再取固定局部度规、固定$n$的杆剪切，令
$S_{ij}=\sum_A\partial_i\phi^A\partial_j\phi^A$，
$\operatorname{tr}\delta S=0$。两来源的无迹应力响应分别为

$$
\delta T^{\rm mass}_{ij}=2nc f_{\mathcal B}\delta S_{ij}
=-d\,\delta S_{ij},\qquad
\delta T^{\rm vol}_{ij}=2d\,\delta S_{ij}.
\tag{21}
$$

它们的差为$-3d\,\delta S_{ij}$。这是独立非齐次响应自由，不表示原有限噪声和菜单已经能分辨该差。

## 7. 验证、方向选择和保留范围

[material_check.py](material_check.py)登记一个作者组：

- 从二次作用的有限多项式系数提取$A,D,K,C$，核式(8)；用有理符号检查两个根异号、线性质量路线的$y_+>1$及有限增长上界。
- 在原$k^2=1/10$、若干有理$a,c$、正完整惯性下核(18)—(21)和体积型系数正性；示例只是代数/主部诊断，不代替原物理资源。
- 输出保存值中的根及增长率小数只作显示，连续符号结论由解析式及精确有理符号承担；不运行全PDE、量子波包或设备实验。

首次 --write-results 排他保存[material_results.json](material_results.json)，默认只读重算匹配。绑定准入、本文/代码及冻结父作用、候选与1267；不绑定同时编写的读取分项。

**方向结论：**暂停把数密度乘几何函数的质量读门作为已合格全空间材料；线性读$\mathcal B$也不能靠固定数流测试自动过关。若继续有限任务，须给明确带宽/时窗/泄漏与联合来源预算。可优先检查体积型杆系数候选，但它是公开的新材料延拓，不能静默替换原质量作用。读取分项如果只给共同齐次量子概率，本轮只承认该数学读取，不把它升级为完整局部材料实现。本轮不修复全部PDE、不新增认知公理、不宣布认知纲领或有效GR失败。
