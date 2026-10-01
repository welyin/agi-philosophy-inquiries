# 第569轮研究：共同Higgs、残余Gauss与全局商群上的磁作用补全

日期：2026-10-01。接续568；八组复算、十六式。代码joint_quotient_gauge_completion.py，结果joint_quotient_gauge_completion_results.json。正式完成以research_round_569_checks.json为准。

## 0. 结论及其对条件合并的意义

本轮将568的同一Higgs、singlet势和传播模型，与531—547已经给定的超荷、全局群候选及同尺度耦合连接起来。**这是一项有附加输入的共同实现，不是减少全部独立公理数的证明。**

得到两个实际接口结果：

1. 加入超荷后，原SU(2)完全Higgs化的约化不能直接照搬。完整图中留下一个电磁Gauss约束；正确约化得到一份有质量中性纵模，不产生额外光子纵模。
2. 如果选择旧表示允许的最大忠实Z6商群，568直接分写的基本磁迹不能下降为该商群上的函数。可以用旧物质表示的正权字符和加正补项，构造全局良定义、逐面只有单位类取零、且具有任意给定正二次系数的磁函数。

第二项消除了这一明确拼接障碍，但引入磁函数形状选择。相同二次传播不决定完整作用、量子相或全局真空。覆盖群支线仍成立，最大商群没有被531强迫选出。

## 1. 来源、去重及适用范围

复用[531](research_note_531.md)的荷、正权指标与中心核，[533](research_note_533.md)的连接归一，[544](research_note_544.md)第三个保存物质状态，[547](research_note_547.md)的已有中性质量，[553](research_note_553.md)的跨尺度限制，以及[566](research_note_566.md)—[568](research_note_568.md)的原物质、Gauss和传播。

[Tong，Line Operators in the Standard Model，arXiv:1705.01853](https://arxiv.org/abs/1705.01853)明确区分全局群的中心商选择及其线算符和theta周期。这里沿用其全局形式问题，不把论文2017年关于当时实验的表述当作2026年数据结论。本轮的具体有限图磁函数和约化按下文直接证明，不声称发现新标准模型质量公式。

|层次|本轮地位|
|---|---|
|认知动机|不同部门须共用同一物质、约束、传播和资源对象|
|既有模型输入|给定图和面胞，正势L、u，h★>0，既有荷表示及运行轨道|
|新增选择|色与超荷链路、全局商分支、正磁函数、经典系数匹配处方|
|解析证明|残余Gauss辛约化、磁函数下降障碍、正补全存在性和逐面零集合|
|数值验证|同尺度数值、完整图谱、非线性Hessian、独立链路lift与节点规范变换|
|仍未建立|手征费米子动力学、量子连续极限、动态几何连接、维数及统一完成|

## 2. 同一荷单位与有限模型

用整数Q=6Y表示周期2pi的圆群，连接z=exp(i alpha)。沿用531的N=3、非零超荷和Majorana分支，LH Weyl表示及Higgs为

$$
Q_L:(3,2)_1,\quad u^c:(\bar3,1)_{-4},\quad d^c:(\bar3,1)_2,
\quad L:(1,2)_{-3},\quad e^c:(1,1)_6,\quad \nu^c:(1,1)_0,
\quad X:(1,2)_3;
\quad G={SU(3)\times SU(2)\times U(1)_Q\over\langle\zeta\rangle},
\quad\zeta=(\omega_3 I_3,-I_2,e^{i\pi/3}).
\tag{1}
$$

本轮动力学只包含X（等价于原四实分量，通常物理H=X/√2）、s和三种规范链路。式(1)费米子表示用来定义磁函数，**没有因此实现手征费米子传播、Yukawa或完整反常量子化**。全部表示核是对角Z6；仅Higgs在SU(2)×U(1)中的核也是Z6，但不带颜色而包含全部旧电弱物质时共同核只有(-I,-1)生成的Z2，三个命题不可混同。

令颜色链路为C_e、弱链路为W_e，U(1)链路为z_e。统一Hamiltonian候选是

$$
H=\sum_i\left[{|p_{X,i}|^2+p_{s,i}^2\over2v}+{v\over4}d_i^TLd_i\right]
+\sum_e(b_c E_c^2+b_w E_w^2+b_0 E_0^2)
+{k_h\over2}\sum_{e:i\to j}|X_i-z_e^3W_eX_j|^2
+{k_s\over2}\sum_e(s_i-s_j)^2+\sum_f V_f(C_f,W_f,z_f),
\quad d_i=(|X_i|^2-u_h,s_i^2-u_s).
\tag{2}
$$

采用tr(t_at_b)=delta_ab/2；弱t3=diag(1,-1)/2。对体积v=epsilon^d的规则格图，一组明示经典匹配是

$$
k_s=k_h=c^2\epsilon^{d-2},\quad
(b_c,b_w,b_0)={1\over2\epsilon^{d-2}}(g_c^2,g_w^2,g_Y^2/36),\quad
\kappa={k_h\over v},\quad
K_i={\kappa\over2b_i}>0\ (i=c,w,0).
\tag{3}
$$

K_i是磁面势在原群角中的Hessian，不是吸收耦合后的字段系数。特别是U(1)采用g_Y/6，不能套用SU(2)半角系数。把旧一圈MS参数代入(3)是**经典匹配处方**，未计算MS到格点的量子匹配。规则格图、d、v、c与尺度族仍是输入。

## 3. 完整有限图中的中性约化

取常真空X=h★n、n=(0,1)，s=s★，平坦链路。D的行是源+1、靶-1，面边关联矩阵记作mathcal C，满足mathcal C D=0。中性角使用X=h★exp(i t3 phi)n、W=exp(i t3 theta)，只作真空附近经典二次展开。

定义a=theta-6alpha-Dphi、p=pi_w、e=pi_0+6pi_w，则正则一形式和Gauss为

$$
p_\phi\,d\phi+\pi_w\,d\theta+\pi_0\,d\alpha
=p\,da+e\,d\alpha+(p_\phi+D^T\pi_w)d\phi,
\quad G_w=p_\phi+D^T\pi_w=0,
\quad G_0=-6p_\phi+D^T\pi_0=0
\ \Longrightarrow\ D^Te=0.
\tag{4}
$$

因此消去第一约束后，alpha仍按alpha+Dgamma等价。全非线性情况下，非零Higgs在SU(2)×U(1)的稳定子是alpha↦(diag(exp(3i alpha),exp(-3i alpha)),exp(i alpha))；亦不能把旧SU(2)的边商张量积原封搬入。

令mu=k_h h★²/4。二次约化动能为

$$
T_{\rm red}={2\over vh_\star^2}\|D^Tp\|^2+b_w\|p\|^2+b_0\|e-6p\|^2,
\qquad D^Te=0,
\tag{5}
$$

二次势为

$$
V_{\rm red}={\mu\over2}\|a\|^2
+{K_w\over2}\|\mathcal C a+6\mathcal C\alpha\|^2
+{K_0\over2}\|\mathcal C\alpha\|^2.
\tag{6}
$$

证明只用(4)的正则变换及mathcal C D=0。非零梯度模中e纵分量被约束为零，故只有一个中性有质量纵模。横向规范化字段中，质量矩阵及零方向为

$$
M^2_{wY}={c^2h_\star^2\over4}
\begin{pmatrix}g_w^2&-g_wg_Y\\-g_wg_Y&g_Y^2\end{pmatrix},
\quad A_\gamma\ \parallel\ (g_Y,g_w),
\quad m_Z^2={c^2h_\star^2\over4}(g_w^2+g_Y^2),
\quad m_W^2={c^2h_\star^2\over4}g_w^2.
\tag{7}
$$

这重现547已有结构，不作为新质量定理。式(3)与原标量匹配一起给每个非零格图动量

$$
\lambda=4\sum_{j=1}^{d}\sin^2(\epsilon p_j/2):\qquad
\omega_{Z,L/T}^2=m_Z^2+\kappa\lambda,\quad
\omega_{W,L/T}^2=m_W^2+\kappa\lambda,\quad
\omega_{\gamma,T}^2=\kappa\lambda,\quad
\omega_{h,s}^2=m_\pm^2+\kappa\lambda.
\tag{8}
$$

在d≥2的周期格图中，电弱加双径向共有3d+(d-1)+2=4d+1个非零动量局部模式；d=3为13，d=2为9。颜色若另计则在线性层有8(d-1)个横向模，仍不是手征物质。p=0存在d个电磁谐和零频方向，不能按有正规Gaussian真空的非紧致零频振子处理。全局残余群及总荷的高阶约束未被线性计数签收。

## 4. 全局商选择暴露的磁拼接障碍

商群中的每条链路允许独立改lift：乘以任意zeta次幂。节点规范不变不等于这种lift不变。568的SU(2)基本迹与独立圆磁项若直接搬到商群，则同一个单位面holonomy的两个lift给

$$
V_{\rm cover}=\beta_w(1-\operatorname{ReTr}W/2)+\beta_0(1-\operatorname{Re}z):
\quad V_{\rm cover}(I)=0,\qquad
V_{\rm cover}(\zeta)=2\beta_w+\beta_0/2>0.
\tag{9}
$$

因此它不是G上的函数。这只排除“最大商群＋不变更原基本磁迹”的接法；覆盖群、其它中心子商、其它磁作用都未被排除。

一种简单的中心不敏感替换是

$$
V_{\rm blind}=a_c(9-|\operatorname{Tr}C|^2)
+a_w(4-|\operatorname{Tr}W|^2)+a_0(1-\operatorname{Re}z^6).
\tag{10}
$$

各项非负，单位元Hessian分别为3a_c、2a_w、36a_0；故能配任意正K。但其共同零集合在覆盖群为Z3×Z2×Z6的36个中心点；除对角Z6后仍剩6个逐面零类。单凭二次主部相同不能把(10)当作全局等价替换。

## 5. 正权字符补全定理

对任意w_q,w_l>0，采用531实际表示的正权字符和

$$
\chi_R(C,W,z)=w_q[\operatorname{Tr}C\operatorname{Tr}W\,z
+\operatorname{Tr}\bar C(z^{-4}+z^2)]
+w_l[\operatorname{Tr}W\,z^{-3}+z^6+1],
\qquad D_R=12w_q+4w_l.
\tag{11}
$$

实数权不是新的非整数粒子重数。这是对真实酉表示字符的正组合。各表示均对zeta不变，且每项dim-ReTr非负。单位元处三块指标为

$$
I_c=2w_q,\qquad I_w={3w_q+w_l\over2},\qquad I_Q=66w_q+54w_l;
\qquad \operatorname{Hess}(D_R-\operatorname{Re}\chi_R)
=\operatorname{diag}(I_c I_8,I_w I_3,I_Q).
\tag{12}
$$

证明由各表示生成元迹求和；混合块含非Abelian生成元的零迹，因此消失。给定任意正(K_c,K_w,K_0)，取

$$
0<\delta<\min\{K_c/I_c,K_w/I_w,K_0/I_Q\},\quad
a_c={K_c-\delta I_c\over3},\quad
a_w={K_w-\delta I_w\over2},\quad
a_0={K_0-\delta I_Q\over36}.
\tag{13}
$$

定义

$$
V_f=\delta(D_R-\operatorname{Re}\chi_R)
+a_c(9-|\operatorname{Tr}C|^2)+a_w(4-|\operatorname{Tr}W|^2)
+a_0(1-\operatorname{Re}z^6).
\tag{14}
$$

**命题569-A。** 式(14)是G上的光滑、非负、共轭不变面势，其单位元Hessian恰为所给(K_c,K_w,K_0)，且每个商群面只有单位类取零。

证明：严格正系数及上述迹界给非负和下降；(12)—(13)直接给Hessian。若首项为零，每个实际酉表示块的ReTr均须达到维数，因其所有特征值在单位圆上，该块只能是恒等。e^c迫使z^6=1；L迫使W=z^3I2；d^c迫使C=z²I3，剩余块一致。因此

$$
V_f(C,W,z)=0
\iff (C,W,z)\in\langle\zeta\rangle
\iff [C,W,z]=1_G.
\tag{15}
$$

逐面单位不排除周期空间的非平凡平坦全局holonomy，更不证明全量子理论唯一真空。delta及正补项属于新增磁作用输入；**不能称为531的“只用同一个加权迹”条件T在任意尺度的实现**。553旧限制完整保留。

## 6. 有限量子载体存在，但不等于连续量子连接

有限图上可以取

$$
\mathcal H=L^2((\mathbb R^4\times\mathbb R)^V\times G^E),
\qquad \mathcal H_{\rm phys}=\mathcal H^{G^V},
\qquad H_{\rm phys}=H_{\rm Friedrichs}\big|_{\mathcal H_{\rm phys}}.
\tag{16}
$$

在紧支撑光滑物质函数和光滑群函数的稠密域上，式(2)的动能二次型非负，物质四次势及边势非负；每个有限面势(14)有界。各项构成标准非负Schrödinger型闭形式的闭包，给Friedrichs自伴算子。紧群平均与演化相容，故Gauss不变子空间约化该算子。群商Casimir及Higgs作用良定义：zeta在z³W中作用为+I。

此有限存在论证不要求预置非正规“经典真空量子态”。但它也不证明存在连续量子相、规范不变粒子极点、受控RG、手征费米子或动态时空。场景中的面胞、图和全局群仍由建模输入指定。

在光滑指定字段的经典尺度层面，(3)、(12)—(14)给与568同类的曲率平方二次主项；不同delta保留相同主项却改变有限holonomy势。本轮没有新增网格扫描或把这个Taylor事实当量子收敛。

## 7. 八组可复算检查及数值

沿用544同一第三个状态mu=173.34 GeV和u=18.0918442410，旧运行给(g_c²,g_w²,g_Y²)=(1.37210090779,0.420371362122,0.128720796917)。取c=epsilon=v=k_h=k_s=1，h★²=0.442814340447；物质势与径向质量不重选。这些间隙使用原共同质量单位，不能直接当作实测GeV²。

|检查|结果|
|---|---|
|同尺度荷与系数|b=(0.686050454,0.210185681,0.001787789)，K=(0.728809371,2.378849013,279.675086406)；m_W²=0.046536617，m_Z²=0.060786471|
|完整图辛形式、能量、Gauss及方程|最大残差2.49e-14|
|完整中性物理图谱|二维28坐标、三维136坐标；与保留约束的谱最大差3.84e-11|
|非线性Higgs边势＋真实磁补全的Hessian|12变量、步长2e-4，最大差1.28e-4；非解析证明的替代品|
|中心lift及零类|错误原磁项在同一商单位给0或158.8683353；简单伴随替换6个零类；正权补全只1个|
|全12生成元Hessian|8色、3弱、1圆；最大有限差分误差7.19e-5，混合块数值0|
|非线性节点规范变换／逐边独立lift|残差分别7.11e-15、2.85e-14；正性和全零集合由解析证明保证|
|剩余形状自由度|三种delta均精确同K；同一有限holonomy势相差0.192660618，故完整作用仍不唯一|

w_q=1、w_l=2.32321410578取自旧历史匹配，只作本正函数的权。选delta为允许上界的一半，得到delta=0.182202342831，(a_c,a_w,a_0)=(0.121468228554,0.946948986281,6.799772192180)。当前低尺度旧单迹关系缺差为1.60395118871，补项没有把它抹成零。

## 8. 回填、物理解释与下一项

本轮连接C02/C14的Gauss、C15/C16的全局表示、C05/C06/C10的共同传播及C17/C20的同尺度物质参数。它没有证明这些条件彼此等价，也没有总体减少输入数：具体不兼容接法被排除，一族可用的正补全被构造，磁形状自由度和全局选择留下。

接下来优先核**共同物质、参考与动态几何的联合约束数据**：548—552的经典曲率／参考候选，与此处新增规范应力、同一h＋s来源是否能放在一组满足完整约束的数据中。先回查旧参考和Einstein约束构造，不能仅通过“把导数换成协变导数”新增一轮，也不把预置引力作用后的相容性称为引力生成。具体子任务须由回查后真正未闭合的接口确定。

不继续扫描本轮delta、局部斜率或旧中性质量；它们现已足以说明构造和限制。旧557—566数值依赖原H，不能因同名物质就直接签收加项后的装置性能。历史反例、冻结报告和统一目标均保持；本轮完成不代表三维、GR、完整标准模型或统一完成。


**复算与审查：** [代码](joint_quotient_gauge_completion.py)、[结果](joint_quotient_gauge_completion_results.json)、[核验](research_round_569_checks.json)、[条件账](unified_physics_condition_ledger_569.md)、[独立预审](round569_drafts/electroweak_pre_review.txt)、[独立终审](round569_drafts/final_review.txt)。原草稿和入口试算均保留。
