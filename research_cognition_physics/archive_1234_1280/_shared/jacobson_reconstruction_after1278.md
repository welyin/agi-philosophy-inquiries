# 1278后：现行CT／HG角色到有效Einstein方程的成熟定理接续

2026-10-17。非编号理论接续；新增科学轮次、实验和核验组均为0。依据[正式采用记录](ct_hg_adoption_after1278.md)及[机器记录](ct_hg_adoption_record_after1278.json)，CT、HG1、HG2、HG3现已采用。旧[1233重建接口](../../archive_1219_1233/_shared/gravity_reconstruction_interface_after1233.md)登记的条件定理可据此直接接续；旧稿中的“未采用”保留历史含义。本稿不修改1278或其它冻结科学证据，不宣布共同实现已经完成。

**结论。** 在下述共同有效域内，HG1—3与标准几何恒等式给出Einstein主阶方程；CT在其原低阶正则类中给出相同动力学系数的独立一致性检查。面积熵密度与宇宙学常数的数值仍未确定，有限过程保留显式余项。这是成熟重建的当前公理接续，不是本项目重新发现Jacobson定理。

## 1. 前提归属与约定

|前提|身份及范围|
|---|---|
|CT、HG1—3|现行已采用角色；按采用记录的共同域及防循环条款使用，不从旧TR或普通记录一致性偷偷补出。|
|共同几何和物质|明确有效域输入：连通四维光滑Lorentz域，同一实际钟／位置／信号字典，同一对称完整非几何应力$T_{ab}$，以及$\nabla^aT_{ab}=0$。来源、控制、记录、相互作用及功／荷分账均属于同一过程。|
|局部平衡边界族|HG的适用域：每点、每个所需零方向都有合格局部因果边界，参考及观察者归一相容；不是只检查一个面片。|
|受控有效描述|正则性允许局部展开，源涨落、近似boost、有限分辨率及非平衡项在主阶可忽略或另有预算。量子应用中取有效平均来源须另有资格，不由公理命名自动取得。|
|Raychaudhuri、Bianchi及零锥代数引理|数学恒等式／引理，不新增认知公理。|

主文取$c=\hbar=k_B=1$，度规签名为$(-,+,+,+)$，曲率号约定以式(3)的Raychaudhuri方程为准。$s>0$为HG2在该连接域内事前固定的面积熵密度；$\kappa_\chi>0$是观察者boost归一，$\kappa_g$是引力耦合，二者不同。

独立几何确定终端二维小面元$\Sigma_0$及其过去侧无扭率、无焦散的零测地线束。未来指向的仿射切向量为$k^a$，仿射参数取$\lambda\in[-\ell,0]$，在终端满足

$$
\theta(0)=0,\qquad \sigma_{ab}(0)=0,\qquad
\chi^a=-\kappa_\chi\lambda k^a,\qquad
 d\Sigma^a=k^a\,d\lambda\,dA .
\tag{1}
$$

局部boost等式按主阶使用；有限偏差另列。仿射参数不等于固有时。原实际钟须先给观察者归一；正比例重标$k,\lambda$不改变加权积分。缩放$\chi$会共同缩放温度和热量，不改变二者之比。

## 2. HG给出的共同热学表达

HG1、HG2给同一边界、同一参考及观察者的

$$
T_\chi=\frac{\kappa_\chi}{2\pi},\qquad
\frac{\delta Q_\chi}{T_\chi}
=2\pi\int_{\mathcal H}(-\lambda)T_{ab}k^ak^b\,d\lambda\,dA,
\qquad \delta S_{\rm geom}=s\,\delta A .
\tag{2}
$$

几何面积由钟尺／传播字典先行定义，$s$不逐过程拟合。HG2认证的是几何量的热学角色，不要求它等于裸区域von Neumann熵或记录熵。HG1中的热量是该观察者的完整穿界能流，不要求入射激发本身先为热态，也不要求先实现两热库接触；工作、荷及参考变化仍须按同一合同处理。

HG3要求式(2)在局部平衡邻域满足可逆主阶Clausius关系。这正是新增的桥接要求；一般能量守恒、Gibbs形式或第二定律不单独推出它。以上读法与[Jacobson 1995原文式(1)—(2)及局部平衡说明](https://arxiv.org/html/gr-qc/9504004v2)一致。本稿只承接该主阶路线，不把它扩成任意有限过程精确可逆。

## 3. 聚焦推出全部零方向曲率关系

复用[312的重建接口](../../archive_301_341/research_note_312.md)和[316的有限几何工具](../../archive_301_341/research_note_316.md)。四维无扭率零束满足

$$
\frac{d\theta}{d\lambda}
=-\frac12\theta^2-\sigma_{ab}\sigma^{ab}-R_{ab}k^ak^b,
\qquad \frac{d(dA)}{d\lambda}=\theta\,dA .
\tag{3}
$$

在终点$p$固定方向$k$。由式(1)和正则性，沿该束

$$
\theta(\lambda)=-\lambda R_{kk}(p)+O(\lambda^2),
\qquad
\delta A=\frac{\ell^2}{2}A_0R_{kk}(p)+o(\ell^2A_0).
\tag{4}
$$

这里同时取横向面元足够局部，使跨面元的系数变化进入余项；$A_0$是终端面积。类似地，式(2)给

$$
\frac{\delta Q_\chi}{T_\chi}
=2\pi\frac{\ell^2}{2}A_0T_{kk}(p)+o(\ell^2A_0).
\tag{5}
$$

HG3的主阶关系意味着这两式在除以$\ell^2A_0/2$后余项消失，因而

$$
\boxed{R_{ab}k^ak^b=\kappa_gT_{ab}k^ak^b,
\qquad \kappa_g=\frac{2\pi}{s}.}
\tag{6}
$$

它对有效域每点的全部合格零方向成立。此处的局部极限是有效主阶系数的提取，不是授予可执行无限加速度、无限精度或无限次实验。若只取得有限尺度／有限方向预算，则保留相应残差，不能宣称式(6)精确普遍成立。

## 4. 零锥代数与协变守恒闭合

令$D_{ab}=R_{ab}-\kappa_gT_{ab}$。在一点的正交标架内，对所有$|\mathbf n|=1$取零向量$k=(1,\mathbf n)$。式(6)成为

$$
D_{00}+2D_{0i}n_i+D_{ij}n_in_j=0 .
$$

分别取$\mathbf n$与$-\mathbf n$得$D_{0i}=0$；余下空间二次型在单位球面为常数$-D_{00}$，故$D_{ij}=-D_{00}\delta_{ij}$。因此存在标量$f$使

$$
R_{ab}-\kappa_gT_{ab}=f g_{ab}.
\tag{7}
$$

$s$在本域固定，所以$\kappa_g$也是常数。由完整来源守恒和缩并Bianchi恒等式，

$$
\nabla_bf=\frac12\nabla_bR,
\qquad f=\frac12R-\Lambda ,
$$

其中$\Lambda$在每个连通域内为常数。于是

$$
\boxed{G_{ab}+\Lambda g_{ab}=\frac{2\pi}{s}T_{ab},
\qquad G_N=\frac1{4s}.}
\tag{8}
$$

这是有效Einstein方程；$G_{ab}=R_{ab}-Rg_{ab}/2$，而$G_N$是Newton常数。恢复$\hbar$、仍取$c=k_B=1$时，右侧系数为$2\pi/(\hbar s)$，$G_N=1/(4\hbar s)$。它不确定$s$、$\Lambda$、物质参数或初边值，也不把Weyl曲率的初始／辐射自由说成重建失败。

若$s$随位置变化，式(7)的散度多出$T_{ab}\nabla^a\kappa_g$，上述闭合不再原样成立。非最小曲率耦合、高阶熵密度和非平衡项必须使用对应的完整来源与熵修正；不能把它们静默塞入常数面积主项。

## 5. 与已采用CT的系数核对

CT本身不要求所有材料都属于同一种Hamiltonian类别。在其原声明类别中，另有独立正空间度规$h_{ij}$与共轭$\pi^{ij}$、局域二次动量／至二阶空间导数主部、最小物质分离和合格五片菜单。记该类系数为$a_g,b_g,\lambda_g$，避免与面积或仿射参数混淆：

$$
\mathcal H_g=\frac{a_g}{\sqrt h}
(\pi^{ij}\pi_{ij}-\lambda_g\pi^2)
-b_g\sqrt h\,{}^{(3)}R+2\Lambda_H\sqrt h .
\tag{9}
$$

[1188](../../archive_1161_1200/1188/proof.md)、[1235](../1235/proof.md)及[1245完整材料接续](../1245/proof.md)可直接复用：CT的独立切片资格与动力身份要求给$\lambda_g=1/2$、$a_gb_g=1$；完整场方程还保原基片约束、全部法向／动量约束及其演化保持。**这些在1245原材料与类别内已有完整条件证明，并非仍缺动量约束**：1245 §§5—6给一般三维约束闭合、合格初片及平滑窗口内保持；1246—1247的内部程序、控制和记忆扩展已按同一最小材料资格承接，见[原联合定理§4](gr_joint_stage_theorem_draft_after1251.md)。新候选只需按其实际变化检查能否继承，不重开这份旧证明。对应作用为

$$
I_g=\frac1{a_g}\int\sqrt{-g}(R-2\Lambda)\,d^4x,
\qquad 8\pi G_N=\frac{a_g}{2}.
$$

所以在同一物理来源、同一正则归一与HG共同有效域内，必须有

$$
\boxed{a_g=2\kappa_g=\frac{4\pi}{s},\quad
b_g=\frac1{2\kappa_g}=\frac{s}{4\pi},\quad
\Lambda_H=\frac{\Lambda}{2\kappa_g}.}
\tag{10}
$$

相应Wald系数为$s_W=4\pi/a_g=s$。这里的$\Lambda$就是物理宇宙学常数，$\Lambda=a_g\Lambda_H$。这是两路线的共同归一核对，不用$T=G/\kappa_g$反定义来源，也不以系数相等倒定义独立几何熵。CT没有选出$a_g$的绝对数值，HG没有选出$s$的绝对数值；式(10)只使它们在同一对象中不再独立。

若先采用CT选后场方程，再推出该对象的Clausius一致性，这属于正向热学后果；不能把同一结论回代为第二次独立生成Einstein的证据。HG路线式(2)—(8)没有调用CT的系数定理，两条逻辑依赖因此可清楚分开。

## 6. 有限窗口：独立修正与残差

不重新证明或计算316／[1250](../1250/proof.md)的有限结果。对同一过去片、固定生成线标签和切口，令$dV_{\mathcal H}=d\lambda\,dA$，定义

$$
E_{\rm geom}=\ell\int_{\Sigma_0}\theta_0\,dA
+\int_{\mathcal H}\left[(2\lambda+\ell)R_{kk}
+(\lambda+\ell)\left(\sigma^2-\frac{\theta^2}{2}\right)\right]dV_{\mathcal H} .
\tag{11}
$$

由式(3)和两次分部积分精确得

$$
\Delta A-E_{\rm geom}
=\int_{\mathcal H}(-\lambda)R_{kk}\,dV_{\mathcal H}.
$$

因此在精确式(2)字典下，有限检验量为

$$
\boxed{\mathcal R_{\mathcal H}
:=s(\Delta A-E_{\rm geom})-\frac{\delta Q_\chi}{T_\chi}
=s\int_{\mathcal H}(-\lambda)
(R_{kk}-\kappa_gT_{kk})\,dV_{\mathcal H}.}
\tag{12}
$$

$E_{\rm geom}$完全由独立几何计算，不能事后定义为面积与热量之差。膨胀和剪切平方的符号不同；末端不精确平衡时，式(11)首项不可删除。有限面积权与boost权不同，也不能混同。有效Einstein关系精确成立时$\mathcal R_{\mathcal H}=0$，但一般仍有$s\Delta A-\delta Q/T=sE_{\rm geom}$，不是任意有限过程的裸Clausius等式。

实际检验还须加入参考／温标、近似boost、完整来源、几何定位、熵主项修正及读数的独立误差。若这些给$|\mathcal R_{\mathcal H}|\le\epsilon_{\mathcal H}$，所得只是该权重的有限曲率—来源残差界；从有限积分到点态／全方向方程还需原[1243](../1243/proof.md)一类覆盖与稳定性前提。不得把极限定理当成有限精度下零误差的实验证书。

## 7. 已完成的接续与共同实现边界

本稿完成“已采用角色＋共同有效域前提$\Rightarrow$有效Einstein主阶方程”的成熟定理接续及CT系数核对。其数学结果、假设和剩余常数均已明确；它不要求重造一般热设备，也不重新开放四律证明。

共同实现仍须核对同一边界、参考、实际钟、独立面积和完整来源。1250已有选后有限余额；1270—1272已有各自真实量子过程与来源预算；1278已有新候选内裸能源交换与独立钟／传播响应。它们的域、准备、边界和误差不能无证拼接。特别是1250大域面片不自动成为1270小块的局部仪器，有限均匀来源的数学积分不自动给局部量子参考。

下一步在既有共同模型上作有限资格与预测核验；不再把CT／HG是否采用列为待决问题，也不因这些角色已采用而宣布实际共同实现、强关系唯一生成猜想或全部目标完成。本文不改变既有科学计数。
