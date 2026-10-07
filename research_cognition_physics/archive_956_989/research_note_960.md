# 960：原生材料的电磁色散接口与有限有效域

日期：2026-10-07。承接[959](research_note_959.md)及[计算前决定](960/drafts/material_propagation_decision.md)。代码：[material_dispersion_bridge.py](960/material_dispersion_bridge.py)；[结果](960/material_dispersion_bridge_results.json)；[交付核验](960/research_round_960_checks.json)。

## 1. 本轮回答的是采用价值

**956—958的电子材料有一条可复用的成熟电磁接口：它的真实电荷跃迁谱，给出的非延迟二阶能量恰为958条件记录相位的领先系数；引入同一Maxwell传播核后，频率响应、延迟能量和距离力共用这些系数。** 这支持继续采用原生材料，而无需另发明一种记录传播作用。本轮没有建立完整实时QED过程，也没有证明完整共同模型存在。整体目标未完成。

这是已有London／Casimir–Polder理论在本项目实际材料上的匹配，不是新的物理相互作用发现。959精确静电父表示、952同核来源变分、941与953内能／协变接口均保留，不重复证明；372、修订382—384、386或425、522—523的空间接口保持原范围。

**停止决定：** 不接着优化色散力、轨道模型、延迟积分、转子或光子器件；这个候选接口有采用价值，已经足够改变取舍。完整父模型的接合仍应从共同有效域整体判断，不能变成“先修完这个双二聚体的一切问题”。

## 2. 输入和认识论层次

|层次|本轮内容|
|---|---|
|认知动机|同一材料的内部差异能够经同一物理作用形成记录，并参与其来源及反作用；复用H2/H3，不新增认知公理|
|已有解析输入|956材料谱和关系编码、958静电条件相位、959保源静电父表示|
|额外物理输入|3+1真空Maxwell传播、固定零温定态、平行电偶极且与分离轴垂直、长度／电荷／传播速率的一次共同匹配|
|本轮推导|实际材料极化率、领先静电匹配、延迟核的能量与力界、高频加权尾界|
|数值核验|六维及36维真实谱、积分双精度网格对照、独立来源差分、旧几何不匹配量|
|未领取|3+1的认知必然性、真实材料全频响应、原958实时误差运输、完整SM／引力／六协议共同实现|

这里没有要求微观连续，也没有设定物理最小长度。连续Maxwell核是选取的有效描述。材料空间方向与内部自旋SU(2)不是同一件事。

## 3. 从原材料得到真实响应

取956的两轨道、两电子六维材料，保留其CAR符号约定。记

$$
H_m=-v\,{\rm hop}+UD,\quad Q=\frac{n_1-n_2}{2},\quad
\Omega=\sqrt{U^2+16v^2},\quad J=\frac{\Omega-U}{2},\quad
d=\frac{1-U/\Omega}{2},\quad \Delta=U+J.
\tag{1}
$$

单态基态能量为$-J$，$Q|g\rangle=\sqrt d\,|d_-\rangle$，而$H_m|d_-\rangle=U|d_-\rangle$。三个低能三重态均被$Q$湮灭。因此，对外部广义电场耦合$-f_{\rm ext}Q$的虚频线性响应，在低部门是

$$
\alpha_Q(i\xi)=
\sum_{n\ne g}\frac{2(E_n-E_g)|\langle n|Q|g\rangle|^2}
{(E_n-E_g)^2+\xi^2}
=\frac{2\Delta d}{\Delta^2+\xi^2},
\qquad \widehat\alpha_Q=\alpha_Q P_s.
\tag{2}
$$

这不是用$PQP=0$误判材料没有电响应：虚跃迁经过高电荷态。三重态无电偶极响应只针对这个两轨道有效模型；不能推广为真实三重态物质的普遍性质。物理偶极为$p=eaQ\,\mathbf u$，物理极化率相应乘$(ea)^2\mathbf u\mathbf u^{\mathsf T}$。完整轨道匹配可能贡献额外响应，不能由本六维谱自动排除。

将两份材料用958的$\kappa Q_AQ_B$相连，普通定态二阶微扰和虚频积分独立给出

$$
\chi_2=\frac{\kappa^2d^2}{2\Delta}
=\frac{\kappa^2}{2\pi}\int_0^\infty
\alpha_Q(i\xi)^2\,d\xi,\qquad
H_{\rm cond}^{(2)}=-\chi_2P_s^AP_s^B.
\tag{3}
$$

所以，原记录相位的领先系数就是此材料的非延迟色散能。计算在固定自旋扇区最低态及真空定态意义下成立；电偶极保持自旋扇区，未假定未知相干输入已经处于完整场—物质的联合定态。

原$U=1,v=.1,\kappa=.2$下，$\chi_2=2.4629285339836\times10^{-5}$，六维极化率公式与谱求和最大差$9.72\times10^{-17}$，36维微扰与积分一致。958精确静电$\chi=2.4866347759903\times10^{-5}$仍保留；两者不相等，旧读取时刻乘其差约$0.0299503$，不能把“领先系数相同”写成“原长时间通道已相同”。

## 4. 一份传播核同时给能量和力

采用零温弱电偶极多重散射的首项。所用框架见[Milton等，2015，式(2.16)](https://arxiv.org/abs/1502.06129)；远区各向异性系数可与[Milton等，2012，式(2.1)](https://arxiv.org/abs/1111.4224)交叉核对。这里只引入已知Maxwell描述，不宣称由认知公理导出。

取$\mathbf u_A=\mathbf u_B\perp\widehat{\mathbf R}$，令$C=e^2/(4\pi\epsilon_0)$、$\kappa_{\rm dip}=Ca^2/R^3$。Maxwell虚频电偶极核在此方向为$\kappa_{\rm dip}f(x)$，其中$f(x)=e^{-x}(1+x+x^2)$，$x=\xi R/c$。它也可由Yukawa Green函数的Hessian与频率项直接取得。于是

$$
\chi_{\rm EM}^{(2)}(R)
=\frac{\kappa_{\rm dip}^2}{2\pi}
\int_0^\infty\alpha_Q(i\xi)^2 f(\xi R/c)^2\,d\xi
=\chi_{\rm L}(R)\,F(r),
\quad
F(r)=\frac4\pi\int_0^\infty\frac{f(ry)^2}{(1+y^2)^2}\,dy,
\quad r=\frac{R\Delta}{c},
\quad\chi_{\rm L}=\frac{\kappa_{\rm dip}^2d^2}{2\Delta}.
\tag{4}
$$

本式是电偶极、弱跨体耦合的**领先定态**式。式(3)中的$\kappa$只有经共同几何匹配后才能换为$\kappa_{\rm dip}$，不能只在数值上相等就宣布原几何已保留。

保持$U,v,C,a,c$固定，同一相互作用能$E_{\rm int}=-\chi_{\rm EM}^{(2)}$给沿增大$R$方向的力：

$$
\mathcal F_R=-\partial_R E_{\rm int}
=\frac{\chi_{\rm L}}R[-6F(r)+rF'(r)],
\qquad
\frac{\mathcal F_R}{\mathcal F_{\rm L}}=F(r)-\frac{rF'(r)}6,
\qquad \mathcal F_{\rm L}=-\frac{6\chi_{\rm L}}R.
\tag{5}
$$

只微分$R^{-6}$会漏掉传播核的同阶来源。式(5)连接条件能量与保留阶机械力；因只依赖相对位置，两物体的静态力相反。这不是完整动态能源平衡或Einstein应力张量的证书。

远区由$\int_0^\infty e^{-2x}(1+x+x^2)^2dx=13/4$得到$F(r)\sim13/(\pi r)$，故从$R^{-6}$过渡到$R^{-7}$。该系数对应本各向异性方向，不能错用各向同性的23。$r=.01$时$F=1.000092682772$；$r=1$时$F=1.074685219129$，所以不应预设“延迟只会单调减小吸引能”。这些是定态比较，不是超光速通信或时空维数证据。

## 5. 有效范围所需的界，及不必追求的全UV

对$x\ge0$，$0\le f(x)\le3/e$，$f'(x)=x(1-x)e^{-x}$且$|f'(x)|\le x$。所以$|f(x)-1|\le x^2/2$。结合$\int_0^\infty y^2(1+y^2)^{-2}dy=\pi/4$，直接得到

$$
|F(r)-1|\le B_0r^2,\quad |rF'(r)|\le B_1r^2,\quad
\left|\frac{\mathcal F_R}{\mathcal F_{\rm L}}-1\right|
\le\left(B_0+\frac{B_1}{6}\right)r^2,\qquad
B_0=\frac{1+3/e}{2},\quad B_1=\frac6e.
\tag{6}
$$

这是对所选响应模型的解析界，不靠截断高频积分证明。在$r\le.01$，领先能量与领先力的非延迟相对误差分别不超过$1.05182\times10^{-4}$和$1.41970\times10^{-4}$。它们不包含高阶跨体耦合、遗漏轨道、有限尺寸、制备或实时辐射误差。

若响应只在$\xi\le\Lambda$得到匹配，高频不必逐频精确求解。以$\Lambda=Y\Delta$，在尾部明示假设$|\alpha_Q(i\xi)|\le2\Delta d/\xi^2$，则相对于$\chi_{\rm L}$和$|\mathcal F_{\rm L}|$有

$$
\epsilon_{\rm tail,E}\le\frac{12}{\pi e^2Y^3},\qquad
\epsilon_{\rm tail,F}\le
\frac{12}{\pi e^2Y^3}+\frac{4r^2}{\pi eY}.
\tag{7}
$$

证明：能量尾部用$f^2\le9/e^2$与$\int_Y^\infty y^{-4}dy$；来源项用$|2f(x)xf'(x)|\le(6/e)x^2$与$\int_Y^\infty y^{-2}dy$，再按式(5)除6。这里还要求尾部响应在改变$R$时无额外未计来源，符合本固定孤立材料比较。

在$Y=100,r\le.01$，单份模型尾界分别$5.16943\times10^{-7}$和$9.85342\times10^{-7}$。若比较两份低频已匹配且都满足此界的完成，尾部差用两界之和。**此尾部假设不是实际SM材料已满足的实验证明；未知高轨道也可能改变低频极化率，仍要做一次低能匹配。** 本结论是需要什么有限信息的充分条件，不能冒充完整物理模型已认证。

积分的256／512点对照和来源独立差分仅为复算校验，不冒充严格积分误差证明；式(6)—(7)才是本轮使用的解析控制。代码另登记远区$x>40$的指数尾界。

## 6. 一个影响有限预测的真实边界

958可选的四点几何使用

$$
\kappa_{\rm site}(R)=2C\!\left(\frac1R-\frac1{\sqrt{R^2+a^2}}\right),
\qquad
1-\frac{3a^2}{4R^2}\le
\frac{\kappa_{\rm site}}{\kappa_{\rm dip}}\le1
\quad(0\le a/R\le1).
\tag{8}
$$

由$(1+z)^{-1/2}$的Taylor余项可证。旧$R=3,a=1,C=5.84604989415$给$\kappa_{\rm site}=.2$、$\kappa_{\rm dip}=.21652036645$；其领先能量比$1.17202672719$，相差约17.20%。这影响有限预测，因此不能直接继承；但旧四点几何本来就是额外选择，**修到任意多极精度不成为整个纲领的必经门槛**。可以选择真正满足电偶极近似的声明域，或暂不采用该几何对应的预测。

## 7. 本轮的全局取舍

|事项|判断与处置|
|---|---|
|原生材料是否有物理接口价值|有：同一谱接上成熟电磁响应，保留作为优先接口|
|材料记录与场来源是否仍可任意独立选|在本领先定态接口不可；必须共享$\Delta,d,\kappa_{\rm dip}$及其来源|
|原958的长时间全未知输入界是否已进完整场论|没有；定态、匹配阶及准备不同，不自动运输|
|所有高能／精细外推是否必需|不必；本有限观测用匹配域与加权尾界即可描述|
|完整共同模型是否完成|没有；全局恢复表所有部门继续保留，不以本色散接口代替|
|是否继续修这个模型|停止；不研究更多偶极排列、高阶色散、仪器寿命或通量网络|

[961入口](961/drafts/STATUS.md)继续选择完整父有效描述的同一保留阶及实际共同菜单。允许把成熟物理、有限匹配参数和准备权限明确作为输入，不要求重造QED／GR或制造每一个仪器；同时，不把“完整共同实现存在”本身作为已满足前提。本轮减少的是一个独立传播交互的需要，不是减少3+1、Maxwell或Einstein作用等物理输入。

## 8. 复算

使用现有Python与NumPy，无下载、无图像检验：

~~~powershell
$env:OPENBLAS_NUM_THREADS='1'
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics\archive_956_989\960\material_dispersion_bridge.py'
& 'C:\Users\admin\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -B -X utf8 'research_cognition_physics\archive_956_989\960\verify_round960.py'
~~~

默认只读核对已保存结果。正式960／累计3745；统计不代表全局完成百分比。
