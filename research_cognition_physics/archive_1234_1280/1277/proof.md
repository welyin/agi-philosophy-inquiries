# 1277：原自治等待的参考准备界

2026-10-17。按[准入](admission.md)，本轮只研究[1237原完整过程](../1237/proof.md)及[1274的精确边缘通道](../1274/proof.md)。新结果是原等待菜单的全时间有限界；不新增门、材料、实际调度器或公理。验收状态以独立审阅与主线发布为准。

## 1. 同一父过程及三种不同任务

原边界空间为$\mathcal H_B=\mathcal H_4\oplus\mathcal H_6$，维数分别为2和5；$P_4,P_6$是对应投影。原面积、边界能源和固定参考为
$$
\widehat A=4aP_4+6aP_6,\qquad H_B=4eP_4+6eP_6,\qquad
r_*:=r_{B,\beta_*}=\frac{P_4}{4}+\frac{P_6}{10},\qquad
\beta_*=\frac{\log(5/2)}{2e}.
\tag{1}
$$
这里$a,e,g>0$及原正谱条件$4e>\sqrt2g$保持不变。$r_*$是1274已证明的**完整父Hamiltonian Gibbs参考的边缘**；本轮不把该数学身份当作实际准备许可。

工作库存、720维相干旗标、三站程序及任意未知输入伙伴$R$均按原合同保留。只为描述边缘取偏迹，不实际擦除完整输出。设$W:\mathbb C^2\to\mathbb C^5$为原等距，$\sigma_4=I_2/2,\sigma_6=I_5/5$。由1274式(5)—(6)，对全部$X\in M_2$，
$$
\Phi_u(X)=u^2X\ \oplus\left[2u(1-u)WXW^\dagger+(1-u)^2\operatorname{Tr}(X)\sigma_6\right],
\quad u=\frac{1+\cos(\sqrt2gt)}2\in[0,1].
\tag{2}
$$
任意$u\in[0,1]$都在原$t\in[0,t_*]$内达到，$t_*=\pi/(\sqrt2g)$；原全时间轨道因此由同一区间覆盖。式(2)张量恒等作用于$R$，不是中途选择程序分支。

本轮分开三个问题：固定等待时间能否准备$r_*$；若**另行声明**输入独立的随机等待数学合同，已知$\sigma_4$能否匹配该边缘；同一随机等待类能否对任意未知输入及伙伴实现参考替换通道。最后一个是更强任务，不是所有TR参考的普遍前提。

## 2. 固定等待的精确最优边缘误差

记$D_{\rm tr}(\rho,\sigma)=\tfrac12\|\rho-\sigma\|_1$。原指定输入的七维谱为
$$
\operatorname{spec}\Phi_u(\sigma_4)=
\left\{\frac{u^2}{2}[\times2],\ \frac{1+3u-4u^2}{5}[\times2],\ \frac{(1-u)^2}{5}[\times3]\right\}.
$$
故全部时间的误差由同一标量函数精确给出：
$$
D(u)=\frac12\left|u^2-\frac12\right|
+\frac1{10}|1+6u-8u^2|+\frac3{20}|1-4u+2u^2|.
\tag{3}
$$
令$r=1-1/\sqrt2,\ q=1/\sqrt2,\ s=(3+\sqrt{17})/8$，则$0<r<q<s<1$，并有
$$
D(u)=
\begin{cases}
\frac12-u^2,&0\le u\le r,\\
(1+6u-8u^2)/5,&r\le u\le q,\\
-3/10+6u/5-3u^2/5,&q\le u\le s,\\
u^2-\frac12,&s\le u\le1.
\end{cases}
\tag{4}
$$
第一段递减，$D(r)=\sqrt2-1$；第二段为凹函数，最小值在端点；第三段导数$6(1-u)/5>0$；第四段递增。因此
$$
\boxed{\inf_{t\in\mathbb R}D_{\rm tr}(\Phi_t(\sigma_4),r_*)
=\frac35(\sqrt2-1)=:d_*\simeq0.248528137424.}
\tag{5}
$$
唯一的$u$最优点为$q$，相应时间可按周期重复。

此时$u^2=1/2$，原面积的两个结果$4a,6a$各占$1/2$，**整个面积分布**已与$r_*$相同，$H_B$分布也相同。仍失配的是五维扇区内部：$W\mathbb C^2$上的概率比参考多$d_*$，其三维正交补相应不足。该投影是数学区分效果，不因此获得原实际读取权限。面积标定相同本身不足以认证参考准备。

### 2.1 任意原二维输入及伙伴

任意$V\in U(2)$可数学地延拓为
$\widetilde V=V\oplus(WVW^\dagger+P_{W^\perp})$。式(2)协变，且$\widetilde V r_*\widetilde V^\dagger=r_*$。对四个Pauli共轭做数学平均得到$\sigma_4$，由迹距离凸性和酉不变性，
$$
D_{\rm tr}(\Phi_u(\rho),r_*)\ge D_{\rm tr}(\Phi_u(\sigma_4),r_*)\ge d_*
\tag{6}
$$
对原输入码内任意密度矩阵$\rho$成立。该平均仅是证明，不新增Pauli实际准备。对任意$\rho_{BR}$，其$B$边缘仍满足式(6)；若比较完整输出边缘$BR$与$r_*\otimes\rho_R$，偏迹收缩也给同一下界。完整原父过程没有失去未知输入的信息，它仍可存在于工作、旗标、程序及其关联中。

## 3. 随机等待数学扩展：已知输入可以匹配

另声明输入独立的概率测度$\mu(dt)$，并仅在数学上定义$\Phi_\mu=\int\Phi_t\,d\mu$。有限维使该积分良好定义，任意连续或离散时序均包括。若
$\alpha=\int u^2d\mu,\ b=\int2u(1-u)d\mu,\ c=\int(1-u)^2d\mu$，则
$$
\Phi_\mu(X)=\alpha X\oplus[bWXW^\dagger+c\operatorname{Tr}(X)\sigma_6],
\qquad\alpha+b+c=1.
\tag{7}
$$
取等权端点$\mu_*=\tfrac12\delta_0+\tfrac12\delta_{t_*}$，有
$$
\Phi_{\mu_*}(X)=\frac12X\oplus\frac12\operatorname{Tr}(X)\sigma_6,
\qquad\Phi_{\mu_*}(\sigma_4)=r_*.
\tag{8}
$$
所以固定时界不延伸为所有混合时序的已知输入不可能结论。

保留时间标签$J$时，对已知$\sigma_4$，
$$
\omega_{JB}=\frac12|0\rangle\langle0|\otimes\sigma_4+
\frac12|1\rangle\langle1|\otimes\sigma_6,
\qquad I(J:B)=\log2.
\tag{9}
$$
两分支在$B$上支持正交，因此$S(B)=\log2+\tfrac12\log2+\tfrac12\log5$而条件熵减去后留下$\log2$。匹配边缘没有消除这份记录关联，也不把原全部$BWFC$后态变为完整Gibbs参考。

式(7)—(9)未认证实际随机调度器。若要把不同运行时长放在共同最终物理切面，选择、真实等待、停止或保持、钟、记录与控制来源均须沿既有许可兑现；不能把活跃钟或完整状态无声冻结。原自治过程在$2t_*$返回原输入，而非自动永久保留终点参考。本轮不构造这些额外权限。

## 4. 未知输入与伙伴无关准备：随机等待的精确通道界

定义参考替换通道$\mathcal R_*(X)=r_*\operatorname{Tr}X$；半diamond距离按所有输入与任意旁挂伙伴的最坏输出半迹距理解。它要求$B$不再携带原输入信息并与$R$解关联，$R$的原边缘保持；完整辅助仍可保存原输入及其关联。

以二维$R$取归一Bell输入$|\Omega\rangle$。由式(7)，实际输出的$P_4\otimes R$块恰为$\alpha|\Omega\rangle\langle\Omega|$。该四维块内的有界效果
$$
F=(P_4\otimes I_R)-|\Omega\rangle\langle\Omega|
\tag{10}
$$
是秩3投影，实际概率恒为0；目标$r_*\otimes I_R/2$在此块为$I_4/8$，概率为$3/8$。故对**任意输入独立$\mu$**，
$$
\frac12\|\Phi_\mu-\mathcal R_*\|_\diamond\ge\frac38.
\tag{11}
$$
式(10)只是范数下界见证，不授联合Bell读取权限，也不要求原未调用工厂先制备Bell态才能成立数学界。

等权端点达到此界。记$\iota_4$为向四腿扇区的等距嵌入，$\mathcal D_2(X)=I_2\operatorname{Tr}X/2$，则
$$
\Phi_{\mu_*}-\mathcal R_*=\frac12\iota_4\circ(\operatorname{id}_2-\mathcal D_2).
\tag{12}
$$
任意纯输入与伙伴的Schmidt秩至多2，等距化后写为
$|\psi_p\rangle=\sqrt p|00\rangle+\sqrt{1-p}|11\rangle$。恒等与完全退极化输出差的四个本征值为
$$
-\frac p2,\quad-\frac{1-p}{2},\quad
\frac{1\pm\sqrt{1+12p(1-p)}}4.
$$
因此其半迹距为
$$
\frac{1+\sqrt{1+12p(1-p)}}4\le\frac34,
\tag{13}
$$
在$p=1/2$取等号。混合输入经纯化、再偏迹收缩即受同界；该证明覆盖任意有限伙伴，未限制为无关联输入。式(12)的$1/2$系数给
$$
\boxed{\inf_{\mu\text{输入独立}}\frac12\|\Phi_\mu-\mathcal R_*\|_\diamond=\frac38.}
\tag{14}
$$
这是更强的统一参考准备任务。它不否定式(8)对已知输入的精确边缘匹配，也不是所有TR任务必须完成的要求；输入依赖的时序、附加测量或新辅助转换不在本比较类内。

## 5. 有限误差与资源范围

若固定时输出和目标参考分别有已独立给定的半迹距误差$\delta,\zeta$，三角不等式给下界$d_*-\delta-\zeta$。若式(14)中实际通道和目标准备各有统一半diamond误差$\eta_a,\eta_t$，则下界为$3/8-\eta_a-\eta_t$。例如两项各$10^{-3}$时，分别仍大于$0.2465$与等于保守界$373/1000$。这些是条件误差运输，不是新仪器精度认证；单个输入的迹距误差不能替代统一通道合同。

本轮完整父过程和原工作/程序/旗标来源保持。没有通过删除辅助来取得真实闭合演化的不可逆性；边缘迹距离只说明所指边缘任务尚未实现。面积仍是原Casimir角色，$\beta_*$仍是选用的匹配参考参数；没有新增同观察者温标、物理穿界流或几何动力学。

## 6. 查重与可复算范围

1274证明整条边缘轨迹、结构量与完整Gibbs边缘身份；1275检验给定完整平衡参考后的两时过程，均未求式(5)的全时最优边缘准备界或式(14)的随机等待统一准备界。本轮针对[最新参考选择筛查](../_admission/thermal_reference_selection_after1276.md)的可达性条件；不重证一般Gibbs分类。

[check.py](check.py)核式(3)—(4)的有理系数和根区间、标准等距坐标下式(2)的Kraus矩阵、全输入数学平均、最优面积分布、端点记录，以及Bell见证和Schmidt谱公式；这是1组新增有限诊断。标准坐标与原$W$由五维酉基变换等价，不重新模拟1237载体或声称新材料实现。全时间、全输入与全部概率测度的量词由上述解析证明承担，有限采样不代替证明。

核验器验证[来源基线](source_baseline.json)，绑定本文与代码，默认只读比对[保存结果](results.json)，仅首次显式`--write`排他创建结果。旧文件不改。本轮只关闭声明的等待菜单和更强随机等待准备任务，完成后停止调度/参考器件扩建；不是用户完整治理猜想、全部实际许可或GR路线的反证。
