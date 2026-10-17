# 记录选区的同对象接续：有限性、信息预算与非平衡边界余额

2026-10-17。非编号解析接续；没有新增数值实验或仪器。本文的原对象是589/598/603的完整有限图Gauss材料、617切分及724原字段读取，全部旧作用、CAR、规范群与原相互作用保留。不是1274有限SU(2)载体的物理映射。

## 1. 共同输入与三个不同的“区域”角色

固定原有限图、正背景几何及全部耦合，令$H$为原完整Hamiltonian。取$c$足够大，至少覆盖637的形式比较常数$C_*$，使$A_H=H+cI\ge1$且该轮的正比较能源$K_0\le A_H$。移位仅供分析，不从物理能源或来源扣除。由[603至637的既有热迹与能源域](../../../archive_629_652/research_note_637.md)，
$$
\sigma_\beta=\frac{e^{-\beta H}}{Z},\qquad \beta>0
\tag{1}
$$
是忠实正常态，具有有限能源及熵。本文条件性指定它为输入；正常Gibbs态存在不等于实际热化或制备权限已经实现。

[724](../../../archive_702_741/research_note_724.md)使用651同一原物质壁
$$
f_v=s_v-\lambda|\phi_{H,v}|,
$$
其中$\lambda$为原先固定的选壁参数。$f_v$在原目标内有界，弱导数及原能源形式已获证；不要求Higgs轴处的经典光滑速度。

本稿区分：

1. 原连续经典关系壁及其因果类型，需要651的额外经典范围；
2. 本次实际有限标签所选的原图节点集合；
3. 该节点集合在617匹配描述中的区域代数及约化态。

本文直接使用第2、3项，不把它们自动等同于第1项的精确连续边界，也不自动认证Lorentz面积、实际同观察者时钟或跨界热流。背景几何仍为该有效材料的输入；未另外增加几何Hamiltonian，也未从物质选区反推出几何动力学。

## 2. 原完整读取、有限记录与注能

设分辨率$\sigma>0$固定。每节点只实施一次724的Gaussian字段读取：
$$
K_{y_v}(f_v)=(2\pi\sigma^2)^{-1/4}
e^{-(y_v-f_v)^2/(4\sigma^2)},\qquad
K_{\mathbf y}=\prod_vK_{y_v}(f_v).
$$
以事前固定的有限区间分割读数，每个标签向量$\ell$对应一个积分域$I_\ell$。包括全部尾部区间、全部结果：
$$
\mathcal I_\ell(\rho)=\int_{I_\ell}K_{\mathbf y}\rho K_{\mathbf y}^\dagger\,d\mathbf y,
\qquad
\sum_\ell\mathcal I_\ell=\mathcal M,\qquad
p_\ell=\operatorname{Tr}\mathcal I_\ell(\sigma_\beta),
$$
$$
\rho_\ell=\mathcal I_\ell(\sigma_\beta)/p_\ell\quad(p_\ell>0),
\qquad
\bar\rho=\sum_\ell p_\ell\rho_\ell=\mathcal M(\sigma_\beta).
\tag{2}
$$
区间积分后通常不是单个Lüders Kraus后态。未公开的连续读数正确求和，有限记录均保留；不将未知输入预先分类或免费重测。仪器可张量恒等作用于被动伙伴，但只保证非选择伙伴边缘不变，系统关联和条件伙伴可改变。以下熵预算针对所列材料与记录，不额外声称任意伙伴联合相对熵的同一上界。

724式(6)—(7)在原完整$H$上给
$$
\operatorname{Tr}H(\bar\rho-\sigma_\beta)
=\operatorname{Tr}(\sigma_\beta D_{\rm meas}),\qquad
D_{\rm meas}=\sum_v\frac{\hbar^2}{8w_v\sigma^2}\Gamma_{f_v},
$$
$$
0\le D_{\rm meas}\le B_{\rm meas}I,\qquad
B_{\rm meas}=
\frac{\hbar^2M(1+\lambda^2)}{8\sigma^2}\sum_vw_v^{-1}<\infty.
\tag{3}
$$
$w_v$是原正节点体积权，$\Gamma_{f_v}$是原目标逆度量中的梯度平方。式(3)保完整质量、规范与相互作用，并非只算自由字段的一项；它仍是原材料系统的测量反作用账，不是自主装置、控制、指针和实际四动量来源的全套证书。这些权限和代价不得免费附送。

由$A_H\ge1$及(3)，
$$
\sum_\ell p_\ell\operatorname{Tr}(A_H\rho_\ell)
=\operatorname{Tr}(A_H\sigma_\beta)+\operatorname{Tr}(\sigma_\beta D_{\rm meas})
\le E_{\rm in}+B_{\rm meas}<\infty.
\tag{4}
$$
每个非零概率分支有有限形式能源，界可随$p_\ell^{-1}$变大；不授稀有分支的统一条件上界。

## 3. 非选择熵与保记录的相对熵预算

这里只为证明熵性质使用一个数学分解，不另授随机控制器。令
$$
U_{\mathbf k}=\exp\!\left(i\sum_vk_vf_v\right),\qquad
d\nu_\sigma(\mathbf k)=
\prod_v\sqrt{\frac{2\sigma^2}{\pi}}\,e^{-2\sigma^2k_v^2}\,dk_v.
$$
在共同字段谱表示中，两种写法都将矩阵核乘以
$\exp[-\sum_v(f_v-f'_v)^2/(8\sigma^2)]$，所以作为原空间上的正常通道，
$$
\mathcal M(\rho)=\int U_{\mathbf k}\rho U_{\mathbf k}^\dagger\,d\nu_\sigma(\mathbf k).
\tag{5}
$$
$\nu_\sigma$归一，各$k_v$均值为零、方差为$1/(4\sigma^2)$。这些规范不变有界酉保原Gauss空间；724的弱梯度界保证每个固定$\mathbf k$保原形式域，能源增长至多为$\mathbf k$的二次函数。高斯平均能源由式(3)有限。对随机酉混合使用积分形式的熵凹性（由有限混合和正能源谱截断延拓），结合式(4)及637的有限熵控制，得到
$$
S(\bar\rho)\ge S(\sigma_\beta).
\tag{6}
$$
它是638二结果随机酉证明在724 Gaussian通道上的直接接续，不是新的普遍熵定理。

设$H(p)=-\sum_\ell p_\ell\log p_\ell$。记录数有限，因而$H(p)<\infty$。若每节点至多$m+1$个标签，则$H(p)\le |V|\log(m+1)$。以同一真实概率构造两份比较态
$$
\Omega=\bigoplus_\ell p_\ell\rho_\ell,\qquad
\Xi=\bigoplus_\ell p_\ell\sigma_\beta.
\tag{7}
$$
$\Xi$只是用于比较的参考直和态，不是“同一次读取作用于$\sigma_\beta$的输出”，也不声称实际能够无成本恢复$\sigma_\beta$。否则初态与自身的零相对熵会被误用于当前有代价的读取。

Gibbs对数及有限熵给
$$
D(\bar\rho\Vert\sigma_\beta)
=\beta\operatorname{Tr}(\sigma_\beta D_{\rm meas})
-[S(\bar\rho)-S(\sigma_\beta)]
\le\beta B_{\rm meas}.
\tag{8}
$$
又有精确经典—量子分解
$$
\sum_\ell p_\ell D(\rho_\ell\Vert\sigma_\beta)
=D(\bar\rho\Vert\sigma_\beta)+
\left[S(\bar\rho)-\sum_\ell p_\ell S(\rho_\ell)\right].
$$
方括号为Holevo量，处于$[0,H(p)]$，故
$$
\boxed{\sum_\ell p_\ell D(\rho_\ell\Vert\sigma_\beta)
\le\beta B_{\rm meas}+H(p).}
\tag{9}
$$
所有相对熵均以同一$H,\beta,Z$定义；各$\rho_\ell$不必为Gibbs态。有限性还可直接由
$D(\rho_\ell\Vert\sigma_\beta)=\beta\operatorname{Tr}(H\rho_\ell)+\log Z-S(\rho_\ell)$
和(4)核验。$H(p)$是信息资源项，不是能源；不把记录熵乘除温度后擅自称为已支付的功。

## 4. 由同一记录选区，参考与后态共同切分

事前固定一个有限规则$\ell\mapsto A_\ell\subset V$，例如选择指定标签的节点并集；补集记为$A_\ell^c$。区域可以不连通。每个区域只切原图中跨两侧的边。令
$$
J_\ell:\mathcal H_{\rm Gauss}\longrightarrow
\mathcal H_{A_\ell}^{\rm ext}\widehat\otimes
\mathcal H_{A_\ell^c}^{\rm ext}
$$
为617的保全部规范匹配资料等距，$J_\ell^\dagger J_\ell=I$。帽号提示原CAR分次次序；状态总偶，沿原次序处理部分迹。其值域是匹配子空间，不授两侧任意独立准备。定义正常CPTP区域映射
$$
R_{A_\ell}(\rho)=\operatorname{Tr}_{A_\ell^c}
(J_\ell\rho J_\ell^\dagger).
\tag{10}
$$
对每个记录，必须对真实后态与同一未测量参考用相同的$J_\ell$和部分迹。对式(7)同时使用记录控制的直和通道并保记录，数据处理给
$$
\boxed{
\sum_\ell p_\ell
D(R_{A_\ell}\rho_\ell\Vert R_{A_\ell}\sigma_\beta)
\le \sum_\ell p_\ell D(\rho_\ell\Vert\sigma_\beta)
\le \beta B_{\rm meas}+H(p).
}
\tag{11}
$$
参考区域随公开标签改变，但不是逐分支另拟合温度、Hamiltonian或权重。所用$p_\ell$只是式(7)共同的经典比较权重，不是新物理准备。

637为每个切分构造$K_{A_\ell}$，其热迹有限，并给
$$
S(R_{A_\ell}\rho_\ell)
\le \eta\,\operatorname{Tr}(A_H\rho_\ell)
+\log Z_{A_\ell}(\eta)<\infty,\qquad \eta>0.
\tag{12}
$$
有限图只有有限个节点子集，可对这里实际出现的有限区域取热迹上界最大值；结合(4)得到加权区域熵有限。规范中心标签熵、表示维数项及块内熵因此各自有限。没有用有限总热迹直接替代这份区域域证明。

允许在记录之后按旧许可等待同一$H$的有限时间，并将等待纳入$\rho_\ell$；此时全局能源、熵、对同一$\sigma_\beta$的相对熵不变，式(11)仍成立。这里不重新读取、增加控制或改变选区规则。记录决定区域、区域随后如何做真实连续运动，是不同任务。

## 5. 同一条件后态的保中心分解与一般非平衡余额

对固定$\ell$，按1221把原物理空间写成切口标签$\lambda$的直和：
$$
\mathcal H_{\rm phys}\simeq
\bigoplus_\lambda(\mathcal H_{a,\lambda}\widehat\otimes\mathcal H_{b,\lambda}),
\qquad
J_\ell|_\lambda=I_{ab}\widehat\otimes|\Omega_\lambda\rangle.
$$
$|\Omega_\lambda\rangle$是维数$d_\lambda$的强制配对，
$d_\lambda=\prod_{e\in\partial A_\ell}d_{\lambda_e}$；
它不同于记录标签$\ell$。令
$$
q_{\lambda|\ell}=\operatorname{Tr}(P_\lambda\rho_\ell),\quad
\rho_{\lambda|\ell}=P_\lambda\rho_\ell P_\lambda/q_{\lambda|\ell},
$$
$$
\tau_\ell:=\tau[\rho_\ell]
=\bigoplus_{q_{\lambda|\ell}>0}
q_{\lambda|\ell}
(\rho_{a,\lambda|\ell}\widehat\otimes\rho_{b,\lambda|\ell}).
\tag{13}
$$
零权块取零。$\tau_\ell$保原两侧边缘和匹配标签，属于同一物理空间；该处方非线性，不是可免费实施的去关联通道。1221的双向形式控制确保$\tau_\ell$也有有限完整$H$能源。其熵、下述全部相对熵和中心项有限，不涉及无穷减无穷。

1221对任意正常总偶有限能源态已证明
$$
\mathcal C_\ell:=D(\rho_\ell\Vert\tau_\ell)
=S(\tau_\ell)-S(\rho_\ell)\ge0,
$$
$$
I(A_\ell:A_\ell^c)_{J_\ell\rho_\ell J_\ell^\dagger}
=H(q_{\cdot|\ell})+
2\sum_\lambda q_{\lambda|\ell}\log d_\lambda+\mathcal C_\ell.
\tag{14}
$$
这不能把全部互信息都当成可消去关联或热熵。

将原全部$H$按实际支撑分为区域内项和完整跨界簇项：
$$
H=H_{A_\ell}+H_{A_\ell^c}+H_{\partial\ell}.
$$
依原分账，完整电星留源节点侧，跨界梯度、磁星及CAR跳跃均保留。因为$\tau_\ell$保两侧边缘，
$$
\Delta E_{\partial\ell}:=
\operatorname{Tr}[(\tau_\ell-\rho_\ell)H]
=\operatorname{Tr}[(\tau_\ell-\rho_\ell)H_{\partial\ell}].
$$
展开同一Gibbs参考的两个相对熵，得到
$$
\boxed{
\beta\Delta E_{\partial\ell}
=\mathcal C_\ell+
D(\tau_\ell\Vert\sigma_\beta)-D(\rho_\ell\Vert\sigma_\beta).
}
\tag{15}
$$
这是1221一般匹配比较与Gibbs对数的直接推论。对$\rho_\ell=\sigma_\beta$才退回原双相对熵热身份。一般$\Delta E_{\partial\ell}$可有符号，既不是已发生的跨界热流，也不是已经实现的最小功。

式(9)控制了不能省略的非平衡资源；例如精确地有
$$
\sum_\ell p_\ell\mathcal C_\ell
\le\beta\sum_\ell p_\ell\Delta E_{\partial\ell}
+\beta B_{\rm meas}+H(p).
\tag{16}
$$
这是保留所有标签的余额，不声称每个条件分支有与$p_\ell$无关的统一界。各项加权有限，因为1221的形式比较常数可在固定图有限个切分中取最大值。

## 6. 独立面积、测量来源与区域变化的不同含义

在1225共同几何域内，固定
$$
h=w^{2/3}G,\qquad \det G=1,\qquad
a_i=\sqrt{(\operatorname{cof}h)_{ii}}>0,
$$
由该几何独立赋予记录选区面积
$$
A_{\partial,\ell}=\sum_{e\in\partial A_\ell}a_{i(e)}.
\tag{17}
$$
可以选取724一般固定图范围与1225的125节点周期图的交集，并保持全部原材料。对一个选定区域，每条跨区域与补集的边在式(17)中只计一次；若另将多个区域面积相加，共用内部界面会被双方分别计入，不能把这称为新增外边界面积。式(17)不借熵定义面积。但1225校准使用另一些事前准备和效果；它们的实际误差、短时及高阶矩条件不自动授予724条件后态中的新读取。

1226的规模一致热上界使用完整平移Gibbs参考；$\rho_\ell$通常不满足这一前提，且选区与局部能源相关。因此不能从式(14)或(16)删去非平衡项，再用原$e_\beta$无条件替换逐分支资源，也不能据此宣布普适面积熵系数。

724的$J_\ell$在固定图和标签下不依赖外部几何参数。故原$H$及$G_\gamma=\partial_\gamma H$可按同一$J_\ell$运输。若真实条件后态、阈值或几何一起改变，必须另外保态导数、分支概率导数和724已有的读数域端点项；725证明同一次位置读取对法向平方也有非零反作用。这里没有新推无穷维熵导数，不把记录域端点通量认作时空$T_{\mu\nu}$通量。

固定$h$时，标签选择不同可以改变式(17)，但度规没有因此改变。旧651的连续类时壁、完整边界变分和在壳来源，不能仅由这些离散后态自动领取；728的另一个费米基态带资格也不补上这一身份。自主测量装置、真实世界管及实际热边界仍按各自已证范围登记，不借本接续新增权限。

## 7. 准入价值与明确停止线

可直接使用的旧证据为：

- [724](../../../archive_702_741/research_note_724.md) §§2—6：原壁、实际Gaussian区间instrument、完整材料注能、记录选区及匹配、移动阈值来源；
- [725](../../../archive_702_741/research_note_725.md) §§3—6：同后态法向平方、能源域和来源运输；
- [637](../../../archive_629_652/research_note_637.md) §§3—7、[638](../../../archive_629_652/research_note_638.md) §§3—4：任意有限能源区域熵及同参考条件信息预算；
- [1221](../../../archive_1219_1233/1221/proof.md) §§2—6：任意后态的匹配比较、有限形式能源、保中心互信息及同来源；
- [1225](../../../archive_1219_1233/1225/proof.md)、[1226](../../../archive_1219_1233/1226/proof.md)：独立胞元面积字典，以及仅限相应平衡域的规模一致上界；
- [644](../../../archive_629_652/research_note_644.md)：固定区域同参考复制与来源已有；不把其Rényi指数2身份当成本次普通熵或移动区域导数。

本次仅把这些证书接在同一次记录后态上，得到式(11)、(15)及相应有限性，不声称原创的一般数学定理。没有计算新的有限矩阵、热谱或控制器，也不重复旧科学核验。

这是一份共同区域—参考—资源的条件正接口，不是$S=sA$、$\delta Q=T\delta S$、几何反作用或强生成的验收。它也不是1274—1275边界载体的空间实现。下一研究动作若只重复测量注能、记录概率或已有面积上界，应停止；真正的身份检验须使用本文同一后态和不可省略的非平衡余额。
