# 致密物质到潮汐波形：同一状态方程与恒星分支的成熟接口

2026-10-09。**成熟采用；科学新增0、经验新增0、新认知公理0，不占1061。作者交付，待独立审阅。** 见[准入与去重](../_admission/after1060_dense_matter/selection.md)、[原始来源定位](neutron_star_tidal/sources.md)。

## 1. 共同对象与增加的覆盖

本件连接同一份致密物质描述中的总能源、静态支撑、有限物体的四极响应和双星波形：

$$
(\mathrm{EOS},p_c,\text{物质状态与稳定支})
\longmapsto (m_r,p,R,m,y_R,k_2,\lambda,\Lambda)
\longmapsto\widetilde\Lambda\longmapsto\delta\Psi(f).
\tag{1}
$$

这里EOS是状态方程，$p_c$是中央压。采用GR、冷的中微子透明、电中性且β平衡物质背景，以及静态、球对称、无旋转的barotropic完美流体扰动支；后者表示本支的扰动闭合取$\delta\epsilon=(d\epsilon/dp)\delta p$。这些均为明示物理输入。EOS包含相互作用与多体状态信息，[P981](../../archive_956_989/981/drafts/common_parent_contract_v1.md)的群、物种和参数名称并不自动给出该函数；[孤立核子来源匹配](nucleon_source_adoption.md)也不能代替致密物质匹配。门户、新组分或相变若改变本任务中的应力，须在同一EOS与扰动字典内处理。

新增覆盖是M3B／M4／M5中这条宏观共同映射。选定物质模型、状态和稳定支后，不能把支撑压力、引力源密度、Love响应与波形潮汐系数分别从不相容的模型取值，再称同一物体。它没有选择唯一EOS，未把GR或该物质状态变成认知推论，也未完成整个里程碑。

## 2. 从同一应力到有限物体的静态响应

本节先用$G=c=1$，$\epsilon$为**总能量密度**，不是只计重子的静质量密度。TOV方程是

$$
\frac{dm_r}{dr}=4\pi r^2\epsilon,\qquad
\frac{dp}{dr}=-(\epsilon+p)\frac{m_r+4\pi r^3p}{r(r-2m_r)}.
\tag{2}
$$

取$m_r(0)=0,p(0)=p_c,p(R)=0,m=m_r(R)$，$R$为面积半径。光滑段上的静态偶宇称$l=2$扰动可写为

$$
\begin{split}
H''+\left[\frac2r+B\left(\frac{2m_r}{r^2}+4\pi r(p-\epsilon)\right)\right]H'
+\left[-\frac{6B}{r^2}+4\pi B\left(5\epsilon+9p+(\epsilon+p)\frac{d\epsilon}{dp}\right)-(\nu')^2\right]H=0,\\
B=(1-2m_r/r)^{-1},\qquad
\nu'=\frac{2(m_r+4\pi r^3p)}{r(r-2m_r)}.
\end{split}
\tag{3}
$$

正则解$H\sim r^2$与外部潮汐解匹配给$y_R=RH'(R)/H(R)$。其任意总体幅度消去。式(2)—(3)采用[Hinderer等0911.3535，§II式(5)—(13)](https://arxiv.org/html/0911.3535)，式(3)采用与[0711.2420v4式(15)](https://arxiv.org/html/0711.2420v4)相同的$g_{tt}=-e^\nu$记法。

令$C=m/R$。使用**已更正**的0711.2420v4式(23)，不是原未更正版本：

$$
k_2=\frac{\frac85 C^5(1-2C)^2[2+2C(y_R-1)-y_R]}{D(C,y_R)},
\tag{4}
$$

$$
\begin{split}
D(C,y)={}&2C[6-3y+3C(5y-8)]\\
&+4C^3[13-11y+C(3y-2)+2C^2(1+y)]\\
&+3(1-2C)^2[2-y+2C(y-1)]\log(1-2C).
\end{split}
$$

恢复SI单位后，所用对称无迹四极和外部潮汐张量约定为

$$
Q_{ij}=-\lambda_{\rm SI}\mathcal E_{ij},\quad
\lambda_{\rm SI}=\frac{2k_2R^5}{3G},\quad
\Lambda=\frac23k_2\left(\frac{Gm}{c^2R}\right)^{-5}
=\frac{c^{10}\lambda_{\rm SI}}{G^4m^5}.
\tag{5}
$$

$\lambda_{\rm SI}$的单位为$\mathrm{kg\,m^2\,s^2}$，几何单位的$\lambda=G\lambda_{\rm SI}$有长度五次方量纲；$\Lambda$无量纲。式(5)的$\mathcal E$是该规范的潮汐张量，不能任意改符号却保留同一$Q$。

**界面与分支。** 自束缚物体表面能密度不为零时，应使用外侧匹配值；0911.3535式(15)给$y_R=RH'(R^-)/H(R)-4\pi R^3\epsilon(R^-)/m$。内部密度跃变也须用相应界面条件，不能把光滑式(3)直接穿过跃变。强一阶转变可有分离的稳定混合星支，见[Alford—Han—Prakash，1302.4732 §II—III](https://arxiv.org/html/1302.4732)。故式(1)保留中央压、稳定支及相关热／组成史，不把同一EOS无条件写成全局单值$\lambda(m)$，也不把某个典型EOS族中的单调关系扩大到所有相变。

## 3. 同一响应进入轨道与辐射

对声明的早期、准圆、非旋转领先潮汐近似，记$M=m_1+m_2,\eta=m_1m_2/M^2$，并在几何单位取$v=(\pi Mf)^{1/3}$。同一对星体给

$$
\widetilde\Lambda=\frac{16}{13}
\frac{(m_1+12m_2)m_1^4\Lambda_1+(m_2+12m_1)m_2^4\Lambda_2}{M^5},
\tag{6}
$$

$$
\delta\Psi(f)=-\frac{117}{256\eta}\widetilde\Lambda v^5
=\frac{3}{128\eta v^5}\left[-\frac{39}{2}\widetilde\Lambda v^{10}\right].
\tag{7}
$$

这是[Flanagan—Hinderer，0709.1915式(5)—(10)](https://arxiv.org/html/0709.1915)的成熟结果换用无量纲记法：原式的加权$\widetilde\lambda$满足$\widetilde\Lambda=32\widetilde\lambda/M^5$。相位约定为原文的$\widetilde h\propto e^{i\Psi}$。该系数同时包含潮汐改变轨道能量和总四极辐射的作用，不能另选一个与式(5)无关的辐射潮汐系数。

若$f$为探测器频率，式(7)用$v=[\pi G(1+z)M_{\rm source}f/c^3]^{1/3}$；EOS端用源系质量，$\Lambda$保持无量纲。红移／距离的资料角色仍须声明。

**有效域。** 静态$\lambda$只是零频保留系数。将它用于式(7)需要慢变且远离相关模共振、微小形变及所列轨道展开；不声称该最低阶表达式足以拟合整个观测频带。冷β平衡背景的EOS并不在任意反应时间下自动决定扰动：组成冻结、分层、超流、粘弹性、温度与自旋修正或动态潮汐若不可忽略，其响应必须重新从同一物质模型匹配。此处不认证完整迟致响应、噪声核、热并合余迹或统一数值余项，也不继承原文特定模型的400 Hz误差估计为普遍界。

## 4. 一份实际资料的限定采用

采用[LIGO/Virgo，1805.11581v2 §II—IV](https://arxiv.org/html/1805.11581v2)的已发表推断，不重拟合。其低自旋$\chi\le0.05$、共同EOS与原波形／标定合同下，公布$\Lambda_{1.4}=190^{+390}_{-120}$；谱EOS并额外要求$M_{\max}\ge1.97M_\odot$时，两个半径的边缘结果均为$11.9^{+1.4}_{-1.4}\,\mathrm{km}$。均为论文的90%可信区间，不是硬界。

前者采用EOS不敏感关系及其误差处理，后者使用谱EOS和额外脉冲星质量要求；不能拼成彼此独立的三次测量。原分析使用23—2048 Hz数据和PhenomPNRT模型，式(7)只解释共同参数入口，不替代该模型。论文自身指出，半径下界不能单独当成无条件发现潮汐效应的证据。

这仍是已用于[多信使传播](multimessenger_propagation_adoption.md)及[标准汽笛](standard_siren_amplitude_adoption.md)讨论的**同一GW170817事件**。没有新增独立实验、重建似然或相乘显著性；共同EOS与脉冲星约束是输入／推断条件，不是同一数据拟合后另获的独立验证。数值只表示固定2018版本的条件经验支持，不称截至当前的全球最新EOS结论。

## 5. 验收边界与停止线

验收的是式(1)的共同字典：在声明EOS、状态、稳定支、GR及近似阶内，结构与潮汐相位须共同变化。若某拼接不能满足该字典，先判断是否更换了EOS、支、扰动闭合或引力作用；不能直接宣告认知原则或全部P981失败。

没有计算新EOS、选择一条数值恒星序列、拟合GW资料或构造全探测器。未由宏观完美流体得到任意未知量子输入／被动参考的完整测后态，也未把观测后验变成有限时间量子仪器证书。所增加的是相容有效描述族的一条成熟跨部门接口；不删旧自由、不要求所有剩余输入先消失，不降低[ROADMAP](../../ROADMAP.md)验收要求。

[小型代数检查](neutron_star_tidal/check.py)只核式(5)—(7)的单位、系数与交换对称性；[结果](neutron_star_tidal/results.json)不是EOS预测或观测复算。
