# 中子星组成闭合到动态潮汐：条件性成熟接口

2026-10-09。**成熟采用，科学新增0、经验新增0、新认知公理0；不占新轮。** 作者交付，待独立审阅。[来源及版本](neutron_star_dynamics/sources.md) · [有限代数核验](neutron_star_dynamics/check.py) · [结果](neutron_star_dynamics/results.json)。

## 1. 接回哪个共同对象

[既有潮汐采用](neutron_star_tidal_adoption.md)已连接同一EOS、中央压与稳定支的TOV背景、静态Love数和领先绝热波形，并明示静态EOS不自动给任意频率的扰动闭合。本件把这一边界接到**实际自引力流体的组成扰动及其引力多极**：

$$
\{\epsilon(n,x_p),\text{平衡支、扰动组成／反应合同}\}
\longrightarrow\{\delta p,\delta\epsilon,\xi,\delta g\}
\longrightarrow\widehat K_2(\omega)
\longrightarrow I^{ij}(\omega)\longrightarrow\text{轨道潮汐作用}.
\tag{1}
$$

这是多个有清楚近似域的成熟映射；下文不把相对论静态定理、牛顿动态算例与相对论响应框架拼成已经复算的同一全阶星体。P981的SM群与物种名称仍不自动提供致密物质函数、组成、温度、反应率或稳定支。

去重结论：[1009 §4—6](../../archive_1009_1043/research_note_1009.md)已给真实六维材料的静态纤维、有限动态差及最坏误差；[1000 §3—4](../../archive_990_1008/research_note_1000.md)已分开热力学与输运并核双向能源交换；[982](../../archive_956_989/research_note_982.md)已区分材料平直谱与潮汐矩。本件不再把“静态不足”计作新定理，增加的是上述星体组成—引力响应的明确采用入口。

## 2. 同背景、同静态Love为何仍未指定组成动态

采用非旋转、球对称、无弹性及无超流相对运动的完美流体。Penner等的相对论静态电潮汐方程要求

$$
\epsilon'\delta p=p'\delta\epsilon.
\tag{2}
$$

背景处于β平衡。令$c_e^2=dp/d\epsilon$为沿平衡曲线的导数，$c_f^2=(\partial p/\partial\epsilon)_{x_p}$为冻结组成导数。冻结合同$\Delta x_p=0$给$\Delta p=c_f^2\Delta\epsilon$。结合式(2)及$\Delta=\delta+\xi^r\partial_r$，得到

$$
(c_f^2-c_e^2)\Delta\epsilon=0.
\tag{3}
$$

在两导数不同处，静态解满足$\Delta\epsilon=0$；求静态Love的度规方程仍由同一背景给出。故该范围内，组成分层不改变静态电Love。导数相等处式(3)只是恒等式，不强迫$\Delta\epsilon=0$。这不是弹性壳、超流、旋转或任意动态响应的同一性结论。[Penner等，1107.0669v1，§IV.A式(68)、(70)—(76)](https://arxiv.org/pdf/1107.0669v1)。

只固定一条平衡EOS，通常只固定多变量物质函数沿平衡曲线的资料；不能据此免费指定横向导数。反过来，若完整物质函数和反应合同都已给定，$c_f$、组成响应和潮汐核须从它们共同求出，不能在保持该父对象不变时任意另选。

## 3. 有限频率：反应合同与真实自引力扰动

Andersson—Pnigouras在牛顿理论中采用$\rho=m_Bn$，并以$\beta=\mu_n-\mu_p-\mu_e$参数化偏离平衡。保留其原频域记法：

$$
\Delta p=c_e^2\Delta\rho+p_\beta\Delta\beta,\qquad
\Delta\beta=B\,R_A(\omega)\Delta\rho,\qquad
R_A(\omega)=\frac1{1+iA/\omega},\quad t_R=|A|^{-1},
\tag{4}
$$

其中$p_\beta=(\partial p/\partial\beta)_\rho$、$B=(\partial\beta/\partial\rho)_{x_p}$，$A$含反应率。这里的$\rho$是牛顿质量密度，不替换上一节的相对论总能量密度。我们只运输式(4)的模长与两极限，不把其复数符号未经Fourier约定转换搬入下一节。

对$y=|\omega|t_R>0$，直接代数给

$$
|R_A|=\frac{y}{\sqrt{1+y^2}}\le y,
\qquad |R_A-1|=\frac1{\sqrt{1+y^2}}\le y^{-1}.
\tag{5}
$$

因此局部压强系数相对平衡／冻结闭合的误差分别不超过$|p_\beta B|y$和$|p_\beta B|/y$。这是本构系数界，**不是恒星响应或波形误差界**。固定有限$t_R$再令$\omega\to0$回到平衡；先取无限慢反应得到冻结闭合，不能把两极限混同。

同一论文实际求解连续性、Euler与扰动Poisson方程；$\delta\rho$既响应潮汐驱动也产生$\delta\Phi$，不是外加一份无来源的振子。其同一$n=1$平衡多方星在不同扰动$\Gamma_1$下具有不同模态／有限频响应，静态模态和却收敛到同一Love数。常数$\Gamma_1$是作者选定的简化组成闭合，不是SM核物质已匹配的结果。[1906.08982v1，§IV式(47)、(57)—(62)、(69)—(75)，表I—IV及式(102)](https://arxiv.org/pdf/1906.08982v1)。

## 4. 怎样接回同一几何与轨道

相对论动态框架中采用$G=c=1$、$\partial_t\mapsto-i\omega$。Hegade等给出同一应力源与流体方程

$$
G_{\mu\nu}=8\pi(T_{\mu\nu}+\mathbb S_{\mu\nu}),\qquad
\nabla_\mu T^{\mu\nu}=-\nabla_\mu\mathbb S^{\mu\nu},
\tag{6}
$$

以及同一星体多极响应

$$
\widehat I_A^{ij}(\omega)
=-\frac23 R_A^5\widehat K_2(\omega)\widehat{\mathcal E}_A^{ij}(\omega).
\tag{7}
$$

式(7)的多极同时进入伴星轨道的多极力；不能为“响应”和“反作用”各选一个核。静态极限对应既有采用的$\lambda=2R^5k_2/(3G)$。原框架求流体与度规扰动并作内／外区匹配，允许组成分层与扰动耗散源；背景上$\mathbb S=0$，并限制$u^\mu\mathbb S_{\mu\nu}=0$。其流体约化及$\Delta\epsilon/(\epsilon+p)=\Delta n/n$使用这一空间应力、小偏离分支，不涵盖任意能源交换源。[2403.03254v3，§III式(1)—(11)、§IV式(13)—(16)](https://arxiv.org/html/2403.03254v3)。

式(6)展示所采用近似中的共同守恒位置，不等于已经给出反应热、中微子输运及所有二阶热源的完整模型。若用有限反应率或黏性描述能量损失，必须由同一材料模型补其接收能源的部门；本件不把局部复压强单独认证为闭合热史。

动态核也不能仅凭$k_2(0)$在全频段任意拟合。选择物质状态、反应与边界后，它须由同一扰动问题决定。含低频g模时，零频附近的普通Taylor展开及统一余项可能失效；不能从“频率小”自动领取整个组成族的$O(\omega^2)$波形界。这里不复制8PN数值相位、数值EOS或观测后验。

## 5. 本次核验、保留自由与停止线

[小型检查](neutron_star_dynamics/check.py)只核式(5)、静态兼容的代数及原模态表达式的非零极点权重；表II的一个已发表g模数值只作转录／归一校准，不当作精确频率、完整模态和或新恒星计算。默认运行只读对照结果。

本件减少的是拼接自由：同一已声明物质闭合必须同时承担压强扰动、自引力多极和轨道作用。它没有由认知原则选出EOS、反应率或组成，也没有宣称在完整父模型固定后另有可任意选择的动态系数。未认证有限观测窗的总误差、完整量子记录／参考、真实检测器、合并期或完整GR波形。

2026年的新匹配论文[2606.19446v2](https://arxiv.org/html/2606.19446v2)已完成[独立PDF与HTML来源审计](../_admission/dynamical_tidal_source_after1062/source_audit.md)：按其印刷定义，式(57)、(65)—(66)的running符号及式(63)—(64)的量纲衔接不一致。式(65)—(71)的所核内部代数可自洽，但不能原样接到式(57)、(66)；详见审计中的两种β记号与核验范围。这是项目核对，不是官方勘误或对整篇理论的否定。本批**仍不采用完整8PN系数或Fisher预报**，上述已核成熟接口不受此局部问题影响。

结论：新增的是中子星共同动态接口的计0覆盖，未形成新的有限任务科学轮。后续只有在具体频带、同一物质／反应／热源和误差合同下需要比较实际预测，才推进相应计算；不为重复静态非识别再造星体模型。
