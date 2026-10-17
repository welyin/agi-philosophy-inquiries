# 1267后最小终末几何读口：有限Ramsey指针的方向筛查

2026-10-16。非编号、未准入；只保存方向判断，不新增科学组，不修改1267或导航。接[终末来源筛查](../../../archive_1234_1280/_admission/terminal_source_screen_after1267.md)及[当前目标](../../../archive_1234_1280/_shared/goal.md)。

**判断：值得保留一个比三元Lüders更小、且可独立变分来源的候选，但今天没有完成其实际联合证书。** 最小目标只保原$\kappa=3$、$B=1,2$两程序的$a(2)$几何分离；不要求复制所有$\rho/P$分箱、不继承整个$\kappa$族在$t>2$的量子证书。真正净项是一个以原杆构造的局部尺度交互，其新增压力、钟功和自由演化可以具体计算；一般Ramsey、弱扰动或简并标签本身不计新发现。

## 1. 查重与可复用接口

- [1240](../../../archive_1234_1280/1240/proof.md)：接收后实际事件、连续钟与来源不能被重置；它没有本尺度指针交互。
- [1245—1247](../../../archive_1234_1280/1247/proof.md)及[1249正则字典](../../../archive_1234_1280/1249/canonical.md)：原杆、内部钟、记忆、完整来源及齐次强度量子化已经给定。
- [1251后内部K指针](../../../archive_1234_1280/_admission/local_collective_record_connection_after1251.md)§2：同标签指针、反作用、钟补偿与正性模板已有，不能把一般QND再计新轮；其几何无关内部质量不能直接替代下文的杆依赖质量。
- [1251后密度读取独审](../../../archive_1234_1280/_admission/local_density_reader_scope_review_after1251.md)：协变表达式也可能改变共同主部。不能仅因下文作用是局部标量就继承原CC/SR权限。
- [1259](../../../archive_1234_1280/1259/proof.md)、[1266](../../../archive_1234_1280/1266/proof.md)：原内部旗标和真实连续门的Gaussian/Duhamel工具可复用；其无新几何压力及$J+V$恒定结论不适用于尺度读取。
- [1267几何分项](../../../archive_1234_1280/1267/geometry_proof.md)：已给两程序$a(2)$差及另一个固定参数族的逆稳定；本筛查不改它的POVM或声称新Ramsey已实现该POVM。新读口只要求恢复同一$a(2)$标定量。
- [保留记录—场余项筛查](../../../archive_1234_1280/_admission/retained_record_field_residual_screen_after1251.md)：完整混合传播、初片和来源误差不能用几何子块代替。其一般伴随身份不重复推导。

## 2. 一项独立于Einstein方程的来源候选

在原三杆的正空间梯度域定义

$$
\mathcal B_\phi=\sum_{A=1}^3g^{\mu\nu}
\partial_\mu\phi^A\partial_\nu\phi^A>0,\qquad
q_{\rm rod}=-\frac12\log\frac{\mathcal B_\phi}{3k^2}.
$$

原齐次绕数$\phi^A=kx^A$给$q_{\rm rod}=Q=\log a$，不使用$T=G/\kappa$或材料密度反演。选择实有界$f(q)$，在原轨迹管道等于$e^q$；可复用原$Q$方向$C^{40}$平台，区间外明确延拓。其支持远离$\mathcal B_\phi=0$，故可将相应函数在其它$\mathcal B_\phi$域平滑置零。有效域及完成属于新输入，不能从坐标命名取得。

最小记录为简并二级系统，$P=\sigma_z$、初态$|+x\rangle$，基线承载于已计源材料的明确简并扩张；记录不是零质量。终末另一个有限标签环境也需登记，不能少算去相干承担者。以下采用原强度作用单位：

$$
h_{\rm int}=\varepsilon\alpha\,g(\theta)f(q_{\rm rod})P,
\qquad \int g\,d\theta=1,\qquad P^2=I.
\tag{1}
$$

同一原钟触发$g$，不另引外时间控制。物理相似尺度和$\varepsilon=\hbar/(V_0s^2)$必须整体运输；不能只把$\hbar$调小。若用实际总Hamiltonian表示，应按原作用倍率恢复总能源，不能把强度项直接称每个标签的物理质量。

在一个明确世界线质量密度扩张中，记局部强度$\eta=\varepsilon\alpha$、质量交互$V=\eta g fP$、数密度$n$。固定原矢量数流密度及全部材料变量，对逆度规独立变分得到

$$
T_{\mu\nu}^{\rm int}
=n\eta gP\left[
f(q_{\rm rod})U_\mu U_\nu
-\frac{f'(q_{\rm rod})}{\mathcal B_\phi}
 \sum_A\partial_\mu\phi^A\partial_\nu\phi^A
\right].
\tag{2}
$$

这是$-\int\sqrt{-g}\,nV$的变分：$\delta(\sqrt{-g}n)=
\frac12\sqrt{-g}nU_\mu U_\nu\delta g^{\mu\nu}$，
$\delta q_{\rm rod}=-(2\mathcal B_\phi)^{-1}\sum_A
\partial_\mu\phi^A\partial_\nu\phi^A\delta g^{\mu\nu}$。
在原齐次共动域，

$$
\Delta\rho=n\eta g fP,\qquad
\Delta P_{\rm iso}=-\frac{n\eta g f'P}{3}.
\tag{3}
$$

这补出了裸$f(Q)$门所没有的明确独立压力。它尚不证明局部量子杆值等于整体$\widehat Q$；该同材料量子字典须单列，不能把全体标签设成同一个算符。

钟方程仍给$\dot J=-\eta g'fP$，但$Q$运动且$P$保持时，

$$
\frac d{d\theta}(J+V)=\eta g f'(Q)\dot Q\,P .
\tag{4}
$$

因此几何读门通常有净功，脉冲首尾为零不保证$J$恢复，亦不删除式(3)压力。由同一变分，杆方程和几何动量也受反作用。这个新导数耦合的有效主部/有限因果预算须核对，不把它自动纳入旧无导数内部门的资格。

## 3. 不冻结Q的最短联合比较

在同一有效量子表示中明确选择

$$
H_{\rm read}(t)=H_c(t)\otimes I+
\varepsilon\alpha g(t)\widehat f\otimes\sigma_z.
\tag{5}
$$

$\widehat f$为新增有界尺度读口的具名量子化；它与原$a_c$效果不必完全相同。需要证明它在同一真实准备上的标定差，而非直接授予旧概率。先只用原$\kappa=3$两程序的保持窗；在新门窗不加1266另一旗标测量。

比较态保留全部原自由演化：

$$
\Psi_{\rm cmp}(t)=U_0(t)\Psi_0\otimes
\exp\!\left[-i\alpha\sigma_z
 \int_2^t g(s)f(Q_{\rm cl}(s))\,ds\right]|+x\rangle.
$$

由同一生成元的Duhamel身份，

$$
\|\Psi_{\rm read}(t)-\Psi_{\rm cmp}(t)\|
\le|\alpha|\int_2^t |g(s)|
\|(\widehat f-f(Q_{\rm cl}(s)))U_0(s)\Psi_0\|\,ds.
\tag{6}
$$

式(6)已包含自由演化与指针引起的量子反作用，不是把几何冻结后另算一个门。1266的Gaussian向量残差和同管道的Lipschitz/方差界可直接用于右侧；尚须对本$f$完成一次具体绑定，不能只凭$\varepsilon$小宣布通过。对有界来源$B$，全态误差贡献至多$2\|B\|$乘式(6)，另加(2)的新来源及实际匹配误差。

脉冲结束后的$\sigma_y$测量在比较态上有

$$
p_+=\frac12\left[1+\sin(2\alpha\bar a)\right],
\qquad
\bar a=\int g(t)a_{\rm cl}(t)\,dt.
\tag{7}
$$

这是短窗平均尺度；必须显式比较到固定$a(2)$，不能把后来状态或平均值换名为原采样。新指针终末投影在$g=0$且其基线简并时不会非选择地改变原系统边缘；新增记录/控制的实际完整来源仍须说明。

## 4. 非空的事前不等式，不冒领已实现

可只预定$\alpha=1/16$、门窗$[2,2+\tau]$、$\tau=10^{-6}$，不扫描参数。旧经典保持管在这个小窗可用$a\in[2,3]$、$|\dot a|\le4$（应在准入时用原端点与$H<11/10$绑定）。则

$$
\left|\bar a-a(2)\right|\le4\tau,\qquad
\frac{d}{da}\frac{1+\sin(2\alpha a)}2
=\alpha\cos(2\alpha a)\ge\frac1{32}.
$$

原两程序$a(2)$差$>683/10^7$，所以理想固定时刻概率差$>683/(32\cdot10^7)$，约$2.1344\times10^{-6}$。两支短窗偏差合计$\le\tau/2=5\times10^{-7}$，仍留正预算。后继若能分别把同态量子比较、实际终末匹配及统计误差控制为每支$10^{-7},10^{-7},2\times10^{-7}$，总扣除后仍为正。有限Bernoulli重复可给所选统计精度；这不是实际系统偏差或来源预算的证明。

此处相位强度是$\lambda_{\rm int}/\varepsilon=\alpha$（物理总作用图相应为耦合作用量$/\hbar$），不能独立令物理交互为零而保持读取。取原归一$C^2$脉冲时$\|g\|_\infty=35/(16\tau)$，源和钟瞬时偏移按$\varepsilon|\alpha|/\tau$控制，累计反作用按$\varepsilon|\alpha|$控制。有限读数、有限门时长和小来源要同时满足，例如

$$
\begin{split}
&2\delta_{\rm q}+2\delta_{\rm sys}+2\delta_{\rm stat}
  +2|\alpha|\,4\tau < \frac{683}{32\cdot10^7},\\
&n_{\max}\varepsilon|\alpha|\frac{35}{16\tau}
 \max(\|f\|,\|f'\|/3)<\delta_{T,{\rm allow}},\\
&2\|B\|\delta_{\rm q}
 +\delta_{\rm new\,source}(B)+\delta_{\rm matching}(B)
 <\delta_{B,{\rm allow}} .
\end{split}
\tag{8}
$$

这里$\delta_{\rm q}$必须由式(6)计算，不是可任选误差；后两行的容差须与同一有限场方程测试（如选用弱式测试，保其原范围）的归一量纲匹配，不把概率与曲率误差相加。原极小$\varepsilon$为此提供可能余量，不代替来源到场测试的实际常数。

门的采样窗可以很短，但结果只能在合法未来事件释放。不能从归一$\tau$小推出全块瞬时访问，也不能从延后可访问推出内部已经实现同一个$\widehat f$；局部杆—整体块表示及共同钟触发须具名登记。

## 5. 准入判断和停止线

本次确有一个未被旧K门填上的连接：(2)—(4)直接给几何依赖读取的独立应力和净功，(6)给同一自由过程与新读门的有限比较对象。它比继续复制三元Lüders来源更小，且不依赖把来源定义成Einstein张量。

**下一项最多只核式(6)、(8)对这一份原准备/窗口/来源字典的共同成立。** 若有效作用还不能给实际局部杆与$\widehat f$的材料匹配，或新杆导数耦合的有限域因果误差吞掉预算，就停止这一接法；不扩建全部PDE、无界指针或通用记录工厂。

本稿未给新门的数值证书、原生世界管实现、完整环境物理账或源到场余项的新数值。只保旧同模型场余项接口为后继输入；尤其不把有限$a$分离直接变成任意时空的Einstein重建，也不重写1267已经成立的条件结论。当前保存方向筛查后停止，不预约1268。
