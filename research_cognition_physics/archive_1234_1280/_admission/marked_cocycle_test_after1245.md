# 指定参考共轭的物理标记检验

2026-10-15。执行[上一项准入问题](reference_cocycle_and_reencoding_scope_after1245.md)。本页是非编号解析范围检验，不是新引力定理，不运行重复矩阵样例。CT仍为未采用猜想；125份冻结科学文件保持。

## 1. 检验结果

**指定的Connes参考共轭保持原材料代数，但通常不逐一保持辅助能源变量或联合生成元。** 随共轭运输完整物理字典能保预测；额外固定某个裸能源算子则可能只禁止改变参考。这两种要求均不能单独选出引力响应。

这比“尚未检查物理标记”更具体，但没有补出独立几何身份，不据此预约1246。它不否定指定共轭定理，也不否定未来可由独立物理字典获得选择。

## 2. 精确共轭与有限标准表示

取同一标准表示下的循环分离参考$\Omega,\Psi$，记$h_\Omega=-\log\Delta_\Omega$，$U=u'_{\Psi|\Omega}(P)$；撇号表示交换代数的cocycle，不是求导。[Witten§3.2，式39—44](https://arxiv.org/html/2112.12828v5#S3.SS2)给定共轭保持$\mathcal A$及$P$，而联合生成元的指数满足

$$
Ue^{is(h_\Omega+X)}U^\dagger
=\Delta_{\Omega|\Psi}^{-is}e^{isX}
=e^{is(h_\Psi+X)}u_{\Psi|\Omega}(s).
$$

因此不能直接把运输后的联合生成元写成$h_\Psi+X$。无限维生成元求导须另有共同可微核；以下有限材料标准表示避开这项未证资格，但不冒充Type III引力模型。

取归一、全秩有限矩阵$\rho,\sigma$；材料Hilbert空间为Hilbert–Schmidt矩阵空间，$L_bZ=bZ$、$R_bZ=Zb$，故$R_bR_c=R_{cb}$。在其与$L^2(\mathbb R,dp)$的张量积上取

$$
P=p,\quad X=i\partial_p,\quad
v(p)=\sigma^{-ip}\rho^{ip},\quad U(p)=R_{v(p)},\quad
\ell_\rho=\log\rho,\quad\ell_\sigma=\log\sigma.
$$

$v$及其导数有界，故$U$保持$D(X)=H^1$。乘法求导直接给

$$
UPU^\dagger=P,\qquad UL_bU^\dagger=L_b,
$$

$$
UXU^\dagger=X+R_{a(P)},\qquad
a(p)=-i v(p)^\dagger v'(p)
=\ell_\rho-v(p)^\dagger\ell_\sigma v(p).
$$

这里$v'=-i\ell_\sigma v+iv\ell_\rho$，已明确右乘顺序，不把左cocycle误认作保持原左代数的共轭。令$h_\sigma=-L_{\ell_\sigma}+R_{\ell_\sigma}$，则

$$
U(h_\sigma+X)U^\dagger
=h_\rho+X+L_{\ell_\rho-\ell_\sigma},
\qquad
U(X+R_{\ell_\sigma})U^\dagger=X+R_{\ell_\rho}.
$$

这些是对全部忠实有限参考成立的解析身份，无须有限采样代替证明。

## 3. 固定哪个标记，决定额外合同是什么

|合同|在上述有限模型内的结论|
|---|---|
|仅保持原左代数与$P$|任意$\rho,\sigma$均可|
|要求同一个裸$X$保持|$a(0)=\ell_\rho-\ell_\sigma=0$，故$\rho=\sigma$|
|允许裸$X$相差常数$cI$|在$p=0$得$\rho=e^c\sigma$；归一化迫使$c=0$|
|运输补偿后的$X+R_{\ell}$|任意忠实参考对均满足，未限制引力参数|

这不是说物理钟不能换参考。Witten的$X$来自一侧集体能源，区域代数使用联合量$h+X$，不能直接把裸$X$当作本项目实际正能钟。JSS使用的观察者能源$q$又位于另一种dressing表示。[Witten§3.1—3.2](https://arxiv.org/html/2112.12828v5#S3.SS1)、[JSS§2](https://arxiv.org/html/2306.01837v3#S2)

例如，仅在JSS的同一$\beta$及其规范对、模流身份前提下，写$V_\psi=e^{iH_\psi p}$、$h_\psi=\beta H_\psi$、$X=-\beta q$、$P=-p/\beta$，则

$$
V_\psi^\dagger qV_\psi=q-H_\psi
=-\frac{h_\psi+X}{\beta}.
$$

这是两种数学表示的条件字典，不证明任一辅助变量已由本项目内部钟实现，更不证明正能切口、局部完整应力或几何响应。不得把两篇文献中的能源符号直接合并。

共轭同时作用于工作状态时，原左材料读数保持；**更换参考向量不是将实际准备$\sigma$制备成$\rho$**。乘积形准备一般也不保持，不能删除新增表示中的关联。[Witten§3.3](https://arxiv.org/html/2112.12828v5#S3.SS3)

## 4. 合格的物理身份检验仍是什么

先在同一有限任务$m$中规定原准备、程序钟读时、末指针及许可未来。两侧由独立的原物理字典给$\iota_\Omega(m),\iota_\Psi(m)$，然后检验

$$
\iota_\Psi(m)=U\iota_\Omega(m)U^\dagger
$$

及相应状态、后继通道的对应。若先用右边定义左边，只得到表示协变；若未经字典便要求同一个裸算子固定，则增加了不必要的合同。这个检验可由理论效果或相关函数给出，不要求先完成全部装置或主动施加模流。

若响应参数$\gamma$已在原过程存在，逐个$\gamma$共同运输全部对象不会自动消去它。对于有界可微算子族，或另已认证共同可微域的无界族，$U_\gamma$也变时必须保留

$$
\partial_\gamma H'_\gamma
=U_\gamma(\partial_\gamma H_\gamma)U_\gamma^\dagger
+[(\partial_\gamma U_\gamma)U_\gamma^\dagger,H'_\gamma],
\qquad H'_\gamma=U_\gamma H_\gamma U_\gamma^\dagger.
$$

遗漏第二项会制造假的来源不匹配。本式只论静态参数族的表示微分；不赠予含时变换下的物理Hamiltonian或实际干预权限。

## 5. 停止线与复用

本项筛查已经回答“指定cocycle本身是否补出选择”：**尚未提供第二份独立的实际几何字典；本页有限诊断中，固定裸能源只约束参考。** 因而停止为此增加谱例、参考维数或交叉积装置；无理由重复1189的抽象同构反例。

并行回查了“完整联合体作为外层参与者”这一替代。它可按[认知方法论](../../../猜想/从认知问题到社会机制再到物理检验的研究方法论.md)理解为合并组织不能遗失既有责任，但不能以社会语言重报旧身份：[1178](../../archive_1161_1200/research_note_1178.md)已选择非中央钟耦合；[1185](../../archive_1161_1200/research_note_1185.md)已有微观结合能、完整交织来源与原钟有限记录。D-P和CC也已保留完整等待及分层打包中的交互，不能再次证明普通结合能应入账。

更关键的是，[585](../../archive_585_628/research_note_585.md)已区分共同总演化与局域lapse来源。统一外层钟率只检验总生成元，不能自动决定非均匀的局域源。[1189后来源审查](../../archive_1161_1200/_admission/after1189_source_conversion/README.md)与[1199后审查](../../archive_1161_1200/_admission/after1199_source_realization/README.md)已登记同一未完问题。故本次也不把“整块参与”列为刚发现的新机制。

重新准入必须具体给出：同一原几何承担者的独立局域响应$J_G[N]$、两份真实接入/展开程序，以及对一项仍自由响应的可区分预测。不能把总$H$的重分组当成两个实际程序，也不能直接把已选ADM或面积律放进$J_G$。这是一项剩余机制问题；没有缺失授权或已核实的外部阻碍，不因暂未找到候选就宣告GR或认知纲领失败。

本次三个子代理分别核原文与字典、有限标准表示、物理标记语义。独立审阅见[同名审阅](marked_cocycle_review_after1245.md)。不改变原阶段目标，不新增角色采用或科学核验组；三项GR共同交付仍待验收。
