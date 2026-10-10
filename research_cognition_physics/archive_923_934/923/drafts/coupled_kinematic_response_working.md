# 923工作推导：用同一有限过程合并运动学误差的传播交叉项

接[922背景读口推导](../../../archive_894_922/922/drafts/background_derivative_proof.md)与[联合误差系统](../../../archive_894_922/922/drafts/joint_error_transport_working.md)。这是下一项实施对象，不另计正式轮次，尚无新数值结果。

## 1. 一个相对发展方程包含全部交叉项

保留同一原背景、受迫方向、接收与测试模式，记Z=(Y,u,z,A)，

$$
 \mathbb F(Z)=\big(F(Y),-iH_D(Y)u,-iH_D(Y)z,DF(Y)A+S_\chi(Y,u)\big).
$$

原有限近似轨迹为Zhat(t)。定义W_ε(0)=0，并计算

$$
 \dot W_\epsilon=\mathbb F(\widehat Z+W_\epsilon)-\mathbb F(\widehat Z)
 +\epsilon(f_Y,0,0,f_A).
$$

f_Y是922背景运动学残差的负投影，f_A是921切向残差的负投影。它们在本分项比较中是固定的原轨迹函数；这是孤立误差分量的运输，不声称改变后的轨迹仍有同一精确残差。

减去相同的F(Zhat)确保即使Zhat不是精确解，W_0也恒为零。此处不能直接重新发展一个邻近背景，再把原轨迹的离散差混进所求方向导数。

有限维光滑条件下对实参数ε微分，(b,η,η_z,C)=∂εW|0满足

$$
 \dot b=L b+f_Y,\qquad
 \dot\eta=-iH_D\eta-iDH_D[b]u,\qquad
 \dot\eta_z=-iH_D\eta_z-iDH_D[b]z,
$$


$$
 \dot C=L C+D^2F[b,A]+D_Y S_\chi[b]+D_uS_\chi[\eta]+f_A.
$$

该式直接由链式规则得到；完整保背景改变传播算符及接收应力的作用。Dirac源依赖旋量及共轭，必须取实微分，不能嵌套普通complex-step将共轭项遗漏。

## 2. 同一原权重的模式变化可直接微分

固定坐标点令n=z†βu，d=u†βu，w=n/d，则

$$
 \delta n=\eta_z^\dagger\beta u+z^\dagger\beta\eta,\quad
 \delta d=\eta^\dagger\beta u+u^\dagger\beta\eta,\quad
 \delta w=(\delta n-w\delta d)/d.
$$

对时空坐标再求导得到δ(∂w)。用同一有限Fourier/Hermite模式及其导数实施，不能随意拼接值和导数。914的Ward配对对(w,∂w)线性，因此模式变化项可直接配对；922已计算的移动材料点效应不得重复计入。

总读口首阶为K(C)+D_Y K[b]+K_{δw}(A)。用一个共同有限变形下的原读口差分核验此分解，可检查交叉项的符号与重复计账。本轮不能凭单独字段小就删去任何本阶源。

## 3. 实施与停止界线

复用N17、原T、原参考/材料路径、原准备和原u,z。需保存原全部正则A与u,z的插值，不能由配置导数擅自重造canonical数据。先作一次有限联合计算与独立幅度检查，不开单项收敛扫描。

本项完成也只合并了一个确定的运动学误差分量：初值/几何约束、其它残差、连续重建与配点算符差、真实物理投影传播、完整872反馈及Q/Q_eff字典仍另列。若希望改成有限模型本身作为共同候选，需要证明其来源/观测匹配，不能默认上述比较误差就是该候选的物理错误。

所有工作继续按声明有限域的有效描述验收；没有无截断完成或物理最小尺度要求。
