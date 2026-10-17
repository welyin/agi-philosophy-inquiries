# QET方向筛查：真实取能不等于双向一阶几何选择

2026-10-16。非编号、0新增科学组／公理。独立核读Hotta，[*Controlled Hawking Process by Quantum Energy Teleportation*，0907.1378v3，§§1—3](https://arxiv.org/html/0907.1378)。

**原文范围。** §2预置CGHS作用；高维解释使用Einstein球对称背景。§3采用近视界平直场、外供测量切换及已知记录后的反馈；提取能量并非热Hawking辐射。式(15)—(16)的测量后源与结果比特无关。准确的反馈归一是式(31)前的$E_B=-\theta\eta+\xi\theta^2$，以及式(34)的$E_B^{\min}=-\eta^2/(4\xi)$；写成$\xi\theta^2+2\eta\theta$须另定义$\eta$，不可混用。

## 本次独立尺度判断

固定光滑紧支形状$f,p$、非零$\xi=\int(p')^2$及原文$T>x_2-x_1$，仅令测量包$\lambda(x)=r f(x)$。将其代入原文式(13)、(16)、(33)，记固定有限常数$A_1=\int(f')^2$、$C=(2/\pi)\int_0^\infty\omega|\widetilde f(\omega)|^2d\omega$及$\eta_1=\eta'(0)$，直接得到

$$
E_A(r)=A_1r^2,\qquad
\eta(r)=\eta_1r e^{-Cr^2},\qquad
\theta_*(r)=\frac{\eta_1r e^{-Cr^2}}{2\xi},
\qquad
E_B^{\min}(r)=-\frac{\eta_1^2r^2e^{-2Cr^2}}{4\xi}.
$$

因此$r$正反变化均不产生该热能源的一阶斜率。更强地，原文每个条件分支的测量后密度都是$r^2(f')^2$；任何事前固定的有限线性热流权重也从二阶开始。此结论限原菜单及固定形状，不是全部QET的分类。

若故意固定$\theta_0\ne0$，$E_B$对$r$可有线性项，但$r=0$已产生$\xi\theta_0^2>0$的反馈注能，不能同时声称仍在未扰动真空处饱和。联合小参数式为

$$
E_A+E_B=A_1r^2+\xi\theta^2-\eta_1r\theta+O(r^3\theta);
$$

在$(r,\theta)=(0,0)$没有一阶项。用$r=\sqrt{t}$重新标记只引入单向或不可微路径，不补[317](../../../archive_301_341/research_note_317.md)的双向线性前件。

**取舍。** 本文给出的真实反馈说明记录能够成为可用资源，比单纯条件重估更强；但几何响应已依赖选定引力方程，所报能量又是上述二阶量，不能提供独立的317系数选择或Jacobson热温标。外供切换／反馈的模型权限与源也未自动接到本项目。保留为[单侧过程审计](one_sided_thermal_process_after1270.md)的机制参考，停止以本菜单争取双向一阶选择；不修装置、不新增轮次。
