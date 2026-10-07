# 第1006轮：常规超导与整体认知操作中的集体相干

2026-10-08。接[1005](research_note_1005.md)与[入口](1006/drafts/STATUS.md)。正式1006／累计3786；新增物理试验组0。本轮是成熟机制采用及最小解析审计，**整体目标未完成**。

## 1. 本轮回答什么

整体假说目前已解释普通材料、有限供给与辐射反馈，但P08尚未具体解释常规超导。认知组织可以同时包含可重复读取的经典记录和宏观量子有序吗？有复序参量、配对或零电阻，分别够不够说明磁屏蔽和实际相位比较？

**结论：它们可以在指定变量和工作域中相容，但不能互相替代。** 常规配对作用、集体相位刚性、电磁响应和选择性记录有不同物理职责；共同解释不能只把它们都称作“认同”。本轮采用成熟BCS、近临界GL、London及弱结机制，给出这些职责如何衔接，不新增认知公理。

## 2. 旧成果直接复用

|旧结果|复用内容|本轮不冒领的结论|
|---|---|---|
|[355](../archive_342_369/research_note_355.md)|中性GP凝聚态、振幅／相位与Bogoliubov接口及有效域|中性声学几何不是带电Meissner效应或GR|
|[1000](1000/material_phase_transport_adoption_v1.md)|相中慢变量、热—力配对及材料来源|一个热弹性例不说明超导配对|
|[1004供给](1004/resource_preparation_adoption_v1.md)|关系参考须有载体与共同访问|时间平移资源定理不能原样改成粒子数／电荷定理|
|[1005](1005/radiation_formation_adoption_v1.md)|辐射再吸收和热／动量反馈|冷光源或超导体不是无限参考与无成本装置|

171是超导量子比特实验平台引用；657的BCS Gram是费米配对计算，960的London是色散力。关键词相同不等于已经解释了本轮现象。旧空间372、修订382—384、386或425、522—523保持。

## 3. 机制与最小解析关系

[机制补充](1006/conventional_superconductivity_adoption_v1.md)完整给出参与者、物理采用和适用边界。常规链为：电子—晶格作用与屏蔽排斥竞争 → 合适通道的集体配对 → 准粒子谱与相位刚性 → 规范不变的电磁响应。不是由预置GL负系数反向证明认知必然产生配对。

以 $\psi=f e^{i\theta}$、配对电荷 $q=-2e$ 为约定，局域静态自由能中的梯度项是 $f^2(\hbar\nabla\theta-q\boldsymbol A)^2/(2m_*)$。因此

$$
\boldsymbol p_s=\hbar\nabla\theta-q\boldsymbol A,
\qquad \boldsymbol j_s=\frac{qf^2}{m_*}\boldsymbol p_s.
$$

同时变换 $\boldsymbol A\mapsto\boldsymbol A+\nabla\chi$、$\theta\mapsto\theta+q\chi/\hbar$，$\boldsymbol p_s$ 不变。真实刚性针对规范不变的流动关系，不是对纯规范重命名收费。

常幅、无涡旋、静态局域Meissner支上，对电流取旋度并联立同一Maxwell磁场：

$$
\nabla^2\boldsymbol B=\lambda^{-2}\boldsymbol B,
\qquad \lambda^{-2}=\frac{\mu_0q^2f_0^2}{m_*}>0.
$$

相应弱场半空间边界给指数穿透。它说明磁屏蔽依赖集体响应及边界，而非单凭零电阻；理想导体的 $\boldsymbol E=0$ 只约束磁场的时间变化。静态自由能不提供达到该支的动力学、热噪声或仪器寿命。

实际弱连接读到的关系是

$$
\varphi=\theta_R-\theta_L-\frac q\hbar\int_L^R\boldsymbol A\cdot d\boldsymbol l,
\qquad E_J(\varphi)=-\mathcal E_J\cos\varphi.
$$

此处 $\mathcal E_J>0$ 是采用的结耦合幅值，不是免费获取的认知能力。固定总粒子数可以保留两侧的关系关联；绝对裸相位不是观测对象。关系参考不要求外部绝对标准，但仍要求真实弱连接、路径和制备。

## 4. 对整体假说的修订

1. **集体有序不等于整体纯态。** 有限温度的配对关联可与准粒子、晶格和环境的热运动并存。
2. **经典记录不等于所有相干消失。** 记录可以由其他变量和时间窗承担；但若探测区分待保留叠加，仍会回扰它。
3. **维持某个序参量不等于保护任意量子信息。** 逻辑态可失相干而材料仍超导；反之，超导有序本身也有温度、场强和相位滑移等限制。
4. **零电阻不等于无成本生命周期。** 被动持久电流不必持续做功；实际制备、低温环境、调控、读出和复位仍须有能源／熵去向。

这些补强P08及P02／P04／P09的相容解释；不把有限温度的吸收混成静态直流电阻，不把配对相位资源混成1004的时间相干。数学957、候选可替换地位和应用目标保持。

## 5. 验收与停止线

[审计代码](1006/superconductivity_adoption_audit.py)与[结果](1006/superconductivity_adoption_results.json)记录采用证据、四项解释义务和未核事项；[交付核验](1006/research_round_1006_checks.json)保护历史并核对入口。没有新增物理模拟，累计3786保持。

宏观有序机制现已具体补充；现实材料Tc、完整微观形成、噪声与寿命、非常规相、真实器件以及共同引力来源仍未在本轮验证。不存在必须先完成这些细节才能采用上述解释的证据。

接[1007](1007/drafts/STATUS.md)补中微子产生—传播—探测及介质味转换；直接复用539的质量接口，不重算质量拟合或W995产额。

## 文献及读取范围

- [BCS 1957原文](https://journals.aps.org/pr/pdf/10.1103/PhysRev.108.1175)：核读全文相关节，采用常规配对与电磁响应机制。
- [Gor'kov 1959官方索引](https://www.jetp.ras.ru/cgi-bin/e/index/e/9/6/p1364?a=list)：索引／摘要已核，PDF抓取失败；用于近Tc微观到GL联系，不声称本项目重证该推导。
- [Josephson 1962原文](https://www.physics.umd.edu/courses/Phys798C/AnlageFall25/Josephson-1962-Possible%20New%20Effects.pdf)：已核两区域数目关联及隧穿流。
- [Mattis—Bardeen 1958](https://journals.aps.org/pr/abstract/10.1103/PhysRev.111.412)：核读出版社摘要，采用超导态仍有频率依赖散射／吸收的机制范围。
