# 两臂协方差与joint_tate的关系

这里比较的是用同一批观察数据得到的两个估计量，不是两个不可同时观察的潜在结局的相关性。当前模拟的目标是真正的总体TATE，target协变量的抽样波动属于其不确定性。

## iid并不意味着两个统计量独立

若不同人的O_i独立同分布，那么Cov(mean f_1(O_i),mean f_0(O_i))等于Cov(f_1(O),f_0(O))/n。不同人的交叉项消失，同一人的两项贡献仍可能相关。

已知真实nuisance函数时，target AIPW的中心化得分为

    psi_a = m_a(X)-mu_a + I(A=a)/e_a(X) * (Y-m_a(X)).

条件残差均值为零，且两个arm的残差乘积为零，因此

    Cov(psi_1,psi_0) = Cov(m_1(X),m_0(X)).

二者都使用同一批target患者的协变量。这个结果不需要识别Cov(Y(1),Y(0)|X)。若只分析两套互不重叠的独立样本原始均值，结论可能不同；若把X视为固定并改为条件/样本特定estimand，也是在改变当前推断问题。

## 有协方差不等于两个优化必然不同

一般的联合目标是Var(mu_agg,1)+Var(mu_agg,0)-2Cov(mu_agg,1,mu_agg,0)，加上分别作用于两臂source权重的惩罚。协方差依赖两套eta时，分别优化两个方差通常不等价于联合优化。

一个重要特殊情形是：target及所有保留source arms的OR人口极限都正确且可迁移。此时target的两臂聚合得分分别是

    G_t^a = m_a^*(X) + w_0^a R_t,a,

source得分是eta_j^a R_sj,a。每个R的条件均值为零，同一个人的两臂残差乘积为零。因此source的总体交叉协方差为零，target部分为Cov_t(m_1^*,m_0^*)，与eta无关。整个均值估计的交叉协方差为这个值除以n_t；其N倍极限为该协方差除以pi_t。

在这个特殊情形下，联合与分别优化的oracle目标只差一个常数，最优权重可以相同。它并没有使两臂协方差变成零。有限样本的经验协方差仍有随机误差；在OR错设分支，候选预测和target anchor的协方差一般也没有上述常数结构。保留joint_tate覆盖这些更一般的情况，不应声称其总比separate_arms有明显有限样本改善。

## 已完成的独立核验

- C1-C3、K2/8、rho0/1.5的12份真实保存拟合，共120个outer训练问题：把完整的两臂候选协方差矩阵转换为代码的increment表示，验证协方差、Wald惩罚及完整目标；最大目标差约1.9e-14。
- 独立按site写出的完整TATE重建与原估计、固定权重方差及报告的eta敏感性方差相符。
- 当前bounded DGP中，真实nuisance下target AIPW两臂均值协方差在n_t=1000时约为6.065e-6（C1/C3）、5.933e-6（C2），尽管患者iid。
- 一个一般正定协方差例子产生joint eta=(0.2051,0.5385)，而分别优化得到(1/3,1/3)；它是数学反例，不是当前DGP的总体结果。

证据保存在scratch `implementation/r11/enar_joint_tate_revision_20260924_v1/formula_checks_v2/`。稿件据此保留joint_tate、明确说明特殊等价情形，并将效率保证写为TATE的联合oracle方差界。
