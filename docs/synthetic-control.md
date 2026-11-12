# 合成对照与空间 placebo 诊断

全部厂区、周值与结果均为合成资料，不代表任何真实项目的效果。实现位于
`analysis_data/synthetic.py`，命令为：

```sh
python -m analysis_data.cli analyze data/normal.json
python -m analysis_data.cli analyze data/hand-example.json
python -m analysis_data.cli analyze data/balanced-panel.json
```

结果默认写到标准输出；`--output` 可写入其他路径，但拒绝覆盖输入文件，
并在输出中保留输入的 `source_sha256`。

## 面板与协议

- 只支持一个处理厂区、一次处理、时间轴为连续整数周，且面板平衡无缺失。
- `treated` 是唯一受处理厂区；`sites` 中其余厂区都视为始终未处理的 donor。
- `treatment_week` 及以后为 post；之前为 pre。
- pre 的最后 `validation_weeks` 个连续周是留出段，更早的周是训练段；两段
  均不能为空。留出段只用于选 lambda，最终权重在**全部 pre 周**上重拟合。
- 不做 z-score、中心化或任何标准化，也不加入截距。
- 处理后的观测、donor 值和差距只在权重冻结后计算。修改任何 post 数值都
  不会改变 lambda 选择、标准化方式或权重；只会改变 post 差距和 post-MSPE。

## 权重目标与 KKT

对训练周 `t=1..T`，令 `X` 为 donor 周值矩阵，`y` 为处理厂区周值，权重为
`w`。优化目标是

```text
f(w) = (1/T) ||y - X w||^2 + lambda * sum_j w_j^2
subject to w_j >= 0, sum_j w_j = 1.
```

梯度为

```text
g = (2/T) X' (X w - y) + 2 lambda w.
```

输出报告目标值、完整梯度、等式乘子 `nu`、`mu = g + nu`（项目符号约定为
`L = f + nu(sum w - 1) - mu'w`）、等式残差、下界违反量、支持集平稳性残差、
非支持集对偶可行性违反量和互补残差 `|mu_j w_j|`。`max_residual` 是这些
量的最大值。求解器使用支持集/有效集方法与单纯形投影，不对 KKT 矩阵直接
求逆；重复列、完全共线列或奇异方向由对称特征分解的最小范数解处理（奇异
方向被截断），仍返回非负且和为 1 的合法权重。

`solver.termination_reason` 为 `kkt_converged` 才表示已满足配置的 KKT
容差；`budget_exhausted`（预算耗尽）和 `stalled_without_dual_violation`
都明确标记为 `optimal: false`，不会被当作已求得最优。

## lambda 选择与同分规则

`lambdas` 是预先给定、非负的候选集合。对每个 lambda：

1. 只在训练段拟合权重；
2. 用该权重计算留出段的均方误差；
3. 取留出 MSE 最小者；若某候选的 MSE 与最小值的绝对差
   `<= tie_tolerance`（默认 `1e-12`），视为同分，取其中**较大**的 lambda；
4. 用选中的 lambda 在全部 pre 周重拟合，得到最终权重。

`lambda_selection.candidates` 保留每个候选的训练 MSE、留出 MSE、目标值、
KKT 残差与终止原因，便于复查；`fit.fit_weeks` 明确列出最终重拟合使用的
全部 pre 周。

## 主结果与失配提示

`weekly` 每周给出 `observed`、`synthetic`、`gap = observed - synthetic`
和 `period`（pre/post）。`fit.post_average_gap` 是 post 周差距的算术平均；
`fit.pre_msp`、`fit.post_msp` 是两个时期的均方差距。主统计为

```text
same_period_ratio = post-MSPE / pre-MSPE
```

零分母按项目约定处理，绝不在 JSON 中写 `NaN`/`Infinity`：

- pre-MSPE 与 post-MSPE 都为 0：`ratio_kind = "zero_over_zero"`，比率定义为 `0`；
- pre-MSPE 为 0、post-MSPE 为正：`ratio_kind = "positive_over_zero"`，
  比率显示为字符串 `"+infinity"`，异常不会被除零丢掉；
- 其他情况给出有限数值。

若 `pre-RMSPE > pre_rmspe_max`，结果照常输出，但 `diagnostics` 给出
`mismatch_warning`，并且不会用 post 结果去重新挑选 donor 或 lambda。

合成对照只按可观测的 donor 序列构造反事实；它**不能自动消除未观察混杂**。
本项目的输出是拟合权重下的对照差距，不是投入有效性的因果证明。

## 空间 placebo 诊断

按同一协议，依次把每个 donor 当作伪处理厂区：

- 伪处理的 donor 池是其余 donor；**原受处理厂区永不进入任何 placebo
  donor 池**，也不作为 placebo 被评估。
- 对每个 pseudo 重新执行相同的训练/留出选 lambda、全 pre 重拟合。
- `pre_rmspe <= pre_rmspe_max` 的 placebo 才进入可比集合；门槛是唯一的
  可比性筛选条件。
- 主统计仍是 `post-MSPE / pre-MSPE`，零分母沿用上面的约定。
- 尾部秩比例为

  ```text
  (1 + #(可比 placebo 且 ratio >= treated ratio))
  / (1 + #(可比 placebo))
  ```

  并列（ratio 相等）计入分子；输出 `included_denominator`、
  `included_numerator`、`tail_rank_proportion`、`tied_placebos`，以及
  `excluded` 中每个被排除 donor 的原因（例如
  `pre_rmspe_exceeds_threshold`）。

若 treated 自己的 pre-RMSPE 超过门槛，则不报告可比秩，只输出失配诊断
（`rank.reportable = false`）。这是**空间 placebo 诊断**，不是真实随机
分配下的精确 p 值，`not_a_randomization_p_value` 恒为 `true`。

## 小型权重手算例

资料：[`data/hand-example.json`](../data/hand-example.json)。周 `0..5`，
`treatment_week = 4`，`validation_weeks = 1`，所以训练周为 `0..2`，留出周
为 `3`，最终重拟合用 pre 周 `0..3`。处理厂区恒为 2：

| week | treated | donor-1 | donor-2 |
|---:|---:|---:|---:|
| 0 | 2 | 1 | 3 |
| 1 | 2 | 2 | 2 |
| 2 | 2 | 3 | 1 |
| 3 | 2 | 4 | 0 |
| 4 (post) | 2 | 5 | -1 |
| 5 (post) | 2 | 6 | -2 |

约束为 `w1 + w2 = 1`、`w1,w2 >= 0`。第 0 周给出

```text
w1 * 1 + w2 * 3 = 2
w1 + w2 = 1
=> w2 = 1/2, w1 = 1/2.
```

验证第 1 周：`0.5*2 + 0.5*2 = 2`；第 2 周：`0.5*3 + 0.5*1 = 2`；第 3 周：
`0.5*4 + 0.5*0 = 2`。因此 pre-MSE 为 0。两个候选 lambda 的留出 MSE 都是
0，按同分规则选择较大的 `lambda = 0.01`。post 周合成值仍为 2，差距为 0，
主统计是 `zero_over_zero`。两个 pseudo placebo 的 pre-RMSPE 都是
`sqrt(6) ≈ 2.449`，低于门槛 3，因此秩分母为 `1 + 2 = 3`；本结果不是 p 值。

## 平衡面板样例

[`data/balanced-panel.json`](../data/balanced-panel.json) 是 12 周、5 个
厂区（1 个 treated + 4 个 donor）的平衡周面板，训练周 `0..5`、留出周
`6..7`、post 周 `8..11`。它用于展示完整的 `lambda_selection`、每周差距
曲线和 4 个 donor 的空间 placebo 轮换。所有厂区与数值均为合成。
