# 厂区周面板

weeks连续整数，treatment_week及以后为post。pre最后validation_weeks个连续周是留出，之前为训练，两段均非空；同分以均方误差绝对差<=配置tie_tolerance判定，再选较大lambda。默认tie_tolerance=1e-12。不做z-score或处理后数据标准化。donor-a/c故意共线。

单纯形约束w>=0、sum w=1。目标为均方误差+lambda||w||²，KKT残差应基于此确切尺度。报告梯度、等式乘子、非负约束余量及互补残差，或等价投影梯度证据并有独立检查方法。pre-MSPE为全部pre周重拟合后的均方误差。placebo筛选仅看pre-RMSPE<=预设门槛。比值pre=post=0时定义0，pre=0且post>0定义+infinity（JSON用明确字符串，不写NaN/Infinity非法常量）。秩比例=(1+通过门槛且ratio>=treated ratio的placebo数)/(1+通过门槛placebo数)，如果treated自己pre超门槛则不报告可比秩，只给失配诊断。此比例不是随机化p值。

domain只验证平衡面板及算各厂区全期均值，没有拟合或对照诊断。

## 合成分析资料字段

输入至少包含：

- `weeks`：连续整数列表，非空。
- `sites`：非空、非重复厂区字符串列表。
- `rows`：每项 `{site, week, value}`；必须覆盖每个厂区/周的笛卡尔积，
  不允许缺失、重复身份或非有限值（布尔值也不接受）。
- `treated`：`sites` 中的唯一受处理厂区；其余厂区视为始终未处理的 donor。
- `treatment_week`：面板中的一周；之前为 pre，之后（含该周）为 post，
  两段均非空。
- `validation_weeks`：正整数，小于 pre 长度；pre 末尾该数量的连续周为
  留出段。
- `lambdas`：非空、有限、非负候选列表。
- `pre_rmspe_max`：有限、非负的 placebo 可比门槛。
- `solver`：可选；`tolerance` 默认 `1e-8` 且必须为正，
  `max_iterations` 默认 `20000` 且必须为正整数。
- `tie_tolerance`：可选，默认 `1e-12`，非负。
- `zero_tolerance`：可选，默认 `1e-12`，非负；用于按同一约定判断
  “pre=0”和“post=0”。

处理后的资料不参与标准化、lambda 选择或权重拟合。修改任何 post 数值只
会改变 post 合成值、差距和 post-MSPE，不会改变权重或选中的 lambda。

## 输出口径

`analysis_data.synthetic.analyze(value, source_sha256=None)` 返回：

- `protocol`：treated、donors、pre/training/holdout/post 周、lambda 候选、
  容差、求解配置，以及“无截距、无标准化、单纯形约束”的声明。
- `lambda_selection.candidates`：每个 lambda 的训练 MSE、留出 MSE、目标、
  KKT 残差、终止原因和训练段权重。
- `fit`：选中 lambda、全部 pre 周重拟合权重、目标、梯度、等式乘子、
  `mu`、KKT 各项残差、终止原因、pre/post MSPE、post 平均差距和主统计。
- `weekly`：每周 `observed`、`synthetic`、`gap`、`period`。
- `diagnostics`：pre 拟合是否在门槛内、失配提示、非因果声明。
- `spatial_placebo`：每个 donor 的伪处理结果、donor 池、排除原因、
  入选分母/分子、并列计入的尾部秩比例和“非随机化 p 值”声明。

预算耗尽（`termination_reason = "budget_exhausted"`）或停滞
（`"stalled_without_dual_violation"`）时，`optimal` 必须为 `false`。

零分母约定：pre-MSPE 与 post-MSPE 都为 0 时比率定义为 0；pre-MSPE 为 0
且 post-MSPE 为正时，JSON 中比率写字符串 `"+infinity"`，不写非法常量。
