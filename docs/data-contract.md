# 厂区周面板

weeks连续整数，treatment_week及以后为post。pre最后validation_weeks个连续周是留出，之前为训练，两段均非空；同分以均方误差绝对差<=配置tie_tolerance判定，再选较大lambda。默认tie_tolerance=1e-12。不做z-score或处理后数据标准化。donor-a/c故意共线。

单纯形约束w>=0、sum w=1。目标为均方误差+lambda||w||²，KKT残差应基于此确切尺度。报告梯度、等式乘子、非负约束余量及互补残差，或等价投影梯度证据并有独立检查方法。pre-MSPE为全部pre周重拟合后的均方误差。placebo筛选仅看pre-RMSPE<=预设门槛。比值pre=post=0时定义0，pre=0且post>0定义+infinity（JSON用明确字符串，不写NaN/Infinity非法常量）。秩比例=(1+通过门槛且ratio>=treated ratio的placebo数)/(1+通过门槛placebo数)，如果treated自己pre超门槛则不报告可比秩，只给失配诊断。此比例不是随机化p值。

domain只验证平衡面板及算各厂区全期均值，没有拟合或对照诊断。
