# 工业试点效果的合成对照与时间外验证

这是现场机器人应用分析的合成资料项目。`data`保存业务记录，
`analysis_data.io`负责严格JSON读取和原始字节指纹，`analysis_data.domain`
提供已有的基础汇总，`analysis_data.synthetic`提供合成对照、留出周选
lambda、每周差距曲线和空间placebo诊断。字段与数学口径见
[资料约定](docs/data-contract.md)，方法与手算例见
[合成对照说明](docs/synthetic-control.md)。所有厂区、设备、收益与结果均为
合成，不代表任何实际项目的效果。

Python 3.11及以上版本可直接运行，当前无第三方Python依赖。项目根目录中的
常用命令：

```sh
python -m unittest discover -s tests -v
python -m compileall -q analysis_data tests
python -m analysis_data.cli inspect data/normal.json
python -m analysis_data.cli analyze data/normal.json
python -m analysis_data.cli analyze data/hand-example.json
python -m analysis_data.cli analyze data/balanced-panel.json --output /tmp/scm-result.json
```

`inspect`展示原始资料摘要与基础统计。`analyze`在不覆盖输入的前提下输出
可复查的donor权重、lambda候选比较、KKT证据、每周观测/合成/差距、post
平均差距、空间placebo轮换与排除原因；输出保留输入文件指纹。合成对照
不能自动消除未观察混杂，空间placebo诊断也不是精确随机化p值。

边界资料与normal采用同样的JSON结构，invalid前缀文件保存不合法输入。
`data/hand-example.json`是小型权重手算例，`data/balanced-panel.json`是
12周、5厂区的平衡面板样例；二者均为合成资料。

容器可执行 `docker build -t field-analysis .` 后用
`docker run --rm field-analysis python -m analysis_data.cli analyze data/normal.json`
查看结果，或
`docker run --rm field-analysis python -m unittest discover -s tests -v`
运行检查。运行期间无外网或账号要求。
