# 图表复现代码

本目录仅包含生成文章各图所需的 Python 代码及依赖清单。无需安装 R，也不需要外部数据。

## 图号与代码对应表

| 文章图号 | Python 文件 | 图片内容 | 输出文件（位于 figures/） |
| --- | --- | --- | --- |
| 图 1 | `figure_01_binary_accuracy.py` | 二元信号：正确最终决策概率与隐私预算 | `figure_01.pdf/png` |
| 图 2(a)、2(b) | `figure_02_binary_stopping_time.py` | 二元信号：报告数停止时间、日历停止时间 | `figure_02a.pdf/png`、`figure_02b.pdf/png` |
| 图 3 | `figure_03_binary_objective.py` | 二元信号：三个信号准确率下的平台目标函数 | `figure_03.pdf/png` |
| 图 4 | `figure_04_signal_threshold.py` | 决策阈值附近的阶跃决策函数与高斯信号分布示意 | `figure_04.pdf/png` |
| 图 5 | `figure_05_smooth_response.py` | 信号密度、平滑随机响应和行动翻转概率三个面板 | `figure_05.pdf/png` |
| 图 6(a) | `figure_06a_continuous_report_time.py` | 连续信号：报告数停止时间，比较动态规划、模拟与渐近近似 | `figure_06a.pdf/png` |
| 图 6(b) | `figure_06b_continuous_calendar_time.py` | 连续信号：纳入参与率后的日历停止时间 | `figure_06b.pdf/png` |
| 图 7 | `figure_07_privacy_regimes.py` | 五种隐私机制下的对数似然比轨迹 | `figure_07.pdf/png` |
| 附录图 J.1 | `figure_J1_likelihood_evolution.py` | 正报告递推轨迹及精确递推与渐近递推的比较 | `figure_J1.pdf/png` |

## 安装与运行

已验证环境：Python 3.13；依赖版本见 `requirements.txt`。建议使用独立虚拟环境。

```sh
python -m pip install -r requirements.txt
python figure_01_binary_accuracy.py
python figure_02_binary_stopping_time.py
python figure_03_binary_objective.py
python figure_04_signal_threshold.py
python figure_05_smooth_response.py
python figure_06a_continuous_report_time.py
python figure_06b_continuous_calendar_time.py
python figure_07_privacy_regimes.py
python figure_J1_likelihood_evolution.py
```

每个文件可独立运行，输出目录固定为脚本所在目录下的 `figures/`。程序自动保存 PDF 和 PNG，
不弹出绘图窗口。图 2、图 6 的两个面板分别输出；图 5 和附录图 J.1 各自输出合并面板。
压缩包和仓库仅收录代码，不包含预生成图片。

## 参数与转换说明

- 图 1–3 的信号准确率为 0.55、0.7、0.9；图 3 的等待成本为 0.01，隐私容忍度服从 U[0,5]。
- 图 6 的默认参数为 B=10、sigma=1，25 个隐私预算取值、300 条模拟路径、50,000 次报告上限。
  未停止路径不会被丢弃后用于计算无条件均值；程序会给出警告及停止比例。
- 图 7 默认 10,000 步，每组 5 条路径，sigma=1；固定隐私预算为 0.1、0.5、1，
  异质预算服从 U[0,1]，另含无隐私基准。保留原程序的随机数序列、路径选取和高亮规则。
  `python figure_07_privacy_regimes.py --quick` 仅用于小规模检查，单独输出 `figure_07_quick`。
- 图 1、图 5、图 7、附录图 J.1 从原绘图计算转换为 Python。逐点数值比较中，
  图 1、5 最大绝对误差小于 1e-14；图 7 全部 25 条、每条 10,000 步的轨迹误差小于 2e-13；
  附录图 J.1 的误差小于 3e-9。字体及图形库的渲染细节可能不同。
- 图 4 未找到原独立源码，现按文章图示重绘，保留阶跃函数、高斯曲线、阈值与邻近信号标记；
  它是示意图，曲线高度和横向位置用于排版，不应从图中读取估计参数。
- 附录图 J.1 保留原图使用的振幅 `1/(1+exp(epsilon))`，正文图 5–7 使用 `1/2`。
  附录图展示指定正报告历史下的确定性递推，不代替完整随机学习过程。

所有默认绘图入口均已实际运行。图 7 的随机序列兼容逻辑直接由 Python 执行，运行时不调用 R。
