# 图表复现代码

本目录仅包含生成文章各图所需的 Python 代码及依赖清单，不需要外部数据。

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

## 参数设置

- **图 1**：信号准确率 p=0.55、0.7、0.9；级联阈值 k=2、3、4、5；隐私预算在 [0,5] 上取 1,000 个等距点。
- **图 2**：信号准确率 p=0.55、0.7、0.9；隐私预算在 [0.05,4.5] 上取 3,000 个等距点。
  图 2(a) 计算报告数停止时间；图 2(b) 使用隐私容忍度 U[0,5]，上述预算范围内的参与率为 1-epsilon/5。
- **图 3**：信号准确率 p=0.55、0.7、0.9；等待成本 c=0.01；隐私容忍度 U[0,5]；
  隐私预算在 [0.1,4.95] 上取 3,000 个等距点。
- **图 4**：示意图的高斯曲线均值为 1、标准差为 1.8，决策阈值为 -1.6；
  曲线峰值缩放为 0.85，阶跃函数高度为 1。这些参数用于示意图排版，不是模拟估计值。
- **图 5**：两种状态的高斯信号均值为 -1、1，标准差 sigma=1；隐私预算 epsilon=1；
  决策阈值为 -3；信号网格为 [-5,5]，步长 0.01；翻转概率为 0.5*exp(-abs(s+3))。
- **图 6**：信念阈值 B=10，信号标准差 sigma=1，真实状态为 +1；隐私预算在 [0.1,2] 上取 25 个等距点。
  每个预算模拟 300 条路径，报告次数上限为 50,000，随机种子为 123；动态规划网格为 401 点，
  最大迭代次数为 100,000，收敛容差为 1e-8。图 6(b) 使用隐私容忍度 U[0,5]。
- **图 7**：五组模拟统一使用信号标准差 **sigma=1**、真实状态 +1、初始对数似然比 0；
  每组 5 条路径，每条 10,000 步。三组固定隐私预算为 0.1、0.5、1，异质预算服从 U[0,1]，
  另含无隐私基准；随机响应振幅为 1/2。异质组种子为 7–11，三组固定预算的种子分别为
  17–21、37–41、47–51，无隐私组种子为 57–61。每组高亮时间均值最接近组内总体均值的路径。
  `--quick` 使用 100 步、每组 2 条路径，其他参数不变，单独输出 `figure_07_quick`。
- **附录图 J.1**：信号标准差 sigma=1，初始对数似然比为 6，递推长度为 200；
  隐私预算为 0.1、0.5、1，无隐私基准使用 epsilon=100；随机响应振幅为 1/(1+exp(epsilon))。
  右侧面板在 epsilon=1 下比较精确递推与渐近递推。
