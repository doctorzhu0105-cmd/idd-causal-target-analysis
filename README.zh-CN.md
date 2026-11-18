# 椎间盘疾病相关 cis-MR 可复现代码

本仓库保留公开的相关 cis-MR 代码与可再分发结果。MASTER 1.5.0 是 v1.2.0 目录分析的历史内部标签，不是 R1 修订版发布号；[v1.3.0](https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases/tag/v1.3.0) 已于 2026 年 11 月 17 日发布，版本 DOI 为 [10.5281/zenodo.24312891](https://doi.org/10.5281/zenodo.24312891)，R1 系列 concept DOI 为 [10.5281/zenodo.24312890](https://doi.org/10.5281/zenodo.24312890)。v1.3.1 为维护性更新，仅调整核验程序名称、命令行入口和发布文档，估计方法和科研结果不变。发布日期以 [GitHub Releases](https://github.com/doctorzhu0105-cmd/idd-causal-target-analysis/releases) 为准；引用时请在上述 R1 系列中选择对应版本的 Zenodo 记录，使用该记录给出的版本 DOI。

R1 软件将 6,087 条历史完整矩阵比较与 1,068 条依赖 C3 COMPUTED_UNVERIFIED 输入的条件性标量重建分开，合计 7,155 / 7,284。其余 129 条不完整，另有 150 个成员信息缺失的配置，二者不是同一分母。1,068 条逐键派生数据的再分发条件尚未核清，未随软件分发。四个确定性数值核验程序是既有结果的独立核验器，不是原始分析生成器。详细范围见 revision_r1/README.md。

v1.2.0 在 v1.1.0 correlated cis-MR 复现层基础上增加：

- 2,479 个基因、5 个 LD 阈值的 Effective K 与联合暴露信息；
- 3 个 FinnGen 结局的 5 阈值广义 IVW 梯度；
- 暴露信息、精度和广义杠杆的配对变化；
- 匹配参数化及分结局的信息-精度关联；
- AKR1C1/AKR1C2 在全目录中的位置；
- K=1 二阶标准误敏感性和全局 BH/BY 排名传播；
- 修订 Figure 1 的确定性源表与生成脚本。

这些结果是配置依赖的描述性分析，不是新估计器、估计器验证、独立复现或治疗靶点验证。

仓库不包含 eQTLGen、FinnGen、1000 Genomes、UKB 的原始文件、缓存或个体数据。公开的逐基因表仅含模型派生结果，不含 SNP、等位基因、坐标、源数据效应量、EAF 或参与者记录。

离线验证：

```bash
python -B -m unittest discover -s tests -p "test_*.py"
python -B -m unittest discover -s revision_r1/tests -p "test_*.py"
python -B tools/verify_package.py --root .
```

外部输入、完整重算接口和边界见 `docs/DATA_ACCESS.md` 与 `docs/REPRODUCIBILITY.md`。
