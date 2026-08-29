# 椎间盘疾病相关 cis-MR 可复现代码

本仓库是科学 MASTER 1.5.0 候选对应的公开代码与可再分发派生结果层，不是稿件权威树。

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
python tools/verify_package.py --root .
```

外部输入、完整重算接口和边界见 `docs/DATA_ACCESS.md` 与 `docs/REPRODUCIBILITY.md`。
