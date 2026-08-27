# 椎间盘疾病相关 cis-MR 可复现代码

本仓库是与锁定科学 MASTER 1.4.0 同步的公开代码和可复现层。当前主线是配置依赖的 correlated cis-MR：signed LD、广义 IVW、广义 Q、暴露递增方向的 GLS MR-Egger、使用删减 Omega 且不重新剪枝的 LOO、有效维度/杠杆/联合信息诊断、多重检验族敏感性、同一 1000G EUR 面板的 LD bootstrap、LD shrinkage、coloc 先验敏感性，以及中性的 13 行证据矩阵与图件源数据。

边界如下：

- 仓库不是稿件权威树，不给出 A/B/C/D 层级或“已验证靶点”结论。
- 7,284 行全局 BH 是正式多重检验参照；13 行只是 K>=2 的下游证据敏感性子集。
- 200 次 LD bootstrap 只反映同一 503 人参考面板的抽样不确定性，不是独立面板复现。
- 已推迟已知真值模拟；不得写成已经完成。
- 不公开 eQTLGen/FinnGen/UKB 缓存、1000G PLINK 原始文件、个体级数据或再分发权利不明确的提取物。
- v1.0.0 旧模块保留用于历史复现，属于次级/历史范围。

离线检查：

```bash
python -B -m unittest discover -s tests -p "test_*.py"
python tools/verify_package.py --root .
```

外部输入和执行顺序见 `docs/DATA_ACCESS.md` 与 `docs/REPRODUCIBILITY.md`。当前只是未发布的 MASTER 1.4.0 同步候选，没有创建 tag、GitHub Release 或新的 Zenodo 版本。
