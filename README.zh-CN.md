# 椎间盘疾病因果靶点分析代码

本目录是用于 GitHub 公开的代码包，覆盖 bulk DEG/WGCNA、druggable-genome cis-eQTL MR、coloc、SMR/HEIDI、Steiger、区域方向分析、单细胞定位、药物证据、docking 适用性、靶点排序和图表派生。

本包不包含论文、作者联系方式、投稿文件、原始或处理后源数据、API 缓存、第三方可执行文件以及当前数值结果。完整运行前需按 `docs/DATA_ACCESS.md` 从原始提供方获取数据，并遵守相应条款。

上传前运行：

```bash
python -B -m unittest discover -s tests -p "test_*.py"
python tools/verify_package.py --root .
```

代码采用 MIT License。真实 GitHub URL、release tag 和归档 DOI 只能在仓库实际发布和归档后填写，当前包不会预先编造这些信息。
