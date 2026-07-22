# GitHub 上传清单

1. 新建空仓库，建议名称为 `idd-causal-target-analysis`。
2. 上传本目录内部的全部文件，不上传外层日期目录或 ZIP 本身。
3. 在仓库根目录运行 `python tools/verify_package.py --root .`，结果必须为 `PASS`。
4. 确认 `LICENSE`、`manifest.tsv`、`config/`、`environment/`、`py/`、`R/` 和 `tests/` 均存在。
5. 不额外上传本机原始数据、处理结果、日志、论文、第三方软件或 API 缓存。
6. 仓库真实创建后记录 URL；需要论文引用时，再固定 release tag 并通过实际归档生成 DOI。
7. 将真实 URL 和 DOI 发回后，再更新文章的数据与代码可得性声明。
