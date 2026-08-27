# 同步候选复核清单

本文件只描述未来作者批准后的操作；Stage 7A 不执行 push、tag、GitHub Release 或 Zenodo 更新。

1. 在 clean clone 的候选分支运行 `python -B -m unittest discover -s tests -p "test_*.py"`。
2. 运行 `python tools/verify_package.py --root .`，结果必须为 `PASS`。
3. 核对 `RELEASE_STATUS.json` 仍为 `UNRELEASED_MASTER_1.4.0_SYNC_CANDIDATE`，且三个公开版本字段均为 `null`。
4. 核对 `provenance/`、`results/master_1_4/`、当前 correlated cis-MR 代码、环境锁和 manifest 均在候选 commit 内。
5. 确认没有原始/受限数据、UKB cache、1000G PLINK 二进制、作者信息、稿件治理文件或凭据。
6. 仅在 GPT/作者明确批准后，将该普通前向 commit 推送至现有公开仓库；不得改写 `main` 历史或删除 `v1.0.0`。
7. 后续 release/tag/Zenodo 版本号须另行决定。本候选未预选版本号，也不得把 MASTER `1.4.0` 自动当作公开 release 号。
