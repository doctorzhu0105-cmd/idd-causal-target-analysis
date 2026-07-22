#!/usr/bin/env python3
"""Render the IDD target dossier from frozen-style tabular artifacts."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path, default: Mapping[str, object] | None = None) -> dict[str, object]:
    if not path.exists():
        return dict(default or {})
    return json.loads(path.read_text(encoding="utf-8"))


def fmt(value: object, digits: int = 6) -> str:
    if value is None or str(value).strip() in {"", "nan", "NA", "not_estimable"}:
        return "not_estimable"
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return str(value)


def yes_no(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return "是"
    if text in {"false", "0", "no"}:
        return "否"
    return "未估计"


def best_gtex(rows: Iterable[Mapping[str, str]], gene: str) -> Mapping[str, str] | None:
    candidates = [
        row
        for row in rows
        if row.get("gene_symbol") == gene and row.get("status") == "estimable"
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda row: float(row.get("pp4") or "-1"))


def target_section(
    row: Mapping[str, str],
    gtex_rows: list[dict[str, str]],
    docking: Mapping[str, str],
    mr_rows: list[dict[str, str]],
    coloc_rows: list[dict[str, str]],
    smr_rows: list[dict[str, str]],
    reverse_rows: list[dict[str, str]],
    sensitivity_rows: list[dict[str, str]],
    loo_rows: list[dict[str, str]],
) -> str:
    gene = row["gene_symbol"]
    gtex = best_gtex(gtex_rows, gene)
    gtex_text = (
        f"{gtex.get('tissue')} / {gtex.get('outcome_id')} / PP4={fmt(gtex.get('pp4'))}"
        if gtex
        else "not_estimable"
    )
    docking_row = next((item for item in docking.get("rows", []) if item.get("gene_symbol") == gene), None)
    docking_text = docking_row.get("reason", "not_estimable") if docking_row else "not_estimable"
    outcomes = []
    outcome_ids = sorted(
        {
            item.get("outcome_id", "")
            for item in mr_rows + coloc_rows + smr_rows + reverse_rows
            if item.get("gene_symbol") == gene and item.get("outcome_id")
        }
    )
    for outcome_id in outcome_ids:
        mr = next(
            (item for item in mr_rows if item.get("gene_symbol") == gene and item.get("outcome_id") == outcome_id),
            {},
        )
        coloc = next(
            (item for item in coloc_rows if item.get("gene_symbol") == gene and item.get("outcome_id") == outcome_id),
            {},
        )
        smr = next(
            (item for item in smr_rows if item.get("gene_symbol") == gene and item.get("outcome_id") == outcome_id),
            {},
        )
        reverse = next(
            (item for item in reverse_rows if item.get("gene_symbol") == gene and item.get("outcome_id") == outcome_id),
            {},
        )
        outcomes.append(
            f"  - `{outcome_id}`：MR beta={fmt(mr.get('beta'))}，P={fmt(mr.get('pvalue'))}，FDR={fmt(mr.get('fdr'))}；coloc PP4={fmt(coloc.get('pp4'))}；SMR P={fmt(smr.get('p_smr'))}，HEIDI P={fmt(smr.get('p_heidi'))}；反向 MR FDR={fmt(reverse.get('fdr')) if reverse.get('status') == 'estimable' else 'not_estimable'}。"
        )
    best_outcome = row.get("best_outcome_id", "")
    egger = next(
        (
            item
            for item in sensitivity_rows
            if item.get("gene_symbol") == gene
            and item.get("outcome_id") == best_outcome
            and item.get("direction") == "expression_to_disease"
            and item.get("method") == "mr_egger"
        ),
        {},
    )
    gene_loo = [
        item
        for item in loo_rows
        if item.get("gene_symbol") == gene
        and item.get("outcome_id") == best_outcome
        and item.get("direction") == "expression_to_disease"
        and item.get("status") == "estimable"
    ]
    worst_loo = max(gene_loo, key=lambda item: float(item.get("pvalue") or "-1"), default={})
    sensitivity_text = (
        f"Cochran Q P={fmt(egger.get('q_pvalue'))}；Egger intercept P={fmt(egger.get('egger_intercept_pvalue'))}，FDR={fmt(egger.get('egger_intercept_fdr'))}；"
        f"LOO 最弱一行排除 {worst_loo.get('excluded_snp', 'not_estimable')} 后 P={fmt(worst_loo.get('pvalue'))}。"
        if egger.get("status") == "estimable" or worst_loo
        else "单工具或工具不足，weighted median、MR-Egger 与 LOO 不可估计；SMR 与 Wald ratio 不能视为独立工具重复。"
    )
    return "\n".join(
        [
            f"### {row['priority_role']}：{gene}",
            f"- 综合等级：{row.get('evidence_tier', 'not_estimable')}；排序分数：{fmt(row.get('score'))}；证据模态数：{row.get('evidence_modalities', 'not_estimable')}。",
            f"- 组织 DEG：{yes_no(row.get('bulk_deg_pass'))}；WGCNA 选定模块：{yes_no(row.get('wgcna_module_pass'))}。该组优先靶点没有同时通过这两项组织发现证据，不能写成组织-遗传闭环。",
            f"- 最佳结局：`{row.get('best_outcome_id', 'not_estimable')}`；MR beta={fmt(row.get('mr_beta'))}，95% CI [{fmt(row.get('mr_ci_lower'))}, {fmt(row.get('mr_ci_upper'))}]；OR={fmt(row.get('mr_odds_ratio'))}，OR 95% CI [{fmt(row.get('mr_or_ci_lower'))}, {fmt(row.get('mr_or_ci_upper'))}]；P={fmt(row.get('mr_pvalue'))}；正式集合 FDR={fmt(row.get('mr_formal_set_fdr'))}；工具数={row.get('mr_n_instruments', 'not_estimable')}。",
            "- 三结局逐行结果（用于避免把最佳结局推广到直接 IDD）：",
            *outcomes,
            f"- 最佳结局敏感性：{sensitivity_text}",
            f"- coloc：共享变异={row.get('coloc_shared_variants', 'not_estimable')}；PP4={fmt(row.get('coloc_pp4'))}。SMR beta={fmt(row.get('smr_beta'))}，P={fmt(row.get('smr_pvalue'))}，FDR={fmt(row.get('smr_fdr'))}；HEIDI P={fmt(row.get('heidi_pvalue'))}，SNP={row.get('heidi_n_snp', 'not_estimable')}。",
            f"- Steiger 方向支持：{yes_no(row.get('steiger_pass'))}；最佳结局反向 MR：`{row.get('reverse_mr_status', 'not_estimable')}`。跨结局反向 MR 信号：{yes_no(row.get('reverse_mr_conflict'))}；信号结局={row.get('reverse_mr_signal_outcomes', 'not_estimable')}；最小 FDR={fmt(row.get('reverse_mr_min_fdr'))}。反向 MR 不可估计不作为阴性证据；若存在信号，需将该靶点标为因果方向冲突。上述遗传结论是 forward-integrated（MR+coloc+SMR/HEIDI），不是完整双向 MR 通过。",
            f"- eQTL：主暴露为 `{row.get('eqtl_source', 'not_estimable')}`；GTEx Catalogue 代理组织最佳结果：{gtex_text}。GTEx 代理组织不是椎间盘组织 eQTL。",
            f"- 单细胞 marker-score 宽泛定位：{row.get('scrna_cell_type', 'not_detected')} / {row.get('scrna_tissue', 'not_estimable')}；表达细胞比例={fmt(row.get('scrna_cell_fraction_expressing'))}；供体数={row.get('scrna_donors', 'not_estimable')}。本数据只作定位，不能作疾病差异效应或高置信细胞注释。",
            f"- MR 推断干预方向：`{row.get('desired_intervention', 'not_estimable')}`。数据库最高调节剂：{row.get('top_modulator_name', 'not_estimable')}（{row.get('top_modulator_source', 'not_estimable')}），方向匹配={yes_no(row.get('top_modulator_direction_aligned'))}。",
            f"- 推荐水凝胶/纳米粒载荷：`{row.get('recommended_local_payload', 'not_estimable')}`。这只是局部制剂可行性假设，不等于盘内滞留、药效或安全性已证明。",
            f"- Docking：`{docking.get('status', 'not_estimable')}`；该靶点判定：{docking_text}。未报告虚构的 Vina 能量。",
        ]
    )


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    priority = read_tsv(root / "artifacts" / "targets" / "priority_targets.tsv")
    gtex_rows = read_tsv(root / "artifacts" / "mr" / "gtex_coloc_results.tsv")
    mr_rows = read_tsv(root / "artifacts" / "mr" / "mr_results.tsv")
    coloc_rows = read_tsv(root / "artifacts" / "mr" / "coloc_results.tsv")
    smr_rows = read_tsv(root / "artifacts" / "mr" / "smr_heidi_results.tsv")
    reverse_rows = read_tsv(root / "artifacts" / "mr" / "reverse_mr_results.tsv")
    sensitivity_rows = read_tsv(root / "artifacts" / "mr" / "mr_sensitivity_results.tsv")
    loo_rows = read_tsv(root / "artifacts" / "mr" / "leave_one_out_results.tsv")
    cluster_annotations = read_tsv(root / "artifacts" / "scrna" / "gse230809_cluster_annotations.tsv")
    docking_rows = read_tsv(root / "artifacts" / "immune_drug" / "docking_assessment.tsv")
    docking_summary = read_json(root / "artifacts" / "immune_drug" / "docking_summary.json")
    docking = {"status": docking_summary.get("status", "not_estimable"), "rows": docking_rows}
    bulk = read_json(root / "logs" / "idd" / "bulk_summary.json")
    gate0 = read_json(root / "logs" / "idd" / "gate0_summary.json")
    scrna = read_json(root / "logs" / "idd" / "scrna_summary.json")
    gtex_summary = read_json(root / "artifacts" / "mr" / "gtex_coloc_summary.json")
    drug = read_json(root / "logs" / "idd" / "drug_evidence_summary.json")
    ranking = read_json(root / "artifacts" / "targets" / "ranking_summary.json")
    analysis_config = read_json(root / "config" / "idd_analysis.json")
    generated_at = str(analysis_config.get("registered_at", "not_estimable"))
    lines = [
        "# IDD/LDH 可成药与可递送因果靶点档案",
        "",
        f"> 生成基准时间：{generated_at}。本档案由已落盘 artifact 自动编排；缺失证据保留为 `not_estimable`。",
        "",
        "## 结论边界",
        "",
        "当前路线为 **GO**：指定 bulk、单细胞与 FinnGen 结局均可获取，且 GSE70362 AF 在预先固定阈值下存在 FDR 控制的真实差异信号。输出目标是后续水凝胶/纳米粒湿实验的候选依据，不是独立新颖性论文，也不声称首创。",
        "",
        "排序主结果为：`NT5E` 主靶；`AKR1C2`、`CCL28`、`P2RY13` 为三个不同位点的备选。`AKR1C1` 与 `AKR1C2` 位于同一 chr10 近邻位点，保留为 locus alternative，不占用独立备选名额。四个优先靶点均未同时通过组织 DEG 与 WGCNA 模块证据，因此不能写成“组织+遗传+单细胞”完整闭环；其核心支持来自 forward MR、coloc、SMR/HEIDI、Steiger 与椎间盘细胞定位的组合。`NT5E` 另有一个可估计的反向 MR 信号，主靶下注前必须做方向冲突复核，具体缺项见下表。",
        "",
        "## 数据与门槛",
        "",
        f"- GSE70362：AF n={gate0.get('primary', {}).get('AF', {}).get('samples', 'not_estimable')}，FDR 特征={gate0.get('primary', {}).get('AF', {}).get('fdr', 'not_estimable')}，FDR+|log2FC| 特征={gate0.get('primary', {}).get('AF', {}).get('fdr_effect', 'not_estimable')}；NP FDR 特征={gate0.get('primary', {}).get('NP', {}).get('fdr', 'not_estimable')}。WGCNA AF power={bulk.get('wgcna_power', [{}])[0].get('power', 'not_estimable')}，选定非灰模块={bulk.get('selected_af_modules', 'not_estimable')}。",
        f"- GSE230809：样本={scrna.get('samples', 'not_estimable')}，QC 前细胞={scrna.get('cells_before_qc', 'not_estimable')}，QC 后={scrna.get('cells_after_qc', 'not_estimable')}，聚类子集={scrna.get('clustered_cells', 'not_estimable')}，cluster={scrna.get('clusters', 'not_estimable')}；定位目标检出={scrna.get('targets_detected', 'not_estimable')}/{scrna.get('targets_requested', 'not_estimable')}。cluster 宽泛注释计数：{', '.join(f'{key}={sum(item.get("broad_cell_type") == key for item in cluster_annotations)}' for key in sorted({item.get('broad_cell_type', '') for item in cluster_annotations if item.get('broad_cell_type')})) or 'not_estimable'}。健康与疾病供体年龄结构不同，因此不作疾病差异效应。",
        f"- GTEx Catalogue：共扫描 {gtex_summary.get('target_count', 'not_estimable')} 个目标（含 1 个 locus alternative）× {gtex_summary.get('dataset_count', 'not_estimable')} 个代理组织 × 3 个结局；{gtex_summary.get('coloc_estimable', 'not_estimable')}/{int(gtex_summary.get('target_count', 0))*int(gtex_summary.get('dataset_count', 0))*3} 个组合可估计，PP4>=0.8 为 {gtex_summary.get('coloc_pp4_at_least_0_8', 'not_estimable')}；代理组织只作补充。",
        f"- 药物数据库：候选调节剂快照由 Open Targets、DGIdb、ChEMBL 生成；已知调节剂候选数={drug.get('known_modulator_candidates', 'not_estimable')}。递送判断仅是局部制剂可行性，不是疗效证据。",
        "- pQTL（deCODE/UKB-PPP）按用户要求延后且不阻断；doublet detection 因 Windows 依赖缺失未运行；MR-PRESSO 未运行；反向 MR 的不可估计项不当作阴性结果。",
        "- 本轮按硬规则未进行新的文献检索，因此不列未经核实的前人引用，也不作“首个/首次”表述。",
        "",
        "## 优先级总表",
        "",
        "| 角色 | 基因 | 等级 | 最佳结局 | MR beta (95% CI) | OR | coloc PP4 | SMR P / HEIDI P | 单细胞定位 | 干预方向 | 推荐载荷 |",
        "|---|---|---|---|---|---:|---:|---|---|---|---|",
    ]
    for row in priority:
        lines.append(
            f"| {row.get('priority_role')} | `{row.get('gene_symbol')}` | {row.get('evidence_tier')} | {row.get('best_outcome_id')} | {fmt(row.get('mr_beta'))} [{fmt(row.get('mr_ci_lower'))}, {fmt(row.get('mr_ci_upper'))}] | {fmt(row.get('mr_odds_ratio'))} | {fmt(row.get('coloc_pp4'))} | {fmt(row.get('smr_pvalue'))} / {fmt(row.get('heidi_pvalue'))} | {row.get('scrna_cell_type')} ({row.get('scrna_tissue')}, {fmt(row.get('scrna_cell_fraction_expressing'))}) | {row.get('desired_intervention')} | `{row.get('recommended_local_payload')}` |"
        )
    lines.extend(["", "## 靶点档案", ""])
    for row in priority:
        lines.extend(
            [
                target_section(
                    row,
                    gtex_rows,
                    docking,
                    mr_rows,
                    coloc_rows,
                    smr_rows,
                    reverse_rows,
                    sensitivity_rows,
                    loo_rows,
                ),
                "",
            ]
        )
    lines.extend(
        [
            "## 同位点替代与未入选候选",
            "",
            "- `AKR1C1`：与 `AKR1C2` 同一 chr10 近邻位点，MR/coloc/SMR/HEIDI 证据可用，但不作为独立备选；其胞内表达增强载荷与 AKR1C2 具有竞争关系，湿实验应先做选择性和位点拆分。",
            "- 组织支持但无完整遗传支持的候选（例如 `HYAL1`、`MAP4K4`、`MGST1`）只能作为 exploratory tissue-supported targets，不进入本次因果优先级。",
            "",
            "## Docking 判定",
            "",
            f"结构门结果为 **{docking_summary.get('status', 'not_estimable')}**：{docking_summary.get('reason', 'not_estimable')}。因此没有生成 Vina affinity、没有把预测结构当作实验结构，也没有用方向相反的抑制剂 docking 冒充治疗支持。",
            "",
            "## 湿实验优先预测",
            "",
            "1. 在退变刺激的 NP/AF 细胞体系中，用 NT5E mRNA 纳米粒或可验证的重组蛋白局部递送，先检测 CD73/腺苷轴、基质合成与分解指标；这是 MR 保护方向的可检验预测，不是已证实疗效。",
            "2. 对 AKR1C2（并以 AKR1C1 作同位点替代）做表达增强而非抑制剂验证，先完成同工酶选择性、细胞内递送和剂量-反应分离。",
            "3. CCL28 与 P2RY13 分别按抑制方向设计 siRNA/局部载荷；P2RY13 在当前单细胞中接近检测下限，必须先用独立 qPCR/IHC/蛋白检测确认椎间盘表达，再决定是否下注。",
            "",
            "## 复现与冻结",
            "",
            "默认入口为 `make all`，旧骨质疏松路线只保留为显式 legacy targets。完成本档案后由 `py/idd_freeze_core.py` 对白名单 artifact 做字节数和 SHA-256 清单冻结；冻结目录只允许 `--verify-only` 校验，不允许下游脚本回写。",
            "",
        ]
    )
    output = root / "HANDOFF" / "IDD_TARGET_DOSSIER.md"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
