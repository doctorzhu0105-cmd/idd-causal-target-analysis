"""Allele-aligned LD matrix reader for 1000 Genomes EUR PLINK binary filesets.

Avoids a PLINK install: the .bed file is SNP-major with a fixed row size, so the
genotypes for a handful of cis variants can be read by seeking directly to the
rows named in the .bim.

PLINK1 .bed layout
------------------
bytes 0-2   magic 0x6c 0x1b, mode byte 0x01 (SNP-major)
then        ceil(n_samples / 4) bytes per variant, 2 bits per sample, low bits first

2-bit codes, in terms of the .bim A1 allele:
    00 -> 2 copies of A1
    01 -> missing
    10 -> 1 copy of A1 (heterozygous)
    11 -> 0 copies of A1
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# index by the 2-bit code -> dosage of the .bim A1 allele
_CODE_TO_A1_DOSAGE = np.array([2.0, np.nan, 1.0, 0.0], dtype=np.float64)

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def load_bim(bim_path: str) -> pd.DataFrame:
    """Read a .bim, adding the 0-based row index needed to seek into the .bed."""
    bim = pd.read_csv(
        bim_path,
        sep="\t",
        header=None,
        names=["chr", "snp", "cm", "pos", "a1", "a2"],
        dtype={"chr": str, "snp": str, "a1": str, "a2": str},
    )
    bim["bed_row"] = np.arange(len(bim))
    return bim


def read_dosages(bed_path: str, bed_rows, n_samples: int) -> np.ndarray:
    """Return an (n_variants, n_samples) matrix of .bim-A1 dosages."""
    bytes_per_variant = (n_samples + 3) // 4
    out = np.empty((len(bed_rows), n_samples), dtype=np.float64)

    with open(bed_path, "rb") as handle:
        magic = handle.read(3)
        if magic[:2] != b"\x6c\x1b":
            raise ValueError(f"{bed_path} is not a PLINK .bed file")
        if magic[2] != 1:
            raise ValueError("only SNP-major .bed files are supported")

        for i, row in enumerate(bed_rows):
            handle.seek(3 + int(row) * bytes_per_variant)
            raw = np.frombuffer(handle.read(bytes_per_variant), dtype=np.uint8)
            # unpack 4 samples per byte, low-order bit pair first
            codes = np.empty(bytes_per_variant * 4, dtype=np.uint8)
            codes[0::4] = raw & 0b11
            codes[1::4] = (raw >> 2) & 0b11
            codes[2::4] = (raw >> 4) & 0b11
            codes[3::4] = (raw >> 6) & 0b11
            out[i] = _CODE_TO_A1_DOSAGE[codes[:n_samples]]

    return out


def _alleles_match(effect: str, other: str, a1: str, a2: str):
    """Resolve orientation of (effect, other) against .bim (a1, a2).

    Returns +1 if the effect allele is a1, -1 if it is a2, and None if the
    variant cannot be matched even after strand flipping.
    """
    effect, other, a1, a2 = effect.upper(), other.upper(), a1.upper(), a2.upper()
    if {effect, other} == {a1, a2}:
        return 1 if effect == a1 else -1

    if effect in _COMPLEMENT and other in _COMPLEMENT:
        flipped_effect, flipped_other = _COMPLEMENT[effect], _COMPLEMENT[other]
        if {flipped_effect, flipped_other} == {a1, a2}:
            return 1 if flipped_effect == a1 else -1

    return None


def is_palindromic(effect: str, other: str) -> bool:
    return _COMPLEMENT.get(effect.upper()) == other.upper()


def ld_matrix(
    variants: pd.DataFrame,
    bim: pd.DataFrame,
    bed_path: str,
    n_samples: int,
    snp_col: str = "snp",
    effect_col: str = "effect_allele",
    other_col: str = "other_allele",
):
    """Correlation matrix for `variants`, signed to their own effect alleles.

    `variants` must carry rsID and effect/other alleles. Rows absent from the
    reference panel, or whose alleles cannot be reconciled with it, are dropped.
    Returns (kept_variants, correlation_matrix).
    """
    merged = variants.merge(
        bim[["snp", "a1", "a2", "bed_row"]], left_on=snp_col, right_on="snp", how="inner"
    )
    if merged.empty:
        return merged, np.zeros((0, 0))

    orientations = [
        _alleles_match(r[effect_col], r[other_col], r["a1"], r["a2"])
        for _, r in merged.iterrows()
    ]
    merged = merged.assign(ld_orientation=orientations)
    merged = merged[merged.ld_orientation.notna()].reset_index(drop=True)
    if merged.empty:
        return merged, np.zeros((0, 0))

    dosages = read_dosages(bed_path, merged.bed_row.to_numpy(), n_samples)

    # mean-impute missing calls so a few no-calls do not drop whole samples
    for i in range(dosages.shape[0]):
        missing = np.isnan(dosages[i])
        if missing.any():
            dosages[i, missing] = np.nanmean(dosages[i])

    # flip to each variant's own effect allele before correlating
    signs = merged.ld_orientation.to_numpy(dtype=float)
    dosages = dosages * signs[:, None] + (1 - signs)[:, None]

    rho = np.corrcoef(dosages)
    if rho.ndim == 0:
        rho = rho.reshape(1, 1)
    rho = np.nan_to_num(rho, nan=0.0)
    np.fill_diagonal(rho, 1.0)

    return merged, rho
