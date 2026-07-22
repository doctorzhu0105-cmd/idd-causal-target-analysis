#!/usr/bin/env python3
"""Download the registered GSE230809 single-cell inputs with resume support."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import requests


SMALL_INPUTS = {
    "GSE230809_family.soft.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230809/soft/GSE230809_family.soft.gz"
    ),
    "GSE229711_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE229nnn/GSE229711/matrix/GSE229711_series_matrix.txt.gz"
    ),
    "GSE230808_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230808/matrix/GSE230808_series_matrix.txt.gz"
    ),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path, expected_bytes: int | None = None) -> None:
    if path.exists() and (expected_bytes is None or path.stat().st_size == expected_bytes):
        return
    partial = path.with_suffix(path.suffix + ".part")
    start = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={start}-"} if start else {}
    with requests.get(url, headers=headers, stream=True, timeout=(30, 120)) as response:
        response.raise_for_status()
        if start and response.status_code != 206:
            start = 0
            partial.unlink(missing_ok=True)
        mode = "ab" if start else "wb"
        with partial.open(mode) as handle:
            for block in response.iter_content(chunk_size=8 * 1024 * 1024):
                if block:
                    handle.write(block)
    if expected_bytes is not None and partial.stat().st_size != expected_bytes:
        raise RuntimeError(
            f"Incomplete download for {path.name}: {partial.stat().st_size} != {expected_bytes}"
        )
    partial.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.repo_root.resolve()
    sources = json.loads((root / "config" / "idd_sources.json").read_text(encoding="utf-8"))
    contract = sources["single_cell"]
    output_dir = root / "data" / "raw" / "idd" / "GSE230809"
    output_dir.mkdir(parents=True, exist_ok=True)
    for name, url in SMALL_INPUTS.items():
        download(url, output_dir / name)
    raw_tar = output_dir / "GSE230809_RAW.tar"
    download(contract["url"], raw_tar, int(contract["expected_content_length"]))
    manifest = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and not path.name.endswith(".part"):
            manifest.append(
                {
                    "path": path.relative_to(root).as_posix(),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    log_path = root / "logs" / "idd" / "scrna_source_manifest.json"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(
        json.dumps({"files": manifest}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({"file_count": len(manifest), "raw_tar_bytes": raw_tar.stat().st_size}))


if __name__ == "__main__":
    main()
