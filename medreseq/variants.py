"""Parse genomic variants from VCF or TSV files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int    # 1-based
    ref: str
    alt: str

    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}>{self.alt}"


def _parse_vcf(path: Path) -> Iterator[Variant]:
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            chrom, pos, _id, ref, alt_field = (
                parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            )
            alt = alt_field.split(',')[0]   # first allele for multi-allelic
            if alt in ('.', '*'):           # missing / spanning deletion
                continue
            yield Variant(chrom=chrom, pos=pos, ref=ref.upper(), alt=alt.upper())


def _parse_tsv(path: Path) -> Iterator[Variant]:
    with open(path) as fh:
        raw_header = fh.readline().strip()
        if not raw_header:
            return
        header = [c.lower().strip() for c in raw_header.split('\t')]

        chrom_col = next(
            (i for i, h in enumerate(header) if h in ('chrom', 'chr', 'chromosome')), None
        )
        pos_col = next(
            (i for i, h in enumerate(header) if h in ('pos', 'position', 'start')), None
        )
        ref_col = next(
            (i for i, h in enumerate(header) if h in ('ref', 'reference', 'ref_allele')), None
        )
        alt_col = next(
            (i for i, h in enumerate(header) if h in ('alt', 'alternate', 'alt_allele')), None
        )

        if any(c is None for c in (chrom_col, pos_col, ref_col, alt_col)):
            raise ValueError(
                f"TSV missing required columns (chrom/pos/ref/alt). Found columns: {header}"
            )

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(chrom_col, pos_col, ref_col, alt_col):  # type: ignore[arg-type]
                continue
            yield Variant(
                chrom=parts[chrom_col].strip(),         # type: ignore[index]
                pos=int(parts[pos_col]),                # type: ignore[index]
                ref=parts[ref_col].strip().upper(),     # type: ignore[index]
                alt=parts[alt_col].strip().upper(),     # type: ignore[index]
            )


def parse_variants(path: Path) -> Iterator[Variant]:
    """Auto-detect VCF or TSV format and yield Variant records.

    VCF is detected by .vcf extension or ##fileformat/# CHROM header lines.
    TSV must have a header row with columns: chrom/chr, pos, ref, alt
    (flexible synonyms accepted).
    """
    suffix = path.suffix.lower()
    if suffix == '.vcf':
        return _parse_vcf(path)
    if suffix in ('.tsv', '.txt', '.csv'):
        return _parse_tsv(path)
    # Sniff first line when extension is ambiguous
    with open(path) as fh:
        first = fh.readline()
    if first.startswith('##fileformat=VCF') or first.startswith('#CHROM'):
        return _parse_vcf(path)
    return _parse_tsv(path)
