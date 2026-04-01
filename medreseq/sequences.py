"""Reference genome sequence retrieval via pysam/htslib."""

from __future__ import annotations

from pathlib import Path

import pysam

from .variants import Variant


class ReferenceGenome:
    """Thin RAII wrapper around a FASTA file with .fai index."""

    def __init__(self, fasta_path: Path) -> None:
        self.path = Path(fasta_path)
        if not self.path.exists():
            raise FileNotFoundError(f"FASTA not found: {self.path}")
        fai = Path(str(self.path) + '.fai')
        if not fai.exists():
            raise FileNotFoundError(
                f"FASTA index not found: {fai}\n"
                f"Create it with: samtools faidx {self.path}"
            )
        self._fasta = pysam.FastaFile(str(self.path))

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch sequence at 0-based half-open [start, end). Returns uppercase."""
        length = self._fasta.get_reference_length(chrom)
        return self._fasta.fetch(chrom, max(0, start), min(end, length)).upper()

    def fetch_around(self, variant: Variant, flank: int) -> tuple[str, int]:
        """Fetch `flank` bp on each side of the variant.

        Returns (sequence, target_offset) where target_offset is the 0-based
        position of the variant within the returned sequence (may be < flank if
        the window is clamped at the start of the chromosome).
        """
        start_0 = variant.pos - 1 - flank   # 1-based pos → 0-based, then subtract flank
        end_0   = variant.pos - 1 + len(variant.ref) + flank
        clamped = max(0, start_0)
        seq     = self.fetch(variant.chrom, start_0, end_0)
        offset  = (variant.pos - 1) - clamped
        return seq, offset

    @property
    def references(self) -> list[str]:
        return list(self._fasta.references)

    def __enter__(self) -> ReferenceGenome:
        return self

    def __exit__(self, *_: object) -> None:
        self._fasta.close()

    def close(self) -> None:
        self._fasta.close()
