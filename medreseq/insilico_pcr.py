"""In silico PCR: validate primer pairs against a reference sequence window."""

from __future__ import annotations

from dataclasses import dataclass, field

from .primer_design import PrimerPair

_COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT)[::-1]


@dataclass
class PCRResult:
    found: bool
    product_size: int | None = None
    product_seq: str | None = None
    off_target_count: int = 0   # additional products beyond the primary one
    notes: list[str] = field(default_factory=list)


def _find_all(haystack: str, needle: str) -> list[int]:
    """Return all start positions of needle in haystack."""
    positions: list[int] = []
    start = 0
    while True:
        idx = haystack.find(needle, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def run_insilico_pcr(
    reference_window: str,
    pair: PrimerPair,
    min_product: int = 50,
    max_product: int = 2000,
) -> PCRResult:
    """Simulate PCR within a reference sequence window.

    Searches for exact matches of the forward primer and the reverse complement
    of the reverse primer. Products are validated by orientation (forward 5' of
    reverse) and size.

    `reference_window` should span at least the expected amplicon, e.g. the
    region ±pcr_window bp around the target variant.
    """
    ref = reference_window.upper()
    fwd = pair.left_seq.upper()
    rev_rc = reverse_complement(pair.right_seq)  # how the rev primer appears on the fwd strand

    fwd_hits = _find_all(ref, fwd)
    rev_hits = _find_all(ref, rev_rc)

    products: list[tuple[int, int, int]] = []   # (start, end, size)
    for fp in fwd_hits:
        for rp in rev_hits:
            end = rp + len(rev_rc)
            size = end - fp
            if min_product <= size <= max_product:
                products.append((fp, end, size))

    if not products:
        notes: list[str] = []
        if not fwd_hits:
            notes.append("forward primer not found in search window")
        if not rev_hits:
            notes.append("reverse primer (RC) not found in search window")
        if fwd_hits and rev_hits:
            notes.append(
                "primers found but no valid product — check orientation, "
                "product size range, or increase --pcr-window"
            )
        return PCRResult(found=False, notes=notes)

    # Prefer the product closest to the centre of the window (most likely on-target)
    mid = len(ref) // 2
    best = min(products, key=lambda p: abs((p[0] + p[1]) // 2 - mid))

    return PCRResult(
        found=True,
        product_size=best[2],
        product_seq=ref[best[0]:best[1]],
        off_target_count=max(0, len(products) - 1),
    )
