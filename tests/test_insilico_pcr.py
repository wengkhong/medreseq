"""Tests for in silico PCR (pure Python — no pysam/primer3 required)."""

from __future__ import annotations

import pytest

from medreseq.insilico_pcr import PCRResult, reverse_complement, run_insilico_pcr
from medreseq.primer_design import PrimerPair


# ── Helpers ───────────────────────────────────────────────────────────────────

def _pair(left: str, right: str) -> PrimerPair:
    """Build a minimal PrimerPair for PCR testing."""
    return PrimerPair(
        rank=0, left_seq=left, right_seq=right,
        left_tm=60.0, right_tm=60.0,
        left_gc=50.0, right_gc=50.0,
        product_size=0, penalty=0.0,
    )


FWD = "ACGTACGTACGTACGTACGT"   # 20 bp
REV = "GCATGCATGCATGCATGCAT"   # 20 bp
REV_RC = reverse_complement(REV)


# ── reverse_complement ────────────────────────────────────────────────────────

def test_reverse_complement_basic() -> None:
    assert reverse_complement("ACGT") == "ACGT"


def test_reverse_complement_asymmetric() -> None:
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement("GCGC") == "GCGC"


def test_reverse_complement_n() -> None:
    assert reverse_complement("ACGTN") == "NACGT"


def test_reverse_complement_involution() -> None:
    seq = "ATCGATCGATCG"
    assert reverse_complement(reverse_complement(seq)) == seq


# ── run_insilico_pcr — found ──────────────────────────────────────────────────

def test_primers_found_correct_size() -> None:
    gap = "T" * 200
    ref = "A" * 100 + FWD + gap + REV_RC + "A" * 100
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert result.found
    expected = len(FWD) + len(gap) + len(REV_RC)   # 20 + 200 + 20 = 240
    assert result.product_size == expected
    assert result.product_seq is not None
    assert len(result.product_seq) == expected


def test_product_seq_starts_with_fwd_ends_with_rev_rc() -> None:
    gap = "C" * 150
    ref = "G" * 50 + FWD + gap + REV_RC + "G" * 50
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=50, max_product=500)

    assert result.found
    assert result.product_seq.startswith(FWD)
    assert result.product_seq.endswith(REV_RC)


def test_no_off_targets_when_primers_appear_once() -> None:
    gap = "T" * 200
    ref = "A" * 100 + FWD + gap + REV_RC + "A" * 100
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert result.off_target_count == 0


def test_off_targets_detected_when_primers_appear_twice() -> None:
    gap = "T" * 200
    unit = FWD + gap + REV_RC
    ref = unit + "A" * 50 + unit    # two amplifiable copies
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert result.found
    assert result.off_target_count >= 1


# ── run_insilico_pcr — not found ─────────────────────────────────────────────

def test_primers_not_found_in_empty_window() -> None:
    result = run_insilico_pcr("A" * 1000, _pair(FWD, REV))
    assert not result.found
    assert any("forward primer" in n for n in result.notes)
    assert any("reverse primer" in n for n in result.notes)


def test_wrong_orientation_gives_no_product() -> None:
    # rev_rc appears BEFORE fwd → no valid (fwd < rev) product
    gap = "T" * 200
    ref = "A" * 50 + REV_RC + gap + FWD + "A" * 50
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert not result.found
    assert result.notes  # should explain why


def test_product_too_small_not_reported() -> None:
    gap = "T" * 5     # product = 20 + 5 + 20 = 45 bp
    ref = "A" * 50 + FWD + gap + REV_RC + "A" * 50
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert not result.found


def test_product_too_large_not_reported() -> None:
    gap = "T" * 2000   # product = 20 + 2000 + 20 = 2040 bp
    ref = FWD + gap + REV_RC
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert not result.found


# ── Centre-picking heuristic ─────────────────────────────────────────────────

def test_on_target_product_preferred_over_off_target() -> None:
    """When two valid products exist the one nearest the window centre is picked."""
    gap = "T" * 200
    unit = FWD + gap + REV_RC   # 240 bp
    # Place off-target at start, on-target in the middle
    ref = unit + "A" * 500 + unit + "A" * 500
    # Centre of ref ≈ position 870; second unit starts at 740 → closer to centre
    result = run_insilico_pcr(ref, _pair(FWD, REV), min_product=100, max_product=500)

    assert result.found
    assert result.product_size == len(unit)
