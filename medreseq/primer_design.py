"""Primer design via the primer3-py library."""

from __future__ import annotations

from dataclasses import dataclass

import primer3


@dataclass(frozen=True)
class PrimerPair:
    rank: int           # 0-based rank returned by primer3 (lower = better penalty)
    left_seq: str
    right_seq: str
    left_tm: float
    right_tm: float
    left_gc: float
    right_gc: float
    product_size: int
    penalty: float


def design_primers(
    sequence: str,
    target_offset: int,
    target_len: int,
    *,
    num_return: int = 3,
    product_size_range: tuple[int, int] = (200, 800),
    min_tm: float = 57.0,
    opt_tm: float = 60.0,
    max_tm: float = 63.0,
) -> list[PrimerPair]:
    """Design primer pairs that flank a target region within `sequence`.

    Args:
        sequence:         Template sequence (uppercase recommended).
        target_offset:    0-based start of the target region within `sequence`.
        target_len:       Length of the target region (e.g. len(ref_allele)).
        num_return:       Maximum number of primer pairs to return.
        product_size_range: (min, max) amplicon size in bp.
        min_tm / opt_tm / max_tm: Melting temperature constraints.

    Returns:
        List of PrimerPair objects ordered by primer3 penalty (best first).
        Returns an empty list if primer3 cannot design primers.
    """
    result = primer3.design_primers(
        seq_args={
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_TARGET': [target_offset, max(target_len, 1)],
        },
        global_args={
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_NUM_RETURN': num_return,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_MIN_TM': min_tm,
            'PRIMER_OPT_TM': opt_tm,
            'PRIMER_MAX_TM': max_tm,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_PRODUCT_SIZE_RANGE': [list(product_size_range)],
        },
    )

    pairs: list[PrimerPair] = []
    i = 0
    while f'PRIMER_LEFT_{i}_SEQUENCE' in result:
        pairs.append(PrimerPair(
            rank=i,
            left_seq=result[f'PRIMER_LEFT_{i}_SEQUENCE'],
            right_seq=result[f'PRIMER_RIGHT_{i}_SEQUENCE'],
            left_tm=result[f'PRIMER_LEFT_{i}_TM'],
            right_tm=result[f'PRIMER_RIGHT_{i}_TM'],
            left_gc=result[f'PRIMER_LEFT_{i}_GC_PERCENT'],
            right_gc=result[f'PRIMER_RIGHT_{i}_GC_PERCENT'],
            product_size=result[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            penalty=result[f'PRIMER_PAIR_{i}_PENALTY'],
        ))
        i += 1
    return pairs
