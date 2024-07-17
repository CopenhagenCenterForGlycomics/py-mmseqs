"""Functions to manipulate DNA sequences."""

COMPLEMENTS = {"T": "A", "C": "G", "A": "T", "G": "C"}


def complement(sequence):
    """Return the complement of a sequence.

    :param sequence: a sequence to complement
    :type sequence: str
    """
    return "".join([COMPLEMENTS[base] for base in sequence])


def reverse_complement(sequence):
    """Return the reverse-complement of a sequence.

    :param sequence: a sequence to reverse-complement
    :type sequence: str
    """
    return complement(sequence)[::-1]


def subsequence(sequence, left_coord, right_coord, fwd_strand=True):
    """Return a sub-sequence of a sequence of nucleotides.

    If `fwd_strand` is False, the reverse strand subsequence will be
    returned instead of the forward strand.

    :param sequence: the sequence to slice
    :type sequence: str
    :param left_coord: the smaller coordinate (1-based) to slice from
    :type left_coord: int
    :param right_coord: the larger coordinate (1-based) to slice to
    :type right_coord: int
    :param fwd_strand: indicate whether to subsequence the top strand
    :type fwd_strand: bool
    :return: subseq
    """
    subseq = ""

    try:
        if fwd_strand:
            subseq = sequence[left_coord - 1:right_coord]
        else:
            subseq = reverse_complement(sequence[left_coord - 1:right_coord])

    except IndexError:
        pass

    return subseq


__all__ = ["complement", "reverse_complement", "subsequence"]
