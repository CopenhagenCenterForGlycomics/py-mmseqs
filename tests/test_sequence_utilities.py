import unittest

from py_mmseqs.sequence_utilities import *


class Test(unittest.TestCase):
    def test_complement(self):
        """Verify that the complement of a sequence is correct."""
        self.assertEqual(complement("TCAG"), "AGTC")

    def test_complement_2(self):
        """Verify that complement is round-trip stable."""
        self.assertEqual(complement(complement("TCAG")), "TCAG")

    def test_complement_3(self):
        """Verify that complement with non-TCAG characters raises a KeyError."""
        with self.assertRaises(KeyError):
            complement("TCAGX")

    def test_reverse_complement(self):
        """Verify that the reverse-complement of a sequence is correct."""
        self.assertEqual(reverse_complement("TCAG"), "CTGA")

    def test_reverse_complement_2(self):
        """Verify that reverse-complement is round-trip stable."""
        self.assertEqual(reverse_complement(reverse_complement("TCAG")), "TCAG")

    def test_reverse_complement_3(self):
        """Verify that reverse-complement with non-TCAG characters raises a
        KeyError."""
        with self.assertRaises(KeyError):
            reverse_complement("TCAGX")

    def test_subsequence(self):
        """Verify that subsequences are correctly extracted."""
        self.assertEqual(subsequence(
            "TCAGTTCCAAGGTTTCCCAAAGGG", 1, 20),
                         "TCAGTTCCAAGGTTTCCCAA")

    def test_subsequence_2(self):
        """Verify that reverse-strand subsequences are correctly extracted."""
        self.assertEqual(subsequence(
            "TCAGTTCCAAGGTTTCCCAAAGGG", 1, 20, fwd_strand=False),
                         "TTGGGAAACCTTGGAACTGA")


if __name__ == '__main__':
    unittest.main()
