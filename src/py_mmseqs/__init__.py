"""Classes and utility functions for interfacing with MMseqs2 and parsing its
outputs."""

from py_mmseqs.mmseqs2 import MMseqs2, MMseqs2Error
from py_mmseqs.parsers import parse_clusters, parse_tsv


__version__ = '0.1.0'
__all__ = ['MMseqs2', 'MMseqs2Error', 'parse_clusters', 'parse_tsv']
