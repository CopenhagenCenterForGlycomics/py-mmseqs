"""Classes and functions to facilitate parsing MMseqs2 outputs."""

import csv
import re
import sys

from py_mmseqs.sequence_utilities import reverse_complement

# Make sure we can handle CSV files with wide fields - this is necessary for
# MMseqs2TableIterator to work properly if the 'qseq'/'tseq'/'qaln'/'taln'
# fields are present for long sequences
csv.field_size_limit(sys.maxsize)

CIGAR_REGEX = re.compile(r"(\d+)([A-Z])")


class MMseqs2HSP:
    """A class to represent a single high-scoring pair (HSP) from MMseqs2."""
    def __init__(self, **kwargs):
        """Initialize HSP.

        Not all parameters will always be used, depending on the type
        of search and the output columns requested.

        :param kwargs: keyword arguments
        """
        # Initialize query attributes as None
        self.query = None           # query sequence ID
        self.qheader = None         # full query header
        self.qseq = None            # query sequence
        self.qlen = None            # query sequence length
        self.qaln = None            # aligned query sequence
        self.qframe = None          # query reading frame
        self.qset = None            # the chunk name (for large query sets)
        self.qsetid = None          # the chunk ID (for large query sets)
        self.qorfstart = None       # query ORF start position
        self.qorfend = None         # query ORF end position

        # Initialize target attributes as None
        self.target = None          # target sequence ID
        self.theader = None         # full target header
        self.tseq = None            # target sequence
        self.tlen = None            # target sequence length
        self.taln = None            # aligned target sequence
        self.tframe = None          # target reading frame
        self.tset = None            # the db chunk name (for large target sets)
        self.tsetid = None          # the db chunk ID (for large target sets)
        self.torfstart = None       # target ORF start position
        self.torfend = None         # target ORF end position

        # Initialize alignment attributes as None
        self.nident = None          # number of identical residues
        self.pident = None          # percentage of identical residues
        self.fident = None          # fraction of identical residues
        self.qcov = None            # query coverage
        self.tcov = None            # target coverage
        self.alnlen = None          # alignment length
        self.mismatch = None        # number of mismatches
        self.gapopen = None         # number of gap openings
        self.qstart = None          # query start position
        self.qend = None            # query end position
        self.tstart = None          # target start position
        self.tend = None            # target end position
        self.raw = None             # raw alignment score
        self.bits = None            # bit score
        self.evalue = None          # E-value
        self.cigar = None           # CIGAR (alignment) string
        self.strand = None          # orientation of alignment

        # Initialize taxonomy attributes as None
        self.taxid = None           # taxonomy ID
        self.taxname = None         # taxonomy (common) name
        self.taxlineage = None      # taxonomy (complete) lineage

        # Set attributes from keyword arguments
        for keyword, argument in kwargs.items():
            # Convert string arguments to integers or floats if possible
            try:
                argument = int(argument)
            except ValueError:
                try:
                    argument = float(argument)
                except ValueError:
                    pass

            # Convert to None if argument is empty string
            if argument == "":
                argument = None

            # Set attribute
            setattr(self, keyword, argument)

        # Infer any attributes not explicitly set from keyword arguments
        self.gap_fill()

    def gap_fill(self):
        """Infer any attributes not explicitly set from keyword arguments.

        Inferences may be made for the following attributes:

        - query sequence ID
        - target sequence ID
        - orientation of alignment
        - query sequence length
        - target sequence length
        - query alignment
        - target alignment
        - number of identical residues
        - percentage of identical residues
        - fraction of identical residues
        - query coverage
        - target coverage
        - alignment length
        """
        # Infer query and target IDs from headers up to first whitespace
        if self.query is None and self.qheader is not None:
            self.query = self.qheader.split()[0]
        if self.target is None and self.theader is not None:
            self.target = self.theader.split()[0]

        # Infer orientation of alignment from query and target positions
        if (self.strand is None and self.qstart is not None and
                self.qend is not None):
            # Alignments are forward unless query start greater than query end
            if self.qstart > self.qend:
                self.strand = "Plus / Minus"
            else:
                self.strand = "Plus / Plus"

        # Infer query and target sequence lengths from sequences
        if self.qlen is None and self.qseq is not None:
            # Query length is the length of the query sequence
            self.qlen = len(self.qseq)
        if self.tlen is None and self.tseq is not None:
            # Target length is the length of the target sequence
            self.tlen = len(self.tseq)

        # REMAINING INFERENCES RELY ON PRESENCE OF CIGAR STRING
        # ------------------------------------------------------
        if self.cigar is None:
            # Cannot infer any remaining attributes without CIGAR string
            return

        # Query & target alignments, and alignment length are in CIGAR string
        if self.qaln is None or self.taln is None:
            self.qaln, self.taln = self._decode_cigar()

        # Alignment length is the length of the query/target alignment
        if len(self.qaln) != len(self.taln):
            print(repr(self))
            print(len(self.qaln))
            print(len(self.taln))
            raise ValueError("query/target alignments are different lengths")
        if self.alnlen is None:
            self.alnlen = len(self.qaln)

        # Infer number of identical residues from query and target alignments
        if self.nident is None:
            # If we got here, qaln and taln should be defined
            self.nident = self._count_nident()

        # Infer fraction of identical residues from number of identical
        # residues and alignment length
        if self.fident is None:
            # If we got here, nident and alnlen should be defined
            self.fident = self.nident / self.alnlen

        # Infer percentage of identical residues from fraction of identical
        # residues
        if self.pident is None:
            # If we got here, fident should be defined
            self.pident = self.fident * 100

        # Infer query coverage from alignment length and query length
        if self.qcov is None:
            # If we got here, alnlen and qlen should be defined
            self.qcov = min([1.0, self.alnlen / self.qlen])

        # Infer target coverage from alignment length and target length
        if self.tcov is None:
            # If we got here, alnlen and tlen should be defined
            self.tcov = min([1.0, self.alnlen / self.tlen])

    def _decode_cigar(self):
        """Use CIGAR string to decode query and target alignments."""
        if not self.cigar:
            return None, None
        elif not self.qseq or not self.qstart or not self.qend:
            return None, None
        elif not self.tseq or not self.tstart or not self.tend:
            return None, None

        # If alignment is on minus strand, re-orient the query sequence,
        # reverse-complement the target sequence, and reverse the cigar string
        if self.strand == "Plus / Minus":
            cigar = "".join(reversed(
                ["".join(x) for x in re.findall(CIGAR_REGEX, self.cigar)]))
            tseq = reverse_complement(self.tseq[self.tstart - 1:self.tend])
            qseq = self.qseq[self.qend - 1:self.qstart]
        else:
            cigar = self.cigar
            qseq = self.qseq[self.qstart - 1:self.qend]
            tseq = self.tseq[self.tstart - 1:self.tend]

        # Since we've sliced the sequences, neither start has any offset
        qstart, tstart = 0, 0

        # Initialize query and target alignments
        qaln, taln = "", ""
        for length, opcode in re.findall(CIGAR_REGEX, cigar):
            length = int(length)
            # Query matches target
            if opcode == "M":
                qaln += qseq[qstart:qstart + length]
                taln += tseq[tstart:tstart + length]
                qstart += length
                tstart += length
            # Query has insertion relative to target
            elif opcode == "I":
                qaln += qseq[qstart:qstart + length]
                taln += "-" * length
                qstart += length
            # Query has deletion relative to target
            elif opcode == "D":
                qaln += "-" * length
                taln += tseq[tstart:tstart + length]
                tstart += length
            # MMseqs2 only uses M, I, and D operations - no soft/hard clipping
            # or other operations
            else:
                raise ValueError("Unrecognized CIGAR operation: %s" % opcode)

        return qaln, taln

    def _count_nident(self):
        """Count the number of identical residues in the alignment. Only
        available if query and target alignments are defined."""
        if self.qaln is None or self.taln is None:
            return None

        nident = 0
        for qres, tres in zip(self.qaln, self.taln):
            if qres == tres:
                nident += 1

        return nident

    def as_tsv(self, keys=None):
        """Return a tab-separated representation of the HSP. If no keys
        are specified, all attributes will be included.

        :param keys: list of keys to include in the TSV [default: all]
        :type keys: list of str
        :return: tab-separated representation of the HSP
        :rtype: str
        """
        if keys is None:
            keys = self.__dict__.keys()

        return "\t".join([str(getattr(self, key)) for key in keys])

    def pretty_alignment(self):
        """Return string representation of HSP."""
        s = ""
        s += f"Query:    {self.query} ({self.qlen} bp)\n"
        s += f"Target:   {self.target} ({self.tlen} bp)\n"
        s += f"Score:    {self.bits} bits\n"
        s += f"Expect:   {self.evalue}\n"
        s += f"Identity: {self.nident}/{self.alnlen} ({self.pident:.1f}%)\n"
        s += f"Strand:   {self.strand}\n"
        s += f"\n\n"

        # Calculate the width of the position column
        width, max_len = 1, max([self.qlen, self.tlen])
        while max_len / 10 > 1:
            max_len /= 10
            width += 1

        qaln, taln = self.qaln, self.taln
        qstart, tstart = self.qstart, self.tstart
        # Alignment needs to compensate for strand and offset
        if self.strand == "Plus / Minus":
            # Reverse the query start and end positions
            qstart = self.qend

            # Reverse the target start and end positions
            tstart = self.tend

        # Now generate the alignment in chunks of 80 characters at a time
        for i in range(0, len(qaln), 80):
            qslice = qaln[i:i + 80]
            tslice = taln[i:i + 80]
            mslice = "".join(["|" if q == t else " " for q, t in
                              zip(qslice, tslice)])

            # Calculate the query and target start and end positions of this
            # alignment block
            qgaps = qslice.count("-")
            tgaps = tslice.count("-")

            qend = qstart + (len(qslice) - qgaps)
            if self.strand == "Plus / Minus":
                tend = tstart + (tgaps - len(tslice))
            else:
                tend = tstart + (len(tslice) - tgaps)

            # Print the alignment block
            s += f"Query  {qstart:{width}} {qslice} {qend - 1:{width}}\n"
            s += f"       {' ' * width} {mslice}\n"
            if self.strand == "Plus / Plus":
                s += f"Target {tstart:{width}} {tslice} {tend - 1:{width}}\n"
            else:
                s += f"Target {tstart:{width}} {tslice} {tend + 1:{width}}\n"

            s += f"\n"

            qstart, tstart = qend, tend

        return s

    def __str__(self):
        """Return string representation of HSP."""
        try:
            return self.pretty_alignment()
        except KeyError or AttributeError:
            return repr(self)

    def __repr__(self):
        """Return string representation of HSP."""
        return f"HSP({self.__dict__})"

    def __bool__(self):
        """Return True if any attribute is not None."""
        return all([v is None for v in self.__dict__.values()])


class MMseqs2TableIterator:
    """Iterator for MMseqs2 TSV output files."""
    def __init__(self, filename, header=None):
        """Initialize the table iterator.

        :param filename: path to the TSV file
        :type filename: os.PathLike
        :param header: list of column names [default: None]
        :type header: list of str
        """
        if not header:
            with open(filename) as f:
                header = f.readline().strip().split("\t")

        self.header = header
        self.filename = filename

    def __iter__(self):
        """Iterate over rows in tabular output file."""
        with open(self.filename) as f:
            reader = csv.DictReader(f, delimiter="\t", fieldnames=self.header)
            next(reader)  # Skip header
            for row in reader:
                row: dict[str, str]  # Type hint for PyCharm
                try:
                    yield MMseqs2HSP(**row)
                except KeyError:
                    yield MMseqs2HSP()


class MMseqs2Cluster:
    """A class to represent a single cluster from MMseqs2."""
    def __init__(self, headers=None, sequences=None):
        self._headers = []
        self._sequences = []
        self._representative = None

        if not headers:
            headers = []
        if not sequences:
            sequences = []

        if len(headers) != len(sequences):
            raise ValueError("headers and sequences must be co-linear lists")

        [self.add_sequence(h, s) for h, s in zip(headers, sequences)]

    def add_sequence(self, header, sequence):
        """Add a (header, sequence) pair to the cluster. Headers are tested
        for uniqueness."""
        if header in self:
            raise ValueError(f"'{header}' is already present in the cluster")

        self._headers.append(header)
        self._sequences.append(sequence)

    @property
    def representative(self):
        """Return the representative sequence of the cluster."""
        return self._representative

    @representative.setter
    def representative(self, value):
        """Set the representative sequence of the cluster."""
        if value not in self:
            raise ValueError(f"the representative must be a member of the "
                             f"cluster; '{value}' not present")

        self._representative = value

    @property
    def headers(self):
        """Return the headers of the cluster."""
        return self._headers

    @property
    def sequences(self):
        """Return the sequences of the cluster."""
        return self._sequences

    def __iter__(self):
        """Iterate over (sequence, header) pairs in the cluster."""
        return iter(zip(self._headers, self._sequences))

    def __str__(self):
        """Return a string (FASTA) representation of the cluster."""
        return "".join(f">{h}\n{s}\n" for h, s in self)

    def __len__(self):
        """Return the size of the cluster."""
        return len(self._headers)

    def __bool__(self):
        """Return True if the cluster is non-empty."""
        return len(self) > 0

    def __contains__(self, item):
        """Test if the cluster contains a header."""
        return item in self._headers


class MMseqs2ClusterIterator:
    """Iterator for MMseqs2 FASTA-like cluster output files."""
    def __init__(self, filename):
        self.filename = filename

    def __iter__(self):
        headers, sequences = list(), list()
        with open(self.filename, "r") as cluster_reader:
            # Keep a reference to the previous and current lines, so we can
            # check for the presence of the duplicate header used to indicate
            # the start of a new cluster (and its representative)
            prior = cluster_reader.readline()
            current = cluster_reader.readline()

            # While loop to iterate until EOF
            while current:
                # Current line is a header
                if current.startswith(">"):
                    # Previous line is a header and sequences are present
                    if prior.startswith(">") and sequences:
                        # Pop the new cluster's representative from the headers
                        headers.pop(-1)

                        # Old cluster is complete, so yield it
                        yield MMseqs2Cluster(headers, sequences)
                        headers, sequences = list(), list()

                    # Current line is a header, so add it to the headers
                    headers.append(current[1:].rstrip())
                else:
                    # Current line is a sequence, so add it to the sequences
                    sequences.append(current.rstrip())

                prior, current = current, cluster_reader.readline()

            # Yield the last cluster
            yield MMseqs2Cluster(headers, sequences)


def parse_tsv(tsv_file, header=None):
    """Parse an MMseqs2 TSV file and return a list of MMseqs2HSP objects."""
    return MMseqs2TableIterator(tsv_file, header=header)


def parse_clusters(cluster_file):
    """Parse an MMseqs2 cluster file and return a list of MMseqs2Cluster
    objects."""
    return MMseqs2ClusterIterator(cluster_file)
