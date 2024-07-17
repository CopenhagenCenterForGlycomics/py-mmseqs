"""Unit tests for the py_mmseqs.mmseqs2 module."""

import pathlib
import shutil
import unittest

from py_mmseqs.mmseqs2 import MMseqs2, MMseqs2Error, _cast_list_objs_to_paths


class TestCastListObjsToPaths(unittest.TestCase):
    """Most important tests suite - almost all MMseqs2 class methods
    require this function to work as expected."""
    def test_cast_list_objs_to_paths_1(self):
        """Verify that we can cast a list of pathlib.Path objects to a list
        of pathlib.Path objects."""
        obj_list = [pathlib.Path("/path/to/file1"),
                    pathlib.Path("/path/to/file2"),
                    pathlib.Path("/path/to/file3")]
        path_list = _cast_list_objs_to_paths(obj_list)

        # We should get a list back
        self.assertIsInstance(path_list, list)

        # The list should be the same length as the input list
        self.assertEqual(len(path_list), len(obj_list))

        # All elements in the returned list should be pathlib.Path objects
        for path in path_list:
            self.assertIsInstance(path, pathlib.Path)

    def test_cast_list_objs_to_paths_2(self):
        """Verify that we can cast a list of pathlike strings to a list
        of pathlib.Path objects."""
        obj_list = ["/path/to/file1", "/path/to/file2", "/path/to/file3"]
        path_list = _cast_list_objs_to_paths(obj_list)

        # We should get a list back
        self.assertIsInstance(path_list, list)

        # The list should be the same length as the input list
        self.assertEqual(len(path_list), len(obj_list))

        # All elements in the returned list should be pathlib.Path objects
        for path in path_list:
            self.assertIsInstance(path, pathlib.Path)

    def test_cast_list_objs_to_paths_3(self):
        """Verify that we can cast a list of strings that don't look like
        paths to a list of pathlib.Path objects."""
        obj_list = ["file1", "file2", "file3"]
        path_list = _cast_list_objs_to_paths(obj_list)

        # We should get a list back
        self.assertIsInstance(path_list, list)

        # The list should be the same length as the input list
        self.assertEqual(len(path_list), len(obj_list))

        # All elements in the returned list should be pathlib.Path objects
        for path in path_list:
            self.assertIsInstance(path, pathlib.Path)

    def test_cast_list_objs_to_paths_4(self):
        """Verify that a TypeError is raised when non-pathlike objects are
        passed to the function."""
        obj_list = [1, 2, 3]
        with self.assertRaises(TypeError):
            _cast_list_objs_to_paths(obj_list)


# Test easy pipelines - easy-linsearch omitted for now
class TestEasySearchAminoAcid(unittest.TestCase):
    """Test the 'easy-search' pipeline with amino acid sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several amino acid FASTA files to tests with, and one
        # compressed FASTA file
        self.aa_fasta1 = self.test_data_dir.joinpath("aa_1.fasta")
        self.aa_fasta2 = self.test_data_dir.joinpath("aa_2.fasta")
        self.aa_fasta3 = self.test_data_dir.joinpath("aa_3.fasta")
        self.aa_fasta_gz = self.test_data_dir.joinpath("aa_all.fasta.gz")

        # Verify that all tests FASTA files exist
        for test_file in (self.aa_fasta1, self.aa_fasta2, self.aa_fasta3,
                          self.aa_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define a temporary directory for MMseqs2 outputs
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define search output files
        self.search_outfile = self.tmpdir.joinpath("test_easy_search_aa.tsv")

    def test_one_query(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        amino acid query file against an amino acid target file with no
        keyword arguments."""
        stdout = self.mmseqs.easy_search(self.aa_fasta1, self.aa_fasta2,
                                         self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_two_query(self):
        """Verify that we can run the 'easy-search' pipeline with two
        amino acid query files against an amino acid target file with no
        keyword arguments."""
        stdout = self.mmseqs.easy_search(
            [self.aa_fasta1, self.aa_fasta2], self.aa_fasta3,
            self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_valid_kwargs(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        amino acid query file against an amino acid target file with valid
        keyword arguments."""
        stdout = self.mmseqs.easy_search(self.aa_fasta1, self.aa_fasta2,
                                         self.search_outfile, self.tmpdir,
                                         min_seq_id=0.5, e=0.01, c=0.5)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_invalid_kwargs(self):
        """Verify that the 'easy-search' pipeline with a single amino acid
        query file against an amino acid target file with invalid keyword
        arguments raises a KeyError."""
        with self.assertRaises(KeyError):
            self.mmseqs.easy_search(self.aa_fasta1, self.aa_fasta2,
                                    self.search_outfile, self.tmpdir,
                                    max_seq_id=0.5, evalue=0.01, coverage=0.5)

        self.assertFalse(self.search_outfile.exists())

    def test_valid_kw_bad_type(self):
        """Verify that the 'easy-search' pipeline with a valid keyword
        argument but the wrong type raises a TypeError."""
        with self.assertRaises(TypeError):
            self.mmseqs.easy_search(self.aa_fasta1, self.aa_fasta2,
                                    self.search_outfile, self.tmpdir,
                                    min_seq_id="0.5", e=0.01, c=0.5)

        self.assertFalse(self.search_outfile.exists())

    def test_gzipped_target(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        amino acid query file against a gzipped amino acid target file."""
        stdout = self.mmseqs.easy_search(self.aa_fasta1, self.aa_fasta_gz,
                                         self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestEasySearchNucleotide(unittest.TestCase):
    """Test the 'easy-search' pipeline with nucleotide sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several nucleotide FASTA files to tests with, and one
        # compressed FASTA file
        self.nt_fasta1 = self.test_data_dir.joinpath("nt_1.fasta")
        self.nt_fasta2 = self.test_data_dir.joinpath("nt_2.fasta")
        self.nt_fasta3 = self.test_data_dir.joinpath("nt_3.fasta")
        self.nt_fasta_gz = self.test_data_dir.joinpath("nt_all.fasta.gz")

        # Verify that all tests FASTA files exist
        for test_file in (self.nt_fasta1, self.nt_fasta2, self.nt_fasta3,
                          self.nt_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define a temporary directory for MMseqs2 outputs
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define search output files
        self.search_outfile = self.tmpdir.joinpath("test_easy_search_nt.tsv")

    def test_one_query(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        nucleotide query file against a nucleotide target file with the
        search-type specified."""
        stdout = self.mmseqs.easy_search(self.nt_fasta1, self.nt_fasta2,
                                         self.search_outfile, self.tmpdir,
                                         search_type=3)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_one_query_no_search_type(self):
        """Verify that the 'easy-search' pipeline with a single nucleotide
        query file against a nucleotide target file with no keyword
        arguments raises an MMseqs2Error."""
        with self.assertRaises(MMseqs2Error):
            self.mmseqs.easy_search(self.nt_fasta1, self.nt_fasta2,
                                    self.search_outfile, self.tmpdir)

        self.assertFalse(self.search_outfile.exists())

    def test_two_query(self):
        """Verify that we can run the 'easy-search' pipeline with two
        nucleotide query files against a nucleotide target file with no
        keyword arguments."""
        stdout = self.mmseqs.easy_search(
            [self.nt_fasta1, self.nt_fasta2], self.nt_fasta3,
            self.search_outfile, self.tmpdir, search_type=3)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_valid_kwargs(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        nucleotide query file against a nucleotide target file with valid
        keyword arguments."""
        stdout = self.mmseqs.easy_search(self.nt_fasta1, self.nt_fasta2,
                                         self.search_outfile, self.tmpdir,
                                         min_seq_id=0.5, e=0.01, c=0.5,
                                         search_type=3)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_invalid_kwargs(self):
        """Verify that the 'easy-search' pipeline with a single nucleotide
        query file against a nucleotide target file with invalid keyword
        arguments raises a KeyError."""
        with self.assertRaises(KeyError):
            self.mmseqs.easy_search(self.nt_fasta1, self.nt_fasta2,
                                    self.search_outfile, self.tmpdir,
                                    max_seq_id=0.5, evalue=0.01, coverage=0.5)

        self.assertFalse(self.search_outfile.exists())

    def test_valid_kw_bad_type(self):
        """Verify that the 'easy-search' pipeline with a valid keyword
        argument but the wrong type raises a TypeError."""
        with self.assertRaises(TypeError):
            self.mmseqs.easy_search(self.nt_fasta1, self.nt_fasta2,
                                    self.search_outfile, self.tmpdir,
                                    min_seq_id="0.5", e=0.01, c=0.5)

        self.assertFalse(self.search_outfile.exists())

    def test_gzipped_target(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        nucleotide query file against a gzipped nucleotide target file."""
        stdout = self.mmseqs.easy_search(self.nt_fasta1, self.nt_fasta_gz,
                                         self.search_outfile, self.tmpdir,
                                         search_type=3)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestEasySearchMixed(unittest.TestCase):
    """Test the 'easy-search' pipeline with amino acid and nucleotide
    sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several amino acid and nucleotide FASTA files to tests with
        self.aa_fasta1 = self.test_data_dir.joinpath("aa_1.fasta")
        self.aa_fasta2 = self.test_data_dir.joinpath("aa_2.fasta")
        self.nt_fasta1 = self.test_data_dir.joinpath("nt_1.fasta")
        self.nt_fasta2 = self.test_data_dir.joinpath("nt_2.fasta")

        # Verify that all tests FASTA files exist
        for test_file in (
                self.aa_fasta1, self.aa_fasta2, self.nt_fasta1, self.nt_fasta2):
            self.assertTrue(test_file.exists())

        # We also have a pair of gzip-compressed FASTA files
        self.aa_fasta_gz = self.test_data_dir.joinpath("aa_all.fasta.gz")
        self.nt_fasta_gz = self.test_data_dir.joinpath("nt_all.fasta.gz")

        # Define and create a temporary directory for MMseqs2 output
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define search output files
        self.search_outfile = self.tmpdir.joinpath("search_mixed.tsv")

    def test_one_nt_query_aa_target(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        nucleotide query file against an amino acid target file."""
        stdout = self.mmseqs.easy_search(self.nt_fasta1, self.aa_fasta1,
                                         self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_two_nt_query_aa_target(self):
        """Verify that we can run the 'easy-search' pipeline with multiple
        nucleotide query files against an amino acid target file."""
        stdout = self.mmseqs.easy_search(
            [self.nt_fasta1, self.nt_fasta2], self.aa_fasta1,
            self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_one_aa_query_nt_target(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        amino acid query file against a nucleotide target file."""
        stdout = self.mmseqs.easy_search(self.aa_fasta1, self.nt_fasta1,
                                         self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_two_aa_query_nt_target(self):
        """Verify that we can run the 'easy-search' pipeline with multiple
        amino acid query files against a nucleotide target file."""
        stdout = self.mmseqs.easy_search(
            [self.aa_fasta1, self.aa_fasta2], self.nt_fasta1,
            self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_nt_query_aa_target_gz(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        nucleotide query file against a gzipped amino acid target file."""
        stdout = self.mmseqs.easy_search(self.nt_fasta1, self.aa_fasta_gz,
                                         self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def test_aa_query_nt_target_gz(self):
        """Verify that we can run the 'easy-search' pipeline with a single
        amino acid query file against a gzipped nucleotide target file."""
        stdout = self.mmseqs.easy_search(self.aa_fasta1, self.nt_fasta_gz,
                                         self.search_outfile, self.tmpdir)

        self.assertTrue(self.search_outfile.exists())
        self.assertIsInstance(stdout, str)

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestEasyClusterAminoAcid(unittest.TestCase):
    """Test the 'easy-cluster' pipeline with amino acid sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several amino acid FASTA files to tests with, and one
        # compressed FASTA file
        self.aa_fasta1 = self.test_data_dir.joinpath("aa_1.fasta")
        self.aa_fasta2 = self.test_data_dir.joinpath("aa_2.fasta")
        self.aa_fasta3 = self.test_data_dir.joinpath("aa_3.fasta")
        self.aa_fasta_gz = self.test_data_dir.joinpath("aa_all.fasta.gz")

        # Verify that all tests FASTA files exist
        for test_file in (self.aa_fasta1, self.aa_fasta2, self.aa_fasta3,
                          self.aa_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define a temporary directory for MMseqs2 outputs
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define cluster output file names
        self.cluster_prefix = self.tmpdir.joinpath("test_easy_cluster_aa")
        suffices = ("_cluster.tsv", "_all_seqs.fasta", "_rep_seq.fasta")
        self.cluster_outfiles = (pathlib.Path(str(self.cluster_prefix) + s) for
                                 s in suffices)

    def test_one_query(self):
        """Verify that we can run the 'easy-cluster' pipeline with a single
        amino acid query file."""
        stdout = self.mmseqs.easy_cluster(self.aa_fasta1, self.cluster_prefix,
                                          self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_two_query(self):
        """Verify that we can run the 'easy-cluster' pipeline with two
        amino acid query files."""
        stdout = self.mmseqs.easy_cluster([self.aa_fasta1, self.aa_fasta2],
                                          self.cluster_prefix, self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_gzipped_query(self):
        """Verify that we can run the 'easy-cluster' pipeline with a gzipped
        amino acid query file."""
        stdout = self.mmseqs.easy_cluster(self.aa_fasta_gz, self.cluster_prefix,
                                          self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_valid_kwargs(self):
        """Verify that we can run the 'easy-cluster' pipeline with a single
        amino acid query file and valid keyword arguments."""
        stdout = self.mmseqs.easy_cluster(self.aa_fasta1, self.cluster_prefix,
                                          self.tmpdir, min_seq_id=0.35,
                                          e=0.01, c=0.75)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_invalid_kwargs(self):
        """Verify that the 'easy-cluster' pipeline with a single amino acid
        query file and invalid keyword arguments raises a KeyError."""
        with self.assertRaises(KeyError):
            self.mmseqs.easy_cluster(self.aa_fasta1, self.cluster_prefix,
                                     self.tmpdir, max_seq_id=0.35, evalue=0.01,
                                     coverage=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def test_valid_kw_bad_type(self):
        """Verify that the 'easy-cluster' pipeline with a valid keyword
        argument but the wrong type raises a TypeError."""
        with self.assertRaises(TypeError):
            self.mmseqs.easy_cluster(self.aa_fasta1, self.cluster_prefix,
                                     self.tmpdir, min_seq_id="0.35",
                                     e=0.01, c=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestEasyClusterNucleotide(unittest.TestCase):
    """Test the 'easy-cluster' pipeline with nucleotide sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several nucleotide FASTA files to tests with, and one
        # compressed FASTA file
        self.nt_fasta1 = self.test_data_dir.joinpath("nt_1.fasta")
        self.nt_fasta2 = self.test_data_dir.joinpath("nt_2.fasta")
        self.nt_fasta3 = self.test_data_dir.joinpath("nt_3.fasta")
        self.nt_fasta_gz = self.test_data_dir.joinpath("nt_all.fasta.gz")

        # Verify that all tests FASTA files exist
        for test_file in (self.nt_fasta1, self.nt_fasta2, self.nt_fasta3,
                          self.nt_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define a temporary directory for MMseqs2 outputs
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define cluster output file names
        self.cluster_prefix = self.tmpdir.joinpath("test_easy_cluster_aa")
        suffices = ("_cluster.tsv", "_all_seqs.fasta", "_rep_seq.fasta")
        self.cluster_outfiles = (pathlib.Path(str(self.cluster_prefix) + s) for
                                 s in suffices)

    def test_one_query(self):
        """Verify that we can run the 'easy-cluster' pipeline with a single
        nucleotide query file."""
        stdout = self.mmseqs.easy_cluster(self.nt_fasta1, self.cluster_prefix,
                                          self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_two_query(self):
        """Verify that we can run the 'easy-cluster' pipeline with two
        nucleotide query files."""
        stdout = self.mmseqs.easy_cluster([self.nt_fasta1, self.nt_fasta2],
                                          self.cluster_prefix, self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_gzipped_query(self):
        """Verify that we can run the 'easy-cluster' pipeline with a gzipped
        nucleotide query file."""
        stdout = self.mmseqs.easy_cluster(self.nt_fasta_gz, self.cluster_prefix,
                                          self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_valid_kwargs(self):
        """Verify that we can run the 'easy-cluster' pipeline with a single
        nucleotide query file and valid keyword arguments."""
        stdout = self.mmseqs.easy_cluster(self.nt_fasta1, self.cluster_prefix,
                                          self.tmpdir, min_seq_id=0.35,
                                          e=0.01, c=0.75)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_invalid_kwargs(self):
        """Verify that the 'easy-cluster' pipeline with a single nucleotide
        query file and invalid keyword arguments raises a KeyError."""
        with self.assertRaises(KeyError):
            self.mmseqs.easy_cluster(self.nt_fasta1, self.cluster_prefix,
                                     self.tmpdir, max_seq_id=0.35, evalue=0.01,
                                     coverage=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def test_valid_kw_bad_type(self):
        """Verify that the 'easy-cluster' pipeline with a valid keyword
        argument but the wrong type raises a TypeError."""
        with self.assertRaises(TypeError):
            self.mmseqs.easy_cluster(self.nt_fasta1, self.cluster_prefix,
                                     self.tmpdir, min_seq_id="0.35",
                                     e=0.01, c=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestEasyLinclustAminoAcid(unittest.TestCase):
    """Test the 'easy-linclust' pipeline with amino acid sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several amino acid FASTA files to tests with, and one
        # compressed FASTA file
        self.aa_fasta1 = self.test_data_dir.joinpath("aa_1.fasta")
        self.aa_fasta2 = self.test_data_dir.joinpath("aa_2.fasta")
        self.aa_fasta3 = self.test_data_dir.joinpath("aa_3.fasta")
        self.aa_fasta_gz = self.test_data_dir.joinpath("aa_all.fasta.gz")

        # Verify that all tests FASTA files exist
        for test_file in (self.aa_fasta1, self.aa_fasta2, self.aa_fasta3,
                          self.aa_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define a temporary directory for MMseqs2 outputs
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define cluster output file names
        self.cluster_prefix = self.tmpdir.joinpath("test_easy_linclust_aa")
        suffices = ("_cluster.tsv", "_all_seqs.fasta", "_rep_seq.fasta")
        self.cluster_outfiles = (pathlib.Path(str(self.cluster_prefix) + s) for
                                 s in suffices)

    def test_one_query(self):
        """Verify that we can run the 'easy-linclust' pipeline with a single
        amino acid query file."""
        stdout = self.mmseqs.easy_linclust(self.aa_fasta1, self.cluster_prefix,
                                           self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_two_query(self):
        """Verify that we can run the 'easy-linclust' pipeline with two
        amino acid query files."""
        stdout = self.mmseqs.easy_linclust([self.aa_fasta1, self.aa_fasta2],
                                           self.cluster_prefix, self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_gzipped_query(self):
        """Verify that we can run the 'easy-linclust' pipeline with a gzipped
        amino acid query file."""
        stdout = self.mmseqs.easy_linclust(self.aa_fasta_gz, self.cluster_prefix,
                                           self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_valid_kwargs(self):
        """Verify that we can run the 'easy-linclust' pipeline with a single
        amino acid query file and valid keyword arguments."""
        stdout = self.mmseqs.easy_linclust(self.aa_fasta1, self.cluster_prefix,
                                           self.tmpdir, min_seq_id=0.9, c=0.75)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_invalid_kwargs(self):
        """Verify that the 'easy-linclust' pipeline with a single amino acid
        query file and invalid keyword arguments raises a KeyError."""
        with self.assertRaises(KeyError):
            self.mmseqs.easy_linclust(self.aa_fasta1, self.cluster_prefix,
                                      self.tmpdir, max_seq_id=0.9,
                                      coverage=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def test_valid_kw_bad_type(self):
        """Verify that the 'easy-linclust' pipeline with a valid keyword
        argument but the wrong type raises a TypeError."""
        with self.assertRaises(TypeError):
            self.mmseqs.easy_linclust(self.aa_fasta1, self.cluster_prefix,
                                      self.tmpdir, min_seq_id="0.9", c=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestEasyLinclustNucleotide(unittest.TestCase):
    """Test the 'easy-linclust' pipeline with nucleotide sequences."""
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several nucleotide FASTA files to tests with, and one
        # compressed FASTA file
        self.nt_fasta1 = self.test_data_dir.joinpath("nt_1.fasta")
        self.nt_fasta2 = self.test_data_dir.joinpath("nt_2.fasta")
        self.nt_fasta3 = self.test_data_dir.joinpath("nt_3.fasta")
        self.nt_fasta_gz = self.test_data_dir.joinpath("nt_all.fasta.gz")

        # Verify that all tests FASTA files exist
        for test_file in (self.nt_fasta1, self.nt_fasta2, self.nt_fasta3,
                          self.nt_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define a temporary directory for MMseqs2 outputs
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")

        # Define cluster output file names
        self.cluster_prefix = self.tmpdir.joinpath("test_easy_linclust_nt")
        suffices = ("_cluster.tsv", "_all_seqs.fasta", "_rep_seq.fasta")
        self.cluster_outfiles = (pathlib.Path(str(self.cluster_prefix) + s) for
                                 s in suffices)

    def test_one_query(self):
        """Verify that we can run the 'easy-linclust' pipeline with a single
        nucleotide query file."""
        stdout = self.mmseqs.easy_linclust(self.nt_fasta1, self.cluster_prefix,
                                           self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_two_query(self):
        """Verify that we can run the 'easy-linclust' pipeline with two
        nucleotide query files."""
        stdout = self.mmseqs.easy_linclust([self.nt_fasta1, self.nt_fasta2],
                                           self.cluster_prefix, self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_gzipped_query(self):
        """Verify that we can run the 'easy-linclust' pipeline with a gzipped
        nucleotide query file."""
        stdout = self.mmseqs.easy_linclust(
            self.nt_fasta_gz, self.cluster_prefix, self.tmpdir)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_valid_kwargs(self):
        """Verify that we can run the 'easy-linclust' pipeline with a single
        nucleotide query file and valid keyword arguments."""
        stdout = self.mmseqs.easy_linclust(self.nt_fasta1, self.cluster_prefix,
                                           self.tmpdir, min_seq_id=0.9,
                                           e=0.01, c=0.75)

        self.assertIsInstance(stdout, str)
        for outfile in self.cluster_outfiles:
            self.assertTrue(outfile.exists())

    def test_invalid_kwargs(self):
        """Verify that the 'easy-linclust' pipeline with a single nucleotide
        query file and invalid keyword arguments raises a KeyError."""
        with self.assertRaises(KeyError):
            self.mmseqs.easy_linclust(self.nt_fasta1, self.cluster_prefix,
                                      self.tmpdir, max_seq_id=0.9, evalue=0.01,
                                      coverage=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def test_valid_kw_bad_type(self):
        """Verify that the 'easy-linclust' pipeline with a valid keyword
        argument but the wrong type raises a TypeError."""
        with self.assertRaises(TypeError):
            self.mmseqs.easy_linclust(self.nt_fasta1, self.cluster_prefix,
                                      self.tmpdir, min_seq_id="0.9",
                                      e=0.01, c=0.75)

        for outfile in self.cluster_outfiles:
            self.assertFalse(outfile.exists())

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


class TestMMseqs2(unittest.TestCase):
    def setUp(self):
        """Set up the testing context."""
        self.mmseqs = MMseqs2()

        # Test data files are stored here
        self.test_data_dir = pathlib.Path(__file__).resolve().parent.joinpath(
            "data")

        # We have several amino acid and nucleotide FASTA files to tests with
        self.aa_fasta1 = self.test_data_dir.joinpath("aa_1.fasta")
        self.aa_fasta2 = self.test_data_dir.joinpath("aa_2.fasta")
        self.aa_fasta3 = self.test_data_dir.joinpath("aa_3.fasta")
        self.nt_fasta1 = self.test_data_dir.joinpath("nt_1.fasta")
        self.nt_fasta2 = self.test_data_dir.joinpath("nt_2.fasta")
        self.nt_fasta3 = self.test_data_dir.joinpath("nt_3.fasta")

        # Verify that all tests FASTA files exist
        for test_file in (self.aa_fasta1, self.aa_fasta2, self.aa_fasta3,
                          self.nt_fasta1, self.nt_fasta2, self.nt_fasta3):
            self.assertTrue(test_file.exists())

        # We also have a pair of gzip-compressed FASTA files
        self.aa_fasta_gz = self.test_data_dir.joinpath("aa_all.fasta.gz")
        self.nt_fasta_gz = self.test_data_dir.joinpath("nt_all.fasta.gz")

        # Verify that the tests compressed FASTA files exist
        for test_file in (self.aa_fasta_gz, self.nt_fasta_gz):
            self.assertTrue(test_file.exists())

        # Define and create a temporary directory for MMseqs2 output
        self.tmpdir = self.test_data_dir.parent.joinpath("mmseqs_tmp")
        # self.tmpdir.mkdir(parents=True, exist_ok=True)

        # Define search output files
        self.search_nt1_nt2 = self.tmpdir.joinpath("search_nt1_nt2.tsv")
        self.search_nt1_nt2_nt3 = self.tmpdir.joinpath("search_nt1_nt2_nt3.tsv")
        self.search_aa1_aa2 = self.tmpdir.joinpath("search_aa1_aa2.tsv")
        self.search_aa1_aa2_aa3 = self.tmpdir.joinpath("search_aa1_aa2_aa3.tsv")

        # Define cluster output files
        self.cluster_aa1_aa2 = self.tmpdir.joinpath("cluster_aa1_aa2")
        self.cluster_aa1_aa2_aa3 = self.tmpdir.joinpath("cluster_aa1_aa2_aa3")

    def test_easy_taxonomy(self):
        self.fail()

    def test_easy_rbh(self):
        self.fail()

    def test_search(self):
        self.fail()

    def test_linsearch(self):
        self.fail()

    def test_map(self):
        self.fail()

    def test_rbh(self):
        self.fail()

    def test_linclust(self):
        self.fail()

    def test_cluster(self):
        self.fail()

    def test_clusterupdate(self):
        self.fail()

    def test_taxonomy(self):
        self.fail()

    def test_databases(self):
        self.fail()

    def test_createdb(self):
        self.fail()

    def test_createindex(self):
        self.fail()

    def test_createlinindex(self):
        self.fail()

    def test_convertmsa(self):
        self.fail()

    def test_tsv2db(self):
        self.fail()

    def test_tar2db(self):
        self.fail()

    def test_db2tar(self):
        self.fail()

    def test_msa2profile(self):
        self.fail()

    def test_compress(self):
        self.fail()

    def test_decompress(self):
        self.fail()

    def test_rmdb(self):
        self.fail()

    def test_mvdb(self):
        self.fail()

    def test_cpdb(self):
        self.fail()

    def test_lndb(self):
        self.fail()

    def test_aliasdb(self):
        self.fail()

    def test_unpackdb(self):
        self.fail()

    def test_touchdb(self):
        self.fail()

    def test_createsubdb(self):
        self.fail()

    def test_concatdbs(self):
        self.fail()

    def test_splitdb(self):
        self.fail()

    def test_mergedbs(self):
        self.fail()

    def test_subtractdbs(self):
        self.fail()

    def test_convertalis(self):
        self.fail()

    def test_createtsv(self):
        self.fail()

    def test_convert2fasta(self):
        self.fail()

    def test_result2flat(self):
        self.fail()

    def test_createseqfiledb(self):
        self.fail()

    def test_taxonomyreport(self):
        self.fail()

    def test_extractorfs(self):
        self.fail()

    def test_extractframes(self):
        self.fail()

    def test_orftocontig(self):
        self.fail()

    def test_reverseseq(self):
        self.fail()

    def test_translatenucs(self):
        self.fail()

    def test_translateaa(self):
        self.fail()

    def test_splitsequence(self):
        self.fail()

    def test_masksequence(self):
        self.fail()

    def test_extractalignedregion(self):
        self.fail()

    def test_swapresults(self):
        self.fail()

    def test_result2rbh(self):
        self.fail()

    def test_result2msa(self):
        self.fail()

    def test_result2dnamsa(self):
        self.fail()

    def test_result2stats(self):
        self.fail()

    def test_filterresult(self):
        self.fail()

    def test_offsetalignment(self):
        self.fail()

    def test_proteinaln2nucl(self):
        self.fail()

    def test_result2repseq(self):
        self.fail()

    def test_sortresult(self):
        self.fail()

    def test_summarizealis(self):
        self.fail()

    def test_summarizeresult(self):
        self.fail()

    def test_createtaxdb(self):
        self.fail()

    def test_createbintaxonomy(self):
        self.fail()

    def test_createbintaxmapping(self):
        self.fail()

    def test_addtaxonomy(self):
        self.fail()

    def test_filtertaxdb(self):
        self.fail()

    def test_filtertaxseqdb(self):
        self.fail()

    def test_aggregatetax(self):
        self.fail()

    def test_aggregatetaxweights(self):
        self.fail()

    def test_lcaalign(self):
        self.fail()

    def test_lca(self):
        self.fail()

    def test_majoritylca(self):
        self.fail()

    def test_multihitdb(self):
        self.fail()

    def test_multihitsearch(self):
        self.fail()

    def test_besthitperset(self):
        self.fail()

    def test_combinepvalperset(self):
        self.fail()

    def test_mergeresultsbyset(self):
        self.fail()

    def test_prefilter(self):
        self.fail()

    def test_ungappedprefilter(self):
        self.fail()

    def test_gappedprefilter(self):
        self.fail()

    def test_kmermatcher(self):
        self.fail()

    def test_kmersearch(self):
        self.fail()

    def test_align(self):
        self.fail()

    def test_alignall(self):
        self.fail()

    def test_transitivealign(self):
        self.fail()

    def test_rescorediagonal(self):
        self.fail()

    def test_alignbykmer(self):
        self.fail()

    def test_clust(self):
        self.fail()

    def test_clusthash(self):
        self.fail()

    def test_mergeclusters(self):
        self.fail()

    def test_result2profile(self):
        self.fail()

    def test_msa2result(self):
        self.fail()

    def test_sequence2profile(self):
        self.fail()

    def test_profile2pssm(self):
        self.fail()

    def test_profile2neff(self):
        self.fail()

    def test_profile2consensus(self):
        self.fail()

    def test_profile2repseq(self):
        self.fail()

    def test_convertprofiledb(self):
        self.fail()

    def test_enrich(self):
        self.fail()

    def test_result2pp(self):
        self.fail()

    def test_profile2cs(self):
        self.fail()

    def test_tsv2exprofiledb(self):
        self.fail()

    def test_convertca3m(self):
        self.fail()

    def test_expandaln(self):
        self.fail()

    def test_expand2profile(self):
        self.fail()

    def test_setextendeddbtype(self):
        self.fail()

    def test_view(self):
        self.fail()

    def test_apply(self):
        self.fail()

    def test_filterdb(self):
        self.fail()

    def test_swapdb(self):
        self.fail()

    def test_prefixid(self):
        self.fail()

    def test_suffixid(self):
        self.fail()

    def test_renamedbkeys(self):
        self.fail()

    def test_diffseqdbs(self):
        self.fail()

    def test_summarizetabs(self):
        self.fail()

    def test_gff2db(self):
        self.fail()

    def test_maskbygff(self):
        self.fail()

    def test_convertkb(self):
        self.fail()

    def test_summarizeheaders(self):
        self.fail()

    def test_nrtotaxmapping(self):
        self.fail()

    def test_extractdomains(self):
        self.fail()

    def test_countkmer(self):
        self.fail()

    def test_help(self):
        """Verify that invoking the help method with no arguments returns a
        string."""
        self.assertIsInstance(self.mmseqs.help(), str)

    def test_help_valid_pipeline(self):
        """Verify that invoking the help method with a valid pipeline returns
        a string."""
        self.assertIsInstance(self.mmseqs.help("easy-search"), str)

    def test_help_invalid_pipeline(self):
        """Verify that invoking the help method with an invalid pipeline raises
        an MMseqs2Error."""
        with self.assertRaises(MMseqs2Error):
            self.mmseqs.help("foo")

    def test_version(self):
        """Test that the version method returns a string."""
        self.assertIsInstance(self.mmseqs.version(), str)

    def test_version_matches_api_version(self):
        """Test that the version method and api_version match."""
        self.assertEqual(self.mmseqs.version().split(".")[0],
                         str(self.mmseqs.api_version))

    def tearDown(self):
        """Clean up the testing context."""
        if self.tmpdir.exists():
            shutil.rmtree(self.tmpdir)


if __name__ == '__main__':
    unittest.main()
