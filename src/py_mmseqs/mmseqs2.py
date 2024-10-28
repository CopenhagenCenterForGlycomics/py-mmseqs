"""A Python wrapper around MMseqs2."""

import json
import logging
import os
import pathlib
import subprocess as sp


def _cast_list_objs_to_paths(path_list):
    """Cast a list of objects to fully resolved pathlib.Path objects.

    :param path_list: a list of objects to cast to pathlib.Path
    :type path_list: list[os.PathLike]
    """
    for i, path in enumerate(path_list):
        try:
            path_list[i] = pathlib.Path(path).expanduser().resolve()
        except TypeError:
            raise TypeError(f"unable to cast '{path}' to pathlib.Path")

    return path_list


class MMseqs2Error(Exception):
    """Raised when MMseqs2 pipelines finish with non-zero exit codes."""
    pass


class MMseqs2:
    """This class is used to invoke MMseqs2 pipelines from Python.

    The API is designed to be as close as possible to the MMseqs2
    command-line interface, with each pipeline implemented as a class
    method. It should be perfectly compatible with the following major
    releases of MMseqs2: 13.45111, 14.7e284, 15.6f452.

    Every pipeline has some number of required positional arguments,
    followed by optional keyword arguments. Optional arguments can be
    passed as keyword arguments with leading dashes removed and internal
    dashes replaced with underscores (e.g., -e 0.5 becomes e=0.5;
    --min-seq-id 0.3 becomes min_seq_id=0.3).

    For help with a specific pipeline, call the help() method with the
    name of the pipeline as an argument. For general help, call help()
    with no arguments. The version() method returns the version string
    from MMseqs2.

    If MMseqs2 is not globally available (i.e., not in the PATH), the
    caller must supply the path to the MMseqs2 binary when creating an
    instance of this class. The API version is determined using the
    version() method, and is used to validate the keyword arguments
    passed to each pipeline. If a keyword argument is not valid for the
    specified pipeline, a KeyError will be raised. If the value passed
    with a valid keyword argument is not of the correct type, a
    TypeError will be raised.

    If a formatted MMseqs2 command returns a non-zero exit code, an
    MMseqs2Error will be raised, with the MMseqs2 error message
    provided. Otherwise, the stdout of the command will be returned as
    a string to the caller."""
    def __init__(self, binary_path='mmseqs'):
        """Initialize the MMseqs2 wrapper.

        :param binary_path: path to the MMseqs2 binary (only needed if
            not in PATH)
        :type binary_path: str | os.PathLike
        :raises FileNotFoundError: if the MMseqs2 binary is not found
        :raises MMseqs2Error: if the installed MMseqs2 version is not
            supported
        """
        self.binary_path = binary_path
        try:
            # Conveniently, the 'mmseqs version' command is the same for all
            # recent versions of MMseqs2 (maybe all versions?)
            self.api_version = int(self.version().split('-')[0].split(".")[0])
        except FileNotFoundError:
            raise FileNotFoundError("MMseqs2 binary not found; please ensure "
                                    "it is installed and in your PATH, or "
                                    "provide the path when creating an MMseqs2 "
                                    "instance.")

        if not 13 <= self.api_version <= 15:
            raise MMseqs2Error(f"installed MMseqs2 version {self.api_version} "
                               f"is not supported; please install one of the "
                               f"major releases 13, 14, or 15.")

        # Find the path to api_XX_argtypes.json
        whereami = pathlib.Path(__file__).resolve().parent
        api_file = whereami / f"api_{self.api_version}_argtypes.json"
        with open(api_file, 'r') as json_file:
            self.arg_types = json.load(json_file)

    # Easy workflows for plain text input/output
    def easy_search(self, query, target, output, tmpdir, **kwargs):
        """Sensitive homology search.

        Query sequence file(s) are searched against a target sequence
        file, with results written to an output file in TSV format.

        Note: for some keyword arguments (e.g., min_seq_id or c), the
        supplied value must be in the range [0, 1]. The responsibility
        is left to the caller to ensure that supplied values are within
        the correct range. If a value outside the range is supplied, an
        error will be raised by MMseqs2.

        :param query: path to query sequence file(s) in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike | list[os.PathLike]
        :param target: path to target sequence file in FASTA/Q format
        :type target: os.PathLike
        :param output: path to desired output file (TSV)
        :type output: os.PathLike
        :param tmpdir: path to temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises TypeError: if query is not an instance of os.PathLike or a
            list of os.PathLike instances
        :raises ValueError: if query is a list and any element cannot be cast
            to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        # If query not a list, make one; all positional args must be PathLike
        if not isinstance(query, list):
            args = [query]
        else:
            args = query
        args.extend([target, output, tmpdir])

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('easy-search', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('easy-search', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def easy_linsearch(self, query, target, output, tmpdir, **kwargs):
        """Fast (linear time), less sensitive homology search.

        Query sequence file(s) are searched against a target sequence
        file, with results written to an output file in TSV format.

        Note: for some keyword arguments (e.g., min_seq_id or c), the
        supplied value must be in the range [0, 1]. The responsibility
        is left to the caller to ensure that supplied values are within
        the correct range. If a value outside the range is supplied, an
        error will be raised by MMseqs2.

        Note: for mixed-type searches with amino-acid query and nucleotide
        target files, you must use 'search_type=2' to avoid an error.

        :param query: path to query sequence file(s) in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike | list[os.PathLike]
        :param target: path to target sequence file in FASTA/Q format
        :type target: os.PathLike
        :param output: path to desired output file (TSV)
        :type output: os.PathLike
        :param tmpdir: path to temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises TypeError: if query is not an instance of os.PathLike or a
            list of os.PathLike instances
        :raises ValueError: if query is a list and any element cannot be cast
            to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        # @milot-mirdita suggests this pipeline is experimental
        raise NotImplementedError("easy-linsearch pipeline is experimental "
                                  "and not supported")

        # If query not a list, make one; all positional args must be PathLike
        if not isinstance(query, list):
            args = [query]
        else:
            args = query
        args.extend([target, output, tmpdir])

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('easy-linsearch', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('easy-linsearch', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def easy_cluster(self, query, output_prefix, tmpdir, **kwargs):
        """Sensitive clustering.

        Sequences from one or more query files are pooled and clustered, with
        results written to three output files:

          - A two-column TSV file with the name <output_prefix>_cluster.tsv
          - A FASTA file with the name <output_prefix>_rep_seq.fasta
          - A FASTA file with the name <output_prefix>_all_seqs.fasta

        The first file maps cluster representative IDs in the first column
        to cluster member IDs in the second column (one pair per line).

        The second file contains the representative sequences for each cluster,
        and the third file contains all clustered sequences in a FASTA-like
        format.

        Note: for some keyword arguments (e.g., min_seq_id or c), the
        supplied value must be in the range [0, 1]. The responsibility
        is left to the caller to ensure that supplied values are within
        the correct range. If a value outside the range is supplied, an
        error will be raised by MMseqs2.

        Clusters are most easily generated and extracted for use by doing
        something like:

        >>> from py_mmseqs.parsers import parse_clusters
        >>> queries = [pathlib.Path("seqs1.fasta"), pathlib.Path("seqs2.fasta")]
        >>> prefix = pathlib.Path("easy_clusters")
        >>> tmp = pathlib.Path("tmp")
        >>> mmseqs = MMseqs2()
        >>> stdout = mmseqs.easy_cluster(queries, prefix, tmp)
        >>> clu_out = pathlib.Path(str(prefix) + "_all_seqs.fasta")
        >>> clusters = list(parse_clusters(clu_out))

        :param query: path to one or more sequence files in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike | list[os.PathLike]
        :param output_prefix: desired output file path prefix (outfile will
            be a two-column TSV file with name <prefix>_cluster.tsv)
        :type output_prefix: os.PathLike
        :param tmpdir: path to temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises TypeError: if query is not an instance of os.PathLike or a
            list of os.PathLike instances
        :raises ValueError: if query is a list and any element cannot be cast
            to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        # If query not a list, make one; all positional args must be PathLike
        if not isinstance(query, list):
            args = [query]
        else:
            args = query
        args.extend([output_prefix, tmpdir])

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('easy-cluster', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('easy-cluster', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def easy_linclust(self, query, output_prefix, tmpdir, **kwargs):
        """Fast (linear time), less sensitive clustering.

        Sequences from one or more query files are pooled and clustered, with
        results written to three output files:

          - A two-column TSV file with the name <output_prefix>_cluster.tsv
          - A FASTA file with the name <output_prefix>_rep_seq.fasta
          - A FASTA file with the name <output_prefix>_all_seqs.fasta

        The first file maps cluster representative IDs in the first column
        to cluster member IDs in the second column (one pair per line).

        The second file contains the representative sequences for each cluster,
        and the third file contains all clustered sequences in a FASTA-like
        format.

        Note: for some keyword arguments (e.g., min_seq_id or c), the
        supplied value must be in the range [0, 1]. The responsibility
        is left to the caller to ensure that supplied values are within
        the correct range. If a value outside the range is supplied, an
        error will be raised by MMseqs2.

        Clusters are most easily generated and extracted for use by doing
        something like:

        >>> from py_mmseqs.parsers import parse_clusters
        >>> queries = [pathlib.Path("seqs1.fasta"), pathlib.Path("seqs2.fasta")]
        >>> prefix = pathlib.Path("easy_clusters")
        >>> tmp = pathlib.Path("tmp")
        >>> mmseqs = MMseqs2()
        >>> stdout = mmseqs.easy_linclust(queries, prefix, tmp)
        >>> clu_out = pathlib.Path(str(prefix) + "_all_seqs.fasta")
        >>> clusters = list(parse_clusters(clu_out))

        :param query: path to one or more sequence files in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike | list[os.PathLike]
        :param output_prefix: desired output file path prefix (outfile will
            be a two-column TSV file with name <prefix>_cluster.tsv)
        :type output_prefix: os.PathLike
        :param tmpdir: path to temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises TypeError: if query is not an instance of os.PathLike or a
            list of os.PathLike instances
        :raises ValueError: if query is a list and any element cannot be cast
            to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        # If query not a list, make one; all positional args must be PathLike
        if not isinstance(query, list):
            args = [query]
        else:
            args = query
        args.extend([output_prefix, tmpdir])

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('easy-linclust', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('easy-linclust', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def easy_taxonomy(self, query, target_db, output_prefix, tmpdir, **kwargs):
        """Taxonomic classification.

        :param query: path to query sequence file(s) in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike | list[os.PathLike]
        :param target_db: path to the target sequence database generated with
            `createdb` pipeline AND annotated with taxonomy information
        :type target_db: os.PathLike
        :param output_prefix: desired output file path prefix (several output
            files will be generated with this prefix)
        :type output_prefix: os.PathLike
        :param tmpdir: path to temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises TypeError: if query is not an instance of os.PathLike or a
            list of os.PathLike instances
        :raises ValueError: if query is a list and any element cannot be cast
            to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        # If query not a list, make one; all positional args must be PathLike
        if not isinstance(query, list):
            args = [query]
        else:
            args = query
        args.extend([target_db, output_prefix, tmpdir])

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('easy-taxonomy', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('easy-taxonomy', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def easy_rbh(self, query, target, output, tmpdir, **kwargs):
        """Find reciprocal best hits between .

        :param query: path to a query sequence file in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike
        :param target: path to a target sequence file in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type target: os.PathLike
        :param output: path to desired output file (TSV)
        :type output: os.PathLike
        :param tmpdir: path to temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query, target, output, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('easy-rbh', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('easy-rbh', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Main workflows for database input/output
    def search(self, query_db, target_db, alignment_db, tmpdir, **kwargs):
        """Sensitive homology search.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param target_db: path to a target sequence database generated with
            `createdb` pipeline
        :type target_db: os.PathLike
        :param alignment_db: path to the output alignment database to be
            generated
        :type alignment_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, alignment_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('search', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('search', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def linsearch(self, query_db, target_db, alignment_db, tmpdir, **kwargs):
        """Fast, less sensitive homology search.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param target_db: path to a target sequence database generated with
            `createdb` pipeline
        :type target_db: os.PathLike
        :param alignment_db: path to the output alignment database to be
            generated
        :type alignment_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, alignment_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('linsearch', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('linsearch', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def map(self, query_db, target_db, alignment_db, tmpdir, **kwargs):
        """Map nearly identical sequences.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param target_db: path to a target sequence database generated with
            `createdb` pipeline
        :type target_db: os.PathLike
        :param alignment_db: path to the output alignment database to be
            generated
        :type alignment_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, alignment_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('map', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('map', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def rbh(self, query_db, target_db, alignment_db, tmpdir, **kwargs):
        """Reciprocal best hit search.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param target_db: path to a target sequence database generated with
            `createdb` pipeline
        :type target_db: os.PathLike
        :param alignment_db: path to the output alignment database to be
            generated
        :type alignment_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, alignment_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('rbh', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('rbh', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def linclust(self, query_db, cluster_db, tmpdir, **kwargs):
        """Fast, less sensitive clustering.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param cluster_db: path to the output cluster database to be
            generated
        :type cluster_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, cluster_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('linclust', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('linclust', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def cluster(self, query_db, cluster_db, tmpdir, **kwargs):
        """Slower, sensitive clustering.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param cluster_db: path to the output cluster database to be
            generated
        :type cluster_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, cluster_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('cluster', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('cluster', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def clusterupdate(self, old_query_db, new_query_db, old_cluster_db,
                      new_mapped_db, new_cluster_db, tmpdir, **kwargs):
        """Update previous clustering with new sequences.

        :param old_query_db: path to the original query sequence database
            generated with `createdb` pipeline
        :type old_query_db: os.PathLike
        :param new_query_db: path to the new query sequence database
            generated with `createdb` pipeline (sequences have been added
            or removed from the original database)
        :type new_query_db: os.PathLike
        :param old_cluster_db: path to the original cluster database
            generated with `cluster` or `linclust` pipeline
        :type old_cluster_db: os.PathLike
        :param new_mapped_db: path to a new sequence database to be generated
            (will contain all unique sequences from the original and new query
            databases and can be used as `old_query_db` in a subsequent call)
        :type new_mapped_db: os.PathLike
        :param new_cluster_db: path to the output cluster database to be
            generated
        :type new_cluster_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [old_query_db, new_query_db, old_cluster_db, new_mapped_db,
                new_cluster_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('clusterupdate', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('clusterupdate', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def taxonomy(self, query_db, target_db, taxa_db, tmpdir, **kwargs):
        """Taxonomic classification.

        :param query_db: path to a query sequence database generated with
            `createdb` pipeline
        :type query_db: os.PathLike
        :param target_db: path to a target sequence database generated with
            `createdb` pipeline AND annotated with taxonomy information
        :type target_db: os.PathLike
        :param taxa_db: path to the output mapped taxa database to be generated
        :type taxa_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, taxa_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('taxonomy', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('taxonomy', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Input database creation
    def databases(self, name, target_db, tmpdir):
        """List and download hosted databases.

        :param name: name of the database to download
        :type name: str
        :param target_db: path to the target database to be created
        :type target_db: os.PathLike
        :param tmpdir: path to a temporary directory to use
        :type tmpdir: os.PathLike
        """
        args = [name, target_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('databases', *args)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createdb(self, query, target_db, **kwargs):
        """Convert FASTA/Q file(s) to a sequence DB.

        :param query: path to query sequence file(s) in FASTA/Q format,
            optionally compressed with gzip or bzip2
        :type query: os.PathLike | list[os.PathLike]
        :param target_db: path to the target sequence database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises TypeError: if query is not an instance of os.PathLike or a
            list of os.PathLike instances
        :raises ValueError: if query is a list and any element cannot be cast
            to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        # If query not a list, make one; all positional args must be PathLike
        if not isinstance(query, list):
            args = [query]
        else:
            args = query

        args.append(target_db)

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createindex(self, target_db, tmpdir, **kwargs):
        """Store precomputed index on disk to reduce search overhead.

        :param target_db: path to the target sequence database to be indexed
        :type target_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [target_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createindex', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createindex', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createlinindex(self, target_db, tmpdir, **kwargs):
        """Create linsearch index.

        :param target_db: path to the target sequence database to be indexed
        :type target_db: os.PathLike
        :param tmpdir: path to a temporary directory
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [target_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createlinindex', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createlinindex', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def convertmsa(self, msa_file, target_db, **kwargs):
        """Convert Stockholm/PFAM MSA file to an MSA DB.

        :param msa_file: path to the MSA file in Stockholm/PFAM format,
            optionally compressed with gzip
        :type msa_file: os.PathLike
        :param target_db: path to the target MSA database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [msa_file, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('convertmsa', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('convertmsa', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def tsv2db(self, tsv_file, target_db, **kwargs):
        """Convert a TSV file to any DB.

        :param tsv_file: path to the TSV file to convert to a database
        :type tsv_file: os.PathLike
        :param target_db: path to the target database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [tsv_file, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('tsv2db', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('tsv2db', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def tar2db(self, tar_file, target_db, **kwargs):
        """Convert content of tar archives to any DB.

        :param tar_file: path to the tar archive file to convert to a database,
            optionally compressed with gzip
        :type tar_file: os.PathLike
        :param target_db: path to the target database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [tar_file, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('tar2db', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('tar2db', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def db2tar(self, db, tar_file, **kwargs):
        """Archive contents of a DB to a tar archive.

        :param db: path to the database to be archived
        :type db: os.PathLike
        :param tar_file: path to the tar archive file to create, optionally
            compressed with gzip
        :type tar_file: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        if self.api_version < 14:
            raise MMseqs2Error("'db2tar' was not introduced until MMseqs2 "
                               "version 14")

        args = [db, tar_file]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('db2tar', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('db2tar', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def msa2profile(self, msa_db, profile_db, **kwargs):
        """Convert an MSA DB to a profile DB.

        :param msa_db: path to the MSA database to convert to a profile
            database
        :type msa_db: os.PathLike
        :param profile_db: path to the target profile database to be created
        :type profile_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [msa_db, profile_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('msa2profile', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('msa2profile', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Handle databases on storage and memory
    def compress(self, db, compressed_db, **kwargs):
        """Compress DB entries.

        :param db: path to the database to compress
        :type db: os.PathLike
        :param compressed_db: path to the compressed database to be created
        :type compressed_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [db, compressed_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('compress', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('compress', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def decompress(self, compressed_db, db, **kwargs):
        """Decompress DB entries.

        :param compressed_db: path to the compressed database to decompress
        :type compressed_db: os.PathLike
        :param db: path to the decompressed database to be created
        :type db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [compressed_db, db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('decompress', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('decompress', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def rmdb(self, db, **kwargs):
        """Remove a DB.

        :param db: path to the database to remove
        :type db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if db cannot be cast to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('rmdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('rmdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def mvdb(self, source_db, target_db, **kwargs):
        """Move a DB.

        :param source_db: path to the source database to move
        :type source_db: os.PathLike
        :param target_db: path to which the database should be moved
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('mvdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('mvdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def cpdb(self, source_db, target_db, **kwargs):
        """Copy a DB.

        :param source_db: path to the source database to copy
        :type source_db: os.PathLike
        :param target_db: path to which the database should be copied
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('cpdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('cpdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def lndb(self, source_db, target_db, **kwargs):
        """Symlink a DB.

        :param source_db: path to the source database to symlink
        :type source_db: os.PathLike
        :param target_db: path to which the database should be symlinked
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('lndb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('lndb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def aliasdb(self, source_db, target_db, **kwargs):
        """Create relative symlink of DB to another name in the same folder.

        :param source_db: path to the source database to alias
        :type source_db: os.PathLike
        :param target_db: path to which the database should be aliased
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        if self.api_version < 14:
            raise MMseqs2Error("'aliasdb' was not introduced until MMseqs2 "
                               "version 14")

        args = [source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('aliasdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('aliasdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def unpackdb(self, db, outdir, **kwargs):
        """Unpack a DB into separate files.

        :param db: path to the database to unpack
        :type db: os.PathLike
        :param outdir: path to the directory in which to unpack the database
        :type outdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [db, outdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('unpackdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('unpackdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def touchdb(self, db, **kwargs):
        """Preload DB into memory (page cache).

        :param db: path to the database to preload
        :type db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if db cannot be cast to a pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('touchdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('touchdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Unite and intersect databases
    def createsubdb(self, subset_file, source_db, target_db, **kwargs):
        """Create a subset of a DB from list of DB keys.

        :param subset_file: path to the file containing the subset of keys
            to extract from the source database
        :type subset_file: os.PathLike
        :param source_db: path to the source database from which to extract
            the subset of keys
        :type source_db: os.PathLike
        :param target_db: path to the target database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [subset_file, source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createsubdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createsubdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def concatdbs(self, db1, db2, target_db, **kwargs):
        """Concatenate two DBs, giving new IDs to entries from 2nd DB.

        :param db1: path to the first database to be concatenated
        :type db1: os.PathLike
        :param db2: path to the second database to be concatenated
        :type db2: os.PathLike
        :param target_db: path to the target database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [db1, db2, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('concatdbs', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('concatdbs', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def splitdb(self, source_db, target_db, **kwargs):
        """Split DB into subsets.

        :param source_db: path to the source database to split
        :type source_db: os.PathLike
        :param target_db: path to the target database to be created
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('splitdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('splitdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def mergedbs(self):
        """Merge entries from multiple DBs."""
        # TODO: Implement mergedbs after determining the correct signature
        raise NotImplementedError("mergedbs is not yet implemented")

    def subtractdbs(self, result_db_left, result_db_right, result_db, **kwargs):
        """Remove all entries from first DB occurring in second DB by
        key.

        :param result_db_left: path to a result database to subtract hits from
        :type result_db_left: os.PathLike
        :param result_db_right: path to a result database to use as a filter
        :type result_db_right: os.PathLike
        :param result_db: path to the target result database to be created
        :type result_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [result_db_left, result_db_right, result_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('subtractdbs', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('subtractdbs', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Format conversion for downstream processing
    def convertalis(self, query_db, target_db, alignment_db, alignment_file,
                    **kwargs):
        """Convert alignment DB to BLAST-tab, SAM or custom format.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param alignment_db: path to the alignment database to be converted
        :type alignment_db: os.PathLike
        :param alignment_file: path to the output alignment file to create
        :type alignment_file: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, alignment_db, alignment_file]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('convertalis', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('convertalis', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createtsv(self, query_db, target_db, result_db, tsv_file, **kwargs):
        """Convert result DB to tab-separated flat file.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to convert to TSV
        :type result_db: os.PathLike
        :param tsv_file: path to the output TSV file to create
        :type tsv_file: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, tsv_file]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createtsv', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createtsv', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def convert2fasta(self, sequence_db, fasta_file, **kwargs):
        """Convert sequence DB to FASTA format.

        :param sequence_db: path to the sequence database to convert to FASTA
        :type sequence_db: os.PathLike
        :param fasta_file: path to the output FASTA file to create
        :type fasta_file: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, fasta_file]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('convert2fasta', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('convert2fasta', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def result2flat(self, query_db, target_db, result_db, fasta_db, **kwargs):
        """Create flat file by adding FASTA headers to DB entries.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to convert to flat
        :type result_db: os.PathLike
        :param fasta_db: path to the FASTA database to use for headers
        :type fasta_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, fasta_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('result2flat', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('result2flat', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createseqfiledb(self, sequence_db, result_db, fasta_db, **kwargs):
        """Create a DB of unaligned FASTA entries.

        :param sequence_db: path to the sequence database to convert to FASTA
        :type sequence_db: os.PathLike
        :param result_db: path to the result database to read clusters from
        :type result_db: os.PathLike
        :param fasta_db: path to the unaligned FASTA database to create
        :type fasta_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, result_db, fasta_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createseqfiledb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createseqfiledb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def taxonomyreport(self, sequence_db, result_db, report_file, **kwargs):
        """Create a taxonomy report in Kraken or Krona format.

        :param sequence_db: path to a taxonomy-labeled sequence database
        :type sequence_db: os.PathLike
        :param result_db: path to the result database to use for taxonomy
            assignment
        :type result_db: os.PathLike
        :param report_file: path to the output report file to create
        :type report_file: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, result_db, report_file]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('taxonomyreport', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('taxonomyreport', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Sequence manipulation/transformation
    def extractorfs(self, sequence_db, orf_db, **kwargs):
        """Six-frame extraction of open reading frames.

        :param sequence_db: path to the sequence database to extract ORFs from
        :type sequence_db: os.PathLike
        :param orf_db: path to the output ORF database to create
        :type orf_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, orf_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('extractorfs', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('extractorfs', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def extractframes(self, sequence_db, frame_db, **kwargs):
        """Extract frames from a nucleotide sequence DB.

        :param sequence_db: path to the sequence database to extract frames from
        :type sequence_db: os.PathLike
        :param frame_db: path to the output frame database to create
        :type frame_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, frame_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('extractframes', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('extractframes', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def orftocontig(self, contig_db, orf_db, alignment_db, **kwargs):
        """Write ORF locations in alignment format.

        :param contig_db: path to a contig sequence database to map orfs against
        :type contig_db: os.PathLike
        :param orf_db: path to the ORF sequence database to map against contigs
        :type orf_db: os.PathLike
        :param alignment_db: path to the output alignment database to create
        :type alignment_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [contig_db, orf_db, alignment_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('orftocontig', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('orftocontig', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def reverseseq(self, sequence_db, rev_seq_db, **kwargs):
        """Reverse (without complement) sequences.

        :param sequence_db: path to the sequence database whose sequences are
            to be reversed
        :type sequence_db: os.PathLike
        :param rev_seq_db: path to the reversed sequence database to create
        :type rev_seq_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, rev_seq_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('reverseseq', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('reverseseq', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def translatenucs(self, sequence_db, translate_db, **kwargs):
        """Translate nucleotides to proteins.

        :param sequence_db: path to the sequence database whose nucleotides
            are to be translated
        :type sequence_db: os.PathLike
        :param translate_db: path to the translated protein sequence database
            to create
        :type translate_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, translate_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('translatenucs', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('translatenucs', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def translateaa(self, sequence_db, translate_db, **kwargs):
        """Translate proteins to lexicographically lowest codons.

        NOTE: this is not round-trip safe with `translatenucs`, as use
        of lexicographically lowest codons may result in a different
        nucleotide sequence than the original.
        """
        # TODO: Implement translateaa after determining the correct signature
        raise NotImplementedError("translateaa is not yet implemented")

    def splitsequence(self, source_db, target_db, **kwargs):
        """Split sequences by length.

        :param source_db: path to the source sequence database whose
            sequences are to be split on length
        :type source_db: os.PathLike
        :param target_db: path to the target split-sequence database to create
        :type target_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [source_db, target_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('splitsequence', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('splitsequence', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def masksequence(self, sequence_db, mask_db, **kwargs):
        """Soft mask sequence DB using tantan.

        :param sequence_db: path to the sequence database to mask
        :type sequence_db: os.PathLike
        :param mask_db: path to the masked sequence database to create
        :type mask_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, mask_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('masksequence', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('masksequence', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def extractalignedregion(self, query_db, target_db, result_db,
                             sequence_db, **kwargs):
        """Extract aligned sequence region from query database.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to extract aligned regions
            from
        :type result_db: os.PathLike
        :param sequence_db: path to the output sequence database to create
        :type sequence_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, sequence_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('extractalignedregion', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('extractalignedregion', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Result manipulation
    def swapresults(self, query_db, target_db, result_db, swapped_db, **kwargs):
        """Transpose prefilter/alignment DB.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to swap query/target
            data for
        :type result_db: os.PathLike
        :param swapped_db: path to the output swapped result database to create
        :type swapped_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, swapped_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('swapresults', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('swapresults', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def result2rbh(self, result_db, rbh_db, **kwargs):
        """Filter a merged result DB to retain only reciprocal best
        hits.

        :param result_db: path to the result database to filter for RBHs
        :type result_db: os.PathLike
        :param rbh_db: path to the filtered RBH result database to create
        :type rbh_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [result_db, rbh_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('result2rbh', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('result2rbh', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def result2msa(self, query_db, target_db, result_db, msa_db, **kwargs):
        """Compute MSA DB from a result DB.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to compute MSAs from
        :type result_db: os.PathLike
        :param msa_db: path to the output MSA database to create
        :type msa_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, msa_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('result2msa', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('result2msa', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def result2dnamsa(self, query_db, target_db, result_db, msa_db, **kwargs):
        """Compute MSA DB without insertions in the query for DNA
        sequences.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to compute MSAs from
        :type result_db: os.PathLike
        :param msa_db: path to the output MSA database to create
        :type msa_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, msa_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('result2dnamsa', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('result2dnamsa', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def result2stats(self, query_db, target_db, result_db, stats_db, **kwargs):
        """Compute statistics for each entry in a DB.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to compute statistics for
        :type result_db: os.PathLike
        :param stats_db: path to the output statistics database to create
        :type stats_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, stats_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('result2stats', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('result2stats', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def filterresult(self, query_db, target_db, result_db, filter_db, **kwargs):
        """Pairwise alignment result filter.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param result_db: path to the result database to filter
        :type result_db: os.PathLike
        :param filter_db: path to the output filtered result database to create
        :type filter_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, target_db, result_db, filter_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('filterresult', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('filterresult', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def offsetalignment(self, query_db, query_orf_db, target_db,
                        target_orf_db, align_db, offset_db, **kwargs):
        """Offset alignment by ORF start position.

        :param query_db: path to the query sequence database that was used
        :type query_db: os.PathLike
        :param query_orf_db: path to the query ORF database that was used
        :type query_orf_db: os.PathLike
        :param target_db: path to the target sequence database that was used
        :type target_db: os.PathLike
        :param target_orf_db: path to the target ORF database that was used
        :type target_orf_db: os.PathLike
        :param align_db: path to the alignment database to offset
        :type align_db: os.PathLike
        :param offset_db: path to the output offset alignment database to create
        :type offset_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [query_db, query_orf_db, target_db, target_orf_db,
                align_db, offset_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('offsetalignment', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('offsetalignment', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def proteinaln2nucl(self, nucl_query_db, nucl_target_db, aa_query_db,
                        aa_target_db, aa_align_db, nucl_align_db, **kwargs):
        """Transform protein alignments to nucleotide alignments.

        :param nucl_query_db: path to the nucleotide query sequence database
        :type nucl_query_db: os.PathLike
        :param nucl_target_db: path to the nucleotide target sequence database
        :type nucl_target_db: os.PathLike
        :param aa_query_db: path to the protein query sequence database
        :type aa_query_db: os.PathLike
        :param aa_target_db: path to the protein target sequence database
        :type aa_target_db: os.PathLike
        :param aa_align_db: path to the protein alignment database
        :type aa_align_db: os.PathLike
        :param nucl_align_db: path to the output nucleotide alignment database
        :type nucl_align_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [nucl_query_db, nucl_target_db, aa_query_db, aa_target_db,
                aa_align_db, nucl_align_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('proteinaln2nucl', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('proteinaln2nucl', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def result2repseq(self, sequence_db, result_db, repseq_db, **kwargs):
        """Get representative sequences from result DB.

        :param sequence_db: path to the sequence database that was used
        :type sequence_db: os.PathLike
        :param result_db: path to the result database to extract representative
            sequences from
        :type result_db: os.PathLike
        :param repseq_db: path to the output representative sequence database
            to create
        :type repseq_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, result_db, repseq_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('result2repseq', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('result2repseq', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def sortresult(self, result_db, sorted_db, **kwargs):
        """Sort a result DB in the same order as the prefilter or
        align module.

        :param result_db: path to the result database to sort
        :type result_db: os.PathLike
        :param sorted_db: path to the output sorted result database to create
        :type sorted_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [result_db, sorted_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('sortresult', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('sortresult', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def summarizealis(self, alignment_db, summary_db, **kwargs):
        """Summarize alignment result to one row (uniq. cov., cov.,
        avg. seq. id.).

        :param alignment_db: path to the alignment database to summarize
        :type alignment_db: os.PathLike
        :param summary_db: path to the output summary database to create
        :type summary_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [alignment_db, summary_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('summarizealis', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('summarizealis', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def summarizeresult(self, alignment_db, summary_db, **kwargs):
        """Extract annotations from alignment DB.

        :param alignment_db: path to the alignment database to summarize
        :type alignment_db: os.PathLike
        :param summary_db: path to the output summary database to create
        :type summary_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [alignment_db, summary_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('summarizeresult', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('summarizeresult', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # Taxonomy assignment
    def createtaxdb(self, sequence_db, tmpdir, **kwargs):
        """Add taxonomic labels to sequence DB.

        :param sequence_db: path to the sequence database to add taxonomic
            labels to
        :type sequence_db: os.PathLike
        :param tmpdir: path to the temporary directory to use
        :type tmpdir: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db, tmpdir]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createtaxdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createtaxdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createbintaxonomy(self, names_dmp, nodes_dmp, merged_dmp, tax_file,
                          **kwargs):
        """Create binary taxonomy from NCBI input.

        :param names_dmp: path to the names.dmp file
        :type names_dmp: os.PathLike
        :param nodes_dmp: path to the nodes.dmp file
        :type nodes_dmp: os.PathLike
        :param merged_dmp: path to the merged.dmp file
        :type merged_dmp: os.PathLike
        :param tax_file: path to the output taxonomy file to create
        :type tax_file: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [names_dmp, nodes_dmp, merged_dmp, tax_file]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('createbintaxonomy', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('createbintaxonomy', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def createbintaxmapping(self):
        """Create binary taxonomy mapping from tabular taxonomy mapping."""
        if self.api_version < 14:
            raise MMseqs2Error("'createbintaxmapping' was not introduced until "
                               "MMseqs2 version 14")

        # TODO: implement createbintaxmapping after verifying signature
        raise NotImplementedError("createbintaxmapping is not yet implemented")

    def addtaxonomy(self, target_db, result_db, tax_result_db, **kwargs):
        """Retroactively add taxonomic labels to a result DB.

        :param target_db: path to the target sequence database that was used
            in search - must have taxonomy annotations for this to work
        :type target_db: os.PathLike
        :param result_db: path to the result database to add taxonomic labels to
        :type result_db: os.PathLike
        :param tax_result_db: path to the output taxonomy-added result
            database to create
        :type tax_result_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [target_db, result_db, tax_result_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('addtaxonomy', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('addtaxonomy', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def filtertaxdb(self, target_db, tax_db, filtered_db, **kwargs):
        """Filter a taxonomy-aware result database."""
        raise NotImplementedError("filtertaxdb is not yet implemented")

    def filtertaxseqdb(self, target_db, filtered_db, **kwargs):
        """Filter a taxonomy-aware sequence database.

        :param target_db: path to a taxonomy-annotated sequence database
            to be filtered
        :type target_db: os.PathLike
        :param filtered_db: path to the output taxonomy-filtered sequence
            database
        :type filtered_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if any positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [target_db, filtered_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('filtertaxseqdb', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('filtertaxseqdb', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    def aggregatetax(self, tax_seq_db, set_to_seq_map, tax_res_per_seq_db,
                     tax_res_per_set_db):
        """Aggregate multiple taxon labels to a single label."""
        raise NotImplementedError("aggregatetax is not yet implemented")

    def aggregatetaxweights(self, tax_seq_db, set_to_seq_map,
                            tax_res_per_seq_db, tax_aln_res_per_seq_db,
                            tax_res_per_set_db):
        """Aggregate multiple taxon labels to a single label."""
        raise NotImplementedError("aggregatetaxweights is not yet implemented")

    def lcaalign(self, query_db, target_db, result_db, alignment_db, **kwargs):
        """Efficient gapped alignment for lca computation.

        :param query_db:
        :type query_db:
        :param target_db:
        :type target_db:
        :param result_db:
        :type result_db:
        :param alignment_db:
        :type alignment_db:
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        """
        raise NotImplementedError("lcaalign is not yet implemented")

    def lca(self):
        """Compute the lowest common ancestor."""
        raise NotImplementedError("lca is not yet implemented")

    def majoritylca(self):
        """Compute the lowest common ancestor using majority voting."""
        raise NotImplementedError("majoritylca is not yet implemented")

    # Multi-hit search
    def multihitdb(self):
        """Create sequence DB for multi hit searches."""
        raise NotImplementedError("multihitdb is not yet implemented")

    def multihitsearch(self):
        """Search with a grouped set of sequences against another
        grouped set."""
        raise NotImplementedError("multihitsearch is not yet implemented")

    def besthitperset(self):
        """For each set of sequences compute the best element and
        update p-value."""
        raise NotImplementedError("besthitperset is not yet implemented")

    def combinepvalperset(self):
        """For each set compute the combined p-value."""
        raise NotImplementedError("combinepvalperset is not yet implemented")

    def mergeresultsbyset(self):
        """Merge results from multiple ORFs back to their respective
         contig."""
        raise NotImplementedError("mergeresultsbyset is not yet implemented")

    # Prefiltering
    def prefilter(self):
        """Double consecutive diagonal k-mer search."""
        raise NotImplementedError("prefilter is not yet implemented")

    def ungappedprefilter(self):
        """Optimal diagonal score search."""
        raise NotImplementedError("ungappedprefilter is not yet implemented")

    def gappedprefilter(self):
        """Optimal Smith-Waterman-based prefiltering (slow)."""
        if self.api_version < 15:
            raise MMseqs2Error("'gappedprefilter' was not introduced until "
                               "MMseqs2 version 15")

        raise NotImplementedError("gappedprefilter is not yet implemented")

    def kmermatcher(self):
        """Find bottom-m-hashed k-mer matches within sequence DB."""
        raise NotImplementedError("kmermatcher is not yet implemented")

    def kmersearch(self):
        """Find bottom-m-hashed k-mer matches between target and query
        DB."""
        raise NotImplementedError("kmersearch is not yet implemented")

    # Alignment
    def align(self):
        """Optimal gapped local alignment."""
        raise NotImplementedError("align is not yet implemented")

    def alignall(self, sequence_db, result_db, alignment_db, **kwargs):
        """Within-result all-vs-all gapped local alignment."""
        raise NotImplementedError("alignall is not yet implemented")

    def transitivealign(self):
        """Transfer alignments via transitivity."""
        raise NotImplementedError("transitivealign is not yet implemented")

    def rescorediagonal(self):
        """Compute sequence identity for diagonal."""
        raise NotImplementedError("rescorediagonal is not yet implemented")

    def alignbykmer(self):
        """Heuristic gapped local k-mer based alignment."""
        raise NotImplementedError("alignbykmer is not yet implemented")

    # Clustering
    def clust(self):
        """Cluster result by Set-Cover/Connected-Component/Greedy-
        Incremental."""
        raise NotImplementedError("clust is not yet implemented")

    def clusthash(self):
        """Hash-based clustering of equal length sequences."""
        raise NotImplementedError("clusthash is not yet implemented")

    def mergeclusters(self):
        """Merge multiple cascaded clustering steps."""
        raise NotImplementedError("mergeclusters is not yet implemented")

    # Profile databases
    def result2profile(self):
        """Compute profile DB from a result DB."""
        raise NotImplementedError("result2profile is not yet implemented")

    def msa2result(self):
        """Convert an MSA DB to a profile DB."""
        raise NotImplementedError("msa2result is not yet implemented")

    def sequence2profile(self):
        """Turn a sequence into profile by adding context specific pseudo
        counts."""
        if self.api_version < 14:
            raise MMseqs2Error("'sequence2profile' was not introduced until "
                               "MMseqs2 version 14")

        raise NotImplementedError("sequence2profile is not yet implemented")

    def profile2pssm(self):
        """Convert a profile DB to a tab-separated PSSM file."""
        raise NotImplementedError("profile2pssm is not yet implemented")

    def profile2neff(self):
        """Convert a profile DB to a tab-separated list of Neff scores."""
        if self.api_version < 15:
            raise MMseqs2Error("'profile2neff' was not introduced until "
                               "MMseqs2 version 15")

        raise NotImplementedError("profile2neff is not yet implemented")

    def profile2consensus(self):
        """Extract consensus sequence DB from a profile DB."""
        raise NotImplementedError("profile2consensus is not yet implemented")

    def profile2repseq(self):
        """Extract representative sequence DB from a profile DB."""
        raise NotImplementedError("profile2repseq is not yet implemented")

    def convertprofiledb(self):
        """Convert a HH-suite HHM DB to a profile DB."""
        raise NotImplementedError("convertprofiledb is not yet implemented")

    # Profile-profile databases
    def enrich(self):
        """Boost diversity of search result."""
        if self.api_version > 13:
            raise MMseqs2Error("'enrich' pipeline was removed from MMseqs2 "
                               "after version 13")

        raise NotImplementedError("enrich is not yet implemented")

    def result2pp(self):
        """Merge two profile DBs by shared hits."""
        if self.api_version > 13:
            raise MMseqs2Error("'result2pp' was removed from MMseqs2 after "
                               "version 13")

        raise NotImplementedError("result2pp is not yet implemented")

    def profile2cs(self):
        """Convert a profile DB into a column state sequence DB."""
        if self.api_version > 13:
            raise MMseqs2Error("'profile2cs' was removed from MMseqs2 after "
                               "version 13")

        raise NotImplementedError("profile2cs is not yet implemented")

    def tsv2exprofiledb(self):
        """Create an expandable profile db from TSV files."""
        if self.api_version < 14:
            raise MMseqs2Error("'tsv2exprofiledb' was not introduced until "
                               "MMseqs2 version 14")

        raise NotImplementedError("tsv2exprofiledb is not yet implemented")

    def convertca3m(self):
        """Convert a cA3M DB to a result DB."""
        raise NotImplementedError("convertca3m is not yet implemented")

    def expandaln(self):
        """Expand an alignment result based on another."""
        raise NotImplementedError("expandaln is not yet implemented")

    def expand2profile(self):
        """Expand an alignment result based on another and create a
        profile."""
        raise NotImplementedError("expand2profile is not yet implemented")

    # Utility modules to manipulate DBs
    def setextendeddbtype(self):
        """Write an extended DB."""
        if self.api_version < 15:
            raise MMseqs2Error("'setextendeddbtype' was not introduced until "
                               "MMseqs2 version 15")

        raise NotImplementedError("setextendeddbtype is not yet implemented")

    def view(self):
        """Print DB entries given in --id-list to stdout."""
        raise NotImplementedError("view is not yet implemented")

    def apply(self):
        """Execute given program on each DB entry."""
        raise NotImplementedError("apply is not yet implemented")

    def filterdb(self):
        """DB filtering by given conditions."""
        raise NotImplementedError("filterdb is not yet implemented")

    def swapdb(self):
        """Transpose DB with integer values in first column."""
        raise NotImplementedError("swapdb is not yet implemented")

    def prefixid(self):
        """For each entry in a DB prepend the entry key to the entry
        itself."""
        raise NotImplementedError("prefixid is not yet implemented")

    def suffixid(self):
        """For each entry in a DB append the entry key to the entry
        itself."""
        raise NotImplementedError("suffixid is not yet implemented")

    def renamedbkeys(self):
        """Create a new DB with original keys renamed."""
        raise NotImplementedError("renamedbkeys is not yet implemented")

    # Special-purpose utilities
    def diffseqdbs(self):
        """Compute diff of two sequence DBs."""
        raise NotImplementedError("diffseqdbs is not yet implemented")

    def summarizetabs(self):
        """Extract annotations from HHblits BLAST-tab-formatted
        results."""
        raise NotImplementedError("summarizetabs is not yet implemented")

    def gff2db(self):
        """Extract regions from a sequence database based on a GFF3
        file."""
        raise NotImplementedError("gff2db is not yet implemented")

    def maskbygff(self):
        """Mask out sequence regions in a sequence DB by features
        selected from a GFF3 file."""
        raise NotImplementedError("maskbygff is not yet implemented")

    def convertkb(self):
        """Convert UniProtKB data to a DB."""
        raise NotImplementedError("convertkb is not yet implemented")

    def summarizeheaders(self):
        """Summarize FASTA headers of result DB."""
        raise NotImplementedError("summarizeheaders is not yet implemented")

    def nrtotaxmapping(self):
        """Create taxonomy mapping for NR database."""
        raise NotImplementedError("nrtotaxmapping is not yet implemented")

    def extractdomains(self, alignment_db, msa_db, domain_db, **kwargs):
        """Extract the highest scoring alignment regions for each sequence.

        :param alignment_db: path to the alignment database to extract domains
            from
        :type alignment_db: os.PathLike
        :param msa_db: path to the MSA database to extract domains from
        :type msa_db: os.PathLike

        """
        # TODO: Implement extractdomains after determining the correct signature
        raise NotImplementedError("extractdomains is not yet implemented")

    def countkmer(self, sequence_db, **kwargs):
        """Count k-mers.

        :param sequence_db: path to the sequence database to count k-mers in
        :type sequence_db: os.PathLike
        :param kwargs: additional arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: stdout from the MMseqs2 command
        :rtype: str
        :raises ValueError: if the positional argument cannot be cast to a
            pathlib.Path
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
            code
        """
        args = [sequence_db]

        # Cast all positional arguments to fully resolved pathlib.Path objects
        args = tuple(str(arg) for arg in _cast_list_objs_to_paths(args))

        # Validate kwargs - don't catch errors here, let them go to the caller
        self._validate_kwargs('countkmer', **kwargs)

        # Run the MMseqs2 command; if non-zero exit code occurs, MMseqs2Error
        # will be raised by _run(); otherwise, return the stdout
        process = self._run('countkmer', *args, **kwargs)

        # If the process completes successfully, return the stdout
        return process.stdout.rstrip()

    # General purpose utilities
    def help(self, pipeline=None):
        """Return help message for a specific pipeline or for all pipelines.

        :param pipeline: name of the pipeline to get help for
        :type pipeline: str
        """
        if not pipeline:
            return self._run('--help').stdout.rstrip()

        if "_" in pipeline:
            pipeline = pipeline.replace("_", "-")

        return self._run(pipeline, '--help').stdout.rstrip()

    def version(self):
        """Return the version of MMseqs2.

        :return: version string from MMseqs2
        :rtype: str
        """
        return self._run('version').stdout.rstrip()

    # The following methods should never be called directly by the user,
    # and are the only methods that don't return strings.
    def _validate_kwargs(self, pipeline, **kwargs):
        """Check that the given keyword arguments are valid for the
        specified pipeline, and that their values are of the correct
        type.

        No validation is done to ensure that the values themselves are
        valid for the specified pipeline (e.g., an invalid header column
        name specified with `--format-output` of the `search` pipeline).

        :param pipeline: name of the pipeline to validate arguments for
        :type pipeline: str
        :param kwargs: keyword arguments to validate
        :type kwargs: dict[str, str | int | float | bool | None]
        :return: True if all kwargs are valid and of the correct type
        :rtype: bool
        :raises KeyError: if any kwarg is not valid for the specified
            pipeline
        :raises TypeError: if any kwarg value is not of the correct type
        """
        # Get the help text for the specified pipeline, and extract the valid
        # arguments from it
        pipeline_help = self.help(pipeline)
        valid_args = set()
        for line in pipeline_help.split('\n'):
            if line.startswith(' -'):
                valid_args.add(line.split()[0])

        # Check that all passed kwargs are valid and of the correct type
        for kw_arg, kw_value in kwargs.items():
            # Assume arguments starting with "-" or "--" are pre-formatted
            if kw_arg.startswith("-") or kw_arg.startswith("--"):
                real_kw_arg = kw_arg
            # Format any other arguments
            # E.g., "a" -> "-a"; "--do_thing" -> "--do-thing"
            elif len(kw_arg) == 1:
                real_kw_arg = f'-{kw_arg.replace("_", "-")}'
            else:
                real_kw_arg = f'--{kw_arg.replace("_", "-")}'

            if real_kw_arg not in valid_args:
                # When raising error, report the kw_arg as it was passed
                raise KeyError(f"{pipeline} does not accept argument {kw_arg}")

            # Check that the value is of the correct type
            if not isinstance(kw_value, eval(self.arg_types[real_kw_arg])):
                raise TypeError(f"value passed with '{kw_arg}' must be of type"
                                f" {self.arg_types[real_kw_arg]}, not "
                                f"{type(kw_value)}")

        # If we got here without raising an error, we think all kwargs are valid
        return True

    def _run(self, pipeline, *args, **kwargs):
        """Run MMseqs2 with the given command and arguments.

        :param pipeline: name of the MMseqs2 pipeline to run
        :type pipeline: str
        :param args: positional arguments to pass to the MMseqs2 command
        :type args: str
        :param kwargs: keyword arguments to pass to the MMseqs2 command
        :type kwargs: str | int | float | bool
        :return: a CompletedProcess object from `subprocess.run()`
        :rtype: sp.CompletedProcess
        :raises FileNotFoundError: if the MMseqs2 binary is not found
        :raises MMseqs2Error: if the MMseqs2 command returns a non-zero exit
        """
        command = [self.binary_path, pipeline]
        for arg in args:
            if isinstance(arg, list):
                command.extend(arg)
            else:
                command.append(str(arg))

        for kw_arg, kw_value in kwargs.items():
            if len(kw_arg) == 1:
                command.append(f'-{kw_arg.replace("_", "-")}')
            else:
                command.append(f'--{kw_arg.replace("_", "-")}')

            # If the value is a boolean, only include the flag if it's True
            if isinstance(kw_value, bool):
                if kw_value:
                    continue
                else:
                    command.pop()

            # If the value is not a boolean, include it in the command
            command.append(str(kw_value))

        # Log the command to be run - assumes logging context has been set up
        logging.info(f"Running MMseqs2 command: {' '.join(command)}")

        try:
            # check=True raises CalledProcessError if the exit code is non-zero
            # text=True will capture stdout and stderr as strings
            return sp.run(command, stdout=sp.PIPE, stderr=sp.PIPE,
                          check=True, text=True)
        except sp.CalledProcessError as e:
            message = ". ".join(e.stdout.rstrip().split("\n")[-3:-1])
            raise MMseqs2Error(message)
