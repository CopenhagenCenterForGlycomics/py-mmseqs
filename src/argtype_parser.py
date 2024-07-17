"""Parse MMseqs2 parameters from Parameters.cpp to a JSON file. Used to
generate api_xx_argtypes.json files for new versions of MMseqs2."""

import argparse
import json
import pathlib
import sys


def parse_args():
    """Parse commandline arguments."""
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=pathlib.Path,
                   help='path to Parameters.cpp to be read')
    p.add_argument('dest', type=pathlib.Path,
                   help='path to JSON file to be written')

    return p.parse_args()


def main():
    """Commandline entrypoint."""
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    args = parse_args()

    params = {}
    with open(args.source, 'r') as param_reader:
        for line in param_reader:
            if line.startswith('        PARAM_'):
                line = line.split(", ")

                # Name of the parameter
                name = line[1]

                # Find the type declaration
                for elem in line:
                    if elem.startswith('typeid'):
                        type_decl = elem[7:-1]
                        break

                # There are several types that don't translate in Python
                # NOTE: this is not comprehensive, so the output file will need
                # to be checked and manually corrected
                if type_decl == "std::string":
                    type_decl = "str"
                elif type_decl == "double":
                    type_decl = "float"
                elif type_decl == "MultiParam<int>":
                    type_decl = "int"
                elif type_decl == "MultiParam<char*>":
                    type_decl = "str"
                elif type_decl == "MultiParam<float>":
                    type_decl = "float"
                elif type_decl == "ByteParser":
                    type_decl = "str"

                params[name] = type_decl

    with open(args.dest, 'w') as param_writer:
        json.dump(params, param_writer, indent=4)


if __name__ == "__main__":
    main()