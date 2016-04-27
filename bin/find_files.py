#!/usr/bin/env python

import argparse
import os
import tinytools as tt
import geoio


def setup():
    ### parse command line input ###
    descrip = """Create the spectral conversion files requested by the arguments."""

    parser = argparse.ArgumentParser(description=descrip)
    parser.add_argument("-input", "--input",
                        required=True,
                        nargs='+',
                        help="Input can be of four forms:"
                             "1. A list of file strings to operate on (this "
                                "will simply be passed through)."
                             "2. A string or list of strings with glob type "
                                "search patterns."
                             "3. A directory that will be searched for "
                                ".TIL files"
                             "4. A file to parse for a list of file strings")

    return parser.parse_args()

def main(args):
    ### Find files to run on ###
    # Figure out what the input format is and find files
    input = args.input
    flist = []

    # print(args.input)
    if len(input) > 1:
        # if isinstance(input,list):
        # Then this is a list of file strings or a list of glob patterns
        for x in input:
            # This decided if it is a file or a glob and expands as needed
            flist.extend(return_from_single_string(x))
    # elif isinstance(input,str):
    else:
        input = input[0]
        print(input)
        # If just a string, then it could be a file, glob, dir to search, or a
        # a file to read.
        if os.path.isdir(input):
            flist.extend(tt.files.search(input, '*.TIL', depth=0))
        elif os.path.isfile(input):
            e = os.path.splitext(input)[-1]
            if tt.files.filter(e, '.txt', case_sensitive=False):
                # Then this is a file of file strings
                print('*** NOT IMPLEMENTED YET ***')
        else:
            flist.extend(return_from_single_string(input))

    return flist

### Repeated string operation ###
def return_from_single_string(x):
    if os.path.isfile(x):
        # Then this is a file so just add it to the list
        return (x)
    else:
        # Then this is probably a glob pattern
        (path, base) = os.path.split(x)
        if not path:
            path = os.getcwd()
        return tt.files.search(path, base, depth=0)


if __name__ == "__main__":
    args = setup()
    main(args)
