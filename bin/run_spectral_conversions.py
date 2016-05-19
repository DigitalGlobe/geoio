#!/usr/bin/env python

import argparse
import os

import geoio.dg
import tinytools as tt
import geoio


def setup():
    ### parse command line input ###
    descrip = """Create the spectral conversion files requested by the arguments."""

    parser = argparse.ArgumentParser(description=descrip)
    parser.add_argument("-input", "--input",
                        required=True,
                        nargs='+',
                        help="A list of file strings to operate on.")
    parser.add_argument("-spectral_ops", "--spectral_ops",
                        required=True,
                        choices=['radiance', 'toa', 'DGAComp'],
                        nargs='+',
                        help="Which spectral operations to perform on the "
                             "input files.")
    parser.add_argument("-outdir", "--outdir",
                        required=False,
                        nargs=1,
                        help="*** NOT IMPLEMENTED YET ***"
                             "Location to write the output files if not the "
                             "default rename at the input location.")
    parser.add_argument("-test_only", "--test_only",
                        required=False,
                        action="store_true",
                        default=False,
                        help="This flag will turn off all spectral operations "
                             "and will return the files that would have "
                             "been operated on if the flag was not passed.")
    parser.add_argument("-clear_all_spectral_files",
                        "--clear_all_spectral_files",
                        required=False,
                        action="store_true",
                        default=False,
                        help="This will delete ALL the spectral files "
                             "associated with "
                             "the files passed in.  This will short circuit"
                             "all spectral operations.  It is STRONGLY "
                             "suggested "
                             "that a test_only check is performed before "
                             "running to prevent unwanted data removal.")

    return parser.parse_args()

    assert (parser.spectral_ops or parser.clear_all_spectral_files), \
        "Either a spectral operation of a deletion of files needs to " \
        "be requested."


def main(args):
    ### Input files ###
    flist = args.input

    ### Output Run Info ###
    # Print the files that are going to be manipulated and, if test_only, exit
    print('')
    print("The input files are:")
    for x in flist:
        print(x)

    ### Validation ###
    # Verify the files all open correctly - you wouldn't want to have stuff crash
    # part of the way through a long run.  So, this will proivde a bit
    # of validation.
    print('')
    print("Verifying the files have the appropriate meta data...")
    for x in flist:
        i1 = geoio.dg.DGImage(x)
        assert i1.meta_dg.IMD
        assert i1.meta.abscalfactor
        assert i1.meta.effbandwidth
        del i1
    print("... they do!")


    # Exit here if the user requested test_only
    if args.test_only:
        return

    ### The Heavy Lifting ###
    # Run the spectral conversions
    if args.clear_all_spectral_files:
        print('')
        print("Deleting all spectral files...")
        for x in flist:
            i1 = geoio.dg.DGImage(x)
            i1.delete_all_spectral_files(test_only=False)
            del i1

    if args.spectral_ops:
        print('')
        print("Each file will be converted to {}:  ".format(args.spectral_ops))
        print('Running the spectral conversions...')
        for x in flist:
            print('')
            print("Staring on file {}:".format(x))
            i1 = geoio.dg.DGImage(x)
            if 'radiance' in args.spectral_ops:
                print("Running at sensor radiance...")
                i1.create_at_sensor_rad_files()
            if 'toa' in args.spectral_ops:
                print("Running toa reflectance...")
                i1.create_toa_ref_files()
            if 'DGAComp' in args.spectral_ops:
                print("Running DGAComp reflectance...")
                i1.create_dgacomp_ref_files()
            del i1

    print('*** Run Complete ***')


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
