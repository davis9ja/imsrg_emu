##############################################################
# Make the parser that parses the command line arguments for #
# the emulator.                                              #
#                                                            #
# Author: Jacob Davison                                      #
# Date:   05/05/2022                                         #
##############################################################

import argparse

def make_argparser():
    parser = argparse.ArgumentParser(description='Emulate dynamical system using DMD. Intended for IMSRG but not exclusive.')
    subparsers = parser.add_subparsers(title='emulation type', description='valid emulation types', help='emulation type', dest='emu_method')
    
    parser_std = subparsers.add_parser('standard', help='run a standard DMD emulation on the data')
    parser_std.add_argument('dataPath', help='path/to/data/file; CSV with rows of snapshots')

    parser_par = subparsers.add_parser('parametric', help='run a parametric DMD emulation on the data')
    parser_par.add_argument('dataPath', help='path/to/data/list; text list of files, each CSV with rows of snapshots')
    parser_par.add_argument('paramList', help='path/to/param/file; CSV of param for each training sample in dataPath')
    parser_par.add_argument('testParam', help='param set to test with parametric DMD (must fall within paramList range)')
    parser_par.add_argument('--testPath', '-tp', default=None, help='/path/to/test/data; CSV with rows of snapshots of data to test (mostly for plotting)')
    parser_par.add_argument('--emuType', '-e', required=True, choices=['rKOI', 'rEPI'], help='Choice of parametric emulator type. Choices are reduced Koopman Interpolation (rKOI) and reduced Eigenpair Interpolation (rEPI).')

    for subparser in [parser_std, parser_par]:
        subparser.add_argument('-N', '--nobs', type=int, default=10, help='number of snapshots per DMD operator')
        subparser.add_argument('-E', '--exact', type=bool, default=False, help='True or False: compute DMD operator exactly instead of reduced')
        subparser.add_argument('-T', '--trunc', type=int, default=6, help='SVD truncation rank in reduced DMD')
        subparser.add_argument('-t', '--tol', type=float, default=None, help='SVD singular value tolerance (if None, default to --trunc)')
        subparser.add_argument('--t0', type=float, default=0.0, help='starting point for emulation')
        subparser.add_argument('--t1', type=float, default=10.0, help='ending point for emulation')
        subparser.add_argument('--dt', type=float, default=0.1, help='emulation step width')
        subparser.add_argument('--plot', type=bool, default=False, help='make plots (default directory: ./plots/)')

    return vars(parser.parse_args())
