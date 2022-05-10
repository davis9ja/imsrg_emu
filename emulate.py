##########################################################
# Main executable for running an emulator interactively. #
#                                                        #
# Author: Jacob Davison                                  #
# Date:   05/05/2022                                     #
##########################################################

import sys
import argparse
import pprint
import numpy as np

import dmd_rkoi as drk
import dmd_std as dst
from imsrg_emu.utils.get_log_data import get_log_data
from imsrg_emu.utils.make_argparser import make_argparser
from imsrg_emu.utils.make_plots import make_energy_plots, make_correlation_plots

args = make_argparser()

print("Args from command line:")
pprint.pprint(args)
print()

dmd = None
test_data = None

if args['emu_method'] == 'standard':

    print("Reading single flow data from ", args['dataPath'])
    #data_matrix = get_log_data(args['dataPath'])
    data_matrix = np.loadtxt(args['dataPath'], delimiter=',', comments="#").T    
    test_data = data_matrix

    print("Fitting standard DMD emulator")
    dmd = dst.DMD_STD()
    dmd.fit(data_matrix, args['nobs'], r=args['trunc'], enforce_physics=True)

elif args['emu_method'] == 'parametric':

    data_list = []
    with open(args['dataPath'], 'r') as f:
        for line in f:
            data_matrix = np.loadtxt(line.replace("\n",""), delimiter=',', comments="#").T
            data_list.append(data_matrix)

    with open(args['paramList'], 'r') as f:
        params = np.asarray(f.readlines(), dtype=np.float64)

    if args['testPath'] is not None:
        test_data = np.loadtxt(args['testPath'], delimiter=',', comments='#').T

    if args['emuType'] == 'rKOI':
        print("Fitting rKOI DMD emulator")
        dmd = drk.DMD_rKOI()
        dmd.fit(data_list, params, args['nobs'], r=args['trunc'])
        dmd.interp_dmd(args['testParam'])

    print(params)
    for data in data_list:
        print('{:0.16f} {}'.format(data[0,0], data.shape))
print("Printing results...")

s_range = np.arange(args['t0'], args['t1']+args['dt'], args['dt'])
pred = dmd.predict(s_range, args['dt'])

print("{:<10s} | {:<10s}".format("s", "E"))
print("-----------------------")
for i,s in enumerate(s_range):
    print("{:10.7f} | {:10.7f}".format(s, pred[0,i]))

if args['plot']:
    assert pred.shape == test_data.shape, "matrices shape {}, {} are not compatible; PLOT requires that predict data and test data are same shape".format(pred.shape, test_data.shape)
    
    make_energy_plots(s_range, pred[0,:], test_data[0,:], args)
    make_correlation_plots(s_range, pred, test_data, args)
