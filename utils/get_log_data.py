###############################################################
# Get log data from the C++ IMSRG code for the pairing model. #
#                                                             #
# Author: Jacob Davison                                       #
# Date:   05/05/2022                                          #
###############################################################

import numpy as np

def get_log_data(data_path):
    """Read .log.imsrg data from TCIMSRG.
    """
    with open(data_path, 'r') as f:
        
        lines = f.readlines()

        lines_trunc = lines[7:]
        lines_trunc = lines_trunc[:-1]

        data_matrix = np.array([np.array(line.split(','))[0:-1] for line in lines_trunc], dtype=np.float64).T


    return data_matrix
