##############################################
# Implementation of the standard DMD method. #
#                                            #
# Author: Jacob Davison                      #
# Date:   05/05/2022                         #
##############################################

import scipy.linalg as la
import numpy as np

from sklearn.utils.extmath import randomized_svd
from sklearn.decomposition import TruncatedSVD

from imsrg_emu.utils.get_log_data import get_log_data

class DMD_STD(object):
    """Standard implementation of the reduced DMD method. Brunton et al. 2021 (arXiv:2102.12086v2)
    """
    def __init__(self):
        """Class initializer.
        """
        
        self._phi = None
        self._eigs = None
        self._b = None

    @property
    def phi(self):
        return self._phi

    @property
    def eigs(self):
        return self._eigs
    
    @property
    def b(self):
        return self._b


    def fit(self, data, nobs, exact=False, r=3, randomize=False, n_components=10, enforce_physics=False):
        """Build the DMD operator.
    
        Arguments:
        
        data -- matrix of snapshot columns, of the evolving dynamical system
        nobs -- number of observations to build the DMD operator
        
        Keyword arguments:
        
        exact -- construct the exact DMD operator, given by X*pinv(X) (default: False)
        r -- truncation rank of SVD on measurement space
        randomize -- compute the randomized SVD on the measurement space (default: False)
        n_components -- number of components to sample from the randomized SVD
        enforce_physics -- enforce physical constraints on the DMD eigenvalues (default: False)        
        """
    
        X,Xp = data[:,:nobs-1], data[:,1:nobs]

        H0 = X[:,0]

        if not exact:
            
            # Compute economy SVD with truncation
            if not randomize:
                U,s,Vh = la.svd(X, full_matrices=False)

                s = s[0:r]
                U = U[:,0:r]
                Vh = Vh[0:r,:]
                sigma = np.diag(s)

            # Compute randomized SVD sampled n_components
            else:
                U,s,Vh = randomized_svd(X, n_components, n_oversamples=2*X.shape[1]-n_components)
                sigma = np.diag(s)
            
            # Compute DMD operator
            A = U.conj().T@Xp@Vh.conj().T@la.inv(sigma)

            # Compute eigendecomposition
            w,v = la.eig(A)

            # Compute DMD modes
            phi = Xp@Vh.conj().T@la.inv(sigma)@v*np.reciprocal(w)
            
        else:

            # Compute the exact DMD
            A = Xp@la.pinv(X)

            # Compute eigendecomposition
            w,phi = la.eig(A)
            
        # Compute DMD amplitudes
        b = la.lstsq(phi, H0, lapack_driver='gelsd')[0]

        if enforce_physics:
            idx_max = np.argmax(w)
            w[idx_max] = 1
            
            # drop imaginary comps
            w = np.real(w) 
        
        # Set class attributes
        self._phi = phi
        self._eigs = w
        self._b = b
        
    def predict(self, s_range, ds):
        """Emulate the dynamical system over the specified range.

        Arguments:

        s -- dynamical variable range (list of numbers to evaluate DMD expansion)
        ds -- stepwidth for continuous time

        Returns:
        
        reconstructed_data -- matrix of reconstructed snaphots from the DMD operator built by fit()
        """

        assert self.phi is not None, "Build DMD operator first via fit()"

        pred_list = []
        for s in s_range:
            Xs = self.phi@np.diag(np.exp(np.log(self.eigs)/ds*s))@self.b
            pred_list.append(Xs)

        reconstructed_data = np.real(np.array(pred_list).T)

        return reconstructed_data

# if __name__ == "__main__":
    
#     data_matrix = get_log_data('/mnt/home/daviso53/Research/tcimsrg/build/flow/HS08-1.00-0.50-0.10-0.00-20.00-0.05.log.imsrg')
#     print(data_matrix.shape)

#     dmd = DMD_STD()
#     dmd.fit(data_matrix, 20, r=6, enforce_physics=True)
    
#     print(dmd.phi.shape, dmd.eigs, dmd.b)

#     s_range = np.arange(0, 20.05, 0.05)

#     pred = dmd.predict(s_range, 0.05)
    
#     print(pred[0,-1])

