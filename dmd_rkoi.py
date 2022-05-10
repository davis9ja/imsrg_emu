#######################################################################
# Implemenation of the reduced Koopman Operator Interpolation method  #
# for parametric DMD.                                                 #
#                                                                     #
# Author: Jacob Davison                                               #
# Date:   05/05/2022                                                  #
#######################################################################

import numpy as np
import scipy.linalg as la
import scipy.interpolate
from imsrg_emu.utils.get_log_data import get_log_data

class DMD_rKOI(object):
    """
    Reduced Koopman Operator Interpolation for parametric DMD. Huhn et al. 2022 (arXiv:2204.12006v1)
    """
    
    def __init__(self):

        self._Ar_training = None
        self._Ur_training = None
        self._b_training = None

        self._AI = None
        self._UI = None
        self._bI = None

        self._Phi_p = None
        self._eigs_p = None
        self._b_p = None

        self._Ar_shape = None
        self._Ur_shape = None

    @property
    def Ar_training(self):
        return self._Ar_training

    @property
    def Ur_training(self):
        return self._Ur_training

    @property
    def b_training(self):
        return self._b_training

    @property
    def AI(self):
        return self._AI

    @property
    def UI(self):
        return self._UI

    @property
    def bI(self):
        return self._bI

    @property
    def Phi_p(self):
        return self._Phi_p

    @property
    def eigs_p(self):
        return self._eigs_p

    @property
    def b_p(self):
        return self._b_p

    def fit(self, data_list, parameters, nobs_t, r=6):
        """Fit the interpolators to build the parametric DMD system.
        
        Arguments:
        
        data_list -- list of numpy snapshot matrices (list of standard DMD inputs)
        parameters -- list of parameters corresponding to each input in data_list
        nobs_t -- number of observations to use per DMD input
        r -- truncation rank of SVD in each DMD input
        """

        Ar_training = []
        Ur_training = []
        b_training = []

        for data in data_list:
            X = data[:, :nobs_t-1]
            Xp = data[:, 1:nobs_t]

            U,sigma,Vt = la.svd(X, full_matrices=False)

            if isinstance(r, int):
                Ur = U[:, :r]
                sr = sigma[:r]
                Vtr = Vt[:r, :]
            elif isinstance(r, float):
                keep_idx = np.argwhere(sigma > r)[:,0]
                print(keep_idx.shape, keep_idx)

                Ur = U[:, keep_idx]
                sr = sigma[keep_idx]
                Vtr = Vt[keep_idx, :]
            print(Ur.shape, sr.shape, Vtr.shape)
            Ar = Ur.conj().T@Xp@Vtr.conj().T@np.diag(np.reciprocal(sr))#la.inv(np.diag(sr))
            Ar_training.append(np.reshape(Ar,(-1,)))
            Ur_training.append(np.reshape(Ur,(-1,)))

            # eigendecomp
            w,v =la.eig(Ar)

            # compute phi (or, try the other way in the Data Driven Science book)
            # which is to evaluate in column space of Xp
            phi = Ur@v
            # phi = Xp@Vtr.conj().T@la.inv(np.diag(sr))@v*np.reciprocal(w)

            # compute mode amplitudes
            b = la.lstsq(phi, X[:,0], lapack_driver='gelsd')[0]

            b_training.append(b)


            self._Ar_shape = Ar.shape
            self._Ur_shape = Ur.shape

        self._Ar_training = np.asarray(Ar_training).T
        self._Ur_training = np.asarray(Ur_training).T
        self._b_training = np.asarray(b_training).T

        self._AI = scipy.interpolate.interp1d(parameters, self.Ar_training, axis=-1)
        self._UI = scipy.interpolate.interp1d(parameters, self.Ur_training, axis=-1)
        self._bI = scipy.interpolate.interp1d(parameters, self.b_training, axis=-1)


    def interp_dmd(self, param_pred):
        """Predict DMD operator for specified parameter via interpolators in fit().

        Arguments:

        param_pred -- parameter to predict (must be within training range)
        """

        Ar_pred = np.reshape(self.AI([param_pred]), self._Ar_shape)

        w_pred, v_pred = la.eig(Ar_pred)
        
        Ur_pred = np.reshape(self.UI([param_pred]), self._Ur_shape)

        b_pred = self.bI([param_pred])[:,0]

        phi_pred = Ur_pred@v_pred
        # w_pred[0] = 1

        # -- enforce physical constraints

        # real eigs
        w_pred = np.real(w_pred)

        # positive eigs     
        idx = np.argwhere(w_pred > 0)[:,0]
        w_pred = w_pred[idx]
        phi_pred = phi_pred[:,idx]
        b_pred = b_pred[idx]

        # set "background"
        sorted_idx = np.argsort(w_pred)[::-1]
        w_pred = w_pred[sorted_idx]
        phi_pred = phi_pred[:,sorted_idx]
        b_pred = b_pred[sorted_idx]

        if w_pred[0] > 1:
            w_pred[0] = 1

        self._Phi_p = phi_pred
        self._eigs_p = w_pred
        self._b_p = b_pred

    def predict(self, s_range, ds):
        """Emulate the dynamical system over the specified range, for the specified parametric realization (in interp_dmd).

        Arguments:
        
        s_range -- list of dynamical variables to evaluate in DMD expansion
        ds -- step width

        Returns:
        
        reconstructed_data -- matrix of reconstructed snaphots from the DMD operator interpolated from fit()
        """

        assert self.Ar_training is not None, "Build DMD interpolators from fit()"
        assert self.Phi_p is not None, "Interpolate the DMD operator for param_test via interp_dmd()"

        pred_list = []
        for s in s_range:
            Xs = self.Phi_p@np.diag(np.exp(np.log(self.eigs_p)/ds*s))@self.b_p
            pred_list.append(Xs)

        reconstructed_data = np.real(np.array(pred_list).T)

        return reconstructed_data
        

# if __name__ == "__main__":
#     data_matrix1 = get_log_data('/mnt/home/daviso53/Research/tcimsrg/build/experiment/HS08-1.00--1.00-0.00-0.00-20.00-0.05.log.imsrg')
#     data_matrix2 = get_log_data('/mnt/home/daviso53/Research/tcimsrg/build/experiment/HS08-1.00-1.00-0.00-0.00-20.00-0.05.log.imsrg')

#     data_list = [data_matrix1, data_matrix2]
#     params = np.array([-1, 1])
    
#     dmd_rkoi = DMD_rKOI()
#     dmd_rkoi.fit(data_list, params, 100, 6)
#     dmd_rkoi.interp_dmd(0.5)

#     print(dmd_rkoi.Ar_training.shape)
#     print(dmd_rkoi.Ur_training.shape)
#     print(dmd_rkoi.eigs_p.shape)
#     print(dmd_rkoi.Phi_p.shape)
#     print(dmd_rkoi.b_p.shape)

#     s_range = np.arange(0, 20.05, 0.05)
#     recon_data = dmd_rkoi.predict(s_range, 0.05)
    
#     print(recon_data[0,-1])
