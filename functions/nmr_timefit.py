"""NMR time fit class."""
import sys
from typing import Optional

sys.path.append("..")
  
#matplotlib.use("TkAgg")
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import least_squares

from functions.nmr_mix import NMR_Mix
  

class NMR_TimeFit(NMR_Mix):
    """Class to fit time domain FIDs to a series of exponentially decaying components.

    Attributes:
        ydata (np.ndarray): Time domain data.
        tdata (np.ndarray): Time points in seconds.
        area (np.ndarray): Area of each component.
        freq (np.ndarray): Resonance frequency of each component.
        phase (np.ndarray): Phase of each component.
        fwhmL (np.ndarray): Lorentzian FWHM of each component.
        fwhmG (np.ndarray): Gaussian FWHM of each component.
        method (str): Method for fitting, either "voigt" or "lorentzian".
        line_broadening (float): Line broadening in Hz.
        zeropad_size (Optional[int]): Zero padding size.
        sample_time (float): Dwell time in seconds.
        spectral_signal (np.ndarray): Spectral signal.
        f (np.ndarray): Frequency points in Hz.
    """

    def __init__(
        self,
        ydata: np.ndarray,
        tdata: np.ndarray,
        area: np.ndarray,
        freq: np.ndarray,
        phase: np.ndarray,
        fwhmL: np.ndarray = np.array([]),
        fwhmG: np.ndarray = np.array([]),
        method: str = "voigt",
        line_broadening: float = 0,
        zeropad_size: Optional[int] = None,
    ):
        """Initialize NMR_TimeFit class.

        Args:
            ydata (np.ndarray): Time domain data.
            tdata (np.ndarray): Time points in seconds.
            area (np.ndarray): Area guess of each component.
            freq (np.ndarray): Resonance frequency guess of each component.
            phase (np.ndarray): Phase guess of each component.
            fwhm (np.ndarray): FWHM guess of each component.
            fwhmL (np.ndarray): Lorentzian FWHM guess of each component.
            fwhmG (np.ndarray): Gaussian FWHM guess of each component.
            method (str): Fitting method, either "voigt" or "lorentzian".
            line_broadening (float): Line broadening in Hz.
            zeropad_size (Optional[int]): Zero padding size.
        """
        super().__init__(
            area=area,
            freq=freq,
            fwhmL=fwhmL,
            fwhmG=fwhmG,
            method=method,
            phase=phase,
        )
        self.line_broadening = line_broadening
        self.tdata = tdata
        if not zeropad_size:
            self.zeropad_size = np.size(tdata)
        else:
            self.zeropad_size = zeropad_size
        # apply line broadening on the time domain signal
        self.ydata = np.multiply(
            ydata, np.exp(-np.pi * self.line_broadening * self.tdata)
        )
        # calculate dwell time from delta tdata
        self.sample_time = self.tdata[1] - self.tdata[0]
        self.spectral_signal = self.sample_time * np.fft.fftshift(
            np.fft.fft(self.ydata, self.zeropad_size)
        )
        self.f = np.linspace(-0.5, 0.5, self.zeropad_size + 1) / self.sample_time
        # take out last sample to have the right number of samples
        self.f = self.f[:-1]
        self.sort_freq()

    def calc_time_fit_residual(self, bounds):
        """Fit the time domain signal using least square curve fitting.

        Running trust region reflection algorithm.
        Args:
            bounds (list): Bounds for the fitting parameters.
        """
        fun = self.get_residual_time_function
        if self.method == "voigt":
            x0 = np.array(
                [self.area, self.freq, self.fwhmL, self.fwhmG, self.phase]
            ).flatten()
        else:
            x0 = np.array([self.area, self.freq, self.fwhm, self.phase]).flatten()
        # curve fitting using trust region reflection algorithm
        ''' 
        fit_result = least_squares(
            fun=fun,
            x0=x0,
            method="trf", # trf lm
            ftol=1e-9,
            xtol=1e-10,
            bounds=bounds,
            max_nfev=30000
        )
        '''      
        fit_result = least_squares(
            fun=fun,
            x0=x0,
            bounds=bounds,
            method='trf', 
            max_nfev=30000,
            ftol=1E-9,
            xtol=1E-15,
            diff_step=None,  # None for central difference approximation
            verbose=0  # 0 for 'off', 1 for 'iter', 2 for 'final'           
        )

        # resolving the fitting results
        fit_param = fit_result["x"]
        n_fre = int(np.size(fit_param) / self.ncomp)
        fit_param = np.reshape(fit_param, (self.ncomp, n_fre))
        fit_freq = fit_param[1, :]
        # check for aliased frequency
        halfBW = 0.5 * (np.amax(self.f) - np.amin(self.f))
        alias_index = np.where(abs(fit_freq) > halfBW)
        n_alias = np.size(alias_index)
        while n_alias > 0:
            for k in range(0, n_alias):
                idx = alias_index[k]
                while fit_freq[idx] < -halfBW:
                    fit_freq[idx] = fit_freq[idx] + 2 * halfBW
                while fit_freq[idx] > halfBW:
                    fit_freq[idx] = fit_freq[idx] - 2 * halfBW

            # fit again after the alias frequency is solved
            x0 = fit_param[:]
            x0[1, :] = fit_freq  # replace the frequency to be within the range
            x0 = x0.flatten()
            fit_result = least_squares(
                fun=fun,
                x0=x0,
                jac="3-point",
                bounds=(-np.inf, np.inf),
                method="trf",
            )

            fit_param = fit_result["x"].reshape([5, 3])
            fit_freq = fit_param[1, :]

            # check for aliased frequency
            alias_index = np.where(abs(fit_freq) > halfBW)
            n_alias = np.size(alias_index)
        ''' 
        # parsing the fitting results
        fit_vec = np.multiply(
            fit_param[0, :], np.exp(1j * np.pi * fit_param[-1, :] / 180.0)
        )
        fit_area = abs(fit_vec)
        fit_freq = fit_param[1, :]
        fit_phase = np.arctan2(np.imag(fit_vec), np.real(fit_vec)) * 180.0 / np.pi
        fit_param[0, :] = fit_area
        fit_param[1, :] = fit_freq
        fit_param[4, :] = fit_phase
        
        if self.method == "voigt":
            fit_fwhmL = fit_param[2, :]
            fit_fwhmG = fit_param[3, :]

            self.set_components(
                area=fit_area,
                freq=fit_freq,
                fwhmL=fit_fwhmL,
                fwhmG=fit_fwhmG,
                phase=fit_phase,
            )
                   
        elif self.method == "lorentzian":
            fit_fwhm = fit_param[2, :]
            self.set_components(
                area=fit_area, freq=fit_freq, fwhm=fit_fwhm, phase=fit_phase
            )
        else:
            raise ValueError("Unknown fitting method.")
        '''
        return fit_param


    def get_residual_time_function(self, x: np.ndarray):
        """Calculate the residual of fitting.

        Args:
            x (np.ndarray): Fitting parameters of shape [area, freq, fwhmL, fwhmG,
             phase]
        """
        if self.method == "voigt":
            x = np.reshape(x, (5, int(np.size(x) / 5)))
            tmpNMRMix = NMR_Mix(
                area=x[0, :],
                freq=x[1, :],
                fwhmL=x[2, :],
                fwhmG=x[3, :],
                phase=x[4, :],
                method="voigt",
            )
        else:
            raise ValueError("Only voigt method is supported for time domain fitting.")

        complex_fit_time = tmpNMRMix.get_time_function(self.tdata)
        fit_sig = np.array([np.real(complex_fit_time), np.imag(complex_fit_time)])
        truth_sig = np.array([np.real(self.ydata), np.imag(self.ydata)])
        residual = (truth_sig - fit_sig).flatten()
        return residual

    def fit_time_signal_residual(
        self,
        bounds: tuple = (-np.inf, np.inf),
    ):
        """Fit the time domain signal using least square curve fitting.

        Calls scipy.optimize.least_squares to fit the time domain signal.
        Also store the fitting results in the class.

        Args:
            bounds (tuple, optional): Bounds for the fitting parameters.
                Defaults to (-np.inf, np.inf).
        """
        fit_param = self.calc_time_fit_residual(bounds)
        # parsing the fitting results
        fit_vec = np.multiply(
            fit_param[0, :], np.exp(1j * np.pi * fit_param[-1, :] / 180.0)
        )
        fit_area = abs(fit_vec)
        fit_freq = fit_param[1, :]
        fit_phase = np.arctan2(np.imag(fit_vec), np.real(fit_vec)) * 180.0 / np.pi

        if self.method == "voigt":
            fit_fwhmL = fit_param[2, :]
            fit_fwhmG = fit_param[3, :]

            self.set_components(
                area=fit_area,
                freq=fit_freq,
                fwhmL=fit_fwhmL,
                fwhmG=fit_fwhmG,
                phase=fit_phase,
            )
        elif self.method == "lorentzian":
            fit_fwhm = fit_param[2, :]
            self.set_components(
                area=fit_area, freq=fit_freq, fwhm=fit_fwhm, phase=fit_phase
            )
        else:
            raise ValueError("Unknown fitting method.")

    def plot_time_spect_fit(self):
        """Plot the time domain and spectral domain fitting results."""
        complex_fit_time = self.get_time_function(tdata=self.tdata)
        plt.figure(figsize=(15, 5))
        plt.subplot(221)

        ax1 = plt.subplot(1, 3, 1)
        ax1.plot(self.tdata, abs(self.ydata))
        ax1.plot(self.tdata, abs(complex_fit_time))

        ax1.legend(["broad time sig", "fit time sig"])

        # calculate fit spectral signal
        complex_fit_spect = self.sample_time * np.fft.fftshift(
            np.fft.fft(complex_fit_time, self.zeropad_size)
        )

        ax2 = plt.subplot(1, 3, 2)
        ax2.plot(self.f, abs(self.spectral_signal), "*-")
        ax2.plot(self.f, abs(complex_fit_spect))
        ax2.set_xlim((-10000, 10000))
        ax2.legend(["spectral sig", "fit spect sig"])
        plt.show()
