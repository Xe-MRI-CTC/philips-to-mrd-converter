"""NMR mix class."""

import numpy as np


class NMR_Mix:
    """Base Class for curve fitting for spectroscopy.

    A class that represents a series of exponentially decaying components.
    The class knows how to calculate time or spectral domain signals. The
    key assumption is that all components experience perfectly exponential
    decay.

        Attributes:
            area (np.ndarray): Area of each component.
            freq (np.ndarray): Frequency of each component.
            phase (np.ndarray): Phase of each component.
            fwhmL (np.ndarray): Lorentzian FWHM of each component.
            fwhmG (np.ndarray): Gaussian FWHM of each component.
            method (str): Method for fitting, only supports "voigt".
            ncomp (int): Number of components.
    """

    def __init__(
        self,
        area: np.ndarray,
        freq: np.ndarray,
        phase: np.ndarray,
        fwhmL: np.ndarray = np.array([]),
        fwhmG: np.ndarray = np.array([]),
        method="voigt",
    ):
        """Initialize NMR_Mix class."""
        self.area = np.array([area]).flatten()
        self.freq = np.array([freq]).flatten()
        self.phase = np.array([phase]).flatten()
        self.fwhmL = np.array([fwhmL]).flatten()
        self.fwhmG = np.array([fwhmG]).flatten()
        self.method = method
        self.ncomp = 5

    def get_time_function(self, tdata: np.ndarray):
        """Get time function for given time points.

        Available for both voigt and lorentzian fitting.

        Args:
            t (np.ndarray): Time points in seconds.
        Returns: Time domain signal.
        """
        n_fre = int(np.size(self.freq))
        assert n_fre == 3, "Number of components must be 3."
        time_sig = np.zeros(np.shape(tdata))
        if self.method == "voigt":
            time_sig = time_sig + self.area[0] * np.exp(
                1j * np.pi / 180.0 * self.phase[0]
                + 1j * 2 * np.pi * tdata * self.freq[0]
            ) * np.exp(-np.pi * tdata * self.fwhmL[0])
            for k in range(1, n_fre):
                time_sig = time_sig + self.area[k] * np.exp(
                    1j * np.pi / 180.0 * self.phase[k]
                    + 1j * 2 * np.pi * tdata * self.freq[k]
                ) * np.exp(-(tdata**2) * 4 * np.log(2) * self.fwhmG[k] ** 2) * np.exp(
                    -np.pi * tdata * self.fwhmL[k]
                )
        else:
            raise ValueError("Method must be either 'voigt' or 'lorentzian'.")
        return time_sig

    def get_init_params(self):
        """Get initial parameters for fitting."""
        return np.concatenate(
            (self.area, self.freq, self.phase, self.fwhmL, self.fwhmG)
        ).flatten()

    def sort_freq(self):
        """Sort components according to resonance frequency in descending order."""
        sort_index = np.argsort(-self.freq)
        self.freq = self.freq[sort_index]
        self.area = self.area[sort_index]
        self.phase = self.phase[sort_index]

        if self.method == "voigt":
            self.fwhmG = self.fwhmG[sort_index]
            self.fwhmL = self.fwhmL[sort_index]
        else:
            self.fwhm = self.fwhm[sort_index]

    def set_components(
        self,
        area: np.ndarray,
        freq: np.ndarray,
        phase: np.ndarray,
        fwhm: np.ndarray = np.array([]),
        fwhmL: np.ndarray = np.array([]),
        fwhmG: np.ndarray = np.array([]),
    ):
        """Set components and sort frequencies in descending order."""
        self.area = np.array([area]).flatten()
        self.freq = np.array([freq]).flatten()
        self.fwhm = np.array([fwhm]).flatten()
        self.phase = np.array([phase]).flatten()
        self.fwhmG = np.array([fwhmG]).flatten()
        self.fwhmL = np.array([fwhmL]).flatten()
        self.sort_freq()

    def calcComponentSpectralDomainSignal(self, f):
        # Calculates the spectral domain signal from the individual components
        # of the NMR_mix at the given spectral frequencies (f is in Hz).
        # Note: This function returns the spectrum for each individual component.
        # The overall "mix" spectrum can be obtained with calcSpectralDomainSignal.

        nFreqSamples = len(f)
        '''
        nComp = 1 #len(self.area)
        componentSpectralDomainSignal = np.zeros((nFreqSamples, nComp), dtype=complex)

        for iComp in range(nComp):
            componentSpectralDomainSignal[:, iComp] = \
                self.area[iComp] * np.exp(1j * np.pi / 180 * self.phase[iComp]) / \
                (1j * 2 * np.pi * (f - self.freq[iComp]) + np.pi * self.fwhmL[iComp])
        '''
        nComp = 1
        componentSpectralDomainSignal = np.zeros((nFreqSamples, nComp), dtype=complex)
        componentSpectralDomainSignal = \
            self.area * np.exp(1j * np.pi / 180 * self.phase) / \
            (1j * 2 * np.pi * (f - self.freq) + np.pi * self.fwhmL)

        return componentSpectralDomainSignal

    def calcComponentSpectralDomainSignal_dis(self, f):
        # Calculates the spectral domain signal from the individual components
        # of the NMR_mix at the given spectral frequencies (f is in Hz).
        # Note: This function returns the spectrum for each individual component.
        # The overall "mix" spectrum can be obtained with calcSpectralDomainSignal.
        nFreqSamples = len(f)
        nComp = len(self.area)
        componentSpectralDomainSignal = np.zeros((nFreqSamples, nComp), dtype=complex)

        for iComp in range(nComp):
            componentSpectralDomainSignal[:, iComp] = \
                self.area[iComp] * np.exp(1j * np.pi / 180 * self.phase[iComp]) / \
                (1j * 2 * np.pi * (f - self.freq[iComp]) + np.pi * self.fwhmL[iComp])

        return componentSpectralDomainSignal

    def calcComponentTimeDomainSignal(self, t):
        # Calculates the time domain signal from the individual components
        # of the NMR mix at the given time points (t is in seconds).
        # Note: This function returns the time signal for each individual
        # component. The overall "mix" signal can be obtained with
        # calcTimeDomainSignal.

        nTimePoints = len(t)
        nComp = len(self.area)
        componentTimeDomainSignal = np.zeros((nTimePoints, nComp), dtype=complex)

        # First component
        componentTimeDomainSignal[:, 0] = self.area[0] * \
            np.exp(1j * np.pi / 180 * self.phase[0] + 1j * 2 * np.pi * t * self.freq[0]) \
            * np.exp(-np.pi * t * self.fwhmL[0])

        # Remaining components
        for iComp in range(1, nComp):
            componentTimeDomainSignal[:, iComp] = self.area[iComp] * \
                np.exp(1j * np.pi / 180 * self.phase[iComp] + 1j * 2 * np.pi * t * self.freq[iComp]) \
                * np.exp(-t**2 * 4 * np.log(2) * (self.fwhmG[iComp])**2) \
                * np.exp(-np.pi * t * self.fwhmL[iComp])

        return componentTimeDomainSignal

    def calcSpectralDomainSignal(self, f):
        # Calculates the spectral domain signal from the mix of NMR components
        # at the given spectral frequencies (f is in Hz).
        # Note: This function returns the net spectrum from the sum of all components
        # in the NMR_Mix. The individual component spectra can be obtained with
        # calcComponentSpectralDomainSignal.

        nFreqSamples = len(f)
        nComp = len(self.area)
        spectralDomainSignal = np.zeros(nFreqSamples, dtype=complex)

        for iComp in range(nComp):
            spectralDomainSignal += self.area[iComp] * np.exp(1j * np.pi / 180 * self.phase[iComp]) / \
                (1j * 2 * np.pi * (f - self.freq[iComp]) + np.pi * self.fwhmL[iComp])

        return spectralDomainSignal
