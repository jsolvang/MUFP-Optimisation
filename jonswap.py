import numpy as np
import matplotlib.pyplot as plt


def jonswap(Hs, Tp, gammaJS, f, plot):
    global axs
    fp = 1/Tp
    s = np.ones(len(f))
    A = np.ones(len(f))
    for ii in np.linspace(0,len(f)-1,len(f)).astype(int):
        if f[ii] <= fp:
            sigma = 0.07
        elif f[ii] > fp:
            sigma = 0.09
        s[ii] = 0.3125*np.square(Hs)*Tp*np.power(np.divide(f[ii],fp),-5)*np.exp(-1.25*np.power(np.divide(f[ii],fp),-4))*((1-0.287*np.log(gammaJS)))*np.power(gammaJS,np.exp(-0.5*np.square(np.divide(np.divide(f[ii],fp)-1,sigma))))
        A[ii] = np.sqrt(2*s[ii]*(1/3600))

    if plot == 1:
        plt.rcParams["figure.figsize"] = (15, 10)
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(f, s, linewidth=2, markersize=5)
        axs[0].set(xlabel='Frequency [Hz]', ylabel='Spectral Density [W/Hz]')
        axs[0].grid(b=True, which='both', axis='both')
        axs[1].plot(f, A, linewidth=2, markersize=5)
        axs[1].set(xlabel='Frequency [Hz]', ylabel='Amplitude [m]')
        axs[1].grid(b=True, which='both', axis='both')
        plt.tight_layout()
    return s, A