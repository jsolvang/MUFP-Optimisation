import numpy as np
import matplotlib.pyplot as plt


def jonswap(Hs, Tp, df, gammaJS, f, plot):
    df = df / (2*np.pi)         # Converting from rads-1 to hz
    f = f / (2 * np.pi)         # Converting from rads-1 to hz

    fp = 1/Tp
    s = np.ones(len(f))
    A = np.ones(len(f))
    for ii in np.linspace(0,len(f)-1,len(f)).astype(int):
        if f[ii] <= fp:
            sigma = 0.07
        elif f[ii] > fp:
            sigma = 0.09
        s[ii] = 0.3125*np.square(Hs)*Tp*np.power(np.divide(f[ii],fp),-5)*np.exp(-1.25*np.power(np.divide(f[ii],fp),-4))*((1-0.287*np.log(gammaJS)))*np.power(gammaJS,np.exp(-0.5*np.square(np.divide(np.divide(f[ii],fp)-1,sigma))))

    # Converting to rad/s
    s = np.divide(s,2*np.pi)
    f = f*2*np.pi
    df_rads = df*2*np.pi
    A = np.sqrt(2*s*df_rads)

    s[s < 0.01 * max(s)] = 0
    A[A < 0.01 * max(A)] = 0



    if plot == 1:
        plt.rcParams["figure.figsize"] = (15, 10)
        fig, axs = plt.subplots(1, 2)
        axs[0].plot(f, s, linewidth=2, markersize=5)
        axs[0].set(xlabel='Frequency [rad/s]', ylabel='Spectral Density [M^2/rad/s]')
        axs[0].grid(b=True, which='both', axis='both')
        axs[1].plot(f, A, linewidth=2, markersize=5)
        axs[1].set(xlabel='Frequency [rad/s]', ylabel='Amplitude [m]')
        axs[1].grid(b=True, which='both', axis='both')
        axs[0].set(xlim=(0, 5))
        axs[1].set(xlim=(0, 5))
        plt.tight_layout()
    return s, A