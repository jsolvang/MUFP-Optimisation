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
    #df_rads = df*2*np.pi
    #A = np.sqrt(2*s*df_rads)

    #s[s < 0.01 * max(s)] = 0
    #A[A < 0.01 * max(A)] = 0.01 * max(A)



    if plot == 1:
        plt.rcParams['axes.grid'] = True
        plt.rcParams["figure.figsize"] = (15, 5)
        plt.rcParams.update({'font.size': 15})
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['lines.markersize'] = 18
        plt.plot(f, s, linewidth=2, markersize=5)
        plt.xlabel('Frequency [rad/s]')
        plt.ylabel('Spectral Density [$m^2/rad/s$]')
        plt.grid(b=True, which='both', axis='both')
        plt.xlim(0, 5)
        plt.tight_layout()
    return s, A