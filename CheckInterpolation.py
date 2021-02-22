import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
import pickle


class InterpolateParameters:
    def __init__(self, pull, interp_at):
        results = []

        for i in np.linspace(0, len(pull[:, 1]) - 1, 3).astype(int):
            file_loc = r'C:\Users\Joar\Documents\1_Education\NTNU\pickle_files'
            file_name = "\sim_x_%f_y_%f_D%f" % (pull[i, 0], pull[i, 1], pull[i, 2])
            file_path = file_loc + file_name
            infile = open(file_path, 'rb')
            results.append(pickle.load(infile))
            infile.close()

        interp_ADDEDMASS = np.zeros(shape=(int(results[1].numwavelengths), 6, 6))
        interp_DAMPING = np.zeros(shape=(int(results[1].numwavelengths), 6, 6))
        interp_WAVEEX = np.zeros(shape=(2, int(results[1].numwavelengths), 6, 6))

        for nn in np.linspace(0, int(results[1].numwavelengths) - 1, int(results[1].numwavelengths)).astype(int):
            for mm in np.linspace(0, 5, 6).astype(int):
                for kk in np.linspace(0, 5, 6).astype(int):
                    interp_ADDEDMASS[nn, mm, kk] = np.interp(interp_at[1],
                                                             [pull[0, interp_at[0]], pull[2, interp_at[0]]],
                                                             [results[0].ADDEDMASS[nn, mm, kk],
                                                              results[2].ADDEDMASS[nn, mm, kk]])
                    interp_DAMPING[nn, mm, kk] = np.interp(interp_at[1],
                                                           [pull[0, interp_at[0]], pull[2, interp_at[0]]],
                                                           [results[0].DAMPING[nn, mm, kk],
                                                            results[2].DAMPING[nn, mm, kk]])

        for nn in np.linspace(0, int(results[0].numwavelengths) - 1, int(results[0].numwavelengths)).astype(int):
            for mm in np.linspace(0, 5, 6).astype(int):
                for kk in np.linspace(0, 3, 4).astype(int):
                    interp_WAVEEX[0, nn, mm, kk] = np.interp(interp_at[1],
                                                             [pull[0, interp_at[0]], pull[2, interp_at[0]]],
                                                             [results[0].WAVEEX[0, nn, mm, kk],
                                                              results[2].WAVEEX[0, nn, mm, kk]])
                    interp_WAVEEX[1, nn, mm, kk] = np.interp(interp_at[1],
                                                             [pull[0, interp_at[0]], pull[2, interp_at[0]]],
                                                             [results[0].WAVEEX[1, nn, mm, kk],
                                                              results[2].WAVEEX[1, nn, mm, kk]])
        self.results = results
        self.ADDEDMASS = interp_ADDEDMASS
        self.DAMPING = interp_DAMPING
        self.WAVEEX = interp_WAVEEX

    def _plot_interpolation(self, pull, interp_at):

        #################################################################################
        results = self.results

        ADDEDMASS_yaxis = [['Added Mass [KG] (DOF:Surge)'],
                            ['Added Mass [KG] (DOF:Sway)'],
                            ['Added Mass [KG] (DOF:Heave)'],
                            ['Added Mass [KG-m] (DOF:Roll)'],
                            ['Added Mass [KG-m] (DOF:Pitch)'],
                            ['Added Mass [KG-m] (DOF:Yaw)']]

        plt.rcParams["figure.figsize"] = (20, 20)
        fig, axs = plt.subplots(6, 1)
        for ii in np.linspace(0, 5, 6).astype(int):
            axs[ii].plot(results[0].wave_disc[:, 4], results[0].ADDEDMASS[:, ii, ii], '-+m',
                         label="sim_x_%f_y_%f_D%f" % (pull[0, 0], pull[0, 1], pull[0, 2]),
                         linewidth=1, markersize=5)
            axs[ii].plot(results[1].wave_disc[:, 4], results[1].ADDEDMASS[:, ii, ii], '-or',
                         label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[1, 2]),
                         linewidth=2, markersize=5)
            axs[ii].plot(results[2].wave_disc[:, 4], results[2].ADDEDMASS[:, ii, ii], '-og',
                         label="sim_x_%f_y_%f_D%f" % (pull[2, 0], pull[2, 1], pull[2, 2]),
                         linewidth=1, markersize=5)
            axs[ii].plot(results[1].wave_disc[:, 4], self.ADDEDMASS[:, ii, ii], '--ob',
                         label="Interpolated at %f" % interp_at[1],
                         linewidth=2, markersize=5)
            axs[ii].legend(loc="upper right")
            axs[ii].set(ylabel=ADDEDMASS_yaxis[ii])
            axs[ii].grid(b=True, which='both', axis='both')

        #################################################################################

        plt.rcParams["figure.figsize"] = (20, 20)
        fig, axs = plt.subplots(6, 1)
        for ii in np.linspace(0, 5, 6).astype(int):
            axs[ii].plot(results[1].wave_disc[:, 4], results[1].ADDEDMASS[:, ii, ii], '-or',
                         label="sim_x_%f_y_%f_D%f" % (pull[0, 1], pull[1, 1], pull[2, 1]),
                         linewidth=2, markersize=5)
            axs[ii].plot(results[1].wave_disc[:, 4], self.ADDEDMASS[:, ii, ii], '--ob',
                         label="Interpolated at %f" % interp_at[1],
                         linewidth=2, markersize=5)
            axs[ii].legend(loc="upper right")
            axs[ii].set(ylabel=ADDEDMASS_yaxis[ii])
            axs[ii].grid(b=True, which='both', axis='both')

        #################################################################################

        Damping_yaxis = [['Damping [kg/s] (DOF:Surge)'],
                            ['Damping [kg/s] (DOF:Sway)'],
                            ['Damping [kg/s] (DOF:Heave)'],
                            ['Damping [kg-m/s] (DOF:Roll)'],
                            ['Damping [kg-m/s] (DOF:Pitch)'],
                            ['Damping [kg-m/s] (DOF:Yaw)']]

        plt.rcParams["figure.figsize"] = (20, 20)
        fig, axs = plt.subplots(6, 1)
        for ii in np.linspace(0, 5, 6).astype(int):
            axs[ii].plot(results[0].wave_disc[:, 4], results[0].DAMPING[:, ii, ii], '-+m',
                         label="sim_x_%f_y_%f_D%f" % (pull[0, 0], pull[0, 1], pull[0, 2]),
                         linewidth=2, markersize=5)
            axs[ii].plot(results[1].wave_disc[:, 4], results[1].DAMPING[:, ii, ii], '-or',
                         label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[1, 2]),
                         linewidth=2, markersize=5)
            axs[ii].plot(results[2].wave_disc[:, 4], results[2].DAMPING[:, ii, ii], '-og',
                         label="sim_x_%f_y_%f_D%f" % (pull[2, 0], pull[2, 1], pull[2, 2]),
                         linewidth=2, markersize=5)
            axs[ii].plot(results[2].wave_disc[:, 4], self.DAMPING[:, ii, ii], '--ob',
                         label="Interpolated at %f" % interp_at[1],
                         linewidth=2, markersize=5)
            axs[ii].legend(loc="upper right")
            axs[ii].set(ylabel=Damping_yaxis[ii])
            axs[ii].grid(b=True, which='both', axis='both')

        ##################################################################################

        plt.rcParams["figure.figsize"] = (20, 20)
        fig, axs = plt.subplots(6, 1)
        for ii in np.linspace(0, 5, 6).astype(int):
            axs[ii].plot(results[1].wave_disc[:, 4], results[1].DAMPING[:, ii, ii], '-or',
                         label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[1, 2]),
                         linewidth=2, markersize=5)
            axs[ii].plot(results[2].wave_disc[:, 4], self.DAMPING[:, ii, ii], '--ob',
                         label="Interpolated at %f" % interp_at[1],
                         linewidth=2, markersize=5)
            axs[ii].legend(loc="upper right")
            axs[ii].set(ylabel=Damping_yaxis[ii])
            axs[ii].grid(b=True, which='both', axis='both')

        ####################################################################################################

        fig, axs = plt.subplots(2, 2)
        for jj in np.linspace(0, 2, 3).astype(int):
            axs[0, 0].plot(results[0].wave_disc[:, 4], results[0].WAVEEX[0, :, jj, 2], '-+g',
                           label="sim_x_%f_y_%f_D%f" % (pull[0, 0], pull[0, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 0].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[0, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 0].plot(results[2].wave_disc[:, 4], results[2].WAVEEX[0, :, jj, 2], '-+m',
                           label="sim_x_%f_y_%f_D%f" % (pull[2, 0], pull[2, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 0].plot(results[2].wave_disc[:, 4], self.WAVEEX[0, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[0, 0].legend(loc="upper right")
            axs[0, 0].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force [N]')
            axs[0, 0].grid(b=True, which='both', axis='both')

        for jj in np.linspace(3, 5, 3).astype(int):
            axs[0, 1].plot(results[0].wave_disc[:, 4], results[0].WAVEEX[0, :, jj, 2], '-+g',
                           label="sim_x_%f_y_%f_D%f" % (pull[0, 0], pull[0, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 1].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[0, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 1].plot(results[2].wave_disc[:, 4], results[2].WAVEEX[0, :, jj, 2], '-+m',
                           label="sim_x_%f_y_%f_D%f" % (pull[2, 0], pull[2, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 1].plot(results[2].wave_disc[:, 4], self.WAVEEX[0, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            #axs[0, 1].legend(loc="upper right")
            axs[0, 1].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force [N]')
            axs[0, 1].grid(b=True, which='both', axis='both')

        for jj in np.linspace(0, 2, 3).astype(int):
            axs[1, 0].plot(results[0].wave_disc[:, 4], results[0].WAVEEX[1, :, jj, 2], '-+g',
                           label="sim_x_%f_y_%f_D%f" % (pull[0, 0], pull[0, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 0].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[1, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 0].plot(results[2].wave_disc[:, 4], results[2].WAVEEX[1, :, jj, 2], '-+m',
                           label="sim_x_%f_y_%f_D%f" % (pull[2, 0], pull[2, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 0].plot(results[2].wave_disc[:, 4], self.WAVEEX[1, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[1, 0].legend(loc="upper right")
            axs[1, 0].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force [N]')
            axs[1, 0].grid(b=True, which='both', axis='both')

        for jj in np.linspace(3, 5, 3).astype(int):
            axs[1, 1].plot(results[0].wave_disc[:, 4], results[0].WAVEEX[1, :, jj, 2], '-+g',
                           label="sim_x_%f_y_%f_D%f" % (pull[0, 0], pull[0, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 1].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[1, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 1].plot(results[2].wave_disc[:, 4], results[2].WAVEEX[1, :, jj, 2], '-+m',
                           label="sim_x_%f_y_%f_D%f" % (pull[2, 0], pull[2, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 1].plot(results[2].wave_disc[:, 4], self.WAVEEX[1, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[1, 1].legend(loc="upper right")
            axs[1, 1].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force [N]')
            axs[1, 1].grid(b=True, which='both', axis='both')

        ######################################################################################

        fig, axs = plt.subplots(2, 2)
        for jj in np.linspace(0, 2, 3).astype(int):
            axs[0, 0].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[0, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 0].plot(results[2].wave_disc[:, 4], self.WAVEEX[0, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[0, 0].legend(loc="upper right")
            axs[0, 0].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force [N]')
            axs[0, 0].grid(b=True, which='both', axis='both')

        for jj in np.linspace(3, 5, 3).astype(int):
            axs[0, 1].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[0, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[0, 1].plot(results[2].wave_disc[:, 4], self.WAVEEX[0, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[0, 1].legend(loc="upper right")
            axs[0, 1].set(xlabel='omega, rad/s', ylabel='0 Deg Exciting Force [N]')
            axs[0, 1].grid(b=True, which='both', axis='both')

        for jj in np.linspace(0, 2, 3).astype(int):
            axs[1, 0].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[1, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 0].plot(results[2].wave_disc[:, 4], self.WAVEEX[1, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[1, 0].legend(loc="upper right")
            axs[1, 0].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force [N]')
            axs[1, 0].grid(b=True, which='both', axis='both')

        for jj in np.linspace(3, 5, 3).astype(int):
            axs[1, 1].plot(results[1].wave_disc[:, 4], results[1].WAVEEX[1, :, jj, 2], '-+b',
                           label="sim_x_%f_y_%f_D%f" % (pull[1, 0], pull[1, 1], pull[0, 2]),
                           linewidth=2, markersize=5)
            axs[1, 1].plot(results[2].wave_disc[:, 4], self.WAVEEX[1, :, jj, 2], '-+r',
                           label="Interpolated at %f" % interp_at[1],
                           linewidth=2, markersize=5)
            axs[1, 1].legend(loc="upper right")
            axs[1, 1].set(xlabel='omega, rad/s', ylabel='90 Deg Exciting Force [N]')
            axs[1, 1].grid(b=True, which='both', axis='both')
