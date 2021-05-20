import struct
import numpy as np

# Path + file name
path = r"C:\Users\Joar\Documents\1_Education\NTNU\OneDrive - NTNU\Thesis\Modelling\SIMA\FixedRotor5MW\Initial\151-20210417101939"
# Total simulation time and storage time step
simtime = 600.
dt = 0.1
transient = 100

# Element number in storage list ()
n_beam = 1

# ---------------------------------------------RIFLEX results---------------------------------------------#
# Read .bin file
with open(path + r'\sima_elmfor.bin', mode='rb') as file:
    fileContent = file.read()

# Reshape results in file to a 2D array 'A', with one result component per row and one time step per column
cols = len(np.asarray(struct.unpack("f" * (len(fileContent) // 4), fileContent))) / (simtime / dt)
A = np.reshape(np.asarray(struct.unpack("f" * (len(fileContent) // 4), fileContent)), (int(cols), int(simtime / dt)),
               order='F')

# Time vector
t = A[1]

# NB! Forces from RIFLEX are in kN!

# Axial force
N = A[2 + 10 * (n_beam - 1)] * 1000.0

# Torsional moment
Mx = A[3 + 10 * (n_beam - 1)] * 1000.0

# Mom. about local y-axis, end 1 and 2
My_end1 = A[4 + 10 * (n_beam - 1)] * 1000.0
My_end2 = A[5 + 10 * (n_beam - 1)] * 1000.0

# Mom. about local z-axis, end 1 and 2
Mz_end1 = A[6 + 10 * (n_beam - 1)] * 1000.0
Mz_end2 = A[7 + 10 * (n_beam - 1)] * 1000.0

# Shear force in local y-direction, end 1 and 2
Qy_end1 = A[8 + 10 * (n_beam - 1)] * 1000.0
Qy_end2 = A[9 + 10 * (n_beam - 1)] * 1000.0

# Shear force in local z-direction, end 1 and 2
Qz_end1 = A[10 + 10 * (n_beam - 1)] * 1000.0
Qz_end2 = A[11 + 10 * (n_beam - 1)] * 1000.0

#---------------------------------------------SIMO results---------------------------------------------#

# Read .bin file
with open(path + r'\results.tda', mode='rb') as file:
    CC = file.read()

# Body number
n_body = 1

# Reshape results in file to a 2D array 'A', with one result component per row and one time step per column
cols = len(np.asarray(struct.unpack("f" * (len(fileContent) // 4), fileContent))) / (simtime / dt)
A = np.transpose(np.reshape(np.asarray(struct.unpack("f" * (len(fileContent) // 4), fileContent)), (int(simtime / dt), int(cols)), order='F'))

# Time vector
t = A[1]

# Translations (m) and rotations (deg) for body
Tx = A[4+8*(n_body-1)]
Ty = A[5+8*(n_body-1)]
Tz = A[6+8*(n_body-1)]
Rx = A[7+8*(n_body-1)]
Ry = A[8+8*(n_body-1)]
Rz = A[9+8*(n_body-1)]

# ---------------------------------------------Wind Turbine Results---------------------------------------------#
with open(path + r'\sima_witurb.bin', mode='rb') as file:
    CC = file.read()

# Element number in storage list ()
n_beam = 1

cols = len(np.asarray(struct.unpack("f" * (len(CC) // 4), CC))) / (simtime / dt)
CC = np.transpose(
    np.reshape(np.asarray(struct.unpack("f" * (len(CC) // 4), CC)), (int(cols), int(simtime / dt)), order='F'))

time = CC[:, 1]
omega = CC[:, 2] * np.pi / 180  # convert from deg/s to rad/s
genTq = CC[:, 4]
genPwr = CC[:, 5]
azimuth = CC[:, 6]
HubWindX = CC[:, 7]
HubWindY = CC[:, 8]
HubWindZ = CC[:, 9]
AeroForceX = CC[:, 10] *1e3  # converting to N
AeroMomX = CC[:, 13] *1e3 # converting to N
Bl1Pitch = CC[:, 16]
Bl2Pitch = CC[:, 17]
Bl3Pitch = CC[:, 18]