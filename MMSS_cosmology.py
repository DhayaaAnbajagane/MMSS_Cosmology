#Load all required packages

import numpy as np #The holy grail of scientific computing!
import matplotlib.pyplot as plt #The only reason I can call myself an "artist"

#Set plotting params
import matplotlib as mpl
mpl.rcParams['xtick.direction'], mpl.rcParams['ytick.direction'] = 'in', 'in'
mpl.rcParams['xtick.major.size'], mpl.rcParams['xtick.minor.size'] = 14, 8
mpl.rcParams['xtick.major.width'], mpl.rcParams['xtick.minor.width'] = 1.2, 0.8
mpl.rcParams['xtick.major.pad'], mpl.rcParams['xtick.minor.pad'] = 10, 10
mpl.rcParams['ytick.major.size'], mpl.rcParams['ytick.minor.size'] = 14, 8
mpl.rcParams['ytick.major.width'], mpl.rcParams['ytick.minor.width'] = 1.2, 0.8

#Command to compute Omega
def Omega_X(a, Omega_0, w):
    '''
    Computes the fraction of energy in the Universe due to a component X.

    Input
    --------
    a: numpy array
        All the values of the scale factor, a, at which we compute Omega_X

    Omega_0: float
        The value of Omega_X at present (a = 1)

    w: float, int, or function of a
        The equation of state w = pressure/energy of the component X.
        In this project, we make it a function of scale factor a

    Output
    --------

    numpy array:
        The value of Omega_X at all values in the array a
    '''

    #Initialize numpy array to store our output
    factor = np.zeros(a.size)

    #Iterate/loop over all input values of a
    for i, a_now in enumerate(a):

        #For each a_now, create an array of values from a = a_now to a = 1
        present_epoch_to_a_now = np.geomspace(a_now, 1, 1000, endpoint=True)

        #Integrate
        factor[i] = np.exp(3*np.trapz(1/present_epoch_to_a_now*(1 + w(present_epoch_to_a_now)), present_epoch_to_a_now))

    return Omega_0*factor

#Compute Hubble constant given list of Omegas
def Hubble(list_of_Omegas):

    H0 = 70 #Fiducial value of hubble constant in km/s/Mpc

    return H0*np.sqrt(np.sum(list_of_Omegas, axis = 0))

def Our_Universe(a_start = 1e-7, a_end = 1e2, Omega_m = 0.3, Omega_r = 2e-5, Omega_L = 0.7, w_m = 0, w_r = 1/3, w_Lambda = -1):

    a = np.geomspace(a_start, a_end, 500, endpoint=True)

    Omega_M = Omega_X(a, Omega_m,  lambda x: w_m)
    Omega_R = Omega_X(a, Omega_r, lambda x: w_r)
    Omega_L = Omega_X(a, Omega_L,  lambda x: w_Lambda)

    Sum = Omega_m + Omega_r + Omega_L

    print("-----------------------")
    print("Our Initial Conditions")
    print("-----------------------")
    print("Omega_m:%0.2f"%(Omega_m/Sum))
    print("Omega_r:%0.2f"%(Omega_r/Sum))
    print("Omega_L:%0.2f"%(Omega_L/Sum))

    plt.rc('xtick',labelsize=22)
    plt.rc('ytick',labelsize=22)

    fig, axes = plt.subplots(2, 1, figsize = (12,14), sharex=True)
    plt.subplots_adjust(hspace = 0.05)
    axes[1].set_yscale('log')
    [ax.set_xscale('log') for ax in axes]

    axes[0].plot(a, Omega_M/Normalize, lw = 2, label = r'$\Omega_m,\, w = ' + str(w_m) + '$')
    axes[0].plot(a, Omega_R/Normalize, lw = 2, label = r'$\Omega_r,\, w = ' + str(w_r) + '$')
    axes[0].plot(a, Omega_L/Normalize, lw = 2, label = r'$\Omega_\Lambda,\, w = ' + str(w_Lambda) + '$')

    axes[0].legend(fontsize = 20)
    axes[0].set_ylabel(r'$\Omega_i(a)$', size = 30)
    axes[0].grid()

    list_of_Omegas = np.array([Omega_M, Omega_R, Omega_L])
    axes[1].plot(a, Hubble(list_of_Omegas), lw = 3, alpha = 1, color = 'C6')

    axes[1].legend(fontsize = 24)
    axes[1].set_ylabel(r'$H(a)$', size = 30)
    axes[1].set_xlabel(r'$a$', size = 30)
    axes[1].grid()

def New_Universe(w_new, Omega_new, a_start = 1e-7, a_end = 1e2, Omega_m = 0.3, Omega_r = 2e-5, Omega_L = 0.7, w_m = 0, w_r = 1/3, w_Lambda = -1):

    a = np.geomspace(a_start, a_end, 500, endpoint=True)

    Omega_M = Omega_X(a, Omega_m,  lambda x: w_m)
    Omega_R = Omega_X(a, Omega_r,  lambda x: w_r)
    Omega_L = Omega_X(a, Omega_L,  lambda x: w_Lambda)
    Omega_new = Omega_X(a, Omega_new, w_new)

    Normalize = (Omega_M + Omega_R + Omega_L + Omega_new)

    Sum = Omega_m + Omega_r + Omega_L + Omega_new

    print("-----------------------")
    print("Our Initial Conditions")
    print("-----------------------")
    print("Omega_m:%0.2f"%(Omega_m/Sum))
    print("Omega_r:%0.2f"%(Omega_r/Sum))
    print("Omega_L:%0.2f"%(Omega_L/Sum))
    print("Omega_new:%0.2f"%(Omega_new/Sum))

    plt.rc('xtick',labelsize=22)
    plt.rc('ytick',labelsize=22)

    fig, axes = plt.subplots(2, 1, figsize = (12,14), sharex=True)
    plt.subplots_adjust(hspace = 0.05)
    axes[1].set_yscale('log')
    [ax.set_xscale('log') for ax in axes]

    axes[0].plot(a, Omega_M/Normalize, lw = 2, label = r'$\Omega_m,\, w = ' + str(w_m) + '$')
    axes[0].plot(a, Omega_R/Normalize, lw = 2, label = r'$\Omega_r,\, w = ' + str(w_r) + '$')
    axes[0].plot(a, Omega_L/Normalize, lw = 2, label = r'$\Omega_\Lambda,\, w = ' + str(w_Lambda) + '$')
    axes[0].plot(a, Omega_new/Normalize, lw = 2, label = r'$\Omega_{\rm new}$')

    axes[0].legend(fontsize = 20)
    axes[0].set_ylabel(r'$\Omega_i(a)$', size = 30)
    axes[0].grid()

    list_of_Omegas = np.array([Omega_M, Omega_R, Omega_L])
    axes[1].plot(a, Hubble(list_of_Omegas), lw = 3, label = r'Normal', alpha = 0.3, color = 'C4')

    list_of_Omegas = np.array([Omega_M, Omega_R, Omega_L, Omega_new])
    axes[1].plot(a, Hubble(list_of_Omegas), lw = 3, label = r'New!', color = 'C6')

    axes[1].legend(fontsize = 24)
    axes[1].set_ylabel(r'$H(a)$', size = 30)
    axes[1].set_xlabel(r'$a$', size = 30)
    axes[1].grid()
