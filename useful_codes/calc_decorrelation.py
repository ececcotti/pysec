#! /usr/bin/env python

from argparse import ArgumentParser
import sys, os
import math
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
from scipy import special


parser = ArgumentParser(description="Estimate peak response to a point source at a given distance from the phase center due to frequency and/or time averaging")
parser.add_argument('-b', dest="bsl", type=float, metavar='<max_bsl>', help='maximum baseline lenght', required=True)
parser.add_argument('-t', nargs='+', dest="time", type=float, metavar='<delta_t>', help='time interval for averaging')
parser.add_argument('-f', nargs='+', dest="freq", type=float, metavar='<delta_f>', help='frequency interval for averaging')
parser.add_argument('-i', dest="freq_0", type=float, metavar='<f_0>', help='central frequency of the interval')
parser.add_argument('-r', dest="offset", type=float, metavar='<theta>', help='offset for which you want the attenuation factor')

def freq_smearing(f0, delta_f, b, r): #units: time=sec, freq=Hz, dist=deg
    c = 3.e8
    theta_b = 1.410 * c / f0 / b
    R = 1.0645 * theta_b * f0 / (r * delta_f) * special.erf(0.8326 * r * delta_f / (theta_b * f0))
    return R

def time_smearing(f0, delta_t, b, r):
    c = 3.e8
    omega_e = 7.27e-5
    theta_b = 1.410 * c / f0 / b
    R = 1.0645 * theta_b / (r * omega_e * delta_t) * special.erf(0.8326 * r * omega_e * delta_t / theta_b)
    return R




def main(argv):
    args=parser.parse_args(argv)

    mpl.style.use('default')
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True

    fig, ax = plt.subplots(1, 2, sharey=True, figsize=(12,5))
    fig.subplots_adjust(wspace=0)

    for delta_f in args.freq:
        Rb = []
        for rx in r:
            Rb.append(freq_smearing(args.freq_0, delta_f, args.bsl, math.radians(rx)))
    
        ax[0].plot(r, Rb, label=r"$\Delta\nu = $%.1f kHz"%(delta_f/1.e3))
        ax[0].set_title('Frequency smearing (baseline = %i km)' % (args.bsl/1.e3)),

    for delta_t in args.time:
        Ra = []
        for rx in r:
            Ra.append(time_smearing(args.freq_0, delta_t, args.bsl, math.radians(rx)))
    
        ax[1].set_title('Time smearing (baseline = %i km)' % (args.bsl/1.e3))
        ax[1].plot(r, Ra, label=r"$\Delta t = $%i s"%delta_t)

    ax[0].set_ylabel('Attenuation factor')

    for i in range(2):
        ax[i].set_ylim(0, 1.1)
        ax[i].set_xlabel('Offset from phase center (deg)')
        ax[i].grid(color='silver', alpha=0.5)
        ax[i].legend()

    #plt.show()
    plt.savefig('smaering1.pdf')
    plt.close()

    if args.offset:
        print ('\nSelected offset =', args.offset, 'deg')

        print ('-> FREQUENCY SMEARING')
        for delta_f in args.freq:
            myRb = freq_smearing(args.freq_0, delta_f, args.bsl, math.radians(args.offset))
            print ('Delta f = %s kHz'%str(delta_f/1.e3), '-> attenuation factor = %.3f '%myRb, '(peak brightness reduced by %.1f%%)'% (np.abs((1-myRb)*100)))

        print ('-> TIME SMEARING')
        for delta_t in args.time:
            myRa = time_smearing(args.freq_0, delta_t, args.bsl, math.radians(args.offset))
            print ('Delta t = %s s'%str(delta_t), '-> attenuation factor = %.3f '%myRa, '(peak brightness reduced by %.1f%%)'% (np.abs((1-myRa)*100)))


if __name__ == '__main__':
    main(sys.argv[1:]) 
