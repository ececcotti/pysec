#! /usr/bin/env python

import os,sys
import numpy as np 
from argparse import ArgumentParser
from casacore import tables as tb 

parser = ArgumentParser(description="Flag channels in MS files")
parser.add_argument('ms_files', nargs='+', help='input MS files')
parser.add_argument('-c', nargs='+', dest="ch", type=str, help="channels to be flagged (more than one channel allowed). If you want to flag a range of channels, use '-' to separate values. Example: '-c 1 2 10-12' will flag channels 1, 2, 10, 11, and 12", def
ault=None)
parser.add_argument('-f', nargs=2, dest="freqs", metavar=('<f1>', '<f2>'), type=float, help="frequency range (in MHz) to be flagged. Using a range helps the flagging to find channels, but if you want to flag a single channel and exactly know it frequency val
ue, you can set the same value for <f1> and <f2>", default=None)


def main(argv):
    args=parser.parse_args(argv)

    for ms in args.ms_files:
        with tb.table(ms + '::SPECTRAL_WINDOW') as spw:
            all_freq = spw.getcol('CHAN_FREQ')[0,:] / 1.e6
        

        if args.ch:
            all_chs = [x.strip() for x in args.ch]
            chs = []
            for ch in all_chs:
                ch = ch.split('-')
                if len(ch) > 1: 
                    chs_range = np.arange(int(ch[0]), int(ch[1])+1, 1)
                    chs_range = list(chs_range)
                    for i in chs_range: chs.append(i)
                else: chs.append(int(ch[0]))

        if args.freqs:
            if args.freqs[0] == args.freqs[1]:
                chs = np.where(all_freq == args.freqs[0])[0]

            else:
                chs = np.where((all_freq >= args.freqs[0]) & (all_freq <= args.freqs[1]))[0]
            
        print ('selected channels:')
        print (chs)

        with tb.table(ms, readonly=False) as t:
            flags = t.getcol('FLAG')

            #print(flags.shape, flags[:,chs,:])
            flags[:,chs,:] = True
            
            t.putcol('FLAG', flags)

        






if __name__ == '__main__':
    main(sys.argv[1:])
