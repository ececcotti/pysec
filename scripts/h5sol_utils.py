import optparse
import os,sys
import numpy as np 
import h5py
import matplotlib.pyplot as plt

def freqUnit(fu='MHz'):
    fu = fu.lower()
    if fu=='ghz': return 1.e9
    elif fu=='mhz': return 1.e6
    elif fu=='khz': return 1.e3
    elif fu=='hz': return 1.
    else: raise ValueError("Invalid frequency unit (only GHz, MHz, kHz, Hz accepted)")

def ant2coord(ants):
    ant_coord = []
    for i in range(ants.shape[0]):
        ant_coord.append(ants[i][1])

    ant_coord = np.array(ant_coord)
    ant_coord -= ant_coord[0]

    return ant_coord

def antennaGains(amp, phase, avg_time=False, avg_frequency=False, avg_antenna=False, avg_polarization=False):

    avg = np.array([avg_time, avg_frequency, avg_antenna, avg_polarization])

    if np.any(avg) == False:
        return amp[:] * np.exp(1j * phase[:])
    else: 
        wTrue = tuple(np.where(avg==True)[0])
        #return np.nanmean(amp, axis=wTrue) * np.exp(1j * np.nanmean(phase, axis=wTrue))
        return np.nanmean(amp[:] * np.exp(1j * phase[:]), axis=wTrue)

def baselineGains(ant_gains, ant_coord, ant_axis, avg_time=False, avg_frequency=False, avg_antenna=False, avg_polarization=False, nolofar=False):

    nant = ant_coord.shape[0]

    bll = np.zeros(int(nant*(nant-1)/2), dtype=np.float64)
    #blgain = np.zeros((int(nant*(nant-1)/2), npol), dtype=np.complex128)

    avg = np.array([avg_time, avg_frequency, avg_antenna, avg_polarization])

    ant_count = 0
    bl_count = 0
    for i in range(nant):
        for j in range(ant_count,nant):
            if i != j: 
                bll[bl_count] = np.sqrt(np.sum(np.square(ant_coord[i] - ant_coord[j]))) 
                if np.any(avg) == False:
                    blgain.append(again[i,:] * np.conj(again[j,:])
                bl_count += 1
        ant_count +=1
#    w_cc = []; w_cr = []; w_rr = []
#    acount = 0
#    bcount = 0
#    for i in range(nant):
#        for j in range(acount,nant):
#            if i != j: 
#                bll[bcount] = np.sqrt(np.sum(np.square(ants_coord[i] - ants_coord[j]))) 
#                #blgain[bcount, :] = again[i,:] * np.conj(again[j,:])
#                blgain[bcount, :] = np.nanmean(again[:,:,i,:] * np.conj(again[:,:,j,:]), axis=(0,1))

#                # index of CS-CS, RS-RS, CS-RS baselines
#                st1 = ants[i][0].decode('utf-8')[0]
#                st2 = ants[j][0].decode('utf-8')[0]
#                if st1 == 'C' and st2 == 'C': w_cc.append(bcount)
#                elif st1 == 'R' and st2 == 'R': w_rr.append(bcount)
#                else: w_cr.append(bcount)

#                bcount += 1
#        acount += 1

    w_cc = np.asarray(w_cc)
    w_rr = np.asarray(w_rr)
    w_cr = np.asarray(w_cr)


if __name__=="__main__":

    o = optparse.OptionParser(usage='%prog [options] *.h5')
    o.add_option('--nsol', nargs=2, dest='nsol', default=('000','000'), type="str", help='Numbers of solution set and subtable (default: 000 000)')
    o.add_option('--frqu', dest='frqu', default='MHz', type="str", help='Frequency unit for plot (default: MHz)')
    o.add_option('--nolofar', dest='nolofar', default=True, action='store_false', help='Disable functionalities which work only with LOFAR data')

    opts,args = o.parse_args(sys.argv[1:])

    solfile = h5py.File(args[0], 'r')
    sols = solfile['sol' + opts.nsol[0]]

    amp = sols['amplitude' + opts.nsol[1]]['val']
    phase = sols['phase' + opts.nsol[1]]['val']

    frqu = freqUnit(opts.frqu)

    frq_range = [np.min(sols['amplitude' + opts.nsol[1]]['freq'])/frqu, np.max(sols['amplitude' + opts.nsol[1]]['freq'])/frqu]
    nt = amp.shape[0]
    nf = amp.shape[1]
    nant = amp.shape[2]
    npol = amp.shape[3]

    ant_gains = antennaGains(amp,phase)
    ant_coord = ant2coord(sols['antenna'])

    """
    ants = sols['antenna']
    ants_coord = []
    for i in range(ants.shape[0]):
        ants_coord.append(ants[i][1])

    ants_coord = np.array(ants_coord)
    ants_coord -= ants_coord[0]

    bll = np.zeros(int(nant*(nant-1)/2), dtype=np.float64)
    blgain = np.zeros((int(nant*(nant-1)/2), npol), dtype=np.complex128)

    w_cc = []; w_cr = []; w_rr = []
    acount = 0
    bcount = 0
    for i in range(nant):
        for j in range(acount,nant):
            if i != j: 
                bll[bcount] = np.sqrt(np.sum(np.square(ants_coord[i] - ants_coord[j]))) 
                #blgain[bcount, :] = again[i,:] * np.conj(again[j,:])
                blgain[bcount, :] = np.nanmean(again[:,:,i,:] * np.conj(again[:,:,j,:]), axis=(0,1))

                # index of CS-CS, RS-RS, CS-RS baselines
                st1 = ants[i][0].decode('utf-8')[0]
                st2 = ants[j][0].decode('utf-8')[0]
                if st1 == 'C' and st2 == 'C': w_cc.append(bcount)
                elif st1 == 'R' and st2 == 'R': w_rr.append(bcount)
                else: w_cr.append(bcount)

                bcount += 1
        acount += 1

    w_cc = np.asarray(w_cc)
    w_rr = np.asarray(w_rr)
    w_cr = np.asarray(w_cr)
    
    fitdeg = 4
    fitp = np.polynomial.polynomial.Polynomial.fit(x=bll, y=np.abs(blgain[:,0]), deg=fitdeg)
    xfit = np.arange(np.min(bll), np.max(bll), np.max(bll)/100)


    plt.figure(figsize=(8,5))
    plt.plot(bll[w_cc], np.abs(blgain[w_cc,0]), '.', mew=0, label='CS-CS baseline', alpha=0.75)
    plt.plot(bll[w_rr], np.abs(blgain[w_rr,0]), '.', mew=0, label='RS-RS baseline', alpha=0.75)
    plt.plot(bll[w_cr], np.abs(blgain[w_cr,0]), '.', mew=0, label='CS-RS baseline', alpha=0.75)
    #plt.plot(xfit, fitp(xfit), '-', color='C3', label='polyfit deg=4')
    plt.axvline(5000, ls='--', color='C3')
    plt.xlabel('Baseline length (m)')
    plt.ylabel('Gain amplitude')
    plt.title('solutions at %.1f-%.1f MHz' % (frq_range[0], frq_range[1]))
    plt.grid(color='silver')
    plt.legend()
    plt.tight_layout()
    plt.savefig('bslgainsss_%i-%iMHz.pdf'%(frq_range[0], frq_range[1]))
    plt.close()
    """