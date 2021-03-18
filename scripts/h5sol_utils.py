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

def antennaGains(amp, phase, avg_time=False, avg_frequency=False, avg_polarization=False):

    ant_axis = 2
    avg_antenna = False
    avg = np.array([avg_time, avg_frequency, avg_antenna, avg_polarization])

    if np.any(avg) == False:
        return amp[:] * np.exp(1j * phase[:]), ant_axis
    else: 
        wTrue = np.where(avg==True)[0]
        if wTrue.size == 1:
            if wTrue[0] == 0 or wTrue[0] == 1: ant_axis =1
        elif wTrue.size > 1:
            if wTrue[0] == 0 and wTrue[1] == 1: ant_axis = 0
        return np.nanmean(amp, axis=tuple(wTrue)) * np.exp(-1j * np.nanmean(phase, axis=tuple(wTrue))), ant_axis
        #return np.nanmean(amp[:] * np.exp(1j * phase[:]), axis=wTrue)

def baselineGains(ant_gains, ant_coord, ant_axis, avg_time=False, avg_frequency=False, avg_polarization=False, ants=None):

    nant = ant_coord.shape[0]

    bll = np.zeros(int(nant*(nant-1)/2), dtype=np.float64)

    shape_blgain = []
    for i in range(len(ant_gains.shape)):
        if i == ant_axis: # substitute the antenna axis with the baseline axis N(N-1)/2
            shape_blgain.append(int(nant*(nant-1)/2.))
        else: shape_blgain.append(ant_gains.shape[i])
    
    bl_gains = np.zeros(shape_blgain, dtype=np.complex128)

    avg_antenna = False
    avg = np.array([avg_time, avg_frequency, avg_antenna, avg_polarization])

    w_cc = []; w_cr = []; w_rr = []
    ant_count = 0
    bl_count = 0
    slc = [slice(None)] * bl_gains.ndim
    for i in range(ant_gains.shape[ant_axis]):
        for j in range(ant_count,ant_gains.shape[ant_axis]):
            if i != j: 
                bll[bl_count] = np.sqrt(np.sum(np.square(ant_coord[i] - ant_coord[j]))) 
                if np.any(avg) == False:
                    slc[ant_axis] = bl_count
                    bl_gains[tuple(slc)] = ant_gains.take(i, axis=ant_axis) * ant_gains.take(j, axis=ant_axis).conj()

                st1 = ants[i][0].decode('utf-8')[0]
                st2 = ants[j][0].decode('utf-8')[0]
                if st1 == 'C' and st2 == 'C': w_cc.append(bl_count)
                elif st1 == 'R' and st2 == 'R': w_rr.append(bl_count)
                else: w_cr.append(bl_count)
        
                bl_count += 1
        ant_count += 1

    w_cc = np.asarray(w_cc)
    w_rr = np.asarray(w_rr)
    w_cr = np.asarray(w_cr)

    return bll, bl_gains, [w_cc, w_rr, w_cr]


if __name__=="__main__":

    o = optparse.OptionParser(usage='%prog [options] *.h5')
    o.add_option('--nsol', nargs=2, dest='nsol', default=('000','000'), type="str", help='Numbers of solution set and subtable (default: 000 000)')
    o.add_option('--frqu', dest='frqu', default='MHz', type="str", help='Frequency unit for plot (default: MHz)')
    o.add_option('--ag_avg', nargs=3, dest='ag_avg', default=(1,1,0), type="int", 
        help='Average antenna gains in time, frequency and/or polarization. Use 0 (False) and 1 (True) to average one or more axis, which are given in the previous order (default: 1 1 0)')

    opts,args = o.parse_args(sys.argv[1:])

    solfile = h5py.File(args[0], 'r')
    sols = solfile['sol' + opts.nsol[0]]

    amp = sols['amplitude' + opts.nsol[1]]['val']
    phase = sols['phase' + opts.nsol[1]]['val']

    frqu = freqUnit(opts.frqu)

    frq_range = [np.min(sols['amplitude' + opts.nsol[1]]['freq'])/frqu, np.max(sols['amplitude' + opts.nsol[1]]['freq'])/frqu]
    nt, nf, nant, npol = [amp.shape[i] for i in range(amp.ndim)]

    # antenna gains
    ant_gains, ant_axis = antennaGains(amp, phase, avg_time=opts.ag_avg[0], avg_frequency=opts.ag_avg[1], avg_polarization=opts.ag_avg[2])
    ant_coord = ant2coord(sols['antenna'])
    
    # baseline gains
    bll, bl_gains, w = baselineGains(ant_gains, ant_coord, ant_axis=ant_axis, ants=sols['antenna'])

    """
    fitdeg = 4
    fitp = np.polynomial.polynomial.Polynomial.fit(x=bll, y=np.abs(blgain[:,0]), deg=fitdeg)
    xfit = np.arange(np.min(bll), np.max(bll), np.max(bll)/100)

    """
    plt.figure(figsize=(8,5))
    plt.plot(bll[w[0]], np.abs(bl_gains[w[0],0]), '.', mew=0, label='CS-CS baseline', alpha=0.75)
    plt.plot(bll[w[1]], np.abs(bl_gains[w[1],0]), '.', mew=0, label='RS-RS baseline', alpha=0.75)
    plt.plot(bll[w[2]], np.abs(bl_gains[w[2],0]), '.', mew=0, label='CS-RS baseline', alpha=0.75)
    #plt.plot(xfit, fitp(xfit), '-', color='C3', label='polyfit deg=4')
    plt.axvline(5000, ls='--', color='C3')
    plt.xlabel('Baseline length (m)')
    plt.ylabel('Gain amplitude')
    plt.title('solutions at %.1f-%.1f MHz' % (frq_range[0], frq_range[1]))
    plt.grid(color='silver')
    plt.legend()
    plt.tight_layout()
    plt.show()
    #plt.savefig('bslgainsss_%i-%iMHz.pdf'%(frq_range[0], frq_range[1]))
    plt.close()
    