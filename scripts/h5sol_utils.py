import optparse
import os,sys
import numpy as np 
import h5py
import matplotlib.pyplot as plt
import matplotlib as mpl

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
        return amp * np.exp(1j * phase), ant_axis
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

def plotBaselineGains(bl_length, bl_gains, f_range, pol=0, w=None, plot_size=(8,5), save_pdf=False, save_png=False, direction=None):

    if pol not in [0,1,2,3]: raise ValueError("Invalid polarization (only 0, 1, 2, 3 accepted)")

    title = 'solutions at %.1f-%.1f MHz' % (f_range[0], f_range[1])
    outfile = 'bslgains_%i-%iMHz' % (f_range[0], f_range[1])

    if direction != None:
        title += ', dir %s' % direction
        outfile += '-%s' % direction

    plt.figure(figsize=plot_size)
    if w != None:
        plt.plot(bll[w[0]], np.abs(bl_gains[w[0],pol]), '.', mew=0, label='CS-CS baseline', alpha=0.75)
        plt.plot(bll[w[1]], np.abs(bl_gains[w[1],pol]), '.', mew=0, label='RS-RS baseline', alpha=0.75)
        plt.plot(bll[w[2]], np.abs(bl_gains[w[2],pol]), '.', mew=0, label='CS-RS baseline', alpha=0.75)
    else: 
        plt.plot(bll, np.abs(bl_gains[:,pol]), '.', mew=0, label='CS-CS baseline', alpha=0.75)

    plt.axvline(5000, ls='--', color='C3')
    plt.xlabel('Baseline length (m)')
    plt.ylabel('Gain amplitude')
    plt.grid(color='silver')
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    if save_pdf:
        plt.savefig(outfile + '.pdf')
    if save_png:
        plt.savefig(outfile + '.png', dpi=200)

    plt.show()
    plt.close()


if __name__=="__main__":

    o = optparse.OptionParser(usage='%prog [options] *.h5')
    o.add_option('--nsol', nargs=2, dest='nsol', default=('000','000'), type="str", help='Numbers of solution set and subtable (default: 000 000)')
    o.add_option('--ddsol', dest='ddsol', default=False, action='store_true', help="Enable selection of direction subtables for direction-dependent solutions")
    o.add_option('--dir', dest='dir', default=0, type="int", help='Direction subtable (default: 0)')
    o.add_option('--ag_avg', nargs=3, dest='ag_avg', default=(1,1,0), type="int", 
        help='Average antenna gains in time, frequency and/or polarization. Use 0 (False) and 1 (True) to average one or more axis, which are given in the previous order (default: 1 1 0)')
    o.add_option('--bg_avg', nargs=3, dest='bg_avg', default=(0,0,0), type="int", 
        help='Average baseline gains in time, frequency and/or polarization. Use 0 (False) and 1 (True) to average one or more axis, which are given in the previous order (default: 0 0 0)')
    
    g = optparse.OptionGroup(o, 'Plotting Options')
    g.add_option('--frqu', dest='frqu', default='MHz', type="str", help='Frequency unit for plot (default: MHz)')
    g.add_option('--ppol', dest='ppol', default=0, type="int", help='Select polarization to plot: 0=xx, 1=xy, 2=yx, 3=yy (default: 0)')
    g.add_option('--plotb', dest='plotb', default=False, action='store_true', help="Plot baseline gains")
    g.add_option('--psize', nargs=2, dest='psize', default=(8,5), type="int", help='Size (x y) of the selected plot (default: 8 5)')
    g.add_option('--pdf', dest='save_pdf', default=False, action='store_true', help='Save plot as pdf')
    g.add_option('--png', dest='save_png', default=False, action='store_true', help='Save plot as png')
    o.add_option_group(g)

    opts,args = o.parse_args(sys.argv[1:])

    solfile = h5py.File(args[0], 'r')
    sols = solfile['sol' + opts.nsol[0]]

    amp = sols['amplitude' + opts.nsol[1]]['val']
    phase = sols['phase' + opts.nsol[1]]['val']

    frqu = freqUnit(opts.frqu)

    frq_range = [np.min(sols['amplitude' + opts.nsol[1]]['freq'])/frqu, np.max(sols['amplitude' + opts.nsol[1]]['freq'])/frqu]
    
    slc_sol = [slice(None)] * amp.ndim
    if opts.ddsol: 
        nt, nf, nant, ndir, npol = [amp.shape[i] for i in range(amp.ndim)]
        slc_sol[3] = opts.dir
        dir_name = sols['amplitude' + opts.nsol[1]]['dir'][opts.dir].decode('utf-8')[1:-1]
    else: 
        nt, nf, nant, npol = [amp.shape[i] for i in range(amp.ndim)]
        dir_name = None

    # antenna gains
    print ('Generate antenna gains')
    ant_gains, ant_axis = antennaGains(amp[tuple(slc_sol)], phase[tuple(slc_sol)], 
            avg_time=opts.ag_avg[0], avg_frequency=opts.ag_avg[1], avg_polarization=opts.ag_avg[2])
    ant_coord = ant2coord(sols['antenna'])
    
    # baseline gains
    print ('Generate baseline gains')
    bll, bl_gains, w = baselineGains(ant_gains, ant_coord, ant_axis=ant_axis, 
        avg_time=opts.bg_avg[0], avg_frequency=opts.bg_avg[1], avg_polarization=opts.bg_avg[2], ants=sols['antenna'])

    # polt baseline gains
    print ('Plot')
    plotBaselineGains(bll, bl_gains, frq_range, pol=opts.ppol, w=w, plot_size=opts.psize, save_pdf=opts.save_pdf, save_png=opts.save_png, direction=dir_name)
    
    """
    clock = sols['clock' + opts.nsol[1]]['val'][:,:,0]
    time = (sols['amplitude' + opts.nsol[1]]['time'][:]-sols['amplitude' + opts.nsol[1]]['time'][0])/3600.
    for i in range(clock.shape[1]):
        wr = sols['antenna'][i][0].decode('utf-8')[0]
        if wr == 'R': 
            firstR = i
            break

    nRstation = len(clock[0,firstR:-1])
    data = range(0, nRstation+1)
    norm = mpl.colors.Normalize(vmin=min(data), vmax=max(data), clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet_r)
    node_color = [(r, g, b) for r, g, b, a in mapper.to_rgba(data)]

    plt.figure(figsize=(10,5))
    count = 0
    for i in range(clock.shape[1]):
        wr = sols['antenna'][i][0].decode('utf-8')[0]
        if wr == 'R':
            plt.plot(time, clock[:,i]*1.e9, '-', color=node_color[count], lw=1)
            count += 1

    plt.xlabel('Time (h)')
    plt.ylabel('Clock (ns)')
    plt.title('clock delay at %.1f-%.1f MHz' % (frq_range[0], frq_range[1]))
    plt.savefig('clock_%i-%iMHz.pdf'%(frq_range[0], frq_range[1]))
    plt.close()
    """