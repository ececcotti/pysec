from casacore import tables as tb
import optparse
import os,sys
import numpy as np
import matplotlib.pyplot as plt

o = optparse.OptionParser()
o.set_usage('python data_inspection.py [options] MS_file')
o.set_description(__doc__)
o.add_option('--info', dest='info', action='store_true', help='Print MS info')
o.add_option('--bsl', dest='bsl', action='store_true', help='Print baselines info')
o.add_option('--bsort', dest='bsort', action='store_true', help='Sort rows by baseline')
o.add_option('--layout', dest='layout', action='store_true', help='Plot antenna layout')
o.add_option('--core', dest='core', action='store_true', help='Plot only core station layout (required --layout option)')

opts, args = o.parse_args(sys.argv[1:])

if opts.core and not opts.layout:
        print ('WARNING: You have to select --layout option too to plot core station layout')
        print ('--> continue without plotting antenna layout \n')


for filename in args:
        print ('Reading ',filename)
        path,srcFile = os.path.split(os.path.realpath(filename))
        filename_split = filename.split('/')
        filename_ms = filename.split('/')[-1]

        if opts.info:
                os.system('msoverview in=' + filename + ' verbose=T')
                print ('--> MS file structure')
                print (tb.tablestructure(filename))

        t = tb.table(filename)
        uvw = t.getcol('UVW')

        if opts.bsl:
                u = uvw[:,0]
                v = uvw[:,1]
                w = uvw[:,2]
                uvLength = np.sqrt(u**2 + v**2 + w**2)
                max_bsl = np.max(uvLength)
                min_bsl = np.min(uvLength[uvLength > 0])
                print ('Maximum baseline length =', max_bsl, 'm')
                print ('Minimum baseline length =', min_bsl, 'm')

                usf_bsl_0 = uvLength[uvLength <= (min_bsl + 0.8)]

        if opts.layout:
                t_ant = tb.table(filename + '/ANTENNA')[:]
                ant_pos = []
                for ant in t_ant:
                        ant_name = ant['NAME']
                        if opts.core and str(ant_name[0:2]) == 'CS':
                                ant_pos_0 = ant['POSITION']
                                ant_pos.append(ant_pos_0)
                        elif not opts.core:
                                ant_pos_0 = ant['POSITION']
                                ant_pos.append(ant_pos_0)
                        print (ant_name, ant_pos_0[1], ant_pos_0[0])

                ant_pos = np.array(ant_pos)
                ant_pos_x = ant_pos[:,1]
                ant_pos_y = -ant_pos[:,0]
                ant_pos_z = ant_pos[:,2]
                z_diff = (np.max(ant_pos_z)-np.min(ant_pos_z))/np.max(ant_pos_z) * 100
                print ('The maximum difference along the z-axis is %.2f' % z_diff + '%')
                if z_diff < 1.0: 
                        print ('--> negligible: safe 2D plot')
                else:
                        print ('--> not negligible: work in progress for 3D plot (2D plot made)')

                ant_pos_x_min = np.min(ant_pos_x)
                ant_pos_y_min = np.min(ant_pos_y)
                ant_pos_x = ant_pos_x - ant_pos_x_min
                ant_pos_y = ant_pos_y - ant_pos_y_min
                plt.plot(ant_pos_x, ant_pos_y, 'b.', mec='none', alpha=0.7)
                plt.xlabel('x position')
                plt.ylabel('y position')
                plt.xlim(-np.max(ant_pos_x)*0.02, np.max(ant_pos_x)*1.02)
                plt.ylim(-np.max(ant_pos_y)*0.02, np.max(ant_pos_y)*1.02)
                plt.show()
