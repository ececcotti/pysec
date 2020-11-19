import optparse
import os,sys
import numpy as np
from casacore import tables as tb

o = optparse.OptionParser(usage='%prog [options] *.MS')
o.add_option('-s', '--save',
            dest='save', action='store_true',  
            help='save the FLAG column as .npy file')
o.add_option('-r', '--rewrite',
            dest='rewrite', action='store_true',
            help='rewrite the FLAG column using the .npy file')
o.add_option('-f', '--fname', 
            dest='fname', type='str', default=None,
            help='name of the .npy file containing flags, with no extension (default ONLY in saving mode: MS_name.flag_v00 ; if the file already exists, a new version number will be used)')

opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    ms_ext = filename.split('.')[-1]
    ms_name = filename.split('.' + ms_ext)[-2]

    if opts.save and opts.rewrite:
        raise Exception('Only one mode at a time can be used: choose -s OR -r')
    if not opts.fname and opts.rewrite:
        raise Exception('-f/--fname must be specified in rewriting mode')

    fversion = 0
    if opts.fname == None: 
        fname_pattern = ms_name + '.flag_v'
        opts.fname = fname_pattern + '%02d'%fversion

    fname_check = os.path.exists(opts.fname + '.npy')
    while fname_check == True and opts.save:
        fversion += 1
        opts.fname = fname_pattern + '%02d'%fversion
        fname_check = os.path.exists(opts.fname + '.npy')
    
    t = tb.table(filename, readonly=opts.save)
    if opts.save: 
        print('Saving FLAG column in ' + opts.fname + '.npy')
        flag = t.getcol('FLAG')
        flagfile = open(opts.fname + '.npy','wb')
        np.save(flagfile, flag)
    if opts.rewrite:
        print('Overwriting FLAG column with ' + opts.fname + '.npy')
        flagfile = open(opts.fname + '.npy','rb')
        flag = np.load(flagfile)
        t.putcol('FLAG', flag)

    flagfile.close()
    t.close()

