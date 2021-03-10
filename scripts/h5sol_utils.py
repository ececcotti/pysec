import optparse
import os,sys
import numpy as np 
import h5py
import matplotlib.pyplot as plt

o = optparse.OptionParser(usage='%prog [options] *.h5')

opts,args = o.parse_args(sys.argv[1:])

solfile = h5py.File(args[0], 'r')
sols = solfile['sol000']

amp = sols['amplitude000']['val']
phase = sols['phase000']['val']

frq_range = [np.min(sols['amplitude000']['freq'])/1.e6, np.max(sols['amplitude000']['freq'])/1.e6]
nt = amp.shape[0]
nf = amp.shape[1]
nant = amp.shape[2]
npol = amp.shape[3]


again = np.zeros((nant, npol), dtype=np.complex128)
#phase_avg = np.zeros((nant, npol), dtype=np.float64)
for a in range(nant):
    for p in range(npol):
        again[a,p] = np.nanmean(amp[:,:,a,p]) * np.exp(1j * np.nanmean(phase[:,:,a,p]))

print (solfile)


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
            blgain[bcount, :] = again[i,:] * np.conj(again[j,:])

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
plt.legend()
plt.tight_layout()
plt.savefig('bslgains_%i-%iMHz.pdf'%(frq_range[0], frq_range[1]))
plt.close()


