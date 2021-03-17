import numpy as np

def remove_auto_correlations(data_array, axes=(0, 1)):
    """Remove the auto-corrlation term from input array.
    Takes an N x N array, removes the diagonal components, and returns a flattened N(N-1) dimenion in place of the array.
    If uses on a M dimensional array, returns an M-1 array.
    Parameters
    ----------
    data_array : array
        Array shaped like (Nbls, Nbls, Ntimes, Nfreqs). Removes same baseline diagonal along the specifed axes.
    axes : tuple of int, length 2
        axes over which the diagonal will be removed.
    Returns
    -------
    data_out : array with the same type as `data_array`.
        (Nbls * (Nbls-1), Ntimes, Nfreqs) array.
        if input has pols: (Npols, Nbls * (Nbls -1), Ntimes, Nfreqs)
    Raises
    ------
    ValueError
        If axes is not a length 2 tuple.
        If axes are not adjecent (e.g. axes=(2,7)).
        If axes do not have the same shape.
    """
    n_inds = data_array.shape[axes[0]]
    indices = np.logical_not(np.tril(np.ones(n_inds, dtype=bool)))
    # move the axes so axes[0] is the 0th axis and axis 1 is the 1th
    data_array = np.moveaxis(data_array, axes[0], 0)
    data_array = np.moveaxis(data_array, axes[1], 1)
    data_out = data_array[indices]
    data_out = np.moveaxis(data_out, 0, axes[0])
    return data_out

def cross_multiply_array(array1,axis=0):
    """Cross multiply the arrays along the given axis.
    Cross multiplies along axis and computes array_1.conj() * array_2
    if axis has length M then a new axis of size M will be inserted directly succeeding the original.
    Parameters
    ----------
    array_1 : array_like
        N-dimensional array_like
    array_2 : array_like, optional
        N-dimenional array.
        Defaults to copy of array_1
    axis : int
        Axis along which to cross multiply
    Returns
    -------
    cross_array : array_like
        N+1 Dimensional array
    Raises
    ------
    ValueError
        If input arrays have different shapes.
    """

    cross_array = np.expand_dims(array1, axis=axis).conj() * np.expand_dims(
        array1, axis=axis + 1)
    return cross_array

a = np.array([[[[1.+1.j, 2.+2.j, 3], [1.+4.j, 3.+0.j, 6+1j], [1.,3., 0+2j]], 
    [[2.+4.j, 5.+0.j, 3], [2.+5.j, 1.+4.j, 8], [0.+3*1j,3., 4]]],[[[1.+1.j, 2.+2.j, 3], [1.+4.j, 3.+0.j, 6+1j], [1.,3., 0+2j]], 
    [[2.+4.j, 5.+0.j, 3], [2.+5.j, 1.+4.j, 8], [0.+3*1j,3., 4]]]])

aaxis=2 # antenna axis
xx = np.zeros(a.shape)
yy = []
for i in range(len(a.shape)):
    if i == aaxis: # substitute the antenna axis with the baseline axis N(N-1)/2
        print ('ciao')
        yy.append(int(a.shape[i]*(a.shape[i]-1)/2.))
    else: yy.append(a.shape[i])
print(a.shape)
print(tuple(yy)) 

xx = np.zeros(yy, dtype=np.complex128)  # generate the empty baseline gains

c = cross_multiply_array(a, axis=aaxis)

d = remove_auto_correlations(c, axes=(aaxis,aaxis+1))
print (d)

# that's work!!!
acount = 0
bcount = 0
slc = [slice(None)] * xx.ndim
print (slc)
for i in range(a.shape[aaxis]):
    for j in range(acount,a.shape[aaxis]):
        if i != j: 
            slc[aaxis] = bcount
            xx[tuple(slc)] = a.take(i, axis=aaxis) * a.take(j, axis=aaxis).conj()
        
            bcount += 1
    acount += 1

print(xx)


