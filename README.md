# astro-pysEC

A collection of Python scripts useful for astronomical and radio-interferometric purposes.

## Installation
The repository can be downloaded with
```
git clone https://github.com/ececcotti/astro-pysEC
```
The following Python packages are required to run all the scripts:
-  [Astropy](https://docs.astropy.org/en/stable/index.html) 4.1 or later
-  [Matplotlib](https://matplotlib.org/) 2.0 or later
-  [NumPy](https://numpy.org/) 1.16 or later
-  [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html) 1.1.4 or later
-  [python-casacore](https://pypi.org/project/python-casacore/) 2.1.0 or later
-  [regions](https://pypi.org/project/regions/) 0.4 or later

These can be installed on Ubuntu with
```
pip install astropy matplotlib numpy python-casacore pandas regions
```
If you do not need all the scripts, you can check what are the required packages in each one. The scripts are developed using Python 3.8, but some of them might also work with Python 2.7 &ndash; verify by checking the scripts description.

## Scripts description
-  `flagManager.py`: it saves or replaces the FLAG column in Measurement Set files using a numpy array (Python 2.7 or later).

-  `plot_fits.py`: it generates plots of a FITS image using WCS coordinates. The output image can be saved in PDF or PNG format. Other plotting options can be listed with `python plot_fits.py -h` (Python 3.8 or later).

-  `reg_info.py`: it gets useful information (e.g. sum, mean and rms) from one or more regions of an image. Multiple FITS images are accepted, but only one region file. However, it can be a collection of more than one region; in this case, individual and global information will be given. DS9 (REG), CASA (CRTF) and FITS region file formats are all allowed (Python 3.8 or later).
