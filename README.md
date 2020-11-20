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
-  [python-casacore](https://pypi.org/project/python-casacore/) 2.1.0 or later

These can be installed on Ubuntu with
```
pip install astropy matplotlib numpy python-casacore
```
If you do not need all the scripts, you can check what are the required packages in each one. The scripts are developed using Python 3.8, but some of them might also work with Python 2.7 &ndash; verify by checking the scripts description.

## Scripts description
-  `flagManager.py`: it saves or replaces the FLAG column in Measurement Set files using a numpy array (Python 2.7 or later).

-  `plot_fits.py`: it generates plots of a fits image using WCS coordinates. The output image can be saved in pdf or png format. Other plotting options can be listed with `python plot_fits.py -h` (Python 3.8).

