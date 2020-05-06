# apogeesb2
Examine CCFs in the APOGEE apstar files to search for double lined spectroscopic binaries

## Required packages:
* Astropy
* Gausspy (https://github.com/gausspy/gausspy)

## Installation
python3 setup.py install

## Example
```
apogeesb2 /path/apogee/spectro/redux/r13/stars/apo25m/105-45
```

##
```
Identify SB2s in the APOGEE spectra via gaussian fitting of the CCF

positional arguments:
  directory             Directory containing apstar fits files

optional arguments:
  -h, --help            show this help message and exit
  --usepath USEPATH     Instead of scanning directory, use txt file with full
                        paths to all apstar files (Default False)
  --out OUT             Name of the fits table to save identified SB2
                        properties (Default sb2s.fits)
  --saveall SAVEALL     Save deconvolution for all sources (Default False)
  --outall OUTALL       Name of the fits table to save identified SB2
                        properties (Default all_deconvolutions.fits)
  --makeccf MAKECCF     Make CCF from the data (default True), or read
                        existing pickle
  --ccfs CCFS           Name where to dump the pickle file containing ccfs
                        (Default ccfs.pickle)
  --meta META           Name where to dump the pickle file containing
                        corresponding metadata (Default ccfs_meta.pickle)
  --deconvolve DECONVOLVE
                        Deconvolve CCFs from scratch (default True), or read
                        existing pickle
  --deconvol DECONVOL   Name where to dump the pickle file containing raw
                        deconvolution (Default ccfs_decomposed.pickle)
  --deletetemp DELETETEMP
                        Delete temporary files (Default True)
  --makeplots MAKEPLOTS
                        Generate CCF plots (Default True)
  --plotdir PLOTDIR     Folder to output the figures (Default plots)
  --width WIDTH         Width from the central peak to scan, in km/s (Default
                        400)
  --alpha ALPHA         Default alpha parameter (Default 1.5)
  --offset OFFSET       Default offset to the continuum (Default 0)
```
