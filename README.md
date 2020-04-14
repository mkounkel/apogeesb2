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
positional arguments:
  directory             Directory containing apstar fits files

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             Name of the fits table to save identified SB2
                        properties (Default sb2s.fits)
  --saveall SAVEALL     Save deconvolution for all sources (Default False)
  --outall OUTALL       Name of the fits table to save identified SB2
                        properties (Default all_deconvolutions.fits)
  --ccfs CCFS           Name where to dump the pickle file containing ccfs
                        (Default ccfs.pickle)
  --meta META           Name where to dump the pickle file containing
                        corresponding metadata (Default ccfs_meta.pickle)
  --deconvol DECONVOL   Name where to dump the pickle file containing raw
                        deconvolution (Default ccfs_decomposed.pickle)
  --deletetemp DELETETEMP
                        Delete temporary files (Default True)
```
