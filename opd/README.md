# PhoSim OPD Validation

```
python checkRayPosition.py output/output.fits.gz zemax/chiefRay.txt zemax/dbg.txt
```

This generates a plot [position.pdf](https://github.com/enhsin/phosim_tool/blob/master/opd/position.pdf) showing the difference of the ray position at each optical surface from PhoSim and Zemax. The field angle is 1.7 degree and the wavelength is 770nm. The general ray ID is 80*OPD_SCREEN_SIZE + 15.

```
python checkRayPosition.py output/output127.fits.gz lsst/optics_3.txt
```

This generates a plot [layout.pdf](https://github.com/enhsin/phosim_tool/blob/master/opd/layout.pdf) showing the telescope layout including the entrance pupil and the reference sphere. The ID of this particular general ray is 15*OPD_SCREEN_SIZE + 127.

```
python plotOPD.py output/opd.fits.gz zemax/OPD_0_1U.TXT
```

This generates [opd.pdf](https://github.com/enhsin/phosim_tool/blob/master/opd/opd.pdf) and [diff.fits.gz](https://github.com/enhsin/phosim_tool/blob/master/opd/diff.fits.gz). 'opd.pdf' plots PhoSim and Zemax OPD maps and their differences. 'diff.fits.gz' is the difference in wave numbers.


* input/ - input catalog and command files for PhoSim

* output/ - PhoSim opd map and event file (storing ray positions)

* zemax/ - Zemax results from Bo Xin

* SEDs/, lsst/, raytrace/ - relevant source and data files for this simulation

