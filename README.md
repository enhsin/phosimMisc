# PhoSimMisc

The Photon Simulator ([PhoSim](https://bitbucket.org/phosim/phosim_release)) is a set of fast photon Monte Carlo codes used to calculate the physics of the atmosphere and a telescope & camera in order to simulate realistic astronomical images. This repository contains miscellaneous scripts to run or validate PhoSim.

* opd/ - validating PhoSim optical path difference (OPD) with Zemax.
* opd_1.2_0.8/ - validating OPD at field angle = (1.2, 0.8). Testing against various n_silica.
* sensitivityMatrix/ - making sensitivity matrix.
* xyPosition.py - conversion between sky and image coordinates
* ZMXtoPhosim.py - convert Zemax file to PhoSim optics_x.txt
* refractive.py - calculate refractive index of silica with Sellmeier equation

