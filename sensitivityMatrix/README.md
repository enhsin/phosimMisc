#PhoSim OPD maps for sensitivity matrix

Example:

```
python runOPD.py -s 'M2' -d 0.2 -n 6 --zernike
```

This adds a 0.2 micron perturbation to M2, characterized by a Noll index = 6 Zernike polynomial (vertical astigmatism).

```
python runOPD.py -k 1 -c 0 -r 4
```

This generates an opd map for the 4th row and the 0th column (indices start from 0) element in linearity_table_bending_short.txt (M2 -0.5 mm y-decenter) for field point #1.  The field angle for #1-#35 can be found in [fieldpoints.png](https://github.com/enhsin/phosimMisc/blob/master/sensitivityMatrix/fieldpoints.png) from Bo Xin.

```
python runOPD.py -r -1 --all
```

This generates opd maps for all field points under perfect optics (set row to -1).
