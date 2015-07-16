#PhoSim OPD maps for sensitivity matrix

Example:

```
python runOPD.py -s 'M2' -d 0.2 -n 6 --zernike
```

This adds a 0.2 micron perturbation to M2, a Noll index = 6 Zernike polynomial (vertical astigmatism).

```
python runOPD.py -k 1 -c 0 -r 4
```

This generates an opd map for the 4th row and the 0th column (indices start from 0) element in linearity_table_bending_short.txt (M2 -0.5 mm y-decenter) for field point #1.  The field angle for all field points can be found in [fieldpoints.png](https://github.com/enhsin/phosimMisc/blob/master/sensitivityMatrix/fieldpoints.png) (from Bo).

```
python runOPD.py -r -1 --all
```

This generates opd maps for all field points under perfect optics (set row to -1).

```
python runOPD.py -k 1 --fea='M1_b20_-0.50_grid.DAT' -s 'M1'
```

M1 bending data 'M1_b20_-0.50_grid.DAT' for field point #1.

```
python runOPD.py -k 1 --fea='M2_b1_-0.25_gridz.txt,M2_b1_-0.25_grid.DAT' -s 'M2'
```

M2 bending data 'M2_b1_-0.25_grid.DAT' and 1-28 Zernikes 'M2_b1_-0.25_gridz.txt' for field point #1. Input files are separated by comma. Zernike coefficient file should be listed first.


