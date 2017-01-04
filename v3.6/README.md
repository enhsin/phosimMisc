## Running PhoSim on NERSC

#### dependencies
```
module load python/2.7.9   (can add this to ~/.bashrc.ext)
scp -r username@conte.rcac.purdue.edu:/depot/lsst/apps/zlib-1.2.8 .
scp -r username@conte.rcac.purdue.edu:/depot/lsst/apps/include .    (cfitsio and fftw)
scp -r username@conte.rcac.purdue.edu:/depot/lsst/apps/lib .
```
### install phosim
```
git clone https://bitbucket.org/phosim/phosim_release.git
```
### generate flats
```
python phosim.py examples/flats/flatu_instcat_0 -t 8 -g condor --checkpoint=0
python tools/cluster_submit.py dag_99992000.dag -w $SCRATCH/phosim_release/work -o $SCRATCH/phosim_release/output  (need to specify the full path)
```

### SBATCH setting
```
#SBATCH -L SCRATCH
#SBATCH -p shared
#SBATCH -t 24:00:00   (single r-flat needs ~ 14 hrs with 8 cores)
#SBATCH -N 1
#SBATCH -n 8 
#SBATCH --mem=10GB
```
