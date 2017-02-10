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
#SBATCH -p shared
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --mem=10GB
#SBATCH -A yourAccountNumber  
#SBATCH -n 8
#SBATCH -J raytrace
#SBATCH -o raytrace_99999999_R22_S11_E000.log
#SBATCH -d afterok:3846182

cd yourWorkDir
cat trimcatalog_99999999_R22_S11.pars >> raytrace_99999999_R22_S11_E000_0.pars
yourPhosimDir/bin/raytrace < raytrace_99999999_R22_S11_E000_0.pars 
```

### PBS setting
```
#!/bin/bash -l
#PBS -q standby
#PBS -l walltime=4:00:00
#PBS -l mem=10GB
#PBS -l naccesspolicy=shared
#PBS -l nodes=1:ppn=4
#PBS -W depend=afterok:5074467

cd yourWorkDir
yourPhosimDir/bin/raytrace < raytrace_99992006_R33_S22_E001_1.pars \
> yourWorkDir/logs/raytrace_99992006_R33_S22_E001.log
```
