1. General Information
login

```bash
haina@login.htcf.wustl.edu
```
for analysis storage:
```bash
cd lts/mtblab/baldridge1/ANALYSIS/haina
```
for scratch:
```bash
cd /scratch/mtblab/haina
```
for sequence storage:
```bash
cd lts/mtblab/baldridge1/PERMANENT_FILE_STORAGE
```
create job:

```bash
#!/bin/bash

#SBATCH --cpus-per-task=2
#SBATCH --mem=256G
#SBATCH -J  fastqc082924
#SBATCH --output=slurm-%x.%j.out  # %j gives job id, %x gives job name
#SBATCH --error=slurm-%x.%j.err  # Errors will be written in this file
#SBATCH --mail-user=haina  # Email address to receive notifications about submitted job
#SBATCH --mail-type=ALL

eval "$(conda shell.bash hook)"

ROOT_DIR="/lts/mtblab/baldridge1/PERMANENT_FILE_STORAGE/WGS/Raw_data/"
OUTPUT_DIR=$SCRATCHDIR/fastqoutput

cd $ROOT_DIR

conda activate magsprocess

# define the root firectory where raw data is located


#create the output directory is it doesnot exist
mkdir -p "$OUTPUT_DIR"

# find all .fastqc.bz2 files and run fastqc on each

find "$ROOT_DIR" -type f -name "*.fastq.bz2" | while read -r file; do
    echo "Running FastQC on $file"
    fastqc -o "$OUTPUT_DIR" "$file"

done
```

before submit yoru job
```bash
chmod +x job.bash # to make your job executable
```



