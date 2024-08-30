sequence cleaning info

sequence are quality controlled using fastqc (FastQC v0.12.1) and multiqc (version 1.24.1)

```bash
find "$ROOT_DIR" -type f -name "*.fastq.bz2" | while read -r file; do
    echo "Running FastQC on $file"
    fastqc -o "$OUTPUT_DIR" "$file" -t 4
```

demultiplexed raw fastq sequences are processed using BBDuk (sourceforge.net/projects/bbmap/; BBMap version bbmap-39.08-0) to quality trim, remove
Illumina adapters. Trimming parameters are set to a k-mer length of 19 and a minimum Phred quality score of 25. Reads with a minimum average Phred quality score below 23 and
length shorter than 50 bp after trimming are discarded. 

```bash
bbduk.sh in="$FILE" out="$OUTPUT_FILE" ref="$ADAPTER_PATH" ktrim=r k=19 mink=11 hdist=1 qtrim=rl trimq=25 maq=23 minlength=50
```
