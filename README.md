# FusionScan
Gene fusion detection and visualization

# Download
Get latest (may be not stable)
```shell
# download use http
https://github.com/OpenGene/FusionScan/archive/master.zip

# or download use git
git clone https://github.com/OpenGene/FusionScan.git
```
Get the stable releases  
https://github.com/OpenGene/FusionScan/releases/latest

# Build
FusionScan only depends on `libz`, which is always available on Linux or Mac systems. If your system has no `libz`, install it first.
```shell
cd FusionScan
make
```

# Usage
You should input following files to run fusionscan
* the reference genome fasta file, specified by `-r` or `--ref=`
* the fusion setting file, specified by `-f` or `--fusion=`
* the fastq file(s), specified by `-1` or `--read1=` for single-end data. If dealing with pair-end data, specify the read2 file by `-2` or `--read2=`
* use `-h` or `--html=` to specify the file name of HTML report
* the plain text resultis directly printed to STDOUT, you can use `>` to pipe it to a file
```shell
fusionscan -r hg19.fasta -f genes/cancer.hg19.csv -1 R1.fq -2 R2.fq -h report.html > result
```
