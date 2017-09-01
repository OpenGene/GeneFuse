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
You should provide following arguments to run fusionscan
* the reference genome fasta file, specified by `-r` or `--ref=`
* the fusion setting file, specified by `-f` or `--fusion=`
* the fastq file(s), specified by `-1` or `--read1=` for single-end data. If dealing with pair-end data, specify the read2 file by `-2` or `--read2=`
* use `-h` or `--html=` to specify the file name of HTML report
* the plain text result is directly printed to STDOUT, you can pipe it to a file using a `>`

## Example
```shell
fusionscan -r hg19.fasta -f genes/cancer.hg19.csv -1 R1.fq -2 R2.fq -h report.html > result
```

## Reference genome
The reference genome should be a single whole FASTA file containg all chromosome data. This file shouldn't be compressed. For human data, typicall `hg19/GRch37` or `hg38/GRch38` assembly is used, which can be downloaded from following sites:
* `hg19/GRch37`: ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
* `hg38/GRch38`: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
