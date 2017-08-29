# FusionScan
Being coded... NOT ready to be used

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
```shell
fusionscan -r hg19.fasta -f genes/cancer.hg19.csv -1 R1.fq -2 R2.fq > result
```
