[![install with conda](
https://anaconda.org/bioconda/genefuse/badges/version.svg)](https://anaconda.org/bioconda/genefuse)
# GeneFuse
A tool to detect and visualize target gene fusions by scanning FASTQ files directly. This tool accepts FASTQ files and reference genome as input, and outputs detected fusion results in TEXT, JSON and HTML formats.

# Take a quick glance of the informative report
* Sample HTML report: http://opengene.org/GeneFuse/report.html
* Sample JSON report: http://opengene.org/GeneFuse/report.json
* Dataset for testing: http://opengene.org/dataset.html  Please download the paired-end FASTQ files for GeneFuse testing (Illumina platform)

# Get genefuse program
## install with Bioconda
[![install with conda](
https://anaconda.org/bioconda/genefuse/badges/version.svg)](https://anaconda.org/bioconda/genefuse)
```shell
conda install -c bioconda genefuse
```
## download binary
This binary is only for Linux systems, http://opengene.org/GeneFuse/genefuse
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/GeneFuse/genefuse
chmod a+x ./genefuse
```
## or compile from source
```shell
# get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/genefuse.git

# build
cd genefuse
make

# Install
sudo make install
```

# Usage
You should provide following arguments to run genefuse
* the reference genome fasta file, specified by `-r` or `--ref=`
* the fusion setting file, specified by `-f` or `--fusion=`
* the fastq file(s), specified by `-1` or `--read1=` for single-end data. If dealing with pair-end data, specify the read2 file by `-2` or `--read2=`
* use `-h` or `--html=` to specify the file name of HTML report
* use `-j` or `--json=` to specify the file name of JSON report
* the plain text result is directly printed to STDOUT, you can pipe it to a file using a `>`

## Example
```shell
genefuse -r hg19.fasta -f genes/druggable.hg19.csv -1 genefuse.R1.fq.gz -2 genefuse.R2.fq.gz -h report.html > result
```

## Reference genome
The reference genome should be a single whole FASTA file containg all chromosome data. This file shouldn't be compressed. For human data, typicall `hg19/GRch37` or `hg38/GRch38` assembly is used, which can be downloaded from following sites:
* `hg19/GRch37`: ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
* `hg38/GRch38`: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz  Remember to decompress hg38.fa.gz since it is gzipped and is not supported currently.

## Fusion file
The fusion file is a list of coordinated target genes together with their exons. A sample is:
```CSV
>EML4_ENST00000318522.5,chr2:42396490-42559688
1,42396490,42396776
2,42472645,42472827
3,42483641,42483770
4,42488261,42488434
5,42490318,42490446
...

>ALK_ENST00000389048.3,chr2:29415640-30144432
1,30142859,30144432
2,29940444,29940563
3,29917716,29917880
4,29754781,29754982
5,29606598,29606725
...
```
The coordination system should be consistent with the reference genome.  
### Fusion files provided in this package
Four fusion files are provided with `genefuse`:
1. `genes/druggable.hg19.csv`: all druggable fusion genes based on `hg19/GRch37` reference assembly.
2. `genes/druggable.hg38.csv`: all druggable fusion genes based on `hg38/GRch38` reference assembly.
3. `genes/cancer.hg19.csv`: all COSMIC curated fusion genes (http://cancer.sanger.ac.uk/cosmic/fusion) based on `hg19/GRch37` reference assembly.
4. `genes/cancer.hg38.csv`: all COSMIC curated fusion genes (http://cancer.sanger.ac.uk/cosmic/fusion) based on `hg38/GRch38` reference assembly.

Notes:
* `genefuse` runs almost ~5X faster with `druggable` genes than `cancer` genes, since `druggable` genes are only a small subset of `cancer` genes. Use this one if you only care about the fusion related personalized medicine for cancers.
* The `cancer` genes should be enough for most cancer related studies, since all COSMIC curated fusion genes are included.
* If you want to create a custom gene list, please follow the instructions given on next section.
### Create a fusion file based on hg19 or hg38
If you'd like to create a custom fusion file, you can use `scripts/make_fusion_genes.py`   
As the script uses `refFlat.txt` file to determine genomic coordinates of exons, you need to download a `refFlat.txt` file from UCSC Genome Browser in advance. Of course, the choice of using either hg19 or hg38 is up to you.

As for the input gene list file, all genes should be listed in separate lines.  By default, the longest transcript will be used. However, you can specify a different transcript by adding the transcript ID to the end of a gene. The gene and its transcript should be separated by a tab or a space. Please note that each gene should be the HGNC official gene symbol, and each transcript should be NCBI RefSeq transcript ID. 

An example of gene list file:

```
BRCA2	NM_000059
FAM155A
IRS2
```

When both input gene list file (`gene_list.txt`) and `refFlat.txt` file are prepared, you can use following command to generate a user-defined fusion file (`fusion.csv`):

```shell
python3 scripts/make_fusion_genes.py gene_list.txt -r /path/to/refflat -o fusion.csv
```

# HTML report
GeneFuse can generate very informative and interactive HTML pages to visualize the fusions with following information:
* the fusion genes, along with their transcripts.
* the inferred break point with reference genome coordinations.
* the inferred fusion protein, with all exons and the transcription direction.
* the supporting reads, with all bases colorized according to their quality scores.
* the number of supporting reads, and how many of them are unique (the rest may be duplications)
## A HTML report example
![image](http://www.opengene.org/GeneFuse/eml4alk.png)  
See the HTML page of this picture: http://opengene.org/GeneFuse/report.html

# All options
```
options:
  -1, --read1       read1 file name (string)
  -2, --read2       read2 file name (string [=])
  -f, --fusion      fusion file name, in CSV format (string)
  -r, --ref         reference fasta file name (string)
  -u, --unique      specify the least supporting read number is required to report a fusion, default is 2 (int [=2])
  -d, --deletion    specify the least deletion length of a intra-gene deletion to report, default is 50 (int [=50])
  -h, --html        file name to store HTML report, default is genefuse.html (string [=genefuse.html])
  -j, --json        file name to store JSON report, default is genefuse.json (string [=genefuse.json])
  -t, --thread      worker thread number, default is 4 (int [=4])
  -?, --help        print this message
```

# Cite GeneFuse
If you used GeneFuse in you work, you can cite it as: 

Shifu Chen, Ming Liu, Tanxiao Huang, Wenting Liao, Mingyan Xu and Jia Gu. GeneFuse: detection and visualization of target gene fusions from DNA sequencing data. International Journal of Biological Sciences, 2018; 14(8): 843-848. doi: 10.7150/ijbs.24626
